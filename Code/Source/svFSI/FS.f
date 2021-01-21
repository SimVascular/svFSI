!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     These routines initialize function spaces.
!
!--------------------------------------------------------------------

!     Initialize function spaces over a mesh
      SUBROUTINE INITFSMSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) insd

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      IF (lM%lFib) insd = 1

      ALLOCATE(lM%fs(lM%nFs))

!     The first set of basis is inherited directly from mesh basis
      lM%fs(1)%lShpF = lM%lShpF
      lM%fs(1)%eType = lM%eType
      lM%fs(1)%eNoN  = lM%eNoN
      lM%fs(1)%nG    = lM%nG

      CALL ALLOCFS(lM%fs(1), insd)

      IF (lM%eType .NE. eType_NRB) THEN
         lM%fs(1)%w   = lM%w
         lM%fs(1)%xi  = lM%xi
         lM%fs(1)%xib = lM%xib
         lM%fs(1)%N   = lM%N
         lM%fs(1)%Nb  = lM%Nb
         lM%fs(1)%Nx  = lM%Nx
      END IF

!     Sets Taylor-Hood basis if invoked by user (fluid, ustruct, FSI)
      IF (lM%nFs .EQ. 2) THEN
!        Select Taylor-Hood element
         CALL SETTHOODFS(lM%fs(2), lM%fs(1)%eType)

!        Initialize the function space
         CALL INITFS(lM%fs(2), insd)
      END IF

      RETURN
      END SUBROUTINE INITFSMSH
!--------------------------------------------------------------------
!     Initialize function spaces over a face
      SUBROUTINE INITFSFACE(lM, lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) g, insd

      insd = nsd-1
      IF (lM%lShl) insd = insd - 1
      IF (lM%lFib) insd = 0

      lFa%nFs = lM%nFs
      ALLOCATE(lFa%fs(lFa%nFs))

!     The first set of basis is inherited directly from face basis
      lFa%fs(1)%lShpF = lM%lShpF
      lFa%fs(1)%eType = lFa%eType
      lFa%fs(1)%eNoN  = lFa%eNoN
      lFa%fs(1)%nG    = lFa%nG

      CALL ALLOCFS(lFa%fs(1), insd)

      IF (lFa%eType .NE. eType_NRB) THEN
         lFa%fs(1)%w   = lFa%w
         lFa%fs(1)%xi  = lFa%xi
         lFa%fs(1)%N   = lFa%N
         lFa%fs(1)%Nx  = lFa%Nx
         CALL GETNNBNDS(lFa%fs(1)%eType, lFa%fs(1)%eNoN, lFa%fs(1)%xib,
     2      lFa%fs(1)%Nb)
      END IF

!     Sets Taylor-Hood basis if invoked by user (fluid, ustruct, FSI)
      IF (lFa%nFs .EQ. 2) THEN
!        Select Taylor-Hood element
         CALL SETTHOODFS(lFa%fs(2), lFa%fs(1)%eType)

!        Allocate arrays
         CALL ALLOCFS(lFa%fs(2), insd)

!        Get Gauss points and shape functions
         CALL GETGIP(insd, lFa%fs(2)%eType, lFa%fs(2)%nG, lFa%fs(2)%w,
     2      lFa%fs(2)%xi)
         DO g=1, lFa%fs(2)%nG
            CALL GETGNN(insd, lFa%fs(2)%eType, lFa%fs(2)%eNoN,
     2         lFa%fs(2)%xi(:,g), lFa%fs(2)%N(:,g), lFa%fs(2)%Nx(:,:,g))
         END DO
         CALL GETNNBNDS(lFa%fs(2)%eType, lFa%fs(2)%eNoN, lFa%fs(2)%xib,
     2      lFa%fs(2)%Nb)
      END IF

      RETURN
      END SUBROUTINE INITFSFACE
!--------------------------------------------------------------------
      SUBROUTINE INITFS(fs, insd)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fsType), INTENT(INOUT) :: fs
      INTEGER(KIND=IKIND), INTENT(IN) :: insd

      INTEGER g

!     Allocate arrays
      CALL ALLOCFS(fs, insd)

!     Get Gauss points and shape functions
      CALL GETGIP(insd, fs%eType, fs%nG, fs%w, fs%xi)
      DO g=1, fs%nG
         CALL GETGNN(insd, fs%eType, fs%eNoN, fs%xi(:,g), fs%N(:,g),
     2      fs%Nx(:,:,g))
      END DO
      CALL GETNNBNDS(fs%eType, fs%eNoN, fs%xib, fs%Nb)

      RETURN
      END SUBROUTINE INITFS
!--------------------------------------------------------------------
!     Allocates arrays within the function space type. Assumes that
!     fs%eNoN and fs%nG are already defined
      SUBROUTINE ALLOCFS(fs, insd)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fsType), INTENT(INOUT) :: fs
      INTEGER(KIND=IKIND), INTENT(IN) :: insd

      INTEGER(KIND=IKIND) nG, eNoN

      nG   = fs%nG
      eNoN = fs%eNoN

      ALLOCATE(fs%w(nG), fs%xi(insd,nG), fs%xib(2,nsd), fs%N(eNoN,nG),
     2   fs%Nb(2,eNoN), fs%Nx(insd,eNoN,nG))

      IF (fs%eType .EQ. eType_NRB) THEN
         IF (insd .EQ. 1) THEN
            ALLOCATE(fs%Nxx(1,eNoN,nG))
         ELSE IF (insd .EQ. 2) THEN
            ALLOCATE(fs%Nxx(3,eNoN,nG))
         ELSE IF (insd .EQ. 3) THEN
            ALLOCATE(fs%Nxx(6,eNoN,nG))
         END IF
      END IF

      RETURN
      END SUBROUTINE ALLOCFS
!####################################################################
!     Sets Tayloor-Hood basis for a parent element type
      SUBROUTINE SETTHOODFS(fs, eType)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fsType), INTENT(INOUT) :: fs
      INTEGER(KIND=IKIND), INTENT(IN) :: eType

      SELECT CASE (eType)
      CASE (eType_TET10)
         fs%eType = eType_TET4
         fs%lShpF = .TRUE.
         fs%eNoN  = 4
         fs%nG    = 4

      CASE (eType_HEX20)
         fs%eType = eType_HEX8
         fs%lShpF = .TRUE.
         fs%eNoN  = 8
         fs%nG    = 8

      CASE (eType_HEX27)
         fs%eType = eType_HEX8
         fs%lShpF = .TRUE.
         fs%eNoN  = 8
         fs%nG    = 8

      CASE (eType_TRI6)
         fs%eType = eType_TRI3
         fs%lShpF = .TRUE.
         fs%eNoN  = 3
         fs%nG    = 3

      CASE (eType_QUD8)
         fs%eType = eType_QUD4
         fs%lShpF = .FALSE.
         fs%eNoN  = 4
         fs%nG    = 4

      CASE (eType_QUD9)
         fs%eType = eType_QUD4
         fs%lShpF = .FALSE.
         fs%eNoN  = 4
         fs%nG    = 4

      CASE (eType_LIN2)
         fs%eType = eType_LIN1
         fs%lShpF = .TRUE.
         fs%eNoN  = 2
         fs%nG    = 2

      CASE DEFAULT
         err = " Cannot choose Taylor-Hood basis"
      END SELECT

      RETURN
      END SUBROUTINE SETTHOODFS
!--------------------------------------------------------------------
      SUBROUTINE GETTHOODFS(fs, lM, lStab, iOpt)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(fsType), INTENT(OUT) :: fs(2)
      TYPE(mshType), INTENT(IN) :: lM
      LOGICAL, INTENT(IN) :: lStab
      INTEGER, INTENT(IN) :: iOpt

      INTEGER i, g

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))
      IF (lStab) THEN
         DO i=1, 2
            fs(i)%nG    = lM%fs(1)%nG
            fs(i)%eType = lM%fs(1)%eType
            fs(i)%lShpF = lM%fs(1)%lShpF
            fs(i)%eNoN  = lM%fs(1)%eNoN
            CALL ALLOCFS(fs(i), nsd)
            fs(i)%w  = lM%fs(1)%w
            fs(i)%xi = lM%fs(1)%xi
            fs(i)%N  = lM%fs(1)%N
            fs(i)%Nx = lM%fs(1)%Nx
         END DO
      ELSE
         IF (iOpt .EQ. 1) THEN
            fs(1)%nG    = lM%fs(1)%nG
            fs(1)%eType = lM%fs(1)%eType
            fs(1)%lShpF = lM%fs(1)%lShpF
            fs(1)%eNoN  = lM%fs(1)%eNoN
            CALL ALLOCFS(fs(1), nsd)
            fs(1)%w  = lM%fs(1)%w
            fs(1)%xi = lM%fs(1)%xi
            fs(1)%N  = lM%fs(1)%N
            fs(1)%Nx = lM%fs(1)%Nx

            fs(2)%nG    = lM%fs(1)%nG
            fs(2)%eType = lM%fs(2)%eType
            fs(2)%lShpF = lM%fs(2)%lShpF
            fs(2)%eNoN  = lM%fs(2)%eNoN
            CALL ALLOCFS(fs(2), nsd)
            fs(2)%w(:)    = lM%fs(1)%w(:)
            fs(2)%xi(:,:) = lM%fs(1)%xi(:,:)
            DO g=1, fs(2)%nG
               CALL GETGNN(nsd, fs(2)%eType, fs(2)%eNoN, fs(2)%xi(:,g),
     2            fs(2)%N(:,g), fs(2)%Nx(:,:,g))
            END DO
         ELSE IF (iOpt .EQ. 2) THEN
            fs(2)%nG    = lM%fs(2)%nG
            fs(2)%eType = lM%fs(2)%eType
            fs(2)%lShpF = lM%fs(2)%lShpF
            fs(2)%eNoN  = lM%fs(2)%eNoN
            CALL ALLOCFS(fs(2), nsd)
            fs(2)%w  = lM%fs(2)%w
            fs(2)%xi = lM%fs(2)%xi
            fs(2)%N  = lM%fs(2)%N
            fs(2)%Nx = lM%fs(2)%Nx

            fs(1)%nG    = lM%fs(2)%nG
            fs(1)%eType = lM%fs(1)%eType
            fs(1)%lShpF = lM%fs(1)%lShpF
            fs(1)%eNoN  = lM%fs(1)%eNoN
            CALL ALLOCFS(fs(1), nsd)
            fs(1)%w(:)    = lM%fs(2)%w(:)
            fs(1)%xi(:,:) = lM%fs(2)%xi(:,:)
            DO g=1, fs(1)%nG
               CALL GETGNN(nsd, fs(1)%eType, fs(1)%eNoN, fs(1)%xi(:,g),
     2            fs(1)%N(:,g), fs(1)%Nx(:,:,g))
            END DO
         END IF
      END IF

      RETURN
      END SUBROUTINE GETTHOODFS
!####################################################################
      SUBROUTINE THOOD_ValRC()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: THflag
      INTEGER(KIND=IKIND) :: a, i, c, s, e, Ac, iM

      INTEGER, ALLOCATABLE :: eNds(:)

      IF ((eq(cEq)%phys .NE. phys_stokes) .AND.
     2    (eq(cEq)%phys .NE. phys_ustruct)) RETURN

      THflag = .FALSE.
      DO iM=1, nMsh
         IF (msh(iM)%nFs .EQ. 2) THEN
            THflag = .TRUE.
            EXIT
         END IF
      END DO

      IF (THflag) THEN
         ALLOCATE(eNds(tnNo))
         eNds(:) = 0
         DO iM=1, nMsh
            IF (msh(iM)%nFs .EQ. 1) CYCLE
            i = msh(iM)%fs(2)%eNoN
            DO e=1, msh(iM)%nEl
               DO a=i+1, msh(iM)%fs(1)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  eNds(Ac) = 1
               END DO
            END DO
         END DO

         DO a=1, tnNo
            IF (eNds(a) .EQ. 1) THEN
               R(nsd+1,a) = 0._RKIND
               s = (nsd+1)*(nsd+1)
               DO i=rowPtr(a), rowPtr(a+1)-1
                  c = colPtr(i)
                  IF (c .EQ. a) THEN
                     Val(s,i) = 1._RKIND
                  ELSE
                     Val(s,i) = 0._RKIND
                  END IF
               END DO
            END IF
         END DO
         DEALLOCATE(eNds)
      END IF

      RETURN
      END SUBROUTINE THOOD_ValRC
!####################################################################


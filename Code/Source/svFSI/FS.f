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

      INTEGER(KIND=IKIND) g, insd

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      IF (lM%lFib) insd = 1

!     The first set of basis is inherited directly from mesh basis
      lM%fs(1)%lShpF = lM%lShpF
      lM%fs(1)%eType = lM%eType
      lM%fs(1)%eNoN  = lM%eNoN
      lM%fs(1)%nG    = lM%nG

      CALL ALLOCFS(lM%fs(1), insd)

      IF (lM%eType .NE. eType_NRB) THEN
         lM%fs(1)%w  = lM%w
         lM%fs(1)%xi = lM%xi
         lM%fs(1)%N  = lM%N
         lM%fs(1)%Nx = lM%Nx
      END IF

!     Sets Taylor-Hood basis if invoked by user (fluid, ustruct, FSI)
      IF (lM%nFs .EQ. 2) THEN
!        Select Taylor-Hood element
         CALL SETTHOODFS(lM%fs(2), lM%fs(1)%eType)

!        Allocate arrays
         CALL ALLOCFS(lM%fs(2), insd)

!        Get Gauss points and shape functions
         CALL GETGIP(insd, lM%fs(2)%eType, lM%fs(2)%nG, lM%fs(2)%w,
     2      lM%fs(2)%xi)
         DO g=1, lM%fs(2)%nG
            CALL GETGNN(insd, lM%fs(2)%eType, lM%fs(2)%eNoN,
     2         lM%fs(2)%xi(:,g), lM%fs(2)%N(:,g), lM%fs(2)%Nx(:,:,g))
         END DO
      END IF

      RETURN
      END SUBROUTINE INITFSMSH
!--------------------------------------------------------------------
!     Initialize function spaces over a face
      SUBROUTINE INITFSFACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) g, iM, insd

      iM   = lFa%iM
      insd = nsd-1
      IF (msh(iM)%lShl) insd = insd - 1
      IF (msh(iM)%lFib) insd = 0

!     The first set of basis is inherited directly from face basis
      lFa%fs(1)%lShpF = msh(iM)%lShpF
      lFa%fs(1)%eType = lFa%eType
      lFa%fs(1)%eNoN  = lFa%eNoN
      lFa%fs(1)%nG    = lFa%nG

      CALL ALLOCFS(lFa%fs(1), insd)

      IF (lFa%eType .NE. eType_NRB) THEN
         lFa%fs(1)%w  = lFa%w
         lFa%fs(1)%xi = lFa%xi
         lFa%fs(1)%N  = lFa%N
         lFa%fs(1)%Nx = lFa%Nx
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
      END IF

      RETURN
      END SUBROUTINE INITFSFACE
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

      ALLOCATE(fs%w(nG), fs%xi(insd,nG), fs%N(eNoN,nG),
     2   fs%Nx(insd,eNoN,nG))

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
!--------------------------------------------------------------------
!     Sets Tayloor-Hood basis for a parent element type
      SUBROUTINE SETTHOODFS(fs, eType)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fsType), INTENT(INOUT) :: fs
      INTEGER(KIND=IKIND), INTENT(IN) :: eType

      SELECT CASE (eType)
      CASE (eType_QUD)
         fs%eType = eType_LIN
         fs%lShpF = .TRUE.
         fs%eNoN  = 2
         fs%nG    = 2

      CASE (eType_QTR)
         fs%eType = eType_TRI
         fs%lShpF = .TRUE.
         fs%eNoN  = 3
         fs%nG    = 3

      CASE (eType_BIQ)
         fs%eType = eType_BIL
         fs%lShpF = .FALSE.
         fs%eNoN  = 4
         fs%nG    = 4

      CASE (eType_QTE)
         fs%eType = eType_TET
         fs%lShpF = .TRUE.
         fs%eNoN  = 4
         fs%nG    = 4

      CASE DEFAULT
         err = " Cannot choose Taylor-Hood basis"
      END SELECT

      RETURN
      END SUBROUTINE SETTHOODFS
!####################################################################


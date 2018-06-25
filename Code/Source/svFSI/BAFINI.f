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
!     This routine initializes required structure for boundaries,
!     faces, those that interface with FSILS and cplBC.
!
!--------------------------------------------------------------------

      SUBROUTINE BAFINI()

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER iM, iFa, iEq, i, a, iBc, lsPtr

      INTEGER, ALLOCATABLE :: gNodes(:)

      DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            msh(iM)%fa(iFa)%iM = iM
            CALL FACEINI(msh(iM)%fa(iFa))
         END DO
      END DO

      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            CALL BCINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa))
         END DO
      END DO

!     cplBC faces are initialized here
      iEq = 1
      IF (ALLOCATED(cplBC%fa)) DEALLOCATE(cplBC%fa)
      IF (ALLOCATED(cplBC%xn)) DEALLOCATE(cplBC%xn)
      ALLOCATE(cplBC%fa(cplBC%nFa), cplBC%xn(cplBC%nX))
      IF (cplBC%coupled) THEN
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_cpl)) THEN
               i = eq(iEq)%bc(iBc)%cplBCPtr
               cplBC%fa(i)%name = TRIM(msh(iM)%fa(iFa)%name)
               cplBC%fa(i)%y    = 0D0
               IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) THEN
                  cplBC%fa(i)%bGrp = cplBC_Dir
               ELSE IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
                  cplBC%fa(i)%bGrp = cplBC_Neu
                  IF (cplBC%schm .NE. cplBC_E) eq(iEq)%bc(iBc)%bType=
     2               IBSET(eq(iEq)%bc(iBc)%bType,bType_res)
               ELSE
                  err = "Not a compatible cplBC_type"
               END IF
            END IF
         END DO
         IF (cplBC%schm .NE. cplBC_E) CALL CALCDERCPLBC
      END IF

!     Setting up FSILS
      lsPtr = 0
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            CALL FSILSINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), lsPtr)
         END DO
      END DO

      IF (mvMsh) THEN
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct)) i = i + 1
         END DO
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct)) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         CALL FSILS_BC_CREATE(lhs, nFacesLS, i, nsd, BC_TYPE_Dir,
     &                        gNodes)
      END IF

      RETURN
      END SUBROUTINE BAFINI

!--------------------------------------------------------------------
!     Initializing faces
      SUBROUTINE FACEINI(lFa)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(faceType), INTENT(INOUT) :: lFa

      LOGICAL flag
      INTEGER e, a, Ac, iM, g
      REAL(KIND=8) tmp, area, n(nsd)

      REAL(KIND=8), ALLOCATABLE :: s(:), sV(:,:)

!     Calculating the center of the face, diameter and its area
      iM = lFa%iM
      ALLOCATE(s(tnNo), sV(nsd,tnNo))
      IF (ALLOCATED(lFa%nV)) DEALLOCATE(lFa%nV)
      ALLOCATE(lFa%nV(nsd,lFa%nNo))
      s        = 1D0
      sV       = 0D0
      area     = Integ(lFa,s)
      lFa%area = area
!     Making sure area is not zero, since it will cause issues later on
      IF (ISZERO(area)) THEN
         IF (cm%mas()) wrn = "<"//TRIM(lFa%name)//"> area is zero"
      END IF
      DO e=1, lFa%nEl
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)
         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, n)
            DO a=1, lFa%eNoN
               Ac       = lFa%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + n*lFa%N(a,g)*lFa%w(g)
            END DO
         END DO
      END DO

      flag = .TRUE.
      CALL COMMU(sV)
      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         tmp = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(tmp)) THEN
            IF (flag) THEN
               wrn = "Skipping normal calculation of node "//a//
     2            " in face <"//TRIM(lFa%name)//">"
               flag = .FALSE.
            END IF
            lFa%nV(:,a) = 0D0
            lFa%nV(1,a) = 1D0
            CYCLE
         END IF
         lFa%nV(:,a) = sV(:,Ac)/tmp
      END DO

      RETURN
      END SUBROUTINE FACEINI
!--------------------------------------------------------------------
      SUBROUTINE BCINI(lBc, lFa)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER iM, iFa, jFa, i, a, b, Ac, j, ierr
      REAL(KIND=8) tmp, nV(nsd), center(nsd), maxN

      INTEGER, ALLOCATABLE :: gNodes(:), sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = "Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      iM  = lFa%iM
      iFa = lBc%iFa
      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lFa%nNo))

      ALLOCATE(s(tnNo), sCount(cm%np()), disp(cm%np()))
      s = 0D0
      IF (BTEST(lBc%bType,bType_flat)) THEN
!     Just a constant value for Flat profile
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = 1D0
         END DO
      ELSE IF (BTEST(lBc%bType,bType_para)) THEN
!     Here is the method that is used for imposing parabolic profile:
!     1- Find the coordinate of the points on the boundary 2- find unit
!     vector from center to each of points on the boundary: ew
!     3- maximize ew(i).e where e is the unit vector from current
!     point to the center 4- Use the point i as the diam here
         DO i=1, nsd
            center(i) = Integ(lFa, x, i)/lFa%area
         END DO
         ALLOCATE(gNodes(tnNo), sVl(nsd,lFa%nNo), sV(nsd,tnNo))
!     gNodes is one if a node located on the boundary (beside iFa)
         gNodes = 0
         DO jFa=1, msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, msh(iM)%fa(jFa)%nNo
               Ac         = msh(iM)%fa(jFa)%gN(a)
               gNodes(Ac) = 1
            END DO
         END DO
!     "j" is a counter for the number of nodes that are located on the
!     boundary of lFa and sVl contains the list of their coordinates
         j   = 0
         sVl = 0D0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (gNodes(Ac) .EQ. 1) THEN
               j = j + 1
               sVl(:,j) = x(:,Ac)
            END IF
         END DO
         IF (.NOT.cm%seq()) THEN
!     Getting the length data that is going to be received at each proc
            i = j*nsd
            CALL MPI_ALLGATHER(i, 1, mpint, sCount, 1, mpint, cm%com(),
     2         ierr)
            disp(1) = 0
            DO i=2, cm%np()
               disp(i) = disp(i-1) + sCount(i-1)
            END DO
            CALL MPI_ALLGATHERV(sVl, j*nsd, mpreal, sV, sCount, disp,
     2         mpreal, cm%com(), ierr)
            j = SUM(sCount)/nsd
         ELSE
            sV(:,1:j) = sVl(:,1:j)
         END IF
         IF (cm%mas() .AND. j.EQ.0) err = "No perimeter"//
     2      " found for face "//lFa%name
!     sVl will keep the normal unit vector from center to perimeter
         DEALLOCATE(sVl)
         ALLOCATE(sVl(nsd,j))
         DO a=1, j
            sV(:,a)  = sV(:,a) - center
            sVl(:,a) = sV(:,a)/SQRT(NORM(sV(:,a)))
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            nV = x(:,Ac) - center
            maxN = NORM(nV, sVl(:,1))
            i = 1
            DO b=2, j
               tmp = NORM(nV, sVl(:,b))
               IF (tmp .GT. maxN) THEN
                  maxN = tmp
                  i = b
               END IF
            END DO
            s(Ac) = 1D0 - NORM(nV)/NORM(sV(:,i))
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ud)) THEN
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = lBc%gx(a)
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ddep)) THEN
         s = 1D0
      END IF

!     Now correcting the inlet BC for the inlet ring
      IF (BTEST(lBc%bType,bType_zp)) THEN
         DO jFa=1, msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, msh(iM)%fa(jFa)%nNo
               Ac    = msh(iM)%fa(jFa)%gN(a)
               s(Ac) = 0D0
            END DO
         END DO
      END IF

!     Normalizing the profile for flux
      tmp = 1D0
      IF (BTEST(lBc%bType,bType_flx)) THEN
         tmp = Integ(lFa, s)
         IF (ISZERO(tmp)) THEN
            tmp = 1D0
            wrn = "Using face <"//TRIM(lFa%name)//
     2         "> to impose BC led to no non-zero node."
         END IF
      END IF

      DO a=1, lFa%nNo
         Ac        = lFa%gN(a)
         lBc%gx(a) = s(Ac)/tmp
      END DO

      RETURN
      END SUBROUTINE BCINI
!--------------------------------------------------------------------
      SUBROUTINE FSILSINI(lBc, lFa, lsPtr)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lsPtr
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER a, e, Ac, g, iM, i
      REAL(KIND=8) n(nsd)
      LOGICAL :: eDrn

      INTEGER, ALLOCATABLE :: gNodes(:)
      REAL(KIND=8), ALLOCATABLE :: sV(:,:), sVl(:,:)

      iM  = lFa%iM
      ALLOCATE(sVl(nsd,lFa%nNo), sV(nsd,tnNo), gNodes(lFa%nNo))
      DO a=1, lFa%nNo
         gNodes(a) = lFa%gN(a)
      END DO

      IF (BTEST(lBc%bType,bType_Dir)) THEN
         lsPtr     = lsPtr + 1
         lBc%lsPtr = lsPtr
         sVl = 0D0
         eDrn = .FALSE.
         DO i=1, nsd
            IF (lBc%eDrn(i) .NE. 0) THEN
               eDrn = .TRUE.
               EXIT
            END IF
         END DO
         IF (eDrn) THEN
            sVl = 1D0
            DO i=1, nsd
               IF (lBc%eDrn(i) .NE. 0) sVl(i,:) = 0D0
            END DO
         END IF
         CALL FSILS_BC_CREATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Dir,
     2      gNodes, sVl)
      ELSE IF (BTEST(lBc%bType,bType_Neu)) THEN
         IF (BTEST(lBc%bType,bType_res)) THEN
            sV = 0D0
            DO e=1, lFa%nEl
               IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM),lFa,e)
               DO g=1, lFa%nG
                  CALL GNNB(lFa, e, g, n)
                  DO a=1, lFa%eNoN
                     Ac = lFa%IEN(a,e)
                     sV(:,Ac) = sV(:,Ac) + lFa%N(a,g)*lFa%w(g)*n
                  END DO
               END DO
            END DO
            DO a=1, lFa%nNo
               Ac       = lFa%gN(a)
               sVl(:,a) = sV(:,Ac)
            END DO
            lsPtr     = lsPtr + 1
            lBc%lsPtr = lsPtr
            CALL FSILS_BC_CREATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Neu,
     2         gNodes, sVl)
         ELSE
            lBc%lsPtr = 0
         END IF
      ELSE
         err = "Unxpected bType in FSILSINI"
      END IF

      RETURN
      END SUBROUTINE FSILSINI

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

!     Compute face normals and area
      DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            msh(iM)%fa(iFa)%iM = iM
            IF (msh(iM)%lFib) CYCLE
            CALL FACEINI(msh(iM)%fa(iFa))
         END DO
         IF (msh(iM)%lShl) CALL SHLINI(msh(iM))
      END DO

!     Initialize face BC profile
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            CALL BCINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa))
            IF (msh(iM)%lShl) THEN
               CALL SHLBCINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), msh(iM))
            END IF
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
            IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_cpl) .OR.
     2          BTEST(eq(iEq)%bc(iBc)%bType,bType_RCR)) THEN
               i = eq(iEq)%bc(iBc)%cplBCPtr
               cplBC%fa(i)%name = TRIM(msh(iM)%fa(iFa)%name)
               cplBC%fa(i)%y    = 0D0
               IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) THEN
                  cplBC%fa(i)%bGrp = cplBC_Dir
               ELSE IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
                  cplBC%fa(i)%bGrp = cplBC_Neu
                  IF (cplBC%schm .NE. cplBC_E) eq(iEq)%bc(iBc)%bType=
     2               IBSET(eq(iEq)%bc(iBc)%bType,bType_res)

!                 Copy RCR structure from bc() to cplBC()
                  cplBC%fa(i)%RCR%Rp = eq(iEq)%bc(iBc)%RCR%Rp
                  cplBC%fa(i)%RCR%C  = eq(iEq)%bc(iBc)%RCR%C
                  cplBC%fa(i)%RCR%Rd = eq(iEq)%bc(iBc)%RCR%Rd
                  cplBC%fa(i)%RCR%Pd = eq(iEq)%bc(iBc)%RCR%Pd
               ELSE
                  err = "Not a compatible cplBC_type"
               END IF
            END IF
         END DO
         IF (cplBC%useGenBC) CALL genBC_Integ_X('I')
         IF (cplBC%schm .NE. cplBC_E) CALL CALCDERCPLBC()
      END IF

!     Setting up FSILS
      lsPtr = 0
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            eq(iEq)%bc(iBc)%lsPtr = 0
            CALL FSILSINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), lsPtr)
         END DO
      END DO

      IF (ibLSptr .NE. 0) THEN
         SELECT CASE (ib%mthd)
         CASE (ibMthd_SSM)
            i = SUM(iblank(:))
            ALLOCATE(gNodes(i))
            i = 0
            DO a=1, tnNo
               IF (iblank(a) .EQ. 1) THEN
                  i = i + 1
                  gNodes(i) = a
               END IF
            END DO

         CASE DEFAULT
            err = " Invalid IB method (BAFINI)"
         END SELECT

         lsPtr = ibLSptr
         CALL FSILS_BC_CREATE(lhs, lsPtr, i, nsd, BC_TYPE_Dir, gNodes)
         DEALLOCATE(gNodes)
      END IF

      IF (mvMsh) THEN
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct) .OR.
     2          ISDOMAIN(1, a, phys_vms_struct) .OR.
     3          ISDOMAIN(1, a, phys_lElas)) i = i + 1
         END DO
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct) .OR.
     2          ISDOMAIN(1, a, phys_vms_struct) .OR.
     3          ISDOMAIN(1, a, phys_lElas)) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         lsPtr = nFacesLS
         CALL FSILS_BC_CREATE(lhs, lsPtr, i, nsd, BC_TYPE_Dir, gNodes)
         DEALLOCATE(gNodes)
      END IF

      IF (cmmInit) THEN
         i = SUM(cmmBdry)
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (cmmBdry(a) .EQ. 1) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         lsPtr = nFacesLS
         CALL FSILS_BC_CREATE(lhs, lsPtr, i, nsd, BC_TYPE_Dir, gNodes)
         DEALLOCATE(gNodes)
      END IF

      RETURN
      END SUBROUTINE BAFINI
!####################################################################
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
      std = " Area of face <"//TRIM(lFa%name)//"> is "//STR(area)
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

      IF (BTEST(lBc%bType,bType_Robin) .AND. .NOT.dFlag) err =
     2    "Robin BC can be set for a displacement-based eqn only"

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
!####################################################################
      SUBROUTINE FSILSINI(lBc, lFa, lsPtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lsPtr
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER a, e, Ac, g, iM, i, nNo
      REAL(KIND=8) n(nsd)
      LOGICAL :: eDrn

      INTEGER, ALLOCATABLE :: gNodes(:)
      REAL(KIND=8), ALLOCATABLE :: sV(:,:), sVl(:,:)

      iM  = lFa%iM
      nNo = lFa%nNo
      ALLOCATE(sVl(nsd,nNo), sV(nsd,tnNo), gNodes(nNo))
      DO a=1, nNo
         gNodes(a) = lFa%gN(a)
      END DO

      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (lBc%weakDir) THEN
            lBc%lsPtr = 0
         ELSE
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
     2         gNodes, sVl)
         END IF
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
      ELSE IF (BTEST(lBc%bType,bType_trac)) THEN
         lBc%lsPtr = 0
      ELSE IF (BTEST(lBc%bType,bType_CMM)) THEN
         lsPtr     = lsPtr + 1
         lBc%lsPtr = lsPtr

         nNo = 0
         DO a=1, lFa%nNo
            IF (ISZERO(lBc%gx(a))) nNo = nNo + 1
         END DO
         DEALLOCATE(gNodes, sVl)
         ALLOCATE(sVl(nsd,nNo), gNodes(nNo))
         sVl  = 0D0

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

         nNo = 0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (ISZERO(lBc%gx(a))) THEN
               nNo = nNo + 1
               gNodes(nNo) = Ac
            END IF
         END DO

         CALL FSILS_BC_CREATE(lhs, lsPtr, nNo, nsd, BC_TYPE_Dir, gNodes,
     2      sVl)
      ELSE
         err = "Unxpected bType in FSILSINI"
      END IF

      DEALLOCATE(sVl, sV, gNodes)

      RETURN
      END SUBROUTINE FSILSINI
!####################################################################
!     Compute shell extended IEN for triangular elements. Reqd. to
!     resolve bending moments
      SUBROUTINE SETSHLXIEN(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER :: a, b, e, f, Ac, Bc, nEl, eNoN, ep(2,3)

      INTEGER, ALLOCATABLE :: incN(:)

      eNoN = lM%eNoN
      nEl  = lM%nEl

      ep   = RESHAPE((/2,3,3,1,1,2/), SHAPE(ep))
      ALLOCATE(incN(eNoN), lM%eIEN(eNoN,nEl), lM%sbc(eNoN,nEl))

      lM%eIEN = 0
      lM%sbc  = 0
      DO e=1, nEl
         DO a=1, eNoN
            Ac = lM%IEN(ep(1,a),e)
            Bc = lM%IEN(ep(2,a),e)
            DO f=1, nEl
               IF (e .EQ. f) CYCLE
               incN = 0
               DO b=1, eNoN
                  IF (lM%IEN(b,f).EQ.Ac .OR. lM%IEN(b,f).EQ.Bc)
     2               incN(b) = incN(b) + 1
               END DO
               IF (SUM(incN) .EQ. 2) THEN
                  DO b=1, eNoN
                     IF (incN(b) .EQ. 0) THEN
                        lM%eIEN(a,e) = lM%IEN(b,f)
                        EXIT
                     END IF
                  END DO
                  EXIT
               END IF
            END DO
            IF (lM%eIEN(a,e) .EQ. 0) THEN
               lM%sbc(a,e) = IBSET(lM%sbc(a,e), bType_free)
            END IF
         END DO
      END DO
      DEALLOCATE(incN)

      RETURN
      END SUBROUTINE SETSHLXIEN
!--------------------------------------------------------------------
!     Initializing shell normals and area
      SUBROUTINE SHLINI(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER :: a, e, g, Ac, nNo, nEl, eNoN
      LOGICAL :: flag
      REAL(KIND=8) :: Jac, area, nV(nsd), tmpR(nsd,nsd-1)

      REAL(KIND=8), ALLOCATABLE :: xl(:,:), sV(:,:)

      nNo  = lM%nNo
      nEl  = lM%nEl
      eNoN = lM%eNoN

!     Compute shell director (normal)
      ALLOCATE(xl(nsd,eNoN), sV(nsd,tnNo), lM%nV(nsd,nNo))
      sV   = 0D0
      area = 0D0
      DO e=1, nEl
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            xl(:,a) = x(:,Ac)
         END DO

         DO g=1, lM%nG
            CALL GNNS(eNoN, lM%Nx(:,:,g), xl, nV, tmpR, tmpR)
            Jac = SQRT(NORM(nV))
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + lM%w(g)*lM%N(a,g)*nV
               area     = area     + lM%w(g)*lM%N(a,g)*Jac
            END DO
         END DO
      END DO
      area = cm%reduce(area)
      CALL COMMU(sV)

      std = " Area of the shell surface <"//TRIM(lM%name)//"> is "//
     2   STR(area)

      flag = .TRUE.
      DO a=1, nNo
         Ac  = lM%gN(a)
         Jac = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(Jac)) THEN
            IF (flag) THEN
               wrn = "Skipping normal calculation of node "//a//
     2            " in mesh <"//TRIM(lM%name)//">"
               flag = .FALSE.
            END IF
            lM%nV(:,a) = 0D0
            lM%nV(1,a) = 1D0
            CYCLE
         END IF
         lM%nV(:,a) = sV(:,Ac)/Jac
      END DO

      DEALLOCATE(xl, sV)

      RETURN
      END SUBROUTINE SHLINI
!--------------------------------------------------------------------
!     Initializing shell boundary condition variables
      SUBROUTINE SHLBCINI(lBc, lFa, lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER :: a, b, e, Ac, Bc, Ec
      LOGICAL :: bFlag

      IF (lFa%eType .EQ. eType_NRB) RETURN
      DO e=1, lFa%nEl
         Ec = lFa%gE(e)
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,Ec)
            bflag = .FALSE.
            DO b=1, lFa%eNoN
               Bc = lFa%IEN(b,e)
               IF (Ac .EQ. Bc) THEN
                  bFlag = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT.bFlag) THEN
               IF (.NOT.BTEST(lM%sbc(a,Ec),bType_free)) err =
     2            "BC detected on a non-boundary shell element. "//
     3            "Correction needed"
               lM%sbc(a,Ec) = IBCLR(lM%sbc(a,Ec), bType_free)
               IF (BTEST(lBc%bType,bType_free)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_free)
               ELSE IF (BTEST(lBc%bType,bType_fix)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_fix)
               ELSE IF (BTEST(lBc%bType,bType_hing)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_hing)
               ELSE IF (BTEST(lBc%bType,bType_symm)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_symm)
               END IF
               EXIT
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE SHLBCINI
!####################################################################

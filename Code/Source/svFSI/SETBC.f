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
!     In this subroutine Dirichlet and Nuemann BC are updated based
!     on the specified flow rate or the type of outlet BC.
!
!--------------------------------------------------------------------

      SUBROUTINE SETBCDIR(lA, lY, lD)

      USE COMMOD

      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: lA(tDof, tnNo), lY(tDof, tnNo),
     2   lD(tDof, tnNo)

      INTEGER iFa, a, Ac, iEq, iBc, s, e, iM, nNo, lDof, i
      LOGICAL :: eDir(maxnsd)

      REAL(KIND=8), ALLOCATABLE :: tmpA(:,:), tmpY(:,:)

      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            IF (.NOT.BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) CYCLE
            s = eq(iEq)%s
            e = eq(iEq)%e
            IF (eq(iEq)%dof .EQ. nsd+1) e = e - 1
            eDir = .FALSE.
            lDof = 0
            DO i=1, nsd
               IF (eq(iEq)%bc(iBc)%eDrn(i) .NE. 0) THEN
                  eDir(i) = .TRUE.
                  lDof = lDof + 1
               END IF
            END DO
            IF (lDof .EQ. 0) lDof = e - s + 1
            iFa  = eq(iEq)%bc(iBc)%iFa
            iM   = eq(iEq)%bc(iBc)%iM
            nNo  = msh(iM)%fa(iFa)%nNo
            IF (ALLOCATED(tmpA)) DEALLOCATE(tmpA, tmpY)
            ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
            CALL SETBCDIRL(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), tmpA, tmpY,
     2         lDof)
            IF (ANY(eDir)) THEN
               DO a=1, msh(iM)%fa(iFa)%nNo
                  Ac = msh(iM)%fa(iFa)%gN(a)
                  IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_impD)) THEN
                     DO i=1, nsd
                        lDof = 0
                        IF (eDir(i)) THEN
                           lDof = lDof + 1
                           lY(i,Ac) = tmpA(lDof,a)
                           lD(i,Ac) = tmpY(lDof,a)
                        END IF
                     END DO
                  ELSE
                     DO i=1, nsd
                        lDof = 0
                        IF (eDir(i)) THEN
                           lDof = lDof + 1
                           lA(i,Ac) = tmpA(lDof,a)
                           lY(i,Ac) = tmpY(lDof,a)
                        END IF
                     END DO
                  END IF
               END DO
            ELSE
               DO a=1, msh(iM)%fa(iFa)%nNo
                  Ac = msh(iM)%fa(iFa)%gN(a)
                  IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_impD)) THEN
                     lY(s:e,Ac) = tmpA(:,a)
                     lD(s:e,Ac) = tmpY(:,a)
                  ELSE
                     lA(s:e,Ac) = tmpA(:,a)
                     lY(s:e,Ac) = tmpY(:,a)
                  END IF
               END DO
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE SETBCDIR
!--------------------------------------------------------------------
      SUBROUTINE SETBCDIRL(lBc, lFa, lA, lY, lDof)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: lDof
      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(INOUT) :: lA(lDof,lFa%nNo), lY(lDof,lFa%nNo)

      INTEGER :: a, i
      REAL(KIND=8) :: dirY, dirA, nV(nsd)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (lDof .NE. lBc%gm%dof) err = "Inconsistent DOF"
         CALL IGBC(lBc%gm, lY, lA)
         RETURN
      ELSE IF (BTEST(lBc%bType,bType_ustd)) THEN
         CALL IFFT(lBc%gt, dirY, dirA)
      ELSE ! std / cpl
         dirA = 0D0
         dirY = lBc%g
      END IF
      IF (lDof .EQ. nsd) THEN
         DO a=1, lFa%nNo
            nV      = lFa%nV(:,a)
            lA(:,a) = dirA*lBc%gx(a)*nV
            lY(:,a) = dirY*lBc%gx(a)*nV
         END DO
      ELSE
         DO a=1, lFa%nNo
            DO i=1, lDof
               lA(i,a) = dirA*lBc%gx(a)
               lY(i,a) = dirY*lBc%gx(a)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE SETBCDIRL

!####################################################################
!     Here for the outlets
      SUBROUTINE SETBCNEU(Yg, Dg)

      USE COMMOD

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER iFa, iBc, iM

      DO iBc=1, eq(cEq)%nBc
         IF (.NOT.BTEST(eq(cEq)%bc(iBc)%bType,bType_Neu)) CYCLE
         iFa = eq(cEq)%bc(iBc)%iFa
         iM  = eq(cEq)%bc(iBc)%iM
         CALL SETBCNEUL(eq(cEq)%bc(iBc), msh(iM)%fa(iFa), Yg, Dg)
      END DO

      RETURN
      END SUBROUTINE SETBCNEU
!--------------------------------------------------------------------
      SUBROUTINE SETBCNEUL(lBc, lFa, Yg, Dg)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER a, Ac, e, s
      REAL(KIND=8) Q, h, tmp

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: hg(:), tmpA(:), yl(:,:), hl(:),
     2   dl(:,:)
!     Geting the contribution of Neu BC
      IF (BTEST(lBc%bType,bType_cpl)) THEN
         h = lBc%g
      ELSE
         IF (BTEST(lBc%bType,bType_gen)) THEN
!     Using "hg" as a temporary variable here
            ALLOCATE(tmpA(lFa%nNo), hg(lFa%nNo))
            CALL IGBC(lBc%gm, tmpA, hg)
            DEALLOCATE(hg)
         ELSE IF (BTEST(lBc%bType,bType_res)) THEN
            Q = Integ(lFa, Yn, eq(cEq)%s, eq(cEq)%s+nsd-1)
            h = Q*lBc%r
         ELSE IF (BTEST(lBc%bType,bType_std)) THEN
            h = lBc%g
         ELSE IF (BTEST(lBc%bType,bType_ustd)) THEN
            CALL IFFT(lBc%gt, h, tmp)
         ELSE
            err = "Correction in SETBCNEU is needed"
         END IF
      END IF

      ALLOCATE(hg(tnNo), yl(tDof,lFa%eNoN), hl(lFa%eNoN), ptr(lFa%eNoN),
     2   dl(tDof,msh(lFa%iM)%eNoN))
!     Transforming it to a unified format
      IF (BTEST(lBc%bType,bType_gen)) THEN
         DO a=1, lFa%nNo
            Ac     = lFa%gN(a)
            hg(Ac) = tmpA(a)
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ddep)) THEN
         s = eq(cEq)%s
         IF (eq(cEq)%dof .EQ. 1) THEN
            DO a=1, lFa%nNo
               Ac     = lFa%gN(a)
               hg(Ac) = -h*Dg(s,Ac)
            END DO
         ELSE IF (eq(cEq)%dof .GE. nsd) THEN
            DO a=1, lFa%nNo
               Ac     = lFa%gN(a)
               hg(Ac) = -h*NORM(Dg(s:s+nsd-1,Ac),lFa%nV(:,a))
            END DO
         ELSE
            err = "Correction in SETBCNEU is needed"
         END IF
      ELSE
         DO a=1, lFa%nNo
            Ac     = lFa%gN(a)
            hg(Ac) = -h*lBc%gx(a)
         END DO
      END IF

!     Constructing LHS/RHS contribution and assembiling them
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac      = lFa%IEN(a,e)
               ptr(a)  = Ac
            yl(:,a) = Yg(:,Ac)
            hl(a)   = hg(Ac)
         END DO
         DO a=1, msh(lFa%iM)%eNoN
            Ac      = msh(lFa%iM)%IEN(a,lFa%gE(e))
            dl(:,a) = Dg(:,Ac)
         END DO
!     Add Neumann BCs contribution to the LHS/RHS
         CALL BCONSTRUCT(lFa, yl, dl, hl, ptr, e)
      END DO

      RETURN
      END SUBROUTINE SETBCNEUL

!####################################################################
!     cplBC is set here
      SUBROUTINE SETBCCPL

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, PARAMETER :: iEq = 1
      INTEGER iFa, ptr, iBc, iM
      REAL(KIND=8) tmp

      IF (cplBC%schm .EQ. cplBC_I) THEN
         CALL CALCDERCPLBC
      ELSE
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            ptr = eq(iEq)%bc(iBc)%cplBCptr
            IF (ptr .NE. 0) THEN
               IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
                  cplBC%fa(ptr)%Qo = Integ(msh(iM)%fa(iFa),Yo,1,nsd)
                  cplBC%fa(ptr)%Qn = Integ(msh(iM)%fa(iFa),Yn,1,nsd)
               ELSE IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) THEN
                  tmp = msh(iM)%fa(iFa)%area
                  cplBC%fa(ptr)%Po = Integ(msh(iM)%fa(iFa),Yo,nsd+1)/tmp
                  cplBC%fa(ptr)%Pn = Integ(msh(iM)%fa(iFa),Yn,nsd+1)/tmp
               END IF
            END IF
         END DO
         CALL cplBC_Integ_X
      END IF

      DO iBc=1, eq(iEq)%nBc
         iFa = eq(iEq)%bc(iBc)%iFa
         ptr = eq(iEq)%bc(iBc)%cplBCptr
         IF (ptr .NE. 0) eq(iEq)%bc(iBc)%g = cplBC%fa(ptr)%y
      END DO

      RETURN
      END SUBROUTINE SETBCCPL

!--------------------------------------------------------------------
!     cplBC derivative is calculated here
      SUBROUTINE CALCDERCPLBC

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, PARAMETER :: iEq = 1
      REAL(KIND=8), PARAMETER :: absTol = 1D-8, relTol = 1D-5

      INTEGER iFa, i, j, ptr, iBc, iM
      REAL(KIND=8) orgQ, orgY, diff, area

      IF (ALL(cplBC%fa%bGrp.EQ.cplBC_Dir)) RETURN

      DO iBc=1, eq(iEq)%nBc
         iFa = eq(iEq)%bc(iBc)%iFa
         iM  = eq(iEq)%bc(iBc)%iM
         ptr = eq(iEq)%bc(iBc)%cplBCptr
         IF (ptr .NE. 0) THEN
            IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
               cplBC%fa(ptr)%Qo = Integ(msh(iM)%fa(iFa),Yo,1,nsd)
               cplBC%fa(ptr)%Qn = Integ(msh(iM)%fa(iFa),Yn,1,nsd)
            ELSE IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) THEN
               area = msh(iM)%fa(iFa)%area
               cplBC%fa(ptr)%Po = Integ(msh(iM)%fa(iFa),Yo,nsd+1)/area
               cplBC%fa(ptr)%Pn = Integ(msh(iM)%fa(iFa),Yn,nsd+1)/area
            END IF
         END IF
      END DO
      CALL cplBC_Integ_X

      j    = 0
      diff = 0D0
      DO iBc=1, eq(iEq)%nBc
         i = eq(iEq)%bc(iBc)%cplBCptr
         IF (i.NE.0 .AND. BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
            diff = diff + cplBC%fa(i)%Qo**2D0
            j = j + 1
         END IF
      END DO
      diff = SQRT(diff/REAL(j,8))
      IF (diff*relTol .LT. absTol) THEN
         diff = absTol
      ELSE
         diff = diff*relTol
      END IF

      DO iBc=1, eq(iEq)%nBc
         i = eq(iEq)%bc(iBc)%cplBCptr
         IF (i.NE.0 .AND. BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
            orgY = cplBC%fa(i)%y
            orgQ = cplBC%fa(i)%Qn
            cplBC%fa(i)%Qn = cplBC%fa(i)%Qn + diff

            CALL cplBC_Integ_X

            eq(iEq)%bc(iBc)%r = (cplBC%fa(i)%y - orgY)/diff

            cplBC%fa(i)%y  = orgY
            cplBC%fa(i)%Qn = orgQ
         END IF
      END DO

      RETURN
      END SUBROUTINE CALCDERCPLBC

!--------------------------------------------------------------------
!     Interface to call 0D code
      SUBROUTINE cplBC_Integ_X

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER fid, iFa, ios
      REAL(KIND=8), ALLOCATABLE :: y(:)

      IF (cm%mas()) THEN
         fid = 1
         OPEN(fid, FILE=cplBC%commuName, FORM='UNFORMATTED')
         WRITE(fid) version
         WRITE(fid) cplBC%nFa
         WRITE(fid) cplBC%nX
         WRITE(fid) dt
         WRITE(fid) time-dt
         WRITE(fid) cplBC%xo
         DO iFa=1, cplBC%nFa
            WRITE(fid) cplBC%fa(iFa)%bGrp
            WRITE(fid) cplBC%fa(iFa)%Qo
            WRITE(fid) cplBC%fa(iFa)%Qn
            WRITE(fid) cplBC%fa(iFa)%Po
            WRITE(fid) cplBC%fa(iFa)%Pn
            WRITE(fid) cplBC%fa(iFa)%name
         END DO
         CLOSE(fid)

         CALL SYSTEM(TRIM(cplBC%binPath)//" "//TRIM(cplBC%commuName))

         OPEN(fid,FILE=cplBC%commuName,STATUS='OLD',FORM='UNFORMATTED')
         READ(fid,IOSTAT=ios) cplBC%xn
         DO iFa=1, cplBC%nFa
            IF (ios .GT. 0) EXIT
            READ(fid,IOSTAT=ios) cplBC%fa(iFa)%y
         END DO
         CLOSE(fid)
         IF (ios .GT. 0) err = "Issue with reading "//
     2      TRIM(cplBC%commuName)//" possibly due to an issue with "//
     3      TRIM(cplBC%binPath)
      END IF

      IF (.NOT.cm%seq()) THEN
         ALLOCATE(y(cplBC%nFa))
         IF (cm%mas()) y = cplBC%fa%y
         CALL cm%bcast(cplBC%xn)
         CALL cm%bcast(y)
         IF (cm%slv()) cplBC%fa%y = y
      END IF

      RETURN
      END SUBROUTINE cplBC_Integ_X


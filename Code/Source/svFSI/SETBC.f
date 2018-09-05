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
            IF(BTEST(eq(iEq)%bc(iBc)%bType,bType_CMM)) THEN
               s = eq(iEq)%s
               e = eq(iEq)%e
               IF (eq(iEq)%dof .EQ. nsd+1) e = e - 1
               iFa  = eq(iEq)%bc(iBc)%iFa
               iM   = eq(iEq)%bc(iBc)%iM
               DO a=1, msh(iM)%fa(iFa)%nNo
                  IF (ISZERO(eq(iEq)%bc(iBc)%gx(a))) THEN
                     Ac = msh(iM)%fa(iFa)%gN(a)
                     lA(s:e,Ac) = 0D0
                     lY(s:e,Ac) = 0D0
                  END IF
               END DO
            END IF ! END bType_CMM

            IF (.NOT.BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) CYCLE
            IF (eq(iEq)%bc(iBc)%weakDir) CYCLE
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
                           lY(s+i-1,Ac) = tmpA(lDof,a)
                           lD(s+i-1,Ac) = tmpY(lDof,a)
                        END IF
                     END DO
                  ELSE
                     DO i=1, nsd
                        lDof = 0
                        IF (eDir(i)) THEN
                           lDof = lDof + 1
                           lA(s+i-1,Ac) = tmpA(lDof,a)
                           lY(s+i-1,Ac) = tmpY(lDof,a)
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

         IF (ibFlag) THEN
            CALL IB_SETBCDIR(ib%Ao, ib%Yo, ib%Uo)
            IF (ib%mthd .EQ. ibMthd_SSM) THEN
               IF (eq(iEq)%phys .EQ. phys_fluid .OR.
     2             eq(iEq)%phys .EQ. phys_FSI) THEN
                  CALL IB_PRJCTU(lY, lD, ib%Yo)
               END IF
            END IF
         END IF
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
      REAL(KIND=8), ALLOCATABLE :: hg(:), tmpA(:), yl(:,:), hl(:)

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

      ALLOCATE(hg(tnNo), yl(tDof,lFa%eNoN), hl(lFa%eNoN), ptr(lFa%eNoN))
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
      IF (eq(cEq)%phys .EQ. phys_ustruct) THEN
         CALL BUSTRUCTNEU(lFa, hg, Dg)
      ELSE
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac      = lFa%IEN(a,e)
               ptr(a)  = Ac
               yl(:,a) = Yg(:,Ac)
               hl(a)   = hg(Ac)
            END DO
!     Add Neumann BCs contribution to the LHS/RHS
            CALL BCONSTRUCT(lFa, yl, hl, ptr, e)
         END DO
      END IF

      RETURN
      END SUBROUTINE SETBCNEUL
!####################################################################
!     Here traction vector is applied for the outlets
      SUBROUTINE SETBCTRAC(Yg)

      USE COMMOD

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo)

      INTEGER iFa, iBc, iM

      DO iBc=1, eq(cEq)%nBc
         IF (.NOT.BTEST(eq(cEq)%bc(iBc)%bType,bType_trac)) CYCLE
         iFa = eq(cEq)%bc(iBc)%iFa
         iM  = eq(cEq)%bc(iBc)%iM
         IF (eq(cEq)%phys .NE. phys_preSt) err = "Traction BC is "//
     2      "applied for prestress equation only"
         CALL SETBCTRACL(eq(cEq)%bc(iBc), msh(iM)%fa(iFa), Yg)
      END DO

      RETURN
      END SUBROUTINE SETBCTRAC
!--------------------------------------------------------------------
      SUBROUTINE SETBCTRACL(lBc, lFa, Yg)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo)

      INTEGER :: a, Ac, e, g, iM, nNo, eNoN, cPhys
      REAL(KIND=8) :: w, nV(nsd), y(tDof), Jac, h(nsd)

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: N(:), yl(:,:), hl(:,:),hg(:,:),
     2   tmpA(:,:), lR(:,:), lK(:,:,:)

      iM   = lFa%iM
      nNo  = lFa%nNo
      eNoN = lFa%eNoN

      IF (.NOT.ALLOCATED(lBc%gm)) err = "Correction needed in SETBCTRAC"

!     Geting the contribution of traction BC
      ALLOCATE(tmpA(nsd,nNo), hg(nsd,nNo))
      CALL IGBC(lBc%gm, tmpA, hg)
      DEALLOCATE(hg)

      ALLOCATE(hg(nsd,tnNo), hl(nsd,eNoN), yl(tDof,eNoN), ptr(eNoN),
     2   N(eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Transforming it to a unified format
      DO a=1, nNo
         Ac       = lFa%gN(a)
         hg(:,Ac) = tmpA(:,a)
      END DO
      DEALLOCATE(tmpA)

!     Constructing LHS/RHS contribution and assembiling them
      DO e=1, lFa%nEl
         DO a=1, eNoN
            Ac      = lFa%IEN(a,e)
            ptr(a)  = Ac
            yl(:,a) = Yg(:,Ac)
            hl(:,a) = hg(:,Ac)
         END DO

!     Add Traction BCs contribution to the LHS/RHS
         lK = 0D0
         lR = 0D0
         cDmn  = DOMAIN(msh(iM), cEq, lFa%gE(e))
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_preSt) CYCLE

         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)

         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g)*Jac
            N   = lFa%N(:,g)

            h = 0D0
            y = 0D0
            DO a=1, eNoN
               h(:) = h(:) + N(a)*hl(:,a)
               y(:) = y(:) + N(a)*yl(:,a)
            END DO

            DO a=1, eNoN
               lR(:,a) = lR(:,a) - w*N(a)*h(:)
            END DO
         END DO

#ifdef WITH_TRILINOS
         IF (useTrilinosAssemAndLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO

      RETURN
      END SUBROUTINE SETBCTRACL
!####################################################################
!     Weak treatment of Dirichlet boundary conditions
      SUBROUTINE SETBCDIRW(Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER :: iBc, iFa, iM

      DO iBc=1, eq(cEq)%nBc
         iM  = eq(cEq)%bc(iBc)%iM
         iFa = eq(cEq)%bc(iBc)%iFa
         IF (.NOT.eq(cEq)%bc(iBc)%weakDir) CYCLE
         CALL SETBCDIRWL(eq(cEq)%bc(iBc), msh(iM), msh(iM)%fa(iFa), Yg,
     2      Dg)
      END DO

      RETURN
      END SUBROUTINE SETBCDIRW
!--------------------------------------------------------------------
      SUBROUTINE SETBCDIRWL(lBc, lM, lFa, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      LOGICAL :: flag, eDir(maxnsd)
      INTEGER :: a, e, i, g, Ac, Ec, ss, ee, lDof, nNo, nEl, nG, eNoN,
     2   eNoNb, cPhys
      REAL(KIND=8) :: w, Jac, xp(nsd), xi(nsd), xi0(nsd), nV(nsd),
     2   ub(nsd), tauB(2), Ks(nsd,nsd)

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: N(:), Nb(:), Nxi(:,:), Nx(:,:),
     2   xl(:,:), xbl(:,:), yl(:,:), ubl(:,:), ubg(:,:), tmpA(:,:),
     3   tmpY(:,:), lR(:,:), lK(:,:,:)

      nNo   = lFa%nNo
      nEl   = lFa%nEl
      nG    = lFa%nG
      eNoNb = lFa%eNoN
      eNoN  = lM%eNoN

      tauB  = lBc%tauB

      ss    = eq(cEq)%s
      ee    = eq(cEq)%e
      IF (eq(cEq)%dof .EQ. nsd+1) ee = ee - 1
      eDir  = .FALSE.
      lDof  = 0
      DO i=1, nsd
         IF (lBc%eDrn(i) .NE. 0) THEN
            eDir(i) = .TRUE.
            lDof = lDof + 1
         END IF
      END DO
      IF (lDof .EQ. 0) lDof = ee - ss + 1

      ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
      CALL SETBCDIRL(lBc, lFa, tmpA, tmpY, lDof)
      IF (BTEST(lBc%bType,bType_impD)) tmpY(:,:) = tmpA(:,:)

      ALLOCATE(ubg(nsd,tnNo))
      ubg = 0D0
      IF (ANY(eDir)) THEN
         DO a=1, nNo
            Ac = lFa%gN(a)
            DO i=1, nsd
               lDof = 0
               IF (eDir(i)) THEN
                  lDof = lDof + 1
                  ubg(i,Ac) = tmpY(lDof,a)
               END IF
            END DO
         END DO
      ELSE
         DO a=1, nNo
            Ac = lFa%gN(a)
            ubg(:,Ac) = tmpY(:,a)
         END DO
      END IF
      DEALLOCATE(tmpA, tmpY)

      ALLOCATE(Nb(eNoNb), xbl(nsd,eNoNb), ubl(nsd,eNoNb))

      ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN), xl(nsd,eNoN),
     2   yl(tDof,eNoN), ptr(eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

      xi0 = 0D0
      DO g=1, lM%nG
         xi0 = xi0 + lM%xi(:,g)
      END DO
      xi0 = xi0 / REAL(lM%nG,KIND=8)

      DO e=1, nEl
         Ec = lFa%gE(e)
         cDmn  = DOMAIN(lM, cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys

         IF (cPhys .NE. phys_fluid) err =
     2      " Weak Dirichlet BC formulated for fluid phys only"

         DO a=1, eNoN
            Ac = lM%IEN(a,Ec)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            IF (mvMsh) xl(:,a) = xl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
            yl(:,a) = Yg(:,Ac)
         END DO

         DO a=1, eNoNb
            Ac = lFa%IEN(a,e)
            xbl(:,a) = x(:,Ac)
            ubl(:,a) = ubg(:,Ac)
            IF (mvMsh) xbl(:,a) = xbl(:,a) + Dg(nsd+2:2*nsd+1,Ac)
         END DO

         lK = 0D0
         lR = 0D0
         DO g=1, nG
            CALL GNNB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g) * Jac
            Nb  = lFa%N(:,g)

            xp = 0D0
            ub = 0D0
            DO a=1, eNoNb
               xp = xp + xbl(:,a)*Nb(a)
               ub = ub + ubl(:,a)*Nb(a)
            END DO

            xi = xi0
            CALL GETXI(lM%eType, eNoN, xl, xp, xi, flag)
            CALL GETGNN(nsd, lM%eType, eNoN, xi, N, Nxi)

            i = 0
            DO a=1, eNoN
               IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
            END DO
            IF (.NOT.flag .OR. i.NE.eNoN) THEN
               err = " Error in computing face shape function"
            END IF

            IF (g.EQ.1 .OR. .NOT.lM%lShpF)
     2         CALL GNN(eNoN, Nxi, xl, Nx, Jac, Ks)

            IF (nsd .EQ. 3) THEN
               CALL BWFLUID3D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
            ELSE
               CALL BWFLUID2D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
            END IF
         END DO

!     Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (useTrilinosAssemAndLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO

      RETURN
      END SUBROUTINE SETBCDIRWL
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
!####################################################################
! Below defines the SET_BC methods for the Coupled Momentum Method (CMM)
      SUBROUTINE SETBCCMM(Ag, Yg, Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER iFa, iBc, iM

      DO iBc=1, eq(cEq)%nBc
          IF(.NOT.BTEST(eq(cEq)%bc(iBc)%bType,bType_CMM)) CYCLE
          iFa = eq(cEq)%bc(iBc)%iFa
          iM = eq(cEq)%bc(iBc)%iM
          CALL SETBCCMML(msh(iM)%fa(iFa), Ag, Yg, Dg)
      END DO

      RETURN
      END SUBROUTINE SETBCCMM
!--------------------------------------------------------------------
!     This defines the pseudo-structural equations to solve on the
!     boundaries of the lateral surfaces for the CMM method. It
!     borrows heavily from the current implementation of the 2D
!     linear elasticity equations, modified for the CMM method
!     It then will add the contributions to the LHS and RHS matrices
      SUBROUTINE SETBCCMML(lFa, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER a, Ac, e, eNoN

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: al(:,:), yl(:,:), dl(:,:), xl(:,:),
     2  pS0l(:,:)

      eNoN = lFa%eNoN
      ALLOCATE(ptr(eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
     2   xl(nsd,eNoN), pS0l(nstd,eNoN))

      pS0l = 0D0
!     Constructing the CMM contributions to the LHS/RHS and
!     assembling them
      DO e=1, lFa%nEl
          DO a=1, lFa%eNoN
              Ac = lFa%IEN(a,e)
              ptr(a)  = Ac
              al(:,a) = Ag(:,Ac)
              yl(:,a) = Yg(:,Ac)
              dl(:,a) = Dg(:,Ac)
              xl(:,a) = x(:,Ac)
              IF(ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
          END DO

!     Add CMM BCs contributions to the LHS/RHS
          CALL CMM_CONSTRUCT(lFa, al, yl, dl, xl, pS0l, ptr, e)
      END DO

      RETURN
      END SUBROUTINE SETBCCMML
!####################################################################

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
!-----------------------------------------------------------------------
!
!     This routine embodies formulation for solving electrophysiology
!     model equation using operator-splitting method.
!
!-----------------------------------------------------------------------

!     This is for solving 3D electrophysiology diffusion equations
      PURE SUBROUTINE CEP3D (eNoN, w, N, Nx, al, yl, fNl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), fNl(nFn*nsd,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER a, b, i, iFn
      REAL(KIND=8) :: T1, amd, wl, Td, Tx(nsd), Diso, Dani(2),
     2   fl(nsd,nFn), D(nsd,nsd), DTx(nsd), DNx(nsd,eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      Dani = eq(cEq)%dmn(cDmn)%cep%Dani
      i    = eq(cEq)%s

      wl = w*T1

      fl = 0D0
      DO a=1, eNoN
         DO iFn=1, nFn
            b = (iFn-1)*nsd
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
            fl(3,iFn) = fl(3,iFn) + N(a)*fNl(b+3,a)
         END DO
      END DO

      D(:,:) = 0D0
      D(1,1) = Diso
      D(2,2) = Diso
      D(3,3) = Diso
      DO iFn=1, nFn
         D(1,1) = D(1,1) + Dani(iFn)*fl(1,iFn)*fl(1,iFn)
         D(1,2) = D(1,2) + Dani(iFn)*fl(1,iFn)*fl(2,iFn)
         D(1,3) = D(1,3) + Dani(iFn)*fl(1,iFn)*fl(3,iFn)

         D(2,1) = D(2,1) + Dani(iFn)*fl(2,iFn)*fl(1,iFn)
         D(2,2) = D(2,2) + Dani(iFn)*fl(2,iFn)*fl(2,iFn)
         D(2,3) = D(2,3) + Dani(iFn)*fl(2,iFn)*fl(3,iFn)

         D(3,1) = D(3,1) + Dani(iFn)*fl(3,iFn)*fl(1,iFn)
         D(3,2) = D(3,2) + Dani(iFn)*fl(3,iFn)*fl(2,iFn)
         D(3,3) = D(3,3) + Dani(iFn)*fl(3,iFn)*fl(3,iFn)
      END DO

      Td = 0D0
      Tx = 0D0
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)

         DNx(1,a) = D(1,1)*Nx(1,a) + D(1,2)*Nx(2,a) + D(1,3)*Nx(3,a)
         DNx(2,a) = D(2,1)*Nx(1,a) + D(2,2)*Nx(2,a) + D(2,3)*Nx(3,a)
         DNx(3,a) = D(3,1)*Nx(1,a) + D(3,2)*Nx(2,a) + D(3,3)*Nx(3,a)
      END DO

      DTx(1) = D(1,1)*Tx(1) + D(1,2)*Tx(2) + D(1,3)*Tx(3)
      DTx(2) = D(2,1)*Tx(1) + D(2,2)*Tx(2) + D(2,3)*Tx(3)
      DTx(3) = D(3,1)*Tx(1) + D(3,2)*Tx(2) + D(3,3)*Tx(3)

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td
     2      + Nx(1,a)*DTx(1) + Nx(2,a)*DTx(2) + Nx(3,a)*DTx(3))

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b) + Nx(3,a)*DNx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP3D
!-----------------------------------------------------------------------
!     This is for solving 2D electrophysiology diffusion equation
      PURE SUBROUTINE CEP2D (eNoN, w, N, Nx, al, yl, fNl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), fNl(nFn*nsd,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER a, b, i, iFn
      REAL(KIND=8) :: T1, amd, wl, Td, Tx(nsd), Diso, Dani(2),
     2   fl(nsd,nFn), D(nsd,nsd), DTx(nsd), DNx(nsd,eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      Dani = eq(cEq)%dmn(cDmn)%cep%Dani
      i    = eq(cEq)%s
      wl   = w*T1

      fl = 0D0
      DO a=1, eNoN
         DO iFn=1, nFn
            b = (iFn-1)*nsd
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
            fl(3,iFn) = fl(3,iFn) + N(a)*fNl(b+3,a)
         END DO
      END DO

      D(:,:) = 0D0
      D(1,1) = Diso
      D(2,2) = Diso
      DO iFn=1, nFn
         D(1,1) = D(1,1) + Dani(iFn)*fl(1,iFn)*fl(1,iFn)
         D(1,2) = D(1,2) + Dani(iFn)*fl(1,iFn)*fl(2,iFn)

         D(2,1) = D(2,1) + Dani(iFn)*fl(2,iFn)*fl(1,iFn)
         D(2,2) = D(2,2) + Dani(iFn)*fl(2,iFn)*fl(2,iFn)
      END DO

      Td = 0D0
      Tx = 0D0
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)

         DNx(1,a) = D(1,1)*Nx(1,a) + D(1,2)*Nx(2,a)
         DNx(2,a) = D(2,1)*Nx(1,a) + D(2,2)*Nx(2,a)
      END DO

      DTx(1) = D(1,1)*Tx(1) + D(1,2)*Tx(2)
      DTx(2) = D(2,1)*Tx(1) + D(2,2)*Tx(2)

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td
     2      + Nx(1,a)*DTx(1) + Nx(2,a)*DTx(2))

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP2D
!-----------------------------------------------------------------------
!     This is for solving 1D electrophysiology diffusion equation
!     for Purkinje fibers
      PURE SUBROUTINE CEP1D (eNoN, insd, w, N, Nx, al, yl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN, insd
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(insd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER a, b, i, iFn
      REAL(KIND=8) :: T1, amd, wl, Td, Tx, Diso, DNx(eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      i    = eq(cEq)%s
      wl   = w*T1

      Td = 0D0
      Tx = 0D0
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)
         Tx = Tx + Nx(1,a)*yl(i,a)
         DNx(a) = Diso*Nx(1,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td + Nx(1,a)*Diso*Tx)

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd + Nx(1,a)*DNx(b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP1D
!-----------------------------------------------------------------------
      PURE SUBROUTINE BCEP (eNoN, w, N, h, lR)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER :: a
      REAL(KIND=8) f

      f = w*h

!     Here the loop is started for constructing left and right hand side
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + N(a)*f
      END DO

      RETURN
      END SUBROUTINE BCEP
!#######################################################################
      SUBROUTINE CEPION(iEq, iDof)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iEq, iDof

      LOGICAL :: IPASS = .TRUE.
      INTEGER :: Ac, iDmn, cPhys, dID, nX, nXmax

      REAL(KIND=8), ALLOCATABLE :: Xl(:), Xion(:,:), sA(:), sF(:,:)

      SAVE IPASS, nXmax, Xion

!     Initialization step
      IF (IPASS) THEN
         IPASS = .FALSE.
!        Determine max. state variables for all domains
         nXmax = 0
         DO iDmn=1, eq(iEq)%nDmn
            cPhys = eq(iEq)%dmn(iDmn)%phys
            IF (cPhys .NE. phys_CEP) CYCLE

            nX = eq(iEq)%dmn(iDmn)%cep%nX
            IF (nX .GT. nXmax) nXmax = nX
         END DO

!        Initialize CEP model state variables
         ALLOCATE(Xion(nxMax,tnNo))
         Xion = 0D0
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sA(tnNo), sF(nXmax,tnNo))
            sA = 0D0
            sF = 0D0
            DO Ac=1, tnNo
               DO iDmn=1, eq(iEq)%nDmn
                  cPhys = eq(iEq)%dmn(iDmn)%phys
                  dID   = eq(iEq)%dmn(iDmn)%Id
                  IF (cPhys.NE.phys_CEP .OR. .NOT.BTEST(dmnId(Ac),dID))
     2                CYCLE
                  nX = eq(iEq)%dmn(iDmn)%cep%nX
                  ALLOCATE(Xl(nX))
                  CALL CEPINIT(eq(iEq)%dmn(iDmn)%cep, nX, Xl)
                  sA(Ac) = sA(Ac) + 1.0D0
                  sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
                  DEALLOCATE(Xl)
               END DO
            END DO
            CALL COMMU(sA)
            CALL COMMU(sF)
            DO Ac=1, tnNo
               IF (.NOT.ISZERO(sA(Ac)))
     2            Xion(:,Ac) = sF(:,Ac)/sA(Ac)
            END DO
            DEALLOCATE(sA, sF)
         ELSE
            DO Ac=1, tnNo
               IF (eq(iEq)%dmn(1)%phys .NE. phys_CEP) CYCLE
               nX = eq(iEq)%dmn(1)%cep%nX
               ALLOCATE(Xl(nX))
               CALL CEPINIT(eq(iEq)%dmn(1)%cep, nX, Xl)
               Xion(1:nX,Ac) = Xl(:)
               DEALLOCATE(Xl)
            END DO
         END IF
      ELSE
!        Copy action potential after diffusion as first state variable
         DO Ac=1, tnNo
            Xion(1,Ac) = Yo(iDof,Ac)
         END DO
      END IF

      IF (ALLOCATED(dmnId)) THEN
         ALLOCATE(sA(tnNo), sF(nXmax,tnNo))
         sA = 0D0
         sF = 0D0
         DO Ac=1, tnNo
            DO iDmn=1, eq(iEq)%nDmn
               cPhys = eq(iEq)%dmn(iDmn)%phys
               dID   = eq(iEq)%dmn(iDmn)%Id
               IF (cPhys.NE.phys_CEP .OR. .NOT.BTEST(dmnId(Ac),dID))
     2             CYCLE
               nX = eq(iEq)%dmn(iDmn)%cep%nX
               ALLOCATE(Xl(nX))
               Xl(:) = Xion(1:nX,Ac)
               CALL CEPINTEG(eq(iEq)%dmn(iDmn)%cep, nX, Xl, time-dt, dt)
               sA(Ac) = sA(Ac) + 1.0D0
               sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
               DEALLOCATE(Xl)
            END DO
         END DO
         CALL COMMU(sA)
         CALL COMMU(sF)
         DO Ac=1, tnNo
            IF (.NOT.ISZERO(sA(Ac)))
     2         Xion(:,Ac) = sF(:,Ac)/sA(Ac)
         END DO
         DEALLOCATE(sA, sF)
      ELSE
         DO Ac=1, tnNo
            IF (eq(iEq)%dmn(1)%phys .NE. phys_CEP) CYCLE
            nX = eq(iEq)%dmn(1)%cep%nX
            ALLOCATE(Xl(nX))
            Xl(:) = Xion(1:nX,Ac)
            CALL CEPINTEG(eq(iEq)%dmn(iDmn)%cep, nX, Xl, time-dt, dt)
            Xion(1:nX,Ac) = Xl(:)
            DEALLOCATE(Xl)
         END DO
      END IF

      DO Ac=1, tnNo
         Yo(iDof,Ac) = Xion(1,Ac)
      END DO

      RETURN
      END SUBROUTINE CEPION
!#######################################################################
      SUBROUTINE CEPINIT(cep, nX, X)
      USE CEPMOD
      IMPLICIT NONE
      TYPE(cepModelType), INTENT(IN) :: cep
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(OUT) :: X(nX)

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         CALL AP_INIT(nX, X)

      CASE (cepModel_FN)
         CALL FN_INIT(nX, X)

      CASE (cepModel_TTP)
         CALL TTP_INIT(nX, X)

      END SELECT

      RETURN
      END SUBROUTINE CEPINIT
!-----------------------------------------------------------------------
      SUBROUTINE CEPINTEG(cep, nX, X, t, dt)
      USE CEPMOD
      USE UTILMOD, ONLY : eps
      IMPLICIT NONE
      TYPE(cepModelType), INTENT(IN) :: cep
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(IN) :: t, dt
      REAL(KIND=8), INTENT(INOUT) :: X(nX)

      REAL(KIND=8) :: Ts, Te, Istim

      INTEGER, ALLOCATABLE :: IPAR(:)
      REAL(KIND=8), ALLOCATABLE :: RPAR(:)

      Ts = cep%Istim%Ts + FLOOR(t/cep%Istim%Tp)
      Te = Ts + cep%Istim%Td
      IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
         Istim = cep%Istim%A
      ELSE
         Istim = 0D0
      END IF

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            CALL AP_INTEGFE(nX, X, t, dt, Istim)
         CASE (tIntType_RK4)
            CALL AP_INTEGRK(nX, X, t, dt, Istim)
         CASE (tIntType_CN2)
            CALL AP_INTEGCN2(nX, X, t, dt, Istim, IPAR, RPAR)
         END SELECT

      CASE (cepModel_FN)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            CALL FN_INTEGFE(nX, X, t, dt, Istim)
         CASE (tIntType_RK4)
            CALL FN_INTEGRK(nX, X, t, dt, Istim)
         CASE (tIntType_CN2)
            CALL FN_INTEGCN2(nX, X, t, dt, Istim, IPAR, RPAR)
         END SELECT

      CASE (cepModel_TTP)
         ALLOCATE(IPAR(2), RPAR(18))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            CALL TTP_INTEGFE(nX, X, t, dt, Istim, RPAR)
         CASE (tIntType_RK4)
            CALL TTP_INTEGRK(nX, X, t, dt, Istim, RPAR)
         CASE (tIntType_CN2)
            CALL TTP_INTEGCN2(nX, X, t, dt, Istim, IPAR, RPAR)
         END SELECT
      END SELECT

      IF (ISNAN(X(1))) THEN
         WRITE(*,'(A)') " NaN occurence. Aborted!"
         CALL STOPSIM()
      END IF

      DEALLOCATE(IPAR, RPAR)

      RETURN
      END SUBROUTINE CEPINTEG
!####################################################################

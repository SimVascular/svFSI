!-----------------------------------------------------------------------
!
!     This is for solving electrophysiology model equation.
!
!-----------------------------------------------------------------------
!#######################################################################
!     This is for solving 3D electrophysiology diffusion equations
      PURE SUBROUTINE CEP3D (eNoN, w, N, Nx, al, yl, fNl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), fNl(nsd,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER i, a, b
      REAL(KIND=8) :: T1, amd, wl, Td, Tx(nsd), Diso, Dani, fg(nsd),
     2   D(nsd,nsd), DTx(nsd), DNx(nsd,eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      Dani = eq(cEq)%dmn(cDmn)%cep%Dani
      i    = eq(cEq)%s

      wl = w*T1

      fg = 0D0
      DO a=1, eNoN
         fg = fg + N(a)*fNl(:,a)
      END DO

      D(1,1) = Diso + Dani*fg(1)*fg(1)
      D(1,2) = Dani*fg(1)*fg(2)
      D(1,3) = Dani*fg(1)*fg(3)

      D(2,1) = Dani*fg(2)*fg(1)
      D(2,2) = Diso + Dani*fg(2)*fg(2)
      D(2,3) = Dani*fg(2)*fg(3)

      D(3,1) = Dani*fg(3)*fg(1)
      D(3,2) = Dani*fg(3)*fg(2)
      D(3,3) = Diso + Dani*fg(3)*fg(3)

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

         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b) + Nx(3,a)*DNx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP3D
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
!-----------------------------------------------------------------------
!     This is for solving 2D electrophysiology diffusion equation
      PURE SUBROUTINE CEP2D (eNoN, w, N, Nx, al, yl, fNl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), fNl(nsd,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER i, a, b
      REAL(KIND=8) :: T1, amd, wl, Td, Tx(nsd), Diso, Dani, fg(nsd),
     2   D(nsd,nsd), DTx(nsd), DNx(nsd,eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      Dani = eq(cEq)%dmn(cDmn)%cep%Dani
      i    = eq(cEq)%s
      wl   = w*T1

      fg = 0D0
      DO a=1, eNoN
         fg = fg + N(a)*fNl(:,a)
      END DO

      D(1,1) = Diso + Dani*fg(1)*fg(1)
      D(1,2) = Dani*fg(1)*fg(2)

      D(2,1) = Dani*fg(2)*fg(1)
      D(2,2) = Diso + Dani*fg(2)*fg(2)

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

         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP2D
!#######################################################################
      SUBROUTINE CEPION(iEq, iDof)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iEq, iDof

      LOGICAL :: IPASS = .TRUE.
      INTEGER :: iM, iDmn, e, a, Ac, cPhys, nX, nXmax

      INTEGER, ALLOCATABLE :: incNd(:,:)
      REAL(KIND=8), ALLOCATABLE :: Xl(:), Xion(:,:), sA(:), sF(:,:)

      SAVE IPASS, nXmax, incNd, Xion

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

!        Determine nodes of a domain where CEP model is solved
         ALLOCATE(incNd(eq(iEq)%nDmn,tnNo))
         incNd = 0
         DO iM=1, nMsh
            DO e=1, msh(iM)%nEl
               iDmn = DOMAIN(msh(iM), iEq, e)
               cPhys = eq(iEq)%dmn(iDmn)%phys
               IF (cPhys .NE. phys_CEP) CYCLE

               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  incNd(iDmn, Ac) = 1
               END DO
            END DO
         END DO

!        Initialize CEP model state variables
         ALLOCATE(sA(tnNo), sF(nXmax,tnNo))
         sA = 0D0
         sF = 0D0
         DO iDmn=1, eq(iEq)%nDmn
            cPhys = eq(iEq)%dmn(iDmn)%phys
            IF (cPhys .NE. phys_CEP) CYCLE

            nX = eq(iEq)%dmn(iDmn)%cep%nX
            ALLOCATE(Xl(nX))
            CALL CEPINIT(eq(iEq)%dmn(iDmn)%cep, nX, Xl)

            DO Ac=1, tnNo
               IF (incNd(iDmn,Ac) .EQ. 0) CYCLE
               sA(Ac) = sA(Ac) + 1.0D0
               sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
            END DO
            DEALLOCATE(Xl)
         END DO

         CALL COMMU(sA)
         CALL COMMU(sF)

         ALLOCATE(Xion(nxMax,tnNo))
         Xion = 0D0
         DO Ac=1, tnNo
            IF (.NOT.ISZERO(sA(Ac))) THEN
               Xion(:,Ac) = sF(:,Ac)/sA(Ac)
            END IF
         END DO
         DEALLOCATE(sA, sF)
      ELSE
!        Copy action potential after diffusion as first state variable
         DO Ac=1, tnNo
            Xion(1,Ac) = Yo(iDof,Ac)
         END DO
      END IF

      ALLOCATE(sA(tnNo), sF(nXmax,tnNo))
      sA = 0D0
      sF = 0D0
      DO iDmn=1, eq(iEq)%nDmn
         cPhys = eq(iEq)%dmn(iDmn)%phys
         IF (cPhys.NE.phys_CEP) CYCLE

         nX = eq(iEq)%dmn(iDmn)%cep%nX
         ALLOCATE(Xl(nX))
         DO Ac=1, tnNo
            IF (incNd(iDmn,Ac) .EQ. 0) CYCLE
            Xl(:) = Xion(1:nX,Ac)
            CALL CEPINTEG(eq(iEq)%dmn(iDmn)%cep, nX, Xl, time-dt, dt)
            sA(Ac)      = sA(Ac) + 1.0D0
            sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
         END DO
         DEALLOCATE(Xl)
      END DO

      CALL COMMU(sA)
      CALL COMMU(sF)

      DO Ac=1, tnNo
         IF (.NOT.ISZERO(sA(Ac))) THEN
            Xion(:,Ac) = sF(:,Ac)/sA(Ac)
         END IF
      END DO
      DEALLOCATE(sA, sF)

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

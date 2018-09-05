!#######################################################################
      MODULE APMOD
      IMPLICIT NONE

      INTERFACE AP_INIT
         MODULE PROCEDURE :: AP_INIT0, AP_INITS, AP_INITV
      END INTERFACE AP_INIT

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE AP_INIT0(nX, X)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)

      INCLUDE "PARAMS_AP.f"

      X(:) = 1D-3
      X(1) = Voffset

      RETURN
      END SUBROUTINE AP_INIT0
!-----------------------------------------------------------------------
      SUBROUTINE AP_INITS(nX, X, X0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: X0

      X(:) = X0

      RETURN
      END SUBROUTINE AP_INITS
!-----------------------------------------------------------------------
      SUBROUTINE AP_INITV(nX, X, X0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: X0(:)

      IF (SIZE(X0,1) .NE. nX) THEN
         STOP "ERROR: inconsistent array size (AP initialization)"
      END IF
      X(:) = X0(:)

      RETURN
      END SUBROUTINE AP_INITV
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE AP_INTEGFE(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=8) :: t, dt, f(nX), fext

      t    = Ts / Tscale
      dt   = Ti / Tscale
      fext = Istim * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
      CALL AP_GETF(nX, t, X, f, fext)
      X = X + dt*f
      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE AP_INTEGRK(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=8) :: t, trk, dt, fext, Xrk(nX,4), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      fext = Istim * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
!     RK4: 1st pass
      trk = t
      Xrk(:,1) = X(:)
      CALL AP_GETF(nX, trk, Xrk(:,1), frk(:,1), fext)

!     RK4: 2nd pass
      trk = t + dt/2.0D0
      Xrk(:,2) = X(:) + dt*frk(:,1)/2.0D0
      CALL AP_GETF(nX, trk, Xrk(:,2), frk(:,2), fext)

!     RK4: 3rd pass
      trk = t + dt/2.0D0
      Xrk(:,3) = X(:) + dt*frk(:,2)/2.0D0
      CALL AP_GETF(nX, trk, Xrk(:,3), frk(:,3), fext)

!     RK4: 4th pass
      trk = t + dt
      Xrk(:,4) = X(:) + dt*frk(:,3)
      CALL AP_GETF(nX, trk, Xrk(:,4), frk(:,4), fext)

      X(:) = X(:) + (dt/6.0D0) * ( frk(:,1) + 2.0D0*frk(:,2) +
     2   2.0D0*frk(:,3) + frk(:,4) )

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGRK
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE AP_INTEGCN2(nX, Xn, Ts, Ti, Istim, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      INTEGER, INTENT(INOUT) :: IPAR(2)
      REAL(KIND=8), INTENT(INOUT) :: Xn(nX), RPAR(2)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=8), PARAMETER :: eps = EPSILON(eps)

      INTEGER :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=8) :: t, dt, fext, atol, rtol, Xk(nX), fn(nX), fk(nX),
     2   rK(nX), Im(nX,nX), JAC(nX,nX), rmsA, rmsR

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      t     = Ts / Tscale
      dt    = Ti / Tscale
      fext  = Istim * Tscale / Vscale
      Xn(1) = (Xn(1) - Voffset)/Vscale
      Im    = MAT_ID(nX)

      CALL AP_GETF(nX, t, Xn, fn, fext)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL AP_GETF(nX, t, Xk, fk, fext)
         rK(:) = Xk(:) - Xn(:) - 0.5D0*dt*(fk(:) + fn(:))

         rmsA = 0D0
         rmsR = 0D0
         DO i=1, nX
            rmsA = rmsA + rK(i)**2.0D0
            rmsR = rmsR + ( rK(i) / (Xk(i)+eps) )**2.0D0
         END DO
         rmsA = SQRT(rmsA/REAL(nX,KIND=8))
         rmsR = SQRT(rmsR/REAL(nX,KIND=8))

         l1   = k .GT. itMax
         l2   = rmsA .LE. atol
         l3   = rmsR .LE. rtol
         IF (l1 .OR. l2 .OR. l3) EXIT

         CALL AP_GETJ(nX, t, Xk, JAC)
         JAC   = Im - 0.5D0*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL AP_GETF(nX, t, Xn, fn, fext)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE AP_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETF(n, t, X, f, fext)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: t, X(n), fext
      REAL(KIND=8), INTENT(OUT) :: f(n)

      INCLUDE "PARAMS_AP.f"

      f(1) = X(1)*(c*(X(1)-alpha)*(1.0D0-X(1)) - X(2)) + fext

      f(2) = (a + mu1*X(2)/(mu2 + X(1))) *
     2       (-X(2) - c*X(1)*(X(1) - b - 1.0D0))

      RETURN
      END SUBROUTINE AP_GETF
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETJ(n, t, X, JAC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: t, X(n)
      REAL(KIND=8), INTENT(OUT) :: JAC(n,n)

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=8) :: n1, n2, n3

      JAC(:,:) = 0.0D0

      n1 = X(1) - alpha
      n2 = 1.0D0 - X(1)
      JAC(1,1) = c * (n1*n2 + X(1)*(n2 - n1)) - X(2)

      JAC(1,2) = -X(1)

      n1 = mu1*X(2)/(mu2 + X(1))
      n2 = n1 / (mu2 + X(1))
      n3 = X(2) + c*X(1)*(X(1) - b - 1.0D0)

      JAC(2,1) = n2*n3 - c*(a+n1)*(2.0D0*X(1) - b - 1.0D0)

      n1 = mu1/(mu2 + X(1))
      n2 = a + n1*X(2)
      n3 = -n3
      JAC(2,2) = n1*n3 - n2

      RETURN
      END SUBROUTINE AP_GETJ
!-----------------------------------------------------------------------
      END MODULE APMOD
!#######################################################################


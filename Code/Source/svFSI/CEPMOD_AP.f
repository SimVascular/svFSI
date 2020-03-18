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
!     This module defines data structures for Aliev-Panfilov cellular
!     activation model for cardiac electrophysiology.
!
!-----------------------------------------------------------------------

      MODULE APMOD
      USE TYPEMOD
      IMPLICIT NONE

      INTERFACE AP_INIT
         MODULE PROCEDURE :: AP_INIT0, AP_INITS, AP_INITV
      END INTERFACE AP_INIT

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE AP_INIT0(nX, X)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)

      INCLUDE "PARAMS_AP.f"

      X(:) = 1.E-3_RKIND
      X(1) = Voffset

      RETURN
      END SUBROUTINE AP_INIT0
!-----------------------------------------------------------------------
      SUBROUTINE AP_INITS(nX, X, X0)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: X0

      X(:) = X0

      RETURN
      END SUBROUTINE AP_INITS
!-----------------------------------------------------------------------
      SUBROUTINE AP_INITV(nX, X, X0)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: X0(:)

      IF (SIZE(X0,1) .NE. nX) THEN
         STOP "ERROR: inconsistent array size (AP initialization)"
      END IF
      X(:) = X0(:)

      RETURN
      END SUBROUTINE AP_INITV
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE AP_INTEGFE(nX, X, Ts, Ti, Istim, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=RKIND) :: t, dt, f(nX), fext, Isac

      t    = Ts / Tscale
      dt   = Ti / Tscale
      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
      CALL AP_GETF(nX, X, f, fext)
      X = X + dt*f
      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE AP_INTEGRK(nX, X, Ts, Ti, Istim, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=RKIND) :: t, dt, dt6, fext, Isac, Xrk(nX), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      dt6  = dt / 6._RKIND

      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale
      X(1) = (X(1) - Voffset)/Vscale

!     RK4: 1st pass
      Xrk  = X
      CALL AP_GETF(nX, Xrk, frk(:,1), fext)

!     RK4: 2nd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,1)
      CALL AP_GETF(nX, Xrk, frk(:,2), fext)

!     RK4: 3rd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,2)
      CALL AP_GETF(nX, Xrk, frk(:,3), fext)

!     RK4: 4th pass
      Xrk  = X + dt*frk(:,3)
      CALL AP_GETF(nX, Xrk, frk(:,4), fext)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGRK
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE AP_INTEGCN2(nX, Xn, Ts, Ti, Istim, Ksac, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      INTEGER(KIND=IKIND), INTENT(INOUT) :: IPAR(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xn(nX), RPAR(2)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INTEGER(KIND=IKIND) :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: t, dt, fext, atol, rtol, Xk(nX), fn(nX),
     2   fk(nX), rK(nX), Im(nX,nX), JAC(nX,nX), rmsA, rmsR, Isac

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      t     = Ts / Tscale
      dt    = Ti / Tscale
      Isac  = Ksac * (Vrest - Xn(1))
      fext  = (Istim + Isac) * Tscale / Vscale
      Xn(1) = (Xn(1) - Voffset)/Vscale
      Im    = MAT_ID(nX)

      CALL AP_GETF(nX, Xn, fn, fext)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL AP_GETF(nX, Xk, fk, fext)
         rK(:) = Xk(:) - Xn(:) - 0.5_RKIND*dt*(fk(:) + fn(:))

         rmsA = 0._RKIND
         rmsR = 0._RKIND
         DO i=1, nX
            rmsA = rmsA + rK(i)**2._RKIND
            rmsR = rmsR + ( rK(i) / (Xk(i)+eps) )**2._RKIND
         END DO
         rmsA = SQRT(rmsA/REAL(nX, KIND=RKIND))
         rmsR = SQRT(rmsR/REAL(nX, KIND=RKIND))

         l1   = k .GT. itMax
         l2   = rmsA .LE. atol
         l3   = rmsR .LE. rtol
         IF (l1 .OR. l2 .OR. l3) EXIT

         CALL AP_GETJ(nX, Xk, JAC, Ksac*Tscale)
         JAC   = Im - 0.5_RKIND*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL AP_GETF(nX, Xn, fn, fext)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE AP_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETF(n, X, f, fext)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), fext
      REAL(KIND=RKIND), INTENT(OUT) :: f(n)

      INCLUDE "PARAMS_AP.f"

      f(1) = X(1)*(c*(X(1)-alpha)*(1._RKIND-X(1)) - X(2)) + fext

      f(2) = (a + mu1*X(2)/(mu2 + X(1))) *
     2       (-X(2) - c*X(1)*(X(1) - b - 1._RKIND))

      RETURN
      END SUBROUTINE AP_GETF
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETJ(n, X, JAC, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), Ksac
      REAL(KIND=RKIND), INTENT(OUT) :: JAC(n,n)

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=RKIND) :: n1, n2, n3

      JAC(:,:) = 0._RKIND

      n1 = X(1) - alpha
      n2 = 1._RKIND - X(1)
      JAC(1,1) = c * (n1*n2 + X(1)*(n2 - n1)) - X(2) - Ksac

      JAC(1,2) = -X(1)

      n1 = mu1*X(2)/(mu2 + X(1))
      n2 = n1 / (mu2 + X(1))
      n3 = X(2) + c*X(1)*(X(1) - b - 1._RKIND)

      JAC(2,1) = n2*n3 - c*(a+n1)*(2._RKIND*X(1) - b - 1._RKIND)

      n1 = mu1/(mu2 + X(1))
      n2 = a + n1*X(2)
      n3 = -n3
      JAC(2,2) = n1*n3 - n2

      RETURN
      END SUBROUTINE AP_GETJ
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model
      SUBROUTINE AP_ACTVSTRS(X, dt, Tact, epsX)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(OUT) :: epsX
      REAL(KIND=RKIND), INTENT(INOUT) :: Tact

      INCLUDE "PARAMS_AP.f"

      REAL(KIND=RKIND) :: nr

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX
      nr   = Tact + epsX*dt*eta_T*(X - Vrest)
      Tact = nr / (1._RKIND + epsX*dt)

      RETURN
      END SUBROUTINE AP_ACTVSTRS
!-----------------------------------------------------------------------
      END MODULE APMOD
!#######################################################################


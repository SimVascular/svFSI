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
!     This module defines data structures for Fitzhugh-Nagumo cellular
!     activation model for cardiac electrophysiology.
!
!-----------------------------------------------------------------------

      MODULE FNMOD
      IMPLICIT NONE

      INTERFACE FN_INIT
         MODULE PROCEDURE :: FN_INIT0, FN_INITS, FN_INITV
      END INTERFACE FN_INIT

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE FN_INIT0(nX, X)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)

      X(:) = 1D-3

      RETURN
      END SUBROUTINE FN_INIT0
!-----------------------------------------------------------------------
      SUBROUTINE FN_INITS(nX, X, X0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: X0

      X(:) = X0

      RETURN
      END SUBROUTINE FN_INITS
!-----------------------------------------------------------------------
      SUBROUTINE FN_INITV(nX, X, X0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: X0(:)

      IF (SIZE(X0,1) .NE. nX) THEN
         STOP "ERROR: inconsistent array size (AP initialization)"
      END IF
      X(:) = X0(:)

      RETURN
      END SUBROUTINE FN_INITV
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE FN_INTEGFE(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_FN.f"

      REAL(KIND=8) :: t, dt, f(nX), fext

      t    = Ts / Tscale
      dt   = Ti / Tscale
      fext = Istim * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
      CALL FN_GETF(nX, t, X, f, fext)
      X(:) = X(:) + dt*f(:)
      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE FN_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE FN_INTEGRK(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(INOUT) :: X(nX)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_FN.f"

      REAL(KIND=8) :: t, trk, dt, fext, Xrk(nX,4), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      fext = Istim * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
!     RK4: 1st pass
      trk = t
      Xrk(:,1) = X(:)
      CALL FN_GETF(nX, trk, Xrk(:,1), frk(:,1), fext)

!     RK4: 2nd pass
      trk = t + dt/2.0D0
      Xrk(:,2) = X(:) + dt*frk(:,1)/2.0D0
      CALL FN_GETF(nX, trk, Xrk(:,2), frk(:,2), fext)

!     RK4: 3rd pass
      trk = t + dt/2.0D0
      Xrk(:,3) = X(:) + dt*frk(:,2)/2.0D0
      CALL FN_GETF(nX, trk, Xrk(:,3), frk(:,3), fext)

!     RK4: 4th pass
      trk = t + dt
      Xrk(:,4) = X(:) + dt*frk(:,3)
      CALL FN_GETF(nX, trk, Xrk(:,4), frk(:,4), fext)

      X(:) = X(:) + (dt/6.0D0) * ( frk(:,1) + 2.0D0*frk(:,2) +
     2   2.0D0*frk(:,3) + frk(:,4) )

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE FN_INTEGRK
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE FN_INTEGCN2(nX, Xn, Ts, Ti, Istim, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      INTEGER, INTENT(INOUT) :: IPAR(2)
      REAL(KIND=8), INTENT(INOUT) :: Xn(nX), RPAR(2)
      REAL(KIND=8), INTENT(IN) :: Ts, Ti, Istim

      INCLUDE "PARAMS_FN.f"

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

      CALL FN_GETF(nX, t, Xn, fn, fext)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL FN_GETF(nX, t, Xk, fk, fext)
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

         CALL FN_GETJ(nX, t, Xk, JAC)
         JAC   = Im - 0.5D0*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL FN_GETF(nX, t, Xn, fn, fext)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE FN_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE FN_GETF(n, t, X, f, fext)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: t, X(n), fext
      REAL(KIND=8), INTENT(OUT) :: f(n)

      INCLUDE "PARAMS_FN.f"

      f(1) = c * ( X(1)*(X(1)-alpha)*(1.0D0-X(1)) - X(2) ) + fext

      f(2) = X(1) - b*X(2) + a

      RETURN
      END SUBROUTINE FN_GETF
!-----------------------------------------------------------------------
      SUBROUTINE FN_GETJ(n, t, X, JAC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(IN) :: t, X(n)
      REAL(KIND=8), INTENT(OUT) :: JAC(n,n)

      INCLUDE "PARAMS_FN.f"

      REAL(KIND=8) :: n1, n2

      n1 = -3.0D0*X(1)**2.0D0
      n2 = 2.0D0*(1.0D0+alpha)*X(1)
      JAC(1,1) = c * (n1 + n2 - 1.0D0)
      JAC(1,2) = -c
      JAC(2,1) = 1.0D0
      JAC(2,2) = -b

      RETURN
      END SUBROUTINE FN_GETJ
!-----------------------------------------------------------------------
      END MODULE FNMOD
!#######################################################################


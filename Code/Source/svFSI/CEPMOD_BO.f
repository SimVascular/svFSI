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
!     This module defines data structures for Bueno-Orovio cellular
!     activation model for cardiac electrophysiology.
!
!--------------------------------------------------------------------

      MODULE BOMOD
      USE TYPEMOD
      IMPLICIT NONE

      PRIVATE :: ISZERO, STEP, DELTA

      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE BO_INIT(nX, X)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)

      INCLUDE "PARAMS_BO.f"

      X(1) = Voffset
      X(2) = 1._RKIND
      X(3) = 1._RKIND
      X(4) = 0._RKIND

      RETURN
      END SUBROUTINE BO_INIT
!--------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE BO_INTEGFE(imyo, nX, X, Ts, Ti, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), RPAR(5)

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: t, dt, f(nX), fext, Isac

      t    = Ts / Tscale
      dt   = Ti / Tscale

      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale

      CALL BO_GETF(imyo, nX, X, f, fext, RPAR)
      X(:) = X(:) + dt*f(:)

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE BO_INTEGFE
!--------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE BO_INTEGRK(imyo, nX, X, Ts, Ti, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), RPAR(5)

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: t, dt, dt6, fext, Isac, Xrk(nX), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      dt6  = dt/6._RKIND

      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale
      X(1) = (X(1) - Voffset)/Vscale

!     RK4: 1st pass
      Xrk  = X
      CALL BO_GETF(imyo, nX, Xrk, frk(:,1), fext, RPAR)

!     RK4: 2nd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,1)
      CALL BO_GETF(imyo, nX, Xrk, frk(:,2), fext, RPAR)

!     RK4: 3rd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,2)
      CALL BO_GETF(imyo, nX, Xrk, frk(:,3), fext, RPAR)

!     RK4: 4th pass
      Xrk  = X + dt*frk(:,3)
      CALL BO_GETF(imyo, nX, Xrk, frk(:,4), fext, RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE BO_INTEGRK
!--------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE BO_INTEGCN2(imyo, nX, Xn, Ts, Ti, Istim, Ksac, IPAR,
     2   RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX
      INTEGER(KIND=IKIND), INTENT(INOUT) :: IPAR(2)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac
      REAL(KIND=RKIND), INTENT(INOUT) :: Xn(nX), RPAR(5)

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INTEGER(KIND=IKIND) :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: t, dt, fext, atol, rtol, Xk(nX), fn(nX),
     2   fK(nX), rK(nX), Im(nX,nX), JAC(nX,nX), rmsA, rmsR, Isac

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      t     = Ts / Tscale
      dt    = Ti / Tscale
      Isac  = Ksac * (Vrest - Xn(1))
      fext  = (Istim + Isac) * Tscale / Vscale

      Xn(1) = (Xn(1) - Voffset)/Vscale
      Im    = MAT_ID(nX)

      CALL BO_GETF(imyo, nX, Xn, fn, fext, RPAR)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL BO_GETF(imyo, nX, Xk, fk, fext, RPAR)
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

         CALL BO_GETJ(imyo, nX, Xk, JAC)
         JAC   = Im - 0.5_RKIND*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL BO_GETF(imyo, nX, Xn, fn, fext, RPAR)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE BO_INTEGCN2
!--------------------------------------------------------------------
      SUBROUTINE BO_GETF(i, n, X, f, fext, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: i, n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), fext
      REAL(KIND=RKIND), INTENT(OUT) :: f(n)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(5)

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: u, v, w, s, H_uv, H_uw, H_umv, H_uo, taum_v,
     2   taum_w, tau_so, tau_s, tau_o, v_inf, w_inf, I_fi, I_so, I_si

!     Create local copies of the 4 state variables
      u = X(1)
      v = X(2)
      w = X(3)
      s = X(4)

!     Define step functions
      H_uv  = STEP(u - theta_v(i))
      H_uw  = STEP(u - theta_w(i))
      H_umv = STEP(u - thetam_v(i))
      H_uo  = STEP(u - theta_o(i))

!     Define additional constants
      taum_v = (1._RKIND-H_umv)*taum_v1(i) + H_umv*taum_v2(i)
      taum_w = taum_w1(i) + 0.5_RKIND*(taum_w2(i)-taum_w1(i))*
     2   (1._RKIND + TANH(km_w(i)*(u-um_w(i))))
      tau_so = tau_so1(i) + 0.5_RKIND*(tau_so2(i)-tau_so1(i))*
     2   (1._RKIND+DTANH(k_so(i)*(u-u_so(i))))
      tau_s  = (1._RKIND-H_uw)*tau_s1(i) + H_uw*tau_s2(i)
      tau_o  = (1._RKIND-H_uo)*tau_o1(i) + H_uo*tau_o2(i)
      v_inf  = (1._RKIND-H_umv)
      w_inf  = (1._RKIND-H_uo)*(1._RKIND - u/tau_winf(i)) +
     2   H_uo*ws_inf(i)

!     Compute RHS of state variable equations
      I_fi = -v*H_uv*(u-theta_v(i))*(u_u(i) - u)/tau_fi(i)
      I_so =  (u-u_o(i))*(1._RKIND-H_uw)/tau_o + H_uw/tau_so
      I_si = -H_uw*w*s/tau_si(i)

      f(1) = -(I_fi + I_so + I_si + fext)

      f(2) = (1._RKIND-H_uv)*(v_inf-v)/taum_v - H_uv*v/taup_v(i)

      f(3) = (1._RKIND-H_uw)*(w_inf-w)/taum_w - H_uw*w/taup_w(i)

      f(4) = (0.5_RKIND*(1._RKIND + TANH(k_s(i)*(u-u_s(i))))-s)/tau_s

      RPAR(3) = I_fi
      RPAR(4) = I_so
      RPAR(5) = I_si

      RETURN
      END SUBROUTINE BO_GETF
!--------------------------------------------------------------------
      SUBROUTINE BO_GETJ(i, n, X, JAC)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: i, n
      REAL(KIND=RKIND), INTENT(IN) :: X(n)
      REAL(KIND=RKIND), INTENT(OUT) :: JAC(n,n)

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: u, v, w, s, H_uv, H_uw, H_umv, H_uo, D_uw,
     2   D_uv, taum_v, taum_w, tau_so, tau_s, tau_o, v_inf, w_inf, n1,
     3   n2, n3

!     Create local variables
      u = X(1)
      v = X(2)
      w = X(3)
      s = X(4)

!     Define step functions
      H_uv  = STEP(u - theta_v(i))
      H_uw  = STEP(u - theta_w(i))
      H_umv = STEP(u - thetam_v(i))
      H_uo  = STEP(u - theta_o(i))

!     Define delta functions
      D_uw  = DELTA(u - theta_w(i))
      D_uv  = DELTA(u - theta_v(i))

!     Define additional constants
      taum_v = (1._RKIND-H_umv)*taum_v1(i) + H_umv*taum_v2(i)
      taum_w = taum_w1(i) + 0.5_RKIND*(taum_w2(i)-taum_w1(i))*
     2   (1._RKIND + TANH(km_w(i)*(u-um_w(i))))
      tau_so = tau_so1(i) + 0.5_RKIND*(tau_so2(i)-tau_so1(i))*
     2   (1._RKIND+DTANH(k_so(i)*(u-u_so(i))))
      tau_s  = (1._RKIND-H_uw)*tau_s1(i) + H_uw*tau_s2(i)
      tau_o  = (1._RKIND-H_uo)*tau_o1(i) + H_uo*tau_o2(i)
      v_inf  = (1._RKIND-H_umv)
      w_inf  = (1._RKIND-H_uo)*(1._RKIND - u/tau_winf(i)) +
     2   H_uo*ws_inf(i)

!     Define Jacobian
      JAC(:,:) = 0._RKIND

      n1 = v*H_uv*(u_u(i) + theta_v(i) - 2._RKIND*u)/tau_fi(i)
      n2 = -(1._RKIND-H_uw)/tau_fi(i)
      n3 = (-1._RKIND/tau_so + (theta_w(i)-u_o(i))/tau_o +
     2   w*s/tau_si(i))*D_uw
      JAC(1,1) = n1 + n2 + n3

      JAC(1,2) = H_uv*(u-theta_v(i))*(u_u(i)-u)/tau_fi(i)

      n1 = H_uw/tau_si(i)
      JAC(1,3) = n1*s
      JAC(1,4) = n1*w

      n1 = -1._RKIND/taum_v
      n2 = -1._RKIND/taup_v(i)
      JAC(2,1) = ((v_inf-v)*n1 + v*n2)*D_uv
      JAC(2,2) = (1._RKIND-H_uv)*n1 + H_uv*n2

      n1 = -1._RKIND/taum_w
      n2 = -1._RKIND/taup_w(i)
      JAC(3,1) = ((w_inf-w)*n1 + w*n2)*D_uw
      JAC(3,3) = (1._RKIND-H_uw)*n1 + H_uw*n2

      n1 = COSH(k_s(i)*(u-u_s(i)))
      n2 = 1._RKIND/(n1*n1)
      n3 = 1._RKIND/tau_s
      JAC(4,1) = 0.5_RKIND*k_s(i)*n2*n3
      JAC(4,4) = -n3

      RETURN
      END SUBROUTINE BO_GETJ
!--------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model
      SUBROUTINE BO_ACTVSTRS(X, dt, Tact, epsX)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(OUT) :: epsX
      REAL(KIND=RKIND), INTENT(INOUT) :: Tact

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: nr

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX
      nr   = Tact + epsX*dt*eta_T*(X - Vrest)
      Tact = nr / (1._RKIND + epsX*dt)

      RETURN
      END SUBROUTINE BO_ACTVSTRS
!--------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and slow inward current variable (s)
      SUBROUTINE BO_ACTVSTRN(c, I4f, dt, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c, I4f, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      INCLUDE "PARAMS_BO.f"

      REAL(KIND=RKIND) :: SL, Fa, rtmp

!     fiber length
      SL = I4f * SL0

!     Sacromere force-length relationship
      IF (SL.GE.SLmin .AND. SL.LE.SLmax) THEN
         SL = 0.5_RKIND*f0 + fc1*COS(SL) + fs1*SIN(SL) +
     2      fc2*COS(2._RKIND*SL) + fs2*SIN(2._RKIND*SL)  +
     3      fc3*COS(3._RKIND*SL) + fs3*SIN(3._RKIND*SL)
      ELSE
         SL = 0._RKIND
      END IF

!     Active force
      Fa   = alFa * (c-c0)*(c-c0) * SL

      rtmp = 2._RKIND*I4f*(1._RKIND/(1._RKIND+gf)**3._RKIND - 1._RKIND)
      gf = gf + dt*(Fa + rtmp)/(mu_C * c * c)

      RETURN
      END SUBROUTINE BO_ACTVSTRN
!--------------------------------------------------------------------
      FUNCTION STEP(r)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: r
      REAL(KIND=RKIND) STEP

      IF (r .LT. 0) THEN
         STEP = 0._RKIND
      ELSE
         STEP = 1._RKIND
      END IF

      RETURN
      END FUNCTION STEP
!--------------------------------------------------------------------
      FUNCTION DELTA(r)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: r
      REAL(KIND=RKIND) DELTA

      DELTA = 0._RKIND
      IF (ISZERO(r)) DELTA = 1._RKIND

      RETURN
      END FUNCTION DELTA
!--------------------------------------------------------------------
      FUNCTION ISZERO(ia)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: ia
      LOGICAL ISZERO

      REAL(KIND=RKIND), PARAMETER :: epsil = EPSILON(epsil)
      REAL(KIND=RKIND) a, b, nrm

      a   = ABS(ia)
      b   = 0._RKIND
      nrm = MAX(a,epsil)

      ISZERO = .FALSE.
      IF ((a-b)/nrm .LT. 1._RKIND*epsil) ISZERO = .TRUE.

      RETURN
      END FUNCTION ISZERO
!--------------------------------------------------------------------
      END MODULE BOMOD
!####################################################################


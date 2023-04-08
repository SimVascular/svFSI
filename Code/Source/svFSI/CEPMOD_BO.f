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
!     activation model for cardiac electrophysiology. Both active stress
!     and active strain models are used for excitation-contraction
!     coupling.
!
!     Reference for Bueno-Orovio electrophysiology model:
!        Bueno-Orovio, A., Cherry, E. M., & Fenton, F. H. (2008).
!        Minimal model for human ventricular action potentials in tissue
!        Journal of Theoretical Biology, 253(3), 544–560.
!        https://doi.org/10.1016/j.jtbi.2008.03.029
!
!     References for active stress/active strain models:
!        Garcia-blanco, E., et al. (2019). A new computational framework
!        for electro-activation in cardiac mechanics. Computer Methods
!        in Applied Mechanics and Engineering.
!        https://doi.org/10.1016/j.cma.2019.01.042
!
!        Simone Rossi, et al. (2012). Orthotropic active strain models
!        for the numerical simulation of cardiac biomechanics.
!        International Journal for Numerical Methods in Biomedical
!        Engineering, 28, 761–788. https://doi.org/10.1002/cnm.2473
!
!--------------------------------------------------------------------

      MODULE BOMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INCLUDE "PARAMS_BO.f"

      PUBLIC :: BO_INIT
      PUBLIC :: BO_READPARFF
      PUBLIC :: BO_INTEGFE
      PUBLIC :: BO_INTEGRK
      PUBLIC :: BO_INTEGCN2
      PUBLIC :: BO_ACTVSTRS_FE
      PUBLIC :: BO_ACTVSTRS_RK
      PUBLIC :: BO_ACTVSTRS_BE
      PUBLIC :: BO_ACTVSTRN_FE
      PUBLIC :: BO_ACTVSTRN_RK
      PUBLIC :: BO_ACTVSTRN_BE

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE BO_INIT(nX, X)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)

      X(1) = Voffset
      X(2) = 1.0_RKIND
      X(3) = 1.0_RKIND
      X(4) = 1.E-5_RKIND

      RETURN
      END SUBROUTINE BO_INIT
!-----------------------------------------------------------------------
      SUBROUTINE BO_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

!     Cellular activation model parameters
      CALL GETRVEC(fid, "u_o", 3, u_o)
      CALL GETRVEC(fid, "u_u", 3, u_u)
      CALL GETRVEC(fid, "theta_v", 3, theta_v)
      CALL GETRVEC(fid, "theta_w", 3, theta_w)
      CALL GETRVEC(fid, "thetam_v", 3, thetam_v)
      CALL GETRVEC(fid, "theta_o", 3, theta_o)
      CALL GETRVEC(fid, "taum_v1", 3, taum_v1)
      CALL GETRVEC(fid, "taum_v2", 3, taum_v2)
      CALL GETRVEC(fid, "taup_v", 3, taup_v)
      CALL GETRVEC(fid, "taum_w1", 3, taum_w1)
      CALL GETRVEC(fid, "taum_w2", 3, taum_w2)
      CALL GETRVEC(fid, "km_w", 3, km_w)
      CALL GETRVEC(fid, "um_w", 3, um_w)
      CALL GETRVEC(fid, "taup_w", 3, taup_w)
      CALL GETRVEC(fid, "tau_fi", 3, tau_fi)
      CALL GETRVEC(fid, "tau_o1", 3, tau_o1)
      CALL GETRVEC(fid, "tau_o2", 3, tau_o2)
      CALL GETRVEC(fid, "tau_so1", 3, tau_so1)
      CALL GETRVEC(fid, "tau_so2", 3, tau_so2)
      CALL GETRVEC(fid, "k_so", 3, k_so)
      CALL GETRVEC(fid, "u_so", 3, u_so)
      CALL GETRVEC(fid, "tau_s1", 3, tau_s1)
      CALL GETRVEC(fid, "tau_s2", 3, tau_s2)
      CALL GETRVEC(fid, "k_s", 3, k_s)
      CALL GETRVEC(fid, "u_s", 3, u_s)
      CALL GETRVEC(fid, "tau_si", 3, tau_si)
      CALL GETRVEC(fid, "tau_winf", 3, tau_winf)
      CALL GETRVEC(fid, "ws_inf", 3, ws_inf)

!     Electromechanics coupling parameters: active stress model
      CALL GETRVAL(fid, "Vrest", Vrest)
      CALL GETRVAL(fid, "Vcrit", Vcrit)
      CALL GETRVAL(fid, "K_T", K_T)
      CALL GETRVAL(fid, "eps_0", eps_0)
      CALL GETRVAL(fid, "eps_i", eps_i)
      CALL GETRVAL(fid, "xi_T", xi_T)

!     Electromechanics coupling parameters: active strain model
      CALL GETRVAL(fid, "alfa", alfa)
      CALL GETRVAL(fid, "c_0", c_0)
      CALL GETRVAL(fid, "mu_c", mu_C)
      CALL GETRVAL(fid, "SL0", SL0)
      CALL GETRVAL(fid, "SLmin", SLmin)
      CALL GETRVAL(fid, "SLmax", SLmax)
      CALL GETRVAL(fid, "f0", f0)
      CALL GETRVAL(fid, "fc1", fc1)
      CALL GETRVAL(fid, "fs1", fs1)
      CALL GETRVAL(fid, "fc2", fc2)
      CALL GETRVAL(fid, "fs2", fs2)
      CALL GETRVAL(fid, "fc3", fc3)
      CALL GETRVAL(fid, "fs3", fs3)

!     Electrophysiology model parameters
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "sV", sV)
      CALL GETRVAL(fid, "rho", rho)

      CLOSE(fid)

      RETURN
      END SUBROUTINE BO_READPARFF
!--------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE BO_INTEGFE(imyo, nX, X, Ts, Ti, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), RPAR(5)

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
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using forward Euler integration
      SUBROUTINE BO_ACTVSTRS_FE(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: epsX, f

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      CALL BO_ASTRS_GETF(X, epsX, Ta, f)
      Ta = Ta + (dt*f)

      RETURN
      END SUBROUTINE BO_ACTVSTRS_FE
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using RK4 integration
      SUBROUTINE BO_ACTVSTRS_RK(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: dt6, epsX, Ta_rk, f_rk(4)

      dt6  = dt / 6._RKIND
      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      Ta_rk = Ta
      CALL BO_ASTRS_GETF(X, epsX, Ta_rk, f_rk(1))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(1))
      CALL BO_ASTRS_GETF(X, epsX, Ta_rk, f_rk(2))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(2))
      CALL BO_ASTRS_GETF(X, epsX, Ta_rk, f_rk(3))

      Ta_rk = Ta + (dt*f_rk(3))
      CALL BO_ASTRS_GETF(X, epsX, Ta_rk, f_rk(4))

      Ta = Ta + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE BO_ACTVSTRS_RK
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using backward Euler integration
      SUBROUTINE BO_ACTVSTRS_BE(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: epsX

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      Ta = (Ta + (dt*epsX*K_T*(X-Vrest))) / (1._RKIND + epsX*dt)

      RETURN
      END SUBROUTINE BO_ACTVSTRS_BE
!-----------------------------------------------------------------------
      SUBROUTINE BO_ASTRS_GETF(X, eX, Ta, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, eX, Ta
      REAL(KIND=RKIND), INTENT(OUT) :: f

      f = eX*(K_T*(X-Vrest) - Ta)

      RETURN
      END SUBROUTINE BO_ASTRS_GETF
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using forward Euler
      SUBROUTINE BO_ACTVSTRN_FE(c_c, dt, I4f, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_c, dt, I4f
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      REAL(KIND=RKIND) :: f

      CALL BO_ASTRN_GETF(c_c, I4f, gf, f)
      gf = gf + (dt*f)

      RETURN
      END SUBROUTINE BO_ACTVSTRN_FE
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using Runge-Kutta
      SUBROUTINE BO_ACTVSTRN_RK(c_c, dt, I4f, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_c, dt, I4f
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      REAL(KIND=RKIND) :: dt6, gf_rk, f_rk(4)

      dt6  = dt / 6._RKIND

      gf_rk = gf
      CALL BO_ASTRN_GETF(c_c, I4f, gf_rk, f_rk(1))

      gf_rk = gf + (0.5_RKIND*dt*f_rk(1))
      CALL BO_ASTRN_GETF(c_c, I4f, gf_rk, f_rk(2))

      gf_rk = gf + (0.5_RKIND*dt*f_rk(2))
      CALL BO_ASTRN_GETF(c_c, I4f, gf_rk, f_rk(3))

      gf_rk = gf + (dt*f_rk(3))
      CALL BO_ASTRN_GETF(c_c, I4f, gf_rk, f_rk(4))

      gf = gf + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE BO_ACTVSTRN_RK
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using backward Euler
      SUBROUTINE BO_ACTVSTRN_BE(c_c, dt, I4f, gf, itMax, atol, rtol)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itMax
      REAL(KIND=RKIND), INTENT(IN) :: c_c, dt, I4f, atol, rtol
      REAL(KIND=RKIND), INTENT(OUT) :: gf

      INTEGER :: k
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: gk, f, rK, Jac

      CALL BO_ASTRN_GETF(c_c, I4f, gf, f)

      k  = 0
      gk = gf
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      DO
         k = k + 1

         CALL BO_ASTRN_GETF(c_c, I4f, gk, f)
         rK = gk - gf - (dt*f)

         l1   = k .GT. itMax
         l2   = ABS(rK) .LE. atol
         l3   = ABS(rK/gk) .LE. rtol
         IF (l1 .OR. l2 .OR. l3) EXIT

         CALL BO_ASTRN_GETJ(c_c, I4f, gk, Jac)
         Jac = 1._RKIND - dt*Jac

         gk  = gk - rK/(Jac + eps)
      END DO
      gf = gk

      IF (l1 .AND. .NOT.(l2) .AND. .NOT.(l3)) THEN
         WRITE(*,'(4X,A)') "Warning: Newton-Raphson failed to "//
     2      "converge (backward Euler)"
      END IF

      RETURN
      END SUBROUTINE BO_ACTVSTRN_BE
!-----------------------------------------------------------------------
      SUBROUTINE BO_ASTRN_GETF(c_c, I4f, gf, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_c, I4f, gf
      REAL(KIND=RKIND), INTENT(OUT) :: f

      INTEGER (KIND=IKIND) i, j
      REAL(KIND=RKIND) :: SL, RFL, mui, F_a, m1, ts

      SL  = SQRT(I4f) * SL0

      RFL = 0._RKIND
      IF (SL.GE.SLmin .AND. SL.LE.SLmax) THEN
         RFL = 0.5_RKIND*f0
     2       + fc1*COS(SL) + fc2*COS(2._RKIND*SL) + fc3*COS(3._RKIND*SL)
     3       + fs1*SIN(SL) + fs2*SIN(2._RKIND*SL) + fs3*SIN(3._RKIND*SL)
      END IF
      F_a  = alfa * RFL * (c_c - c_0)**2

      ts = 0._RKIND
      DO i=1, 5
         j  = (i+1)*(i+2)
         m1 = (-1._RKIND)**REAL(i,KIND=RKIND)
         ts = ts + m1*REAL(j,KIND=RKIND)*(gf**REAL(i,KIND=RKIND))
      END DO

      mui = 1._RKIND/(mu_c * c_c * c_c)
      f   = mui * (F_a + I4f*ts)

      RETURN
      END SUBROUTINE BO_ASTRN_GETF
!-----------------------------------------------------------------------
      SUBROUTINE BO_ASTRN_GETJ(c_c, I4f, gf, Jac)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_c, I4f, gf
      REAL(KIND=RKIND), INTENT(OUT) :: Jac

      INTEGER (KIND=IKIND) i, j
      REAL(KIND=RKIND) :: m1, mui, dts
c      REAL(KIND=RKIND) :: SL, dRFL, dFa

c      SL  = SQRT(I4f) * SL0

c      dRFL = 0._RKIND
c      IF (SL.GE.SLmin .AND. SL.LE.SLmax) THEN
c         dRFL = SL0 *  (fs1*COS(SL) - fc1*SIN(SL)
c     2      + 2._RKIND*(fs2*COS(2._RKIND*SL) - fc2*SIN(2._RKIND*SL))
c     3      + 3._RKIND*(fs3*COS(3._RKIND*SL) - fc3*SIN(3._RKIND*SL)))
c      END IF
c      dFa = alfa * dRFL * (c_c - c_0)**2

      dts = 0._RKIND
      DO i=1, 5
         j   = i*(i+1)*(i+2)
         m1  = (-1._RKIND)**REAL(i,KIND=RKIND)
         dts = dts + m1*REAL(j,KIND=RKIND)*(gf**REAL(i-1,KIND=RKIND))
      END DO

      mui = 1._RKIND/(mu_c * c_c * c_c)
      Jac = mui*I4f*dts

      RETURN
      END SUBROUTINE BO_ASTRN_GETJ
!-----------------------------------------------------------------------
      SUBROUTINE GETRVAL(fileId, skwrd, rVal)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fileId
      CHARACTER(LEN=*), INTENT(IN) :: skwrd
      REAL(KIND=RKIND), INTENT(INOUT) :: rVal

      INTEGER(KIND=IKIND) :: slen, i, ios
      CHARACTER(LEN=stdL) :: sline, scmd, sval

      REWIND(fileId)
      slen = LEN(TRIM(skwrd))
      DO
         READ(fileId,'(A)',END=001) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         DO i=1, slen
            IF (sline(i:i) .EQ. ':') EXIT
         END DO

         IF (i .GE. slen) THEN
            STOP "Error: inconsistent input file format"
         END IF

         scmd = sline(1:i-1)
         sval = sline(i+1:slen)
         sval = ADJUSTC(sval)

!        Remove any trailing comments
         slen = LEN(TRIM(sval))
         DO i=1, slen
            IF (sval(i:i) .EQ. '#') EXIT
         END DO
         sval = TRIM(ADJUSTC(sval(1:i-1)))

         IF (TRIM(skwrd) .EQ. TRIM(scmd)) THEN
            READ(sval,*,IOSTAT=ios) rval
            IF (ios .NE. 0) THEN
               WRITE(*,'(A)') " Error: while reading "//TRIM(skwrd)
               STOP
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

! 001  WRITE(*,'(A)') " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"
!      STOP

      END SUBROUTINE GETRVAL
!-----------------------------------------------------------------------
      SUBROUTINE GETRVEC(fileId, skwrd, nd, rVec)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fileId, nd
      CHARACTER(LEN=*), INTENT(IN) :: skwrd
      REAL(KIND=RKIND), INTENT(INOUT) :: rVec(nd)

      INTEGER(KIND=IKIND) :: slen, i, ios, nt
      CHARACTER(LEN=stdL) :: sline, scmd, sval
      CHARACTER(LEN=stdL), DIMENSION(1024) :: tokList

      REWIND(fileId)
      slen = LEN(TRIM(skwrd))
      DO
         READ(fileId,'(A)',END=001) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         DO i=1, slen
            IF (sline(i:i) .EQ. ':') EXIT
         END DO

         IF (i .GE. slen) THEN
            STOP "Error: inconsistent input file format"
         END IF

         scmd = sline(1:i-1)
         sval = sline(i+1:slen)
         sval = ADJUSTC(sval)

!        Remove any trailing comments
         slen = LEN(TRIM(sval))
         DO i=1, slen
            IF (sval(i:i) .EQ. '#') EXIT
         END DO
         sval = TRIM(ADJUSTC(sval(1:i-1)))

         IF (TRIM(skwrd) .EQ. TRIM(scmd)) THEN
            CALL PARSESTR(sval, tokList, nt)
            IF (nt .NE. nd) THEN
               DO i=1, nt
                  WRITE(*,'(I2,2X,A)') i, TRIM(tokList(i))
               END DO
               WRITE(*,'(A)') " Error: Unexpected token length "//
     2            TRIM(skwrd)
               STOP
            END IF

            DO i=1, nt
               READ(tokList(i),*,IOSTAT=ios) rvec(i)
               IF (ios .NE. 0) THEN
                  WRITE(*,'(A)') " Error: while reading "//TRIM(skwrd)
                  STOP
               END IF
            END DO
            EXIT
         END IF
      END DO

      RETURN

 001  WRITE(*,'(A)') " Error: EOF reached while finding command <"//
     2   TRIM(skwrd)//">"
      STOP

      END SUBROUTINE GETRVEC
!-----------------------------------------------------------------------
!     Removes any leading spaces or tabs
      PURE FUNCTION ADJUSTC(str)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str
      CHARACTER(LEN=LEN(str)) ADJUSTC

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(str)
         IF (str(i:i) .NE. " " .AND. str(i:i) .NE. "  ") EXIT
      END DO
      IF (i .GT. LEN(str)) THEN
         ADJUSTC = ""
      ELSE
         ADJUSTC = str(i:)
      END IF

      RETURN
      END FUNCTION ADJUSTC
!-----------------------------------------------------------------------
      SUBROUTINE PARSESTR(strng, toks, ntoks)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: strng
      CHARACTER(LEN=*), DIMENSION(1024), INTENT(OUT) :: toks
      INTEGER(KIND=IKIND), INTENT(OUT) :: ntoks

      INTEGER(KIND=IKIND) :: i, j, slen
      CHARACTER(LEN=stdL) :: dlmtr, token

      dlmtr = ''
      token = ''

      dlmtr = '< (,=")>'
      ntoks = 1
      slen  = LEN(TRIM(strng))

      ntoks = 0
      i = 0
      DO WHILE (i .LT. slen)
         i = i + 1
         IF (INDEX(dlmtr,strng(i:i)) .NE. 0) CYCLE
         DO j=i+1, slen
            IF (INDEX(dlmtr,strng(j:j)) .NE. 0) EXIT
         END DO
         IF (j .LE. slen) THEN
            ntoks = ntoks + 1
            toks(ntoks) = strng(i:j-1)
            i = j-1
         ELSE
            EXIT
         END IF
      END DO

      RETURN
      END SUBROUTINE PARSESTR
!--------------------------------------------------------------------
      END MODULE BOMOD
!####################################################################


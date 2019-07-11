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
!     This module defines data structures for ten Tusscher-Panfilov
!     epicardial cellular activation model for cardiac electrophysiology
!
!-----------------------------------------------------------------------

      MODULE TTPMOD
      IMPLICIT NONE

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE TTP_INIT(imyo, nX, nG, X, Xg, X0, Xg0)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imyo, nX, nG
      REAL(KINd=8), INTENT(OUT) :: X(nX), Xg(nG)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: X0(nX), Xg0(nG)

      SELECT CASE (imyo)
      CASE (1) ! EPI !
!        Initialize state variables
         X(1)   = -85.23D0      ! V      (units: mV)
         X(2)   =  136.89D0     ! K_i    (units: mM)
         X(3)   =  8.6040D0     ! Na_i   (units: mM)
         X(4)   =  1.26D-4      ! Ca_i   (units: mM)
         X(5)   =  3.6D-4       ! Ca_ss  (units: mM)
         X(6)   =  3.64D0       ! Ca_sr  (units: mM)
         X(7)   =  9.073D-1     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  6.21D-3      ! x_r1   (dimensionless)
         Xg(2)  =  4.712D-1     ! x_r2   (dimensionless)
         Xg(3)  =  9.5D-3       ! x_s    (dimensionless)
         Xg(4)  =  1.72D-3      ! m      (dimensionless)
         Xg(5)  =  7.444D-1     ! h      (dimensionless)
         Xg(6)  =  7.045D-1     ! j      (dimensionless)
         Xg(7)  =  3.373D-5     ! d      (dimensionless)
         Xg(8)  =  7.888D-1     ! f      (dimensionless)
         Xg(9)  =  9.755D-1     ! f_2    (dimensionless)
         Xg(10) =  9.953D-1     ! f_cass (dimensionless)
         Xg(11) =  9.99998D-1   ! s      (dimensionless)
         Xg(12) =  2.42D-8      ! r      (dimensionless)

      CASE (2) ! ENDO !
!        Initialize state variables
         X(1)   = -86.709D0     ! V      (units: mV)
         X(2)   =  138.4D0      ! K_i    (units: mM)
         X(3)   =  10.355D0     ! Na_i   (units: mM)
         X(4)   =  1.30D-4      ! Ca_i   (units: mM)
         X(5)   =  3.6D-4       ! Ca_ss  (units: mM)
         X(6)   =  3.715D0      ! Ca_sr  (units: mM)
         X(7)   =  9.068D-1     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  4.48D-3      ! x_r1   (dimensionless)
         Xg(2)  =  4.76D-1      ! x_r2   (dimensionless)
         Xg(3)  =  8.7D-3       ! x_s    (dimensionless)
         Xg(4)  =  1.55D-3      ! m      (dimensionless)
         Xg(5)  =  7.573D-1     ! h      (dimensionless)
         Xg(6)  =  7.225D-1     ! j      (dimensionless)
         Xg(7)  =  3.164D-5     ! d      (dimensionless)
         Xg(8)  =  8.009D-1     ! f      (dimensionless)
         Xg(9)  =  9.778D-1     ! f_2    (dimensionless)
         Xg(10) =  9.953D-1     ! f_cass (dimensionless)
         Xg(11) =  3.212D-1     ! s      (dimensionless)
         Xg(12) =  2.235D-8     ! r      (dimensionless)

      CASE (3) ! MID-MYO !
!        Initialize state variables
         X(1)   = -85.423D0     ! V      (units: mV)
         X(2)   =  138.52D0     ! K_i    (units: mM)
         X(3)   =  10.132D0     ! Na_i   (units: mM)
         X(4)   =  1.53D-4      ! Ca_i   (units: mM)
         X(5)   =  4.2D-4       ! Ca_ss  (units: mM)
         X(6)   =  4.272D0      ! Ca_sr  (units: mM)
         X(7)   =  8.978D-1     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  1.65D-2      ! x_r1   (dimensionless)
         Xg(2)  =  4.730D-1     ! x_r2   (dimensionless)
         Xg(3)  =  1.74D-2      ! x_s    (dimensionless)
         Xg(4)  =  1.65D-3      ! m      (dimensionless)
         Xg(5)  =  7.490D-1     ! h      (dimensionless)
         Xg(6)  =  6.788D-1     ! j      (dimensionless)
         Xg(7)  =  3.288D-5     ! d      (dimensionless)
         Xg(8)  =  7.026D-1     ! f      (dimensionless)
         Xg(9)  =  9.526D-1     ! f_2    (dimensionless)
         Xg(10) =  9.942D-1     ! f_cass (dimensionless)
         Xg(11) =  9.99998D-1   ! s      (dimensionless)
         Xg(12) =  2.347D-8     ! r      (dimensionless)

      END SELECT

      IF (PRESENT(X0)) THEN
         IF (SIZE(X0,1) .NE. nX) THEN
            STOP "ERROR: inconsistent array size (TTP initialization)"
         END IF
         X(:)  = X0(:)
      END IF

      IF (PRESENT(Xg0)) THEN
         IF (SIZE(Xg0,1) .NE. nG) THEN
            STOP "ERROR: inconsistent array size (TTP initialization)"
         END IF
         Xg(:) = Xg0(:)
      END IF

      RETURN
      END SUBROUTINE TTP_INIT
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE TTP_INTEGFE(imyo, nX, nG, X, Xg, Ts, dt, Istim, Ksac,
     2   RPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imyo, nX, nG
      REAL(KIND=8), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(18)
      REAL(KIND=8), INTENT(IN) :: Ts, dt, Istim, Ksac

      REAL(KIND=8) :: f(nX)

!     Get time derivatives (RHS)
      CALL TTP_GETF(imyo, nX, nG, X, Xg, f, Istim, Ksac, RPAR)

!     Update gating variables
      CALL TTP_UPDATEG(imyo, dt, nX, nG, X, Xg)

!     Update state variables
      X = X + dt*f(:)

      RETURN
      END SUBROUTINE TTP_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE TTP_INTEGRK(imyo, nX, nG, X, Xg, Ts, dt, Istim, Ksac,
     2   RPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imyo, nX, nG
      REAL(KIND=8), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(18)
      REAL(KIND=8), INTENT(IN) :: Ts, dt, Istim, Ksac

      REAL(KIND=8) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6.0D0

!     RK4: 1st pass
      Xrk = X
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac,
     2   RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL TTP_UPDATEG(imyo, 0.5D0*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5D0*dt*frk(:,1)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac,
     2   RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5D0*dt*frk(:,2)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac,
     2   RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL TTP_UPDATEG(imyo, dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac,
     2   RPAR)

      X = X + dt6*(frk(:,1) + 2.0D0*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE TTP_INTEGRK
!-----------------------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE TTP_GETF(i, nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, nX, nG
      REAL(KIND=8), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=8), INTENT(OUT) :: dX(nX)
      REAL(KIND=8), INTENT(INOUT) :: RPAR(18)

      INCLUDE "PARAMS_TTP.f"

      REAL(KIND=8) :: RT, a, b, tau, sq5, e1, e2, e3, e4, n1, n2, d1,
     2   d2, d3, I_sac

!     Local copies of state variables
      V     = X(1)
      K_i   = X(2)
      Na_i  = X(3)
      Ca_i  = X(4)
      Ca_ss = X(5)
      Ca_sr = X(6)
      R_bar = X(7)

!     Local copies of gating variables
      xr1   = Xg(1)
      xr2   = Xg(2)
      xs    = Xg(3)
      m     = Xg(4)
      h     = Xg(5)
      j     = Xg(6)
      d     = Xg(7)
      f     = Xg(8)
      f2    = Xg(9)
      fcass = Xg(10)
      s     = Xg(11)
      r     = Xg(12)

!     Stretch-activated currents
      I_sac = K_sac * (Vrest - V)

!      Diff = 1.0D0 / (1.0D1 * rho * Cm * sV)
      RT   = Rc * Tc / Fc
      E_K  = RT * LOG(K_o/K_i)
      E_Na = RT * LOG(Na_o/Na_i)
      E_Ca = 0.5D0 * RT * LOG(Ca_o/Ca_i)
      E_Ks = RT * LOG( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) )

!     I_Na: Fast sodium current
      I_Na = G_Na * (m**3.0D0) * h * j * (V - E_Na)

!     I_to: transient outward current
      I_to = G_to(i) * r * s * (V - E_K)

!     I_K1: inward rectifier outward current
      e1   = EXP(0.06D0*(V - E_K - 200.0D0))
      e2   = EXP(2.0D-4*(V - E_K + 100.0D0))
      e3   = EXP(0.1D0*(V - E_K - 10.0D0))
      e4   = EXP(-0.5D0*(V - E_K))
      a    = 0.1D0/(1.0D0 + e1)
      b    = (3.0D0*e2 + e3) / (1.0D0 + e4)
      tau  = a / (a + b)
      sq5  = SQRT(K_o/5.4D0)
      I_K1 = G_K1 * sq5 * tau * (V - E_K)

!     I_Kr: rapid delayed rectifier current
      I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K)

!     I_Ks: slow delayed rectifier current
      I_Ks = G_Ks(i) * (xs**2.0D0) * (V - E_Ks)

!     I_CaL: L-type Ca current
      a     = 2.0D0*(V-15.0D0)/RT
      b     = 2.0D0*a*Fc * (0.25D0*Ca_ss*EXP(a) - Ca_o) / (EXP(a)-1.0D0)
      I_CaL = G_CaL * d * f * f2 * fcass * b

!     I_NaCa: Na-Ca exchanger current
      e1     = EXP(gamma*V/RT)
      e2     = EXP((gamma-1.0D0)*V/RT)
      n1     = e1*(Na_i**3.0D0)*Ca_o - e2*(Na_o**3.0D0)*Ca_i*alpha
      d1     = K_mNai**3.0D0 + Na_o**3.0D0
      d2     = K_mCa + Ca_o
      d3     = 1.0D0 + K_sat*e2
      I_NaCa = K_NaCa * n1 / (d1*d2*d3)

!     I_NaK: Na-K pump current
      e1    = EXP(-0.1D0*V/RT)
      e2    = EXP(-V/RT)
      n1    = P_NaK * K_o * Na_i
      d1    = K_o + K_mK
      d2    = Na_i + K_mNa
      d3    = 1.0D0 + 0.1245D0*e1 + 0.0353D0*e2
      I_NaK = n1 / (d1*d2*d3)

!     I_pCa: plateau Ca current
      I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i)

!     I_pK: plateau K current
      I_pK  = G_pK * (V-E_K) / (1.0D0 + EXP((25.0D0-V)/5.98D0))

!     I_bCa: background Ca current
      I_bCa = G_bCa * (V - E_Ca)

!     I_bNa: background Na current
      I_bNa = G_bNa * (V - E_Na)

!     I_leak: Sacroplasmic Reticulum Ca leak current
      I_leak = V_leak * (Ca_sr - Ca_i)

!     I_up: Sacroplasmic Reticulum Ca pump current
      I_up  = Vmax_up / (1.0D0 + (K_up/Ca_i)**2.0D0)

!     I_rel: Ca induced Ca current (CICR)
      k_casr = max_sr - ((max_sr-min_sr)/(1.0D0 + (EC/Ca_sr)**2.0D0) )
      k1     = k1p / k_casr
      O      = k1 * R_bar * (Ca_ss**2.0D0) / (k3 + k1*(Ca_ss**2.0D0))
      I_rel  = V_rel * O * (Ca_sr - Ca_ss)

!     I_xfer: diffusive Ca current between Ca subspae and cytoplasm
      I_xfer = V_xfer * (Ca_ss - Ca_i)

!-----------------------------------------------------------------------
!     Now compute time derivatives
!     dV/dt: rate of change of transmembrane voltage
      dX(1)  = -(I_Na + I_to + I_K1 + I_Kr + I_Ks + I_CaL + I_NaCa +
     2   I_NaK + I_pCa + I_pK + I_bCa + I_bNa  + I_stim) + I_sac

!     dK_i/dt
      dX(2)  = -(Cm/(V_c*Fc)) * (I_K1 + I_to + I_Kr + I_Ks + I_pK -
     2   2D0*I_NaK + I_stim)

!     dNa_i/dt
      dX(3)  = -(Cm/(V_c*Fc)) * (I_Na + I_bNa + 3.0D0*(I_NaK + I_NaCa))

!     dCa_i/dt
      n1     = (I_leak - I_up)*V_sr/V_c + I_xfer
      n2     = -(Cm/(V_c*Fc)) * (I_bCa + I_pCa - 2.0D0*I_Naca) /2.0D0
      d1     = 1.0D0 + K_bufc*Buf_c/(Ca_i + K_bufc)**2.0D0
      dX(4)  = (n1 + n2)/d1

!     dCa_ss: rate of change of Ca_ss
      n1     = (-I_CaL*Cm/(2.0D0*Fc) + I_rel*V_sr - V_c*I_xfer) / V_ss
      d1     = 1.0D0 + K_bufss*Buf_ss/(Ca_ss + K_bufss)**2.0D0
      dX(5)  = n1 / d1

!     dCa_sr: rate of change of Ca_sr
      n1     = I_up - I_leak - I_rel
      d1     = 1.0D0 + K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**2.0D0
      dX(6)  = n1 / d1

!     Rbar: ryanodine receptor
      k2     = k2p * k_casr
      dX(7)  = -k2*Ca_ss*R_bar + k4*(1.0D0 - R_bar)

!     Quantities to be written to file
      RPAR(3)  = I_Na
      RPAR(4)  = I_K1
      RPAR(5)  = I_to
      RPAR(6)  = I_Kr
      RPAR(7)  = I_Ks
      RPAR(8)  = I_CaL
      RPAR(9)  = I_NaCa
      RPAR(10) = I_NaK
      RPAR(11) = I_pCa
      RPAR(12) = I_pK
      RPAR(13) = I_bCa
      RPAR(14) = I_bNa
      RPAR(15) = I_leak
      RPAR(16) = I_up
      RPAR(17) = I_rel
      RPAR(18) = I_xfer

      RETURN
      END SUBROUTINE TTP_GETF
!-----------------------------------------------------------------------
!     Update all the gating variables
      SUBROUTINE TTP_UPDATEG(i, dt, n, nG, X, Xg)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, n, nG
      REAL(KIND=8), INTENT(IN) :: dt, X(n)
      REAL(KIND=8), INTENT(INOUT) :: Xg(nG)

      INCLUDE "PARAMS_TTP.f"

      REAL(KIND=8) :: a, b, c, tau

      V     = X(1)
      Ca_ss = X(5)

      xr1   = Xg(1)
      xr2   = Xg(2)
      xs    = Xg(3)
      m     = Xg(4)
      h     = Xg(5)
      j     = Xg(6)
      d     = Xg(7)
      f     = Xg(8)
      f2    = Xg(9)
      fcass = Xg(10)
      s     = Xg(11)
      r     = Xg(12)

!     xr1: activation gate for I_Kr
      xr1i   = 1.0D0/(1.0D0 + EXP(-(26.0D0+V)/7.0D0))
      a      = 450.0D0/(1.0D0 + EXP(-(45.0D0+V)/10.0D0))
      b      = 6.0D0/(1.0D0 + EXP((30.0D0+V)/11.5D0))
      tau    = a*b
      Xg(1)  = xr1i - (xr1i - xr1)*EXP(-dt/tau)

!     xr2: inactivation gate for I_Kr
      xr2i   = 1.0D0 /(1.0D0 + EXP((88.0D0+V)/24.0D0))
      a      = 3.0D0 /(1.0D0 + EXP(-(60.0D0+V)/20.0D0))
      b      = 1.12D0/(1.0D0 + EXP(-(60.0D0-V)/20.0D0))
      tau    = a*b
      Xg(2)  = xr2i - (xr2i - xr2)*EXP(-dt/tau)

!     xs: activation gate for I_Ks
      xsi    = 1.0D0/(1.0D0 + EXP(-(5.0D0+V)/14.0D0))
      a      = 1400.0D0/SQRT(1.0D0 + EXP((5.0D0-V)/6.0D0))
      b      = 1.0D0/(1D0 + EXP((V-35.0D0)/15.0D0))
      tau    = a*b + 80.0D0
      Xg(3)  = xsi - (xsi - xs)*EXP(-dt/tau)

!     m: activation gate for I_Na
      mi     = 1.0D0/( (1.0D0 + EXP(-(56.86D0+V)/9.03D0))**2.0D0 )
      a      = 1.0D0/(1.0D0 + EXP(-(60.0D0+V)/5.0D0))
      b      = 0.1D0/(1.0D0 + EXP((35.0D0+V)/5.0D0))
     2       + 0.1D0/(1.0D0 + EXP((V-50.0D0)/200.0D0))
      tau    = a*b
      Xg(4)  = mi - (mi - m)*EXP(-dt/tau)

!     h: fast inactivation gate for I_Na
      hi     = 1.0D0/( (1.0D0 + EXP((71.55D0+V)/7.43D0))**2.0D0 )
      IF (V .GE. -40.0D0) THEN
         a   = 0.0D0
         b   = 0.77D0/(0.13D0*(1.0D0 + EXP(-(10.66D0+V)/11.1D0)))
      ELSE
         a   = 5.7D-2*EXP(-(80.0D0+V)/6.8D0)
         b   = 2.7D0*EXP(0.079D0*V) + 310000.0D0*EXP(0.3485D0*V)
      END IF
      tau    = 1.0D0 / (a + b)
      Xg(5)  = hi - (hi - h)*EXP(-dt/tau)

!     j: slow inactivation gate for I_Na
      ji     = 1.0D0/( (1.0D0 + EXP((71.55D0+V)/7.43D0))**2.0D0 )
      IF (V .GE. -40.0D0) THEN
         a   = 0.0D0
         b   = 6D-1*EXP(5.7D-2*V)/(1.0D0 + EXP(-0.1D0*(V+32.0D0)))
      ELSE
         a   = -(25428.0D0*EXP(0.2444D0*V) + 6.948D-6*EXP(-0.04391D0*V))
     2       * (V+37.78D0) / (1.0D0 + EXP(0.311D0*(79.23D0+V)))
         b   = 0.02424D0*EXP(-0.01052D0*V) /
     2         (1.0D0 + EXP(-0.1378D0*(40.14D0+V)))
      END IF
      tau    = 1.0D0 / (a + b)
      Xg(6)  = ji - (ji - j)*EXP(-dt/tau)

!     d: activation gate for I_CaL
      di     = 1.0D0/(1.0D0 + EXP(-(8.0D0+V)/7.5D0))
      a      = 1.4D0/(1.0D0 + EXP(-(35.0D0+V)/13.0D0)) + 0.25D0
      b      = 1.4D0/(1.0D0 + EXP((5.0D0+V)/5.0D0))
      c      = 1.0D0/(1.0D0 + EXP((50.0D0-V)/20.0D0))
      tau    = a*b + c
      Xg(7)  = di - (di - d)*EXP(-dt/tau)

!     f: slow inactivation gate for I_CaL
      fi     = 1.0D0/(1.0D0 + EXP((20.0D0+V)/7.0D0))
      a      = 1102.5D0*EXP(-((V+27.0D0)**2.0D0)/225.0D0)
      b      = 200.0D0/(1.0D0 + EXP((13.0D0-V)/10.0D0))
      c      = 180.0D0/(1.0D0 + EXP((30.0D0+V)/10.0D0)) + 20.0D0
      tau    = a + b + c
      Xg(8)  = fi - (fi - f)*EXP(-dt/tau)

!     f2: fast inactivation gate for I_CaL
      f2i    = 0.67D0/(1.0D0 + EXP((35.0D0+V)/7.0D0)) + 0.33D0
      a      = 562.0D0*EXP(-((27.0D0+V)**2.0D0) /240.0D0)
      b      = 31.0D0/(1.0D0 + EXP((25.0D0-V)/10.0D0))
      c      = 80.0D0/(1.0D0 + EXP((30.0D0+V)/10.0D0))
      tau    = a + b + c
      Xg(9)  = f2i - (f2i - f2)*EXP(-dt/tau)

!     fCass: inactivation gate for I_CaL into subspace
      c      = 1.0D0/(1.0D0 + (Ca_ss/0.05D0)**2.0D0)
      fcassi = 0.6D0*c  + 0.4D0
      tau    = 80.0D0*c + 2.0D0
      Xg(10) = fcassi - (fcassi - fcass)*EXP(-dt/tau)

!     s: inactivation gate for I_to
      IF (i.EQ.1 .OR. i.EQ.3) THEN
         si  = 1.0D0/(1.0D0 + EXP((20.0D0+V)/5.0D0))
         tau = 85.0D0*EXP(-((V+45.0D0)**2.0D0) /320.0D0)
     2       + 5.0D0/(1.0D0+EXP((V-20.0D0)/5.0D0)) + 3.0D0
      ELSE IF (i .EQ. 2) THEN
         si  = 1.0D0/(1.0D0 + EXP((28.0D0+V)/5.0D0))
         tau = 1000.0D0*EXP(-((V+67.0D0)**2.0D0) /1000.0D0) + 8.0D0
      END IF
      Xg(11) = si - (si - s)*EXP(-dt/tau)

!     r: activation gate for I_to
      ri     = 1.0D0/(1.0D0 + EXP((20.0D0-V)/6.0D0))
      tau    = 9.5D0*EXP(-((V+40.0D0)**2.0D0) /1800.0D0) + 0.8D0
      Xg(12) = ri - (ri - r)*EXP(-dt/tau)

      RETURN
      END SUBROUTINE TTP_UPDATEG
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE TTP_INTEGCN2(imyo, nX, nG, Xn, Xg, Ts, dt, Istim,
     2   Ksac, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: imyo, nX, nG
      INTEGER, INTENT(INOUT) :: IPAR(2)
      REAL(KIND=8), INTENT(INOUT) :: Xn(nX), Xg(nG), RPAR(18)
      REAL(KIND=8), INTENT(IN) :: Ts, dt, Istim, Ksac

      REAL(KIND=8), PARAMETER :: eps = EPSILON(eps)

      INTEGER :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=8) :: t, atol, rtol, Xk(nX), fn(nX), fk(nX), rK(nX),
     2   Im(nX,nX), JAC(nX,nX), rmsA, rmsR

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      Im = MAT_ID(nX)
      CALL TTP_GETF(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1

         CALL TTP_GETF(imyo, nX, nG, Xk, Xg, fk, Istim, Ksac, RPAR)
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

         CALL TTP_GETJ(imyo, nX, nG, Xk, Xg, JAC, Ksac)
         JAC = Im - 0.5D0*dt*JAC
         JAC = MAT_INV(JAC, nX)
         rK  = MATMUL(JAC, rK)
         Xk  = Xk - rK
      END DO
      Xn(:) = Xk(:)
      CALL TTP_UPDATEG(imyo, dt, nX, nG, Xn, Xg)
      CALL TTP_GETF(imyo, nX, nG, Xn, Xg, fn, Istim, Ksac, RPAR)

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE TTP_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE TTP_GETJ(i, nX, nG, X, Xg, JAC, Ksac)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i, nX, nG
      REAL(KIND=8), INTENT(IN) :: X(nX), Xg(nG), Ksac
      REAL(KIND=8), INTENT(OUT) :: JAC(nX,nX)

      INCLUDE "PARAMS_TTP.f"

      REAL(KIND=8) :: RT, a, b, c, tau, sq5, e1, e2, e3, e4, n1, n2, d1,
     2   d2, d3

!     Local copies of state variables
      V     = X(1)
      K_i   = X(2)
      Na_i  = X(3)
      Ca_i  = X(4)
      Ca_ss = X(5)
      Ca_sr = X(6)
      R_bar = X(7)

!     Local copies of gating variables
      xr1   = Xg(1)
      xr2   = Xg(2)
      xs    = Xg(3)
      m     = Xg(4)
      h     = Xg(5)
      j     = Xg(6)
      d     = Xg(7)
      f     = Xg(8)
      f2    = Xg(9)
      fcass = Xg(10)
      s     = Xg(11)
      r     = Xg(12)

      RT    = Rc * Tc / Fc
      E_K   = RT * LOG(K_o/K_i)
      E_Na  = RT * LOG(Na_o/Na_i)
      E_Ca  = 0.5D0 * RT * LOG(Ca_o/Ca_i)
      E_Ks  = RT * LOG( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) )

      E_K_Ki   = -RT / K_i
      E_Na_Nai = -RT / Na_i
      E_Ca_Cai = -RT / Ca_i / 2.0D0
      E_Ks_Ki  = -RT / (K_i + p_KNa*Na_i)
      E_Ks_Nai = p_KNa * E_Ks_Ki

!     I_Na: Fast sodium current
      I_Na = G_Na * (m**3.0D0) * h * j * (V - E_Na)
      I_Na_V   = G_Na * (m**3.0D0) * h * j
      I_Na_Nai = I_Na_V * (-E_Na_Nai)

!     I_to: transient outward current
      I_to = G_to(i) * r * s * (V - E_K)
      I_to_V  = G_to(i) * r * s
      I_to_Ki = I_to_V * (-E_K_Ki)

!     I_K1: inward rectifier outward current
      e1   = EXP(0.06D0*(V - E_K - 200.0D0))
      e2   = EXP(2.0D-4*(V - E_K + 100.0D0))
      e3   = EXP(0.1D0*(V - E_K - 10.0D0))
      e4   = EXP(-0.5D0*(V - E_K))
      a    = 0.1D0/(1.0D0 + e1)
      b    = (3.0D0*e2 + e3) / (1.0D0 + e4)
      tau  = a / (a + b)
      sq5  = SQRT(K_o/5.4D0)
      n1   = -0.006D0*e1/(1.0D0 + e1)**2.0D0
      n2   = (6.0D-4*e2 + 0.1D0*e3 + 0.5D0*b*e4)/(1.0D0 + e4)
      n1   = (a + b)*n1 - a*(n1 + n2)
      d1   = (a + b)**2.0D0
      I_K1 = G_K1 * sq5 * tau * (V - E_K)
      I_K1_V  = G_K1 * sq5 * (tau + (V - E_K)*n1/d1)
      I_K1_Ki = I_K1_V * (-E_K_Ki)

!     I_Kr: rapid delayed rectifier current
      I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K)
      I_Kr_V   = G_Kr * sq5 * xr1 * xr2
      I_Kr_Ki  = I_Kr_V * (-E_K_Ki)

!     I_Ks: slow delayed rectifier current
      I_Ks = G_Ks(i) * (xs**2.0D0) * (V - E_Ks)
      I_Ks_V   = G_Ks(i) * (xs**2.0D0)
      I_Ks_Ki  = I_Ks_V * (-E_Ks_Ki)
      I_Ks_Nai = I_Ks_V * (-E_Ks_Nai)

!     I_CaL: L-type Ca current
      a     = 2.0D0*(V-15.0D0)/RT
      b     = (0.25D0*Ca_ss*EXP(a) - Ca_o) / (EXP(a)-1.0D0)
      c     = G_CaL * d * f * f2 * fcass * (2.0D0*a*Fc)
      n1    = (EXP(a)/RT) / (EXP(a)-1.0D0)
      I_CaL = c * b
      I_CaL_V     = I_CaL/(V-15.0D0) + n1*(c*0.5D0*Ca_ss - 2.0D0*I_CaL)
      I_CaL_Cass  = c * 0.25D0 * n1 * RT

!     I_NaCa: Na-Ca exchanger current
      e1     = EXP(gamma*V/RT)
      e2     = EXP((gamma-1.0D0)*V/RT)
      n1     = e1*(Na_i**3.0D0)*Ca_o - e2*(Na_o**3.0D0)*Ca_i*alpha
      d1     = K_mNai**3.0D0 + Na_o**3.0D0
      d2     = K_mCa + Ca_o
      d3     = 1.0D0 + K_sat*e2
      c      = 1.0D0/(d1*d2*d3)
      I_NaCa = K_NaCa * n1 * c

      n1     = K_NaCa * c * ( e1*(Na_i**3.0D0)*Ca_o*(gamma/RT) -
     2    e2*(Na_o**3.0D0)*Ca_i*alpha*((gamma-1.0D0)/RT) )
      n2     = I_NaCa*K_sat*((gamma-1.0D0)/RT)*e2/d3
      I_NaCa_V   = n1 - n2
      I_NaCa_Nai =  K_NaCa * e1 * (3.0D0*Na_i**2.0D0) * Ca_o * c
      I_NaCa_Cai = -K_NaCa * e2 * (Na_o**3.0D0) * alpha * c

!     I_NaK: Na-K pump current
      e1    = EXP(-0.1D0*V/RT)
      e2    = EXP(-V/RT)
      n1    = P_NaK * K_o * Na_i
      d1    = K_o + K_mK
      d2    = Na_i + K_mNa
      d3    = 1.0D0 + 0.1245D0*e1 + 0.0353D0*e2
      I_NaK = n1 / (d1*d2*d3)
      n1    = (0.01245D0*e1 + 0.0353D0*e2)/RT
      I_NaK_V   = I_NaK * n1 / d3
      I_NaK_Nai = I_NaK * K_mNa/(Na_i*d2)

!     I_pCa: plateau Ca current
      d1 = (K_pCa + Ca_i)
      I_pCa = G_pCa * Ca_i / d1
      I_pCa_Cai = G_pCa * K_pCa/(d1*d1)

!     I_pK: plateau K current
      e1 = EXP((25.0D0-V)/5.98D0)
      I_pK  = G_pK * (V-E_K) / (1.0D0 + e1)
      I_pK_V  = (G_pK + I_pK*e1/5.98D0) / (1.0D0+e1)
      I_pK_Ki = G_pK * (-E_K_Ki) / (1.0D0 + e1)

!     I_bCa: background Ca current
      I_bCa = G_bCa * (V - E_Ca)
      I_bCa_V = G_bCa
      I_bCa_Cai = G_bca * (-E_Ca_Cai)

!     I_bNa: background Na current
      I_bNa = G_bNa * (X(1) - E_Na)
      I_bNa_V = G_bNa
      I_bNa_Nai = G_bNa * (-E_Na_Nai)

!     I_leak: Sacroplasmic Reticulum Ca leak current
      I_leak = V_leak * (Ca_sr - Ca_i)
      I_leak_Cai  = -V_leak
      I_leak_Casr =  V_leak

!     I_up: Sacroplasmic Reticulum Ca pump current
      d1 = 1.0D0 + (K_up/Ca_i)**2.0D0
      I_up  = Vmax_up / d1
      I_up_Cai = (I_up / d1) * (2.0D0*K_up**2.0D0 / Ca_i**3.0D0)

!     I_rel: Ca induced Ca current (CICR)
      n1 = max_sr - min_sr
      d1 = 1.0D0 + (EC/Ca_sr)**2.0D0
      k_casr = max_sr - (n1/d1)
      k1     = k1p / k_casr
      n2     = Ca_ss*2.0D0
      d2     = k3 + k1*n2
      O      = k1 * R_bar * n2 / d2
      I_rel  = V_rel * O * (Ca_sr - Ca_ss)

      k_casr_sr = (n1 / (d1**2.0D0)) * (2.0D0*EC**2.0D0 / Ca_sr**3.0D0)
      k1_casr   = -k1p * k_casr_sr / (k_casr**2.0D0)
      O_Cass = 2.0D0 * k3 * O / (Ca_ss * d2)
      O_Casr = k1_casr * n2 * (R_bar - O) / d2
      O_Rbar = k1 * n2 / d2

      I_rel_Cass = V_rel * (O_Cass*(Ca_sr - Ca_ss) - O)
      I_rel_Casr = V_rel * (O_Casr*(Ca_sr - Ca_ss) + O)
      I_rel_Rbar = V_rel * O_Rbar *(Ca_sr - Ca_ss)

!     I_xfer: diffusive Ca current between Ca subspae and cytoplasm
      I_xfer = V_xfer * (Ca_ss - Ca_i)
      I_xfer_Cai  = -V_xfer
      I_xfer_Cass =  V_xfer

!-----------------------------------------------------------------------
!     Compute Jacobian matrix
      JAC = 0D0
      c   = Cm/(V_c*Fc)

!     V
      JAC(1,1)  = -(I_Na_V + I_to_V + I_K1_V + I_Kr_V + I_Ks_V + I_CaL_V
     2   + I_NaCa_V + I_NaK_V + I_pK_V + I_bCa_V + I_bNa_V + Ksac)
      JAC(1,2)  = -(I_to_Ki + I_K1_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki)
      JAC(1,3)  = -(I_Na_Nai + I_Ks_Nai + I_NaCa_Nai + I_NaK_Nai
     2   + I_bNa_Nai)
      JAC(1,4)  = -(I_NaCa_Cai + I_pCa_Cai + I_bCa_Cai)
      JAC(1,5) = -I_CaL_Cass

!     K_i
      JAC(2,1)  = -c * (I_K1_V + I_to_V + I_Kr_V + I_Ks_V + I_pK_V
     2   - 2D0*I_NaK_V )
      JAC(2,2)  = -c * (I_K1_Ki + I_to_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki)
      JAC(2,3)  = -c * (I_Ks_Nai - 2D0*I_NaK_Nai)

!     Na_i
      JAC(3,1)  = -c * (I_Na_V + I_bNa_V + 3.0D0*(I_NaK_V + I_NaCa_V))
      JAC(3,3)  = -c * (I_Na_Nai + I_bNa_Nai + 3.0D0*(I_NaK_Nai +
     2   I_NaCa_Nai))
      JAC(3,4)  = -c * (3.0D0*I_NaCa_Cai)

!     Ca_i
      n1 = (I_leak - I_up)*V_sr/V_c + I_xfer
     2   - c*(I_bCa + I_pCa - 2.0D0*I_NaCa)/2.0D0
      n2 = (I_leak_Cai - I_up_Cai)*V_sr/V_c + I_xfer_Cai
     2   - c*(I_bCa_Cai + I_pCa_Cai - 2.0D0*I_NaCa_Cai)/2.0D0
      d1 = 1.0D0 + K_bufc*Buf_c/(Ca_i + K_bufc)**2.0D0
      d2 = 2.0D0*K_bufc*Buf_c/(Ca_i + K_bufc)**3.0D0
      JAC(4,1)  = -c * (I_bCa_V - 2.0D0*I_NaCa_V) / 2.0D0 / d1
      JAC(4,3)  = c * I_NaCa_Nai / d1
      JAC(4,4)  = (n2 + n1*d2/d1) / d1
      JAC(4,5) = I_xfer_Cass / d1
      JAC(4,6) = (I_leak_Casr*V_sr/V_c) / d1

!     Ca_ss
      a  = Cm/(2.0D0*Fc*V_ss)
      b  = V_sr / V_ss
      c  = V_c / V_ss
      n1 = -a*I_CaL + b*I_rel - c*I_xfer
      n2 = -a*I_CaL_Cass + b*I_rel_Cass - c*I_xfer_Cass
      d1 = 1.0D0 + K_bufss*Buf_ss/(Ca_ss + K_bufss)**2.0D0
      d2 = 2.0D0*K_bufss*Buf_ss/(Ca_ss + K_bufss)**3.0D0
      JAC(5,1)  = -a * I_CaL_V / d1
      JAC(5,4)  = -c * I_xfer_Cai / d1
      JAC(5,5) = (n2 + n1*d2/d1) / d1
      JAC(5,6) = b * I_rel_Casr / d1
      JAC(5,7) = b * I_rel_Rbar / d1

!     Ca_sr
      n1 = I_up - I_leak - I_rel
      n2 = -(I_leak_Casr + I_rel_Casr)
      d1 = 1.0D0 + K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**2.0D0
      d2 = 2.0D0*K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**3.0D0
      JAC(6,4) = (I_up_Cai - I_leak_Cai) / d1
      JAC(6,5) = -I_rel_Cass / d1
      JAC(6,6) = (n2 + n1*d2/d1) / d1
      JAC(6,7) = -I_rel_Rbar / d1

!     Rbar: ryanodine receptor
      k2 = k2p * k_casr
      JAC(7,5) = -k2 * R_bar
      JAC(7,6) = -(k2p * k_casr_sr) * Ca_ss * R_bar
      JAC(7,7) = -(k2*Ca_ss + k4)

      RETURN
      END SUBROUTINE TTP_GETJ
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model
      SUBROUTINE TTP_ACTVSTRS(c_Ca, dt, Tact, epsX)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: c_Ca, dt
      REAL(KIND=8), INTENT(OUT) :: epsX
      REAL(KIND=8), INTENT(INOUT) :: Tact

      INCLUDE "PARAMS_TTP.f"

      REAL(KIND=8) :: nr

      epsX = EXP(-EXP(-xi_T*(c_Ca - Ca_crit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX
      nr   = Tact + epsX*dt*eta_T*(c_Ca - Ca_rest)
      Tact = nr / (1.0D0 + epsX*dt)

      RETURN
      END SUBROUTINE TTP_ACTVSTRS
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration
      SUBROUTINE TTP_ACTVSTRN(c_Ca, I4f, dt, gf)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: c_Ca, I4f, dt
      REAL(KIND=8), INTENT(INOUT) :: gf

      INCLUDE "PARAMS_TTP.f"

      REAL(KIND=8) :: SL, Fa, rtmp

!     fiber length
      SL = I4f * SL0

!     Sacromere force-length relationship
      IF (SL.GE.SLmin .AND. SL.LE.SLmax) THEN
         SL = 0.5D0*f0 + fc1*COS(SL) + fs1*SIN(SL) +
     2      fc2*COS(2.0D0*SL) + fs2*SIN(2.0D0*SL)  +
     3      fc3*COS(3.0D0*SL) + fs3*SIN(3.0D0*SL)
      ELSE
         SL = 0.0D0
      END IF

!     Active force
      Fa   = alFa * (c_Ca-c_Ca0)*(c_Ca-c_Ca0) * SL

      rtmp = 2.0D0*I4f*(1.0D0/(1.0D0+gf)**3.0D0 - 1.0D0)
      gf = gf + dt*(Fa + rtmp)/(mu_Ca * c_Ca * c_Ca)

      RETURN
      END SUBROUTINE TTP_ACTVSTRN
!-----------------------------------------------------------------------
      END MODULE TTPMOD
!#######################################################################

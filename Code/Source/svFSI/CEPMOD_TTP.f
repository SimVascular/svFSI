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
!     This module defines data structures for tenTusscher-Panfilov
!     cellular activation model for cardiac electrophysiology. Both
!     active stress and active strain models are used for excitation-
!     contraction coupling.
!
!     Reference for tenTusscher-Panfilov electrophysiology model:
!        ten Tusscher, et al. (2004). A model for human ventricular
!        tissue. American Journal of Physiology - Heart and Circulatory
!        Physiology, 286(4), H1573–H1589.
!
!        ten Tusscher, K. H. W. J., & Panfilov, A. V. (2006).
!        Alternans and spiral breakup in a human ventricular tissue
!        model. American Journal of Physiology. Heart and Circulatory
!        Physiology, 291(3), H1088–H1100.
!        https://doi.org/10.1152/ajpheart.00109.2006
!
!     References for active stress/active strain models:
!        Garcia-blanco, E., et al. (2019). A new computational framework
!        for electro-activation in cardiac mechanics. Computer Methods
!        in Applied Mechanics and Engineering.
!        https://doi.org/10.1016/j.cma.2019.01.042
!
!        Wong, J., Goketepe, S., & Kuhl, E. (2013). Computational
!        modeling of chemo-electro-mechanical coupling: A novel implicit
!        monolithic finite element approach. International Journal for
!        Numerical Methods in Biomedical Engineering, 29, 1104–1133.
!        https://doi.org/10.1002/cnm.2565
!
!        Simone Rossi, et al. (2012). Orthotropic active strain models
!        for the numerical simulation of cardiac biomechanics.
!        International Journal for Numerical Methods in Biomedical
!        Engineering, 28, 761–788. https://doi.org/10.1002/cnm.2473
!
!-----------------------------------------------------------------------

      MODULE TTPMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INCLUDE "PARAMS_TTP.f"

      PUBLIC :: TTP_INIT
      PUBLIC :: TTP_READPARFF
      PUBLIC :: TTP_INTEGFE
      PUBLIC :: TTP_INTEGRK
      PUBLIC :: TTP_INTEGCN2
      PUBLIC :: TTP_ACTVSTRS_FE
      PUBLIC :: TTP_ACTVSTRS_RK
      PUBLIC :: TTP_ACTVSTRS_BE
      PUBLIC :: TTP_ACTVSTRN_FE
      PUBLIC :: TTP_ACTVSTRN_RK
      PUBLIC :: TTP_ACTVSTRN_BE

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE TTP_INIT(imyo, nX, nG, X, Xg)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

      SELECT CASE (imyo)
      CASE (1) ! EPI !
!        Initialize state variables
         X(1)   = -85.23_RKIND      ! V      (units: mV)
         X(2)   =  136.89_RKIND     ! K_i    (units: mM)
         X(3)   =  8.6040_RKIND     ! Na_i   (units: mM)
         X(4)   =  1.26E-4_RKIND    ! Ca_i   (units: mM)
         X(5)   =  3.6E-4_RKIND     ! Ca_ss  (units: mM)
         X(6)   =  3.64_RKIND       ! Ca_sr  (units: mM)
         X(7)   =  0.9073_RKIND     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  6.21E-3_RKIND    ! x_r1   (dimensionless)
         Xg(2)  =  0.4712_RKIND     ! x_r2   (dimensionless)
         Xg(3)  =  9.5E-3_RKIND     ! x_s    (dimensionless)
         Xg(4)  =  1.72E-3_RKIND    ! m      (dimensionless)
         Xg(5)  =  0.7444_RKIND     ! h      (dimensionless)
         Xg(6)  =  0.7045_RKIND     ! j      (dimensionless)
         Xg(7)  =  3.373E-5_RKIND   ! d      (dimensionless)
         Xg(8)  =  0.7888_RKIND     ! f      (dimensionless)
         Xg(9)  =  0.9755_RKIND     ! f_2    (dimensionless)
         Xg(10) =  0.9953_RKIND     ! f_cass (dimensionless)
         Xg(11) =  0.999998_RKIND   ! s      (dimensionless)
         Xg(12) =  2.42E-8_RKIND    ! r      (dimensionless)

      CASE (2) ! ENDO !
!        Initialize state variables
         X(1)   = -86.709_RKIND     ! V      (units: mV)
         X(2)   =  138.4_RKIND      ! K_i    (units: mM)
         X(3)   =  10.355_RKIND     ! Na_i   (units: mM)
         X(4)   =  1.3E-4_RKIND    ! Ca_i   (units: mM)
         X(5)   =  3.6E-4_RKIND     ! Ca_ss  (units: mM)
         X(6)   =  3.715_RKIND      ! Ca_sr  (units: mM)
         X(7)   =  0.9068_RKIND     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  4.48E-3_RKIND    ! x_r1   (dimensionless)
         Xg(2)  =  0.476_RKIND      ! x_r2   (dimensionless)
         Xg(3)  =  8.7E-3_RKIND     ! x_s    (dimensionless)
         Xg(4)  =  1.55E-3_RKIND    ! m      (dimensionless)
         Xg(5)  =  0.7573_RKIND     ! h      (dimensionless)
         Xg(6)  =  0.7225_RKIND     ! j      (dimensionless)
         Xg(7)  =  3.164E-5_RKIND   ! d      (dimensionless)
         Xg(8)  =  0.8009_RKIND     ! f      (dimensionless)
         Xg(9)  =  0.9778_RKIND     ! f_2    (dimensionless)
         Xg(10) =  0.9953_RKIND     ! f_cass (dimensionless)
         Xg(11) =  0.3212_RKIND     ! s      (dimensionless)
         Xg(12) =  2.235E-8_RKIND   ! r      (dimensionless)

      CASE (3) ! MID-MYO !
!        Initialize state variables
         X(1)   = -85.423_RKIND     ! V      (units: mV)
         X(2)   =  138.52_RKIND     ! K_i    (units: mM)
         X(3)   =  10.132_RKIND     ! Na_i   (units: mM)
         X(4)   =  1.53E-4_RKIND    ! Ca_i   (units: mM)
         X(5)   =  4.2E-4_RKIND     ! Ca_ss  (units: mM)
         X(6)   =  4.272_RKIND      ! Ca_sr  (units: mM)
         X(7)   =  0.8978_RKIND     ! R'     (dimensionless)

!        Initialize gating variables
         Xg(1)  =  1.65E-2_RKIND    ! x_r1   (dimensionless)
         Xg(2)  =  0.4730_RKIND     ! x_r2   (dimensionless)
         Xg(3)  =  1.74E-2_RKIND    ! x_s    (dimensionless)
         Xg(4)  =  1.65E-3_RKIND    ! m      (dimensionless)
         Xg(5)  =  0.7490_RKIND     ! h      (dimensionless)
         Xg(6)  =  0.6788_RKIND     ! j      (dimensionless)
         Xg(7)  =  3.288E-5_RKIND   ! d      (dimensionless)
         Xg(8)  =  0.7026_RKIND     ! f      (dimensionless)
         Xg(9)  =  0.9526_RKIND     ! f_2    (dimensionless)
         Xg(10) =  0.9942_RKIND     ! f_cass (dimensionless)
         Xg(11) =  0.999998_RKIND   ! s      (dimensionless)
         Xg(12) =  2.347E-8_RKIND   ! r      (dimensionless)

      END SELECT

      RETURN
      END SUBROUTINE TTP_INIT
!-----------------------------------------------------------------------
      SUBROUTINE TTP_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))
      CALL GETRVAL(fid, "Rc", Rc)
      CALL GETRVAL(fid, "Tc", Tc)
      CALL GETRVAL(fid, "Fc", Fc)
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "sV", sV)
      CALL GETRVAL(fid, "rho", rho)
      CALL GETRVAL(fid, "V_c", V_c)
      CALL GETRVAL(fid, "V_sr", V_sr)
      CALL GETRVAL(fid, "V_ss", V_ss)
      CALL GETRVAL(fid, "K_o", K_o)
      CALL GETRVAL(fid, "Na_o", Na_o)
      CALL GETRVAL(fid, "Ca_o", Ca_o)
      CALL GETRVAL(fid, "G_Na", G_Na)
      CALL GETRVAL(fid, "G_K1", G_K1)
      CALL GETRVEC(fid, "G_to", 3, G_to)
      CALL GETRVAL(fid, "G_Kr", G_Kr)
      CALL GETRVEC(fid, "G_Ks", 3, G_Ks)
      CALL GETRVAL(fid, "p_KNa", p_KNa)
      CALL GETRVAL(fid, "G_CaL", G_CaL)
      CALL GETRVAL(fid, "K_NaCa", K_NaCa)
      CALL GETRVAL(fid, "gamma", gamma)
      CALL GETRVAL(fid, "K_mCa", K_mCa)
      CALL GETRVAL(fid, "K_mNai", K_mNai)
      CALL GETRVAL(fid, "K_sat", K_sat)
      CALL GETRVAL(fid, "alpha", alpha)
      CALL GETRVAL(fid, "p_NaK", p_NaK)
      CALL GETRVAL(fid, "K_mK", K_mK)
      CALL GETRVAL(fid, "K_mNa", K_mNa)
      CALL GETRVAL(fid, "G_pK", G_pK)
      CALL GETRVAL(fid, "G_pCa", G_pCa)
      CALL GETRVAL(fid, "K_pCa", K_pCa)
      CALL GETRVAL(fid, "G_bNa", G_bNa)
      CALL GETRVAL(fid, "G_bCa", G_bCa)
      CALL GETRVAL(fid, "Vmax_up", Vmax_up)
      CALL GETRVAL(fid, "K_up", K_up)
      CALL GETRVAL(fid, "V_rel", V_rel)
      CALL GETRVAL(fid, "k1p", k1p)
      CALL GETRVAL(fid, "k2p", k2p)
      CALL GETRVAL(fid, "k3", k3)
      CALL GETRVAL(fid, "k4", k4)
      CALL GETRVAL(fid, "EC", EC)
      CALL GETRVAL(fid, "max_sr", max_sr)
      CALL GETRVAL(fid, "min_sr", min_sr)
      CALL GETRVAL(fid, "V_leak", V_leak)
      CALL GETRVAL(fid, "V_xfer", V_xfer)
      CALL GETRVAL(fid, "Buf_c", Buf_c)
      CALL GETRVAL(fid, "K_bufc", K_bufc)
      CALL GETRVAL(fid, "Buf_sr", Buf_sr)
      CALL GETRVAL(fid, "K_bufsr", K_bufsr)
      CALL GETRVAL(fid, "Buf_ss", Buf_ss)
      CALL GETRVAL(fid, "K_bufss", K_bufss)
      CALL GETRVAL(fid, "Vrest", Vrest)

!     Electromechanics coupling parameters: active stress model
      CALL GETRVAL(fid, "Ca_rest", Ca_rest)
      CALL GETRVAL(fid, "Ca_crit", Ca_crit)
      CALL GETRVAL(fid, "K_T", K_T)
      CALL GETRVAL(fid, "eps_0", eps_0)
      CALL GETRVAL(fid, "eps_i", eps_i)
      CALL GETRVAL(fid, "xi_T", xi_T)

!     Electromechanics coupling parameters: active strain model
      CALL GETRVAL(fid, "alfa", alfa)
      CALL GETRVAL(fid, "c_Ca0", c_Ca0)
      CALL GETRVAL(fid, "mu_Ca", mu_Ca)
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

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

      CLOSE(fid)

      RETURN
      END SUBROUTINE TTP_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE TTP_INTEGFE(imyo, nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(18)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: f(nX)

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
      SUBROUTINE TTP_INTEGRK(imyo, nX, nG, X, Xg, dt, Istim, Ksac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX, nG
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), RPAR(18)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      REAL(KIND=RKIND) :: dt6, Xrk(nX), Xgr(nG), frk(nX,4)

      dt6 = dt/6._RKIND

!     RK4: 1st pass
      Xrk = X
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xg, frk(:,1), Istim, Ksac,
     2   RPAR)

!     Update gating variables by half-dt
      Xgr = Xg
      CALL TTP_UPDATEG(imyo, 0.5_RKIND*dt, nX, nG, X, Xgr)

!     RK4: 2nd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,1)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,2), Istim, Ksac,
     2   RPAR)

!     RK4: 3rd pass
      Xrk = X + 0.5_RKIND*dt*frk(:,2)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,3), Istim, Ksac,
     2   RPAR)

!     Update gating variables by full-dt
      Xgr = Xg
      CALL TTP_UPDATEG(imyo, dt, nX, nG, X, Xgr)

!     RK4: 4th pass
      Xrk = X + dt*frk(:,3)
      CALL TTP_GETF(imyo, nX, nG, Xrk, Xgr, frk(:,4), Istim, Ksac,
     2   RPAR)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))
      Xg = Xgr

      RETURN
      END SUBROUTINE TTP_INTEGRK
!-----------------------------------------------------------------------
!     Compute currents and time derivatives of state variables
      SUBROUTINE TTP_GETF(i, nX, nG, X, Xg, dX, I_stim, K_sac, RPAR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: i, nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), I_stim, K_sac
      REAL(KIND=RKIND), INTENT(OUT) :: dX(nX)
      REAL(KIND=RKIND), INTENT(INOUT) :: RPAR(18)

      REAL(KIND=RKIND) :: RT, a, b, tau, sq5, e1, e2, e3, e4, n1, n2,
     2   d1, d2, d3, I_sac

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

!      Diff = 1._RKIND / (1.0D1 * rho * Cm * sV)
      RT   = Rc * Tc / Fc
      E_K  = RT * LOG(K_o/K_i)
      E_Na = RT * LOG(Na_o/Na_i)
      E_Ca = 0.5_RKIND * RT * LOG(Ca_o/Ca_i)
      E_Ks = RT * LOG( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) )

!     I_Na: Fast sodium current
      I_Na = G_Na * (m**3._RKIND) * h * j * (V - E_Na)

!     I_to: transient outward current
      I_to = G_to(i) * r * s * (V - E_K)

!     I_K1: inward rectifier outward current
      e1   = EXP(0.06_RKIND*(V - E_K - 200._RKIND))
      e2   = EXP(2.E-4_RKIND*(V - E_K + 100._RKIND))
      e3   = EXP(0.1_RKIND*(V - E_K - 10._RKIND))
      e4   = EXP(-0.5_RKIND*(V - E_K))
      a    = 0.1_RKIND/(1._RKIND + e1)
      b    = (3._RKIND*e2 + e3) / (1._RKIND + e4)
      tau  = a / (a + b)
      sq5  = SQRT(K_o/5.4_RKIND)
      I_K1 = G_K1 * sq5 * tau * (V - E_K)

!     I_Kr: rapid delayed rectifier current
      I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K)

!     I_Ks: slow delayed rectifier current
      I_Ks = G_Ks(i) * (xs**2._RKIND) * (V - E_Ks)

!     I_CaL: L-type Ca current
      a     = 2._RKIND*(V-15._RKIND)/RT
      b     = 2._RKIND*a*Fc * (0.25_RKIND*Ca_ss*EXP(a) - Ca_o) /
     2   (EXP(a)-1._RKIND)
      I_CaL = G_CaL * d * f * f2 * fcass * b

!     I_NaCa: Na-Ca exchanger current
      e1     = EXP(gamma*V/RT)
      e2     = EXP((gamma-1._RKIND)*V/RT)
      n1     = e1*(Na_i**3._RKIND)*Ca_o - e2*(Na_o**3._RKIND)*Ca_i*alpha
      d1     = K_mNai**3._RKIND + Na_o**3._RKIND
      d2     = K_mCa + Ca_o
      d3     = 1._RKIND + K_sat*e2
      I_NaCa = K_NaCa * n1 / (d1*d2*d3)

!     I_NaK: Na-K pump current
      e1    = EXP(-0.1_RKIND*V/RT)
      e2    = EXP(-V/RT)
      n1    = P_NaK * K_o * Na_i
      d1    = K_o + K_mK
      d2    = Na_i + K_mNa
      d3    = 1._RKIND + 0.1245_RKIND*e1 + 0.0353_RKIND*e2
      I_NaK = n1 / (d1*d2*d3)

!     I_pCa: plateau Ca current
      I_pCa = G_pCa * Ca_i / (K_pCa + Ca_i)

!     I_pK: plateau K current
      I_pK  = G_pK * (V-E_K) /
     2   (1._RKIND + EXP((25._RKIND-V)/5.98_RKIND))

!     I_bCa: background Ca current
      I_bCa = G_bCa * (V - E_Ca)

!     I_bNa: background Na current
      I_bNa = G_bNa * (V - E_Na)

!     I_leak: Sacroplasmic Reticulum Ca leak current
      I_leak = V_leak * (Ca_sr - Ca_i)

!     I_up: Sacroplasmic Reticulum Ca pump current
      I_up  = Vmax_up / (1._RKIND + (K_up/Ca_i)**2._RKIND)

!     I_rel: Ca induced Ca current (CICR)
      k_casr = max_sr - ((max_sr-min_sr)/
     2   (1._RKIND + (EC/Ca_sr)**2._RKIND) )
      k1     = k1p / k_casr
      O      = k1 * R_bar * (Ca_ss**2._RKIND) /
     2   (k3 + k1*(Ca_ss**2._RKIND))
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
     2   2._RKIND*I_NaK + I_stim)

!     dNa_i/dt
      dX(3)  = -(Cm/(V_c*Fc)) * (I_Na + I_bNa +
     2   3._RKIND*(I_NaK + I_NaCa))

!     dCa_i/dt
      n1     = (I_leak - I_up)*V_sr/V_c + I_xfer
      n2     = -(Cm/(V_c*Fc)) * (I_bCa + I_pCa - 2._RKIND*I_Naca)
     2   / 2._RKIND
      d1     = 1._RKIND + K_bufc*Buf_c/(Ca_i + K_bufc)**2._RKIND
      dX(4)  = (n1 + n2)/d1

!     dCa_ss: rate of change of Ca_ss
      n1     = (-I_CaL*Cm/(2._RKIND*Fc) + I_rel*V_sr - V_c*I_xfer)/V_ss
      d1     = 1._RKIND + K_bufss*Buf_ss/(Ca_ss + K_bufss)**2._RKIND
      dX(5)  = n1 / d1

!     dCa_sr: rate of change of Ca_sr
      n1     = I_up - I_leak - I_rel
      d1     = 1._RKIND + K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**2._RKIND
      dX(6)  = n1 / d1

!     Rbar: ryanodine receptor
      k2     = k2p * k_casr
      dX(7)  = -k2*Ca_ss*R_bar + k4*(1._RKIND - R_bar)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: i, n, nG
      REAL(KIND=RKIND), INTENT(IN) :: dt, X(n)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xg(nG)

      REAL(KIND=RKIND) :: a, b, c, tau

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
      xr1i   = 1._RKIND/(1._RKIND + EXP(-(26._RKIND+V)/7._RKIND))
      a      = 450._RKIND/(1._RKIND + EXP(-(45._RKIND+V)/10._RKIND))
      b      = 6._RKIND/(1._RKIND + EXP((30._RKIND+V)/11.5_RKIND))
      tau    = a*b
      Xg(1)  = xr1i - (xr1i - xr1)*EXP(-dt/tau)

!     xr2: inactivation gate for I_Kr
      xr2i   = 1._RKIND /(1._RKIND + EXP((88._RKIND+V)/24._RKIND))
      a      = 3._RKIND /(1._RKIND + EXP(-(60._RKIND+V)/20._RKIND))
      b      = 1.12_RKIND/(1._RKIND + EXP(-(60._RKIND-V)/20._RKIND))
      tau    = a*b
      Xg(2)  = xr2i - (xr2i - xr2)*EXP(-dt/tau)

!     xs: activation gate for I_Ks
      xsi    = 1._RKIND/(1._RKIND + EXP(-(5._RKIND+V)/14._RKIND))
      a      = 1400._RKIND/SQRT(1._RKIND + EXP((5._RKIND-V)/6._RKIND))
      b      = 1._RKIND/(1._RKIND + EXP((V-35._RKIND)/15._RKIND))
      tau    = a*b + 80._RKIND
      Xg(3)  = xsi - (xsi - xs)*EXP(-dt/tau)

!     m: activation gate for I_Na
      mi     = 1._RKIND/( (1._RKIND +
     2   EXP(-(56.86_RKIND+V)/9.03_RKIND))**2._RKIND )
      a      = 1._RKIND/(1._RKIND + EXP(-(60._RKIND+V)/5._RKIND))
      b      = 0.1_RKIND/(1._RKIND + EXP((35._RKIND+V)/5._RKIND))
     2       + 0.1_RKIND/(1._RKIND + EXP((V-50._RKIND)/200._RKIND))
      tau    = a*b
      Xg(4)  = mi - (mi - m)*EXP(-dt/tau)

!     h: fast inactivation gate for I_Na
      hi     = 1._RKIND/( (1._RKIND
     2       + EXP((71.55_RKIND+V)/7.43_RKIND))**2._RKIND )
      IF (V .GE. -40._RKIND) THEN
         a   = 0._RKIND
         b   = 0.77_RKIND/(0.13_RKIND*(1._RKIND
     2       + EXP(-(10.66_RKIND+V)/11.1_RKIND)))
      ELSE
         a   = 5.7E-2_RKIND*EXP(-(80._RKIND+V)/6.8_RKIND)
         b   = 2.7_RKIND*EXP(0.079_RKIND*V)
     2       + 310000._RKIND*EXP(0.3485_RKIND*V)
      END IF
      tau    = 1._RKIND / (a + b)
      Xg(5)  = hi - (hi - h)*EXP(-dt/tau)

!     j: slow inactivation gate for I_Na
      ji     = 1._RKIND/( (1._RKIND
     2       + EXP((71.55_RKIND+V)/7.43_RKIND))**2._RKIND )
      IF (V .GE. -40._RKIND) THEN
         a   = 0._RKIND
         b   = 0.6_RKIND*EXP(5.7E-2_RKIND*V)
     2       / (1._RKIND + EXP(-0.1_RKIND*(V+32._RKIND)))
      ELSE
         a   = -(25428._RKIND*EXP(0.2444_RKIND*V)
     2       + 6.948E-6_RKIND*EXP(-0.04391_RKIND*V)) * (V+37.78_RKIND)
     3       / (1._RKIND + EXP(0.311_RKIND*(79.23_RKIND+V)))
         b   = 0.02424_RKIND*EXP(-0.01052_RKIND*V) /
     2         (1._RKIND + EXP(-0.1378_RKIND*(40.14_RKIND+V)))
      END IF
      tau    = 1._RKIND / (a + b)
      Xg(6)  = ji - (ji - j)*EXP(-dt/tau)

!     d: activation gate for I_CaL
      di     = 1._RKIND/(1._RKIND + EXP(-(8._RKIND+V)/7.5_RKIND))
      a      = 1.4_RKIND/(1._RKIND + EXP(-(35._RKIND+V)/13._RKIND))
     2       + 0.25_RKIND
      b      = 1.4_RKIND/(1._RKIND + EXP((5._RKIND+V)/5._RKIND))
      c      = 1._RKIND/(1._RKIND + EXP((50._RKIND-V)/20._RKIND))
      tau    = a*b + c
      Xg(7)  = di - (di - d)*EXP(-dt/tau)

!     f: slow inactivation gate for I_CaL
      fi     = 1._RKIND/(1._RKIND + EXP((20._RKIND+V)/7._RKIND))
      a      = 1102.5_RKIND*EXP(-((V+27._RKIND)**2._RKIND)/225._RKIND)
      b      = 200._RKIND/(1._RKIND + EXP((13._RKIND-V)/10._RKIND))
      c      = 180._RKIND/(1._RKIND + EXP((30._RKIND+V)/10._RKIND))
     2       + 20._RKIND
      tau    = a + b + c
c!     for spiral wave breakup
c      IF (V .GT. 0._RKIND) tau = tau*2._RKIND
      Xg(8)  = fi - (fi - f)*EXP(-dt/tau)

!     f2: fast inactivation gate for I_CaL
      f2i    = 0.67_RKIND/(1._RKIND + EXP((35._RKIND+V)/7._RKIND))
     2       + 0.33_RKIND
      a      = 562._RKIND*EXP(-((27._RKIND+V)**2._RKIND) /240._RKIND)
      b      = 31._RKIND/(1._RKIND + EXP((25._RKIND-V)/10._RKIND))
      c      = 80._RKIND/(1._RKIND + EXP((30._RKIND+V)/10._RKIND))
      tau    = a + b + c
      Xg(9)  = f2i - (f2i - f2)*EXP(-dt/tau)

!     fCass: inactivation gate for I_CaL into subspace
      c      = 1._RKIND/(1._RKIND + (Ca_ss/0.05_RKIND)**2._RKIND)
      fcassi = 0.6_RKIND*c  + 0.4_RKIND
      tau    = 80._RKIND*c + 2._RKIND
      Xg(10) = fcassi - (fcassi - fcass)*EXP(-dt/tau)

!     s: inactivation gate for I_to
      IF (i.EQ.1 .OR. i.EQ.3) THEN
         si  = 1._RKIND/(1._RKIND + EXP((20._RKIND+V)/5._RKIND))
         tau = 85._RKIND*EXP(-((V+45._RKIND)**2._RKIND) /320._RKIND)
     2       + 5._RKIND/(1._RKIND+EXP((V-20._RKIND)/5._RKIND))
     3       + 3._RKIND
      ELSE IF (i .EQ. 2) THEN
         si  = 1._RKIND/(1._RKIND + EXP((28._RKIND+V)/5._RKIND))
         tau = 1000._RKIND*EXP(-((V+67._RKIND)**2._RKIND) /1000._RKIND)
     2       + 8._RKIND
      END IF
      Xg(11) = si - (si - s)*EXP(-dt/tau)

!     r: activation gate for I_to
      ri     = 1._RKIND/(1._RKIND + EXP((20._RKIND-V)/6._RKIND))
      tau    = 9.5_RKIND*EXP(-((V+40._RKIND)**2._RKIND) /1800._RKIND)
     2       + 0.8_RKIND
      Xg(12) = ri - (ri - r)*EXP(-dt/tau)

      RETURN
      END SUBROUTINE TTP_UPDATEG
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE TTP_INTEGCN2(imyo, nX, nG, Xn, Xg, dt, Istim, Ksac,
     2   IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: imyo, nX, nG
      INTEGER(KIND=IKIND), INTENT(INOUT) :: IPAR(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xn(nX), Xg(nG), RPAR(18)
      REAL(KIND=RKIND), INTENT(IN) :: dt, Istim, Ksac

      INTEGER(KIND=IKIND) :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: atol, rtol, Xk(nX), fn(nX), fk(nX), rK(nX),
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
      DO
         k = k + 1

         CALL TTP_GETF(imyo, nX, nG, Xk, Xg, fk, Istim, Ksac, RPAR)
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

         CALL TTP_GETJ(imyo, nX, nG, Xk, Xg, JAC, Ksac)
         JAC = Im - 0.5_RKIND*dt*JAC
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
      INTEGER(KIND=IKIND), INTENT(IN) :: i, nX, nG
      REAL(KIND=RKIND), INTENT(IN) :: X(nX), Xg(nG), Ksac
      REAL(KIND=RKIND), INTENT(OUT) :: JAC(nX,nX)

      REAL(KIND=RKIND) :: RT, a, b, c, tau, sq5, e1, e2, e3, e4, n1, n2,
     2   d1, d2, d3

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
      E_Ca  = 0.5_RKIND * RT * LOG(Ca_o/Ca_i)
      E_Ks  = RT * LOG( (K_o + p_KNa*Na_o)/(K_i + p_KNa*Na_i) )

      E_K_Ki   = -RT / K_i
      E_Na_Nai = -RT / Na_i
      E_Ca_Cai = -RT / Ca_i / 2._RKIND
      E_Ks_Ki  = -RT / (K_i + p_KNa*Na_i)
      E_Ks_Nai = p_KNa * E_Ks_Ki

!     I_Na: Fast sodium current
      I_Na = G_Na * (m**3._RKIND) * h * j * (V - E_Na)
      I_Na_V   = G_Na * (m**3._RKIND) * h * j
      I_Na_Nai = I_Na_V * (-E_Na_Nai)

!     I_to: transient outward current
      I_to = G_to(i) * r * s * (V - E_K)
      I_to_V  = G_to(i) * r * s
      I_to_Ki = I_to_V * (-E_K_Ki)

!     I_K1: inward rectifier outward current
      e1   = EXP(0.06_RKIND*(V - E_K - 200._RKIND))
      e2   = EXP(2.E-4_RKIND*(V - E_K + 100._RKIND))
      e3   = EXP(0.1_RKIND*(V - E_K - 10._RKIND))
      e4   = EXP(-0.5_RKIND*(V - E_K))
      a    = 0.1_RKIND/(1._RKIND + e1)
      b    = (3._RKIND*e2 + e3) / (1._RKIND + e4)
      tau  = a / (a + b)
      sq5  = SQRT(K_o/5.4_RKIND)
      n1   = -6.E-3_RKIND*e1/(1._RKIND + e1)**2._RKIND
      n2   = (6.E-4_RKIND*e2 + 0.1_RKIND*e3 + 0.5_RKIND*b*e4)
     2     / (1._RKIND + e4)
      n1   = (a + b)*n1 - a*(n1 + n2)
      d1   = (a + b)**2._RKIND
      I_K1 = G_K1 * sq5 * tau * (V - E_K)
      I_K1_V  = G_K1 * sq5 * (tau + (V - E_K)*n1/d1)
      I_K1_Ki = I_K1_V * (-E_K_Ki)

!     I_Kr: rapid delayed rectifier current
      I_Kr = G_Kr * sq5 * xr1 * xr2 * (V - E_K)
      I_Kr_V   = G_Kr * sq5 * xr1 * xr2
      I_Kr_Ki  = I_Kr_V * (-E_K_Ki)

!     I_Ks: slow delayed rectifier current
      I_Ks = G_Ks(i) * (xs**2._RKIND) * (V - E_Ks)
      I_Ks_V   = G_Ks(i) * (xs**2._RKIND)
      I_Ks_Ki  = I_Ks_V * (-E_Ks_Ki)
      I_Ks_Nai = I_Ks_V * (-E_Ks_Nai)

!     I_CaL: L-type Ca current
      a     = 2._RKIND*(V-15._RKIND)/RT
      b     = (0.25_RKIND*Ca_ss*EXP(a) - Ca_o) / (EXP(a)-1._RKIND)
      c     = G_CaL * d * f * f2 * fcass * (2._RKIND*a*Fc)
      n1    = (EXP(a)/RT) / (EXP(a)-1._RKIND)
      I_CaL = c * b
      I_CaL_V = I_CaL/(V-15._RKIND) + n1*(c*0.5_RKIND*Ca_ss
     2        - 2._RKIND*I_CaL)
      I_CaL_Cass  = c * 0.25_RKIND * n1 * RT

!     I_NaCa: Na-Ca exchanger current
      e1     = EXP(gamma*V/RT)
      e2     = EXP((gamma-1._RKIND)*V/RT)
      n1     = e1*(Na_i**3._RKIND)*Ca_o - e2*(Na_o**3._RKIND)*Ca_i*alpha
      d1     = K_mNai**3._RKIND + Na_o**3._RKIND
      d2     = K_mCa + Ca_o
      d3     = 1._RKIND + K_sat*e2
      c      = 1._RKIND/(d1*d2*d3)
      I_NaCa = K_NaCa * n1 * c

      n1     = K_NaCa * c * ( e1*(Na_i**3._RKIND)*Ca_o*(gamma/RT) -
     2    e2*(Na_o**3._RKIND)*Ca_i*alpha*((gamma-1._RKIND)/RT) )
      n2     = I_NaCa*K_sat*((gamma-1._RKIND)/RT)*e2/d3
      I_NaCa_V   = n1 - n2
      I_NaCa_Nai =  K_NaCa * e1 * (3._RKIND*Na_i**2._RKIND) * Ca_o * c
      I_NaCa_Cai = -K_NaCa * e2 * (Na_o**3._RKIND) * alpha * c

!     I_NaK: Na-K pump current
      e1    = EXP(-0.1_RKIND*V/RT)
      e2    = EXP(-V/RT)
      n1    = P_NaK * K_o * Na_i
      d1    = K_o + K_mK
      d2    = Na_i + K_mNa
      d3    = 1._RKIND + 0.1245_RKIND*e1 + 0.0353_RKIND*e2
      I_NaK = n1 / (d1*d2*d3)
      n1    = (0.01245_RKIND*e1 + 0.0353_RKIND*e2)/RT
      I_NaK_V   = I_NaK * n1 / d3
      I_NaK_Nai = I_NaK * K_mNa/(Na_i*d2)

!     I_pCa: plateau Ca current
      d1 = (K_pCa + Ca_i)
      I_pCa = G_pCa * Ca_i / d1
      I_pCa_Cai = G_pCa * K_pCa/(d1*d1)

!     I_pK: plateau K current
      e1 = EXP((25._RKIND-V)/5.98_RKIND)
      I_pK  = G_pK * (V-E_K) / (1._RKIND + e1)
      I_pK_V  = (G_pK + I_pK*e1/5.98_RKIND) / (1._RKIND+e1)
      I_pK_Ki = G_pK * (-E_K_Ki) / (1._RKIND + e1)

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
      d1 = 1._RKIND + (K_up/Ca_i)**2._RKIND
      I_up  = Vmax_up / d1
      I_up_Cai = (I_up / d1) * (2._RKIND*K_up**2._RKIND /Ca_i**3._RKIND)

!     I_rel: Ca induced Ca current (CICR)
      n1 = max_sr - min_sr
      d1 = 1._RKIND + (EC/Ca_sr)**2._RKIND
      k_casr = max_sr - (n1/d1)
      k1     = k1p / k_casr
      n2     = Ca_ss*2._RKIND
      d2     = k3 + k1*n2
      O      = k1 * R_bar * n2 / d2
      I_rel  = V_rel * O * (Ca_sr - Ca_ss)

      k_casr_sr = (n1 / (d1**2._RKIND)) * (2._RKIND*EC**2._RKIND
     2          / Ca_sr**3._RKIND)
      k1_casr   = -k1p * k_casr_sr / (k_casr**2._RKIND)
      O_Cass = 2._RKIND * k3 * O / (Ca_ss * d2)
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
      JAC = 0._RKIND
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
     2   - 2._RKIND*I_NaK_V )
      JAC(2,2)  = -c * (I_K1_Ki + I_to_Ki + I_Kr_Ki + I_Ks_Ki + I_pK_Ki)
      JAC(2,3)  = -c * (I_Ks_Nai - 2._RKIND*I_NaK_Nai)

!     Na_i
      JAC(3,1)  = -c * (I_Na_V + I_bNa_V +
     2   3._RKIND*(I_NaK_V + I_NaCa_V))
      JAC(3,3)  = -c * (I_Na_Nai + I_bNa_Nai + 3._RKIND*(I_NaK_Nai +
     2   I_NaCa_Nai))
      JAC(3,4)  = -c * (3._RKIND*I_NaCa_Cai)

!     Ca_i
      n1 = (I_leak - I_up)*V_sr/V_c + I_xfer
     2   - 0.5_RKIND*c*(I_bCa + I_pCa - 2._RKIND*I_NaCa)
      n2 = (I_leak_Cai - I_up_Cai)*V_sr/V_c + I_xfer_Cai
     2   - 0.5_RKIND*c*(I_bCa_Cai + I_pCa_Cai - 2._RKIND*I_NaCa_Cai)
      d1 = 1._RKIND + K_bufc*Buf_c/(Ca_i + K_bufc)**2._RKIND
      d2 = 2._RKIND*K_bufc*Buf_c/(Ca_i + K_bufc)**3._RKIND
      JAC(4,1)  = -c * (I_bCa_V - 2._RKIND*I_NaCa_V) / 2._RKIND / d1
      JAC(4,3)  = c * I_NaCa_Nai / d1
      JAC(4,4)  = (n2 + n1*d2/d1) / d1
      JAC(4,5) = I_xfer_Cass / d1
      JAC(4,6) = (I_leak_Casr*V_sr/V_c) / d1

!     Ca_ss
      a  = Cm/(2._RKIND*Fc*V_ss)
      b  = V_sr / V_ss
      c  = V_c / V_ss
      n1 = -a*I_CaL + b*I_rel - c*I_xfer
      n2 = -a*I_CaL_Cass + b*I_rel_Cass - c*I_xfer_Cass
      d1 = 1._RKIND + K_bufss*Buf_ss/(Ca_ss + K_bufss)**2._RKIND
      d2 = 2._RKIND*K_bufss*Buf_ss/(Ca_ss + K_bufss)**3._RKIND
      JAC(5,1)  = -a * I_CaL_V / d1
      JAC(5,4)  = -c * I_xfer_Cai / d1
      JAC(5,5) = (n2 + n1*d2/d1) / d1
      JAC(5,6) = b * I_rel_Casr / d1
      JAC(5,7) = b * I_rel_Rbar / d1

!     Ca_sr
      n1 = I_up - I_leak - I_rel
      n2 = -(I_leak_Casr + I_rel_Casr)
      d1 = 1._RKIND + K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**2._RKIND
      d2 = 2._RKIND*K_bufsr*Buf_sr/(Ca_sr + K_bufsr)**3._RKIND
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
!     stress model using forward Euler integration
      SUBROUTINE TTP_ACTVSTRS_FE(c_Ca, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: eps_Ca, f

      eps_Ca = EXP(-EXP(-xi_T*(c_Ca - Ca_crit)))
      eps_Ca = eps_0 + (eps_i - eps_0)*eps_Ca

      CALL TTP_ASTRS_GETF(c_Ca, eps_Ca, Ta, f)
      Ta = Ta + (dt*f)

      RETURN
      END SUBROUTINE TTP_ACTVSTRS_FE
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using RK4 integration
      SUBROUTINE TTP_ACTVSTRS_RK(c_Ca, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: dt6, eps_Ca, Ta_rk, f_rk(4)

      dt6  = dt / 6._RKIND
      eps_Ca = EXP(-EXP(-xi_T*(c_Ca - Ca_crit)))
      eps_Ca = eps_0 + (eps_i - eps_0)*eps_Ca

      Ta_rk = Ta
      CALL TTP_ASTRS_GETF(c_Ca, eps_Ca, Ta_rk, f_rk(1))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(1))
      CALL TTP_ASTRS_GETF(c_Ca, eps_Ca, Ta_rk, f_rk(2))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(2))
      CALL TTP_ASTRS_GETF(c_Ca, eps_Ca, Ta_rk, f_rk(3))

      Ta_rk = Ta + (dt*f_rk(3))
      CALL TTP_ASTRS_GETF(c_Ca, eps_Ca, Ta_rk, f_rk(4))

      Ta = Ta + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE TTP_ACTVSTRS_RK
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using backward Euler integration
      SUBROUTINE TTP_ACTVSTRS_BE(c_Ca, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: eps_Ca

      eps_Ca = EXP(-EXP(-xi_T*(c_Ca - Ca_crit)))
      eps_Ca = eps_0 + (eps_i - eps_0)*eps_Ca

      Ta = (Ta + (dt*eps_Ca*K_T*(c_Ca-Ca_rest)))/(1._RKIND + eps_Ca*dt)

      RETURN
      END SUBROUTINE TTP_ACTVSTRS_BE
!-----------------------------------------------------------------------
      SUBROUTINE TTP_ASTRS_GETF(c_Ca, e_Ca, Ta, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, e_Ca, Ta
      REAL(KIND=RKIND), INTENT(OUT) :: f

      f = e_Ca*(K_T*(c_Ca-Ca_rest) - Ta)

      RETURN
      END SUBROUTINE TTP_ASTRS_GETF
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using forward Euler
      SUBROUTINE TTP_ACTVSTRN_FE(c_Ca, dt, I4f, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt, I4f
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      REAL(KIND=RKIND) :: f

      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf, f)
      gf = gf + (dt*f)

      RETURN
      END SUBROUTINE TTP_ACTVSTRN_FE
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using Runge-Kutta
      SUBROUTINE TTP_ACTVSTRN_RK(c_Ca, dt, I4f, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt, I4f
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      REAL(KIND=RKIND) :: dt6, gf_rk, f_rk(4)

      dt6  = dt / 6._RKIND

      gf_rk = gf
      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf_rk, f_rk(1))

      gf_rk = gf + (0.5_RKIND*dt*f_rk(1))
      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf_rk, f_rk(2))

      gf_rk = gf + (0.5_RKIND*dt*f_rk(2))
      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf_rk, f_rk(3))

      gf_rk = gf + (dt*f_rk(3))
      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf_rk, f_rk(4))

      gf = gf + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE TTP_ACTVSTRN_RK
!-----------------------------------------------------------------------
!     Compute macroscopic fiber strain based on sacromere force-length
!     relationship and calcium concentration using backward Euler
      SUBROUTINE TTP_ACTVSTRN_BE(c_Ca, dt, I4f, gf, itMax, atol, rtol)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: itMax
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, dt, I4f, atol, rtol
      REAL(KIND=RKIND), INTENT(OUT) :: gf

      INTEGER :: k
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: gk, f, rK, Jac

      CALL TTP_ASTRN_GETF(c_Ca, I4f, gf, f)

      k  = 0
      gk = gf
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      DO
         k = k + 1

         CALL TTP_ASTRN_GETF(c_Ca, I4f, gk, f)
         rK = gk - gf - (dt*f)

         l1   = k .GT. itMax
         l2   = ABS(rK) .LE. atol
         l3   = ABS(rK/gk) .LE. rtol
         IF (l1 .OR. l2 .OR. l3) EXIT

         CALL TTP_ASTRN_GETJ(c_Ca, I4f, gk, Jac)
         Jac = 1._RKIND - dt*Jac

         gk  = gk - rK/(Jac + eps)
      END DO
      gf = gk

      IF (l1 .AND. .NOT.(l2) .AND. .NOT.(l3)) THEN
         WRITE(*,'(4X,A)') "Warning: Newton-Raphson failed to "//
     2      "converge (backward Euler)"
      END IF

      RETURN
      END SUBROUTINE TTP_ACTVSTRN_BE
!-----------------------------------------------------------------------
      SUBROUTINE TTP_ASTRN_GETF(c_Ca, I4f, gf, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, I4f, gf
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
      F_a  = alfa * RFL * (c_Ca - c_Ca0)**2

      ts = 0._RKIND
      DO i=1, 5
         j  = (i+1)*(i+2)
         m1 = (-1._RKIND)**REAL(i,KIND=RKIND)
         ts = ts + m1*REAL(j,KIND=RKIND)*(gf**REAL(i,KIND=RKIND))
      END DO

      mui = 1._RKIND/(mu_Ca * c_Ca * c_Ca)
      f   = mui * (F_a + I4f*ts)

      RETURN
      END SUBROUTINE TTP_ASTRN_GETF
!-----------------------------------------------------------------------
      SUBROUTINE TTP_ASTRN_GETJ(c_Ca, I4f, gf, Jac)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: c_Ca, I4f, gf
      REAL(KIND=RKIND), INTENT(OUT) :: Jac

      INTEGER (KIND=IKIND) i, j
      REAL(KIND=RKIND) :: mui, m1, dts
c      REAL(KIND=RKIND) :: SL, dRFL, dFa

c      SL   = SQRT(I4f) * SL0
c      dRFL = 0._RKIND
c      IF (SL.GE.SLmin .AND. SL.LE.SLmax) THEN
c         dRFL = SL0 *  (fs1*COS(SL) - fc1*SIN(SL)
c     2      + 2._RKIND*(fs2*COS(2._RKIND*SL) - fc2*SIN(2._RKIND*SL))
c     3      + 3._RKIND*(fs3*COS(3._RKIND*SL) - fc3*SIN(3._RKIND*SL)))
c      END IF
c      dFa = alfa * dRFL * (c_Ca - c_Ca0)**2

      dts = 0._RKIND
      DO i=1, 5
         j   = i*(i+1)*(i+2)
         m1  = (-1._RKIND)**REAL(i,KIND=RKIND)
         dts = dts + m1*REAL(j,KIND=RKIND)*(gf**REAL(i-1,KIND=RKIND))
      END DO

      mui  = 1._RKIND/(mu_Ca * c_Ca * c_Ca)
      Jac  = mui*I4f*dts

      RETURN
      END SUBROUTINE TTP_ASTRN_GETJ
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

 001  RETURN

! 001  WRITE(*,'(A)') " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"
!      STOP

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
!-----------------------------------------------------------------------
      END MODULE TTPMOD
!#######################################################################

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
!     Parameters for TenTusscher-Panfilov Ventricular Myocyte Model.
!     Parameters are chosen based on below references:
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
!--------------------------------------------------------------------

!     Default model parameters
!     R: Gas constant
      REAL(KIND=RKIND) :: Rc = 8314.472_RKIND      ! units: J/mol/K
!     T: Temperature
      REAL(KIND=RKIND) :: Tc = 310._RKIND          ! units: K
!     F: Faraday constant
      REAL(KIND=RKIND) :: Fc = 96485.3415_RKIND    ! units: C/mmol
!     Cm: Cell capacitance per unit surface area
      REAL(KIND=RKIND) :: Cm = 0.185_RKIND         ! units: uF/cm^{2}
!     sV: Surface to volume ratio
      REAL(KIND=RKIND) :: sV = 0.2_RKIND           ! units: um^{-1}
!     rho: Cellular resistivity
      REAL(KIND=RKIND) :: rho = 162._RKIND         ! units: \Omega-cm
!     V_c: Cytoplasmic volume
      REAL(KIND=RKIND) :: V_c = 16.404E-3_RKIND    ! units: um^{3}
!     V_sr: Sacroplasmic reticulum volume
      REAL(KIND=RKIND) :: V_sr = 1.094E-3_RKIND    ! units: um^{3}
!     V_ss: Subspace volume
      REAL(KIND=RKIND) :: V_ss = 5.468E-5_RKIND    ! units: um^{3}
!     K_o: Extracellular K concentration
      REAL(KIND=RKIND) :: K_o = 5.4_RKIND          ! units: mM
!     Na_o: Extracellular Na concentration
      REAL(KIND=RKIND) :: Na_o = 140._RKIND        ! units: mM
!     Ca_o: Extracellular Ca concentration
      REAL(KIND=RKIND) :: Ca_o = 2._RKIND          ! units: mM
!     G_Na: Maximal I_Na conductance
      REAL(KIND=RKIND) :: G_Na = 14.838_RKIND      ! units: nS/pF
!     G_K1: Maximal I_K1 conductance
      REAL(KIND=RKIND) :: G_K1 = 5.405_RKIND       ! units: nS/pF
!     G_to: Maximal epicardial I_to conductance, units: nS/pF
      REAL(KIND=RKIND) :: G_to(3) =
     2   (/0.294_RKIND, 0.073_RKIND, 0.294_RKIND/)
!     G_Kr: Maximal I_Kr conductance
      REAL(KIND=RKIND) :: G_Kr = 0.153_RKIND       ! units: nS/pF
c!     G_Kr for spiral wave breakup
c      REAL(KIND=RKIND) :: G_Kr = 0.172_RKIND      ! units: nS/pF
!     G_Ks: Maximal epicardial I_Ks conductance, units: nS/pF
      REAL(KIND=RKIND) :: G_Ks(3) =
     2   (/0.392_RKIND, 0.392_RKIND, 0.098_RKIND/)
c!     G_Ks for spiral wave breakup (epi)
c      REAL(KIND=RKIND) :: G_Ks(3) =
c     2   (/0.441_RKIND, 0.392_RKIND, 0.098_RKIND/)
!     p_KNa: Relative I_Ks permeability to Na
      REAL(KIND=RKIND) :: p_KNa = 3.E-2_RKIND      ! dimensionless
!     G_CaL: Maximal I_CaL conductance
      REAL(KIND=RKIND) :: G_CaL = 3.98E-5_RKIND    ! units: cm^{3}/uF/ms
!     K_NaCa: Maximal I_NaCa
      REAL(KIND=RKIND) :: K_NaCa = 1000._RKIND     ! units: pA/pF
!     gamma: Voltage dependent parameter of I_NaCa
      REAL(KIND=RKIND) :: gamma = 0.35_RKIND       ! dimensionless
!     K_mCa: Ca_i half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: K_mCa = 1.38_RKIND       ! units: mM
!     K_mNai: Na_i half-saturation constant for I_NaCa
      REAL(KIND=RKIND) :: K_mNai = 87.5_RKIND      ! units: mM
!     K_sat: Saturation factor for I_NaCa
      REAL(KIND=RKIND) :: K_sat = 0.1_RKIND        ! dimensionless
!     alpha: Factor enhancing outward nature of I_NaCa
      REAL(KIND=RKIND) :: alpha = 2.5_RKIND        ! dimensionless
!     p_NaK: Maximal I_NaK
      REAL(KIND=RKIND) :: p_NaK = 2.724_RKIND      ! units: pA/pF
!     K_mK: K_o half-saturation constant of I_NaK
      REAL(KIND=RKIND) :: K_mK = 1._RKIND          ! units: mM
!     K_mNa: Na_i half-saturation constant of I_NaK
      REAL(KIND=RKIND) :: K_mNa = 40._RKIND        ! units: mM
!     G_pK: Maximal I_pK conductance
      REAL(KIND=RKIND) :: G_pK = 1.46E-2_RKIND     ! units: nS/pF
c!     G_pK for spiral wave breakup
c      REAL(KIND=RKIND) :: G_pK = 2.19E-3_RKIND     ! units: nS/pF
!     G_pCa: Maximal I_pCa conductance
      REAL(KIND=RKIND) :: G_pCa = 0.1238_RKIND     ! units: pA/pF
c!     G_pCa for spiral wave breakup
c      REAL(KIND=RKIND) :: G_pCa = 0.8666_RKIND     ! units: pA/pF
!     K_pCa: Half-saturation constant of I_pCa
      REAL(KIND=RKIND) :: K_pCa = 5.E-4_RKIND      ! units: mM
!     G_bNa: Maximal I_bNa conductance
      REAL(KIND=RKIND) :: G_bNa = 2.9E-4_RKIND     ! units: nS/pF
!     G_bCa: Maximal I_bCa conductance
      REAL(KIND=RKIND) :: G_bCa = 5.92E-4_RKIND    ! units: nS/pF
!     Vmax_up: Maximal I_up conductance
      REAL(KIND=RKIND) :: Vmax_up = 6.375E-3_RKIND ! units: mM/ms
!     K_up: Half-saturation constant of I_up
      REAL(KIND=RKIND) :: K_up = 2.5E-4_RKIND      ! units: mM
!     V_rel: Maximal I_rel conductance
      REAL(KIND=RKIND) :: V_rel = 0.102_RKIND      ! units: mM/ms
!     k1p: R to O and RI to I, I_rel transition rate
      REAL(KIND=RKIND) :: k1p = 0.15_RKIND         ! units: mM^{-2}/ms
!     k2p: O to I and R to RI, I_rel transition rate
      REAL(KIND=RKIND) :: k2p = 4.5E-2_RKIND       ! units: mM^{-1}/ms
!     k3: O to R and I to RI, I_rel transition rate
      REAL(KIND=RKIND) :: k3 = 6.E-2_RKIND         ! units: ms^{-1}
!     k4: I to O and Ri to I, I_rel transition rate
      REAL(KIND=RKIND) :: k4 = 5.E-3_RKIND         ! units: ms^{-1}
!     EC: Ca_sr half-saturation constant of k_casr
      REAL(KIND=RKIND) :: EC = 1.5_RKIND           ! units: mM
!     max_sr: Maximum value of k_casr
      REAL(KIND=RKIND) :: max_sr = 2.5_RKIND       ! dimensionless
!     min_sr: Minimum value of k_casr
      REAL(KIND=RKIND) :: min_sr = 1._RKIND        ! dimensionless
!     V_leak: Maximal I_leak conductance
      REAL(KIND=RKIND) :: V_leak = 3.6E-4_RKIND    ! units: mM/ms
!     V_xfer: Maximal I_xfer conductance
      REAL(KIND=RKIND) :: V_xfer = 3.8E-3_RKIND    ! units: mM/ms
!     Buf_c: Total cytoplasmic buffer concentration
      REAL(KIND=RKIND) :: Buf_c = 0.2_RKIND        ! units: mM
!     K_bufc: Ca_i half-saturation constant for cytplasmic buffer
      REAL(KIND=RKIND) :: K_bufc = 1.E-3_RKIND     ! units: mM
!     Buf_sr: Total sacroplasmic buffer concentration
      REAL(KIND=RKIND) :: Buf_sr = 10._RKIND       ! units: mM
!     K_bufsr: Ca_sr half-saturation constant for subspace buffer
      REAL(KIND=RKIND) :: K_bufsr = 0.3_RKIND      ! units: mM
!     Buf_ss: Total subspace buffer concentration
      REAL(KIND=RKIND) :: Buf_ss = 0.4_RKIND       ! units: mM
!     K_bufss: Ca_ss half-saturation constant for subspace buffer
      REAL(KIND=RKIND) :: K_bufss = 2.5E-4_RKIND   ! units: mM
!     Resting potential
      REAL(KIND=RKIND) :: Vrest = -85.23_RKIND     ! units: mV
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Ca_rest: Resting Ca concentration
      REAL(KIND=RKIND) :: Ca_rest = 5.E-5_RKIND    ! units: mM
!     Ca_crit: Critical Ca concentration
      REAL(KIND=RKIND) :: Ca_crit = 8.E-4_RKIND    ! units: mM
!     K_T: Saturation of concentration
      REAL(KIND=RKIND) :: K_T = 12.5_RKIND         ! units: MPa/mM
!     eps_0: Minimum activation
      REAL(KIND=RKIND) :: eps_0 = 0.1_RKIND        ! units: ms^{-1}
!     eps_i: Maximum activation
      REAL(KIND=RKIND) :: eps_i = 1._RKIND         ! units: ms^{-1}
!     Transition rate
      REAL(KIND=RKIND) :: xi_T = 4.E3_RKIND        ! units: mM^{-1}
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active strain model
!     Active force of sacromere (-mM^{-2})
      REAL(KIND=RKIND) :: alfa = -4.E6_RKIND
!     Resting Ca concentration (mM)
      REAL(KIND=RKIND) :: c_Ca0 = 2.155E-4_RKIND
!     Viscous-type constant (ms-mM^{-2})
      REAL(KIND=RKIND) :: mu_Ca = 5.E9_RKIND

!     Force-length relationship parameters
!     Initial length of sacromeres (um)
      REAL(KIND=RKIND) :: SL0 = 1.95_RKIND
!     Min. length of sacromeres (um)
      REAL(KIND=RKIND) :: SLmin = 1.7_RKIND
!     Max. length of sacromeres (um)
      REAL(KIND=RKIND) :: SLmax = 2.6_RKIND
!     Fourier coefficients
      REAL(KIND=RKIND) :: f0  = -4333.618335582119_RKIND
      REAL(KIND=RKIND) :: fs1 =  2570.395355352195_RKIND
      REAL(KIND=RKIND) :: fs2 =  1329.53611689133_RKIND
      REAL(KIND=RKIND) :: fs3 =  104.943770305116_RKIND
      REAL(KIND=RKIND) :: fc1 = -2051.827278991976_RKIND
      REAL(KIND=RKIND) :: fc2 =  302.216784558222_RKIND
      REAL(KIND=RKIND) :: fc3 =  218.375174229422_RKIND

!-----------------------------------------------------------------------
!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!-----------------------------------------------------------------------
!     Variables
!     Reverse potentials for Na, K, Ca
      REAL(KIND=RKIND) :: E_Na
      REAL(KIND=RKIND) :: E_K
      REAL(KIND=RKIND) :: E_Ca
      REAL(KIND=RKIND) :: E_Ks
!     Cellular transmembrane currents
!     I_Na: Fast sodium current
      REAL(KIND=RKIND) :: I_Na
!     I_K1: inward rectifier outward current
      REAL(KIND=RKIND) :: I_K1
!     I_to: transient outward current
      REAL(KIND=RKIND) :: I_to
!     I_Kr: rapid delayed rectifier current
      REAL(KIND=RKIND) :: I_Kr
!     I_Ks: slow delayed rectifier current
      REAL(KIND=RKIND) :: I_Ks
!     I_CaL: L-type Ca current
      REAL(KIND=RKIND) :: I_CaL
!     I_NaCa: Na-Ca exchanger current
      REAL(KIND=RKIND) :: I_NaCa
!     I_NaK: Na-K pump current
      REAL(KIND=RKIND) :: I_NaK
!     I_pCa: plateau Ca current
      REAL(KIND=RKIND) :: I_pCa
!     I_pK: plateau K current
      REAL(KIND=RKIND) :: I_pK
!     I_bCa: background Ca current
      REAL(KIND=RKIND) :: I_bCa
!     I_lean: background Na current
      REAL(KIND=RKIND) :: I_bNa
!     I_leak: sacroplasmic reticulum Ca leak current
      REAL(KIND=RKIND) :: I_leak
!     I_up: sacroplasmic reticulum Ca pump current
      REAL(KIND=RKIND) :: I_up
!     I_rel: Ca induced Ca release current
      REAL(KIND=RKIND) :: I_rel
!     I_xfer: diffusive Ca current
      REAL(KIND=RKIND) :: I_xfer
!-----------------------------------------------------------------------
!     State variables
      REAL(KIND=RKIND) :: V
      REAL(KIND=RKIND) :: K_i
      REAL(KIND=RKIND) :: Na_i
      REAL(KIND=RKIND) :: Ca_i
      REAL(KIND=RKIND) :: Ca_ss
      REAL(KIND=RKIND) :: Ca_sr
      REAL(KIND=RKIND) :: R_bar

!     Gating variables (runtime, steady state)
      REAL(KIND=RKIND) :: xr1, xr1i
      REAL(KIND=RKIND) :: xr2, xr2i
      REAL(KIND=RKIND) :: xs, xsi
      REAL(KIND=RKIND) :: m, mi
      REAL(KIND=RKIND) :: h, hi
      REAL(KIND=RKIND) :: j, ji
      REAL(KIND=RKIND) :: d, di
      REAL(KIND=RKIND) :: f, fi
      REAL(KIND=RKIND) :: f2, f2i
      REAL(KIND=RKIND) :: fcass, fcassi
      REAL(KIND=RKIND) :: s, si
      REAL(KIND=RKIND) :: r, ri

!     Other variables
      REAL(KIND=RKIND) :: k1
      REAL(KIND=RKIND) :: k2
      REAL(KIND=RKIND) :: k_casr
      REAL(KIND=RKIND) :: O

!     Jacobian variables
      REAL(KIND=RKIND) :: E_Na_Nai, E_K_Ki, E_Ca_Cai, E_Ks_Ki, E_Ks_Nai
      REAL(KIND=RKIND) :: I_Na_V, I_Na_Nai
      REAL(KIND=RKIND) :: I_to_V, I_to_Ki
      REAL(KIND=RKIND) :: I_K1_V, I_K1_Ki
      REAL(KIND=RKIND) :: I_Kr_V, I_Kr_Ki
      REAL(KIND=RKIND) :: I_Ks_V, I_Ks_Ki, I_Ks_Nai
      REAL(KIND=RKIND) :: I_CaL_V, I_CaL_Cass
      REAL(KIND=RKIND) :: I_NaCa_V, I_NaCa_Nai, I_NaCa_Cai
      REAL(KIND=RKIND) :: I_NaK_V, I_NaK_Nai
      REAL(KIND=RKIND) :: I_pCa_Cai
      REAL(KIND=RKIND) :: I_pK_V, I_pK_Ki
      REAL(KIND=RKIND) :: I_bCa_V, I_bCa_Cai
      REAL(KIND=RKIND) :: I_bNa_V, I_bNa_Nai
      REAL(KIND=RKIND) :: I_leak_Cai, I_leak_Casr
      REAL(KIND=RKIND) :: I_up_Cai
      REAL(KIND=RKIND) :: I_rel_Cass, I_rel_Casr, I_rel_Rbar
      REAL(KIND=RKIND) :: I_xfer_Cai, I_xfer_Cass
      REAL(KIND=RKIND) :: k_casr_sr, k1_casr, O_Casr, O_Cass, O_Rbar
!#######################################################################

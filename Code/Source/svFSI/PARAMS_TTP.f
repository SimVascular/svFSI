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
!     Constants for TenTusscher-Panfilov Ventricular Myocyte Model.
!
!--------------------------------------------------------------------

!     Default model parameters
!     R: Gas constant
      REAL(KIND=8) :: Rc = 8314.472D0       ! units: J/mol/K
!     T: Temperature
      REAL(KIND=8) :: Tc = 310.0D0          ! units: K
!     F: Faraday constant
      REAL(KIND=8) :: Fc = 96485.3415D0     ! units: C/mmol
!     Cm: Cell capacitance per unit surface area
      REAL(KIND=8) :: Cm = 0.185D0         ! units: uF/cm^{2}
!     sV: Surface to volume ratio
      REAL(KIND=8) :: sV = 0.2D0            ! units: um^{-1}
!     rho: Cellular resistivity
      REAL(KIND=8) :: rho = 162.0D0        ! units: \Omega-cm
!     V_c: Cytoplasmic volume
      REAL(KIND=8) :: V_c = 16.404D-3      ! units: um^{3}
!     V_sr: Sacroplasmic reticulum volume
      REAL(KIND=8) :: V_sr = 1.094D-3      ! units: um^{3}
!     V_ss: Subspace volume
      REAL(KIND=8) :: V_ss = 5.468D-5      ! units: um^{3}
!     K_o: Extracellular K concentration
      REAL(KIND=8) :: K_o = 5.4D0          ! units: mM
!     Na_o: Extracellular Na concentration
      REAL(KIND=8) :: Na_o = 140.0D0       ! units: mM
!     Ca_o: Extracellular Ca concentration
      REAL(KIND=8) :: Ca_o = 2.0D0         ! units: mM
!     G_Na: Maximal I_Na conductance
      REAL(KIND=8) :: G_Na = 14.838D0      ! units: nS/pF
!     G_K1: Maximal I_K1 conductance
      REAL(KIND=8) :: G_K1 = 5.405D0       ! units: nS/pF
!     I_to: Maximal epicardial I_to conductance
      REAL(KIND=8) :: G_to = 0.294D0       ! units: nS/pF
!     G_Kr: Maximal I_Kr conductance
      REAL(KIND=8) :: G_Kr = 0.153D0       ! units: nS/pF
!     G_Ks: Maximal epicardial I_Ks conductance
      REAL(KIND=8) :: G_Ks = 0.392D0       ! units: nS/pF
!     p_KNa: Relative I_Ks permeability to Na
      REAL(KIND=8) :: p_KNa = 0.03D0       ! dimensionless
!     G_CaL: Maximal I_CaL conductance
      REAL(KIND=8) :: G_CaL = 3.98D-5      ! units: cm^{3}/uF/ms
!     K_NaCa: Maximal I_NaCa
      REAL(KIND=8) :: K_NaCa = 1000.0D0    ! units: pA/pF
!     gamma: Voltage dependent parameter of I_NaCa
      REAL(KIND=8) :: gamma = 0.35D0       ! dimensionless
!     K_mCa: Ca_i half-saturation constant for I_NaCa
      REAL(KIND=8) :: K_mCa = 1.38D0       ! units: mM
!     K_mNai: Na_i half-saturation constant for I_NaCa
      REAL(KIND=8) :: K_mNai = 87.5D0      ! units: mM
!     K_sat: Saturation factor for I_NaCa
      REAL(KIND=8) :: K_sat = 0.1D0        ! dimensionless
!     alpha: Factor enhancing outward nature of I_NaCa
      REAL(KIND=8) :: alpha = 2.5D0        ! dimensionless
!     p_NaK: Maximal I_NaK
      REAL(KIND=8) :: p_NaK = 2.724D0      ! units: pA/pF
!     K_mK: K_o half-saturation constant of I_NaK
      REAL(KIND=8) :: K_mK = 1.0D0         ! units: mM
!     K_mNa: Na_i half-saturation constant of I_NaK
      REAL(KIND=8) :: K_mNa = 40.0D0       ! units: mM
!     G_pK: Maximal I_pK conductance
      REAL(KIND=8) :: G_pK = 1.46D-2       ! units: nS/pF
!     G_pCa: Maximal I_pCa conductance
      REAL(KIND=8) :: G_pCa = 0.1238D0     ! units: pA/pF
!     K_pCa: Half-saturation constant of I_pCa
      REAL(KIND=8) :: K_pCa = 5.0D-4       ! units: mM
!     G_bNa: Maximal I_bNa conductance
      REAL(KIND=8) :: G_bNa = 2.9D-4       ! units: nS/pF
!     G_bCa: Maximal I_bCa conductance
      REAL(KIND=8) :: G_bCa = 5.92D-4      ! units: nS/pF
!     Vmax_up: Maximal I_up conductance
      REAL(KIND=8) :: Vmax_up = 6.375D-3   ! units: mM/ms
!     K_up: Half-saturation constant of I_up
      REAL(KIND=8) :: K_up = 2.5D-4        ! units: mM
!     V_rel: Maximal I_rel conductance
      REAL(KIND=8) :: V_rel = 0.102D0      ! units: mM/ms
!     k1p: R to O and RI to I, I_rel transition rate
      REAL(KIND=8) :: k1p = 0.15D0         ! units: mM^{-2}/ms
!     k2p: O to I and R to RI, I_rel transition rate
      REAL(KIND=8) :: k2p = 0.045D0        ! units: mM^{-1}/ms
!     k3: O to R and I to RI, I_rel transition rate
      REAL(KIND=8) :: k3  = 0.06D0         ! units: ms^{-1}
!     k4: I to O and Ri to I, I_rel transition rate
      REAL(KIND=8) :: k4  = 5.0D-3         ! units: ms^{-1}
!     EC: Ca_sr half-saturation constant of k_casr
      REAL(KIND=8) :: EC  = 1.5D0          ! units: mM
!     max_sr: Maximum value of k_casr
      REAL(KIND=8) :: max_sr = 2.5D0       ! dimensionless
!     min_sr: Minimum value of k_casr
      REAL(KIND=8) :: min_sr = 1.0D0       ! dimensionless
!     V_leak: Maximal I_leak conductance
      REAL(KIND=8) :: V_leak = 3.6D-4      ! units: mM/ms
!     V_xfer: Maximal I_xfer conductance
      REAL(KIND=8) :: V_xfer = 3.8D-3      ! units: mM/ms
!     Buf_c: Total cytoplasmic buffer concentration
      REAL(KIND=8) :: Buf_c = 0.2D0        ! units: mM
!     K_bufc: Ca_i half-saturation constant for cytplasmic buffer
      REAL(KIND=8) :: K_bufc = 1.0D-3      ! units: mM
!     Buf_sr: Total sacroplasmic buffer concentration
      REAL(KIND=8) :: Buf_sr = 10.0D0      ! units: mM
!     K_bufsr: Ca_sr half-saturation constant for subspace buffer
      REAL(KIND=8) :: K_bufsr = 0.3D0      ! units: mM
!     Buf_ss: Total subspace buffer concentration
      REAL(KIND=8) :: Buf_ss = 0.4D0       ! units: mM
!     K_bufss: Ca_ss half-saturation constant for subspace buffer
      REAL(KIND=8) :: K_bufss = 2.5D-4     ! units: mM
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters
!     Ca_rest: Resting Ca concentration
      REAL(KIND=8) :: Ca_rest = 5.0D-5     ! units: mM
!     Ca_crit: Critical Ca concentration
      REAL(KIND=8) :: Ca_crit = 8.0D-4     ! units: mM
!     eta_T: Saturation of concentration
      REAL(KIND=8) :: eta_T = 12.5D0       ! units: MPa/mM
!     eps_0: Minimum activation
      REAL(KIND=8) :: eps_0 = 0.1D0        ! units: ms^{-1}
!     eps_i: Maximum activation
      REAL(KIND=8) :: eps_i = 1.0D0        ! units: ms^{-1}
!     Transition rate
      REAL(KIND=8) :: xi_T = 4.0D3         ! units: mM^{-1}
!-----------------------------------------------------------------------
!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 1.0D0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 1.0D0
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = 0.0D0
!-----------------------------------------------------------------------
!     Variables
!     Reverse potentials for Na, K, Ca
      REAL(KIND=8) :: E_Na
      REAL(KIND=8) :: E_K
      REAL(KIND=8) :: E_Ca
      REAL(KIND=8) :: E_Ks
!     Cellular transmembrane currents
!     I_Na: Fast sodium current
      REAL(KIND=8) :: I_Na
!     I_K1: inward rectifier outward current
      REAL(KIND=8) :: I_K1
!     I_to: transient outward current
      REAL(KIND=8) :: I_to
!     I_Kr: rapid delayed rectifier current
      REAL(KIND=8) :: I_Kr
!     I_Ks: slow delayed rectifier current
      REAL(KIND=8) :: I_Ks
!     I_CaL: L-type Ca current
      REAL(KIND=8) :: I_CaL
!     I_NaCa: Na-Ca exchanger current
      REAL(KIND=8) :: I_NaCa
!     I_NaK: Na-K pump current
      REAL(KIND=8) :: I_NaK
!     I_pCa: plateau Ca current
      REAL(KIND=8) :: I_pCa
!     I_pK: plateau K current
      REAL(KIND=8) :: I_pK
!     I_bCa: background Ca current
      REAL(KIND=8) :: I_bCa
!     I_lean: background Na current
      REAL(KIND=8) :: I_bNa
!     I_leak: sacroplasmic reticulum Ca leak current
      REAL(KIND=8) :: I_leak
!     I_up: sacroplasmic reticulum Ca pump current
      REAL(KIND=8) :: I_up
!     I_rel: Ca induced Ca release current
      REAL(KIND=8) :: I_rel
!     I_xfer: diffusive Ca current
      REAL(KIND=8) :: I_xfer
!-----------------------------------------------------------------------
!     Other auxillary variables
      REAL(KIND=8) :: V
      REAL(KIND=8) :: K_i
      REAL(KIND=8) :: Na_i
      REAL(KIND=8) :: Ca_i
      REAL(KIND=8) :: xr1, xr1i
      REAL(KIND=8) :: xr2, xr2i
      REAL(KIND=8) :: xs, xsi
      REAL(KIND=8) :: m, mi
      REAL(KIND=8) :: h, hi
      REAL(KIND=8) :: j, ji
      REAL(KIND=8) :: Ca_ss
      REAL(KIND=8) :: d, di
      REAL(KIND=8) :: f, fi
      REAL(KIND=8) :: f2, f2i
      REAL(KIND=8) :: fcass, fcassi
      REAL(KIND=8) :: s, si
      REAL(KIND=8) :: r, ri
      REAL(KIND=8) :: Ca_sr
      REAL(KIND=8) :: R_bar
      REAL(KIND=8) :: k1
      REAL(KIND=8) :: k2
      REAL(KIND=8) :: k_casr
      REAL(KIND=8) :: O

!     Jacobian variables
      REAL(KIND=8) :: E_Na_Nai, E_K_Ki, E_Ca_Cai, E_Ks_Ki, E_Ks_Nai
      REAL(KINd=8) :: I_Na_V, I_Na_Nai, I_Na_m, I_Na_h, I_Na_j
      REAL(KIND=8) :: I_to_V, I_to_Ki, I_to_s, I_to_r
      REAL(KIND=8) :: I_K1_V, I_K1_Ki
      REAL(KIND=8) :: I_Kr_V, I_Kr_Ki, I_Kr_xr1, I_Kr_xr2
      REAL(KIND=8) :: I_Ks_V, I_Ks_Ki, I_Ks_Nai, I_Ks_xs
      REAL(KIND=8) :: I_CaL_V, I_CaL_Cass, I_CaL_d, I_CaL_f, I_CaL_f2,
     2   I_CaL_fcass
      REAL(KIND=8) :: I_NaCa_V, I_NaCa_Nai, I_NaCa_Cai
      REAL(KIND=8) :: I_NaK_V, I_NaK_Nai
      REAL(KIND=8) :: I_pCa_Cai
      REAL(KIND=8) :: I_pK_V, I_pK_Ki
      REAL(KIND=8) :: I_bCa_V, I_bCa_Cai
      REAL(KIND=8) :: I_bNa_V, I_bNa_Nai
      REAL(KIND=8) :: I_leak_Cai, I_leak_Casr
      REAL(KIND=8) :: I_up_Cai
      REAL(KIND=8) :: I_rel_Cass, I_rel_Casr, I_rel_Rbar
      REAL(KIND=8) :: I_xfer_Cai, I_xfer_Cass
      REAL(KIND=8) :: k_casr_sr, k1_casr, O_Casr, O_Cass, O_Rbar
!#######################################################################

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
!     Parameters for Bueno-Orovio cellular activation model.
!     Parameters are chosen based on below references.
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

!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale = 85.7_RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = -84._RKIND
!--------------------------------------------------------------------
!     Model parameters (epi, endo, myo)
      REAL(KIND=RKIND) :: u_o(3) =
     2   (/0._RKIND, 0._RKIND, 0._RKIND/)
      REAL(KIND=RKIND) :: u_u(3) =
     2   (/1.55_RKIND, 1.56_RKIND, 1.61_RKIND/)
      REAL(KIND=RKIND) :: theta_v(3) =
     2   (/0.3_RKIND, 0.3_RKIND, 0.3_RKIND/)
      REAL(KIND=RKIND) :: theta_w(3) =
     2   (/0.13_RKIND, 0.13_RKIND, 0.13_RKIND/)
      REAL(KIND=RKIND) :: thetam_v(3) =
     2   (/6.E-3_RKIND, 0.2_RKIND, 0.1_RKIND/)
      REAL(KIND=RKIND) :: theta_o(3) =
     2   (/6.E-3_RKIND, 6.E-3_RKIND, 5.E-3_RKIND/)
      REAL(KIND=RKIND) :: taum_v1(3) =
     2   (/60._RKIND, 75._RKIND, 80._RKIND/)
      REAL(KIND=RKIND) :: taum_v2(3) =
     2   (/1.15E3_RKIND, 10._RKIND, 1.4506_RKIND/)
      REAL(KIND=RKIND) :: taup_v(3) =
     2   (/1.4506_RKIND, 1.4506_RKIND, 1.4506_RKIND/)
      REAL(KIND=RKIND) :: taum_w1(3) =
     2   (/60._RKIND, 6._RKIND, 70._RKIND/)
      REAL(KIND=RKIND) :: taum_w2(3) =
     2   (/15._RKIND, 140._RKIND, 8._RKIND/)
      REAL(KIND=RKIND) :: km_w(3) =
     2   (/65._RKIND, 200._RKIND, 200._RKIND/)
      REAL(KIND=RKIND) :: um_w(3) =
     2   (/3.E-2_RKIND, 1.6E-2_RKIND, 1.6E-2_RKIND/)
      REAL(KIND=RKIND) :: taup_w(3) =
     2   (/200._RKIND, 280._RKIND, 280._RKIND/)
      REAL(KIND=RKIND) :: tau_fi(3) =
     2   (/0.11_RKIND, 0.1_RKIND, 0.078_RKIND/)
      REAL(KIND=RKIND) :: tau_o1(3) =
     2   (/400._RKIND, 470._RKIND, 410._RKIND/)
      REAL(KIND=RKIND) :: tau_o2(3) =
     2   (/6._RKIND, 6._RKIND, 7._RKIND/)
      REAL(KIND=RKIND) :: tau_so1(3) =
     2   (/30.0181_RKIND, 40._RKIND, 91._RKIND/)
      REAL(KIND=RKIND) :: tau_so2(3) =
     2   (/0.9957_RKIND, 1.2_RKIND, 0.8_RKIND/)
      REAL(KIND=RKIND) :: k_so(3) =
     2   (/2.0458_RKIND, 2._RKIND, 2.1_RKIND/)
      REAL(KIND=RKIND) :: u_so(3) =
     2   (/0.65_RKIND, 0.65_RKIND, 0.6_RKIND/)
      REAL(KIND=RKIND) :: tau_s1(3) =
     2   (/2.7342_RKIND, 2.7342_RKIND, 2.7342_RKIND/)
      REAL(KIND=RKIND) :: tau_s2(3) =
     2   (/16._RKIND, 2._RKIND, 2._RKIND/)
      REAL(KIND=RKIND) :: k_s(3) =
     2   (/2.0994_RKIND, 2.0994_RKIND, 2.0994_RKIND/)
      REAL(KIND=RKIND) :: u_s(3) =
     2   (/0.9087_RKIND, 0.9087_RKIND, 0.9087_RKIND/)
      REAL(KIND=RKIND) :: tau_si(3) =
     2   (/1.8875_RKIND, 2.9013_RKIND, 3.3849_RKIND/)
      REAL(KIND=RKIND) :: tau_winf(3) =
     2   (/7.E-2_RKIND, 2.73E-2_RKIND, 1.E-2_RKIND/)
      REAL(KIND=RKIND) :: ws_inf(3) =
     2   (/0.94_RKIND, 0.78_RKIND, 0.5_RKIND/)
!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Resting voltage (mV)
      REAL(KIND=RKIND) :: Vrest = -84._RKIND
!     Critical voltage (mV)
      REAL(KIND=RKIND) :: Vcrit = -30._RKIND
!     Saturation potential (Pa/mV)
      REAL(KIND=RKIND) :: K_T = 5.0E3_RKIND
!     Minimum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_0 = 0.1_RKIND
!     Maximum activation (ms^{-1})
      REAL(KIND=RKIND) :: eps_i = 1._RKIND
!     Transition rate (mV^{-1})
      REAL(KIND=RKIND) :: xi_T = 1._RKIND
!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active strain model
!     Active force of sacromere (-mM^{-2})
      REAL(KIND=RKIND) :: alFa = -4.E+6_RKIND
!     Resting Ca concentration (mM) := slow inward current variable (s)
      REAL(KIND=RKIND) :: c_0 = 2.155E-4_RKIND
!     Viscous-type constant (ms-mM^{-2})
      REAL(KIND=RKIND) :: mu_c = 5.E9_RKIND

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

!--------------------------------------------------------------------
!     Cm: Cell capacitance per unit surface area
      REAL(KIND=RKIND) :: Cm  = 1._RKIND
!     sV: Surface to volume ratio
      REAL(KIND=RKIND) :: sV  = 1._RKIND
!     rho: Cellular resistivity
      REAL(KIND=RKIND) :: rho = 1._RKIND
!####################################################################

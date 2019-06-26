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
!     Parameters used for Bueno-Orovio Ventricular Myocyte Model.
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 85.7D0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 1.0D0
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = -84.0D0
!-----------------------------------------------------------------------
!     Model parameters (epi, endo, myo)
      REAL(KIND=8) :: u_o(3)      = (/0.0D0,     0.0D0,    0.0D0   /)
      REAL(KIND=8) :: u_u(3)      = (/1.55D0,    1.56D0,   1.61D0  /)
      REAL(KIND=8) :: theta_v(3)  = (/0.3D0,     0.3D0,    0.3D0   /)
      REAL(KIND=8) :: theta_w(3)  = (/0.13D0,    0.13D0,   0.13D0  /)
      REAL(KIND=8) :: thetam_v(3) = (/6.0D-3,    0.2D0,    0.1D0   /)
      REAL(KIND=8) :: theta_o(3)  = (/6.0D-3,    6.0D-3,   5.0D-3  /)
      REAL(KIND=8) :: taum_v1(3)  = (/60.0D0,    75.0D0,   80.0D0  /)
      REAL(KIND=8) :: taum_v2(3)  = (/1.15D3,    10.0D0,   1.4506D0/)
      REAL(KIND=8) :: taup_v(3)   = (/1.4506D0,  1.4506D0, 1.4506D0/)
      REAL(KIND=8) :: taum_w1(3)  = (/60.0D0,    6.0D0,    70.0D0  /)
      REAL(KIND=8) :: taum_w2(3)  = (/15.0D0,    140.0D0,  8.0D0   /)
      REAL(KIND=8) :: km_w(3)     = (/65.0D0,    200.0D0,  200.0D0 /)
      REAL(KIND=8) :: um_w(3)     = (/0.03D0,    0.016D0,  0.016D0 /)
      REAL(KIND=8) :: taup_w(3)   = (/200.0D0,   280.0D0,  280.0D0 /)
      REAL(KIND=8) :: tau_fi(3)   = (/0.11D0,    0.1D0,    0.078D0 /)
      REAL(KIND=8) :: tau_o1(3)   = (/400.0D0,   470.0D0,  410.0D0 /)
      REAL(KIND=8) :: tau_o2(3)   = (/6.0D0,     6.0D0,    7.0D0   /)
      REAL(KIND=8) :: tau_so1(3)  = (/30.0181D0, 40.0D0,   91.0D0  /)
      REAL(KIND=8) :: tau_so2(3)  = (/0.9957D0,  1.2D0,    0.8D0   /)
      REAL(KIND=8) :: k_so(3)     = (/2.0458D0,  2.0D0,    2.1D0   /)
      REAL(KIND=8) :: u_so(3)     = (/0.65D0,    0.65D0,   0.6D0   /)
      REAL(KIND=8) :: tau_s1(3)   = (/2.7342D0,  2.7342D0, 2.7342D0/)
      REAL(KIND=8) :: tau_s2(3)   = (/16.0D0,    2.0D0,    2.0D0   /)
      REAL(KIND=8) :: k_s(3)      = (/2.0994D0,  2.0994D0, 2.0994D0/)
      REAL(KIND=8) :: u_s(3)      = (/0.9087D0,  0.9087D0, 0.9087D0/)
      REAL(KIND=8) :: tau_si(3)   = (/1.8875D0,  2.9013D0, 3.3849D0/)
      REAL(KIND=8) :: tau_winf(3) = (/0.07D0,    0.0273D0, 0.01D0  /)
      REAL(KIND=8) :: ws_inf(3)   = (/0.94D0,    0.78D0,   0.5D0   /)
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Resting voltage (mV)
      REAL(KIND=8) :: Vrest = -84.0D0
!     Critical voltage (mV)
      REAL(KIND=8) :: Vcrit = -30.0D0
!     Saturation potential
      REAL(KIND=8) :: eta_T = 0.005D0
!     Minimum activation (ms^{-1})
      REAL(KIND=8) :: eps_0 = 0.1D0
!     Maximum activation (ms^{-1})
      REAL(KIND=8) :: eps_i = 1.0D0
!     Transition rate (mV^{-1})
      REAL(KIND=8) :: xi_T  = 1.0D0
!-----------------------------------------------------------------------
!     Electromechanics coupling parameters: active strain model
!     Active force of sacromere (-mM^{-2})
      REAL(KIND=8) :: alFa = -4.0D6
!     Resting Ca concentration (mM) := slow inward current variable (s)
      REAL(KIND=8) :: c0 = 2.155D-4
!     Viscous-type constant (ms-mM^{-2})
      REAL(KIND=8) :: mu_C = 5.0D6

!     Force-length relationship parameters
!     Initial length of sacromeres (um)
      REAL(KIND=8) :: SL0 = 1.95D0
!     Min. length of sacromeres (um)
      REAL(KIND=8) :: SLmin = 1.7D0
!     Max. length of sacromeres (um)
      REAL(KIND=8) :: SLmax = 2.6D0
!     Fourier coefficients
      REAL(KIND=8) :: f0  = -4333.618335582119D0
      REAL(KIND=8) :: fc1 =  2570.395355352195D0
      REAL(KIND=8) :: fs1 = -2051.827278991976D0
      REAL(KIND=8) :: fc2 =  1329.53611689133D0
      REAL(KIND=8) :: fs2 =  302.216784558222D0
      REAL(KIND=8) :: fc3 =  104.943770305116D0
      REAL(KIND=8) :: fs3 =  218.375174229422D0

!-----------------------------------------------------------------------
!     Cm: Cell capacitance per unit surface area
      REAL(KIND=8) :: Cm  = 1.0D0
!     sV: Surface to volume ratio
      REAL(KIND=8) :: sV  = 1.0D0
!     rho: Cellular resistivity
      REAL(KIND=8) :: rho = 1.0D0
!#######################################################################

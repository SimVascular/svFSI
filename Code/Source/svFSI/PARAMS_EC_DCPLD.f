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
!     Parameters for decoupled and uniformly activated excitation model
!     for excitation-contraction coupling. Parameters are chosen based
!     on below references.
!
!     Reference for active stress model:
!        Pfaller, M. R., et al.(2019). The importance of the pericardium
!        for cardiac biomechanics: from physiology to computational
!        modeling. Biomechanics and Modeling in Mechanobiology,
!        18(2), 503–529. https://doi.org/10.1007/s10237-018-1098-4
!
!     Reference for active strain model:
!        Barbarotta, L., et al.(2018). A transmurally heterogeneous
!        orthotropic activation model for ventricular contraction and
!        its numerical validation. International Journal for Numerical
!        Methods in Biomedical Engineering, 34(12), 1–24.
!        https://doi.org/10.1002/cnm.3137
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active stress model
!     Contractility (Pa)
      REAL(KIND=RKIND) :: sigm0 = 9.0E4_RKIND

!     Min activation rate (ms^-1)
      REAL(KIND=RKIND) :: min_alpha = -0.03_RKIND

!     Max activation rate (ms^-1)
      REAL(KIND=RKIND) :: max_alpha = 0.005_RKIND

!     Onset of systole (ms)
      REAL(KIND=RKIND) :: t_sys = 170._RKIND

!     Onset of diastole (ms)
      REAL(KIND=RKIND) :: t_dia = 484._RKIND

!     Gamma (ms)
      REAL(KIND=RKIND) :: gamma = 5._RKIND
!--------------------------------------------------------------------
!     Electromechanics coupling parameters: active strain model
!     Min fiber shortening
      REAL(KIND=RKIND) :: gf_min = -0.13_RKIND

!     Initial time for activation (ms)
      REAL(KIND=RKIND) :: ta_s = 170._RKIND

!     Duration of activation (ms)
      REAL(KIND=RKIND) :: ta_e = 480._RKIND
!--------------------------------------------------------------------


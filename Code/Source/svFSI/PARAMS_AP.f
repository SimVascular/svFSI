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
!     Parameters used for Aliev-Panfilov Ventricular Myocyte Model.
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 100.0D0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 12.9D0
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = -80.0D0
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=8) :: alpha = 0.01D0
      REAL(KIND=8) :: a     = 0.002D0
      REAL(KIND=8) :: b     = 0.15D0
      REAL(KIND=8) :: c     = 8.0D0
      REAL(KIND=8) :: mu1   = 0.2D0
      REAL(KIND=8) :: mu2   = 0.3D0

!     Cm: Cell capacitance per unit surface area
      REAL(KIND=8) :: Cm  = 1.0d0
!     sV: Surface to volume ratio
      REAL(KIND=8) :: sV  = 1.0D0
!     rho: Cellular resistivity
      REAL(KIND=8) :: rho = 1.0D0
!#######################################################################

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
!     Parameters used for Fitzhugh-Nagumo Myocyte Activation Model.
!
!     Reference for Aliev-Panfilov electrophysiology model:
!        Goktepe, S., & Kuhl, E. (2009). Computational modeling of
!        cardiac electrophysiology: A novel finite element approach.
!        Int. J. Numer. Meth. Engng, 79, 156–178.
!        https://doi.org/10.1002/nme
!
!--------------------------------------------------------------------

!     Scaling factors
!     Voltage scaling
      REAL(KIND=RKIND) :: Vscale  = 1._RKIND
!     Time scaling
      REAL(KIND=RKIND) :: Tscale  = 1._RKIND
!     Voltage offset parameter
      REAL(KIND=RKIND) :: Voffset = 0._RKIND
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=RKIND) :: alpha = -0.5_RKIND
      REAL(KIND=RKIND) :: a = 0._RKIND
      REAL(KIND=RKIND) :: b = -0.6_RKIND
      REAL(KIND=RKIND) :: c = 50._RKIND
!#######################################################################


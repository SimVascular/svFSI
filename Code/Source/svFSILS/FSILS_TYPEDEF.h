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
!     Data types are defined here
!
!--------------------------------------------------------------------

!     Integer, 8 bits
      INTEGER, PARAMETER :: LSIP1  = SELECTED_INT_KIND(2)
!     Integer, 16 bits
      INTEGER, PARAMETER :: LSIP2  = SELECTED_INT_KIND(4)
!     Integer, 32 bits
      INTEGER, PARAMETER :: LSIP4  = SELECTED_INT_KIND(9)
!     Integer, 64 bits
      INTEGER, PARAMETER :: LSIP8  = SELECTED_INT_KIND(18)

!     Real 32 bits (single)
      INTEGER, PARAMETER :: LSRP4  = SELECTED_REAL_KIND(6, 37)
!     Real 64 bits (double)
      INTEGER, PARAMETER :: LSRP8  = SELECTED_REAL_KIND(15, 307)
!     Real 128 bits (long double)
      INTEGER, PARAMETER :: LSRP16 = SELECTED_REAL_KIND(33, 4931)

!     Default integer precision
      INTEGER, PARAMETER :: LSIP = LSIP4
!     Default real precision
      INTEGER, PARAMETER :: LSRP = LSRP8

!####################################################################

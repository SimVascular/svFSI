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
!     This is a header file for the data structure used for 
!     communication with svFSI
!      
!--------------------------------------------------------------------

      INTEGER, PARAMETER :: cplBC_Dir = 66112, cplBC_Neu = 66113,
     2   cplBCVersion = 8
     
      TYPE cplFaceType
         SEQUENCE
         INTEGER bGrp            ! (IN)  geBC_Dir/genBC_Neu
         INTEGER xPtr            ! (USE) pointer to x
         INTEGER :: eqv = 0      ! (USE) internal genBC use
         INTEGER reserved        ! ( - ) reserved for alignment
         REAL(KIND=8) Qo         ! (IN)  flow rate at t
         REAL(KIND=8) Qn         ! (IN)  flow rate at t+dt
         REAL(KIND=8) Po         ! (IN)  pressure at t
         REAL(KIND=8) Pn         ! (IN)  pressure at t+dt
         REAL(KIND=8) y          ! (OUT) imposed flow/pressure
         CHARACTER(LEN=128) name ! (IN)  name of the face
      END TYPE cplFaceType
     

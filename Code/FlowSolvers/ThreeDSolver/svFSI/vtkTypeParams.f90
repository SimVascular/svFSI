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

!**************************************************

        module dataTypeParams
        integer, parameter :: IK1 = selected_int_kind(2)                  ! integer 8 bits
        integer, parameter :: IK2 = selected_int_kind(4)                  ! integer 16 bits
        integer, parameter :: IK4 = selected_int_kind(9)                  ! integer 32 bits
        integer, parameter :: IK8 = selected_int_kind(18)                 ! integer 64 bits
        integer, parameter :: IK = IK4                                    ! default integer type

        integer, parameter :: RK4 = selected_real_kind(6,37)              ! real 32 bits (single)
        integer, parameter :: RK8 = selected_real_kind(15,307)            ! real 64 bits (double)
        integer, parameter :: RK16 = selected_real_kind(33,4931)          ! real 128 bits (long double)
        integer, parameter :: RK = RK8                                    ! default real type (double)

        ! interface procedure to transfer bits (type cast)
        ! from any type to default type
        interface transferBits
            module procedure :: trBitsIK1, trBitsIK1A, &
                                trbitsIK2, trBitsIK2A, &
                                trBitsIK4, trBitsIK4A, &
                                trBitsIK8, trBitsIK8A, &
                                trBitsRK4, trBitsRK4A, &
                                trBitsRK8, trBitsRK8A
        end interface transferBits

        contains

            !==========================================

            subroutine trBitsIK1(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK1), intent(in) :: ikind
            integer(IK), intent(out) :: res
            integer(IK1) :: intK1

            intK1 = transfer(p1,intK1)
            res = int(intK1, kind=IK)
            
            end subroutine trBitsIK1

            !==========================================

            subroutine trBitsIK2(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK2), intent(in) :: ikind
            integer(IK), intent(out) :: res
            integer(IK2) :: intK2

            intK2 = transfer(p1,intK2)
            res = int(intK2, kind=IK)
            
            end subroutine trBitsIK2

            !==========================================

            subroutine trBitsIK4(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK), intent(out) :: res
            integer(IK4) :: intK4

            intK4 = transfer(p1,intK4)
            res = int(intK4, kind=IK)
            
            end subroutine trBitsIK4

            !==========================================

            subroutine trBitsIK8(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK), intent(out) :: res
            integer(IK8) :: intK8

            intK8 = transfer(p1,intK8)
            res = int(intK8, kind=IK)
            
            end subroutine trBitsIK8

            !==========================================

            subroutine trBitsRK4(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            real(RK), intent(out) :: res
            real(RK4) :: realK4

            realK4 = transfer(p1,realK4)
            res = real(realK4, kind=RK)
            
            end subroutine trBitsRK4

            !==========================================

            subroutine trBitsRK8(p1,np1,ikind,res)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            real(RK), intent(out) :: res
            real(RK8) :: realK8

            realK8 = transfer(p1,realK8)
            res = real(realK8, kind=RK)
            
            end subroutine trBitsRK8

            !==========================================

            subroutine trBitsIK1A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK1), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK), intent(out) :: res(nout)
            integer(IK1), dimension(:), allocatable :: intK1

            allocate(intK1(nout)); intK1=0_IK1
            intK1 = transfer(p1,intK1)
            res = int(intK1, kind=IK)
            deallocate(intK1)
            
            end subroutine trBitsIK1A

            !==========================================

            subroutine trBitsIK2A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK2), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK), intent(out) :: res(nout)
            integer(IK2), dimension(:), allocatable :: intK2

            allocate(intK2(nout)); intK2=0_IK2
            intK2 = transfer(p1,intK2)
            res = int(intK2, kind=IK)
            deallocate(intK2)
            
            end subroutine trBitsIK2A

            !==========================================

            subroutine trBitsIK4A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK), intent(out) :: res(nout)
            integer(IK4), dimension(:), allocatable :: intK4

            allocate(intK4(nout)); intK4=0_IK4
            intK4 = transfer(p1,intK4)
            res = int(intK4, kind=IK)
            deallocate(intK4)

            end subroutine trBitsIK4A

            !==========================================

            subroutine trBitsIK8A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            integer(IK), intent(out) :: res(nout)
            integer(IK8), dimension(:), allocatable :: intK8

            allocate(intK8(nout)); intK8=0_IK8
            intK8 = transfer(p1,intK8)
            res = int(intK8, kind=IK)
            deallocate(intK8)
            
            end subroutine trBitsIK8A

            !==========================================

            subroutine trBitsRK4A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK4), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            real(RK), intent(out) :: res(nout)
            real(RK4), dimension(:), allocatable :: realK4

            allocate(realK4(nout)); realK4=0._RK4
            realK4 = transfer(p1,realK4)
            res = real(realK4, kind=RK)
            deallocate(realK4)
            
            end subroutine trBitsRK4A

            !==========================================

            subroutine trBitsRK8A(p1,np1,ikind,res,nout)
            implicit none
            integer(IK4), intent(in) :: np1
            integer(IK1), intent(in) :: p1(np1)
            integer(IK8), intent(in) :: ikind
            integer(IK4), intent(in) :: nout
            real(RK), intent(out) :: res(nout)
            real(RK8), dimension(:), allocatable :: realK8

            allocate(realK8(nout)); realK8=0._RK8
            realK8 = transfer(p1,realK8)
            res = real(realK8, kind=RK)
            deallocate(realK8)
            
            end subroutine trBitsRK8A

            !==========================================

        end module dataTypeParams

!**************************************************

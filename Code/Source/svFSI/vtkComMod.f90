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

      module stdParams
        character(len=8), parameter :: ftab1="(4X,A)"
        character(len=8), parameter :: ftab2="(8X,A)"
        character(len=9), parameter :: ftab3="(12X,A)"
        character(len=9), parameter :: ftab4="(14X,A)"
        character, parameter :: eol=achar(0)
        character, parameter :: newl=achar(10)

      integer, parameter :: stdout=6
      double precision, parameter :: eps=EPSILON(eps)

      integer, parameter :: strL=400
      integer, parameter :: maxToks=30
      integer, parameter :: maxNSD=3
      end module stdParams

!**************************************************

      module genUtils
      use stdParams, only: strL, newl, maxToks, stdout
      implicit none

      interface STR
         module procedure :: ITSTR, RTSTR, DTSTR, NDTSTR
      end interface STR

      contains

         !==========================================

         function fileLineCount(fileName) result(lineCntr)
         implicit none
         character(len=*), intent(in) :: fileName
         integer :: lineCntr, ios

         lineCntr = 0
         open(10,file=trim(fileName))
         do
            read(10,*,iostat=ios)
            if(ios.lt.0) exit
            lineCntr = lineCntr + 1
         end do
         close(10)

         return
         end function fileLineCount

         !==========================================

         subroutine parseString(strng,toks,ntoks)
         implicit none
         character(len=*), intent(in) :: strng
         character(len=*), dimension(maxToks), intent(out) :: toks
         integer, intent(out) :: ntoks

         character(len=strL) :: dlmtr, token

         dlmtr = ''
         token = ''

         dlmtr = '< =">'
         ntoks = 1
         toks(1) = strtok(trim(strng), trim(dlmtr))
         toks(1) = adjustl(toks(1))
         do
            token = strtok(newl,trim(dlmtr))
            if ( token.ne.newl ) then
               ntoks = ntoks+1
               toks(ntoks) = adjustl(token)
            else
               exit
            end if
         end do

         end subroutine parseString

         !==========================================

         character(len=strL) function strtok(strng, dlms)
         implicit none
         character(len=*), intent(in) :: strng
         character(len=*), intent(in) :: dlms

         integer :: ist, iend

         integer, save :: ist0,slen
         character(len=strL), save :: str0

         if ( strng(1:1).ne.newl ) then
            ist0 = 1
            str0 = strng
            slen = len(trim(str0))
         end if

         ist = ist0
         do
            if ( (ist.le.slen) .and. (index(dlms,str0(ist:ist)).ne.0) ) then
               ist = ist+1
            else
               exit
            end if
         end do

         if ( ist.gt.slen ) then
            strtok = newl
            return
         end if

         iend = ist
         do
            if ( (iend.le.slen) .and. (index(dlms,str0(iend:iend)).eq.0) ) then
               iend = iend+1
            else
               exit
            end if
         end do

         strtok = str0(ist:iend-1)
         ist0 = iend+1
         return

         end function strtok

         !==========================================

         function getToken(toks,ntoks,kwrd) result (res)
         implicit none
         integer, intent(in) :: ntoks
         character(len=*), dimension(ntoks), intent(in) :: toks
         character(len=*), intent(in) :: kwrd
         character(len=strL) :: res

         integer :: i

         res=''
         do i=1, ntoks
            if (trim(toks(i)).eq.trim(kwrd) ) then
               res = toks(i)
               return
            end if
         end do

         return

         end function getToken

         !==========================================

         function getTokenValue(toks,ntoks,kwrd) result (res)
         implicit none
         integer, intent(in) :: ntoks
         character(len=*), dimension(ntoks), intent(in) :: toks
         character(len=*), intent(in) :: kwrd
         character(len=strL) :: res

         integer :: i

         res=''
         do i=1, ntoks
            if (trim(toks(i)).eq.trim(kwrd) ) then
               res = toks(i+1)
               return
            end if
         end do

         return

         end function getTokenValue

         !==========================================

         pure function ITSTR(iVal) result(str)
         implicit none
         integer, intent(in) :: iVal
         integer :: n,ist,j,k,slen
         character(len=256) str

         str = ''
         if ( iVal.lt.0 ) then
            str(1:1) = '-'
            ist = 2
         else
            ist = 1
         end if

         slen = 2 - max(0,sign(1,iVal)) + &
            maxval( (/ (min(abs(iVal)/10**n,1)*n,n=1, 9) /) )

         n = iVal
         do j=slen, ist, -1
            k = modulo(abs(n),10) + 1
            str(j:j) = '0123456789'(k:k)
            n = n/10
         end do

         return
         end function ITSTR

         !==========================================

         pure function RTSTR(rVal) result(str)
         implicit none
         integer, parameter :: l=8
         REAL(KIND=4), intent(in) :: rVal
         character(len=l) :: str

         str = NDTSTR(dble(rVal),l)

         return
         end function RTSTR

         !==========================================

         pure function DTSTR(dVal) result(str)
         implicit none
         integer, parameter :: l=8
         double precision, intent(in) :: dVal
         character(len=l) :: str

         str = NDTSTR(dVal,l)

         return
         end function DTSTR

         !==========================================

         pure function NDTSTR(dVal,l) result(str)
         implicit none
         integer, intent(in) :: l
         double precision, intent(in) :: dVal
         character(len=l) :: str

         integer :: i,k,ipos,cnt,ex,abex,nex
         double precision :: absd

         ! check NaN !
         if ( dVal.ne.dVal ) then
            if ( l.ge.3) then
               str = ''
               str(l-2:l) = 'NaN'
            else
               str = 'NaN'(1:l)
            end if
            return
         end if

         ! check infinity !
         absd = abs(dVal)
         if ( absd.gt.huge(absd) ) then
            if ( l.ge.8 ) then
               str = ''
               str(l-7:l) = 'Infinity'
            else
               str = 'Infinity'(1:l)
            end if
            return
         end if

         ! check zero !
         if ( absd.lt.tiny(absd) ) then
            ipos = 1
            if ( l.ge.3 ) then
               str(1:3) = '0.0'
               ipos = 4
            end if
            do i=ipos, l
               str(i:i) = '0'
            end do
            return
         end if

         ! count exp and number !
         ex = 0; nex = 0
         if ( absd.ne.0d0 ) ex = floor(log10(absd)) ! exponent !
         abex = abs(ex)

         ! no of digits in exponent !
         if ( ex.ne.0 ) nex = floor(log10(real(abex,8))) + 1

         cnt = nex+1  ! atleast 1 number before exponent !
         if ( dVal.lt.0d0 ) cnt = cnt+1 ! for - sign of number !
         if ( ex.lt.0 ) cnt = cnt+1 ! for - sign of exponent !
         if ( ex.ne.0 ) cnt = cnt+1 ! for letter E !

         if ( l.lt.cnt ) then   ! check total sum !
            do i=1, l
               str(i:i) = '*'
            end do
            return
         end if

         if ( ex.ne.0 ) then
            do ipos=l, l-nex+1, -1
               k = modulo(abex,10) + 1
               str(ipos:ipos) = '0123456789'(k:k)
               abex = abex/10
            end do
            if ( ex.lt.0 ) then
               str(ipos:ipos) = '-'
               ipos = ipos-1
            end if
            str(ipos:ipos) = 'E'
            ipos = ipos-1
         else
            ipos = l
         end if

         if ( l.gt.cnt ) then
            absd = absd*(1D1**(l-cnt-1-ex))
            do i=ipos, ipos-(l-cnt)+2, -1
               k = idint(modulo(absd,1D1)) + 1
               str(i:i) = '0123456789'(k:k)
               absd = absd/1D1
            end do
            ipos = i
            str(ipos:ipos) = '.'
            ipos = ipos-1
            k = idint(modulo(absd,1D1)) + 1
            str(ipos:ipos) = '0123456789'(k:k)
         else ! l.eq.cnt
            absd = absd*(1D1**(l-cnt-ex))
            k = floor(modulo(absd,1D1)) + 1
            str(ipos:ipos) = '0123456789'(k:k)
         end if
         if ( dVal.lt.0.0d0 ) str(1:1) = '-'

         return
         end function NDTSTR

         !==========================================

   end module genUtils

!**************************************************

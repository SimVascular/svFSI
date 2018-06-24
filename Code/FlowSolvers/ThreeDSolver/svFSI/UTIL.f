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
!     This module contains a set of routine that are handy, such as
!     string manupulation, file handeling, and norm
!
!--------------------------------------------------------------------

      MODULE UTILMOD
      IMPLICIT NONE

!     This is the standard length for all the strings used in this code
      INTEGER, PARAMETER :: stdL = 400
!     This is the standard double number precision
      INTEGER, PARAMETER :: dblPr = 14
!     This is the standard integer number precision
      INTEGER, PARAMETER :: intPr = 10

!     This is for writing binary files
      CHARACTER, PARAMETER :: eol = CHAR(10)

!     \pi value
      REAL(KIND=8), PARAMETER :: pi = 3.1415926535897932384626D0
!     minimum value for a double number that is not considered zero
      REAL(KIND=8), PARAMETER :: eps = EPSILON(eps)

!     This is a general variable that will be used as a stack
      TYPE stackType
!        Maximum length of the stack
         INTEGER :: maxN = 0
!        Current size of stack
         INTEGER :: n = 0
!        Values inside stack
         INTEGER, ALLOCATABLE :: v(:)
      END TYPE stackType

      TYPE queueType
         INTEGER :: n = 0
         INTEGER :: maxN = 0
         INTEGER, ALLOCATABLE :: v(:)
      END TYPE queueType

      TYPE fileType
         CHARACTER(LEN=stdL) :: fname = "DEFAULT_Unspecified"
      CONTAINS
         PROCEDURE :: open => OPENFILE
      END TYPE

      INTERFACE NORM
         MODULE PROCEDURE NORMS, NORMV
      END INTERFACE NORM

      INTERFACE GET
         MODULE PROCEDURE GETINT, GETHEX, GETSTR
      END INTERFACE GET

      INTERFACE PUSHSTACK
         MODULE PROCEDURE PUSHSTACKS, PUSHSTACKV
      END INTERFACE PUSHSTACK

      INTERFACE ENQUEUE
         MODULE PROCEDURE ENQUEUES, ENQUEUEV
      END INTERFACE ENQUEUE

      INTERFACE STR
         MODULE PROCEDURE :: DTSTR, NDTSTR, VDTSTR, RTSTR, NRTSTR,
     2      ITSTR, NITSTR
      END INTERFACE STR

      INTERFACE OPERATOR(//)
         MODULE PROCEDURE :: CNCSL, CNCSI, CNCIS, CNCSR
      END INTERFACE OPERATOR(//)

      INTERFACE CONV
         MODULE PROCEDURE :: CONVI, CONVR
      END INTERFACE CONV

      CONTAINS

!####################################################################
      FUNCTION OPENFILE(f) RESULT(fid)

      IMPLICIT NONE

      CLASS(fileType), INTENT(INOUT) :: f
      INTEGER fid

      LOGICAL flag

      INQUIRE (FILE=TRIM(f%fname), NUMBER=fid)
      IF (fid .NE. -1) RETURN

      INQUIRE (FILE=TRIM(f%fname), EXIST=flag)
      IF (.NOT.flag) THEN
         PRINT *, CLR("File does not exist or can not be opened: "//
     2      TRIM(f%fname))
         STOP
      END IF

      DO fid=11, 1024
         INQUIRE(UNIT=fid, OPENED=flag)
         IF (.NOT.flag) EXIT
      END DO
      OPEN(UNIT=fid, FILE=TRIM(f%fname), STATUS='OLD')

      RETURN
      END FUNCTION OPENFILE

!####################################################################
!     This function will compute second NORM of a vector
      PURE FUNCTION NORMS(U, V)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: U(:)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: V(:)
      REAL(KIND=8) NORMS

      INTEGER i, n

      n = SIZE(U)
      NORMS = 0D0
      IF (PRESENT(V)) THEN
         DO i=1, n
            NORMS = NORMS + U(i)*V(i)
         END DO
      ELSE
         DO i=1, n
            NORMS = NORMS + U(i)*U(i)
         END DO
      END IF

      RETURN
      END FUNCTION NORMS
!--------------------------------------------------------------------
      PURE FUNCTION NORMV(U, V)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: V(:,:)
      REAL(KIND=8) NORMV

      INTEGER m, n, i, j

      m = SIZE(U,1)
      n = SIZE(U,2)
      NORMV = 0D0
      IF (PRESENT(V)) THEN
         SELECT CASE(m)
         CASE(1)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i)
            END DO
         CASE(2)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
            END DO
         CASE(3)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
     2                       + U(3,i)*V(3,i)
            END DO
         CASE(4)
            DO i=1, n
               NORMV = NORMV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
     2                       + U(3,i)*V(3,i) + U(4,i)*V(4,i)
            END DO
         CASE DEFAULT
            DO i=1, n
               NORMV = NORMV + SUM(U(:,i)*V(:,i))
            END DO
         END SELECT
      ELSE
         SELECT CASE(m)
         CASE(1)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i)
            END DO
         CASE(2)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
            END DO
         CASE(3)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
     2                       + U(3,i)*U(3,i)
            END DO
         CASE(4)
            DO i=1, n
               NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
     2                       + U(3,i)*U(3,i) + U(4,i)*U(4,i)
            END DO
         CASE DEFAULT
            DO i=1, n
               DO j=1, m
                  NORMV = NORMV + U(j,i)*U(j,i)
               END DO
            END DO
         END SELECT
      END IF

      RETURN
      END FUNCTION NORMV

!####################################################################
!     This routine does the cross product for a two given vector of
!     V1 and V2.
      PURE FUNCTION CROSS(V) RESULT(U)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: V(:,:)
      REAL(KIND=8) U(SIZE(V,1))

      IF (SIZE(V,1) .EQ. 2) THEN
         U(1) =  V(2,1)
         U(2) = -V(1,1)
      ELSE
         U(1) = V(2,1)*V(3,2) - V(3,1)*V(2,2)
         U(2) = V(3,1)*V(1,2) - V(1,1)*V(3,2)
         U(3) = V(1,1)*V(2,2) - V(2,1)*V(1,2)
      END IF

      RETURN
      END FUNCTION CROSS

!####################################################################
!     This function removes the leading spaces and tabs
      PURE FUNCTION ADJUSTC(str)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: str
      CHARACTER(LEN=LEN(str)) ADJUSTC

      INTEGER i

      DO i=1, LEN(str)
         IF (str(i:i) .NE. " " .AND. str(i:i) .NE. "	") EXIT
      END DO
      IF (i .GT. LEN(str)) THEN
         ADJUSTC = ""
      ELSE
         ADJUSTC = str(i:)
      END IF

      RETURN
      END FUNCTION ADJUSTC
!####################################################################
!     This routine reads a word from "strg" and puts in the "rVal" and
!     removes the readed part from "strg". Valid charecters to read are
!     defined in "valS"
      SUBROUTINE GETSTR(strg, rVal, valS)

      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(INOUT) :: strg
      CHARACTER(LEN=stdL), INTENT(OUT) :: rVal
      CHARACTER(LEN=*), INTENT(IN) :: valS

      LOGICAL flag
      INTEGER i, j, s

      s    = 1
      flag = .FALSE.
      DO i=1, stdL
         IF (flag) THEN
            DO j=1, LEN(valS)
               IF (strg(i:i) .EQ. valS(j:j)) GOTO 005
            END DO
            EXIT
         ELSE
            DO j=1, LEN(valS)
               IF (strg(i:i) .EQ. valS(j:j)) THEN
                  s    = i
                  flag = .TRUE.
               END IF
            END DO
         END IF
 005     CONTINUE
      END DO
      IF (i .GT. stdL) THEN
         rVal = ""
         RETURN
      END IF
      rVal = strg(s:i-1)
      strg = strg(i:)

      RETURN
      END SUBROUTINE GETSTR

!--------------------------------------------------------------------
!     This reads a Hexadecimal number from the begining of "str" and
!     returns the results in "r" and removes the readed part from "str"
      SUBROUTINE GETHEX(str, r)

      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(INOUT) :: str
      INTEGER(KIND=8), INTENT(OUT) :: r

      CHARACTER(LEN=*), PARAMETER :: valS="+-1234567890ABCDEFabcdefZz"
      CHARACTER(LEN=stdL) tmp

      CALL GETSTR(str, tmp, valS)
      READ(tmp,"(Z8)") r

      RETURN
      END SUBROUTINE GETHEX

!--------------------------------------------------------------------
!     This reads an integer from the begining of "str" and returns
!     the results in "r" and removes the readed part from "str"
      SUBROUTINE GETINT(str, r)

      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(INOUT) :: str
      INTEGER, INTENT(OUT) :: r

      CHARACTER(LEN=*), PARAMETER :: valS="+-1234567890"
      CHARACTER(LEN=stdL) tmp

      CALL GETSTR(str, tmp, valS)
      READ(tmp,"(I8)") r

      RETURN
      END SUBROUTINE GETINT

!####################################################################
!     This function will return the current time in sec
      PURE FUNCTION SGN(u)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: u
      INTEGER SGN

      IF (ISZERO(u)) THEN
         SGN = 0
      ELSE IF (u .GT. 0D0) THEN
         SGN = 1
      ELSE
         SGN = -1
      END IF

      RETURN
      END FUNCTION SGN

!####################################################################
!      Returning number of data/words in a string
      PURE FUNCTION CheckNoNumbers(sTmp)

      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(IN) :: sTmp

      LOGICAL isBnk
      INTEGER CheckNoNumbers, i

      isBnk = .TRUE.
      CheckNoNumbers = 0
      DO i=stdL, 1, -1
         IF (isBnk) THEN
            IF (sTmp(i:i) .NE. " ") THEN
               CheckNoNumbers = CheckNoNumbers + 1
               isBnk = .FALSE.
            END IF
         ELSE
            IF (sTmp(i:i) .EQ. " ") THEN
               isBnk = .TRUE.
            END IF
         END IF
      END DO

      RETURN
      END FUNCTION CheckNoNumbers

!####################################################################
!     Pushing and pulling values to stack variables
      PURE SUBROUTINE PUSHSTACKS(stk, iVal)

      IMPLICIT NONE

      TYPE(stackType), INTENT(INOUT) :: stk
      INTEGER, INTENT(IN) :: iVal

      INTEGER, ALLOCATABLE :: tmp(:)

      IF (stk%maxN .EQ. 0) THEN
!     This is a new stack
         stk%maxN = 8
         ALLOCATE(stk%v(stk%maxN))
         stk%n    = 1
         stk%v(1) = iVal
      ELSE
         IF (stk%maxN .LE. stk%n) THEN
!     Stack size is too small
            ALLOCATE(tmp(stk%maxN))
            tmp = stk%v
            DEALLOCATE(stk%v)
            stk%maxN = 4*stk%maxN
            ALLOCATE(stk%v(stk%maxN))
            stk%v(1:stk%n) = tmp
         END IF
         stk%n = stk%n + 1
         stk%v(stk%n) = iVal
      END IF

      RETURN
      END SUBROUTINE PUSHSTACKS
!--------------------------------------------------------------------
      PURE SUBROUTINE PUSHSTACKV(stk, iVal)

      IMPLICIT NONE

      TYPE(stackType), INTENT(INOUT) :: stk
      INTEGER, INTENT(IN) :: iVal(:)

      INTEGER i

      DO i=1, SIZE(iVal)
         CALL PUSHSTACKS(stk, iVal(i))
      END DO

      RETURN
      END SUBROUTINE PUSHSTACKV
!--------------------------------------------------------------------
      FUNCTION PULLSTACK(stk, iVal) RESULT(flag)

      IMPLICIT NONE

      TYPE(stackType), INTENT(INOUT) :: stk
      INTEGER, INTENT(OUT) :: iVal
      LOGICAL flag

      IF (stk%n .EQ. 0) THEN
         iVal = 0
         flag = .FALSE.
      ELSE
         iVal  = stk%v(stk%n)
         stk%n = stk%n - 1
         flag  = .TRUE.
      END IF

      RETURN
      END FUNCTION PULLSTACK

!####################################################################
!     Pushing and pulling values to queue variables
      PURE SUBROUTINE ENQUEUES(que, iVal)

      IMPLICIT NONE

      TYPE(queueType), INTENT(INOUT) :: que
      INTEGER, INTENT(IN) :: iVal

      INTEGER, ALLOCATABLE :: tmp(:)

      IF (que%maxN .EQ. 0) THEN
         que%n = 1
         que%maxN = 8
         ALLOCATE(que%v(que%maxN))
         que%v(1) = iVal
      ELSE
         IF (que%maxN .LE. que%n) THEN
!     Stack size is too small
            ALLOCATE(tmp(que%maxN))
            tmp = que%v
            DEALLOCATE(que%v)
            que%maxN = 4*que%maxN
            ALLOCATE(que%v(que%maxN))
            que%v(1:que%n) = tmp
         END IF
         que%n = que%n + 1
         que%v(que%n) = iVal
      END IF

      RETURN
      END SUBROUTINE ENQUEUES
!--------------------------------------------------------------------
      PURE SUBROUTINE ENQUEUEV(que, iVal)

      IMPLICIT NONE

      TYPE(queueType), INTENT(INOUT) :: que
      INTEGER, INTENT(IN) :: iVal(:)

      INTEGER i

      DO i=1, SIZE(iVal)
         CALL ENQUEUES(que, iVal(i))
      END DO

      RETURN
      END SUBROUTINE ENQUEUEV
!--------------------------------------------------------------------
      FUNCTION DEQUEUE(que, iVal) RESULT(flag)

      IMPLICIT NONE

      TYPE(queueType), INTENT(INOUT) :: que
      INTEGER, INTENT(OUT) :: iVal

      INTEGER i
      LOGICAL flag

      IF (que%n .EQ. 0) THEN
         iVal = 0
         flag = .FALSE.
      ELSE
         iVal  = que%v(1)
         DO i=2, que%n
            que%v(i-1) = que%v(i)
         END DO
         que%n = que%n - 1
         flag  = .TRUE.
      END IF

      RETURN
      END FUNCTION DEQUEUE

!####################################################################
!     These functions produce string from numbers. Following is for
!     doubles
      PURE FUNCTION DTSTR(dVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, PARAMETER :: l = 8

      REAL(KIND=8), INTENT(IN) :: dVal
      CHARACTER(LEN=l) string

      string = NDTSTR(dVal,l)

      RETURN
      END FUNCTION DTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION NDTSTR(dVal,l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l
      REAL(KIND=8), INTENT(IN) :: dVal
      CHARACTER(LEN=l) string

      INTEGER ex, pos, i, j, abex, k, exex
      REAL(KIND=8) absn

      IF (l .EQ. 0) RETURN

      IF (dVal .NE. dVal) THEN
         IF (l .GE. 3) THEN
            string = ""
            string(l-2:l) = "NaN"
         ELSE
            string = "NaN"(1:l)
         END IF
         RETURN
      END IF

      absn = ABS(dVal)
      IF (absn .GT. HUGE(absn)) THEN
         IF (l .GE. 8) THEN
            string = ""
            string(l-7:l) = "Infinity"
         ELSE
            string = "Infinity"(1:l)
         END IF
         RETURN
      END IF

      IF (absn .LT. TINY(absn)) THEN
         pos = 1
         IF (l .GE. 3) THEN
            string(1:3) = "0.0"
            pos = 4
         END IF
         DO i=pos, l
            string(i:i) = "0"
         END DO
         RETURN
      END IF

!     This is the exponent of the number
      IF (absn .NE. 0D0) THEN
         ex = FLOOR(LOG10(absn))
      ELSE
         ex = 0
      END IF
      abex = ABS(ex)
!     How many digits exponent has
      IF (ex .NE. 0) THEN
         exex = FLOOR(LOG10(REAL(abex,8))) + 1
      ELSE
         exex = 0
      END IF

!     Number of digits of exponents and at least one for number
      i = exex + 1
!     Negative sign of the number
      IF (dVal .LT. 0D0) i = i + 1
!     For negative sign on exponent
      IF (ex .LT. 0) i = i + 1
!     That is for 'E'
      IF (ex .NE. 0) i = i + 1
!     If the sum of above is greater than the availale slots, it can
!     not be represnted with in the length. So we replace it with stars
      IF (i .GT. l) THEN
         DO i=1, l
            string(i:i) = "*"
         END DO
         RETURN
      END IF

!     Constructing the exponent first
!     pos is the position of the charecter to be written into string
      IF (ex .NE. 0) THEN
!     Writing the digits first
         DO pos=l, l-exex+1,-1
            k = MODULO(abex,10) + 1
            string(pos:pos) = "0123456789"(k:k)
            abex = abex/10
         END DO
!     Then negative sign if necessary
         IF (ex .LT. 0) THEN
            string(pos:pos) = "-"
            pos = pos - 1
         END IF
         string(pos:pos) = "E"
         pos = pos - 1
      ELSE
         pos = l
      END IF

!     l-i is the useful length remaining, beside the first number
      IF (l-i .GE. 1) THEN
         absn = absn*(1D1**(-ex/2))
         absn = absn*(1D1**(l-i-1-ex+ex/2))
         j = pos
         DO pos=j,j-l+i+2,-1
            k = FLOOR(MODULO(absn,1D1)) + 1
            string(pos:pos) = "0123456789"(k:k)
            absn = absn/1D1
         END DO
         string(pos:pos) = "."
         pos = pos - 1
         k = FLOOR(MODULO(absn,1D1)) + 1
         string(pos:pos) = "0123456789"(k:k)
      ELSE ! l-i .EQ. 0
         absn = absn*(1D1**(l-i-ex))
         k = FLOOR(MODULO(absn,1D1)) + 1
         string(pos:pos) = "0123456789"(k:k)
      END IF
      IF (dVal .LT. 0D0) string(1:1) = "-"

      RETURN
      END FUNCTION NDTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION VDTSTR(dVal) RESULT(string)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: dVal(:)
      CHARACTER(LEN=9*SIZE(dVal)) string

      INTEGER i, n

      n = SIZE(dVal)
      string = ""
      DO i=1, n
         string = TRIM(string)//" "//NDTSTR(dVal(i),8)
      END DO

      RETURN
      END FUNCTION VDTSTR
!--------------------------------------------------------------------
!     This is for real numbers
      PURE FUNCTION RTSTR(rVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, PARAMETER :: l = 7

      REAL, INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal,8),l)

      RETURN
      END FUNCTION RTSTR
!--------------------------------------------------------------------
!     This is for real numbers with specified length
      PURE FUNCTION NRTSTR(rVal,l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l
      REAL, INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal,8), l)

      RETURN
      END FUNCTION NRTSTR
!--------------------------------------------------------------------
!     This is for integers
      PURE FUNCTION ITSTR(iVal) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal
      INTEGER n
      CHARACTER(LEN=2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) string

      INTEGER absn, j, k, is

      absn = ABS(iVal)
      IF (absn .EQ. iVal) THEN
         is = 1
      ELSE
         is = 2
         string(1:1) = "-"
      END IF
      DO j=LEN(string),is,-1
         k = MODULO(absn,10) + 1
         string(j:j) = "0123456789"(k:k)
         absn = absn/10
      END DO

      RETURN
      END FUNCTION ITSTR
!--------------------------------------------------------------------
!     Similar to last one, but with minimum length of l
      PURE FUNCTION NITSTR(iVal, l) RESULT(string)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal, l
      INTEGER n
      CHARACTER(LEN=MAX(2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/)),l)) string

      string = ITSTR(iVal)
      string = ADJUSTR(string)

      RETURN
      END FUNCTION NITSTR

!####################################################################
!     Produces a color
      PURE FUNCTION CLR(iStr,clId) RESULT(oStr)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: iStr
      INTEGER, INTENT(IN), OPTIONAL :: clId
      CHARACTER(LEN=LEN(TRIM(iStr))+9) oStr
!     Colors are: 1: White, 2: Red, 3: Green, 4: Yellow, 5: Blue,
!     6: Magenta, 7: Cyan
      CHARACTER(LEN=2),PARAMETER :: clCdL(7) = (/"29","31","32","33",
     2   "34","35","36"/)

      CHARACTER(LEN=2) clCd

      IF (PRESENT(clId)) THEN
         IF (clId .LE. 0) THEN
            clCd = clCdL(1)
         ELSE
            clCd = clCdL(1+MOD(clId-1,SIZE(clCdL)))
         END IF
      ELSE
         clCd = clCdL(2)
      END IF

      oStr = CHAR(27)//'['//clCd//'m'//TRIM(iStr)//CHAR(27)//'[0m'

      RETURN
      END FUNCTION CLR
!--------------------------------------------------------------------
!     This is for removing color
      PURE FUNCTION RMCLR(iStr) RESULT(oStr)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: iStr
      CHARACTER(LEN=LEN(iStr)) oStr

      INTEGER i, j, l

!     Searching and removing all the color codes
      l    = LEN(TRIM(iStr))
      i    = 0
      j    = 0
      oStr = ""
      DO
         i = i + 1
         IF (i .GT. l-3) EXIT
         IF (iStr(i:i+1) .EQ. CHAR(27)//"[") THEN
            IF (iStr(i+3:i+3) .EQ. "m") THEN
               i = i + 3
               CYCLE
            ELSE IF (i .GT. l-4) THEN
               EXIT
            ELSE IF (iStr(i+4:i+4) .EQ. "m") THEN
               i = i + 4
               CYCLE
            END IF
         END IF
         j = j + 1
         oStr(j:j) = iStr(i:i)
      END DO
      j = j + 1
      oStr(j:j+l-i) = iStr(i:l)

      RETURN
      END FUNCTION RMCLR
!####################################################################
      FUNCTION CNCSL(sVal,lVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      LOGICAL, INTENT(IN) :: lVal

      CHARACTER(LEN=LEN(sVal) + 10) CNCSL

      IF (lVal) THEN
         CNCSL = sVal//CLR("T",3)
      ELSE
         CNCSL = sVal//CLR("F",3)
      END IF

      RETURN
      END FUNCTION CNCSL
!--------------------------------------------------------------------
!     Attaches strings and integer
      FUNCTION CNCSI(sVal,iVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      INTEGER, INTENT(IN) :: iVal

      INTEGER n
      CHARACTER(LEN=LEN(sVal) + 11 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCSI

      CNCSI = sVal//CLR(STR(iVal),3)

      RETURN
      END FUNCTION CNCSI
!--------------------------------------------------------------------
      FUNCTION CNCIS(iVal,sVal)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iVal
      CHARACTER(LEN=*), INTENT(IN) :: sVal

      INTEGER n
      CHARACTER(LEN=LEN(sVal) + 2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCIS

      CNCIS = STR(iVal)//sVal

      RETURN
      END FUNCTION CNCIS
!--------------------------------------------------------------------
      FUNCTION CNCSR(sVal,rVal)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: sVal
      REAL(KIND=8), INTENT(IN) :: rVal

      CHARACTER(LEN=LEN(sVal) + 17) CNCSR

      CNCSR = sVal//CLR(STR(rVal),3)

      RETURN
      END FUNCTION CNCSR

!####################################################################
!     This routine compares two doubles and returns .TRUE. if their
!     difference is less than "eps". If one argument is mising, it
!     will compare against 0D0

      PURE FUNCTION ISZERO(ia, ib)

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: ia
      REAL(KIND=8), INTENT(IN), OPTIONAL :: ib
      LOGICAL ISZERO

      REAL(KIND=8) a, b, tmp, nrm

!     Absolute values are calculated and I make sure "a" is bigger
      a = ABS(ia)
      b = 0D0
      IF (PRESENT(ib)) b = ABS(ib)

      IF (ABS(b) .GT. ABS(a)) THEN
         tmp = a
         a   = b
         b   = tmp
      END IF
      nrm = MAX(a,eps)

      ISZERO = .FALSE.
      IF ((a-b)/nrm .LT. 1D1*eps) ISZERO = .TRUE.

      RETURN
      END FUNCTION ISZERO

!####################################################################

      FUNCTION CPUT()

      IMPLICIT NONE

      INTEGER timeArray(8), i
      INTEGER, PARAMETER::nD(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

      REAL(KIND=8) CPUT

      CALL DATE_AND_TIME (VALUES=timeArray)
!     Year and Month
      timeArray(3) = timeArray(3) + (timeArray(1) - 2010)*365
      DO i=1, timeArray(2) - 1
         timeArray(3) = timeArray(3) + nD(i)
      END DO
!     In order: day, hr, min, sec, msec
      CPUT = timeArray(3)*8.64D4 + timeArray(5)*3.6D3 +
     2   timeArray(6)*6D1 + timeArray(7)*1D0 + timeArray(8)*1D-3

      RETURN
      END FUNCTION CPUT

!####################################################################

      PURE FUNCTION CONVI(s) RESULT(r)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: s(:)
      INTEGER r(SIZE(s))

      INTEGER i, n

      n = SIZE(r)
      DO i=1, n
         CALL MVBITS(s(i),24,8,r(i),0)
         CALL MVBITS(s(i),16,8,r(i),8)
         CALL MVBITS(s(i), 8,8,r(i),16)
         CALL MVBITS(s(i), 0,8,r(i),24)
      END DO

      RETURN
      END FUNCTION CONVI
!--------------------------------------------------------------------
      PURE FUNCTION CONVR(s) RESULT(r)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: s(:)
      REAL(KIND=8) r(SIZE(s))

      INTEGER i, n
      INTEGER(KIND=8) iS, iR

      iR = 0
      n  = SIZE(r)
      DO i=1, n
         iS   = TRANSFER(s(i),iS)
         CALL MVBITS(iS,56,8,iR, 0)
         CALL MVBITS(iS,48,8,iR, 8)
         CALL MVBITS(iS,40,8,iR,16)
         CALL MVBITS(iS,32,8,iR,24)
         CALL MVBITS(iS,24,8,iR,32)
         CALL MVBITS(iS,16,8,iR,40)
         CALL MVBITS(iS, 8,8,iR,48)
         CALL MVBITS(iS, 0,8,iR,56)
         r(i) = TRANSFER(iR,0D0)
      END DO

      RETURN
      END FUNCTION CONVR
!--------------------------------------------------------------------
      SUBROUTINE RSEED(s)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: s

      INTEGER n
      INTEGER, ALLOCATABLE :: seed(:)
      REAL, ALLOCATABLE :: seedR(:)

      CALL RANDOM_SEED(SIZE=n)
      ALLOCATE(seed(n), seedR(n))
      CALL RANDOM_NUMBER(seedR)
      seed = s + INT(n*seedR)
      CALL RANDOM_SEED(PUT=seed)

      RETURN
      END SUBROUTINE RSEED
!--------------------------------------------------------------------
      FUNCTION SEARCHARG(arg, exist) RESULT(r)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: arg
      LOGICAL, INTENT(OUT), OPTIONAL :: exist
      CHARACTER(LEN=stdL) r

      INTEGER nArg, i
      CHARACTER(LEN=stdL) tmp

      nArg = COMMAND_ARGUMENT_COUNT()
      IF (nArg .EQ. 0) STOP "No argument was specified"
      IF (PRESENT(arg)) THEN
         DO i=1, nArg
            CALL GETARG(i, tmp)
            IF (tmp.EQ.arg .OR. tmp.EQ."-"//arg) EXIT
         END DO
         IF (i .GT. nArg) THEN
            IF (PRESENT(exist)) exist = .FALSE.
            r = ""
         ELSE
            IF (PRESENT(exist)) exist = .TRUE.
            IF (i .LT. nArg) THEN
               CALL GETARG(i+1, r)
            ELSE
               r = ""
            END IF
         END IF
      ELSE
         CALL GETARG(nArg, r)
      END IF
      r = ADJUSTL(r)

      RETURN
      END FUNCTION SEARCHARG
!--------------------------------------------------------------------
      SUBROUTINE parseString(strng, toks, ntoks)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(IN) :: strng
      CHARACTER(len=*), DIMENSION(1024), INTENT(OUT) :: toks
      INTEGER, INTENT(out) :: ntoks

      CHARACTER(LEN=stdL) :: dlmtr, token

      dlmtr = ''
      token = ''

      dlmtr = '< (,=")>'
      ntoks = 1
      toks(1) = STRTOK(TRIM(strng), TRIM(dlmtr))
      toks(1) = ADJUSTL(toks(1))
      DO
         token = STRTOK(eol, TRIM(dlmtr))
         IF (token .NE. eol) THEN
            ntoks = ntoks+1
            toks(ntoks) = ADJUSTL(token)
         ELSE
            EXIT
         END IF
      END DO

      RETURN
      END SUBROUTINE parseString
!--------------------------------------------------------------------
      CHARACTER(LEN=stdL) FUNCTION STRTOK(strng, dlms)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: strng
      CHARACTER(len=*), INTENT(in) :: dlms

      INTEGER :: ist, iend
      INTEGER, SAVE :: ist0, slen
      CHARACTER(LEN=stdL), SAVE :: str0

      IF (strng(1:1) .NE. eol) then
         ist0 = 1
         str0 = strng
         slen = LEN(TRIM(str0))
      END IF

      ist = ist0
      DO
         IF (ist.LE.slen .AND. INDEX(dlms,str0(ist:ist)).NE.0) THEN
            ist = ist + 1
         ELSE
            EXIT
         END IF
      END DO

      IF (ist .GT. slen) THEN
         STRTOK = eol
         RETURN
      END IF

      iend = ist
      DO
         IF (iend.LE.slen .AND. INDEX(dlms,str0(iend:iend)).EQ.0) THEN
            iend = iend + 1
         ELSE
            EXIT
         END IF
      END DO

      STRTOK = str0(ist:iend-1)
      ist0 = iend + 1

      RETURN
      END FUNCTION STRTOK
!--------------------------------------------------------------------

      END MODULE UTILMOD

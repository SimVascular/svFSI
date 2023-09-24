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
!     string manupulation, file handeling, and norm.
!
!--------------------------------------------------------------------

      MODULE UTILMOD
      USE TYPEMOD
      IMPLICIT NONE

!     This is the standard length for all the strings used in this code
      INTEGER(KIND=IKIND), PARAMETER :: stdL = 400
!     This is the standard double number precision
      INTEGER(KIND=IKIND), PARAMETER :: dblPr = 14
!     This is the standard INTEGER(KIND=IKIND) number precision
      INTEGER(KIND=IKIND), PARAMETER :: intPr = 10

!     This is for writing binary files
      CHARACTER, PARAMETER :: eol = CHAR(10)

!     \pi value
      REAL(KIND=RKIND), PARAMETER :: pi = 3.1415926535897932384626_RKIND
!     minimum value for a double number that is not considered zero
      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

!     This is a general variable that will be used as a stack
      TYPE stackType
!        Maximum length of the stack
         INTEGER(KIND=IKIND) :: maxN = 0
!        Current size of stack
         INTEGER(KIND=IKIND) :: n = 0
!        Values inside stack
         INTEGER(KIND=IKIND), ALLOCATABLE :: v(:)
      END TYPE stackType

      TYPE queueType
         INTEGER(KIND=IKIND) :: n = 0
         INTEGER(KIND=IKIND) :: maxN = 0
         INTEGER(KIND=IKIND), ALLOCATABLE :: v(:)
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
     2      ITSTR, NITSTR, LLSTR
      END INTERFACE STR

      INTERFACE OPERATOR(//)
         MODULE PROCEDURE :: CNCSL, CNCSI, CNCIS, CNCSR
      END INTERFACE OPERATOR(//)

      INTERFACE SWAP
         MODULE PROCEDURE :: SWAPI, SWAPR
      END INTERFACE SWAP

      CONTAINS

!####################################################################
      FUNCTION OPENFILE(f) RESULT(fid)
      IMPLICIT NONE
      CLASS(fileType), INTENT(INOUT) :: f
      INTEGER(KIND=IKIND) fid

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
!     Computes U.U or U.V
      PURE FUNCTION NORMS(U, V)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: U(:)
      REAL(KIND=RKIND), INTENT(IN), OPTIONAL :: V(:)
      REAL(KIND=RKIND) NORMS

      INTEGER(KIND=IKIND) i, n

      n = SIZE(U)
      NORMS = 0._RKIND
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
      REAL(KIND=RKIND), INTENT(IN) :: U(:,:)
      REAL(KIND=RKIND), INTENT(IN), OPTIONAL :: V(:,:)
      REAL(KIND=RKIND) NORMV

      INTEGER(KIND=IKIND) m, n, i, j

      m = SIZE(U,1)
      n = SIZE(U,2)
      NORMV = 0._RKIND
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
      REAL(KIND=RKIND), INTENT(IN) :: V(:,:)
      REAL(KIND=RKIND) U(SIZE(V,1))

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

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(str)
         IF (str(i:i) .NE. " " .AND. str(i:i) .NE. "  ") EXIT
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
!     removes the read part from "strg". Valid charecters to read are
!     defined in "valS"
      SUBROUTINE GETSTR(strg, rVal, valS)
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(INOUT) :: strg
      CHARACTER(LEN=stdL), INTENT(OUT) :: rVal
      CHARACTER(LEN=*), INTENT(IN) :: valS

      LOGICAL flag
      INTEGER(KIND=IKIND) i, j, s

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
      INTEGER(KIND=IKIND8), INTENT(OUT) :: r

      CHARACTER(LEN=*), PARAMETER :: valS="+-1234567890ABCDEFabcdefZz"
      CHARACTER(LEN=stdL) tmp

      CALL GETSTR(str, tmp, valS)
      READ(tmp,"(Z8)") r

      RETURN
      END SUBROUTINE GETHEX
!--------------------------------------------------------------------
!     This reads an INTEGER from the begining of "str" and returns
!     the results in "r" and removes the readed part from "str"
      SUBROUTINE GETINT(str, r)
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(INOUT) :: str
      INTEGER(KIND=IKIND), INTENT(OUT) :: r

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
      REAL(KIND=RKIND), INTENT(IN) :: u
      INTEGER(KIND=IKIND) SGN

      IF (ISZERO(u)) THEN
         SGN = 0
      ELSE IF (u .GT. 0._RKIND) THEN
         SGN = 1
      ELSE
         SGN = -1
      END IF

      RETURN
      END FUNCTION SGN
!####################################################################
!     Returning number of data/words in a string
      PURE FUNCTION CheckNoNumbers(sTmp)
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: sTmp

      LOGICAL isBnk
      INTEGER(KIND=IKIND) CheckNoNumbers, i

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal

      INTEGER(KIND=IKIND), ALLOCATABLE :: tmp(:)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal(:)

      INTEGER(KIND=IKIND) i

      DO i=1, SIZE(iVal)
         CALL PUSHSTACKS(stk, iVal(i))
      END DO

      RETURN
      END SUBROUTINE PUSHSTACKV
!--------------------------------------------------------------------
      FUNCTION PULLSTACK(stk, iVal) RESULT(flag)
      IMPLICIT NONE
      TYPE(stackType), INTENT(INOUT) :: stk
      INTEGER(KIND=IKIND), INTENT(OUT) :: iVal
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
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal

      LOGICAL flag
      INTEGER(KIND=IKIND) i

      INTEGER(KIND=IKIND), ALLOCATABLE :: tmp(:)

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
!     Check if the new val to be added is already a member of the queue
         flag = .TRUE.
         DO i=1, que%n
            IF (que%v(i) .EQ. iVal) THEN
               flag = .FALSE.
               EXIT
            END IF
         END DO
         IF (.NOT.flag) RETURN
         que%n = que%n + 1
         que%v(que%n) = iVal
      END IF

      RETURN
      END SUBROUTINE ENQUEUES
!--------------------------------------------------------------------
      PURE SUBROUTINE ENQUEUEV(que, iVal)
      IMPLICIT NONE
      TYPE(queueType), INTENT(INOUT) :: que
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal(:)

      INTEGER(KIND=IKIND) i

      DO i=1, SIZE(iVal)
         CALL ENQUEUES(que, iVal(i))
      END DO

      RETURN
      END SUBROUTINE ENQUEUEV
!--------------------------------------------------------------------
      FUNCTION DEQUEUE(que, iVal) RESULT(flag)
      IMPLICIT NONE
      TYPE(queueType), INTENT(INOUT) :: que
      INTEGER(KIND=IKIND), INTENT(OUT) :: iVal

      INTEGER(KIND=IKIND) i
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
      REAL(KIND=RKIND8), INTENT(IN) :: dVal
      CHARACTER(LEN=RKIND8) string

      string = NDTSTR(dVal, RKIND8)

      RETURN
      END FUNCTION DTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION NDTSTR(dVal,l) RESULT(string)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: l
      REAL(KIND=RKIND8), INTENT(IN) :: dVal
      CHARACTER(LEN=l) string

      INTEGER(KIND=IKIND) ex, pos, i, j, abex, k, exex
      REAL(KIND=RKIND8) absn

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
      IF (absn .NE. 0._RKIND8) THEN
         ex = FLOOR(LOG10(absn), KIND=IKIND)
      ELSE
         ex = 0
      END IF
      abex = ABS(ex)
!     How many digits exponent has
      IF (ex .NE. 0) THEN
         exex = FLOOR(LOG10(REAL(abex, KIND=RKIND8)), KIND=IKIND) + 1
      ELSE
         exex = 0
      END IF

!     Number of digits of exponents and at least one for number
      i = exex + 1
!     Negative sign of the number
      IF (dVal .LT. 0._RKIND8) i = i + 1
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
         absn = absn*(10._RKIND8**(REAL(-ex/2, KIND=RKIND8)))
         absn = absn*(10._RKIND8**(REAL(l-i-1-ex+ex/2, KIND=RKIND8)))
         j = pos
         DO pos=j,j-l+i+2,-1
            k = FLOOR(MODULO(absn,10._RKIND8), KIND=IKIND) + 1
            string(pos:pos) = "0123456789"(k:k)
            absn = absn/10._RKIND8
         END DO
         string(pos:pos) = "."
         pos = pos - 1
         k = FLOOR(MODULO(absn,10._RKIND8)) + 1
         string(pos:pos) = "0123456789"(k:k)
      ELSE ! l-i .EQ. 0
         absn = absn*(10._RKIND8**(REAL(l-i-ex, KIND=RKIND)))
         k = FLOOR(MODULO(absn,10._RKIND8), KIND=IKIND) + 1
         string(pos:pos) = "0123456789"(k:k)
      END IF
      IF (dVal .LT. 0._RKIND8) string(1:1) = "-"

      RETURN
      END FUNCTION NDTSTR
!--------------------------------------------------------------------
!     Similar to last one, but with a specified length
      PURE FUNCTION VDTSTR(dVal) RESULT(string)
      IMPLICIT NONE
      REAL(KIND=RKIND8), INTENT(IN) :: dVal(:)
      CHARACTER(LEN=9*SIZE(dVal)) string

      INTEGER(KIND=IKIND) i, n

      n = SIZE(dVal)
      string = ""
      DO i=1, n
         string = TRIM(string)//" "//NDTSTR(dVal(i), RKIND8)
      END DO

      RETURN
      END FUNCTION VDTSTR
!--------------------------------------------------------------------
!     This is for real numbers
      PURE FUNCTION RTSTR(rVal) RESULT(string)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), PARAMETER :: l = 7
      REAL(KIND=RKIND4), INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal, KIND=RKIND8),l)

      RETURN
      END FUNCTION RTSTR
!--------------------------------------------------------------------
!     This is for real numbers with specified length
      PURE FUNCTION NRTSTR(rVal,l) RESULT(string)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: l
      REAL(KIND=RKIND4), INTENT(IN) :: rVal
      CHARACTER(LEN=l) string

      string = NDTSTR(REAL(rVal, KIND=RKIND8), l)

      RETURN
      END FUNCTION NRTSTR
!--------------------------------------------------------------------
!     This is for integers
      PURE FUNCTION ITSTR(iVal) RESULT(string)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal
      INTEGER(KIND=IKIND) n
      CHARACTER(LEN=2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) string

      INTEGER(KIND=IKIND) absn, j, k, is

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal, l
      INTEGER(KIND=IKIND) n
      CHARACTER(LEN=MAX(2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/)),l)) string

      string = ITSTR(iVal)
      string = ADJUSTR(string)

      RETURN
      END FUNCTION NITSTR
!--------------------------------------------------------------------
!     This is for logicals
      PURE FUNCTION LLSTR(flag) RESULT(string)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: flag
      CHARACTER(LEN=1) string

      IF (flag) THEN
         string = 'T'
      ELSE
         string = 'F'
      END IF

      RETURN
      END FUNCTION LLSTR
!####################################################################
!     Produces a color
      PURE FUNCTION CLR(iStr,clId) RESULT(oStr)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: iStr
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: clId
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

      INTEGER(KIND=IKIND) i, j, l

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal

      INTEGER(KIND=IKIND) n
      CHARACTER(LEN=LEN(sVal) + 11 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCSI

      CNCSI = sVal//CLR(STR(iVal),3)

      RETURN
      END FUNCTION CNCSI
!--------------------------------------------------------------------
      FUNCTION CNCIS(iVal,sVal)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iVal
      CHARACTER(LEN=*), INTENT(IN) :: sVal

      INTEGER(KIND=IKIND) n
      CHARACTER(LEN=LEN(sVal) + 2 - MAX(0,SIGN(1,ival)) +
     2   MAXVAL((/(MIN(ABS(iVal)/10**n,1)*n,n=1,9)/))) CNCIS

      CNCIS = STR(iVal)//sVal

      RETURN
      END FUNCTION CNCIS
!--------------------------------------------------------------------
      FUNCTION CNCSR(sVal,rVal)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: sVal
      REAL(KIND=RKIND), INTENT(IN) :: rVal

      CHARACTER(LEN=LEN(sVal) + 17) CNCSR

      CNCSR = sVal//CLR(STR(rVal),3)

      RETURN
      END FUNCTION CNCSR
!####################################################################
!     This routine compares two doubles and returns .TRUE. if their
!     difference is less than "eps". If one argument is mising, it
!     will compare against 0.0
      PURE FUNCTION ISZERO(ia, ib)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: ia
      REAL(KIND=RKIND), INTENT(IN), OPTIONAL :: ib
      LOGICAL ISZERO

      REAL(KIND=RKIND) a, b, tmp, nrm

!     Absolute values are calculated and I make sure "a" is bigger
      a = ABS(ia)
      b = 0._RKIND
      IF (PRESENT(ib)) b = ABS(ib)

      IF (ABS(b) .GT. ABS(a)) THEN
         tmp = a
         a   = b
         b   = tmp
      END IF
      nrm = MAX(a,eps)

      ISZERO = .FALSE.
      IF ((a-b)/nrm .LT. 10._RKIND*eps) ISZERO = .TRUE.

      RETURN
      END FUNCTION ISZERO
!####################################################################
      FUNCTION CPUT()
      IMPLICIT NONE
      REAL(KIND=RKIND) CPUT

      INTEGER(KIND=IKIND) ct, cr

      CALL SYSTEM_CLOCK(COUNT=ct, COUNT_RATE=cr)

      CPUT = REAL(ct,KIND=RKIND)/REAL(cr,KIND=RKIND)

      RETURN
      END FUNCTION CPUT
!####################################################################
      SUBROUTINE RSEED(s)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: s

      INTEGER(KIND=IKIND) n
      INTEGER(KIND=IKIND), ALLOCATABLE :: seed(:)
      REAL, ALLOCATABLE :: seedR(:)

      CALL RANDOM_SEED(SIZE=n)
      ALLOCATE(seed(n), seedR(n))
      CALL RANDOM_NUMBER(seedR)
      seed = s + n*INT(seedR, KIND=IKIND)
      CALL RANDOM_SEED(PUT=seed)

      RETURN
      END SUBROUTINE RSEED
!####################################################################
      FUNCTION SEARCHARG(arg, exst) RESULT(r)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: arg
      LOGICAL, INTENT(OUT), OPTIONAL :: exst
      CHARACTER(LEN=stdL) r

      INTEGER(KIND=IKIND) nArg, i
      CHARACTER(LEN=stdL) tmp

      nArg = COMMAND_ARGUMENT_COUNT()
      IF (nArg .EQ. 0) STOP "No argument was specified"
      IF (PRESENT(arg)) THEN
         DO i=1, nArg
            CALL GETARG(i, tmp)
            IF (tmp.EQ.arg .OR. tmp.EQ."-"//arg) EXIT
         END DO
         IF (i .GT. nArg) THEN
            IF (PRESENT(exst)) exst = .FALSE.
            r = ""
         ELSE
            IF (PRESENT(exst)) exst = .TRUE.
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
!####################################################################
      SUBROUTINE parseString(strng, toks, ntoks)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: strng
      CHARACTER(LEN=*), DIMENSION(250), INTENT(OUT) :: toks
      INTEGER(KIND=IKIND), INTENT(OUT) :: ntoks

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
      CHARACTER(LEN=*), INTENT(IN) :: strng
      CHARACTER(LEN=*), INTENT(IN) :: dlms

      INTEGER(KIND=IKIND) :: ist, iend
      INTEGER(KIND=IKIND), SAVE :: ist0, slen
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
!####################################################################
      SUBROUTINE TO_UPPER(strng)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: strng

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(strng)
         SELECT CASE(strng(i:i))
         CASE("a":"z")
            strng(i:i) = ACHAR(IACHAR(strng(i:i))-32)
         END SELECT
      END DO

      RETURN
      END SUBROUTINE TO_UPPER
!--------------------------------------------------------------------
      SUBROUTINE TO_LOWER(strng)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(INOUT) :: strng

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(strng)
         SELECT CASE(strng(i:i))
         CASE("A":"Z")
            strng(i:i) = ACHAR(IACHAR(strng(i:i))+32)
         END SELECT
      END DO

      RETURN
      END SUBROUTINE TO_LOWER
!####################################################################
      SUBROUTINE SWAPI(a, b)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(INOUT) :: a, b

      INTEGER(KIND=IKIND) c

      c = a
      a = b
      b = c

      RETURN
      END SUBROUTINE SWAPI
!--------------------------------------------------------------------
      SUBROUTINE SWAPR(a, b)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: a, b

      REAL(KIND=RKIND) c

      c = a
      a = b
      b = c

      RETURN
      END SUBROUTINE SWAPR
!####################################################################
      END MODULE UTILMOD
!####################################################################

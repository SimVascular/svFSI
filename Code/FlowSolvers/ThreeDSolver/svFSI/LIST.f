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
!     This module is used for I/O.
!
!--------------------------------------------------------------------

      MODULE LISTMOD
      USE CHNLMOD
      IMPLICIT NONE

      TYPE listType
         PRIVATE
!        Whether this line has been used sofar
         LOGICAL :: used = .FALSE.
!        Length of sub list
         INTEGER :: l = 0
!        Line number associated with this list
         INTEGER line
!        Command
         CHARACTER(LEN=stdL) :: kwd = 'NONE'
!        Value associated with the command
         CHARACTER(LEN=stdL) :: val = 'NONE'
!        Sublist, under current list
         TYPE(listType), POINTER :: sub(:)
         TYPE(ioType), POINTER :: io
      CONTAINS
         PROCEDURE :: PING
         PROCEDURE :: LSRCH
         PROCEDURE, PUBLIC :: srch
         PROCEDURE, PUBLIC :: check => CHECKLIST
         PROCEDURE :: GFLL
         PROCEDURE :: GFLIS
         PROCEDURE :: GFLIV
         PROCEDURE :: GFLRS
         PROCEDURE :: GFLRV
         PROCEDURE :: GFLS
         PROCEDURE :: GFLF
         GENERIC, PUBLIC :: get => GFLL, GFLIS, GFLIV, GFLRS, GFLRV,
     2      GFLS, GFLF
      END TYPE listType

      INTERFACE listType
         MODULE PROCEDURE NEWLIST
      END INTERFACE listType

      CONTAINS

!####################################################################
      FUNCTION NEWLIST(fileName, inIO)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: fileName
      TYPE(ioType), TARGET, INTENT(IN) :: inIO
      TYPE(listType) NEWLIST

      INTEGER, PARAMETER :: maxL = 100

      LOGICAL flag
      INTEGER s, e, l, i, j, A, fNl, lInd(stdL), lvl, fid
      CHARACTER(LEN=stdL) ctmp, sTmp

      INTEGER, ALLOCATABLE :: fInd(:), fLineN(:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: fCon(:)

      NEWLIST%io   => inIO
      NEWLIST%kwd  = fileName
      NEWLIST%used = .TRUE.

!     These are the default values used for the first few lines below
      fid = 1
      INQUIRE (FILE=TRIM(fileName), EXIST=flag)
      IF (.NOT.flag) THEN
         NEWLIST%io%e = "File <"//TRIM(fileName)//"> does not exist"//
     2      " or can not be opened"
      END IF
      OPEN(fid, FILE=fileName, STATUS="OLD")

!     Reading the entire ini file
      fNl = 0
      DO
         READ(fid,"(A)",IOSTAT=i) sTmp
         sTmp = ADJUSTC(sTmp)
         IF (i .GT. 0) THEN
            NEWLIST%io%e = "While reading: "//fileName
         ELSE IF (i .LT. 0) THEN
            EXIT
         END IF
         IF (sTmp .EQ. "End") EXIT
         fNl = fNl + 1
      END DO

!     This for the case a curly brackets is at the middle of a line
      e = 100
      ALLOCATE(fCon(fNl+e), fInd(fNl+e), fLineN(fNl+e))
      fCon = ""
      REWIND(fid)
      A = 0
      lvl = 1
      DO i=1, fNl
         READ(fid,"(A)") sTmp
         sTmp = ADJUSTC(sTmp)
         l = LEN(TRIM(sTmp))
!     Empty and commented lines are ignored
         IF (l.EQ.0 .OR. sTmp(1:1).EQ."#") CYCLE

         e = 0
         ctmp = ""
!     Going one level up or down if "{" or "}" are found
         DO j=1, l
            IF (sTmp(j:j) .EQ. "{") THEN
               lvl = lvl + 1
            ELSE IF (sTmp(j:j) .EQ. "}") THEN
               lvl = lvl - 1
            ELSE
               e = e + 1
               lInd(e) = lvl
               ctmp(e:e) = sTmp(j:j)
            END IF
         END DO
         l = e
         sTmp = ctmp
         IF (l .EQ. 0) CYCLE
!     Cutting a line in pieces if the level is changing in a particular
!     character in the line
         e = lInd(1)
         s = 1
         DO j=1, l
            IF (lInd(j) .NE. e) THEN
               A = A + 1
               fCon(A) = ADJUSTC(sTmp(s:j-1))
               fInd(A) = e
               fLineN(A) = i
               s = j
               e = lInd(j)
            END IF
         END DO
         A = A + 1
         fCon(A) = ADJUSTC(sTmp(s:l))
         fInd(A) = e
         fLineN(A) = i
      END DO
      IF (ANY(fInd(1:A) .LT. 1)) THEN
         NEWLIST%io%e = "Too many } in "//fileName
      END IF
      IF (lvl .NE. 1) NEWLIST%io%e = "Too many { in "//fileName
      CLOSE(fid)
      fNl = A

!     transforming the entire thing into a list
      CALL SETSUBLIST(1, fNl, NEWLIST)

      RETURN
      CONTAINS
!--------------------------------------------------------------------
!     This SUBROUTINE sets the length of a list, allocates and sets
!     value to list%sub based on the lines that are in the same
!     level as first line. To call this, you need to provide the list
!     structure, first line number (s) and last line number (e)
      RECURSIVE SUBROUTINE SETSUBLIST(s, e, list)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: s, e
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL flag
      INTEGER lvl, i, j, k, l
      CHARACTER(LEN=stdL) sTmp

!     Counting the numebr of lines to be allocated
      lvl = fInd(s)
      DO i=s, e
         IF (fInd(i) .EQ. lvl) list%l = list%l + 1
      END DO

      ALLOCATE(list%sub(list%l))
      j = 0
!     We assume this list does not have a sublist.
      flag = .FALSE.
      DO i=s, e
         IF (fInd(i) .EQ. lvl) THEN
            IF (flag) THEN
!     This is the end of sublist
               flag = .FALSE.
               CALL SETSUBLIST(k, i-1, list%sub(j))
            END IF
            j = j + 1
            sTmp = fCon(i)
            DO l=1, LEN(TRIM(sTmp))
               IF (sTmp(l:l) .EQ. ":") EXIT
            END DO
            list%sub(j)%kwd  = sTmp(:l-1)
            list%sub(j)%val  = ADJUSTC(sTmp(l+1:))
            list%sub(j)%line = fLineN(i)
            list%sub(j)%io   => inIO
         ELSE
!     There is a sublist!
            IF (.NOT. flag) k = i
            flag = .TRUE.
         END IF
      END DO
!     This is for the case that index of last line, fInd(e), is
!     different from lvl
      IF (flag) CALL SETSUBLIST(k, i-1, list%sub(j))

      RETURN
      END SUBROUTINE SETSUBLIST
      END FUNCTION NEWLIST

!####################################################################
!     This function search through list for "cmnd" and returns the line
!     number that is found. The index of the searched line can be
!     inputted by iInd. If iInd is provided, a line is not found, error
!     is thrown
      FUNCTION LSRCH(list, cmnd, iInd)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: iInd
      INTEGER LSRCH

      INTEGER i, n, ind

      ind = 1
      IF (PRESENT(iInd)) ind = iInd

!     If the line is not found, -1 is returned
      LSRCH = -1
      n = 0
      DO i=1, list%l
         IF (list%sub(i)%kwd .EQ. cmnd) THEN
            n = n + 1
            IF (n .EQ. ind) THEN
               list%sub(i)%used = .TRUE.
               LSRCH = i
               RETURN
            END IF
         END IF
      END DO
      IF (PRESENT(iInd)) THEN
         list%io%e = TRIM(list%ping(cmnd))//" Command not found"
      END IF

      RETURN
      END FUNCTION LSRCH
!--------------------------------------------------------------------
      FUNCTION SRCH(list, cmnd, ll, ul)

      IMPLICIT NONE

      CLASS(listType), INTENT(IN) :: list
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ll, ul
      INTEGER SRCH

      INTEGER i

      SRCH = 0
      DO i=1, list%l
         IF (list%sub(i)%kwd .EQ. cmnd) SRCH = SRCH + 1
      END DO

      IF (PRESENT(ll)) THEN
         IF (SRCH .LT. ll) list%io%e = TRIM(list%ping(cmnd))//
     2      " COMMAND must be repeated more/equal than "//ll//" time/s"
      END IF
      IF (PRESENT(ul)) THEN
         IF (SRCH .GT. ul) list%io%e = TRIM(list%ping(cmnd))//
     2      " COMMAND must be repeated less/equal than "//ul//" time/s"
      END IF

      RETURN
      END FUNCTION SRCH

!####################################################################
      FUNCTION PING(list, cmnd, subL)

      IMPLICIT NONE

      CLASS(listType), INTENT(IN) :: list
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      TYPE(listType), OPTIONAL :: subL
      CHARACTER(LEN=stdL) PING

      IF (PRESENT(subL)) THEN
         PING = "At LINE "//subL%line//", searched COMMAND <"//
     2      TRIM(cmnd)//"> under <"//TRIM(list%kwd)//": "//
     3      TRIM(list%val)//"> ::"
      ELSE
         PING = "Near LINE "//list%line//" at <"//TRIM(list%kwd)//
     2      ": "//TRIM(list%val)//">, searched COMMAND <"//
     3      TRIM(cmnd)//"> ::"
      END IF

      RETURN
      END FUNCTION PING

!####################################################################
!     This routine checks all the lines of the list to make sure they
!     have been used
      RECURSIVE SUBROUTINE CHECKLIST(list)

      IMPLICIT NONE

      CLASS(listType), INTENT(IN) :: list

      INTEGER i

      IF (.NOT.list%used) list%io%e = "Near line "//list%line//
     2   " keyword <"//TRIM(list%kwd)//"> is not recognized"

      DO i=1, list%l
         CALL CHECKLIST(list%sub(i))
      END DO

      RETURN
      END SUBROUTINE CHECKLIST

!####################################################################
!     This is for reading logical
      FUNCTION GFLL(list, lVal, cmnd, ind)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      LOGICAL, INTENT(OUT) :: lVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(listType), POINTER :: GFLL

      CHARACTER(*), PARAMETER :: true(5)= (/"t   ", "T   ", "1   ",
     2   "True", "true"/), false(5) = (/"f    ", "F    ", "0    ",
     3   "False", "false"/)

      INTEGER i
      CHARACTER(LEN=stdL) c

      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
         GFLL => list%GFLS(c, cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLL => NULL()
            RETURN
         END IF
         GFLL => list%GFLS(c, cmnd, 1)
      END IF

      lVal = .TRUE.
      IF (ANY(c .EQ. true)) THEN
         list%io%d = TRIM(list%ping(cmnd,GFLL))//" Read TRUE"
      ELSE IF (ANY(c .EQ. false)) THEN
         lVal = .FALSE.
         list%io%d = TRIM(list%ping(cmnd,GFLL))//" Read FALSE"
      ELSE
         list%io%e = TRIM(list%ping(cmnd,GFLL))//" Reading error"
      END IF

      RETURN
      END FUNCTION GFLL
!--------------------------------------------------------------------
!     This function return the value of listPtr at a particular line
!     (ind) if ll/lb or ul/ub are specified, the value is checked
!     to be in that specified range. This is for reading integers.
      FUNCTION GFLIS(list, iVal, cmnd, ind, ll, ul)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      INTEGER, INTENT(OUT) :: iVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind, ll, ul
      TYPE(listType), POINTER :: GFLIS

      INTEGER i, ioS

      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLIS => NULL()
            RETURN
         END IF
      END IF
      GFLIS => list%sub(i)

      READ(GFLIS%val, *, IOSTAT=ioS) iVal
      IF (ioS .NE. 0) THEN
         list%io%e = TRIM(list%ping(cmnd,GFLIS))//" Reading error"
      END IF
      IF (PRESENT(ll)) THEN
         IF (iVal .LT. ll) list%io%e =
     2      TRIM(list%ping(cmnd,GFLIS))//" Lower limit is "//ll
      END IF
      IF (PRESENT(ul)) THEN
         IF (iVal .GT. ul) list%io%e =
     2      TRIM(list%ping(cmnd,GFLIS))//" Upper limit is "//ul
      END IF
      list%io%d = TRIM(list%ping(cmnd,GFLIS))//" Read integer value "//
     2   iVal

      RETURN
      END FUNCTION GFLIS
!--------------------------------------------------------------------
!     This is for reading a vector. The results are returned into v
      FUNCTION GFLIV(list, vVal, cmnd, ind)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      INTEGER, INTENT(OUT) :: vVal(:)
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(listType), POINTER :: GFLIV

      INTEGER i, ioS, n, nToks
      CHARACTER(LEN=stdL), DIMENSION(1024) :: tokList

      n = SIZE(vVal)
      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLIV => NULL()
            RETURN
         END IF
      END IF
      GFLIV => list%sub(i)

      CALL parseString(GFLIV%val, tokList, nToks)
      IF (nToks .NE. n) list%io%e = TRIM(list%ping(cmnd,GFLIV))//
     2   " Reading error"

      DO i=1, nToks
         READ(tokList(i), *, IOSTAT=ioS) vVal(i)
         IF (ioS .NE. 0) THEN
            list%io%e = TRIM(list%ping(cmnd,GFLIV))//" Reading error"
         END IF
      END DO

      DO i=1, n
         list%io%d = TRIM(list%ping(cmnd,GFLIV))//" Read vector value "
     2      //vVal(i)
      END DO

      RETURN
      END FUNCTION GFLIV
!--------------------------------------------------------------------
!     This is for reading real numbers
      FUNCTION GFLRS(list, rVal, cmnd, ind, ll, ul, lb, ub)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      REAL(KIND=8), INTENT(OUT) :: rVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      REAL(KIND=8), INTENT(IN), OPTIONAL :: ll, ul, lb, ub
      TYPE(listType), POINTER :: GFLRS

      INTEGER i, ioS

      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLRS => NULL()
            RETURN
         END IF
      END IF
      GFLRS => list%sub(i)

      READ(GFLRS%val, *, IOSTAT=ioS) rVal
      IF (ioS .NE. 0) THEN
         list%io%e = TRIM(list%ping(cmnd,GFLRS))//" Reading error"
      END IF
      IF (PRESENT(ll)) THEN
         IF (rVal .LT. ll) list%io%e =
     2      TRIM(list%ping(cmnd,GFLRS))//" Lower limit is "//ll
      END IF
      IF (PRESENT(ul)) THEN
         IF (rVal .GT. ul) list%io%e =
     2      TRIM(list%ping(cmnd,GFLRS))//" Upper limit is "//ul
      END IF
      IF (PRESENT(lb)) THEN
         IF (rVal .LE. lb) list%io%e =
     2      TRIM(list%ping(cmnd,GFLRS))//" Lower bound is "//lb
      END IF
      IF (PRESENT(ub)) THEN
         IF (rVal .GE. ub) list%io%e =
     2      TRIM(list%ping(cmnd,GFLRS))//" Upper bound is "//ub
      END IF
      list%io%d = TRIM(list%ping(cmnd,GFLRS))//" Read real value "//
     2   rVal

      RETURN
      END FUNCTION GFLRS
!--------------------------------------------------------------------
!     This is for reading a vector. The results are returned into v
      FUNCTION GFLRV(list, vVal, cmnd, ind)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      REAL(KIND=8), INTENT(OUT) :: vVal(:)
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(listType), POINTER :: GFLRV

      INTEGER i, ioS, n, nToks
      CHARACTER(LEN=stdL), DIMENSION(1024) :: tokList

      n = SIZE(vVal)
      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLRV => NULL()
            RETURN
         END IF
      END IF
      GFLRV => list%sub(i)

      CALL parseString(GFLRV%val, tokList, nToks)
      IF (nToks .LT. n) list%io%e = TRIM(list%ping(cmnd,GFLRV))//
     2   " Parse error. Expected no of values not found."

      DO i=1, nToks
         READ(tokList(i), *, IOSTAT=ioS) vVal(i)
         IF (ioS .NE. 0) THEN
            list%io%e = TRIM(list%ping(cmnd,GFLRV))//" Reading error"
         END IF
      END DO

      DO i=1, n
         list%io%d = TRIM(list%ping(cmnd,GFLRV))//" Read vector value "
     2      //vVal(i)
      END DO

      RETURN
      END FUNCTION GFLRV
!--------------------------------------------------------------------
!     This is for reading strings
      FUNCTION GFLS(list, sVal, cmnd, ind)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=stdL), INTENT(OUT) :: sVal
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(listType), POINTER :: GFLS

      INTEGER i, ioS

      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLS => NULL()
            RETURN
         END IF
      END IF
      GFLS => list%sub(i)

      READ(GFLS%val,*,IOSTAT=ioS) sVal
      IF (ioS .NE. 0) THEN
         list%io%e = TRIM(list%ping(cmnd,GFLS))//" Reading error"
      END IF
      list%io%d = TRIM(list%ping(cmnd,GFLS))//" Read char value <"
     2   //TRIM(sVal)//">"

      RETURN
      END FUNCTION GFLS
!--------------------------------------------------------------------
!     This function opens a file and returns the handle to the file
      FUNCTION GFLF(list, file, cmnd, ind)

      IMPLICIT NONE

      CLASS(listType), INTENT(INOUT) :: list
      TYPE(fileType), INTENT(INOUT) :: file
      CHARACTER(LEN=*), INTENT(IN) :: cmnd
      INTEGER, INTENT(IN), OPTIONAL :: ind
      TYPE(listType), POINTER :: GFLF

      INTEGER i
      CHARACTER(LEN=stdL) fName

      IF (PRESENT(ind)) THEN
         i = list%LSRCH(cmnd, ind)
      ELSE
         i = list%LSRCH(cmnd)
         IF (i .LT. 0) THEN
            GFLF => NULL()
            RETURN
         END IF
      END IF
      GFLF => list%sub(i)

      fName = GFLF%val
      file%fname = fName

      list%io%d = TRIM(list%ping(cmnd,GFLF))//
     2   " Opened file from path <"//TRIM(fName)//">"

      RETURN
      END FUNCTION GFLF

!####################################################################
      RECURSIVE SUBROUTINE DESTROYLIST(list)
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list
      INTEGER :: stat, i

      DO i=1, list%l
         IF (list%sub(i)%l /= 0) THEN
            CALL DESTROYLIST(list%sub(i))
         END IF
      END DO

      IF (list%l /= 0) THEN
         DEALLOCATE(list%sub, STAT = stat)
      END IF
      IF ( stat /= 0 ) STOP "*** Trouble deallocating ***"

      END SUBROUTINE DESTROYLIST
!####################################################################

      END MODULE LISTMOD

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

      MODULE CHNLMOD

      USE ISO_FORTRAN_ENV
      USE UTILMOD

      IMPLICIT NONE

!     Channels cases: standard output, error output, warning output,
!     debugging output
      INTEGER, PARAMETER :: CHNL_O = 601
      INTEGER, PARAMETER :: CHNL_E = 602
      INTEGER, PARAMETER :: CHNL_W = 603
      INTEGER, PARAMETER :: CHNL_D = 604

!     Whether to use color in printing outputs
      LOGICAL :: pClr = .TRUE.

!     A general counter for file ID
      INTEGER :: gFID = 314

!     Appended path to all files that are going to be saved
      CHARACTER(LEN=stdL) :: appPath = ""

!     Channel type, used in I/O
      TYPE chnlType
!        Whether it is open to the screen
         LOGICAL :: oTS = .FALSE.
!        Whether it is open to the file
         LOGICAL :: oTF = .FALSE.
!        Channel identifier
         INTEGER id
!        File ID
         INTEGER fId
!        File address
         CHARACTER(LEN=stdL) :: fName = "histor"
!        Channel tag
         CHARACTER(LEN=stdL) :: tag = ""
      CONTAINS
!        Creates a new channel
         PROCEDURE :: new => newChnl
!        Closes the channel
         PROCEDURE :: close => closeChnl
!        To send a string to channel
         PROCEDURE chnlAssign
         GENERIC :: ASSIGNMENT(=) => chnlAssign
      END TYPE chnlType

!     Only to group four channels, in case I rather have them as one
!     variable
      TYPE ioType
!        Standard output
         TYPE(chnlType) :: o
!        Error
         TYPE(chnlType) :: e
!        Warning
         TYPE(chnlType) :: w
!        Debugging
         TYPE(chnlType) :: d
!        Status file
         TYPE(chnlType) :: s
      CONTAINS
!        Opens all as standard channels
         PROCEDURE :: new => newIO
!        Closes the channel
         PROCEDURE :: close => closeIO
      END TYPE ioType

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE chnlAssignChnl, ioAssignIO
      END INTERFACE ASSIGNMENT(=)

      CONTAINS

!####################################################################

      SUBROUTINE newChnl(chnl,id,fName,tag,oTS,oTF)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(OUT) :: chnl
      INTEGER, INTENT(IN) :: id
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fName, tag
      LOGICAL, INTENT(IN), OPTIONAL :: oTS, oTF

      LOGICAL flag

      IF (ALL(id.NE.(/CHNL_O,CHNL_E,CHNL_W,CHNL_D/)))
     2   STOP "Non-existent channel"

      chnl%id = id
      IF (PRESENT(fName)) chnl%fName = ADJUSTL(fName)
      IF (PRESENT(tag)) THEN
         chnl%tag = ADJUSTL(tag)
         IF (chnl%tag .NE. "") chnl%fName =
     2      TRIM(chnl%fName)//"-"//chnl%tag
      END IF

!     By defult outputs are open for error and warning
      IF (PRESENT(oTS)) THEN
         chnl%oTS = oTS
      ELSE
         IF (id.EQ.CHNL_E .OR. id.EQ.CHNL_W) chnl%oTS = .TRUE.
      END IF
      IF (PRESENT(oTF)) THEN
         chnl%oTF = oTF
      ELSE
         IF (id.EQ.CHNL_E .OR. id.EQ.CHNL_W) chnl%oTF = .TRUE.
      END IF
      chnl%fName = TRIM(chnl%fName)//".dat"

      gFID     = gFID + 1
      chnl%fId = gFID
      DO
         INQUIRE(UNIT=chnl%fId, OPENED=flag)
         IF (.NOT.flag) EXIT
         chnl%fId = chnl%fId + 1
      END DO

      RETURN
      END SUBROUTINE newChnl
!--------------------------------------------------------------------
      SUBROUTINE closeChnl(chnl)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl

      LOGICAL flag

      INQUIRE(UNIT=chnl%fId, OPENED=flag)
      IF (flag) CLOSE(chnl%fId)

      RETURN
      END SUBROUTINE closeChnl
!--------------------------------------------------------------------
      SUBROUTINE chnlAssignChnl(s,r)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(OUT) :: s
      TYPE(chnlType), INTENT(IN) :: r

      s%oTS   = r%oTS
      s%oTF   = r%oTF
      s%id    = r%id
      s%fId   = r%fId
      s%fName = r%fName
      s%tag   = r%tag

      RETURN
      END SUBROUTINE chnlAssignChnl
!--------------------------------------------------------------------
      SUBROUTINE ioAssignIO(s,r)
      IMPLICIT NONE
      CLASS(ioType), INTENT(OUT) :: s
      TYPE(ioType), INTENT(IN) :: r

      s%o = r%o
      s%e = r%e
      s%w = r%w
      s%d = r%d
      s%s = r%s

      RETURN
      END SUBROUTINE ioAssignIO

!####################################################################
!     Sending strings to correct channel
      SUBROUTINE chnlAssign(chnl,sInp)
      IMPLICIT NONE
      CLASS(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: sInp

      SELECT CASE(chnl%Id)
      CASE (CHNL_O)
         CALL CHNLOUTPUT(chnl,sInp)
      CASE (CHNL_E)
         CALL CHNLERROR(chnl,sInp)
      CASE (CHNL_W)
         CALL CHNLWARNING(chnl,sInp)
      CASE (CHNL_D)
         CALL CHNLDEBUGGING(chnl,sInp)
      END SELECT

      RETURN
      END SUBROUTINE chnlAssign
!--------------------------------------------------------------------
!     This routine is used for keeping track of what is printed on the
!     screen and history file
      SUBROUTINE CHNLOUTPUT(chnl,isTmp)
      IMPLICIT NONE
      TYPE(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: isTmp

      LOGICAL flag
      CHARACTER(LEN=stdL) sTmp, fName

      IF (.NOT.chnl%oTS .AND. .NOT.chnl%oTF) RETURN

!     Searching and removing all the color codes
      IF (chnl%oTF .OR. (chnl%oTS.AND.(.NOT.pClr))) sTmp = RMCLR(isTmp)

      IF (chnl%oTF) THEN
         INQUIRE(UNIT=chnl%fId, OPENED=flag)
         IF (.NOT.flag) THEN
            IF (appPath .NE. "") CALL SYSTEM("mkdir -p "//TRIM(appPath))
            fName = TRIM(appPath)//TRIM(chnl%fName)
            INQUIRE(FILE=fName, OPENED=flag)
            IF (.NOT.flag) THEN
               OPEN(chnl%fId, FILE=fName, POSITION="APPEND")
            ELSE
               INQUIRE(FILE=fName, NUMBER=chnl%fId)
            END IF
         END IF
         WRITE(chnl%fId,"(A)") TRIM(sTmp)
         CALL FLUSH(chnl%fId)
      END IF

      IF (chnl%oTS) THEN
         IF (chnl%tag .NE. "") THEN
            IF (pClr) THEN
               WRITE(*,"(A)") CLR(TRIM(chnl%tag)//" >>",4)//" "//
     2            TRIM(isTmp)
            ELSE
               WRITE(*,"(A)") TRIM(chnl%tag)//" >> "//TRIM(sTmp)
            END IF
         ELSE
            IF (pClr) THEN
               WRITE(*,"(A)") TRIM(isTmp)
            ELSE
               WRITE(*,"(A)") TRIM(sTmp)
            END IF
         END IF
         CALL FLUSH(OUTPUT_UNIT)
      END IF

      RETURN
      END SUBROUTINE CHNLOUTPUT
!--------------------------------------------------------------------
!     This routine is used for ERROR handling
      SUBROUTINE CHNLERROR(chnl,isTmp)
      IMPLICIT NONE
      TYPE(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: isTmp

      CHARACTER(LEN=LEN(isTmp)) sTmp

      IF (chnl%oTS .OR. chnl%oTF) THEN
         sTmp = RMCLR(isTmp)

         CALL CHNLOUTPUT(chnl,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//
     2      "!!!!!!!!")
         CALL CHNLOUTPUT(chnl,"ERROR occured, see below for more"//
     2      " explanation")
         CALL CHNLOUTPUT(chnl,CLR("ERROR: "//TRIM(ADJUSTL(sTmp))))
         CALL CHNLOUTPUT(chnl,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"//
     2      "!!!!!!!!")
      END IF
      STOP "All processors are forced to stop by a fatal error"

      END SUBROUTINE CHNLERROR
!--------------------------------------------------------------------
!     This routine is used for ERROR handling
      SUBROUTINE CHNLWARNING(chnl,sTmp)
      IMPLICIT NONE
      TYPE(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: sTmp

      IF (chnl%oTS .OR. chnl%oTF) THEN
         IF (pClr) THEN
            CALL CHNLOUTPUT(chnl,CLR(" WARNING:",4)//" "//ADJUSTL(sTmp))
         ELSE
            CALL CHNLOUTPUT(chnl,"!!! WARNING: "//ADJUSTL(sTmp))
         END IF
      END IF

      RETURN
      END SUBROUTINE CHNLWARNING
!--------------------------------------------------------------------
!     This routine is used for ERROR handling
      SUBROUTINE CHNLDEBUGGING(chnl,sTmp)
      IMPLICIT NONE
      TYPE(chnlType), INTENT(INOUT) :: chnl
      CHARACTER(LEN=*), INTENT(IN) :: sTmp

      IF (chnl%oTS .OR. chnl%oTF) THEN
         CALL CHNLOUTPUT(chnl,CLR(" DEBUG-INFO:",3)//" "//ADJUSTL(sTmp))
      END IF

      RETURN
      END SUBROUTINE CHNLDEBUGGING

!####################################################################

      SUBROUTINE newIO(io, iflag)
      IMPLICIT NONE
      CLASS(ioType), INTENT(OUT) :: io
      LOGICAL, INTENT(IN), OPTIONAL :: iflag

      LOGICAL flag

      IF (PRESENT(iflag)) flag = iflag

      CALL io%o%new(CHNL_O,oTS=flag,oTF=flag)
      CALL io%e%new(CHNL_E)
      CALL io%w%new(CHNL_W)
      CALL io%d%new(CHNL_D)
      CALL io%s%new(CHNL_O,oTF=flag)

      RETURN
      END SUBROUTINE newIO
!--------------------------------------------------------------------
      SUBROUTINE closeIO(io)
      IMPLICIT NONE
      CLASS(ioType), INTENT(INOUT) :: io

      CALL io%o%close()
      CALL io%e%close()
      CALL io%w%close()
      CALL io%d%close()
      CALL io%s%close()

      RETURN
      END SUBROUTINE closeIO

      END MODULE CHNLMOD


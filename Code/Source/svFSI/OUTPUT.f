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
!     This routine contains multiple functions that are generally
!     desined to interface with user.
!
!--------------------------------------------------------------------

!     Prepares the output of svFSI to the standard output.
      SUBROUTINE OUTRESULT(timeP, co, iEq)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: co, iEq
      REAL(KIND=RKIND), INTENT(INOUT) :: timeP(3)

      CHARACTER(LEN=stdL) :: sepLine

      INTEGER(KIND=IKIND) fid, i
      REAL(KIND=RKIND) tmp, tmp1, tmp2
      CHARACTER c1, c2
      CHARACTER(LEN=stdL) sOut

      IF (cm%slv()) RETURN

      fid = 1
      tmp = CPUT()

      sepLine = REPEAT("-", 69)

      IF (co .EQ. 1) THEN
         timeP(1) = tmp - timeP(1)
         timeP(2) = 0._RKIND
         std = " "
         std = TRIM(sepLine)
         std = " Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri  "//
     2      "   lsIt   dB  %t"
         IF (nEq .EQ. 1) std = TRIM(sepLine)
         RETURN
      END IF

      IF (nEq.GT.1 .AND. iEq.EQ.1 .AND. eq(iEq)%itr.EQ.1) THEN
         std = TRIM(sepLine)
      END IF

      c1 = " "
      IF (co .EQ. 3) c1 = "s"
      timeP(3) = tmp - timeP(1)
      sOut = " "//eq(iEq)%sym//" "//STR(cTS,5)//"-"//STR(eq(iEq)%itr)//
     2   c1//" "//STR(timeP(3),6)

      IF (ISZERO(eq(iEq)%iNorm)) THEN
         tmp  = 1._RKIND
         tmp1 = 1._RKIND
         tmp2 = 1._RKIND
         i    = 0
      ELSE
         tmp  = eq(iEq)%FSILS%RI%iNorm/eq(iEq)%iNorm
         tmp1 = tmp/eq(iEq)%pNorm
         tmp2 = eq(iEq)%FSILS%RI%fNorm/eq(iEq)%FSILS%RI%iNorm
         i    = INT(20._RKIND*LOG10(tmp1), KIND=IKIND)
      END IF

      IF (i .GT. 20) THEN
         c1 = "!"; c2 = "!"
      ELSE
         c1 = "["; c2 = "]"
      END IF
      sOut = TRIM(sOut)//"  "//c1//STR(i,4)//" "//STR(tmp1,7)//" "//
     2   STR(tmp,7)//" "//STR(tmp2,7)//c2

      IF (ISZERO(timeP(3),timeP(2)))
     2   timeP(3) = (1._RKIND+eps)*timeP(2) + eps
      tmp = 100._RKIND*eq(iEq)%FSILS%RI%callD/(timeP(3) - timeP(2))
      timeP(2) = timeP(3)
      IF (ABS(tmp) .GT. 100._RKIND) tmp = 100._RKIND

      IF (eq(iEq)%FSILS%RI%suc) THEN
         c1 = "["; c2 = "]"
      ELSE
         c1 = "!"; c2 = "!"
      END IF
      sOut = TRIM(sOut)//"  "//c1//STR(eq(iEq)%FSILS%RI%itr,4)//" "//
     2   STR(NINT(eq(iEq)%FSILS%RI%dB, KIND=IKIND),4)//" "//
     3   STR(NINT(tmp, KIND=IKIND),3)//c2

      IF (nEq .GT. 1) THEN
         std = CLR(sOut,iEq)
         IF (sepOutput .AND. iEq.NE.cEq) std = "--"
      ELSE
         std = sOut
      END IF

      RETURN
      END SUBROUTINE OUTRESULT
!####################################################################
      SUBROUTINE WRITERESTART(timeP)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: timeP(3)

      INTEGER(KIND=IKIND) fid, myID
      CHARACTER(LEN=stdL) fName, tmpS

      fid  = 27
      myID = cm%tf()

      fName = TRIM(stFileName)//"_last.bin"
      tmpS  = fName
      IF (.NOT.stFileRepl) THEN
         WRITE(fName,'(I3.3)') cTS
         IF (cTS .GE. 1000) fName = STR(cTS)
         fName = TRIM(stFileName)//"_"//TRIM(fName)//".bin"
      END IF
      IF (cm%mas()) THEN
         OPEN(fid, FILE=TRIM(fName))
         CLOSE(fid, STATUS='DELETE')
      END IF

!     This call is to block all processors
      CALL cm%bcast(fid)

      OPEN(fid, FILE=TRIM(fName), ACCESS='DIRECT', RECL=recLn)
      IF (.NOT.ibFlag) THEN
         IF (dFlag) THEN
!           VMS_STRUCT
            IF (sstEq) THEN
!              Prestress
               IF (pstEq) THEN
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn, pS0, Ad
!              Electromechanics
               ELSE IF (cepEq) THEN
                  IF (.NOT.cem%cpld) err = "Incorrect equation "//
     2               "combination. Cannot write restart files"
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn, Ad, Xion, cem%Ya
               ELSE
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn, Ad
               END IF
            ELSE
!              Prestress
               IF (pstEq) THEN
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn, pS0
!              Electromechanics
               ELSE IF (cepEq) THEN
                  IF (.NOT.cem%cpld) err = "Incorrect equation "//
     2               "combination. Cannot write restart files"
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn, Xion, cem%Ya
               ELSE
                  WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1),
     2               eq%iNorm, cplBC%xn, Yn, An, Dn
               END IF
            END IF
         ELSE
!           Electrophysiology
            IF (cepEq) THEN
               WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1),
     2            eq%iNorm, cplBC%xn, Yn, An, Xion
            ELSE
               WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1),
     2            eq%iNorm, cplBC%xn, Yn, An
            END IF
         END IF
      ELSE
         IF (dFlag) THEN
            IF (pstEq) THEN
               WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1),
     2            eq%iNorm, cplBC%xn, Yn, An, Dn, pS0, ib%An, ib%Yn,
     3            ib%Un, ib%Rfb
            ELSE
               WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1),
     2            eq%iNorm, cplBC%xn, Yn, An, Dn, ib%An, ib%Yn, ib%Un,
     3            ib%Rfb
            END IF
         ELSE
            WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1),
     2         eq%iNorm, cplBC%xn, Yn, An, ib%An, ib%Yn, ib%Un, ib%Rfb
         END IF
      END IF
      CLOSE(fid)

      IF (.NOT.stFileRepl .AND. cm%mas()) THEN
         CALL SYSTEM("ln -f "//TRIM(fName)//" "//TRIM(tmpS))
      END IF

      RETURN
      END SUBROUTINE WRITERESTART
!####################################################################
!     Prints norm of the displacement in the solid domain when being
!     solved for prestress
      SUBROUTINE OUTDNORM()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iEq, iDmn, s, e
      REAL(KIND=RKIND) dnorm, vol
      CHARACTER(LEN=stdL) sepLine

      sepLine = REPEAT('-',67)
      std = TRIM(sepLine)
      DO iEq=1, nEq
         s = eq(iEq)%s
         e = s + nsd - 1
         DO iDmn=1, eq(iEq)%nDmn
            IF ( (eq(iEq)%dmn(iDmn)%phys .NE. phys_struct) .AND.
     2           (eq(iEq)%dmn(iDmn)%phys .NE. phys_vms_struct) .AND.
     3           (eq(iEq)%dmn(iDmn)%phys .NE. phys_lElas)  .AND.
     4           (eq(iEq)%dmn(iDmn)%phys .NE. phys_shell) .AND.
     5           (eq(iEq)%dmn(iDmn)%phys .NE. phys_CMM) ) CYCLE
            dnorm = Integ(eq(iEq)%dmn(iDmn)%Id, Dn, s, e)
            vol = eq(iEq)%dmn(iDmn)%v
            IF (.NOT.ISZERO(vol)) dnorm = dnorm / vol
            std = " Displacement norm on domain "//STR(iDmn)//": "//
     2         STR(dnorm)
         END DO
      END DO
      std = TRIM(sepLine)

      RETURN
      END SUBROUTINE OUTDNORM
!####################################################################
!     This is to find if any exception occures, commenting this out
!     untill fortran2003 is standard in all compilers
      SUBROUTINE EXCEPTIONS
c      USE IEEE_EXCEPTIONS
      USE COMMOD, ONLY: stdL
      IMPLICIT NONE
!     I am assuming nExCk is 5, if a compilation error occured, you need
!     to adjust the following arrays accordingely
c      INTEGER(KIND=IKIND), PARAMETER :: nExCk = SIZE(IEEE_ALL)
c      LOGICAL, PARAMETER :: check(nExCk) =
c     2   (/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.FALSE./)
c      CHARACTER(LEN=stdL), PARAMETER :: ieWarn(nExCk) =
c     2   (/"Overflow", "Divide by zero", "Invalid arithmetic operation",
c     3   "Underflow", "Inexact operation"/)
c
c      LOGICAL, SAVE :: iniSet = .TRUE., sprtFlag(nExCk)
c
c      LOGICAL fg
c      INTEGER(KIND=IKIND) i
c      REAL(KIND=RKIND) r
c
c      IF (iniSet) THEN
c         iniSet = .FALSE.
c         DO i=1, nExCk
c            IF (.NOT.check(i)) CYCLE
c            sprtFlag(i) = IEEE_SUPPORT_FLAG(IEEE_ALL(i), r)
c         END DO
c
c         fg = .FALSE.
c         DO i=1, nExCk
c            IF (.NOT.check(i)) CYCLE
c            IF (sprtFlag(i)) CALL IEEE_SET_HALTING_MODE(IEEE_ALL(i), fg)
c         END DO
c      END IF
c
c      DO i=1, nExCk
c         IF (.NOT.check(i)) CYCLE
c         IF (sprtFlag(i)) THEN
c            CALL IEEE_GET_FLAG(IEEE_ALL(i), fg)
c            IF (fg) THEN
c               CALL WARNING(ieWarn(i))
c               CALL IEEE_SET_FLAG(IEEE_ALL(i), .FALSE.)
c            END IF
c         END IF
c      END DO

      RETURN
      END SUBROUTINE EXCEPTIONS
!####################################################################

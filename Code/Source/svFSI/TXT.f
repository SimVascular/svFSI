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
!     This is to calculate average/flux of variables at the boundaries
!     and print them into the "B_" files.
!
!--------------------------------------------------------------------

      SUBROUTINE TXT(flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: flag

      LOGICAL ltmp, wtn(2), div
      INTEGER(KIND=IKIND) fid, i, l, e, s, a, iOut, iEq, oGrp
      CHARACTER(LEN=stdL) fName(2)

      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)

      fid = 1

!     Writing data related to cplBC
      IF (.NOT.resetSim) THEN
         IF (cplBC%coupled) THEN
            IF (.NOT.flag) THEN
               IF (cplBC%useGenBC) THEN
                  CALL genBC_Integ_X('L')
               ELSE
                  ltmp = ANY(BTEST(eq(1)%bc(:)%bType,bType_RCR))
                  CALL cplBC_Integ_X(ltmp)
               END IF
            END IF

            IF (cm%mas() .AND. .NOT.cplBC%useGenBC) THEN
               IF (flag) THEN
                  INQUIRE(FILE=cplBC%saveName, EXIST=ltmp)
                  IF (cTS.EQ.0 .OR. .NOT.ltmp) THEN
                     OPEN(fid, FILE=cplBC%saveName)
                     CLOSE(fid, STATUS='DELETE')
                  ELSE
                     CALL TRIMFILE(cTS,cplBC%saveName)
                  END IF
               ELSE
                  OPEN(fid, FILE=cplBC%saveName, POSITION='APPEND')
                  DO i=1, cplBC%nX
                     WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC%xo(i)
                  END DO
                  DO i=1, cplBC%nXp
                     WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC%xp(i)
                  END DO
                  WRITE(fid,*)
                  CLOSE(fid)
               END IF
            END IF
         END IF
      END IF ! resetSim

      DO iEq=1, nEq
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE (tmpV(maxnsd,tnNo))
         DO iOut=1, eq(iEq)%nOutput
!     Don't write it when it doesn't suppose to be written
            wtn = eq(iEq)%output(iOut)%wtn(2:3)
            IF (ALL(.NOT.wtn)) CYCLE
            l = eq(iEq)%output(iOut)%l
            s = eq(iEq)%s + eq(iEq)%output(iOut)%o
            e = s + l - 1

            oGrp = eq(iEq)%output(iOut)%grp
            div  = .TRUE.
            SELECT CASE (oGrp)
            CASE (outGrp_NA)
               err = "NA outGrp in TXT"
            CASE (outGrp_A)
               tmpV(1:l,:) = An(s:e,:)
            CASE (outGrp_Y)
               tmpV(1:l,:) = Yn(s:e,:)
            CASE (outGrp_D)
               tmpV(1:l,:) = Dn(s:e,:)
            CASE (outGrp_WSS, outGrp_vort, outGrp_trac)
               CALL ALLPOST(tmpV, Yn, Dn, oGrp, iEq)
               DO a=1, tnNo
                  tmpV(1,a) = SQRT(NORM(tmpV(1:nsd,a)))
               END DO
               l = 1
            CASE (outGrp_eFlx, outGrp_hFlx, outGrp_divV, outGrp_J)
               CALL ALLPOST(tmpV, Yn, Dn, oGrp, iEq)
            CASE (outGrp_absV)
               DO a=1, tnNo
                  tmpV(1:l,a) = Yn(1:nsd,a) - Yn(nsd+2:2*nsd+1,a)
               END DO
            CASE (outGrp_I)
               div = .FALSE.
               DO a=1, tnNo
                  tmpV(1:l,a) = 1._RKIND
               END DO
            CASE DEFAULT
               err = "Undefined output"
            END SELECT

            fName = eq(iEq)%sym//"_"//TRIM(eq(iEq)%output(iOut)%name)
            IF (l .EQ. nsd) THEN
               fName(1) = TRIM(appPath)//"B_"//TRIM(fName(1))//
     2            "_flux.txt"
            ELSE
               fName(1) = TRIM(appPath)//"B_"//TRIM(fName(1))//
     2            "_average.txt"
            END IF
            fName(2) = TRIM(appPath)//"V_"//TRIM(fName(2))//
     2         "_average.txt"

            IF (flag) THEN
               CALL CCTXT(eq(iEq), fName, wtn)
            ELSE
               CALL WTXT(eq(iEq), l, fName, tmpV, wtn, div)
            END IF
         END DO

!        IB outputs
         IF (.NOT.ibFlag) CYCLE
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE (tmpV(maxnsd,ib%tnNo))
         DO iOut=1, eq(iEq)%nOutIB
!     Don't write it when it doesn't suppose to be written
            wtn = eq(iEq)%outIB(iOut)%wtn(2:3)
            IF (ALL(.NOT.wtn)) CYCLE
            l = eq(iEq)%outIB(iOut)%l
            s = eq(iEq)%s + eq(iEq)%outIB(iOut)%o
            e = s + l - 1

            oGrp = eq(iEq)%outIB(iOut)%grp
            div  = .TRUE.
            SELECT CASE (oGrp)
            CASE (outGrp_NA)
               err = "NA outGrp in TXT"
            CASE (outGrp_A)
               DO a=1, ib%tnNo
                  tmpV(1:l,a) = ib%An(s:e,a)/REAL(cm%np(), KIND=RKIND)
               END DO
            CASE (outGrp_Y)
               DO a=1, ib%tnNo
                  tmpV(1:l,a) = ib%Yn(s:e,a)/REAL(cm%np(), KIND=RKIND)
               END DO
            CASE (outGrp_D)
               DO a=1, ib%tnNo
                  tmpV(1:l,a) = ib%Un(s:e,a)/REAL(cm%np(), KIND=RKIND)
               END DO
            CASE (outGrp_I)
               div = .FALSE.
               DO a=1, ib%tnNo
                  tmpV(1:l,a) = 1._RKIND
               END DO
            CASE DEFAULT
               err = "Undefined output"
            END SELECT

            fName = eq(iEq)%sym//"_"//TRIM(eq(iEq)%outIB(iOut)%name)
            IF (l .EQ. nsd) THEN
               fName(1) = TRIM(appPath)//"IB_B_"//TRIM(fName(1))//
     2            "_flux.txt"
            ELSE
               fName(1) = TRIM(appPath)//"IB_B_"//TRIM(fName(1))//
     2            "_average.txt"
            END IF
            fName(2) = TRIM(appPath)//"IB_V_"//TRIM(fName(2))//
     2         "_average.txt"

            IF (flag) THEN
               CALL IB_CCTXT(eq(iEq), fName, wtn)
            ELSE
               CALL IB_WTXT(eq(iEq), l, fName, tmpV, wtn, div)
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE TXT
!####################################################################
!     This is to check/create the txt file
      SUBROUTINE CCTXT(lEq, fName, wtn)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      LOGICAL, INTENT(IN) :: wtn(2)

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      LOGICAL flag
      INTEGER(KIND=IKIND) iM, iFa, fid, iDmn, i
      CHARACTER(LEN=stdL) stmp

      fid = 1
      IF (cm%slv()) RETURN

!     i=1 are the boundary values and i=2 are volume values
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         INQUIRE(FILE=TRIM(fName(i)), EXIST=flag)
         IF (cTS.NE.0 .AND. flag) THEN
            CALL TRIMFILE(cTS+3,fName(i))
            CYCLE
         END IF

         OPEN(fid, FILE=TRIM(fName(i)))
         IF (i .EQ. 1) THEN
            DO iM=1, nMsh
               DO iFa=1, msh(iM)%nFa
                  stmp = msh(iM)%fa(iFa)%name
                  IF (LEN(TRIM(stmp)) .LE. prL) THEN
                     WRITE(fid,'(A)', ADVANCE='NO')
     2                  ADJUSTR(stmp(1:prL))//" "
                  ELSE
                     WRITE(fid,'(A)',ADVANCE='NO') TRIM(stmp)//" "
                  END IF
               END DO
            END DO
            WRITE(fid,*)
            DO iM=1, nMsh
               DO iFa=1, msh(iM)%nFa
                  stmp = STR(msh(iM)%fa(iFa)%area,prL)
                  WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmn
               stmp = "DOMAIN-"//STR(lEq%dmn(iDmn)%Id)
               IF (lEq%dmn(iDmn)%Id .EQ. -1) stmp = "ENTIRE"
               WRITE(fid,'(A)', ADVANCE='NO') ADJUSTR(stmp(1:prL))//" "
            END DO
            WRITE(fid,*)
            DO iDmn=1, lEq%nDmn
               stmp = STR(lEq%dmn(iDmn)%v,prL)
               WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
            END DO
         END IF
         WRITE(fid,*)
         WRITE(fid,*)
         CLOSE(fid)
      END DO

      RETURN
      END SUBROUTINE CCTXT
!--------------------------------------------------------------------
!     Only keeping first "n" lines of fileName
      SUBROUTINE TRIMFILE(n, fileName)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      CHARACTER(LEN=*), INTENT(IN) :: fileName

      INTEGER(KIND=IKIND) i, fin, fout, ios
      CHARACTER(LEN=stdL) sTmp, fTmp

      fin = 1
      fout = 2
      fTmp = TRIM(appPath)//".scratchFile.tmp"
      OPEN(fin, FILE=fileName, STATUS="OLD")
      OPEN(fout, FILE=TRIM(fTmp), STATUS="NEW")

!     since lines might be longer stdL, I am using ADVANCE="NO"
      i = 0
      DO
         READ(fin,"(A)",ADVANCE="NO",END=004,IOSTAT=ios) sTmp
         IF (ios .EQ. 0) THEN
            WRITE(fout,"(A)",ADVANCE="NO") sTmp
         ELSE IF (ios .LT. 0) THEN
            i = i + 1
            WRITE(fout,"(A)") TRIM(sTmp)
         ELSE
            err = "Issue with reading "//fileName
         END IF
         IF (i .GE. n) EXIT
      END DO
      REWIND(fin)
      REWIND(fout)

      i = 0
      DO
         READ(fout,"(A)",ADVANCE="NO",IOSTAT=ios) sTmp
         IF (ios .EQ. 0) THEN
            WRITE(fin,"(A)",ADVANCE="NO") sTmp
         ELSE IF (ios .LT. 0) THEN
            i = i + 1
            WRITE(fin,"(A)") TRIM(sTmp)
         ELSE
            err = "Issue with reading scratch file"
         END IF
         IF (i .GE. n) EXIT
      END DO
      CLOSE(fin)
      CLOSE(fout)
      sTmp = "rm "//TRIM(fTmp)
      CALL SYSTEM(TRIM(sTmp))
      RETURN

 004  wrn = "File "//TRIM(fileName)//" is too short"
      CLOSE(fin)
      CLOSE(fout)
      sTmp = "rm "//TRIM(fTmp)
      CALL SYSTEM(TRIM(sTmp))

      RETURN
      END SUBROUTINE TRIMFILE
!####################################################################
!     This is to write to txt file
      SUBROUTINE WTXT(lEq, m, fName, tmpV, wtn, div)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      LOGICAL, INTENT(IN) :: wtn(2), div
      INTEGER(KIND=IKIND), INTENT(IN) :: m
      REAL(KIND=RKIND), INTENT(IN) :: tmpV(maxnsd,tnNo)
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      INTEGER(KIND=IKIND) iM, iFa, fid, i, iDmn
      REAL(KIND=RKIND) tmp

      fid = 1
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         IF (cm%mas()) OPEN(fid, FILE=TRIM(fName(i)), STATUS='OLD',
     2      POSITION='APPEND')

         IF (i .EQ. 1) THEN
            DO iM=1, nMsh
               DO iFa=1, msh(iM)%nFa
                  IF (m .EQ. 1) THEN
                     IF (div) THEN
                        tmp = msh(iM)%fa(iFa)%area
                        tmp = Integ(msh(iM)%fa(iFa),tmpV,1)/tmp
                     ELSE
                        tmp = Integ(msh(iM)%fa(iFa),tmpV,1)
                     END IF
                  ELSE IF (m .EQ. nsd) THEN
                     tmp = Integ(msh(iM)%fa(iFa),tmpV,1,m)
                  ELSE
                     err = "WTXT only accepts 1 and nsd"
                  END IF
                  IF (cm%mas())
     2               WRITE(fid,'(A)',ADVANCE='NO') STR(tmp,prL)//" "
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmn
               IF (div) THEN
                  tmp = lEq%dmn(iDmn)%v
                  tmp = Integ(lEq%dmn(iDmn)%Id, tmpV, 1, m)/tmp
               ELSE
                  tmp = Integ(lEq%dmn(iDmn)%Id, tmpV, 1, m)
               END IF
               IF (cm%mas())
     2            WRITE(fid,'(A)', ADVANCE='NO') STR(tmp,prL)//" "
            END DO
         END IF
         IF (cm%mas()) THEN
            WRITE(fid,*)
            CLOSE(fid)
         END IF
      END DO

      RETURN
      END SUBROUTINE WTXT
!####################################################################


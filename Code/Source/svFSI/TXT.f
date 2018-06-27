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

      INTEGER fid, i, l, e, s, a, iOut, iEq, oGrp
      CHARACTER(LEN=stdL) fName(2)

      LOGICAL ltmp, wtn(2)
      REAL(KIND=8), ALLOCATABLE :: tmpV(:,:)
      
      fid = 1
      ALLOCATE (tmpV(maxnsd,tnNo))

!     Writing data related to cplBC 
      IF (.NOT.resetSim) THEN     
         IF (cm%mas() .AND. cplBC%coupled) THEN
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
                  WRITE(fid,'(A)',ADVANCE='NO') STR(cplBC%xo(i))//" "
               END DO
               WRITE(fid,*) 
               CLOSE(fid)
            END IF
         END IF
      END IF ! resetSim

      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutput
!     Don't write it when it doesn't suppose to be written
            wtn = eq(iEq)%output(iOut)%wtn(2:3)
            IF (ALL(.NOT.wtn)) CYCLE
            l = eq(iEq)%output(iOut)%l
            s = eq(iEq)%s + eq(iEq)%output(iOut)%o
            e = s + l - 1
            
            oGrp = eq(iEq)%output(iOut)%grp
            SELECT CASE (oGrp)
            CASE (outGrp_NA)
               err = "NA outGrp in TXT"
            CASE (outGrp_A)
               tmpV(1:l,:) = An(s:e,:)
            CASE (outGrp_Y)
               tmpV(1:l,:) = Yn(s:e,:)
            CASE (outGrp_D)
               tmpV(1:l,:) = Dn(s:e,:)
            CASE (outGrp_WSS, outGrp_vort)
               CALL ALLPOST(tmpV, Yn, Dn, oGrp, iEq)
               DO a=1, tnNo
                  tmpV(1,a) = SQRT(NORM(tmpV(1:nsd,a)))
               END DO
               l = 1
            CASE (outGrp_eFlx, outGrp_hFlx)
               CALL ALLPOST(tmpV, Yn, Dn, oGrp, iEq)
            CASE (outGrp_absV)
               DO a=1, tnNo
                  tmpV(1:l,a) = Yn(1:nsd,a) - Yn(nsd+2:2*nsd+1,a)
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
               CALL WTXT(eq(iEq), l, fName, tmpV, wtn)
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE TXT

!--------------------------------------------------------------------
!     This is to check/create the txt file
      SUBROUTINE CCTXT(lEq, fName, wtn)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(eqType), INTENT(IN) :: lEq
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      LOGICAL, INTENT(IN) :: wtn(2)

      INTEGER, PARAMETER :: prL = 10

      LOGICAL flag
      INTEGER iM, iFa, fid, iDmn, i
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
!     This is to write to txt file
      SUBROUTINE WTXT(lEq, m, fName, tmpV, wtn)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(eqType), INTENT(IN) :: lEq
      INTEGER, INTENT(IN) :: m
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      REAL(KIND=8), INTENT(IN) :: tmpV(maxnsd,tnNo)
      LOGICAL, INTENT(IN) :: wtn(2)

      INTEGER, PARAMETER :: prL = 10
      
      INTEGER iM, iFa, fid, i, iDmn
      REAL(KIND=8) tmp

      fid = 1
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         IF (cm%mas()) OPEN(fid, FILE=TRIM(fName(i)), STATUS='OLD', 
     2      POSITION='APPEND')

         IF (i .EQ. 1) THEN
            DO iM=1, nMsh
               DO iFa=1, msh(iM)%nFa
                  IF (m .EQ. 1) THEN
                     tmp = msh(iM)%fa(iFa)%area
                     tmp = Integ(msh(iM)%fa(iFa),tmpV,1)/tmp
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
               tmp = Integ(lEq%dmn(iDmn)%Id, tmpV, 1, m)/lEq%dmn(iDmn)%v
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

!--------------------------------------------------------------------
!     Only keeping first "n" lines of fileName      
      SUBROUTINE TRIMFILE(n, fileName)

      USE COMMOD 
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      CHARACTER(LEN=*), INTENT(IN) :: fileName

      INTEGER i, fin, fout, ios
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


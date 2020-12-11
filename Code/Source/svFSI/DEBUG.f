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
!     Some useful debugging subroutines
!
!--------------------------------------------------------------------

!     Writes solution data and residue to file for debugging. Can also
!     be used with multiple processes
      SUBROUTINE DEBUGR()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a, i, Ac, fid, iM, nNo, gnNo, s, e, ldof
      CHARACTER(LEN=stdL) fName

      REAL(KIND=RKIND), ALLOCATABLE :: lX(:,:), lA(:,:), lY(:,:),
     2   lD(:,:), lR(:,:), gX(:,:), gA(:,:), gY(:,:), gD(:,:), gR(:,:)

      fid = 1889
      WRITE(fName,'(A)') TRIM(appPath)//"dbgR_"//TRIM(eq(cEq)%sym)//
     2   "_"//STR(cTS)//"_"//STR(eq(cEq)%itr)
      IF (cm%mas()) OPEN(fid,FILE=TRIM(fName))

      s = eq(cEq)%s
      e = eq(cEq)%e
      ldof = e-s+1
      DO iM=1, nMsh
         nNo  = msh(iM)%nNo
         gnNo = msh(iM)%gnNo
         ALLOCATE(lX(nsd,nNo), lA(lDof,nNo), lY(lDof,nNo), lD(lDof,nNo),
     2      lR(dof,nNo))
         IF (cm%mas()) THEN
            ALLOCATE(gX(nsd,gnNo), gA(lDof,gnNo), gY(lDof,gnNo),
     2       gD(lDof,tnNo), gR(dof,gnNo))
         ELSE
            ALLOCATE(gX(0,0), gA(0,0), gY(0,0), gD(0,0), gR(0,0))
         END IF
         lX = 0._RKIND
         lA = 0._RKIND
         lY = 0._RKIND
         lD = 0._RKIND
         lR = 0._RKIND
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            lX(:,a) = x(:,Ac)
            lA(:,a) = An(s:e,Ac)
            lY(:,a) = Yn(s:e,Ac)
            lD(:,a) = Dn(s:e,Ac)
            lR(:,a) = R(:,Ac)
         END DO

         gX = GLOBAL(msh(iM), lX)
         gA = GLOBAL(msh(iM), lA)
         gY = GLOBAL(msh(iM), lY)
         gD = GLOBAL(msh(iM), lD)
         gR = GLOBAL(msh(iM), lR)

         IF (cm%mas()) THEN
            WRITE(fid,'(A)') "Mesh: <"//TRIM(msh(iM)%name)//">"
            DO a=1, gnNo
               WRITE(fid,'(A)',ADVANCE='NO') STR(a)
               DO i=1, nsd
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gX(i,a))
               END DO
               DO i=1, dof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gR(i,a))
               END DO
               DO i=1, lDof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gA(i,a))
               END DO
               DO i=1, lDof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gY(i,a))
               END DO
               DO i=1, lDof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gD(i,a))
               END DO
               WRITE(fid,'(A)')
            END DO
         END IF
         DEALLOCATE(lX, lA, lY, lD, lR, gX, gA, gY, gD, gR)
      END DO

      IF (cm%mas()) CLOSE(fid)

      iM = cm%reduce(iM)

      RETURN
      END SUBROUTINE DEBUGR
!--------------------------------------------------------------------
      SUBROUTINE PDEBUGVALR()
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a, b, i, j, k, fid
      CHARACTER(LEN=stdL) fName

      fid = 1256 + cm%tF()
      i   = eq(cEq)%itr
      WRITE(fName,'(A)') TRIM(appPath)//"Val_R_"//STR(cTS)//"_"//
     2   STR(i)//"_"//STR(cm%tF())
      OPEN(fid,FILE=TRIM(fName))
      DO a=1, tnNo
         WRITE(fid,'(A)') REPEAT('=',32)
         WRITE(fid,'(A)',ADVANCE='NO') "Row: "//STR(a)//" grow: "//
     2      STR(ltg(a))//" x: "

         DO i=1, nsd
            WRITE(fid,'(A)',ADVANCE='NO') " "//STR(x(i,a))
         END DO
         WRITE(fid,'(A)')

         WRITE(fid,'(4X,A)',ADVANCE='NO') "R:"
         DO i=1, dof
            IF (i .EQ. 1) THEN
               WRITE(fid,'(1X,1pE25.18)') R(i,a)
            ELSE
               WRITE(fid,'(7X,1pE25.18)') R(i,a)
            END IF
         END DO

         DO i=rowPtr(a), rowPtr(a+1)-1
            b = colPtr(i)
            WRITE(fid,'(4X,A)') REPEAT('-',36)
            WRITE(fid,'(4X,A)') "Col: "//STR(b)//" gcol: "//STR(ltg(b))
            WRITE(fid,'(4X,A)',ADVANCE='NO') "R: "
            DO j=1, dof
               WRITE(fid,'(A)',ADVANCE='NO') " "//STR(R(j,b))
            END DO
            WRITE(fid,'(A)')

            DO k=1, dof
               WRITE(fid,'(4X,A)',ADVANCE='NO')
               DO j=1, dof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//
     2               STR(Val((k-1)*dof+j,i))
               END DO
               WRITE(fid,'(A)')
            END DO
         END DO
      END DO
      CLOSE(fid)

      RETURN
      END SUBROUTINE PDEBUGVALR
!--------------------------------------------------------------------
      SUBROUTINE DEBUGVALR(K, U)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: K(dof*dof,lhs%nnz), U(dof,tnNo)

      INTEGER(KIND=IKIND) a, b, i, Ac, iM, nNo, gnNo, fid
      REAL(KIND=RKIND) :: KU(dof,tnNo)
      CHARACTER(LEN=stdL) fName

      REAL(KIND=RKIND), ALLOCATABLE :: lX(:,:), gX(:,:), lR(:,:),
     2   gR(:,:), lU(:,:), gU(:,:)

      KU(:,:) = 0._RKIND
      SELECT CASE (dof)
      CASE (1)
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               b = colPtr(i)
               KU(1,a) = KU(1,a) + K(1 ,i)*U(1,b) + K(2 ,i)*U(2,b) +
     2                             K(3 ,i)*U(3,b) + K(4 ,i)*U(4,b)
            END DO
         END DO
      CASE (2)
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               b = colPtr(i)
               KU(1,a) = KU(1,a) + K(1 ,i)*U(1,b) + K(2 ,i)*U(2,b) +
     2                             K(3 ,i)*U(3,b) + K(4 ,i)*U(4,b)
               KU(2,a) = KU(2,a) + K(5 ,i)*U(1,b) + K(6 ,i)*U(2,b) +
     2                             K(7 ,i)*U(3,b) + K(8 ,i)*U(4,b)
            END DO
         END DO
      CASE (3)
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               b = colPtr(i)
               KU(1,a) = KU(1,a) + K(1 ,i)*U(1,b) + K(2 ,i)*U(2,b) +
     2                             K(3 ,i)*U(3,b) + K(4 ,i)*U(4,b)
               KU(2,a) = KU(2,a) + K(5 ,i)*U(1,b) + K(6 ,i)*U(2,b) +
     2                             K(7 ,i)*U(3,b) + K(8 ,i)*U(4,b)
               KU(3,a) = KU(3,a) + K(9 ,i)*U(1,b) + K(10,i)*U(2,b) +
     2                             K(11,i)*U(3,b) + K(12,i)*U(4,b)
            END DO
         END DO
      CASE (4)
         DO a=1, tnNo
            DO i=rowPtr(a), rowPtr(a+1)-1
               b = colPtr(i)
               KU(1,a) = KU(1,a) + K(1 ,i)*U(1,b) + K(2 ,i)*U(2,b) +
     2                             K(3 ,i)*U(3,b) + K(4 ,i)*U(4,b)
               KU(2,a) = KU(2,a) + K(5 ,i)*U(1,b) + K(6 ,i)*U(2,b) +
     2                             K(7 ,i)*U(3,b) + K(8 ,i)*U(4,b)
               KU(3,a) = KU(3,a) + K(9 ,i)*U(1,b) + K(10,i)*U(2,b) +
     2                             K(11,i)*U(3,b) + K(12,i)*U(4,b)
               KU(4,a) = KU(4,a) + K(13,i)*U(1,b) + K(14,i)*U(2,b) +
     2                             K(15,i)*U(3,b) + K(16,i)*U(4,b)
            END DO
         END DO
      END SELECT

      CALL COMMU(KU)

      fid = 1256
      i   = eq(cEq)%itr
      WRITE(fName,'(A)') TRIM(appPath)//"dbgValR_"//TRIM(eq(cEq)%sym)//
     2   "_"//STR(cTS)//"_"//STR(i)
      IF (cm%mas()) OPEN(fid,FILE=TRIM(fName))

      DO iM=1, nMsh
         nNo  = msh(iM)%nNo
         gnNo = msh(iM)%gnNo
         ALLOCATE(lX(nsd,nNo), lR(dof,nNo), lU(dof,nNo))
         IF (cm%mas()) THEN
            ALLOCATE(gX(nsd,gnNo), gR(dof,gnNo), gU(dof,gnNo))
         ELSE
            ALLOCATE(gX(0,0), gR(0,0), gU(0,0))
         END IF
         DO a=1, nNo
            Ac = msh(iM)%gN(a)
            lX(:,a) = x (:,Ac)
            lR(:,a) = R (:,Ac)
            lU(:,a) = KU(:,Ac)
         END DO

         gX = GLOBAL(msh(iM), lX)
         gR = GLOBAL(msh(iM), lR)
         gU = GLOBAL(msh(iM), lU)

         IF (cm%mas()) THEN
            WRITE(fid,'(A)') REPEAT('=',36)
            WRITE(fid,'(A)') "Mesh: <"//TRIM(msh(iM)%name)//">"
            DO a=1, gnNo
               WRITE(fid,'(A)',ADVANCE='NO') STR(a)
               DO i=1, nsd
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gX(i,a))
               END DO
               DO i=1, dof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gR(i,a))
               END DO
               DO i=1, dof
                  WRITE(fid,'(A)',ADVANCE='NO') " "//STR(gU(i,a))
               END DO
               WRITE(fid,'(A)')
            END DO
         END IF
         DEALLOCATE(lX, gX, lR, lU, gR, gU)
      END DO

      IF (cm%mas()) CLOSE(fid)

      iM = cm%reduce(iM)

      RETURN
      END SUBROUTINE DEBUGVALR
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix STIFF, only works for serial code
      SUBROUTINE OUTPUTVALR(lVal, lR, Prefix)
      USE TYPEMOD
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: lVal(dof*dof,lhs%nnz), 
     2   lR(dof,lhs%nNo)
      CHARACTER(LEN=*), INTENT(IN) :: Prefix

      INTEGER(KIND=IKIND) ptr, rowN, colN, left, right, fid
      INTEGER(KIND=IKIND) iref, jref, recLn_K, recLn_R
      REAL(KIND=RKIND) :: STIFF(dof*tnNo,dof*tnNo)
      REAL(KIND=RKIND) :: RR(dof*tnNo)
      INQUIRE(IOLENGTH=recLn_K) STIFF
      INQUIRE(IOLENGTH=recLn_R) RR

      PRINT *, "dof=", dof 
      PRINT *, "tnNo=", tnNo

!     Write stiffness matrix in binary
      STIFF = 0._RKIND
      DO rowN=1, tnNo
         iref = dof*(rowN-1)
         left  = rowPtr(rowN)
         right = rowPtr(rowN+1) - 1
         DO ptr=left, right
            colN = colPtr(ptr)
            jref = dof*(colN-1)
            STIFF((iref+1):(iref+dof),(jref+1):(jref+dof)) = 
     2          TRANSPOSE(RESHAPE(lVal(:,ptr),(/dof,dof/)))
         END DO
      END DO
      fid = 59999
      OPEN(UNIT=fid,FILE=TRIM(Prefix)//'_K.bin',STATUS='replace',
     2     FORM='unformatted',ACCESS='direct',recl=recLn_K)    
      WRITE(fid,rec=1) STIFF
      CLOSE(fid)

!     Write right hand side in binary 
      RR = RESHAPE(lR,(/dof*tnNo/))
      fid = 59999
      OPEN(UNIT=fid,FILE=TRIM(Prefix)//'_R.bin',STATUS='replace',
     2     FORM='unformatted',ACCESS='direct',recl=recLn_R)    
      WRITE(fid,rec=1) RR
      CLOSE(fid)

      RETURN
      END SUBROUTINE OUTPUTVALR
!####################################################################

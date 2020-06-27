!--------------------------------------------------------------------
!     Created by Mahdi Esmaily Moghadam
!     contact memt63@gmail.com for reporting the bugs.
!--------------------------------------------------------------------
!
!     UC Copyright Notice
!     This software is Copyright ©2012 The Regents of the University of
!     California. All Rights Reserved.
!
!     Permission to copy and modify this software and its documentation
!     for educational, research and non-profit purposes, without fee,
!     and without a written agreement is hereby granted, provided that
!     the above copyright notice, this paragraph and the following three
!     paragraphs appear in all copies.
!
!     Permission to make commercial use of this software may be obtained
!     by contacting:
!     Technology Transfer Office
!     9500 Gilman Drive, Mail Code 0910
!     University of California
!     La Jolla, CA 92093-0910
!     (858) 534-5815
!     invent@ucsd.edu
!
!     This software program and documentation are copyrighted by The
!     Regents of the University of California. The software program and
!     documentation are supplied "as is", without any accompanying
!     services from The Regents. The Regents does not warrant that the
!     operation of the program will be uninterrupted or error-free. The
!     end-user understands that the program was developed for research
!     purposes and is advised not to rely exclusively on the program for
!     any reason.
!
!     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
!     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
!     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
!     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
!     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
!     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
!     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
!     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
!     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
!
!--------------------------------------------------------------------
!     Jacobi symmetic preconditioner, to precondition both LHS and RHS.
!--------------------------------------------------------------------

      SUBROUTINE PRECONDDIAG(lhs, rowPtr, colPtr, diagPtr, dof, Val, R, &
     &   W)
      INCLUDE "FSILS_STD.h"
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: rowPtr(2,lhs%nNo),              &
     &  colPtr(lhs%nnz), diagPtr(lhs%nNo)
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(INOUT) :: Val(dof*dof,lhs%nnz),           &
     &   R(dof,lhs%nNo)
      REAL(KIND=LSRP), INTENT(OUT) :: W(dof,lhs%nNo)

      INTEGER(KIND=LSIP) nNo, i, a, d, Ac, faIn

      nNo = lhs%nNo

!     Calculating W: W = diag(K)
      SELECT CASE (dof)
      CASE (1)
         DO Ac=1, nNo
            W(1,Ac) = Val(1,diagPtr(Ac))
         END DO
      CASE(2)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            W(1,Ac) = Val(1,d)
            W(2,Ac) = Val(4,d)
         END DO
      CASE(3)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            W(1,Ac) = Val(1,d)
            W(2,Ac) = Val(5,d)
            W(3,Ac) = Val(9,d)
         END DO
      CASE(4)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            W(1,Ac) = Val(1,d)
            W(2,Ac) = Val(6,d)
            W(3,Ac) = Val(11,d)
            W(4,Ac) = Val(16,d)
         END DO
      CASE DEFAULT
         DO Ac=1, nNo
            d = diagPtr(Ac)
            DO i=1, dof
               W(i,Ac) = Val(i*dof-dof+i,d)
            END DO
         END DO
      END SELECT

      CALL FSILS_COMMUV(lhs, dof, W)

!     Accounting for Dirichlet BC and inversing W = W^{-1/2}
      DO Ac=1, nNo
         d = diagPtr(Ac)
         DO i=1, dof
            IF (W(i,Ac) .EQ. 0._LSRP) W(i,Ac) = 1._LSRP
         END DO
      END DO

      W = 1._LSRP/SQRT(ABS(W))
      DO faIn=1, lhs%nFaces
         IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
         i = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               W(1:i,Ac) = W(1:i,Ac)*lhs%face(faIn)%val(1:i,a)
            END DO
         END IF
      END DO

!     Pre-multipling K with W: K = W*K
      CALL PREMUL(lhs, rowPtr, dof, Val, W)

!     Multipling R with W: R = W*R
      R = W*R

!     Now post-multipling K by W: K = K*W
      CALL POSMUL(lhs, rowPtr, colPtr, dof, Val, W)

      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, MIN(lhs%face(faIn)%dof,dof)
                  lhs%face(faIn)%valM(i,a) =                            &
     &               lhs%face(faIn)%val(i,a)*W(i,Ac)
               END DO
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE PRECONDDIAG

!--------------------------------------------------------------------
!     Row and column preconditioner, to precondition both LHS and RHS.
!--------------------------------------------------------------------
      SUBROUTINE PRECONDRNC(lhs, rowPtr, colPtr, diagPtr, dof, Val, R, 
     2   W1, W2)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: rowPtr(2,lhs%nNo),       
     2   colPtr(lhs%nnz), diagPtr(lhs%nNo)
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(INOUT) :: Val(dof*dof,lhs%nnz), 
     2   R(dof,lhs%nNo)
      REAL(KIND=LSRP), INTENT(OUT) :: W1(dof,lhs%nNo), W2(dof,lhs%nNo)

      REAL(KIND=LSRP) :: Wr(dof,lhs%nNo), Wc(dof,lhs%nNo), tol
      INTEGER(KIND=LSIP) nNo, i, j, a, b, d, Ac, faIn
      INTEGER(KIND=LSIP) iter, maxiter, ierr
      LOGICAL flag, gflag(lhs%commu%nTasks)

      nNo     = lhs%nNo
      maxiter = 10
      tol     = 2._LSRP
      iter    = 0
      flag    = .TRUE.
      W1      = 1._LSRP
      W2      = 1._LSRP
      
      !*****************************************************
      ! Apply Dirichlet BC
      !*****************************************************
      Wr = 1._LSRP
      DO faIn=1, lhs%nFaces
         IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
         i = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               Wr(1:i,Ac) = Wr(1:i,Ac)*lhs%face(faIn)%val(1:i,a)
            END DO
         END IF
      END DO
      CALL FSILS_COMMUV(lhs, dof, Wr)
      ! For parallel case, val and Wr can be larger than 1 due to
      ! the addition operator in FSILS_COMMUV. Hence need renormalization.
      Wr = Wr - 0.5_LSRP
      Wr = Wr/ABS(Wr)
      Wr = (Wr + ABS(Wr))*0.5_LSRP
      ! Kill the row and column corresponding to Dirichlet BC
      CALL PREMUL(lhs, rowPtr, dof, Val, Wr)
      R = Wr*R
      CALL POSMUL(lhs, rowPtr, colPtr, dof, Val, Wr)
      ! Set diagnal term to one
      SELECT CASE (dof)
      CASE (1)
         DO Ac=1, nNo
            d        = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._LSRP) + 1._LSRP 
         END DO
      CASE(2)
         DO Ac=1, nNo
            d        = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._LSRP) + 1._LSRP 
            Val(4,d) = Wr(2,Ac)*(Val(4,d)-1._LSRP) + 1._LSRP 
         END DO
      CASE(3)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._LSRP) + 1._LSRP 
            Val(5,d) = Wr(2,Ac)*(Val(5,d)-1._LSRP) + 1._LSRP 
            Val(9,d) = Wr(3,Ac)*(Val(9,d)-1._LSRP) + 1._LSRP 
         END DO
      CASE(4)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            Val(1 ,d) = Wr(1,Ac)*(Val(1 ,d)-1._LSRP) + 1._LSRP 
            Val(6 ,d) = Wr(2,Ac)*(Val(6 ,d)-1._LSRP) + 1._LSRP 
            Val(11,d) = Wr(3,Ac)*(Val(11,d)-1._LSRP) + 1._LSRP 
            Val(16,d) = Wr(4,Ac)*(Val(16,d)-1._LSRP) + 1._LSRP 
         END DO
      CASE DEFAULT
         DO Ac=1, nNo
            d = diagPtr(Ac)
            DO i=1, dof
               Val(i*dof-dof+i,d) = Wr(i,Ac)*(Val(i*dof-dof+i,d)
     2                            - 1._LSRP) + 1._LSRP
            END DO
         END DO
      END SELECT

      !*****************************************************
      ! Row and column scaling
      !*****************************************************
      DO While (flag)
         Wr = 0._LSRP
         Wc = 0._LSRP
         iter = iter + 1
         IF (iter .GE. maxiter) THEN 
            PRINT *, "Warning: maximum iteration number reached"//
     2            "@ SUBROUTINE PRECONDRNC."
            PRINT *, MAXVAL(Wr), MAXVAL(Wc)
            flag = .False.
         END IF

         ! Max norm along row and column
         SELECT CASE (dof)
         CASE (1)
            DO Ac=1, nNo
               a        = rowPtr(1,Ac)
               b        = rowPtr(2,Ac)
               Wr(1,Ac) = MAXVAL(ABS(Val(1,a:b)))

               DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                  a = colPtr(i)
                  Wc(1,a) = MAX(ABS(Val(1,i)),Wc(1,a))
               END DO
            END DO
         CASE(2)
            DO Ac=1, nNo
               a        = rowPtr(1,Ac)
               b        = rowPtr(2,Ac)
               Wr(1,Ac) = MAXVAL(ABS(Val(1:2,a:b)))
               Wr(2,Ac) = MAXVAL(ABS(Val(3:4,a:b)))

               DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                  a = colPtr(i)
                  Wc(1,a) = MAX(MAXVAL(ABS(Val(1:3:2,i))), Wc(1,a))
                  Wc(2,a) = MAX(MAXVAL(ABS(Val(2:4:2,i))), Wc(2,a))
               END DO
            END DO
         CASE(3)
            DO Ac=1, nNo
               a        = rowPtr(1,Ac)
               b        = rowPtr(2,Ac)
               Wr(1,Ac) = MAXVAL(ABS(Val(1:3,a:b)))
               Wr(2,Ac) = MAXVAL(ABS(Val(4:6,a:b)))
               Wr(3,Ac) = MAXVAL(ABS(Val(7:9,a:b)))

               DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                  a = colPtr(i)
                  Wc(1,a) = MAX(MAXVAL(ABS(Val(1:7:3,i))), Wc(1,a))
                  Wc(2,a) = MAX(MAXVAL(ABS(Val(2:8:3,i))), Wc(2,a))
                  Wc(3,a) = MAX(MAXVAL(ABS(Val(3:9:3,i))), Wc(3,a))
               END DO
            END DO
         CASE(4)
            DO Ac=1, nNo
               a        = rowPtr(1,Ac)
               b        = rowPtr(2,Ac)
               Wr(1,Ac) = MAXVAL(ABS(Val(1 :4 ,a:b)))
               Wr(2,Ac) = MAXVAL(ABS(Val(5 :8 ,a:b)))
               Wr(3,Ac) = MAXVAL(ABS(Val(9 :12,a:b)))
               Wr(4,Ac) = MAXVAL(ABS(Val(13:16,a:b)))

               DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                  a = colPtr(i)
                  Wc(1,a) = MAX(MAXVAL(ABS(Val(1:13:4,i))), Wc(1,a))
                  Wc(2,a) = MAX(MAXVAL(ABS(Val(2:14:4,i))), Wc(2,a))
                  Wc(3,a) = MAX(MAXVAL(ABS(Val(3:15:4,i))), Wc(3,a))
                  Wc(4,a) = MAX(MAXVAL(ABS(Val(4:16:4,i))), Wc(4,a))
               END DO
            END DO
         CASE DEFAULT
            DO Ac=1, nNo
               a = rowPtr(1,Ac)
               b = rowPtr(2,Ac)
               DO i=1, dof
                  j = i*dof - dof + 1
                  Wr(i,Ac) = MAXVAL(ABS(Val(j:i*dof,a:b)))
               END DO
               DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                  a = colPtr(i)
                  DO b=1, dof
                     j = dof*(dof-1) + b
                     Wc(b,a) = MAX(MAXVAL(ABS(Val(b:j:dof,i))),Wc(b,a))
                  END DO
               END DO
            END DO
         END SELECT ! dof

         CALL FSILS_COMMUV(lhs, dof, Wr)
         CALL FSILS_COMMUV(lhs, dof, Wc)
         
         Wr = 1._LSRP/SQRT(Wr)
         Wc = 1._LSRP/SQRT(Wc)

         CALL PREMUL(lhs, rowPtr, dof, Val, Wr)
         CALL POSMUL(lhs, rowPtr, colPtr, dof, Val, Wc)

         W1 = W1*Wr
         W2 = W2*Wc

         IF (ABS( 1._LSRP - MAXVAL(Wr) ) .LT. tol .AND. 
     2       ABS( 1._LSRP - MAXVAL(Wc) ) .LT. tol) flag = .False.

         IF (lhs%commu%nTasks .GT. 1) THEN 
            CALL MPI_ALLGATHER(flag, 1, mplog, gflag, 1, mplog, 
     2            lhs%commu%comm, ierr)
            flag = ANY(gflag)
         END IF

      END DO ! do while

      ! Multipling R with Wr: R = Wr*R
      R = W1*R


   !    DO faIn=1, lhs%nFaces
   !       IF (lhs%face(faIn)%coupledFlag) THEN
   !          DO a=1, lhs%face(faIn)%nNo
   !             Ac = lhs%face(faIn)%glob(a)
   !             DO i=1, MIN(lhs%face(faIn)%dof,dof)
   !                lhs%face(faIn)%valM(i,a) =                        
   !   2               lhs%face(faIn)%val(i,a)*W(i,Ac)
   !             END DO
   !          END DO
   !       END IF
   !    END DO

      RETURN
      END SUBROUTINE PRECONDRNC

!--------------------------------------------------------------------
!     Pre-multipling Val with W: Val = W*Val
      SUBROUTINE PREMUL(lhs, rowPtr, dof, Val, W)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: rowPtr(2,lhs%nNo)
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(INOUT) :: Val(dof*dof,lhs%nnz)
      REAL(KIND=LSRP), INTENT(IN) :: W(dof,lhs%nNo)

      INTEGER(KIND=LSIP) nNo, i, j, a, b, Ac
      
      nNo = lhs%nNo

      SELECT CASE (dof)
      CASE (1)
         DO Ac=1, nNo
            a          = rowPtr(1,Ac)
            b          = rowPtr(2,Ac)
            Val(1,a:b) = Val(1,a:b)*W(1,Ac)
         END DO
      CASE(2)
         DO Ac=1, nNo
            a            = rowPtr(1,Ac)
            b            = rowPtr(2,Ac)
            Val(1:2,a:b) = Val(1:2,a:b)*W(1,Ac)
            Val(3:4,a:b) = Val(3:4,a:b)*W(2,Ac)
         END DO
      CASE(3)
         DO Ac=1, nNo
            a            = rowPtr(1,Ac)
            b            = rowPtr(2,Ac)
            Val(1:3,a:b) = Val(1:3,a:b)*W(1,Ac)
            Val(4:6,a:b) = Val(4:6,a:b)*W(2,Ac)
            Val(7:9,a:b) = Val(7:9,a:b)*W(3,Ac)
         END DO
      CASE(4)
         DO Ac=1, nNo
            a              = rowPtr(1,Ac)
            b              = rowPtr(2,Ac)
            Val(1:4,a:b)   = Val(1:4,a:b)*W(1,Ac)
            Val(5:8,a:b)   = Val(5:8,a:b)*W(2,Ac)
            Val(9:12,a:b)  = Val(9:12,a:b)*W(3,Ac)
            Val(13:16,a:b) = Val(13:16,a:b)*W(4,Ac)
         END DO
      CASE DEFAULT
         DO Ac=1, nNo
            a = rowPtr(1,Ac)
            b = rowPtr(2,Ac)
            DO i=1, dof
               j = i*dof - dof + 1
               Val(j:i*dof,a:b) = Val(j:i*dof,a:b)*W(i,Ac)
            END DO
         END DO
      END SELECT

      RETURN
      END SUBROUTINE PREMUL

!--------------------------------------------------------------------
!     Post-multipling Val by W: Val = Val*W
      SUBROUTINE POSMUL(lhs, rowPtr, colPtr, dof, Val, W)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: rowPtr(2,lhs%nNo),
     2                                  colPtr(lhs%nnz)
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(INOUT) :: Val(dof*dof,lhs%nnz)
      REAL(KIND=LSRP), INTENT(IN) :: W(dof,lhs%nNo)

      INTEGER(KIND=LSIP) nNo, i, j, a, b, Ac
      
      nNo = lhs%nNo

      SELECT CASE (dof)
      CASE (1)
         DO Ac=1, nNo
            DO i=rowPtr(1,Ac), rowPtr(2,Ac)
               a = colPtr(i)
               Val(1,i) = Val(1,i)*W(1,a)
            END DO
         END DO
      CASE (2)
         DO Ac=1, nNo
            DO i=rowPtr(1,Ac), rowPtr(2,Ac)
               a = colPtr(i)
               Val(1:3:2,i) = Val(1:3:2,i)*W(1,a)
               Val(2:4:2,i) = Val(2:4:2,i)*W(2,a)
            END DO
         END DO
      CASE (3)
         DO Ac=1, nNo
            DO i=rowPtr(1,Ac), rowPtr(2,Ac)
               a = colPtr(i)
               Val(1:7:3,i) = Val(1:7:3,i)*W(1,a)
               Val(2:8:3,i) = Val(2:8:3,i)*W(2,a)
               Val(3:9:3,i) = Val(3:9:3,i)*W(3,a)
            END DO
         END DO
      CASE (4)
         DO Ac=1, nNo
            DO i=rowPtr(1,Ac), rowPtr(2,Ac)
               a = colPtr(i)
               Val(1:13:4,i) = Val(1:13:4,i)*W(1,a)
               Val(2:14:4,i) = Val(2:14:4,i)*W(2,a)
               Val(3:15:4,i) = Val(3:15:4,i)*W(3,a)
               Val(4:16:4,i) = Val(4:16:4,i)*W(4,a)
            END DO
         END DO
      CASE DEFAULT
         DO Ac=1, nNo
            DO i=rowPtr(1,Ac), rowPtr(2,Ac)
               a = colPtr(i)
               DO b=1, dof
                  j = dof*(dof-1) + b
                  Val(b:j:dof,i) = Val(b:j:dof,i)*W(b,a)
               END DO
            END DO
         END DO
      END SELECT

      RETURN 
      END SUBROUTINE POSMUL

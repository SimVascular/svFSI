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

      SUBROUTINE PRECOND(lhs, rowPtr, colPtr, diagPtr, dof, Val, R, W)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz),        &
     &   diagPtr(lhs%nNo)
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(INOUT) :: Val(dof*dof,lhs%nnz),              &
     &   R(dof,lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: W(dof,lhs%nNo)

      INTEGER nNo, i, j, a, b, d, Ac, faIn
 
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
!     Accounding for Dirichlet BC and inversing W = W^{-1/2}
      
      DO Ac=1, nNo
         d = diagPtr(Ac)
         DO i=1, dof
            IF (W(i,Ac) .EQ. 0D0) W(i,Ac) = 1D0
         END DO
      END DO

      W = 1D0/SQRT(ABS(W))
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

!     Multipling R with W: R = W*R
      R = W*R
      
!     Now post-multipling K by W: K = K*W
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
      END SUBROUTINE PRECOND

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
!     Graduate minimum residual algorithm is implemented here for vector
!     and scaler problems. 
!--------------------------------------------------------------------
      
      SUBROUTINE GMRES(lhs, ls, dof, Val, R, X)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: Val(dof*dof,lhs%nnz), R(dof,lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: X(dof,lhs%nNo)
     
      INTEGER nNo, mynNo, i, j, k, l
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMV, FSILS_DOTV, FSILS_NCDOTV
      REAL(KIND=8) eps, tmp, time
      REAL(KIND=8), ALLOCATABLE :: u(:,:,:), h(:,:), unCondU(:,:), y(:),&
     &   c(:), s(:), err(:)

      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(h(ls%sD+1,ls%sD), u(dof,nNo,ls%sD+1), unCondU(dof,nNo),  &
     &   y(ls%sD), c(ls%sD), s(ls%sD), err(ls%sD+1))

      time   = FSILS_CPUT()
      ls%suc = .FALSE.

      eps = 0D0
      X   = 0D0
      DO l=1, ls%mItr
         IF (l .EQ. 1) THEN
            u(:,:,1) = R
         ELSE
            ls%itr = ls%itr + 1
            CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof,      &
     &         Val, X, u(:,:,1))
            CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, X, u(:,:,1))
            u(:,:,1) = R - u(:,:,1)
         END IF
         IF (ANY(lhs%face%coupledFlag)) THEN
            unCondU = u(:,:,1)
            CALL ADDBCMUL(lhs, BCOP_TYPE_PRE, dof, unCondU, u(:,:,1))
         END IF

         err(1) = FSILS_NORMV(dof, mynNo, lhs%commu, u(:,:,1))
         IF (l .EQ. 1) THEN
            eps       = err(1)
            IF (eps .LE. ls%absTol) THEN
               ls%callD = EPSILON(ls%callD)
               ls%dB    = 0D0
               RETURN
            END IF
            ls%iNorm  = eps
            ls%fNorm  = eps
            eps       = MAX(ls%absTol,ls%relTol*eps)
         END IF
         ls%dB = ls%fNorm
         u(:,:,1) = u(:,:,1)/err(1)
         DO i=1, ls%sD
            ls%itr = ls%itr + 1
            CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof,      &
     &         Val, u(:,:,i), u(:,:,i+1))
            CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, u(:,:,i), u(:,:,i+1))
            IF (ANY(lhs%face%coupledFlag)) THEN
               unCondU = u(:,:,i+1)
               CALL ADDBCMUL(lhs,BCOP_TYPE_PRE,dof, unCondU, u(:,:,i+1))
            END IF
            DO j=1, i+1
               h(j,i) = FSILS_NCDOTV(dof, mynno, u(:,:,j), u(:,:,i+1))
            END DO
            CALL FSILS_BCASTV(i+1, h(:,i), lhs%commu)

            DO j=1, i
               CALL OMPSUMV(dof, nNo, -h(j,i), u(:,:,i+1), u(:,:,j))
               !u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
               h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i)
            END DO
            h(i+1,i) = SQRT(ABS(h(i+1,i)))

            CALL OMPMULV(dof, nNo, 1D0/h(i+1,i), u(:,:,i+1))
            !u(:,:,i+1) = u(:,:,i+1)/h(i+1,i)
            DO j=1, i-1
               tmp      =  c(j)*h(j,i) + s(j)*h(j+1,i)
               h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
               h(j,i)   =  tmp
            END DO
            tmp      = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i)     = h(i,i)/tmp
            s(i)     = h(i+1,i)/tmp
            h(i,i)   = tmp
            h(i+1,i) = 0D0
            err(i+1) = -s(i)*err(i)
            err(i)   =  c(i)*err(i)
            IF (ABS(err(i+1)) .LT. eps) THEN
               ls%suc = .TRUE.
               EXIT
            END IF
         END DO
         IF (i .GT. ls%sD) i = ls%sD

         y = err(1:i)
         DO j=i, 1, -1
            DO k=j+1, i
               y(j) = y(j) - h(j,k)*y(k)
            END DO
            y(j) = y(j)/h(j,j)
         END DO
 
         DO j=1, i
            CALL OMPSUMV(dof, nNo, y(j), X, u(:,:,j))
            !X = X + u(:,:,j)*y(j)
         END DO
         ls%fNorm = ABS(err(i+1))
         IF (ls%suc) EXIT
      END DO

      ls%callD = FSILS_CPUT() - time + ls%callD
      ls%dB    = 1D1*LOG(ls%fNorm/ls%dB)

      RETURN
      END SUBROUTINE GMRES

!====================================================================
      
      SUBROUTINE GMRESS(lhs, ls, Val, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      REAL(KIND=8), INTENT(IN) :: Val(lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
 
      INTEGER nNo, mynNo, i, j, k, l
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMS, FSILS_DOTS, FSILS_NCDOTS
      REAL(KIND=8) eps, tmp
      REAL(KIND=8), ALLOCATABLE :: u(:,:), h(:,:), X(:), y(:), c(:),    &
     &   s(:), err(:)

      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(h(ls%sD+1,ls%sD), u(nNo,ls%sD+1), X(nNo), y(ls%sD),      &
     &   c(ls%sD), s(ls%sD), err(ls%sD+1))
       
      ls%callD  = FSILS_CPUT()
      ls%suc    = .FALSE.
      eps       = FSILS_NORMS(mynNo, lhs%commu, R)
      ls%iNorm  = eps
      ls%fNorm  = eps
      eps       = MAX(ls%absTol,ls%relTol*eps)
      ls%itr    = 0
      X         = 0D0

      IF (ls%iNorm .LE. ls%absTol) THEN
         ls%callD = EPSILON(ls%callD)
         ls%dB    = 0D0
         RETURN
      END IF
      DO l=1, ls%mItr
         ls%dB = ls%fNorm
         ls%itr = ls%itr + 1
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, Val,X,u(:,1))
         
         u(:,1) = R - u(:,1)
         err(1) = FSILS_NORMS(mynNo, lhs%commu, u(:,1))
         u(:,1) = u(:,1)/err(1)
         DO i=1, ls%sD
            ls%itr = ls%itr + 1
            CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, Val,      &
     &         u(:,i), u(:,i+1))
            
            DO j=1, i+1
               h(j,i) = FSILS_NCDOTS(mynno, u(:,j), u(:,i+1))
            END DO
            CALL FSILS_BCASTV(i+1, h(:,i), lhs%commu)

            DO j=1, i
               CALL OMPSUMS(nNo, -h(j,i), u(:,i+1), u(:,j))
               !u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
               h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i)
            END DO
            h(i+1,i) = SQRT(ABS(h(i+1,i)))

            CALL OMPMULS(nNo, 1D0/h(i+1,i), u(:,i+1))
            !u(:,i+1) = u(:,i+1)/h(i+1,i)
            DO j=1, i-1
               tmp      =  c(j)*h(j,i) + s(j)*h(j+1,i)
               h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
               h(j,i)   =  tmp
            END DO
            tmp      = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i)     = h(i,i)/tmp
            s(i)     = h(i+1,i)/tmp
            h(i,i)   = tmp
            h(i+1,i) = 0D0
            err(i+1) = -s(i)*err(i)
            err(i)   =  c(i)*err(i)
            IF (ABS(err(i+1)) .LT. eps) THEN
               ls%suc = .TRUE.
               EXIT
            END IF
         END DO
         IF (i .GT. ls%sD) i = ls%sD

         y = err(1:i)
         DO j=i, 1, -1
            DO k=j+1, i
               y(j) = y(j) - h(j,k)*y(k)
            END DO
            y(j) = y(j)/h(j,j)
         END DO
 
         DO j=1, i
            CALL OMPSUMS(nNo, y(j), X, u(:,j))
            !X = X + u(:,j)*y(j)
         END DO
         ls%fNorm = ABS(err(i+1))
         IF (ls%suc) EXIT
      END DO
      R = X
      ls%callD = FSILS_CPUT() - ls%callD
      ls%dB    = 1D1*LOG(ls%fNorm/ls%dB)

      RETURN
      END SUBROUTINE GMRESS
!====================================================================
      
      SUBROUTINE GMRESV(lhs, ls, dof, Val, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: Val(dof*dof,lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(dof,lhs%nNo)
 
      LOGICAL flag
      INTEGER nNo, mynNo, i, j, k, l
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMV, FSILS_DOTV, FSILS_NCDOTV
      REAL(KIND=8) eps, tmp
      REAL(KIND=8), ALLOCATABLE :: u(:,:,:), h(:,:), X(:,:), y(:), c(:),&
     &   s(:), err(:), unCondU(:,:)

      flag = .FALSE.
      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(h(ls%sD+1,ls%sD), u(dof,nNo,ls%sD+1), X(dof,nNo),        &
     &   y(ls%sD), c(ls%sD), s(ls%sD), err(ls%sD+1), unCondU(dof,nNo))
       
      ls%callD  = FSILS_CPUT()
      ls%suc    = .FALSE.
      eps       = FSILS_NORMV(dof, mynNo, lhs%commu, R)
      ls%iNorm  = eps
      ls%fNorm  = eps
      eps       = MAX(ls%absTol,ls%relTol*eps)
      ls%itr    = 0
      X         = 0D0
      
      CALL BCPRE

      IF (ls%iNorm .LE. ls%absTol) THEN
         ls%callD = EPSILON(ls%callD)
         ls%dB    = 0D0
         RETURN
      END IF

      DO l=1, ls%mItr
         ls%dB = ls%fNorm
         ls%itr = ls%itr + 1
         CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, Val, X, &
     &      u(:,:,1))
         CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, X, u(:,:,1))
 
         u(:,:,1) = R - u(:,:,1)
         IF (ANY(lhs%face%coupledFlag).AND.flag) THEN
            unCondU = u(:,:,1)
            CALL ADDBCMUL(lhs, BCOP_TYPE_PRE, dof, unCondU, u(:,:,1))
         END IF
         err(1)   = FSILS_NORMV(dof, mynNo, lhs%commu, u(:,:,1))
         u(:,:,1) = u(:,:,1)/err(1)
         DO i=1, ls%sD
            ls%itr = ls%itr + 1
            CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, Val, &
     &         u(:,:,i), u(:,:,i+1))
            CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, u(:,:,i), u(:,:,i+1))

            IF (ANY(lhs%face%coupledFlag).AND.flag) THEN
               unCondU = u(:,:,i+1)
               CALL ADDBCMUL(lhs,BCOP_TYPE_PRE,dof, unCondU, u(:,:,i+1))
            END IF
            DO j=1, i+1
               h(j,i) = FSILS_NCDOTV(dof, mynno, u(:,:,j), u(:,:,i+1))
            END DO
            CALL FSILS_BCASTV(i+1, h(:,i), lhs%commu)
            
            DO j=1, i
               CALL OMPSUMV(dof, nNo, -h(j,i), u(:,:,i+1), u(:,:,j))
!              u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
               h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i)
            END DO
            h(i+1,i) = SQRT(ABS(h(i+1,i)))

            CALL OMPMULV(dof, nNo, 1D0/h(i+1,i), u(:,:,i+1))
            !u(:,:,i+1) = u(:,:,i+1)/h(i+1,i)
            DO j=1, i-1
               tmp      =  c(j)*h(j,i) + s(j)*h(j+1,i)
               h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
               h(j,i)   =  tmp
            END DO
            tmp      = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i)     = h(i,i)/tmp
            s(i)     = h(i+1,i)/tmp
            h(i,i)   = tmp
            h(i+1,i) = 0D0
            err(i+1) = -s(i)*err(i)
            err(i)   =  c(i)*err(i)
            IF (ABS(err(i+1)) .LT. eps) THEN
               ls%suc = .TRUE.
               EXIT
            END IF
         END DO
         IF (i .GT. ls%sD) i = ls%sD

         y = err(1:i)
         DO j=i, 1, -1
            DO k=j+1, i
               y(j) = y(j) - h(j,k)*y(k)
            END DO
            y(j) = y(j)/h(j,j)
         END DO
 
         DO j=1, i
            CALL OMPSUMV(dof, nNo, y(j), X, u(:,:,j))
            !X = X + u(:,:,j)*y(j)
         END DO
         ls%fNorm = ABS(err(i+1))
         IF (ls%suc) EXIT
      END DO
      R = X
      ls%callD = FSILS_CPUT() - ls%callD
      ls%dB    = 1D1*LOG(ls%fNorm/ls%dB)

      RETURN

      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE BCPRE

      IMPLICIT NONE

      INTEGER faIn, i, a, Ac, nsd
      REAL(KIND=8) FSILS_NORMV
      REAL(KIND=8), ALLOCATABLE :: v(:,:)

      nsd = dof -  1
      ALLOCATE(v(nsd,nNo))
      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            IF (lhs%face(faIn)%sharedFlag) THEN
               v = 0D0
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     v(i,Ac) = lhs%face(faIn)%valM(i,a)
                  END DO
               END DO
               lhs%face(faIn)%nS = FSILS_NORMV(nsd, mynNo, lhs%commu,   &
     &            v)**2D0
            ELSE
               lhs%face(faIn)%nS = 0D0
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     lhs%face(faIn)%nS = lhs%face(faIn)%nS +            &
     &                  lhs%face(faIn)%valM(i,a)**2D0
                  END DO
               END DO
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE BCPRE

      END SUBROUTINE GMRESV


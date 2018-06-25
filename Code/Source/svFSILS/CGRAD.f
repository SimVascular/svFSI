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
!     Conjugate-gradient algorithm for scaler, vector and Schur
!     complement cases.
!--------------------------------------------------------------------

      SUBROUTINE CGRAD_SCHUR(lhs, ls, dof, D, G, L, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: D(dof,lhs%nnz), G(dof,lhs%nnz),       &
     &   L(lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
 
      INTEGER nNo, mynNo, i
      REAL(KIND=8) errO, err, alpha, eps, time
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMS, FSILS_DOTS
      REAL(KIND=8), ALLOCATABLE :: X(:), P(:), SP(:), DGP(:), GP(:,:),  &
     &   unCondU(:,:)
      
      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(X(nNo), P(nNo), SP(nNo), DGP(nNo), GP(dof,nNo),          &
     &   unCondU(dof,nNo))
 
      time     = FSILS_CPUT()
      ls%suc   = .FALSE.
      ls%iNorm = FSILS_NORMS(mynNo, lhs%commu, R)
      eps      = MAX(ls%absTol,ls%relTol*ls%iNorm)**2D0
      errO     = ls%iNorm*ls%iNorm
      err      = errO
      X        = 0D0
      P        = R

      DO i=1, ls%mItr
         IF (err .LT. eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         errO = err
         CALL FSILS_SPARMULSV(lhs, lhs%rowPtr, lhs%colPtr, dof, G, P,GP)
         IF (ANY(lhs%face%coupledFlag)) THEN
            unCondU = GP
            CALL ADDBCMUL(lhs, BCOP_TYPE_PRE, dof, unCondU, GP)
         END IF
         CALL FSILS_SPARMULVS(lhs, lhs%rowPtr, lhs%colPtr,dof, D,GP,DGP)
         
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, L, P, SP)
       
         CALL OMPSUMS(nNo, -1D0, SP, DGP)
         !SP    = SP - DGP
         alpha = errO/FSILS_DOTS(mynNo, lhs%commu, P, SP)
         CALL OMPSUMS(nNo, alpha, X, P)
         !X     = X + alpha*P
         CALL OMPSUMS(nNo, -alpha, R, SP)
         !R     = R - alpha*SP
         err   = FSILS_NORMS(mynNo, lhs%commu, R)
         err   = err*err
         CALL OMPSUMS(nNo, errO/err, P, R)
         CALL OMPMULS(nNo, err/errO, P)
         !P     = R + err/errO*P
      END DO
      R        = X
      ls%fNorm = SQRT(err)
      ls%callD = FSILS_CPUT() - time + ls%callD
      ls%itr   = ls%itr + i - 1
      IF (errO .LT. EPSILON(errO)) THEN
         ls%dB = 0D0
      ELSE
         ls%dB = 5D0*LOG(err/errO)
      END IF
      
      DEALLOCATE(X, P, SP, DGP, GP)
      
      RETURN
      END SUBROUTINE CGRAD_SCHUR

!====================================================================
       
      SUBROUTINE CGRADS(lhs, ls, K, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      REAL(KIND=8), INTENT(IN) :: K(lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
     
      INTEGER nNo, mynNo, i
      REAL(KIND=8) errO, err, alpha, eps
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMS, FSILS_DOTS
      REAL(KIND=8), ALLOCATABLE :: P(:), KP(:), X(:)
      
      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(P(nNo), KP(nNo), X(nNo))
      
      ls%callD = FSILS_CPUT()
      ls%suc   = .FALSE.
      ls%iNorm = FSILS_NORMS(mynNo, lhs%commu, R)
      eps      = MAX(ls%absTol,ls%relTol*ls%iNorm)**2D0
      errO     = ls%iNorm*ls%iNorm
      err      = errO
      X        = 0D0
      P        = R

      DO i=1, ls%mItr
         IF (err .LT. eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         errO = err
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, K, P, KP)
         alpha = errO/FSILS_DOTS(mynNo, lhs%commu, P, KP)
         CALL OMPSUMS(nNo, alpha, X, P)
         !X     = X + alpha*P
         CALL OMPSUMS(nNo, -alpha, R, KP)
         !R     = R - alpha*KP
         err   = FSILS_NORMS(mynNo, lhs%commu, R)
         err   = err*err
         CALL OMPSUMS(nNo, errO/err, P, R)
         CALL OMPMULS(nNo, err/errO, P)
         !P = R + err/errO*P
      END DO

      R        = X
      ls%itr   = i - 1
      ls%fNorm = SQRT(err)
      ls%callD = FSILS_CPUT() - ls%callD
      IF (errO .LT. EPSILON(errO)) THEN
         ls%dB = 0D0
      ELSE
         ls%dB = 5D0*LOG(err/errO)
      END IF

      RETURN
      END SUBROUTINE CGRADS

!====================================================================
       
      SUBROUTINE CGRADV(lhs, ls, dof, K, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: K(dof*dof,lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(dof,lhs%nNo)
     
      INTEGER nNo, mynNo, i
      REAL(KIND=8) errO, err, alpha, eps
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMV, FSILS_DOTV
      REAL(KIND=8), ALLOCATABLE :: P(:,:), KP(:,:), X(:,:)
      
      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(P(dof,nNo), KP(dof,nNo), X(dof,nNo))
      
      ls%callD = FSILS_CPUT()
      ls%suc   = .FALSE.
      ls%iNorm = FSILS_NORMV(dof, mynNo, lhs%commu, R)
      eps      = MAX(ls%absTol,ls%relTol*ls%iNorm)**2D0
      errO     = ls%iNorm*ls%iNorm
      err      = errO
      X        = 0D0
      P        = R

      DO i=1, ls%mItr
         IF (err .LT. eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         errO = err
         CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, K, P,KP)
         alpha = errO/FSILS_DOTV(dof, mynNo, lhs%commu, P, KP)
         CALL OMPSUMV(dof, nNo, alpha, X, P)
         !X     = X + alpha*P
         CALL OMPSUMV(dof, nNo, -alpha, R, KP)
         !R     = R - alpha*KP
         err   = FSILS_NORMV(dof, mynNo, lhs%commu, R)
         err   = err*err
         CALL OMPSUMV(dof, nNo, errO/err, P, R)
         CALL OMPMULV(dof, nNo, err/errO, P)
         !P = R + err/errO*P
      END DO

      R        = X
      ls%itr   = i - 1
      ls%fNorm = SQRT(err)
      ls%callD = FSILS_CPUT() - ls%callD
      IF (errO .LT. EPSILON(errO)) THEN
         ls%dB = 0D0
      ELSE
         ls%dB = 5D0*LOG(err/errO)
      END IF

      RETURN
      END SUBROUTINE CGRADV



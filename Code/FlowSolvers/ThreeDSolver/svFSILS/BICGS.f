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
!     Biconjugate-gradient algorithm, available for scaler and vectors
!--------------------------------------------------------------------

      SUBROUTINE BICGSS(lhs, ls, K, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      REAL(KIND=8), INTENT(IN) :: K(lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
     
      INTEGER nNo, mynNo, i
      REAL(KIND=8) errO, err, alpha, beta, rho, rhoO, omega, eps
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMS, FSILS_DOTS
      REAL(KIND=8), ALLOCATABLE :: P(:), Rh(:), X(:), V(:), S(:), T(:)

      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(P(nNo), Rh(nNo), X(nNo), V(nNo), S(nNo), T(nNo))
      
      ls%callD = FSILS_CPUT()
      ls%suc   = .FALSE.
      err      = FSILS_NORMS(mynNo, lhs%commu, R)
      ls%iNorm = err
      errO     = err
      eps      = MAX(ls%absTol,ls%relTol*err)
      rho      = err*err
      beta     = rho
      X        = 0D0
      P        = R
      Rh       = R

      DO i=1, ls%mItr
         IF (err .LT. eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, K, P, V)
         alpha = rho/FSILS_DOTS(mynNo, lhs%commu, Rh, V)
         S     = R - alpha*V
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, K, S, T)
         omega = FSILS_NORMS(mynNo, lhs%commu, T)
         omega = FSILS_DOTS(mynNo, lhs%commu, T, S)/(omega*omega)
         X     = X + alpha*P + omega*S
         R     = S - omega*T
         errO  = err
         err   = FSILS_NORMS(mynNo, lhs%commu, R)
         rhoO  = rho
         rho   = FSILS_DOTS(mynNo, lhs%commu, R, Rh)
         beta  = rho*alpha/(rhoO*omega)
         P     = R + beta*(P - omega*V)
      END DO

      R        = X
      ls%itr   = i - 1
      ls%fNorm = err
      ls%callD = FSILS_CPUT() - ls%callD
      IF (errO .LT. EPSILON(errO)) THEN
         ls%dB = 0D0
      ELSE
         ls%dB = 10D0*LOG(err/errO)
      END IF

      RETURN
      END SUBROUTINE BICGSS
!====================================================================
       
      SUBROUTINE BICGSV(lhs, ls, dof, K, R)
      
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: K(dof*dof,lhs%nnz)
      REAL(KIND=8), INTENT(INOUT) :: R(dof,lhs%nNo)
     
      INTEGER nNo, mynNo, i
      REAL(KIND=8) errO, err, alpha, beta, rho, rhoO, omega, eps
      REAL(KIND=8) FSILS_CPUT, FSILS_NORMV, FSILS_DOTV
      REAL(KIND=8), ALLOCATABLE :: P(:,:), Rh(:,:), X(:,:), V(:,:),     &
     &   S(:,:), T(:,:)

      nNo = lhs%nNo
      mynNo = lhs%mynNo

      ALLOCATE(P(dof,nNo), Rh(dof,nNo), X(dof,nNo), V(dof,nNo),         &
     &   S(dof,nNo), T(dof,nNo))
      
      ls%callD = FSILS_CPUT()
      ls%suc   = .FALSE.
      err      = FSILS_NORMV(dof, mynNo, lhs%commu, R)
      errO     = err
      ls%iNorm = err
      eps      = MAX(ls%absTol,ls%relTol*err)
      rho      = err*err
      beta     = rho
      X        = 0D0
      P        = R
      Rh       = R

      DO i=1, ls%mItr
         IF (err .LT. eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, K, P, V)
         alpha = rho/FSILS_DOTV(dof, mynNo, lhs%commu, Rh, V)
         S     = R - alpha*V
         CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, K, S, T)
         omega = FSILS_NORMV(dof, mynNo, lhs%commu, T)
         omega = FSILS_DOTV(dof, mynNo, lhs%commu, T, S)/(omega*omega)
         X     = X + alpha*P + omega*S
         R     = S - omega*T
         errO  = err
         err   = FSILS_NORMV(dof, mynNo, lhs%commu, R)
         rhoO  = rho
         rho   = FSILS_DOTV(dof, mynNo, lhs%commu, R, Rh)
         beta  = rho*alpha/(rhoO*omega)
         P     = R + beta*(P - omega*V)
      END DO

      R        = X
      ls%itr   = i - 1
      ls%fNorm = err
      ls%callD = FSILS_CPUT() - ls%callD
      IF (errO .LT. EPSILON(errO)) THEN
         ls%dB = 0D0
      ELSE
         ls%dB = 10D0*LOG(err/errO)
      END IF

      RETURN
      END SUBROUTINE BICGSV

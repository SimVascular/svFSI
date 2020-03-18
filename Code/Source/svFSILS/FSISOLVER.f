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
!     This routine is mainley intended for solving incompressible NS or
!     FSI equations with a form of AU=R, in which A = [K D;-G L] and
!     G = -D^t
!--------------------------------------------------------------------

      SUBROUTINE FSISOLVER(lhs, ls, dof, Val, R, isS)
      INCLUDE "FSILS_STD.h"
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(IN) :: Val(dof*dof,lhs%nnz)
      REAL(KIND=LSRP), INTENT(INOUT) :: R(dof,lhs%nNo)
      LOGICAL, INTENT(IN) :: isS(lhs%nNo)

      LOGICAL flag, GE
      INTEGER(KIND=LSIP) nNo, mynNo, i, j, k, n, l
      REAL(KIND=LSRP) FSILS_CPUT, FSILS_NORMV, FSILS_DOTV, FSILS_NCDOTV
      REAL(KIND=LSRP) eps, temp

      REAL(KIND=LSRP), ALLOCATABLE :: u(:,:,:), h(:,:), X(:,:), y(:),
     2   c(:), s(:), err(:), unCondU(:,:), v(:,:,:,:), A(:,:), b(:),
     3   xb(:), oldxb(:), tmpU(:,:,:), tmp(:)

      flag  = .FALSE.
      nNo   = lhs%nNo
      mynNo = lhs%mynNo
      n     = ls%sD
      ALLOCATE(h(ls%sD+1,ls%sD), u(dof,nNo,ls%sD+1), X(dof,nNo),        &
     &   y(ls%sD), c(ls%sD), s(ls%sD), err(ls%sD+1), unCondU(dof,nNo),  &
     &   A(2*n,2*n), v(dof,nNo,n,2), b(2*n), xb(2*n), tmpU(dof,nNo,2),  &
     &   tmp(4*n+1), oldxb(2*n))

      ls%callD  = FSILS_CPUT()
      ls%suc    = .FALSE.
      eps       = FSILS_NORMV(dof, mynNo, lhs%commu, R)
      ls%iNorm  = eps
      ls%fNorm  = eps
      eps       = MAX(ls%absTol,ls%relTol*eps)
      ls%itr    = 0
      X         = 0._LSRP
      h         = 0._LSRP

      CALL BCPRE

      IF (ls%iNorm .LE. ls%absTol) THEN
         ls%callD = EPSILON(ls%callD)
         ls%dB    = 0._LSRP
         RETURN
      END IF

      err(1)   = FSILS_NORMV(dof, mynNo, lhs%commu, R)
      u(:,:,1) = R/err(1)
      i = 0
      DO
         i     = i + 1
         ls%dB = ls%fNorm
         CALL SEPARATE(isS, tmpu(:,:,1), tmpu(:,:,2), u(:,:,i))

         l = 0
         DO j=1, 2
            CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, dof, Val, &
     &         tmpu(:,:,j), v(:,:,i,j))
            CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, tmpu(:,:,j),         &
     &         v(:,:,i,j))
            IF (ANY(lhs%face%coupledFlag)) THEN
               unCondU = v(:,:,i,j)
               CALL ADDBCMUL(lhs,BCOP_TYPE_PRE, dof,unCondU, v(:,:,i,j))
            END IF

            DO k=1, i
               l = l + 1
               tmp(l) = FSILS_NCDOTV(dof, mynNo, v(:,:,k,j), v(:,:,i,j))
            END DO
         END DO
         DO k=1, i
            l = l + 1
            tmp(l) = FSILS_NCDOTV(dof, mynNo, v(:,:,k,1), v(:,:,i,2))
         END DO
         DO k=1, i-1
            l = l + 1
            tmp(l) = FSILS_NCDOTV(dof, mynNo, v(:,:,i,1), v(:,:,k,2))
         END DO
         DO j=1, 2
            l = l + 1
            tmp(l) = FSILS_NCDOTV(dof, mynNo, R, v(:,:,i,j))
         END DO
         u(:,:,i+1) = v(:,:,i,1) + v(:,:,i,2)

!     Computing H
         CALL FSILS_BCASTV(l, tmp, lhs%commu)
         l = 0
         DO k=1, i
            l = l + 1
            A(2*k-1,2*i-1) = tmp(l)
            A(2*i-1,2*k-1) = tmp(l)
         END DO
         DO k=1, i
            l = l + 1
            A(2*k,2*i) = tmp(l)
            A(2*i,2*k) = tmp(l)
         END DO
         DO k=1, i
            l = l + 1
            A(2*k-1,2*i) = tmp(l)
            A(2*i,2*k-1) = tmp(l)
         END DO
         DO k=1, i-1
            l = l + 1
            A(2*k,2*i-1) = tmp(l)
            A(2*i-1,2*k) = tmp(l)
         END DO
         DO j=1, 2
            l = l + 1
            b(2*i-2+j) = tmp(l)
         END DO

         xb = b
         IF (GE(2*n, 2*i, A, xB)) THEN
            oldxB = xB
         ELSE
            IF(lhs%commu%masF) PRINT *,"FSILS: Singular matrix detected"
            i  = i - 1
            xB = oldxB
            EXIT
         END IF
         ls%fNorm = ABS(ls%iNorm**2._LSRP - SUM(xB(1:2*i)*B(1:2*i)))
c         print *, sqrt(ls%fnorm)/ls%inorm

         IF(ls%fNorm .LT. eps*eps) THEN
            ls%suc = .TRUE.
            EXIT
         END IF
         IF (i .GE. ls%sD) EXIT

         DO j=1, i+1
         ! You may simplify this using A
            h(j,i) = FSILS_NCDOTV(dof, mynno, u(:,:,j), u(:,:,i+1))
         END DO
            CALL FSILS_BCASTV(i+1, h(:,i), lhs%commu)

         DO j=1, i
            CALL OMPSUMV(dof, nNo, -h(j,i), u(:,:,i+1), u(:,:,j))
!           u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
            h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i)
         END DO
         h(i+1,i) = SQRT(ABS(h(i+1,i)))

         CALL OMPMULV(dof, nNo, 1._LSRP/h(i+1,i), u(:,:,i+1))
         !u(:,:,i+1) = u(:,:,i+1)/h(i+1,i)
         DO j=1, i-1
            temp     =  c(j)*h(j,i) + s(j)*h(j+1,i)
            h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
            h(j,i)   =  temp
         END DO
         temp     = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
         c(i)     = h(i,i)/temp
         s(i)     = h(i+1,i)/temp
         h(i,i)   = temp
         h(i+1,i) = 0._LSRP
         err(i+1) = -s(i)*err(i)
         err(i)   =  c(i)*err(i)
         IF (ABS(err(i+1)) .LT. eps) THEN
            print *, "this is strange!!"
            stop
            ls%suc = .TRUE.
            EXIT
         END IF
      END DO

      R = 0._LSRP
      DO j=1, i
         CALL SEPARATE(isS, tmpu(:,:,1), tmpu(:,:,2), u(:,:,j))
         R = R + xB(2*j-1)*tmpu(:,:,1) + xB(2*j)*tmpu(:,:,2)
      END DO
      ls%itr   = i
      ls%callD = FSILS_CPUT() - ls%callD
      ls%dB    = 5._LSRP*LOG(ls%fNorm/ls%dB)
      ls%fNorm = SQRT(ls%fNorm)

      RETURN
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE SEPARATE(isS, Rf, Rs, Rfi, Rsi)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: isS(nNo)
      REAL(KIND=LSRP), INTENT(OUT) :: Rf(dof,nNo), Rs(dof,nNo)
      REAL(KIND=LSRP), INTENT(IN) :: Rfi(dof,nNo)
      REAL(KIND=LSRP), INTENT(IN), OPTIONAL :: Rsi(dof,nNo)

      INTEGER(KIND=LSIP) Ac

      IF (PRESENT(Rsi)) THEN
         DO Ac=1, nNo
            IF (isS(Ac)) THEN
               Rf(:,Ac) = 0._LSRP
               Rs(:,Ac) = Rfi(:,Ac) + Rsi(:,Ac)
            ELSE
               Rf(:,Ac) = Rfi(:,Ac) + Rsi(:,Ac)
               Rs(:,Ac) = 0._LSRP
            END IF
         END DO
      ELSE
         DO Ac=1, nNo
            IF (isS(Ac)) THEN
               Rf(:,Ac) = 0._LSRP
               Rs(:,Ac) = Rfi(:,Ac)
            ELSE
               Rf(:,Ac) = Rfi(:,Ac)
               Rs(:,Ac) = 0._LSRP
            END IF
         END DO
      END IF

      RETURN
      END SUBROUTINE SEPARATE
!--------------------------------------------------------------------
      SUBROUTINE BCPRE
      IMPLICIT NONE
      INTEGER(KIND=LSIP) faIn, i, a, Ac, nsd
      REAL(KIND=LSRP) FSILS_NORMV
      REAL(KIND=LSRP), ALLOCATABLE :: v(:,:)

      nsd = dof - 1
      ALLOCATE(v(nsd,nNo))
      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            IF (lhs%face(faIn)%sharedFlag) THEN
               v = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     v(i,Ac) = lhs%face(faIn)%valM(i,a)
                  END DO
               END DO
               lhs%face(faIn)%nS = FSILS_NORMV(nsd, mynNo, lhs%commu,   &
     &            v)**2._LSRP
            ELSE
               lhs%face(faIn)%nS = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     lhs%face(faIn)%nS = lhs%face(faIn)%nS +            &
     &                  lhs%face(faIn)%valM(i,a)**2._LSRP
                  END DO
               END DO
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE BCPRE
!--------------------------------------------------------------------
      END SUBROUTINE FSISOLVER
!####################################################################

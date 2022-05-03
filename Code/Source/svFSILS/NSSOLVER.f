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

      SUBROUTINE NSSOLVER(lhs, ls, dof, Val, Ri)
      INCLUDE "FSILS_STD.h"
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_lsType), INTENT(INOUT) :: ls
      INTEGER(KIND=LSIP), INTENT(IN) :: dof
      REAL(KIND=LSRP), INTENT(IN) :: Val(dof*dof,lhs%nnz)
      REAL(KIND=LSRP), INTENT(INOUT) :: Ri(dof,lhs%nNo)

      LOGICAL GE
      INTEGER(KIND=LSIP) nNo, nnz, mynNo, i, j, k, iB, iBB, nB, nsd, c
      REAL(KIND=LSRP) FSILS_CPUT, FSILS_NORMS, FSILS_NORMV, FSILS_DOTS, &
     &   FSILS_DOTV, eps, FSILS_NCDOTS, FSILS_NCDOTV

      REAL(KIND=LSRP), ALLOCATABLE :: U(:,:,:), P(:,:), MU(:,:,:),      &
     &   MP(:,:), A(:,:), tmp(:), tmpG(:), B(:), xB(:), oldxB(:),       &
     &   mK(:,:), mG(:,:), mD(:,:), mL(:), Gt(:,:), Rm(:,:), Rc(:),     &
     &   Rmi(:,:), Rci(:)

      nNo = lhs%nNo
      nnz = lhs%nnz
      mynNo = lhs%mynNo
      nsd = dof - 1
      iB = ls%RI%mItr
      nB = 2*iB
      ALLOCATE(Rm(nsd,nNo), Rc(nNo), Rmi(nsd,nNo), Rci(nNo),            &
     &   U(nsd,nNo,iB), P(nNo,iB), MU(nsd,nNo,nB), MP(nNo,nB),          &
     &   tmp(nB*nB+nB), tmpG(nB*nB+nB), A(nB,nB), B(nB), xB(nB),        &
     &   oldxB(nB))

      Rmi = Ri(1:nsd,:)
      Rci = Ri(dof,:)

      xB          = 0._LSRP
      B           = 0._LSRP
      oldxB       = 0._LSRP
      Rm          = Rmi
      Rc          = Rci
      eps         = SQRT(FSILS_NORMV(nsd,mynNo,lhs%commu,Rm)**2._LSRP   &
     &            +      FSILS_NORMS(    mynNo,lhs%commu,Rc)**2._LSRP)
      ls%RI%iNorm = eps
      ls%RI%fNorm = eps*eps
      ls%CG%callD = 0._LSRP
      ls%GM%callD = 0._LSRP
      ls%CG%itr   = 0
      ls%GM%itr   = 0
      ls%RI%callD = FSILS_CPUT()
      ls%RI%suc   = .FALSE.
      eps         = MAX(ls%RI%absTol,ls%RI%relTol*eps)

      CALL DEPART
      CALL BCPRE

      DO i=1, ls%RI%mItr
         iB  = 2*i - 1
         iBB = 2*i
         ls%RI%dB = ls%RI%fNorm

!        U  = K^-1*Rm
         CALL GMRES(lhs, ls%GM, nsd, mK, Rm, U(:,:,i))
!        P  = D*U
         CALL FSILS_SPARMULVS(lhs, lhs%rowPtr, lhs%colPtr, nsd, mD,     &
     &      U(:,:,i), P(:,i))
!        P  = Rc - P
         P(:,i) = Rc - P(:,i)
!        P  = [L + G^t*G]^-1*P
         CALL CGRAD_SCHUR(lhs, ls%CG, nsd, Gt, mG, mL, P(:,i))
!        MU1 = G*P
         CALL FSILS_SPARMULSV(lhs, lhs%rowPtr, lhs%colPtr, nsd, mG,     &
     &      P(:,i), MU(:,:,iB))
!        MU2 = Rm - G*P
         MU(:,:,iBB) = Rm - MU(:,:,iB)
!        U  = K^-1*[Rm - G*P]
         CALL GMRES(lhs, ls%GM, nsd, mK, MU(:,:,iBB), U(:,:,i))
!        MU2 = K*U
         CALL FSILS_SPARMULVV(lhs, lhs%rowPtr, lhs%colPtr, nsd, mK,     &
     &      U(:,:,i), MU(:,:,iBB))
!        The contribution of coupled BCs is added to the matrix-vector
!        product operation. These coupled BCs contribute through a "resistance"
!        that relates velocity (integrated to give flow rate) to pressure. 
!        See Moghadam et al. 2013 eq. 27.
         CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, nsd, U(:,:,i), MU(:,:,iBB))
!        MP1 = L*P
         CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, mL, P(:,i),  &
     &      MP(:,iB))
!        MP2 = D*U
         CALL FSILS_SPARMULVS(lhs, lhs%rowPtr, lhs%colPtr, nsd, mD,     &
     &      U(:,:,i), MP(:,iBB))

         c = 0
         DO k=iB, iBB
            DO j=1, k
               c = c + 1
               tmp(c) = FSILS_NCDOTV(nsd, mynNo, MU(:,:,j), MU(:,:,k))  &
     &                + FSILS_NCDOTS(     mynNo, MP(:,j),   MP(:,k))
            END DO
            c = c + 1
            tmp(c) = FSILS_NCDOTV(nsd, mynNo, MU(:,:,k), Rmi)           &
     &             + FSILS_NCDOTS(     mynNo,  MP(:,k),  Rci)
         END DO
         IF (lhs%commu%nTasks .GT. 1) THEN
            CALL MPI_ALLREDUCE(tmp, tmpG, c, mpreal, MPI_SUM,           &
     &         lhs%commu%comm, j)
            tmp = tmpG
         END IF

         c = 0
         DO k=iB, iBB
            DO j=1, k
               c      = c + 1
               A(j,k) = tmp(c)
               A(k,j) = tmp(c)
            END DO
            c    = c + 1
            B(k) = tmp(c)
         END DO

         xB = B
         IF (GE(nB, iBB, A, xB)) THEN
            oldxB = xB
         ELSE
            IF(lhs%commu%masF) PRINT *,"FSILS: Singular matrix detected"
            xB = oldxB
            IF (i .GT. 1) THEN
               iB  = iB  - 2
               iBB = iBB - 2
            END IF
            EXIT
         END IF

         ls%RI%fNorm = ls%RI%iNorm**2._LSRP - SUM(xB(1:iBB)*B(1:iBB))
         IF(ls%RI%fNorm .LT. eps*eps) THEN
            ls%RI%suc = .TRUE.
            EXIT
         END IF

         Rm = Rmi - xB(1)*MU(:,:,1)
         Rc = Rci - xB(1)*MP(:,1)
         DO j=2, iBB
            Rm = Rm - xB(j)*MU(:,:,j)
            Rc = Rc - xB(j)*MP(:,j)
         END DO
      END DO
      IF (i .GT. ls%RI%mItr) THEN
         ls%RI%itr = ls%RI%mItr
      ELSE
         ls%RI%itr = i

         Rc = Rci - xB(1)*MP(:,1)
         DO j=2, iBB
            Rc = Rc - xB(j)*MP(:,j)
         END DO
      END IF
      ls%Resc = NINT(100._LSRP*FSILS_NORMS(mynNo, lhs%commu,            &
     &   Rc)**2._LSRP / ls%RI%fNorm, KIND=LSIP)
      ls%Resm = 100 - ls%Resc

      Rmi = xB(2)*U(:,:,1)
      Rci = xB(1)*P(:,1)
      DO i=2, ls%RI%itr
         iB  = 2*i - 1
         iBB = 2*i

         Rmi = Rmi + xB(iBB)*U(:,:,i)
         Rci = Rci + xB(iB)*P(:,i)
      END DO

      ls%RI%callD = FSILS_CPUT() - ls%RI%callD
      ls%RI%dB    = 5._LSRP*LOG(ls%RI%fNorm/ls%RI%dB)

      IF (ls%Resc.LT.0 .OR. ls%Resm.LT.0) THEN
         ls%Resc = 0
         ls%Resm = 0
         ls%RI%db = 0
         ls%RI%fNorm = 0._LSRP
         IF (lhs%commu%masF) THEN
            PRINT "(A)", "Warning: unexpected behavior in FSILS"//      &
     &        " (likely due to the ill-conditioned LHS matrix)"
         END IF
      END IF
      ls%RI%fNorm = SQRT(ls%RI%fNorm)

      Ri(1:nsd,:) = Rmi
      Ri(dof,:) = Rci

      IF (lhs%commu%masF) CALL LOGFILE

      RETURN
      CONTAINS
!--------------------------------------------------------------------
      SUBROUTINE DEPART
      IMPLICIT NONE

      INTEGER(KIND=LSIP) i, j, k, l
      REAL(KIND=LSRP), ALLOCATABLE :: tmp(:)

      ALLOCATE(mK(nsd*nsd,nnz), mG(nsd,nnz), mD(nsd,nnz), mL(nnz),      &
     &   Gt(nsd,nnz), tmp((nsd+1)*(nsd+1)))

      IF (nsd .EQ. 2) THEN
         DO i=1, nnz
            tmp = Val(:,i)

            mK(1,i) = tmp(1)
            mK(2,i) = tmp(2)
            mK(3,i) = tmp(4)
            mK(4,i) = tmp(5)

            mG(1,i) = tmp(3)
            mG(2,i) = tmp(6)

            mD(1,i) = tmp(7)
            mD(2,i) = tmp(8)

            mL(i)   = tmp(9)
         END DO
      ELSE IF(nsd .EQ. 3) THEN
         DO i=1, nnz
            tmp = Val(:,i)

            mK(1,i) = tmp(1)
            mK(2,i) = tmp(2)
            mK(3,i) = tmp(3)
            mK(4,i) = tmp(5)
            mK(5,i) = tmp(6)
            mK(6,i) = tmp(7)
            mK(7,i) = tmp(9)
            mK(8,i) = tmp(10)
            mK(9,i) = tmp(11)

            mG(1,i) = tmp(4)
            mG(2,i) = tmp(8)
            mG(3,i) = tmp(12)

            mD(1,i) = tmp(13)
            mD(2,i) = tmp(14)
            mD(3,i) = tmp(15)

            mL(i)   = tmp(16)
         END DO
      ELSE
         PRINT *, "FSILS: Not defined nsd for DEPART", nsd
         STOP "FSILS: FATAL ERROR"
      END IF

      DO i=1, nNo
         Do j=lhs%rowPtr(1,i), lhs%rowPtr(2,i)
            k = lhs%colPtr(j)
            DO l=lhs%rowPtr(1,k), lhs%rowPtr(2,k)
               IF (lhs%colPtr(l) .EQ. i) THEN
                  Gt(:,l) = -mG(:,j)
                  EXIT
               END IF
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE DEPART
!--------------------------------------------------------------------
      SUBROUTINE BCPRE
      IMPLICIT NONE

      INTEGER(KIND=LSIP) faIn, i, a, Ac
      REAL(KIND=LSRP) FSILS_NORMV
      REAL(KIND=LSRP), ALLOCATABLE :: v(:,:)

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
      SUBROUTINE LOGFILE
      IMPLICIT NONE

      LOGICAL flag
      INTEGER(KIND=LSIP) fid, i, j
      REAL(KIND=LSRP) rtmp
      CHARACTER(LEN=*), PARAMETER :: fName = 'FSILS_NS.log'

      INQUIRE(FILE=fName, EXIST=flag)
      fid = 11232
      OPEN(fid, FILE=fName, POSITION='APPEND')
      IF (.NOT.flag) THEN
         i = 0
         DO j=1, lhs%nFaces
            IF (lhs%face(j)%coupledFlag) i = i + 1
         END DO
         WRITE(fid,*) lhs%gnNo, lhs%commu%nTasks, i
      END IF

      i = 0
      IF (ls%RI%suc) i = i + 100
      IF (ls%GM%suc) i = i + 10
      IF (ls%CG%suc) i = i + 1

      IF (ls%RI%CallD .LT. TINY(ls%RI%CallD)) ls%RI%CallD = 1.E-99_LSRP
      rtmp = (ls%RI%CallD-ls%GM%CallD-ls%CG%CallD)/ls%RI%CallD
      WRITE(fid,"(I4.3,I3,I4,I5,3I4,3ES9.2E2,3I4)")                     &
     &   i, ls%RI%itr, ls%GM%itr, ls%CG%itr,                            &
     &   NINT(100._LSRP*rtmp, KIND=LSIP),                               &
     &   NINT(100._LSRP*ls%GM%callD/ls%RI%CallD, KIND=LSIP),            &
     &   NINT(100._LSRP*ls%CG%callD/ls%RI%CallD, KIND=LSIP),            &
     &   ls%RI%iNorm, ls%RI%fNorm/ls%RI%iNorm, ls%RI%CallD,             &
     &   ls%Resm, ls%Resc, NINT(ls%RI%dB, KIND=LSIP)

      CLOSE(fid)

      RETURN
      END SUBROUTINE LOGFILE
!--------------------------------------------------------------------
      END SUBROUTINE NSSOLVER
!####################################################################

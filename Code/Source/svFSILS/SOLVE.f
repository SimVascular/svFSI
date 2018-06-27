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
!     In this routine, the appropriate LS algorithm is called and
!     the solution is returned
!--------------------------------------------------------------------

      SUBROUTINE FSILS_SOLVE (lhs, ls, dof, Ri, Val, isS, incL, res)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(INOUT) :: Ri(dof,lhs%nNo)
      REAL(KIND=8), INTENT(INOUT) :: Val(dof*dof,lhs%nnz)
      LOGICAL, INTENT(IN), OPTIONAL :: isS(lhs%nNo)
      INTEGER, INTENT(IN), OPTIONAL :: incL(lhs%nFaces)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: res(lhs%nFaces)

      LOGICAL flag
      INTEGER faIn, a, nNo, nnz, nFaces
      REAL(KIND=8), ALLOCATABLE :: R(:,:), W(:,:)
      
      nNo    = lhs%nNo
      nnz    = lhs%nnz
      nFaces = lhs%nFaces

      IF (lhs%nFaces .NE. 0) THEN
         lhs%face%incFlag = .TRUE.
         IF (PRESENT(incL)) THEN
            DO faIn=1, lhs%nFaces
               IF (incL(faIn) .EQ. 0) lhs%face(faIn)%incFlag = .FALSE.
            END DO
         END IF

         flag = ANY(lhs%face%bGrp.EQ.BC_TYPE_Neu)
         IF (.NOT.PRESENT(res) .AND. flag) THEN
            PRINT *, "FSILS: res is required for Neu surfaces"
            STOP "FSILS: FATAL ERROR"
         END IF
         DO faIn=1, lhs%nFaces
            lhs%face(faIn)%coupledFlag = .FALSE.
            IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
            flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
            IF (flag .AND. res(faIn).NE.0D0) THEN
               lhs%face(faIn)%res = res(faIn)
               lhs%face(faIn)%coupledFlag = .TRUE.
            END IF
         END DO
      END IF
      
      ALLOCATE(R(dof,nNo), W(dof,nNo))
      DO a=1, nNo
         R(:,lhs%map(a)) = Ri(:,a)
      END DO

      CALL FSILS_COMMUV(lhs, dof, R)
      CALL PRECOND(lhs, lhs%rowPtr, lhs%colPtr, lhs%diagPtr,dof,Val,R,W)
 
      SELECT CASE (ls%LS_type)
         CASE (LS_TYPE_NS)
            CALL NSSOLVER(lhs, ls, dof, Val, R)
         CASE (LS_TYPE_GMRES)
            IF (dof .EQ. 1) THEN
               CALL GMRESS(lhs, ls%RI, Val, R)
            ELSE
               CALL GMRESV(lhs, ls%RI, dof, Val, R)
            END IF
         CASE (LS_TYPE_CG)
            IF (dof .EQ. 1) THEN
               CALL CGRADS(lhs, ls%RI, Val, R)
            ELSE
               CALL CGRADV(lhs, ls%RI, dof, Val, R)
            END IF
         CASE (LS_TYPE_BICGS)
            IF (dof .EQ. 1) THEN
               CALL BICGSS(lhs, ls%RI, Val, R)
            ELSE
               CALL BICGSV(lhs, ls%RI, dof, Val, R)
            END IF
         CASE DEFAULT
            PRINT *, 'FSILS: LS_type not defined'
            STOP "FSILS: FATAL ERROR"
      END SELECT
      R = R*W

      DO a=1, nNo
         Ri(:,a) = R(:,lhs%map(a))
      END DO

      DEALLOCATE(R, W)

      RETURN
      END SUBROUTINE FSILS_SOLVE


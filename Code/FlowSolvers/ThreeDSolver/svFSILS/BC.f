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
!     The structure of BC input is created here.
!--------------------------------------------------------------------


      SUBROUTINE FSILS_BC_CREATE (lhs, faIn, nNo, dof, BC_type, gNodes, &
     &   Val)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: faIn, nNo, dof
      INTEGER, INTENT(IN) :: BC_type
      INTEGER, INTENT(IN) :: gNodes(nNo)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: Val(dof,nNo)

      INTEGER a, Ac, i
      REAL(KIND=8), ALLOCATABLE :: v(:,:)

      IF (faIn .GT. lhs%nFaces) THEN
         PRINT *, "FSILS: faIn is exceeding lhs structure maximum",     &
     &      " number of face:", lhs%nFaces, " .LT. ", faIn
         STOP "FSILS: FATAL ERROR"
      END IF
      IF (faIn .LE. 0) THEN
         PRINT *, "FSILS: faIn should be greater than zero"
         STOP "FSILS: FATAL ERROR"
      END IF

      IF (lhs%face(faIn)%foC) THEN
         PRINT *, "FSILS: face is not free, you may use FSILS_BC_FREE", &
     &      " to free it"
         STOP "FSILS: FATAL ERROR"
      END IF

      lhs%face(faIn)%nNo  = nNo
      lhs%face(faIn)%dof  = dof
      lhs%face(faIn)%bGrp = BC_type

      ALLOCATE(lhs%face(faIn)%glob(nNo), lhs%face(faIn)%val(dof,nNo),   &
     &   lhs%face(faIn)%valM(dof,nNo))

      DO a=1, nNo
         Ac = lhs%map(gNodes(a))
         lhs%face(faIn)%glob(a) = Ac
      END DO

      IF (PRESENT(Val)) THEN
         DO a=1, nNo
            lhs%face(faIn)%val(:,a) = Val(:,a)
         END DO
      ELSE
         lhs%face(faIn)%val = 0D0
      END IF

      IF (lhs%commu%nTasks .GT. 1) THEN
         a = 0
         IF (lhs%face(faIn)%nNo .NE. 0) a = 1
         CALL MPI_ALLREDUCE(a, Ac, 1, mpint, MPI_SUM, lhs%commu%comm, i)
         IF (Ac .GT. 1) THEN
            lhs%face(faIn)%sharedFlag = .TRUE.
            IF (.NOT.ALLOCATED(v)) ALLOCATE(v(dof,lhs%nNo))
            v = 0D0
            DO a=1, nNo
               Ac = lhs%face(faIn)%glob(a)
               v(:,Ac) = lhs%face(faIn)%val(:,a)
            END DO
            CALL FSILS_COMMUV(lhs, dof, v)

            DO a=1, nNo
               Ac = lhs%face(faIn)%glob(a)
               lhs%face(faIn)%val(:,a) = v(:,Ac)
            END DO
         END IF
      END IF

      RETURN
      END SUBROUTINE FSILS_BC_CREATE

!====================================================================

      SUBROUTINE external_BC_CREATE (lhs, faIn, nNo, dof, BC_type,      &
     &         gNodes, Val)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: faIn, nNo, dof
      INTEGER, INTENT(IN) :: BC_type
      INTEGER, INTENT(IN) :: gNodes(nNo)
      REAL(KIND=8), INTENT(IN), OPTIONAL :: Val(dof,nNo)

      INTEGER a, Ac, i
      REAL(KIND=8), ALLOCATABLE :: v(:,:)

      IF (faIn .GT. lhs%nFaces) THEN
         PRINT *, "FSILS: faIn is exceeding lhs structure maximum",     &
     &      " number of face:", lhs%nFaces, " .LT. ", faIn
         STOP "FSILS: FATAL ERROR"
      END IF
      IF (faIn .LE. 0) THEN
         PRINT *, "FSILS: faIn should be greater than zero"
         STOP "FSILS: FATAL ERROR"
      END IF

      IF (lhs%face(faIn)%foC) THEN
         PRINT *, "FSILS: face is not free, you may use FSILS_BC_FREE", &
     &      " to free it"
         STOP "FSILS: FATAL ERROR"
      END IF

      lhs%face(faIn)%nNo  = nNo
      lhs%face(faIn)%dof  = dof
      lhs%face(faIn)%bGrp = BC_type

      ALLOCATE(lhs%face(faIn)%glob(nNo), lhs%face(faIn)%val(dof,nNo),   &
     &   lhs%face(faIn)%valM(dof,nNo))

      DO a=1, nNo
         Ac = lhs%map(gNodes(a))
         lhs%face(faIn)%glob(a) = Ac
      END DO

      IF (PRESENT(Val)) THEN
         DO a=1, nNo
            lhs%face(faIn)%val(:,a) = Val(:,a)
         END DO
      ELSE
         lhs%face(faIn)%val = 0D0
      END IF

      IF (lhs%commu%nTasks .GT. 1) THEN
         a = 0
         IF (lhs%face(faIn)%nNo .NE. 0) a = 1
         CALL MPI_ALLREDUCE(a, Ac, 1, mpint, MPI_SUM, lhs%commu%comm, i)
         IF (Ac .GT. 1) THEN
            lhs%face(faIn)%sharedFlag = .TRUE.
            IF (.NOT.ALLOCATED(v)) ALLOCATE(v(dof,lhs%nNo))
            v = 0D0
            DO a=1, nNo
               Ac = lhs%face(faIn)%glob(a)
               v(:,Ac) = lhs%face(faIn)%val(:,a)
            END DO
            CALL FSILS_COMMUV(lhs, dof, v)

            DO a=1, nNo
               Ac = lhs%face(faIn)%glob(a)
               lhs%face(faIn)%val(:,a) = v(:,Ac)
            END DO
         END IF
      END IF

      RETURN
      END SUBROUTINE external_BC_CREATE

!====================================================================

      SUBROUTINE FSILS_BC_FREE (lhs, faIn)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: faIn

      IF (.NOT.lhs%face(faIn)%foC) THEN
         PRINT *, 'FSILS: Cannot free a face that is not created yet'
         STOP "FSILS: FATAL ERROR"
      END IF
      lhs%face(faIn)%foC        = .FALSE.
      lhs%face(faIn)%nNo        = 0
      lhs%face(faIn)%bGrp       = BC_TYPE_Dir
      lhs%face(faIn)%res        = 0D0
      lhs%face(faIn)%sharedFlag = .FALSE.

      DEALLOCATE(lhs%face(faIn)%glob, lhs%face(faIn)%val,               &
     &   lhs%face(faIn)%valM)

      RETURN
      END SUBROUTINE FSILS_BC_FREE


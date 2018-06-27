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
!     This creates the communication structure based on the shared
!     nodes on different processors.
!--------------------------------------------------------------------

      SUBROUTINE FSILS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz, gNodes,   &
     &   rowPtr, colPtr, nFaces)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      INTEGER, INTENT(IN) :: gnNo, nNo, nnz
      INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1), colPtr(nnz)
      INTEGER, INTENT(IN) :: nFaces

      INTEGER ip, i, j, a, Ac, ai, s, e, nTasks, tF, maxnNo, ierr,      &
     &   comm, stat(mpsts)

      INTEGER, ALLOCATABLE :: aNodes(:,:), gtlPtr(:), ltg(:), part(:),  &
     &   sCount(:), disp(:)


      IF (lhs%foC) THEN
         PRINT *, "FSILS: LHS is not free You may use FSILS_LHS_FREE",  &
     &      " to free this structure"
         STOP "FSILS: FATAL ERROR"
      END IF

      lhs%foC    = .TRUE.
      lhs%gnNo   = gnNo
      lhs%nNo    = nNo
      lhs%nnz    = nnz
      lhs%commu  = commu
      lhs%nFaces = nFaces

      nTasks = commu%nTasks
      comm   = commu%comm
      tF     = commu%tF

      ALLOCATE (lhs%colPtr(nnz), lhs%rowPtr(2,nNo), lhs%diagPtr(nNo),   &
     &   lhs%map(nNo), lhs%face(nFaces))

!     IF it is called sequentially
      IF (nTasks .EQ. 1) THEN
         DO i=1, nnz
            lhs%colPtr(i) = colPtr(i)
         END DO
         DO Ac=1, nNo
            s = rowPtr(Ac)
            e = rowPtr(Ac+1) - 1
            DO i=s, e
               a = colPtr(i)
               IF (Ac .EQ. a) THEN
                  lhs%diagPtr(Ac) = i
                  EXIT
               END IF
            END DO

            lhs%rowPtr(1,Ac) = s
            lhs%rowPtr(2,Ac) = e

            lhs%map(Ac) = Ac
         END DO

         lhs%mynNo = nNo
         RETURN
      END IF

      CALL MPI_ALLREDUCE (nNo, maxnNo, 1, mpint, MPI_MAX, comm, ierr)

      ALLOCATE(aNodes(maxnNo,nTasks), part(maxnNo), sCount(nTasks),     &
     &   disp(nTasks), gtlPtr(gnNo), ltg(nNo))

      part = 0
      part(1:nNo) = gNodes

      DO i=1, nTasks
         disp(i)   = (i-1)*maxnNo
         sCount(i) = maxnNo
      END DO
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,    &
     &   mpint, comm, ierr)

      gtlPtr = 0
      DO a=1, nNo
         Ac = gNodes(a)
         gtlPtr(Ac) = a
      END DO

!     First including the nodes that belong to processors with higher ID
!     then including the nodes shared by lower processors IDs. shnNo is
!     counter for lower ID and mynNo is counter for higher ID
      lhs%mynNo = nNo
      lhs%shnNo = 0
      DO i=nTasks, 1, -1
!     Will include local nodes later
         IF (i .EQ. tF) CYCLE
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,i)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = gtlPtr(Ac)
            IF (ai .NE. 0) THEN
!     If this node has not been included already
               IF (aNodes(ai,tF) .NE. 0) THEN
                  IF (i .LT. tF) THEN
!     If the processor ID is lower, it is appended to the beginning
                     lhs%shnNo = lhs%shnNo + 1
                     ltg(lhs%shnNo) = Ac
                  ELSE
!     If the processor ID is higher, it is appended to the end
                     ltg(lhs%mynNo) = Ac
                     lhs%mynNo = lhs%mynNo - 1
                  END IF
                  aNodes(ai,tF) = 0
               END IF
            END IF
         END DO
      END DO
!     Now including the local nodes that are left behind
      j = lhs%shnNo + 1
      DO a=1, nNo
         Ac = aNodes(a,tF)
!     If this node has not been included already
         IF (Ac .NE. 0) THEN
            ltg(j) = Ac
            j = j + 1
         END IF
      END DO
      IF (j .NE. lhs%mynNo+1) THEN
         PRINT *, "FSILS: Unexpected behavior", j, lhs%mynNo
         CALL MPI_FINALIZE(ierr)
         STOP "FSILS: FATAL ERROR"
      END IF

!     Having the new ltg pointer, map is constructed
      gtlPtr = 0
      DO a=1, nNo
         Ac = ltg(a)
         gtlPtr(Ac) = a
      END DO
      DO a=1, nNo
         Ac = gNodes(a)
         lhs%map(a) = gtlPtr(Ac)
      END DO

!     Based on the new ordering of the nodes, rowPtr and colPtr are
!     constructed
      DO a=1, nNo
         Ac = lhs%map(a)
         lhs%rowPtr(1,Ac) = rowPtr(a)
         lhs%rowPtr(2,Ac) = rowPtr(a+1) - 1
      END DO
      DO i=1, nnz
         lhs%colPtr(i) = lhs%map(colPtr(i))
      END DO
!     diagPtr points to the diagonal entries of LHS
      DO Ac=1, nNo
         DO i=lhs%rowPtr(1,Ac), lhs%rowPtr(2,Ac)
            a = lhs%colPtr(i)
            IF (Ac .EQ. a) THEN
               lhs%diagPtr(Ac) = i
               EXIT
            END IF
         END DO
      END DO

!     Constructing the communication data structure based on the ltg
      part(1:nNo) = ltg
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,    &
     &   mpint, comm, ierr)

!     This variable keeps track of number of shared nodes
      disp = 0
      lhs%nReq = 0
      DO i=1, nTasks
         IF (i .EQ. tF) CYCLE
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,i)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = gtlPtr(Ac)
            IF (ai .NE. 0) THEN
               disp(i) = disp(i) + 1
            END IF
         END DO
         IF (disp(i) .NE. 0) lhs%nReq = lhs%nReq + 1
      END DO
      ALLOCATE(lhs%cS(lhs%nReq))
!     Now that we know which processor is communicating to which, we can
!     setup the handels and structures
      j = 0
      DO i=1, nTasks
         a = disp(i)
         IF (a .NE. 0) THEN
            j = j + 1
            lhs%cS(j)%iP = i
            lhs%cS(j)%n = a
            ALLOCATE(lhs%cS(j)%ptr(a))
         END IF
      END DO

!     Order of nodes in ptr is based on the node order in porcessor
!     with higher ID. ptr is calculated for tF+1:nTasks and will be
!     sent over.
      DO i=1, lhs%nReq
         iP = lhs%cS(i)%iP
         IF (iP .LT. tF) THEN
            CALL MPI_RECV(lhs%cS(i)%ptr, lhs%cS(i)%n, mpint, iP-1, 1,   &
     &         comm, stat, ierr)

            DO j=1, lhs%cS(i)%n
               lhs%cS(i)%ptr(j) = gtlPtr(lhs%cS(i)%ptr(j))
            END DO
         ELSE
!     This is a counter for the shared nodes
            j = 0
            DO a=1, maxnNo
!     Global node number in processor i at location a
               Ac = aNodes(a,iP)
!     Exit if this is the last node
               IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
               ai = gtlPtr(Ac)
               IF (ai .NE. 0) THEN
                  j = j + 1
!     Just for now global node ID is used. Later on this will be changed
!     to make sure nodes corresponds to each other on both processors
!     and then will be transformed to local node IDs
                  lhs%cS(i)%ptr(j) = Ac
               END IF
            END DO
            CALL MPI_SEND(lhs%cS(i)%ptr, lhs%cS(i)%n, mpint, iP-1, 1,   &
     &         comm, stat, ierr)

            DO j=1, lhs%cS(i)%n
               lhs%cS(i)%ptr(j) = gtlPtr(lhs%cS(i)%ptr(j))
            END DO
         END IF
      END DO

      DEALLOCATE(aNodes, part, gtlPtr, sCount, disp, ltg)

      RETURN
      END SUBROUTINE FSILS_LHS_CREATE

!====================================================================

      SUBROUTINE external_LHS_CREATE(lhs, commu, gnNo, nNo, gNodes,     &
     &   nFaces)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      INTEGER, INTENT(IN) :: gnNo, nNo
      INTEGER, INTENT(IN) :: gNodes(nNo)
      INTEGER, INTENT(IN) :: nFaces

      INTEGER ip, i, j, a, Ac, ai, nTasks, tF, maxnNo, ierr,            &
     &   comm, stat(mpsts)

      INTEGER, ALLOCATABLE :: aNodes(:,:), gtlPtr(:), ltg(:), part(:),  &
     &   sCount(:), disp(:)

      IF (lhs%foC) THEN
         PRINT *, "FSILS: LHS is not free You may use FSILS_LHS_FREE",  &
     &      " to free this structure"
         STOP "FSILS: FATAL ERROR"
      END IF

      lhs%foC    = .TRUE.
      lhs%gnNo   = gnNo
      lhs%nNo    = nNo
      lhs%nnz    = 0
      lhs%commu  = commu
      lhs%nFaces = nFaces

      nTasks = commu%nTasks
      comm   = commu%comm
      tF     = commu%tF

      ALLOCATE (lhs%map(nNo), lhs%face(nFaces))

      CALL MPI_ALLREDUCE (nNo, maxnNo, 1, mpint, MPI_MAX, comm, ierr)

      ALLOCATE(aNodes(maxnNo,nTasks), part(maxnNo), sCount(nTasks),     &
     &   disp(nTasks), gtlPtr(gnNo), ltg(nNo))

      part = 0
      part(1:nNo) = gNodes

      DO i=1, nTasks
         disp(i)   = (i-1)*maxnNo
         sCount(i) = maxnNo
      END DO
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,    &
     &   mpint, comm, ierr)

      gtlPtr = 0
      DO a=1, nNo
         Ac = gNodes(a)
         gtlPtr(Ac) = a
      END DO

!     First including the nodes that belong to processors with higher ID
!     then including the nodes shared by lower processors IDs. shnNo is
!     counter for lower ID and mynNo is counter for higher ID
      lhs%mynNo = nNo
      lhs%shnNo = 0
      DO i=nTasks, 1, -1
!     Will include local nodes later
         IF (i .EQ. tF) CYCLE
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,i)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = gtlPtr(Ac)
            IF (ai .NE. 0) THEN
!     If this node has not been included already
               IF (aNodes(ai,tF) .NE. 0) THEN
                  IF (i .LT. tF) THEN
!     If the processor ID is lower, it is appended to the beginning
                     lhs%shnNo = lhs%shnNo + 1
                     ltg(lhs%shnNo) = Ac
                  ELSE
!     If the processor ID is higher, it is appended to the end
                     ltg(lhs%mynNo) = Ac
                     lhs%mynNo = lhs%mynNo - 1
                  END IF
                  aNodes(ai,tF) = 0
               END IF
            END IF
         END DO
      END DO
!     Now including the local nodes that are left behind
      j = lhs%shnNo + 1
      DO a=1, nNo
         Ac = aNodes(a,tF)
!     If this node has not been included already
         IF (Ac .NE. 0) THEN
            ltg(j) = Ac
            j = j + 1
         END IF
      END DO
      IF (j .NE. lhs%mynNo+1) THEN
         PRINT *, "FSILS: Unexpected behavior", j, lhs%mynNo
         CALL MPI_FINALIZE(ierr)
         STOP "FSILS: FATAL ERROR"
      END IF

!     Having the new ltg pointer, map is constructed
      gtlPtr = 0
      DO a=1, nNo
         Ac = ltg(a)
         gtlPtr(Ac) = a
      END DO
      DO a=1, nNo
         Ac = gNodes(a)
         lhs%map(a) = gtlPtr(Ac)
      END DO

!     Constructing the communication data structure based on the ltg
      part(1:nNo) = ltg
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,    &
     &   mpint, comm, ierr)

!     This variable keeps track of number of shared nodes
      disp = 0
      lhs%nReq = 0
      DO i=1, nTasks
         IF (i .EQ. tF) CYCLE
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,i)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = gtlPtr(Ac)
            IF (ai .NE. 0) THEN
               disp(i) = disp(i) + 1
            END IF
         END DO
         IF (disp(i) .NE. 0) lhs%nReq = lhs%nReq + 1
      END DO
      ALLOCATE(lhs%cS(lhs%nReq))
!     Now that we know which processor is communicating to which, we can
!     setup the handels and structures
      j = 0
      DO i=1, nTasks
         a = disp(i)
         IF (a .NE. 0) THEN
            j = j + 1
            lhs%cS(j)%iP = i
            lhs%cS(j)%n = a
            ALLOCATE(lhs%cS(j)%ptr(a))
         END IF
      END DO

!     Order of nodes in ptr is based on the node order in porcessor
!     with higher ID. ptr is calculated for tF+1:nTasks and will be
!     sent over.
      DO i=1, lhs%nReq
         iP = lhs%cS(i)%iP
         IF (iP .LT. tF) THEN
            CALL MPI_RECV(lhs%cS(i)%ptr, lhs%cS(i)%n, mpint, iP-1, 1,   &
     &         comm, stat, ierr)

            DO j=1, lhs%cS(i)%n
               lhs%cS(i)%ptr(j) = gtlPtr(lhs%cS(i)%ptr(j))
            END DO
         ELSE
!     This is a counter for the shared nodes
            j = 0
            DO a=1, maxnNo
!     Global node number in processor i at location a
               Ac = aNodes(a,iP)
!     Exit if this is the last node
               IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
               ai = gtlPtr(Ac)
               IF (ai .NE. 0) THEN
                  j = j + 1
!     Just for now global node ID is used. Later on this will be changed
!     to make sure nodes corresponds to each other on both processors
!     and then will be transformed to local node IDs
                  lhs%cS(i)%ptr(j) = Ac
               END IF
            END DO
            CALL MPI_SEND(lhs%cS(i)%ptr, lhs%cS(i)%n, mpint, iP-1, 1,   &
     &         comm, stat, ierr)

            DO j=1, lhs%cS(i)%n
               lhs%cS(i)%ptr(j) = gtlPtr(lhs%cS(i)%ptr(j))
            END DO
         END IF
      END DO

      DEALLOCATE(aNodes, part, gtlPtr, sCount, disp, ltg)

      RETURN
      END SUBROUTINE external_LHS_CREATE
!====================================================================

      SUBROUTINE FSILS_LHS_FREE(lhs)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs

      INTEGER faIn, i

      IF (.NOT.lhs%foC) THEN
         PRINT *, 'FSILS: LHS is not created yet to be freed'
         STOP "FSILS: FATAL ERROR"
      END IF

      DO faIn = 1, lhs%nFaces
         IF (lhs%face(faIn)%foC) CALL FSILS_BC_FREE(lhs, faIn)
      END DO

      DO i=1, lhs%nReq
         IF (ALLOCATED(lhs%cS(i)%ptr)) DEALLOCATE(lhs%cS(i)%ptr)
      END DO

      lhs%foC    = .FALSE.
      lhs%gnNo   = 0
      lhs%nNo    = 0
      lhs%nnz    = 0
      lhs%nFaces = 0

      IF (ALLOCATED(lhs%colPtr)) DEALLOCATE(lhs%colPtr)
      IF (ALLOCATED(lhs%rowPtr)) DEALLOCATE(lhs%rowPtr)
      IF (ALLOCATED(lhs%diagPtr)) DEALLOCATE(lhs%diagPtr)
      IF (ALLOCATED(lhs%cS)) DEALLOCATE(lhs%cS)
      DEALLOCATE(lhs%map, lhs%face)

      RETURN
      END SUBROUTINE FSILS_LHS_FREE

!====================================================================

      SUBROUTINE FSILS_LHS_CREATE_C(pLHS, commu, gnNo, nNo, nnz, gNodes,&
     &   rowPtr, colPtr, nFaces)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), POINTER, INTENT(OUT) :: pLHS
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      INTEGER, INTENT(IN) :: gnNo, nNo, nnz
      INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1), colPtr(nnz)
      INTEGER, INTENT(IN) :: nFaces

      TYPE(FSILS_lhsType), TARGET, SAVE :: lhs

      CALL FSILS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz, gNodes, rowPtr, &
     &   colPtr, nFaces)

      pLHS => lhs

      RETURN
      END SUBROUTINE FSILS_LHS_CREATE_C

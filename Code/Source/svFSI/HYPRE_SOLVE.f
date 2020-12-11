!
! Copyright (c) Stanford University, The Regents of the University of
!           California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to DO so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
! Interface to Hypre linear solver package.
!
!--------------------------------------------------------------------
      MODULE HYPREMOD
      USE TYPEMOD

!     Data type for parallel communication
      TYPE hp_cSType
         ! The processor to communicate with
         INTEGER(KIND=IKIND) iP
         ! Number of rows to be commu
         INTEGER(KIND=IKIND) n
         ! Pointer to the row ID for commu
         INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
         ! No. of non-zero elements for commu
         INTEGER(KIND=IKIND) nnz
         ! Pointer to colPtr for commu
         INTEGER(KIND=IKIND), ALLOCATABLE :: cptr(:)
      END TYPE hp_cSType

!     Data type for Hypre Linear Solver related arrays
      TYPE hpType
         LOGICAL :: foC = .FALSE.
         ! Initialized or not
         LOGICAL :: init = .FALSE.
         ! No. of rows in this process
         INTEGER(KIND=IKIND) :: mynNo
         ! No. of none zero entries
         INTEGER(KIND=IKIND) :: nnz
         ! Used to determine if re-initialization is needed 
         INTEGER(KIND=IKIND) :: dof = 0
         ! Global lower ID
         INTEGER(KIND=IKIND) :: ilower
         ! Local to global mapping
         INTEGER(KIND=IKIND), ALLOCATABLE :: ltg(:)
         ! Row pointer, in local index
         INTEGER(KIND=IKIND), ALLOCATABLE :: rowPtr(:,:)
         ! Column pointer, in svFSI global index
         INTEGER(KIND=IKIND), ALLOCATABLE :: colPtr(:)
         ! Diagonal pointer
         INTEGER(KIND=IKIND), ALLOCATABLE :: diagPtr(:)
         ! RHS vector pointer
         REAL(KIND=RKIND), ALLOCATABLE :: R(:,:)
         ! Stiffness matrix pointer
         REAL(KIND=RKIND), ALLOCATABLE :: Val(:,:)
         ! Number of send requests    (USE)
         INTEGER(KIND=IKIND) :: snReq = 0
         ! Number of recv requests    (USE)
         INTEGER(KIND=IKIND) :: rnReq = 0
         ! TYPE(FSILS_commuType) commu
         TYPE(hp_cSType), ALLOCATABLE :: scS(:) ! send
         TYPE(hp_cSType), ALLOCATABLE :: rcS(:) ! recv
      END TYPE hpType

      ! Hypre parameter
      INTEGER, PARAMETER ::   HYPRE_PARCSR=5555
      INTEGER(KIND=IKIND8) :: parcsr_A
      INTEGER(KIND=IKIND8) :: AA, bb, xx
      INTEGER(KIND=IKIND8) :: par_b, par_x
      INTEGER(KIND=IKIND8) :: solver, precond
      REAL(KIND=RKIND) :: reltol

      ! Sparse matrix in terms of individual dof
      INTEGER(KIND=IKIND), ALLOCATABLE :: lcolPtr(:), lrowPtr(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: lVal(:), lR(:)
      ! Mapping between individual dof and node variables
      INTEGER(KIND=IKIND), ALLOCATABLE :: pvmap(:,:)

      ! Hypre Linear Solver data type
      TYPE(hpType) hp

      CONTAINS
         
         SUBROUTINE HYPRE_INITIALIZE(dof,lower,nNo,mpi_comm)
         USE COMMOD, ONLY : dbg

         INTEGER(KIND=IKIND), INTENT(IN) :: dof,lower,nNo,mpi_comm

         INTEGER(KIND=IKIND) :: ierr, local_size, ilower, iupper

         dbg = "Initialize Hypre....."

         local_size = dof*nNo
         ilower = (lower-1)*dof + 1 
         iupper = ilower + local_size - 1

         IF (.NOT. hp%init) CALL HYPRE_INIT(ierr)
         
         ! Special case for multiple equaitons
         IF (hp%init) THEN
            CALL HYPRE_IJMatrixDestroy(AA, ierr)
            CALL HYPRE_IJVectorDestroy(bb, ierr)
            CALL HYPRE_IJVectorDestroy(xx, ierr)
         END IF

         ! Create the matrix.
         ! Note that this is a square matrix, so we indicate the row partition
         ! size twice (since number of rows = number of cols)
         ! This is true even for parallel case.
         CALL HYPRE_IJMatrixCreate(mpi_comm, ilower,
     2      iupper, ilower, iupper, AA, ierr)      
         IF (ierr .NE. 0) 
     2      PRINT *, "HYPRE_IJMatrixCreate", ierr
         ! Set CSR format for A
         CALL HYPRE_IJMatrixSetObjectType(AA, HYPRE_PARCSR, ierr)
         IF (ierr .NE. 0)
     2      PRINT *, "HYPRE_IJMatrixSetObjectType", ierr

         ! Create the rhs and solution
         CALL HYPRE_IJVectorCreate(mpi_comm,
     2      ilower, iupper, bb, ierr)
         CALL HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR, ierr)

         CALL HYPRE_IJVectorCreate(mpi_comm,
     2     ilower, iupper, xx, ierr)
         CALL HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR, ierr)

         END SUBROUTINE HYPRE_INITIALIZE

         !-----------------------------------------------------------

         SUBROUTINE HYPRE_DESTROY()
         USE COMMOD, ONLY : dbg
         INTEGER(KIND=IKIND) ierr
         dbg = "Destroy Hypre ...."

         CALL HYPRE_IJMatrixDestroy(AA, ierr)
         CALL HYPRE_IJVectorDestroy(bb, ierr)
         CALL HYPRE_IJVectorDestroy(xx, ierr)
         CALL HYPRE_Finalize(ierr)

         END SUBROUTINE HYPRE_DESTROY

      END MODULE HYPREMOD  

!--------------------------------------------------------------------

      SUBROUTINE HYPRE_LHS_CREATE(commu, gnNo, nNo, gNodes)
      USE COMMOD
      USE HYPREMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, gnNo
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      INTEGER(KIND=IKIND), INTENT(IN) :: gNodes(nNo)

      INTEGER(KIND=IKIND) :: i, j, k, iP, ptr, is, ie, Ac, a, ai
      INTEGER(KIND=IKIND) :: nTasks, myID, ierr, comm, 
     2                       stat(MPI_STATUS_SIZE)
      INTEGER(KIND=IKIND) :: maxrow, maxnnz, maxnNo, itemp
      INTEGER(KIND=IKIND), ALLOCATABLE :: l2g(:), g2l(:),
     2                                    mynNo(:), colindx(:), part(:),
     3                                    disp(:), aNodes(:,:), 
     4                                    sCount(:), sReq(:), rReq(:),
     5                                    snrow(:,:), rnrow(:,:),
     6                                    scolPtr(:,:), rcolPtr(:,:),
     7                                    nnrow(:), mask(:,:)

      nTasks = commu%nTasks
      comm   = commu%comm
      myID   = commu%tF

      hp%mynNo = lhs%mynNo
      hp%foC   = .TRUE.

!     If it is called sequentially
      IF (nTasks .EQ. 1) THEN
         hp%nnz = lhs%nnz
         ALLOCATE(hp%rowPtr(2,hp%mynNo),hp%colPtr(hp%nnz))
         ALLOCATE(hp%diagPtr(hp%mynNo),hp%ltg(hp%mynNo))
         hp%ilower  = 1
         hp%ltg     = lhs%map
         hp%rowPtr  = lhs%rowPtr
         hp%colPtr  = lhs%colPtr
         hp%diagPtr = lhs%diagPtr
         RETURN
      END IF

      ALLOCATE(mynNo(nTasks))
      mynNo = 0
      CALL MPI_ALLGATHER(lhs%mynNo, 1, mpint, mynNo, 1, mpint, 
     2   comm, ierr)

!     Lower bound of index
      hp%ilower = 1
      DO i = 1, myID-1
         hp%ilower = hp%ilower + mynNo(i)
      END DO

!     Local to global mapping
      ALLOCATE(l2g(nNo),hp%ltg(hp%mynNo))
      DO i = 1, nNo
         l2g(lhs%map(i)) = gNodes(i)
      END DO
      hp%ltg = l2g(1:hp%mynNo)

!     Map colPtr to global index
      ALLOCATE(colindx(lhs%nnz))
      DO i = 1, lhs%nnz
         colindx(i) = l2g(lhs%colPtr(i))
      END DO

!     ------------------------------------------------------
!     Track recv request
!     ------------------------------------------------------
      CALL MPI_ALLREDUCE(nNo, maxnNo, 1, mpint, MPI_MAX, comm, ierr)
      ALLOCATE(aNodes(maxnNo,nTasks),part(maxnNo))
      ALLOCATE(sCount(nTasks),disp(nTasks))
      part = 0
      part(1:nNo) = l2g
      DO i=1, nTasks
         disp(i)   = (i-1)*maxnNo
         sCount(i) = maxnNo
      END DO
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,
     2   mpint, comm, ierr)

!     Global to local mapping 
!     only points owned solely by the current processor
      ALLOCATE(g2l(gnNo))
      g2l = 0
      DO i = 1, hp%mynNo
         g2l(l2g(i)) = i
      END DO

      disp     = 0
      hp%rnReq = 0
      DO i=1, nTasks
         IF (i .EQ. myID) CYCLE
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,i)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = g2l(Ac)
            IF (ai .NE. 0) THEN
               disp(i) = disp(i) + 1
            END IF
         END DO
         IF (disp(i) .NE. 0) hp%rnReq = hp%rnReq + 1
      END DO
      ALLOCATE(hp%rcS(hp%rnReq))
!     Now that we know which processor is sending to currnt processor,
!     we can setup the handles and structures
      hp%rCs(:)%n = 0
      j = 0
      DO i=1, nTasks
         a = disp(i)
         IF (a .NE. 0) THEN
            j = j + 1
            hp%rcS(j)%iP = i
            hp%rcS(j)%n  = a
            ALLOCATE(hp%rcS(j)%ptr(a))
         END IF
      END DO
!     Constract ID of received node
      DO i=1, hp%rnReq
         iP = hp%rcS(i)%iP
!     This is a counter for the shared nodes
         j = 0
         DO a=1, maxnNo
!     Global node number in processor i at location a
            Ac = aNodes(a,iP)
!     Exit if this is the last node
            IF (Ac .EQ. 0) EXIT
!     Corresponding local node in current processor
            ai = g2l(Ac)
            IF (ai .NE. 0) THEN
               j = j + 1
!     Just for now global node ID is used. Later on this will be changed
!     to make sure nodes corresponds to each other on both processors
!     and then will be transformed to local node IDs
               hp%rcS(i)%ptr(j) = Ac
            END IF
         END DO
      END DO

!     ------------------------------------------------------
!     Track send request
!     ------------------------------------------------------
      CALL MPI_ALLREDUCE(hp%rnReq,maxnNo,1,mpint,MPI_MAX,comm,ierr)
      DEALLOCATE(aNodes,part)
      ALLOCATE(aNodes(maxnNo,nTasks),part(maxnNo))
      part = 0
      part(1:hp%rnReq) = hp%rcS(:)%iP
      DO i=1, nTasks
         disp(i)   = (i-1)*maxnNo
         sCount(i) = maxnNo
      END DO
      CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp,
     2   mpint, comm, ierr)

!     No. of send request from current processor
      disp     = 0
      hp%snReq = 0
      j        = 0
      DO i=1, nTasks
         IF (i .EQ. myID) CYCLE
         DO a=1, maxnNo
            Ac = aNodes(a,i)
            IF (Ac .EQ. myID) THEN
               j       = j + 1
               disp(j) = i
            END IF
         END DO
      END DO
      hp%snReq = j
      ALLOCATE(hp%scS(hp%snReq))
      hp%scS(:)%iP = disp(1:hp%snReq)

!     Get No. of outgoing nodes from receiving processors
      hp%sCs(:)%n = 0
      ALLOCATE(rReq(hp%rnReq),sReq(hp%snReq))
      DO i=1, hp%rnReq
         CALL MPI_ISEND(hp%rcS(i)%n,1,mpint,hp%rcS(i)%iP-1,1,
     2                  comm,rReq(i),ierr)
      END DO
      DO i=1, hp%snReq
         CALL MPI_IRECV(hp%scS(i)%n,1,mpint,hp%scS(i)%iP-1,1,
     2                  comm,sReq(i),ierr)
      END DO
      DO i=1, hp%snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO
!     Get outgoing node IDs from receiving processors
      DO i=1, hp%rnReq
         CALL MPI_ISEND(hp%rcS(i)%ptr,hp%rcS(i)%n,mpint,hp%rcS(i)%iP-1,
     2                  1,comm,rReq(i),ierr)
      END DO
      DO i=1, hp%snReq
         ALLOCATE(hp%scS(i)%ptr(hp%scS(i)%n))
         CALL MPI_IRECV(hp%scS(i)%ptr,hp%scS(i)%n,mpint,hp%scS(i)%iP-1,
     2                  1,comm,sReq(i),ierr)
      END DO
      DO i=1, hp%snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO
!     Convert send/recv global node IDs to local node IDs
      g2l = 0
      DO i = 1, nNo
         g2l(l2g(i)) = i
      END DO
      DO i=1, hp%rnReq
         DO j=1,hp%rcS(i)%n
            hp%rcS(i)%ptr(j) = g2l(hp%rcS(i)%ptr(j))
         END DO
      END DO
      DO i=1, hp%snReq
         DO j=1,hp%scS(i)%n
            hp%scS(i)%ptr(j) = g2l(hp%scS(i)%ptr(j))
         END DO
      END DO

!     Determine size of colPtr needs to be communicated.
      itemp = MAXVAL(hp%scS(:)%n)
      CALL MPI_ALLREDUCE(itemp, maxrow, 1, mpint, MPI_MAX, comm, ierr)
      ALLOCATE(snrow(maxrow,hp%snReq))
      snrow = 0
      hp%scS(:)%nnz = 0
      DO i = 1, hp%snReq
         DO j = 1, hp%scS(i)%n
            ptr           = hp%scS(i)%ptr(j)
            snrow(j,i)    = lhs%rowPtr(2,ptr) - lhs%rowPtr(1,ptr) + 1
            hp%scS(i)%nnz = hp%scS(i)%nnz + snrow(j,i)
         END DO
      END DO

!     Build temporary colPtr for communication
      itemp = MAXVAL(hp%scS(:)%nnz)
      CALL MPI_ALLREDUCE(itemp, maxnnz, 1, mpint, MPI_MAX, comm, ierr)
      ALLOCATE(scolPtr(maxnnz,hp%snReq))
      scolPtr = 0
      DO i = 1, hp%snReq
         ALLOCATE(hp%scS(i)%cptr(hp%scS(i)%nnz))
         itemp = 0
         is = 0
         ie = 0
         ! temporary colPtr
         DO j = 1, hp%scS(i)%n
            is  = ie + 1
            ie  = is + snrow(j,i) - 1
            ptr = hp%scS(i)%ptr(j)
            scolPtr(is:ie,i) = 
     2          colindx(lhs%rowPtr(1,ptr):lhs%rowPtr(2,ptr))
            DO k = lhs%rowPtr(1,ptr),lhs%rowPtr(2,ptr)
               itemp = itemp + 1
               hp%scS(i)%cptr(itemp) = k
            END DO
         END DO
      END DO

!     Communicate to update nnz information on receiving process
      ALLOCATE(rnrow(maxrow,hp%rnReq))
      hp%rcS(:)%nnz  = 0
      rnrow = 0
      DO i=1, hp%snReq
         CALL MPI_ISEND(snrow(:,i),maxrow,mpint,hp%scS(i)%iP-1,
     2                  1,comm,sReq(i),ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_IRECV(rnrow(:,i),maxrow,mpint,hp%rcS(i)%iP-1,
     2                  1,comm,rReq(i),ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO
      DO i=1, hp%snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO
      DO i=1, hp%rnReq
         DO j = 1, hp%rcS(i)%n
            hp%rcS(i)%nnz = hp%rcS(i)%nnz + rnrow(j,i)
         END DO
      END DO

!     Communicate colPtr
      ALLOCATE(rcolPtr(maxnnz,hp%rnReq))
      rcolPtr = 0
      DO i=1, hp%snReq
         CALL MPI_ISEND(scolPtr(:,i),maxnnz,mpint,hp%scS(i)%iP-1,
     2                  1,comm,sReq(i),ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_IRECV(rcolPtr(:,i),maxnnz,mpint,hp%rcS(i)%iP-1,
     2                  1,comm,rReq(i),ierr)
      END DO
      DO i=1, hp%rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO
      DO i=1, hp%snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO

!     Search across all hp%rcS(:)%ptr(:) to identify unique rows   
      DEALLOCATE(part)
      ALLOCATE(part(hp%mynNo))
      part = 0
      k    = 0
      DO i=1, hp%rnReq
         DO j=1,hp%rcS(i)%n
            IF (part(hp%rcS(i)%ptr(j)) .EQ. 0) THEN
               k = k + 1
               part(hp%rcS(i)%ptr(j)) = k
            END IF
         END DO
      END DO
      maxnNo = k
      DEALLOCATE(disp)
      ALLOCATE(disp(maxnNo))
      DO i=1, hp%mynNo
         IF (part(i) .GT. 0) THEN
            disp(part(i)) = i
         END IF
      END DO

!     Calculate total No. of non-zero elements in each row
      ALLOCATE(nnrow(hp%mynNo), mask(gnNo,maxnNo))
      mask = 0
      ! Old size
      DO i = 1, hp%mynNo
         nnrow(i) = lhs%rowPtr(2,i) - lhs%rowPtr(1,i) + 1
      END DO
      ! New size
      DO k = 1, maxnNo
         itemp = 0
         ai    = disp(k)
         ! existing colPtr
         DO i = lhs%rowPtr(1,ai), lhs%rowPtr(2,ai)
            itemp = itemp + 1
            mask(colindx(i),k) =  itemp
         END DO
         ! received colPtr
         DO i = 1, hp%rnReq
            is = 0
            ie = 0
            DO j = 1, hp%rcS(i)%n
               is  = ie + 1
               ie  = is + rnrow(j,i) - 1
               IF (hp%rcS(i)%ptr(j) .EQ. ai) THEN
                  DO a = is, ie
                     IF (mask(rcolPtr(a,i),k) .EQ. 0) THEN
                        itemp = itemp + 1
                        mask(rcolPtr(a,i),k) = itemp
                     END IF
                  END DO
               END IF
            END DO
         END DO
         nnrow(ai) = itemp
      END DO

!     Build rowPtr for Hypre
      ALLOCATE(hp%rowPtr(2,hp%mynNo))
      hp%rowPtr(1,1) = 1
      hp%rowPtr(2,1) = nnrow(1)
      DO i = 2, hp%mynNo
         hp%rowPtr(1,i) = hp%rowPtr(2,i-1) + 1 
         hp%rowPtr(2,i) = hp%rowPtr(2,i-1) + nnrow(i)        
      END DO

!     Build ptr to convert received node ID for local node ID
      DO i = 1, hp%rnReq
         ALLOCATE(hp%rcS(i)%cptr(hp%rcS(i)%nnz))
         hp%rcS(i)%cptr = 0
         itemp = 0
         is = 0
         ie = 0
         DO j = 1, hp%rcS(i)%n
            is  = ie + 1
            ie  = is + rnrow(j,i) - 1
            ptr = hp%rcS(i)%ptr(j)
            ai  = part(ptr)
            Ac  = hp%rowPtr(1,ptr)
            DO k = is, ie
               hp%rcS(i)%cptr(k) = mask(rcolPtr(k,i),ai) +
     2                             Ac - 1
            END DO
         END DO
      END DO

!     Build colPtr for Hypre using global index
      hp%nnz = hp%rowPtr(2,hp%mynNo)
      ALLOCATE(hp%colPtr(hp%nnz))
      ! existing colPtr
      DO i = 1, hp%mynNo
         is = lhs%rowPtr(1,i)
         ie = lhs%rowPtr(2,i)
         j  = hp%rowPtr(1,i)
         hp%colPtr(j:(j+ie-is)) = colindx(is:ie)
      END DO
      ! received colPtr
      DO i = 1, hp%rnReq
         DO j = 1, hp%rcS(i)%nnz
            ai = hp%rcS(i)%cptr(j)
            hp%colPtr(ai) = rcolPtr(j,i)
         END DO
      END DO

!     diagPtr points to the diagonal entries of LHS
      ALLOCATE(hp%diagPtr(hp%mynNo))
      DO i=1, hp%mynNo
         DO j=hp%rowPtr(1,i), hp%rowPtr(2,i)
            ai = hp%colPtr(j)
            IF (hp%ltg(i) .EQ. ai) THEN
               hp%diagPtr(i) = j
               EXIT
            END IF
         END DO
      END DO

!     Deallocate
      DEALLOCATE(mynNo, l2g, colindx, aNodes, part, sCount)
      DEALLOCATE(disp, g2l, rReq, sReq, snrow, scolPtr)
      DEALLOCATE(rnrow, rcolPtr, nnrow, mask)

      RETURN
      END SUBROUTINE HYPRE_LHS_CREATE

!--------------------------------------------------------------------

      SUBROUTINE HYPRE_LHS_FREE()
      USE HYPREMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND) :: i

      IF (.NOT.hp%foC) THEN
         PRINT *, 'Hypre solve is not initialized.'
         STOP "HYPRE: FATAL ERROR"
      END IF

      DEALLOCATE(lcolPtr, lrowPtr, lVal, lR, pvmap)
      DEALLOCATE(hp%ltg, hp%rowPtr, hp%colPtr, hp%diagPtr)
      DEALLOCATE(hp%Val, hp%R)
      DO i = 1, hp%snReq
         IF(ALLOCATED(hp%scS(i)%ptr)) DEALLOCATE(hp%scS(i)%ptr)
         IF(ALLOCATED(hp%scS(i)%cptr))DEALLOCATE(hp%scS(i)%cptr)
      END DO
      DO i = 1, hp%rnReq
         IF(ALLOCATED(hp%rcS(i)%ptr)) DEALLOCATE(hp%rcS(i)%ptr)
         IF(ALLOCATED(hp%rcS(i)%cptr))DEALLOCATE(hp%rcS(i)%cptr)
      END DO
      IF(ALLOCATED(hp%scS)) DEALLOCATE(hp%scS)
      IF(ALLOCATED(hp%rcS)) DEALLOCATE(hp%rcS)

      hp%init   = .FALSE.
      hp%foC    = .FALSE.
      hp%mynNo  = 0
      hp%nnz    = 0
      hp%ilower = 0

      RETURN
      END SUBROUTINE HYPRE_LHS_FREE

!--------------------------------------------------------------------

      SUBROUTINE HYPRE_COMMU(Val,R)
      USE CMMOD
      USE HYPREMOD
      USE COMMOD, ONLY : lhs, dof
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(IN) :: Val(dof*dof,lhs%nnz), 
     2                                R(dof,lhs%nNo)

      INTEGER(KIND=IKIND) snReq, rnReq, i, j, k, ierr, is, ie,
     2                    stat(MPI_STATUS_SIZE)

      REAL(KIND=RKIND), ALLOCATABLE :: sB(:,:,:), rB(:,:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: sReq(:), rReq(:)

      IF (ALLOCATED(hp%Val)) DEALLOCATE(hp%Val)
      ALLOCATE(hp%Val(dof*dof,hp%nnz))
      IF (ALLOCATED(hp%R))   DEALLOCATE(hp%R)
      ALLOCATE(hp%R(dof,hp%mynNo))

      IF (lhs%commu%nTasks .EQ. 1) THEN
         hp%Val = Val
         hp%R   = R
         RETURN
      END IF

      DO i = 1, lhs%nNo
         j = lhs%map(i)
         IF (j .LE. hp%mynNo) THEN
            hp%R(:,j) = R(:,i)
         END IF
      END DO

      snReq = hp%snReq
      rnReq = hp%rnReq
      i = MAXVAL(hp%scS(:)%nnz)
      j = MAXVAL(hp%rcS(:)%nnz)
      ALLOCATE(sB(dof*dof,i,snReq), rB(dof*dof,j,rnReq))
      ALLOCATE(rReq(rnReq), sReq(snReq))

      DO i=1, snReq
         DO j=1, hp%scS(i)%nnz
            k = hp%scS(i)%cptr(j)
            sB(:,j,i) = Val(:,k)
         END DO
      END DO

      DO i=1, rnReq
         CALL MPI_IRECV(rB(:,:,i), hp%rcS(i)%nnz*dof*dof, mpreal,
     2      hp%rcS(i)%iP-1, 1, lhs%commu%comm, rReq(i), ierr)
      END DO
      DO i=1, snReq
         CALL MPI_ISEND(sB(:,:,i), hp%scS(i)%nnz*dof*dof, mpreal,
     2      hp%scS(i)%iP-1, 1, lhs%commu%comm, sReq(i), ierr)
      END DO

      DO i=1, rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO

      hp%Val = 0._RKIND
      ! existing Val
      DO i = 1, hp%mynNo
         is = lhs%rowPtr(1,i)
         ie = lhs%rowPtr(2,i)
         j  = hp%rowPtr(1,i)
         DO k = 0, ie-is
            hp%Val(:,j+k) = Val(:,is+k)
         END DO
      END DO
      ! received Val
      DO i=1, rnReq
         DO j = 1, hp%rcS(i)%nnz
            k = hp%rcS(i)%cptr(j)
            hp%Val(:,k) = hp%Val(:,k) + rB(:,j,i)
         END DO
      END DO

      DO i=1, snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO

      RETURN
      END SUBROUTINE HYPRE_COMMU

!--------------------------------------------------------------------
   
      SUBROUTINE HYPRE_SOLVE(Val, Ri, dof, nnz, nNo, lower, 
     2   rowPtr, colPtr, mpi_comm, myID, final_norm, initial_norm,
     3   num_iterations, ls_type, prec_type, success, maxiter, sD)
      USE HYPREMOD
      INCLUDE "CONSTS.f"
      INTEGER(KIND=IKIND), INTENT(IN) :: dof, nnz, nNo, lower
      INTEGER(KIND=IKIND), INTENT(IN) :: mpi_comm, myID
      INTEGER(KIND=IKIND), INTENT(IN) :: rowPtr(2,dof*nNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: colPtr(dof*dof*nnz)
      REAL(KIND=RKIND), INTENT(IN) :: Val(dof*dof*nnz)
      REAL(KIND=RKIND), INTENT(INOUT) :: Ri(dof*nNo)
      REAL(KIND=RKIND), INTENT(OUT) :: final_norm
      REAL(KIND=RKIND), INTENT(IN) :: initial_norm
      INTEGER(KIND=IKIND), INTENT(OUT) :: num_iterations
      INTEGER(KIND=IKIND), INTENT(IN) :: ls_type, prec_type
      INTEGER(KIND=IKIND), INTENT(IN) :: maxiter, sD
      LOGICAL, INTENT(OUT) :: success

      ! Local variable
      INTEGER(KIND=IKIND) :: local_size, nrow
      INTEGER(KIND=IKIND) :: i, ierr
      INTEGER(KIND=IKIND) :: ilower, iupper
      INTEGER(KIND=IKIND) :: rs, re
      
      INTEGER(KIND=IKIND) :: solver_id, precond_id
      INTEGER(KIND=IKIND) :: rows(dof*nNo)
      REAL(KIND=RKIND) :: x_values(dof*nNo)

      ! Hypre has specific ID for preconditioners
      solver_id = ls_type
      precond_id = prec_type - PREC_HYPRE_NONE

      local_size = dof*nNo
      ilower = (lower-1)*dof + 1 
      iupper = ilower + local_size - 1

      ! Initialize A
      CALL HYPRE_IJMatrixInitialize(AA, ierr) 
      IF (ierr .NE. 0)
     2      PRINT *, "HYPRE_IJMatrixInitialize", ierr
      ! Go through local rows and set the matrix entries.
      DO i=1, local_size
         rs = rowPtr(1,i)
         re = rowPtr(2,i)
         nrow = re - rs + 1
         ! Set the values for row i
         CALL HYPRE_IJMatrixSetValues( AA, 1, nrow,
     2      i+ilower-1, colPtr(rs:re), Val(rs:re), ierr)
         IF (ierr .NE. 0) THEN
           PRINT *, myID, "HYPRE_IJMatrixSetValues", ierr
           STOP
         END IF
      END DO 
      CALL HYPRE_IJMatrixAssemble(AA, ierr)
      CALL HYPRE_IJMatrixGetObject(AA, parcsr_A, ierr)

      ! Initialize the rhs and solution
      CALL HYPRE_IJVectorInitialize(bb, ierr)
      CALL HYPRE_IJVectorInitialize(xx, ierr)

      ! Set the rhs values
      DO i = 1, local_size
         x_values(i) = 0._RKIND
         rows(i)     = ilower + i - 1
      END DO
      CALL HYPRE_IJVectorSetValues(
     2     bb, local_size, rows, lR, ierr)
      CALL HYPRE_IJVectorSetValues(
     2     xx, local_size, rows, x_values, ierr)
      CALL HYPRE_IJVectorAssemble(bb, ierr)
      CALL HYPRE_IJVectorAssemble(xx, ierr)
      CALL HYPRE_IJVectorGetObject(bb, par_b, ierr)
      CALL HYPRE_IJVectorGetObject(xx, par_x, ierr)

      ! Debug CZHU
      ! Check the matrix system
      ! call HYPRE_IJVectorPrint(bb, "b.out", ierr)
      ! call HYPRE_IJMatrixPrint(AA, "A.out", ierr)

      ! Set iterative solvers
      SELECT CASE(solver_id)
      CASE (lSolver_GMRES)
         ! IF ( myID .EQ. 0 ) PRINT *, "GMRES solver"
         CALL HYPRE_ParCSRGMRESCreate(mpi_comm, solver, ierr)
         CALL HYPRE_ParCSRGMRESSetKDim(solver, sD, ierr)
         CALL HYPRE_ParCSRGMRESSetMaxIter(solver, maxiter, ierr)
         CALL HYPRE_ParCSRGMRESSetAbsoluteTol(solver, 0.D0, ierr)
         CALL HYPRE_ParCSRGMRESSetTol(solver, reltol, ierr)
         CALL HYPRE_ParCSRGMRESSetPrintLevel(solver, 0, ierr)
         CALL HYPRE_ParCSRGMRESSetLogging(solver, 0, ierr)
      CASE (lSolver_BICGS)
         ! IF ( myID .EQ. 0 ) PRINT *, "BiCGSTAB solver"
         CALL HYPRE_ParCSRBiCGSTABCreate(mpi_comm, solver, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetMaxIter(solver, maxiter, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetATol(solver, 0.D0, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetTol(solver, reltol, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetPrintLev(solver, 0, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetLogging(solver, 0, ierr)
      CASE DEFAULT
         PRINT *, "Hypre module only supports GMRES and BiCGSTAB."
      END SELECT

      ! Set preconditioners
      SELECT CASE(precond_id)
      CASE(0)
         PRINT *, "No preconditioner"
         precond = 0

      CASE(1)
         PRINT *, "Diagonal scaling as preconditioner"
         precond = 0

      CASE(2) 
         PRINT *, "AMG as preconditioner"
         ! Now set up the AMG preconditioner and specify any parameters
         CALL HYPRE_BoomerAMGCreate(precond, ierr)
         IF (ierr .NE. 0) WRITE(*,*) 'BoomerAMGCreate error'

         ! Set some parameters (See Reference Manual for more parameters)
         ! PRINT less solver info since a preconditioner
         CALL HYPRE_BoomerAMGSetPrintLevel(precond, 1, ierr)
         ! Falgout coarsening
         CALL HYPRE_BoomerAMGSetCoarsenType(precond, 0, ierr)
         ! old defaults
         CALL HYPRE_BoomerAMGSetOldDefault(precond, ierr)
         ! SYMMETRIC G-S/Jacobi hybrid relaxation
         CALL HYPRE_BoomerAMGSetRelaxType(precond, 1, ierr)
         ! Sweeeps on each level
         CALL HYPRE_BoomerAMGSetNumSweeps(precond, 1, ierr)
         ! conv. tolerance
         CALL HYPRE_BoomerAMGSetTol(precond, 1d-3, ierr)
         ! DO only one iteration!
         CALL HYPRE_BoomerAMGSetMaxIter(precond, 10, ierr)
    
      CASE(3)
         PRINT *, "PILUT as preconditioner"
         CALL HYPRE_ParCSRPilutCreate(mpi_comm,
     2                                precond, ierr) 
         IF (ierr .NE. 0) WRITE(*,*) 'ParCSRPilutCreate error'

      !   IF (drop_tol .ge. 0.)
      !   1   CALL HYPRE_ParCSRPilutSetDropToleran(precond,
      !   2                                        drop_tol, ierr)

      CASE(6)
         ! PRINT *, "Hypre ILU as preconditioner"

         ! CALL HYPRE_ILUCreate(precond, ierr)
         ! CALL HYPRE_ILUSetType(precond, 0, ierr)
         ! CALL HYPRE_ILUSetLevelOfFill(precond, 1, ierr)
         ! CALL HYPRE_ILUSetLocalReordering(precond, 0, ierr)
         ! ! CALL HYPRE_ILUSetSchurMaxIter()
         ! CALL HYPRE_ILUSetMaxIter(precond, 1000, ierr)

         ! CALL HYPRE_ILUSetPrintLevel(precond, 3, ierr)
         
      CASE(5)
         ! IF ( myID .EQ. 0 ) PRINT *, "Euclid as preconditioner"
         CALL HYPRE_EuclidCreate(mpi_comm,
     2                           precond, ierr) 
         ! Use block Jacobi ILU preconditioning instead of PILU
         CALL HYPRE_EuclidSetBJ(precond, 0, ierr)
         ! Set level k for ILU(k) factorization, default: 1
         CALL HYPRE_EuclidSetLevel(precond, 0, ierr)

         ! CALL HYPRE_EuclidSetRowScale(precond, 1, ierr)
         ! ! Defines a drop tolerance for ILU(k). Default: 0 
         ! CALL HYPRE_EuclidSetSparseA(precond, 1.D-4, ierr)

         ! Summary of runtime settings and timing information is printed to stdout.
         ! CALL HYPRE_EuclidSetStats(precond, 0, ierr)

      CASE DEFAULT
         PRINT *, " Preconditioner not supported in Hypre module."
      END SELECT

      ! Now setup and solve!
      SELECT CASE(solver_id)
      CASE(lSolver_GMRES)
         CALL HYPRE_ParCSRGMRESSetPrecond(solver, precond_id,
     2                                    precond, ierr)
         
         CALL HYPRE_ParCSRGMRESSetup(solver, parcsr_A, par_b,
     2                               par_x, ierr)
         CALL HYPRE_ParCSRGMRESSolve(solver, parcsr_A, par_b,
     2                               par_x, ierr)
         IF (ierr .EQ. 0) success = .TRUE.
         ! Run info - needed logging turned on
         CALL HYPRE_ParCSRGMRESGetNumIteratio(solver, num_iterations,
     2                                        ierr)
         CALL HYPRE_ParCSRGMRESGetFinalRelati(solver, final_norm,
     2                                        ierr)     
         CALL HYPRE_ParCSRGMRESDestroy(solver, ierr)
      
      CASE(lSolver_BICGS)
         CALL HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id,
     2                                       precond, ierr)
         CALL HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b,
     2                                  par_x, ierr)
         CALL HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b,
     2                                  par_x, ierr)
         IF (ierr .EQ. 0) success = .TRUE.
         ! Run info - needed logging turned on
         CALL HYPRE_ParCSRBiCGSTABGetNumIter(solver, num_iterations,
     2                                           ierr)
         CALL HYPRE_ParCSRBiCGSTABGetFinalRel(solver, final_norm,
     2                                           ierr)   
         CALL HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
      END SELECT

      ! Absolute final residual norm
      final_norm = final_norm * initial_norm
      ! Get solution
      CALL HYPRE_IJVectorGetValues(
     2     xx, local_size, rows, x_values, ierr)      

      ! Destroy precond
      SELECT CASE(precond_id)
      CASE(2)
         CALL HYPRE_BoomerAMGDestroy(precond, ierr)
      CASE(3) 
         CALL HYPRE_ParCSRPilutDestroy(precond, ierr)
      CASE(6)
         ! CALL HYPRE_ILUDestroy(precond, ierr)
      CASE(5)
         CALL HYPRE_EuclidDestroy(precond, ierr)
      END SELECT

      Ri = x_values

      ! Debug CZHU
      ! Check the matrix system
      ! call HYPRE_IJVectorPrint(xx, "x.out", ierr)

      RETURN
      END SUBROUTINE HYPRE_SOLVE

!--------------------------------------------------------------------

      SUBROUTINE HYPRE_PRE_SOLVE(dof, commu, ls)
      USE HYPREMOD
      INCLUDE "FSILS_STD.h"

      INTEGER(KIND=IKIND), INTENT(IN) :: dof
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls

      REAL(KIND=RKIND) FSILS_CPUT, FSILS_NORMV

      INTEGER(KIND=IKIND) :: i, j, k, l, i0, j0, rs, re, a

      ls%iNorm = FSILS_NORMV(dof, hp%mynNo, commu, hp%R)
      ls%callD = FSILS_CPUT()
      ls%suc   = .FALSE.

      ! tolerance
      reltol = MAX(ls%abstol/ls%iNorm, ls%reltol)

      ! Convert Val to dof sparse matrix
      ! Pre-vel separation is done through lcolPtr and lrowPtr
      a = 0
      lVal = 0._RKIND
      DO i=1, hp%mynNo
         i0   = (i-1)*dof
         rs   = hp%rowPtr(1,i)
         re   = hp%rowPtr(2,i)
         DO j=1, dof
            j0 = (j-1)*dof
            DO k=rs, re
               DO l=1, dof
                  a = a + 1
                  lVal(a) = hp%Val(j0+l,k)
               END DO
            END DO
         END DO 
      END DO
      lR = 0._RKIND
      DO i = 1, hp%mynNo
         DO j = 1, dof
            lR(pvmap(j,i)) = hp%R(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE HYPRE_PRE_SOLVE

!--------------------------------------------------------------------

      SUBROUTINE HYPRE_POST_SOLVE(Wc, Ri, ls)
      USE HYPREMOD
      USE COMMOD, ONLY : lhs, dof, tnNo, gtnNo
      INCLUDE "FSILS_STD.h"

      ! Debug CZHU, for FSI, is lhs%nNo == tnNo??
      REAL(KIND=RKIND), INTENT(IN)  :: Wc(dof,gtnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Ri(dof,tnNo)
      TYPE(FSILS_subLsType), INTENT(INOUT) :: ls

      INTEGER(KIND=IKIND) :: i, j, k, i0, ierr
      INTEGER(KIND=IKIND) :: nNo, comm, myID, stat(MPI_STATUS_SIZE)
      REAL(KIND=RKIND) :: Rt(dof,tnNo)
      REAL(KIND=RKIND), ALLOCATABLE :: sB(:,:,:), rB(:,:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: sReq(:), rReq(:)
      REAL(KIND=RKIND) FSILS_CPUT

      nNo  = hp%mynNo
      comm = lhs%commu%comm
      myID = lhs%commu%tF
      Ri   = 0._RKIND
      Rt   = 0._RKIND

      ls%callD = FSILS_CPUT() - ls%callD

      IF (ls%fNorm .LT. EPSILON(ls%fNorm)) THEN
         ls%dB = 0._RKIND
      ELSE
         ls%dB = 10._RKIND*LOG(ls%fNorm/ls%iNorm)
      END IF

      ! Reverse pre-vel separation
      DO i=1, nNo
         DO j=1, dof
            i0 = pvmap(j,i)
            hp%R(j,i) = lR(i0)
         END DO
      END DO

      IF (lhs%commu%nTasks .EQ. 1) THEN
         Ri = hp%R*Wc
         RETURN
      END IF

      ! Preconditioning
      DO i=1, nNo
         i0 = hp%ltg(i)
         Rt(:,i) = Wc(:,i0)*hp%R(:,i)
      END DO

      ! Communicate to fill Ri's outside of hp%mynNo.
      ! Processors that receive Val are sending Ri here.
      ALLOCATE(sB(dof,MAXVAL(hp%scS%n),hp%snReq))
      ALLOCATE(rB(dof,MAXVAL(hp%rcS%n),hp%rnReq))
      ALLOCATE(rReq(hp%rnReq),sReq(hp%snReq))

      DO i = 1, hp%rnReq
         DO j = 1, hp%rcS(i)%n
            k = hp%rcS(i)%ptr(j)
            rB(:,j,i) = Rt(:,k)
         END DO
      END DO

      DO i = 1, hp%rnReq
         CALL MPI_ISEND(rB(:,:,i), hp%rcS(i)%n*dof, mpreal, 
     2             hp%rcS(i)%iP-1, 1, comm, rReq(i), ierr)
      END DO
      DO i = 1, hp%snReq
         CALL MPI_IRECV(sB(:,:,i), hp%scS(i)%n*dof, mpreal,
     2             hp%scS(i)%iP-1, 1, comm, sReq(i), ierr)
      END DO
      DO i = 1, hp%snReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO
      DO i = 1, hp%rnReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO

      DO i = 1, hp%snReq
         DO j = 1, hp%scS(i)%n
            k = hp%scS(i)%ptr(j)
            Rt(:,k) = sB(:,j,i)
         END DO
      END DO

      ! lhs mapping
      DO i = 1, tnNo
         Ri(:,i) = Rt(:,lhs%map(i))
      END DO

      DEALLOCATE(rB, sB, rReq, sReq)

      RETURN
      END SUBROUTINE HYPRE_POST_SOLVE

!--------------------------------------------------------------------
!     Build sparse matrix in terms of individual dof
      SUBROUTINE HYPRE_DOFSPARSE(dof, gnNo, comm)
      USE HYPREMOD
      USE CMMOD
      IMPLICIT NONE 
      
      INTEGER(KIND=IKIND), INTENT(IN) :: dof, gnNo, comm
      
      INTEGER(KIND=IKIND) :: g2g(gnNo), temp(gnNo*dof), grmap(gnNo*dof)
      INTEGER(KIND=IKIND) :: i, j, k, a, ai, Ac, ierr
      INTEGER(KIND=IKIND) :: rs, re, ilower, iupper

      ! Special case for multiple equaitons
      IF (hp%init) DEALLOCATE(lrowPtr, lcolPtr, lVal, lR, pvmap)
      hp%dof = dof

      ! Global to global mapping for Hypre solve.
      ! hp%R assumes a new global index hp%ilower.
      ! We need to map current global index to new global index.
      ! g2g(trube global indx) = hypre global indx
      g2g = 0
      temp = 0
      DO i=1, hp%mynNo
         j = hp%ltg(i)
         temp(j) = i + hp%ilower - 1 
      END DO
      CALL MPI_ALLREDUCE(temp(1:gnNo), g2g, gnNo, mpint, MPI_SUM, 
     2   comm, ierr)

      ! Build sparse matrix
      ALLOCATE(lrowPtr(2,dof*hp%mynNo),lcolPtr(dof*dof*hp%nnz))
      ALLOCATE(lVal(dof*dof*hp%nnz), lR(dof*hp%mynNo))
      ALLOCATE(pvmap(dof,hp%mynNo))

      ! Mapping for local pre-vel separation
      pvmap = 0
      a = (dof-1)*hp%mynNo
      DO i = 1, hp%mynNo
         ai = (i-1)*(dof-1)
         DO j = 1, dof-1
            pvmap(j,i) = ai + j
         END DO
         pvmap(dof,i) = a + i
      END DO

      ! Mapping for global pre-vel separation
      ! This is needed since colPtr can exceed hp%mynNo
      ilower = (hp%ilower-1)*dof + 1 
      iupper = ilower + dof*hp%mynNo - 1
      grmap = 0
      temp  = 0
      temp(ilower:iupper) = RESHAPE(pvmap,(/hp%mynNo*dof/)) + ilower - 1
      CALL MPI_ALLREDUCE(temp, grmap, gnNo*dof, mpint, MPI_SUM, 
     2   comm, ierr)

      ! Build sparse matrix, with pres-vel separation
      lrowPtr = 0
      lcolPtr = 0
      DO i = 1, hp%mynNo
         rs = hp%rowPtr(1,i)
         re = hp%rowPtr(2,i)

         a = 0
         DO j = rs, re
            ai = (g2g(hp%colPtr(j))-1)*dof
            DO k = 1, dof
               a = a + 1
               temp(a) = grmap(ai + k)
            END DO
         END DO

         ! a  = (re-rs+1)*dof now
         ai = dof*dof*(rs-1)
         DO j = 1, dof
            Ac = pvmap(j,i)
            lrowPtr(1,Ac) = ai + (j-1)*a + 1
            lrowPtr(2,Ac) = lrowPtr(1,Ac) + a - 1
            lcolPtr(lrowPtr(1,Ac):lrowPtr(2,Ac)) = temp(1:a)
         END DO

      END DO

      RETURN 
      END SUBROUTINE HYPRE_DOFSPARSE

!--------------------------------------------------------------------
!     Row and column preconditioner to precondition both LHS and RHS.
      SUBROUTINE HYPRE_PRECONDRCS(lhs, rowPtr, colPtr, diagPtr,
     2   nNo, gnNo, nnz, dof, Val, R, ltg, W1, W2)

      USE TYPEMOD
      INCLUDE "FSILS_STD.h"
      
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=IKIND), INTENT(IN) :: dof, nNo, gnNo, nnz
      INTEGER(KIND=IKIND), INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz),     
     2                                   diagPtr(nNo), ltg(nNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: Val(dof*dof,nnz), 
     2   R(dof,nNo)
      REAL(KIND=RKIND), INTENT(OUT) :: W1(dof,nNo), W2(dof,gnNo)

      REAL(KIND=RKIND) :: Wr(dof,nNo), Wc(dof,gnNo), temp(dof,gnNo), tol
      INTEGER(KIND=IKIND) i, j, a, b, d, Ac, faIn
      INTEGER(KIND=IKIND) iter, maxiter, ierr
      LOGICAL flag, gflag(lhs%commu%nTasks)

      maxiter = 10
      tol     = 2._RKIND
      iter    = 0
      flag    = .TRUE.
      
      !*****************************************************
      ! Apply Dirichlet BC
      !*****************************************************
      Wr = 1._RKIND
      DO faIn=1, lhs%nFaces
         IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
         i = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               IF (Ac .LE. nNo)
     2            Wr(1:i,Ac) = Wr(1:i,Ac)*lhs%face(faIn)%val(1:i,a)
            END DO
         END IF
      END DO
      ! For parallel case, val and Wr can be larger than 1 due to
      ! the addition operator in FSILS_COMMUV. Hence need renormalization.
      Wr = Wr - 0.5_LSRP
      Wr = Wr/ABS(Wr)
      Wr = (Wr + ABS(Wr))*0.5_LSRP
      ! Kill the row and column corresponding to Dirichlet BC
      temp = 0._RKIND
      DO Ac=1, nNo
         temp(:,ltg(Ac)) = Wr(:,Ac)
      END DO
      CALL MPI_ALLREDUCE(temp, Wc, dof*gnNo, mpreal, MPI_SUM,
     2   lhs%commu%comm, ierr)
      CALL HYPRE_PREMUL(rowPtr, nNo, nnz, dof, Val, Wr)
      CALL HYPRE_POSMUL(rowPtr, colPtr, nNo, gnNo, nnz, dof, Val, Wc)
      R = Wr*R
      ! Set diagnal term to one
      SELECT CASE (dof)
      CASE (1)
         DO Ac=1, nNo
            d        = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._RKIND) + 1._RKIND 
         END DO
      CASE(2)
         DO Ac=1, nNo
            d        = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._RKIND) + 1._RKIND 
            Val(4,d) = Wr(2,Ac)*(Val(4,d)-1._RKIND) + 1._RKIND 
         END DO
      CASE(3)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            Val(1,d) = Wr(1,Ac)*(Val(1,d)-1._RKIND) + 1._RKIND 
            Val(5,d) = Wr(2,Ac)*(Val(5,d)-1._RKIND) + 1._RKIND 
            Val(9,d) = Wr(3,Ac)*(Val(9,d)-1._RKIND) + 1._RKIND 
         END DO
      CASE(4)
         DO Ac=1, nNo
            d       = diagPtr(Ac)
            Val(1 ,d) = Wr(1,Ac)*(Val(1 ,d)-1._RKIND) + 1._RKIND 
            Val(6 ,d) = Wr(2,Ac)*(Val(6 ,d)-1._RKIND) + 1._RKIND 
            Val(11,d) = Wr(3,Ac)*(Val(11,d)-1._RKIND) + 1._RKIND 
            Val(16,d) = Wr(4,Ac)*(Val(16,d)-1._RKIND) + 1._RKIND 
         END DO
      CASE DEFAULT
         DO Ac=1, nNo
            d = diagPtr(Ac)
            DO i=1, dof
               Val(i*dof-dof+i,d) = Wr(i,Ac)*(Val(i*dof-dof+i,d)
     2                            - 1._RKIND) + 1._RKIND
            END DO
         END DO
      END SELECT

      !*****************************************************
      ! Row and column scaling
      !*****************************************************
      W1 = 1._RKIND
      W2 = 1._RKIND
      DO While (flag)
         Wr   = 0._RKIND
         Wc   = 0._RKIND
         temp = 0._RKIND
         iter = iter + 1
         IF (iter .GE. maxiter) THEN 
            PRINT *, "Warning: maximum iteration number reached"//
     2            "@ SUBROUTINE PRECONDRCS."
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

         CALL MPI_ALLREDUCE(Wc, temp, dof*gnNo, mpreal, MPI_MAX,
     2      lhs%commu%comm, ierr)
         Wc = temp

         IF (MAXVAL(ABS(1._RKIND - Wr)) .LT. tol .AND. 
     2       MAXVAL(ABS(1._RKIND - Wc)) .LT. tol) flag = .False.

         Wr = 1._RKIND/SQRT(Wr)
         Wc = 1._RKIND/SQRT(Wc)

         CALL HYPRE_PREMUL(rowPtr, nNo, nnz, dof, Val, Wr)
         CALL HYPRE_POSMUL(rowPtr, colPtr, nNo, gnNo, nnz, dof, Val, Wc)

         W1 = W1*Wr
         W2 = W2*Wc

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
      END SUBROUTINE HYPRE_PRECONDRCS

!--------------------------------------------------------------------
!     Pre-multipling Val with W: Val = W*Val
      SUBROUTINE HYPRE_PREMUL(rowPtr, nNo, nnz, dof, Val, W)

      USE TYPEMOD

      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, nnz, dof
      INTEGER(KIND=IKIND), INTENT(IN) :: rowPtr(2,nNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: Val(dof*dof,nnz)
      REAL(KIND=RKIND), INTENT(IN) :: W(dof,nNo)

      INTEGER(KIND=IKIND) i, j, a, b, Ac
      
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
      END SUBROUTINE HYPRE_PREMUL

!--------------------------------------------------------------------
!     Post-multipling Val by W: Val = Val*W
      SUBROUTINE HYPRE_POSMUL(rowPtr, colPtr, nNo, gnNo, nnz, 
     2   dof, Val, W)

      USE TYPEMOD

      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, gnNo, nnz, dof
      INTEGER(KIND=IKIND), INTENT(IN) :: rowPtr(2,nNo),
     2                                   colPtr(nnz)
      REAL(KIND=RKIND), INTENT(INOUT) :: Val(dof*dof,nnz)
      REAL(KIND=RKIND), INTENT(IN) :: W(dof,gnNo)

      INTEGER(KIND=IKIND) i, j, a, b, Ac
      
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
      END SUBROUTINE HYPRE_POSMUL

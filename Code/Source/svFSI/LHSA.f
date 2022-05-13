!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
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
! permit persons to whom the Software is furnished to do so, subject
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
!     Creating the data structure and assembling LHS sparse matrix.
!
!--------------------------------------------------------------------

      SUBROUTINE LHSA(nnz)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: nnz

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, i, j, rowN, colN, iM, iFa, masN,
     2   mnnzeic

      INTEGER(KIND=IKIND), ALLOCATABLE :: uInd(:,:)

      ALLOCATE(idMap(tnNo))
      DO a=1, tnNo
         idMap(a) = a
      END DO

      mnnzeic = 10*MAXVAL(msh%eNoN)

!     First fill uInd array depending on mesh connectivity as is.
      ALLOCATE (uInd(mnnzeic,tnNo))
      uInd = 0
      DO iM=1, nMsh
!        Treat shell with triangular elements separately
         IF (shlEq .AND. msh(iM)%eType.EQ.eType_TRI3) CYCLE
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               rowN = msh(iM)%IEN(a,e)
               DO b=1, msh(iM)%eNoN
                  colN = msh(iM)%IEN(b,e)
                  CALL ADDCOL(rowN, colN)
               END DO
            END DO
         END DO
      END DO

!     Treat shells with triangular elements here
      DO iM=1, nMsh
         IF (.NOT.shlEq .OR. .NOT.msh(iM)%lShl) CYCLE
         IF (msh(iM)%eType .EQ. eType_NRB) CYCLE
         DO e=1, msh(iM)%nEl
            DO a=1, 2*msh(iM)%eNoN
               IF (a .LE. msh(iM)%eNoN) THEN
                  rowN = msh(iM)%IEN(a,e)
               ELSE
                  rowN = msh(iM)%eIEN(a-msh(iM)%eNoN,e)
               END IF
               IF (rowN .EQ. 0) CYCLE
               DO b=1, 2*msh(iM)%eNoN
                  IF (b .LE. msh(iM)%eNoN) THEN
                     colN = msh(iM)%IEN(b,e)
                  ELSE
                     colN = msh(iM)%eIEN(b-msh(iM)%eNoN,e)
                  END IF
                  IF (colN .EQ. 0) CYCLE
                  CALL ADDCOL(rowN, colN)
               END DO
            END DO
         END DO
      END DO

!     Now reset idMap for undeforming Neumann BC faces. Then insert
!     master node as a column entry in each row for all the slave nodes.
!     This step is performed even for ghost master nodes where the idMap
!     points to the ghost master node.
      flag = .FALSE.
      DO i=1, nEq
         DO j=1, eq(i)%nBc
            iM  = eq(i)%bc(j)%iM
            iFa = eq(i)%bc(j)%iFa
            IF (BTEST(eq(i)%bc(j)%bType, bType_undefNeu)) THEN
               masN = eq(i)%bc(j)%masN
               IF (masN .EQ. 0) CYCLE
               DO a=1, msh(iM)%fa(iFa)%nNo
                  rowN = msh(iM)%fa(iFa)%gN(a)
                  IF (rowN .EQ. masN) CYCLE
                  idMap(rowN) = masN
!                 Insert master to the row if not already present
                  CALL ADDCOL(rowN, masN)
               END DO
               flag = .TRUE.
            END IF
         END DO
      END DO

!     Change uInd if idMap has been changed
      IF (flag) THEN
         DO a=1, tnNo
            rowN = idMap(a)
!           If the mapping is not changed, examine the mapping of the
!           column entries of uInd. Don't do anything if the mapping is
!           unchanged. If the mapping is changed, then add to the column\
!           indices of uInd if the entry is not already present.
            IF (rowN .EQ. a) THEN
               i = 0
               DO
                  i = i + 1
                  b = uInd(i,rowN)
!                 Terminate if end of column entries are reached
                  IF (b .EQ. 0) EXIT
                  colN = idMap(b)
!                 Ignore if the column entry mapping is not changed.
!                 This entry is already present and will be used to
!                 assemble mass (D) matrix
                  IF (b .EQ. colN) CYCLE
!                 As the column entry is now mapped to a new node,
!                 search all column entries and insert the new node if
!                 it is not present. This step is performed to assemble
!                 the Divergence (C) matrix
                  CALL ADDCOL(rowN, colN)
                  IF (i .EQ. mnnzeic) EXIT
               END DO
            ELSE
!           If the row mapping is changed, insert the mapped/unmapped
!           column entries of the old row into the new row if not
!           already present.
               i = 0
               DO
                  i = i + 1
                  b = uInd(i,a)
!                 Terminate if end of column entries are reached
                  IF (b .EQ. 0) EXIT
!                 Add unmapped column to assemble gradient matrix
                  CALL ADDCOL(rowN, b)
!                 If column is mapped, add the mapped column to assemble
!                 stiffness matrix
                  colN = idMap(b)
                  IF (b .NE. colN) CALL ADDCOL(rowN, colN)
                  IF (i .EQ. mnnzeic) EXIT
               END DO
            END IF
         END DO
      END IF

!--------------------------------------------------------------------
!     Finding number of non-zeros in colPtr vector
      nnz = 0
      DO rowN=1, tnNo
         IF (uInd(1,rowN) .EQ. 0) THEN
            err = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               nnz = nnz + 1
            END IF
         END DO
      END DO

!--------------------------------------------------------------------
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE (colPtr(nnz), rowPtr(tnNo+1))
      j  = 1
      rowPtr(1) = 1
      DO rowN=1, tnNo
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE ADDCOL(row, col)
         IMPLICIT NONE
         INTEGER(KIND=IKIND), INTENT(IN) :: row, col

         INTEGER(KIND=IKIND) i, j

         i = 0
         DO
            i = i + 1
            IF (i .EQ. mnnzeic) CALL RESIZ()

!           If current entry is zero, then  fill it up
            IF (uInd(i,row) .EQ. 0) THEN
               uInd(i,row) = col
               EXIT
            END IF

!           If current entry is still smaller, keep going
            IF (col .GT. uInd(i,row)) CYCLE

!           If column entry already exists, exit
            IF (col .EQ. uInd(i,row)) EXIT

!           If we are this point, then then the current entry is bigger.
!           Shift all the entries from here to the end of the list. If
!           list is full, we request a larger list, otherwise we shift
!           and add the item at the current entry position.
            IF (uInd(mnnzeic,row) .NE. 0) CALL RESIZ()
            DO j=mnnzeic, i+1, -1
               uInd(j,row) = uInd(j-1,row)
            END DO
            uInd(i,row) = col
            EXIT
         END DO

         RETURN
         END SUBROUTINE ADDCOL
!--------------------------------------------------------------------
         SUBROUTINE RESIZ()
         IMPLICIT NONE

         INTEGER(KIND=IKIND) n
         INTEGER(KIND=IKIND), ALLOCATABLE :: tmp(:,:)

         n = mnnzeic
         ALLOCATE(tmp(n,tnNo))
         tmp(:,:) = uInd(:,:)
         DEALLOCATE(uInd)
         mnnzeic = n + MAX(5,n/5)
         ALLOCATE(uInd(mnnzeic,tnNo))
         uInd(:,:)   = 0
         uInd(1:n,:) = tmp(:,:)
         DEALLOCATE(tmp)

         RETURN
         END SUBROUTINE RESIZ
!--------------------------------------------------------------------
      END SUBROUTINE LHSA
!####################################################################
!     This subroutine assembels the element stiffness matrix into the
!     global stiffness matrix (Val sparse matrix formatted as a vector)
!     Also assembles residual vector into global residual vector
      SUBROUTINE DOASSEM (d, eqN, lK, lR)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, rowPtr, colPtr, R, Val
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqN(d)
      REAL(KIND=RKIND), INTENT(IN) :: lK(dof*dof,d,d), lR(dof,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqN(a)
         IF (rowN .EQ. 0) CYCLE
         R(:,rowN) = R(:,rowN) + lR(:,a)
         DO b=1, d
            colN = eqN(b)
            IF (colN .EQ. 0) CYCLE
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Val(:,ptr) = Val(:,ptr) + lK(:,a,b)
         END DO
      END DO

      RETURN
      END SUBROUTINE DOASSEM
!####################################################################

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
!     All the routines that require an interface are included here.
!     This Mainly involves small routines with a very well-defined
!     functionality. To use these routines, just add "USE ALLFUN" to
!     your routine.
!
!--------------------------------------------------------------------
      
!####################################################################
!     Matrix and tensor operations
      MODULE MATFUN
      
      IMPLICIT NONE
      
      INTEGER, ALLOCATABLE :: t_ind(:,:)
      
      CONTAINS
!--------------------------------------------------------------------
!     Create a second order identity matrix of rank nd
      FUNCTION MAT_ID(nd) RESULT(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8) :: A(nd,nd)
      
      INTEGER :: i
      
      A = 0D0
      DO i=1, nd
         A(i,i) = 1D0
      END DO
      
      RETURN
      END FUNCTION MAT_ID
!--------------------------------------------------------------------
!     Trace of second order matrix of rank nd
      FUNCTION MAT_TRACE(A, nd) RESULT(trA)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      REAL(KIND=8) :: trA
      
      INTEGER :: i
      
      trA = 0D0
      DO i=1, nd
         trA = trA + A(i,i)
      END DO
      
      RETURN
      END FUNCTION MAT_TRACE
!--------------------------------------------------------------------
!     Create a matrix from outer product of two vectors
      FUNCTION MAT_DYADPROD(u, v, nd) RESULT(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: u(nd), v(nd)
      REAL(KIND=8) :: A(nd,nd)
      
      INTEGER :: i, j
      
      DO j=1, nd
         DO i=1, nd
            A(i,j) = u(i)*v(j)
         END DO
      END DO
      
      RETURN
      END FUNCTION MAT_DYADPROD
!--------------------------------------------------------------------
!     Computes the determinant of a square matrix
      RECURSIVE FUNCTION MAT_DET(A, nd) RESULT(D)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      
      INTEGER :: i, j, n
      REAL(KIND=8) :: D, Am(nd-1,nd-1)
      
      D = 0D0
      IF (nd .EQ. 2) THEN
         D = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      ELSE
         DO i=1, nd
            n = 0
            DO j=1, nd
               IF (i .EQ. j) THEN
                  CYCLE
               ELSE
                  n = n + 1
                  Am(:,n) = A(2:nd,j)
               END IF
            END DO ! j
            D = D + ( (-1D0)**REAL(1+i,KIND=8) * A(1,i) *
     2                 MAT_DET(Am,nd-1) )
         END DO ! i
      END IF ! nd.EQ.2
      
      RETURN
      END FUNCTION MAT_DET
!--------------------------------------------------------------------
!     This function uses LAPACK library to compute eigen values of a 
!     square matrix, A of dimension nd
      FUNCTION MAT_EIG(A, nd)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd)
      COMPLEX*16 :: Amat(nd,nd), MAT_EIG(nd), b(nd), DUMMY(1,1),
     2   WORK(2*nd)
      
      INTEGER :: i, j, iok
      
      Amat = (0D0, 0D0)
      DO j=1, nd
         DO i=1, nd
            Amat(i,j) = CMPLX(A(i,j))
         END DO
      END DO
      
      CALL ZGEEV('N', 'N', nd, Amat, nd, b, DUMMY, 1, DUMMY, 1, WORK,
     2   2*nd, WORK, iok)
      IF (iok .NE. 0) THEN
         WRITE(*,'(A)') "Failed to compute eigen values"
         CALL STOPSIM()
      ELSE
         MAT_EIG(:) = b(:)
      END IF
      
      RETURN
      END FUNCTION MAT_EIG
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix using 
!     Gauss Elimination method
      FUNCTION MAT_INV(A, nd)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(:,:)
      REAL(KIND=8), ALLOCATABLE :: MAT_INV(:,:)
      
      REAL(KIND=8), PARAMETER :: epsil = EPSILON(epsil)
      
      INTEGER :: i, j, k
      REAL(KIND=8) :: d, B(nd,2*nd)
      
      d = MAT_DET(A, nd)
      IF (ABS(d) .LT. 1D2*epsil) THEN
         WRITE(*,'(A)') "Singular matrix detected to compute inverse"
         CALL STOPSIM()
      END IF
      
      ALLOCATE(MAT_INV(nd,nd))
      
      IF (nd .EQ. 2) THEN
         MAT_INV(1,1) =  A(2,2)/d
         MAT_INV(1,2) = -A(1,2)/d
         
         MAT_INV(2,1) = -A(2,1)/d
         MAT_INV(2,2) =  A(1,1)/d
         
         RETURN
      ELSE IF (nd .EQ. 3) THEN
         MAT_INV(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2)) / d
         MAT_INV(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3)) / d
         MAT_INV(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2)) / d
         
         MAT_INV(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3)) / d
         MAT_INV(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1)) / d
         MAT_INV(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3)) / d
         
         MAT_INV(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1)) / d
         MAT_INV(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2)) / d
         MAT_INV(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / d
         
         RETURN
      END IF
      
!     Auxillary matrix
      B = 0D0
      DO i=1, nd
         DO j=1, nd
            B(i,j) = A(i,j)
         END DO
         B(i,nd+i) = 1D0
      END DO
      
!     Pivoting
      DO i=nd, 2, -1
         IF (B(i,1) .GT. B(i-1,1)) THEN
            DO j=1, 2*nd
               d = B(i,j)
               B(i,j) = B(i-1,j)
               B(i-1,j) = d
            END DO
         END IF
      END DO
      
!     Do row-column operations and reduce to diagonal
      DO i=1, nd
         DO j=1, nd
            IF (j .NE. i) THEN
               d = B(j,i)/B(i,i)
               DO k=1, 2*nd
                  B(j,k) = B(j,k) - d*B(i,k)
               END DO
            END IF
         END DO
      END DO
      
!     Unit matrix
      DO i=1, nd
         d = B(i,i)
         DO j=1, 2*nd
            B(i,j) = B(i,j)/d
         END DO
      END DO
      
!     Inverse
      DO i=1, nd
         DO j=1, nd
            MAT_INV(i,j) = B(i,j+nd)
         END DO
      END DO
      
      RETURN
      END FUNCTION MAT_INV
!--------------------------------------------------------------------
!     Initialize tensor index pointer
      SUBROUTINE TEN_INIT(nd)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      
      INTEGER :: ii, nn, i, j, k, l
      
      nn = nd**4
      IF (ALLOCATED(t_ind)) DEALLOCATE(t_ind)
      ALLOCATE(t_ind(4,nn))
      
      ii = 0
      DO l=1, nd
         DO k=1, nd
            DO j=1, nd
               DO i=1, nd
                  ii = ii + 1
                  t_ind(1,ii) = i
                  t_ind(2,ii) = j
                  t_ind(3,ii) = k
                  t_ind(4,ii) = l
               END DO
            END DO
         END DO
      END DO
      
      RETURN
      END SUBROUTINE TEN_INIT
!--------------------------------------------------------------------
!     Create a 4th order order symmetric identity tensor
      FUNCTION TEN_IDs(nd) RESULT(A)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8) :: A(nd,nd,nd,nd)
      
      INTEGER :: i, j
      
      A = 0D0
      DO j=1, nd
         DO i=1, nd
            A(i,j,i,j) = A(i,j,i,j) + 5D-1
            A(i,j,j,i) = A(i,j,j,i) + 5D-1
         END DO
      END DO
      
      RETURN
      END FUNCTION TEN_IDs
!--------------------------------------------------------------------
!     Create a 4th order tensor from outer product of two matrices
      FUNCTION TEN_DYADPROD(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd), B(nd,nd)
      REAL(KIND=8) :: C(nd,nd,nd,nd)
      
      INTEGER :: ii, nn, i, j, k, l
      
      nn = nd**4
      DO ii=1, nn
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = A(i,j) * B(k,l)
      END DO
      
      RETURN
      END FUNCTION TEN_DYADPROD
!--------------------------------------------------------------------
!     Create a 4th order tensor from symmetric outer product of two
!     matrices
      FUNCTION TEN_SYMMPROD(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nd
      REAL(KIND=8), INTENT(IN) :: A(nd,nd), B(nd,nd)
      REAL(KIND=8) :: C(nd,nd,nd,nd)
      
      INTEGER :: ii, nn, i, j, k, l
      
      nn = nd**4
      DO ii=1, nn
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = 5D-1 * ( A(i,k)*B(j,l) + A(i,l)*B(j,k) )
      END DO
      
      RETURN
      END FUNCTION TEN_SYMMPROD
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on first index (2D)
      FUNCTION TEN_MATMULT2D_1(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(2,2), T(2,2,2,2)
      REAL(KIND=8) :: C(2,2,2,2)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 16
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(i,1)*T(1,j,k,l) + A(i,2)*T(2,j,k,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT2D_1
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on second index (2D)
      FUNCTION TEN_MATMULT2D_2(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(2,2), T(2,2,2,2)
      REAL(KIND=8) :: C(2,2,2,2)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 16
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(j,1)*T(i,1,k,l) + A(j,2)*T(i,2,k,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT2D_2
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on third index (2D)
      FUNCTION TEN_MATMULT2D_3(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(2,2), T(2,2,2,2)
      REAL(KIND=8) :: C(2,2,2,2)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 16
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(k,1)*T(i,j,1,l) + A(k,2)*T(i,j,2,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT2D_3
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on fourth index (2D)
      FUNCTION TEN_MATMULT2D_4(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(2,2), T(2,2,2,2)
      REAL(KIND=8) :: C(2,2,2,2)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 16
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(l,1)*T(i,j,k,1) + A(l,2)*T(i,j,k,2)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT2D_4
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on first index (3D)
      FUNCTION TEN_MATMULT3D_1(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(3,3), T(3,3,3,3)
      REAL(KIND=8) :: C(3,3,3,3)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 81
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(i,1)*T(1,j,k,l) + 
     2      A(i,2)*T(2,j,k,l) + A(i,3)*T(3,j,k,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT3D_1
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on second index (3D)
      FUNCTION TEN_MATMULT3D_2(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(3,3), T(3,3,3,3)
      REAL(KIND=8) :: C(3,3,3,3)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 81
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(j,1)*T(i,1,k,l) + 
     2      A(j,2)*T(i,2,k,l) + A(j,3)*T(i,3,k,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT3D_2
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on third index (3D)
      FUNCTION TEN_MATMULT3D_3(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(3,3), T(3,3,3,3)
      REAL(KIND=8) :: C(3,3,3,3)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 81
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(k,1)*T(i,j,1,l) + 
     2      A(k,2)*T(i,j,2,l) + A(k,3)*T(i,j,3,l)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT3D_3
!--------------------------------------------------------------------
!     Matrix-Tensor multiplication with contraction on fourth index (3D)
      FUNCTION TEN_MATMULT3D_4(A, T) RESULT(C)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: A(3,3), T(3,3,3,3)
      REAL(KIND=8) :: C(3,3,3,3)
      
      INTEGER :: ii, i, j, k, l
      
      C  = 0D0
      DO ii=1, 81
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = C(i,j,k,l) + A(l,1)*T(i,j,k,1) + 
     2      A(l,2)*T(i,j,k,2) + A(l,3)*T(i,j,k,3)
      END DO
      
      RETURN
      END FUNCTION TEN_MATMULT3D_4
!--------------------------------------------------------------------
      SUBROUTINE MAT_DISP(A, m, n)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m
      INTEGER, INTENT(IN), OPTIONAL :: n
      REAL(KIND=8), INTENT(IN) :: A(:,:)
      
      INTEGER :: i, j, nr, nc
      
      nr = m
      nc = nr
      IF (PRESENT(n)) nc = n
      
      DO i=1, nr
         WRITE(*,'(3X,A)',ADVANCE='NO')
         DO j=1, nc
            WRITE(*,'(X,1pE15.6,X,A)',ADVANCE='NO') A(i,j), ' | '
         END DO
         WRITE(*,'(A)')
      END DO
      
      RETURN
      END SUBROUTINE MAT_DISP
!--------------------------------------------------------------------
      END MODULE MATFUN
!####################################################################

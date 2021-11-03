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
!     This module defines generally performed matrix and tensor
!     operations.
!
!--------------------------------------------------------------------

!     Matrix and tensor operations
      MODULE MATFUN
      USE TYPEMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND), ALLOCATABLE :: t_ind(:,:)

      PRIVATE :: ISZERO

      INTERFACE TEN_DDOT
         MODULE PROCEDURE TEN_DDOT_3434
      END INTERFACE TEN_DDOT

      CONTAINS
!--------------------------------------------------------------------
!     Create a second order identity matrix of rank nd
      FUNCTION MAT_ID(nd) RESULT(A)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND) :: A(nd,nd)

      INTEGER(KIND=IKIND) :: i

      A = 0._RKIND
      DO i=1, nd
         A(i,i) = 1._RKIND
      END DO

      RETURN
      END FUNCTION MAT_ID
!--------------------------------------------------------------------
!     Trace of second order matrix of rank nd
      FUNCTION MAT_TRACE(A, nd) RESULT(trA)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      REAL(KIND=RKIND) :: trA

      INTEGER(KIND=IKIND) :: i

      trA = 0._RKIND
      DO i=1, nd
         trA = trA + A(i,i)
      END DO

      RETURN
      END FUNCTION MAT_TRACE
!--------------------------------------------------------------------
!     Symmetric part of a matrix, S = (A + A.T)/2
      FUNCTION MAT_SYMM(A, nd) RESULT(S)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      REAL(KIND=RKIND) :: S(nd,nd)

      INTEGER(KIND=IKIND) :: i, j

      S = 0._RKIND
      DO i=1, nd
         DO j=1, nd
            S(i,j) = 0.5_RKIND*(A(i,j) + A(j,i))
         END DO
      END DO

      RETURN
      END FUNCTION MAT_SYMM
!--------------------------------------------------------------------
!     Deviatoric part of a matrix, D = A - (tr.(A)/nd) I
      FUNCTION MAT_DEV(A, nd) RESULT(D)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      REAL(KIND=RKIND) :: D(nd,nd)

      INTEGER(KIND=IKIND) :: trA

      trA = MAT_TRACE(A,nd)
      D   = A - (trA/REAL(nd,KIND=RKIND)) * MAT_ID(nd)

      RETURN
      END FUNCTION MAT_DEV
!--------------------------------------------------------------------
!     Create a matrix from outer product of two vectors
      FUNCTION MAT_DYADPROD(u, v, nd) RESULT(A)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: u(nd), v(nd)
      REAL(KIND=RKIND) :: A(nd,nd)

      INTEGER(KIND=IKIND) :: i, j

      DO j=1, nd
         DO i=1, nd
            A(i,j) = u(i)*v(j)
         END DO
      END DO

      RETURN
      END FUNCTION MAT_DYADPROD
!--------------------------------------------------------------------
!     Create a matrix from symmetric product of two vectors
      FUNCTION MAT_SYMMPROD(u, v, nd) RESULT(A)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: u(nd), v(nd)
      REAL(KIND=RKIND) :: A(nd,nd)

      INTEGER(KIND=IKIND) :: i, j

      DO j=1, nd
         DO i=1, nd
            A(i,j) = 0.5_RKIND*(u(i)*v(j) + u(j)*v(i))
         END DO
      END DO

      RETURN
      END FUNCTION MAT_SYMMPROD
!--------------------------------------------------------------------
!     Double dot product of 2 square matrices
      FUNCTION MAT_DDOT(A, B, nd) RESULT(s)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd), B(nd,nd)
      REAL(KIND=RKIND) :: s

      INTEGER(KIND=IKIND) :: i, j

      s = 0._RKIND
      DO j=1, nd
         DO i=1, nd
            s = s + A(i,j) * B(i,j)
         END DO
      END DO

      RETURN
      END FUNCTION MAT_DDOT
!--------------------------------------------------------------------
!     Computes the determinant of a square matrix
      RECURSIVE FUNCTION MAT_DET(A, nd) RESULT(D)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)

      INTEGER(KIND=IKIND) :: i, j, n
      REAL(KIND=RKIND) :: D, Am(nd-1,nd-1)

      D = 0._RKIND
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
            D = D + ( (-1._RKIND)**REAL(1+i, KIND=RKIND) * A(1,i) *
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
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      COMPLEX(KIND=CXKIND) :: Amat(nd,nd), MAT_EIG(nd), b(nd),
     2   DUMMY(1,1), WORK(2*nd)

      INTEGER(KIND=IKIND) :: i, j, iok

      Amat = (0._RKIND, 0._RKIND)
      DO j=1, nd
         DO i=1, nd
            Amat(i,j) = CMPLX(A(i,j), 0._RKIND, KIND=CXKIND)
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
!     This function computes inverse of a square matrix
      FUNCTION MAT_INV(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      REAL(KIND=RKIND) :: Ainv(nd,nd)

      REAL(KIND=RKIND), PARAMETER :: epsil = EPSILON(epsil)

      INTEGER(KIND=IKIND) :: iok = 0
      REAL(KIND=RKIND) :: d

      IF (nd .EQ. 2) THEN
         d = MAT_DET(A, nd)
         IF (ISZERO(ABS(d))) iok = -1

         Ainv(1,1) =  A(2,2)/d
         Ainv(1,2) = -A(1,2)/d

         Ainv(2,1) = -A(2,1)/d
         Ainv(2,2) =  A(1,1)/d

      ELSE IF (nd .EQ. 3) THEN
         d = MAT_DET(A, nd)
         IF (ISZERO(ABS(d))) iok = -1

         Ainv(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2)) / d
         Ainv(1,2) = (A(1,3)*A(3,2)-A(1,2)*A(3,3)) / d
         Ainv(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2)) / d

         Ainv(2,1) = (A(2,3)*A(3,1)-A(2,1)*A(3,3)) / d
         Ainv(2,2) = (A(1,1)*A(3,3)-A(1,3)*A(3,1)) / d
         Ainv(2,3) = (A(1,3)*A(2,1)-A(1,1)*A(2,3)) / d

         Ainv(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1)) / d
         Ainv(3,2) = (A(1,2)*A(3,1)-A(1,1)*A(3,2)) / d
         Ainv(3,3) = (A(1,1)*A(2,2)-A(1,2)*A(2,1)) / d

      ELSE IF (nd.GT.3 .AND. nd.LT.10) THEN
         d = MAT_DET(A, nd)
         IF (ISZERO(ABS(d))) iok = -1
         Ainv = MAT_INV_GE(A, nd)

      ELSE
         Ainv = MAT_INV_LP(A, nd)

      END IF

      IF (iok .ne. 0) THEN
         WRITE(*,'(A)') "Singular matrix detected to compute inverse"
         WRITE(*,'(A)') "ERROR: Matrix inversion failed"
         CALL STOPSIM()
      END IF

      RETURN
      END FUNCTION MAT_INV
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix using
!     Gauss Elimination method
      FUNCTION MAT_INV_GE(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd)
      REAL(KIND=RKIND) :: Ainv(nd,nd)

      INTEGER(KIND=IKIND) :: i, j, k
      REAL(KIND=RKIND) :: d, B(nd,2*nd)

!     Auxillary matrix
      B = 0._RKIND
      DO i=1, nd
         DO j=1, nd
            B(i,j) = A(i,j)
         END DO
         B(i,nd+i) = 1._RKIND
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
            Ainv(i,j) = B(i,j+nd)
         END DO
      END DO

      RETURN
      END FUNCTION MAT_INV_GE
!--------------------------------------------------------------------
!     This function computes inverse of a square matrix using
!     Lapack functions (DGETRF + DGETRI)
      FUNCTION MAT_INV_LP(A, nd) RESULT(Ainv)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd, nd)
      REAL(KIND=RKIND) :: Ainv(nd,nd)

      INTEGER(KIND=IKIND) iok, IPIV(nd)
      REAL(KIND=RKIND) :: Ad(nd,nd), WORK(2*nd)

      Ad = A
      CALL DGETRF(nd, nd, Ad, nd, IPIV, iok)
      IF (iok .NE. 0) THEN
         WRITE(*,'(A)') "Singular matrix detected to compute inverse"
         WRITE(*,'(A)') "ERROR: Matrix inversion failed"
         CALL STOPSIM()
      END IF

      CALL DGETRI(nd, Ad, nd, IPIV, WORK, 2*nd, iok)
      IF (iok .NE. 0) THEN
         WRITE(*,'(A)') "ERROR: Matrix inversion failed (LAPACK)"
         CALL STOPSIM()
      END IF

      Ainv(:,:) = Ad(:,:)

      RETURN
      END FUNCTION MAT_INV_LP
!--------------------------------------------------------------------
!     Initialize tensor index pointer
      SUBROUTINE TEN_INIT(nd)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

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
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND) :: A(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: i, j

      A = 0._RKIND
      DO j=1, nd
         DO i=1, nd
            A(i,j,i,j) = A(i,j,i,j) + 0.5_RKIND
            A(i,j,j,i) = A(i,j,j,i) + 0.5_RKIND
         END DO
      END DO

      RETURN
      END FUNCTION TEN_IDs
!--------------------------------------------------------------------
!     Create a 4th order tensor from outer product of two matrices
      FUNCTION TEN_DYADPROD(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd), B(nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

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
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd), B(nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

      nn = nd**4
      DO ii=1, nn
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         C(i,j,k,l) = 0.5_RKIND * ( A(i,k)*B(j,l) + A(i,l)*B(j,k) )
      END DO

      RETURN
      END FUNCTION TEN_SYMMPROD
!--------------------------------------------------------------------
!     Transpose of a 4th order tensor [A^T]_ijkl = [A]_klij
      FUNCTION TEN_TRANSPOSE(A, nd) RESULT(B)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd,nd,nd)
      REAL(KIND=RKIND) :: B(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

      nn = nd**4
      DO ii=1, nn
         i = t_ind(1,ii)
         j = t_ind(2,ii)
         k = t_ind(3,ii)
         l = t_ind(4,ii)
         B(i,j,k,l) = A(k,l,i,j)
      END DO

      RETURN
      END FUNCTION TEN_TRANSPOSE
!--------------------------------------------------------------------
!     Double dot product of a 4th order tensor and a 2nd order tensor
!     C_ij = (A_ijkl * B_kl)
      FUNCTION TEN_MDDOT(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd,nd,nd), B(nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd)

      INTEGER(KIND=IKIND) :: i, j

      IF (nd .EQ. 2) THEN
         DO i=1, nd
            DO j=1, nd
               C(i,j) = A(i,j,1,1)*B(1,1) + A(i,j,1,2)*B(1,2)
     2                + A(i,j,2,1)*B(2,1) + A(i,j,2,2)*B(2,2)
            END DO
         END DO
      ELSE
         DO i=1, nd
            DO j=1, nd
               C(i,j) = A(i,j,1,1)*B(1,1) + A(i,j,1,2)*B(1,2)
     2                + A(i,j,1,3)*B(1,3) + A(i,j,2,1)*B(2,1)
     3                + A(i,j,2,2)*B(2,2) + A(i,j,2,3)*B(2,3)
     4                + A(i,j,3,1)*B(3,1) + A(i,j,3,2)*B(3,2)
     5                + A(i,j,3,3)*B(3,3)
            END DO
         END DO
      END IF

      RETURN
      END FUNCTION TEN_MDDOT
!--------------------------------------------------------------------
!     Double dot product of 2 4th order tensors
!     T_ijkl = A_ijmn * B_klmn
      FUNCTION TEN_DDOT_3434(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd,nd,nd), B(nd,nd,nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

      C  = 0._RKIND
      nn = nd**4
      IF (nd .EQ. 2) THEN
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,j,1,1)*B(k,l,1,1)
     2                              + A(i,j,1,2)*B(k,l,1,2)
     3                              + A(i,j,2,1)*B(k,l,2,1)
     4                              + A(i,j,2,2)*B(k,l,2,2)
         END DO
      ELSE
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,j,1,1)*B(k,l,1,1)
     2                              + A(i,j,1,2)*B(k,l,1,2)
     3                              + A(i,j,1,3)*B(k,l,1,3)
     4                              + A(i,j,2,1)*B(k,l,2,1)
     5                              + A(i,j,2,2)*B(k,l,2,2)
     6                              + A(i,j,2,3)*B(k,l,2,3)
     7                              + A(i,j,3,1)*B(k,l,3,1)
     8                              + A(i,j,3,2)*B(k,l,3,2)
     9                              + A(i,j,3,3)*B(k,l,3,3)
         END DO
      END IF

      RETURN
      END FUNCTION TEN_DDOT_3434
!--------------------------------------------------------------------
!     T_ijkl = A_ijmn * B_kmln
      FUNCTION TEN_DDOT_3424(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd,nd,nd), B(nd,nd,nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

      C  = 0._RKIND
      nn = nd**4
      IF (nd .EQ. 2) THEN
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,j,1,1)*B(k,1,l,1)
     2                              + A(i,j,1,2)*B(k,1,l,2)
     3                              + A(i,j,2,1)*B(k,2,l,1)
     4                              + A(i,j,2,2)*B(k,2,l,2)
         END DO
      ELSE
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,j,1,1)*B(k,1,l,1)
     2                              + A(i,j,1,2)*B(k,1,l,2)
     3                              + A(i,j,1,3)*B(k,1,l,3)
     4                              + A(i,j,2,1)*B(k,2,l,1)
     5                              + A(i,j,2,2)*B(k,2,l,2)
     6                              + A(i,j,2,3)*B(k,2,l,3)
     7                              + A(i,j,3,1)*B(k,3,l,1)
     8                              + A(i,j,3,2)*B(k,3,l,2)
     9                              + A(i,j,3,3)*B(k,3,l,3)
         END DO
      END IF

      RETURN
      END FUNCTION TEN_DDOT_3424
!--------------------------------------------------------------------
!     T_ijkl = A_imjn * B_mnkl
      FUNCTION TEN_DDOT_2412(A, B, nd) RESULT(C)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nd
      REAL(KIND=RKIND), INTENT(IN) :: A(nd,nd,nd,nd), B(nd,nd,nd,nd)
      REAL(KIND=RKIND) :: C(nd,nd,nd,nd)

      INTEGER(KIND=IKIND) :: ii, nn, i, j, k, l

      C  = 0._RKIND
      nn = nd**4
      IF (nd .EQ. 2) THEN
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,1,j,1)*B(1,1,k,l)
     2                              + A(i,1,j,2)*B(1,2,k,l)
     3                              + A(i,2,j,1)*B(2,1,k,l)
     4                              + A(i,2,j,2)*B(2,2,k,l)
         END DO
      ELSE
         DO ii=1, nn
            i = t_ind(1,ii)
            j = t_ind(2,ii)
            k = t_ind(3,ii)
            l = t_ind(4,ii)
            C(i,j,k,l) = C(i,j,k,l) + A(i,1,j,1)*B(1,1,k,l)
     2                              + A(i,1,j,2)*B(1,2,k,l)
     3                              + A(i,1,j,3)*B(1,3,k,l)
     4                              + A(i,2,j,1)*B(2,1,k,l)
     5                              + A(i,2,j,2)*B(2,2,k,l)
     6                              + A(i,2,j,3)*B(2,3,k,l)
     7                              + A(i,3,j,1)*B(3,1,k,l)
     8                              + A(i,3,j,2)*B(3,2,k,l)
     9                              + A(i,3,j,3)*B(3,3,k,l)
         END DO
      END IF

      RETURN
      END FUNCTION TEN_DDOT_2412
!--------------------------------------------------------------------
      FUNCTION ISZERO(ia)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: ia
      LOGICAL ISZERO

      REAL(KIND=RKIND), PARAMETER :: epsil = EPSILON(epsil)
      REAL(KIND=RKIND) a, b, nrm

      a   = ABS(ia)
      b   = 0._RKIND
      nrm = MAX(a,epsil)

      ISZERO = .FALSE.
      IF ((a-b)/nrm .LT. 10._RKIND*epsil) ISZERO = .TRUE.

      RETURN
      END FUNCTION ISZERO
!--------------------------------------------------------------------
      END MODULE MATFUN
!####################################################################

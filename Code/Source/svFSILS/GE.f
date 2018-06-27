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
!     This subroutine solve a linear system of equations AX=B using 
!     Gauss elimination and replaces the B with X
!--------------------------------------------------------------------
      
      FUNCTION GE (nV, N, A, B)
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nV, N
      REAL(KIND=8), INTENT(IN) :: A(nV,N)
      REAL(KIND=8), INTENT(INOUT) :: B(N)
      LOGICAL GE

      INTEGER m, ipv, i, j
      REAL(KIND=8) pivot, saveEl
      REAL(KIND=8), ALLOCATABLE :: C(:,:), W(:)

!     Constructing a preconditioner. This is to prevent latter problem
!     with singular matrix
      ALLOCATE(W(N))
      DO i=1, N
         IF (ABS(A(i,i)) .LT. TINY(A(1,1))) THEN
            B  = 0D0
            GE = .FALSE. 
            RETURN
         END IF
         W(i) = 1D0/SQRT(ABS(A(i,i)))
      END DO

      ALLOCATE(C(N,N+1))
      DO i=1, N
         DO j=1, N
            C(i,j) = W(i)*W(j)*A(i,j)
         END DO
         C(i,N+1) = W(i)*B(i)
      END DO

      GE = .TRUE.
      IF (N .LE. 0) THEN
         GE = .FALSE. 
         RETURN
      ELSE IF (N .EQ. 1) THEN
         B(1) = C(1,2)/C(1,1)
         B(1) = B(1)*W(1)
         RETURN
      ELSE IF (N .EQ. 2) THEN
         pivot = C(1,1)*C(2,2) - C(2,1)*C(1,2)
         IF (ABS(pivot) .LT. EPSILON(pivot)) THEN
!     Singular matrix            
            B  = 0D0
            GE = .FALSE.
            RETURN
         END IF
         B(1) = (C(1,3)*C(2,2) - C(2,3)*C(1,2))/pivot
         B(2) = (C(2,3)*C(1,1) - C(1,3)*C(2,1))/pivot
         B(1) = W(1)*B(1)
         B(2) = W(2)*B(2)
         RETURN
      END IF

      DO m=1, N-1
         ipv   = m
         pivot = ABS(C(m,m))
         DO i=m+1, N
            IF (ABS(C(i,m)) .GT. pivot) THEN
               ipv   = i
               pivot = ABS(C(i,m))
            END IF
         END DO
         IF (ABS(pivot) .LT. EPSILON(pivot)) THEN
!     Singular matrix            
            B  = 0D0
            GE = .FALSE.
            RETURN
         END IF
         IF (ipv .NE. m) THEN
            DO j=m, N+1
               saveEl   = C(m,j)
               C(m,j)   = C(ipv,j)
               C(ipv,j) = saveEl
            END DO
         END IF

         DO i=m+1, N
            saveEl = C(i,m)/C(m,m)
            C(i,m) = 0D0
            DO j=m+1, N+1
               C(i,j) = C(i,j) - saveEl*C(m,j)
            END DO
         END DO
      END DO

      DO j=N, 1, -1
         DO i=j+1, N
            C(j,N+1) = C(j,N+1) - C(j,i)*C(i,N+1)
         END DO
         C(j,N+1) = C(j,N+1)/C(j,j)
      END DO
      DO i=1, N
         B(i) = W(i)*C(i,N+1)
      END DO

      DEALLOCATE(C)

      RETURN
      END FUNCTION GE

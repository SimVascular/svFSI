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
!     Product of a sparse matrix and a vector. The matrix might be
!     vector in neither, one or both dimensions.
!--------------------------------------------------------------------

      SUBROUTINE FSILS_SPARMULVV(lhs, rowPtr, colPtr, dof, K, U, KU)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: K(dof*dof,lhs%nnz), U(dof,lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: KU(dof,lhs%nNo)

      INTEGER nNo, i, j, l, col, s, e

      nNo = lhs%nNo

      KU = 0D0
      SELECT CASE (dof)
      CASE (1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               KU(1,i) = KU(1,i) + K(1,j)*U(1,colPtr(j))
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)
               KU(2,i) = KU(2,i) + K(3,j)*U(1,col) + K(4,j)*U(2,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)    &
     &                           + K(3,j)*U(3,col)
               KU(2,i) = KU(2,i) + K(4,j)*U(1,col) + K(5,j)*U(2,col)    &
     &                           + K(6,j)*U(3,col)
               KU(3,i) = KU(3,i) + K(7,j)*U(1,col) + K(8,j)*U(2,col)    &
     &                           + K(9,j)*U(3,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1 ,j)*U(1,col) + K(2 ,j)*U(2,col)  &
     &                           + K(3 ,j)*U(3,col) + K(4 ,j)*U(4,col)
               KU(2,i) = KU(2,i) + K(5 ,j)*U(1,col) + K(6 ,j)*U(2,col)  &
     &                           + K(7 ,j)*U(3,col) + K(8 ,j)*U(4,col)
               KU(3,i) = KU(3,i) + K(9 ,j)*U(1,col) + K(10,j)*U(2,col)  &
     &                           + K(11,j)*U(3,col) + K(12,j)*U(4,col)
               KU(4,i) = KU(4,i) + K(13,j)*U(1,col) + K(14,j)*U(2,col)  &
     &                           + K(15,j)*U(3,col) + K(16,j)*U(4,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, l, e, s, col) 
!$OMP&   SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               DO l=1, dof
                  e = l*dof
                  s = e - dof + 1
                  KU(l,i) = KU(l,i) + SUM(K(s:e,j)*U(:,col))
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      CALL FSILS_COMMUV(lhs, dof, KU)

      RETURN
      END SUBROUTINE FSILS_SPARMULVV

!====================================================================

      SUBROUTINE FSILS_SPARMULVS(lhs, rowPtr, colPtr, dof, K, U, KU)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: K(dof,lhs%nnz), U(dof,lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: KU(lhs%nNo)
      
      INTEGER nNo, i, j, col
      
      nNo = lhs%nNo

      KU = 0D0
      SELECT CASE (dof)
      CASE (1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               KU(i) = KU(i) + K(1,j)*U(1,colPtr(j))
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)        &
     &                       + K(3,j)*U(3,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)        &
     &                       + K(3,j)*U(3,col) + K(4,j)*U(4,col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               KU(i) = KU(i) + SUM(K(:,j)*U(:,colPtr(j)))
            END DO
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      CALL FSILS_COMMUS(lhs, KU)

      RETURN
      END SUBROUTINE FSILS_SPARMULVS

!====================================================================

      SUBROUTINE FSILS_SPARMULSV(lhs, rowPtr, colPtr, dof, K, U, KU)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(IN) :: K(dof,lhs%nnz), U(lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: KU(dof,lhs%nNo)

      INTEGER nNo, i, j, col

      nNo = lhs%nNo
      
      KU = 0D0
      SELECT CASE (dof)
      CASE (1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               KU(1,i) = KU(1,i) + K(1,j)*U(colPtr(j))
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1,j)*U(col)
               KU(2,i) = KU(2,i) + K(2,j)*U(col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1,j)*U(col)
               KU(2,i) = KU(2,i) + K(2,j)*U(col)
               KU(3,i) = KU(3,i) + K(3,j)*U(col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, col) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               col = colPtr(j)
               KU(1,i) = KU(1,i) + K(1,j)*U(col)
               KU(2,i) = KU(2,i) + K(2,j)*U(col)
               KU(3,i) = KU(3,i) + K(3,j)*U(col)
               KU(4,i) = KU(4,i) + K(4,j)*U(col)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
         DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
               KU(:,i) = KU(:,i) + K(:,j)*U(colPtr(j))
            END DO
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      CALL FSILS_COMMUV(lhs, dof, KU)

      RETURN
      END SUBROUTINE FSILS_SPARMULSV

!====================================================================

      SUBROUTINE FSILS_SPARMULSS(lhs, rowPtr, colPtr, K, U, KU)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: rowPtr(2,lhs%nNo), colPtr(lhs%nnz)
      REAL(KIND=8), INTENT(IN) :: K(lhs%nnz), U(lhs%nNo)
      REAL(KIND=8), INTENT(OUT) :: KU(lhs%nNo)

      INTEGER nNo, i, j

      nNo = lhs%nNo

      KU = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j) SCHEDULE(GUIDED)
      DO i=1, nNo
         DO j=rowPtr(1,i), rowPtr(2,i)
            KU(i) = KU(i) + K(j)*U(colPtr(j))
         END DO
      END DO
!$OMP END PARALLEL DO         

      CALL FSILS_COMMUS(lhs, KU)

      RETURN
      END SUBROUTINE FSILS_SPARMULSS


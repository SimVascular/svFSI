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
!     A bunch of operation that benefits from OMP hyperthreading
!--------------------------------------------------------------------
      
      SUBROUTINE OMPSUMS (nNo, r, U, V)
      
      INCLUDE "FSILS_STD.h"

      INTEGER, INTENT(IN) :: nNo
      REAL(KIND=8), INTENT(IN) :: r, V(nNo)
      REAL(KIND=8), INTENT(INOUT) :: U(nNo)
     
      INTEGER i

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
      DO i=1, nNo
         U(i) = U(i) + r*V(i)
      END DO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE OMPSUMS

!====================================================================

      SUBROUTINE OMPSUMV (dof, nNo, r, U, V)
      
      INCLUDE "FSILS_STD.h"

      INTEGER, INTENT(IN) :: dof, nNo
      REAL(KIND=8), INTENT(IN) :: r, V(dof,nNo)
      REAL(KIND=8), INTENT(INOUT) :: U(dof,nNo)
     
      INTEGER i

      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = U(1,i) + r*V(1,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = U(1,i) + r*V(1,i)
            U(2,i) = U(2,i) + r*V(2,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = U(1,i) + r*V(1,i)
            U(2,i) = U(2,i) + r*V(2,i)
            U(3,i) = U(3,i) + r*V(3,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = U(1,i) + r*V(1,i)
            U(2,i) = U(2,i) + r*V(2,i)
            U(3,i) = U(3,i) + r*V(3,i)
            U(4,i) = U(4,i) + r*V(4,i)
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(:,i) = U(:,i) + r*V(:,i)
         END DO
!$OMP END PARALLEL DO
      END SELECT

      RETURN
      END SUBROUTINE OMPSUMV

!====================================================================
      
      SUBROUTINE OMPMULS (nNo, r, U)
      
      INCLUDE "FSILS_STD.h"

      INTEGER, INTENT(IN) :: nNo
      REAL(KIND=8), INTENT(IN) :: r
      REAL(KIND=8), INTENT(INOUT) :: U(nNo)
     
      INTEGER i

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
      DO i=1, nNo
         U(i) = r*U(i)
      END DO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE OMPMULS

!====================================================================

      SUBROUTINE OMPMULV (dof, nNo, r, U)
      
      INCLUDE "FSILS_STD.h"

      INTEGER, INTENT(IN) :: dof, nNo
      REAL(KIND=8), INTENT(IN) :: r
      REAL(KIND=8), INTENT(INOUT) :: U(dof,nNo)
     
      INTEGER i

      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = r*U(1,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = r*U(1,i)
            U(2,i) = r*U(2,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = r*U(1,i)
            U(2,i) = r*U(2,i)
            U(3,i) = r*U(3,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(1,i) = r*U(1,i)
            U(2,i) = r*U(2,i)
            U(3,i) = r*U(3,i)
            U(4,i) = r*U(4,i)
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
         DO i=1, nNo
            U(:,i) = r*U(:,i)
         END DO
!$OMP END PARALLEL DO
      END SELECT

      RETURN
      END SUBROUTINE OMPMULV

!====================================================================


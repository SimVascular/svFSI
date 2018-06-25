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
!     For calculating the norm of a scaler or vector based vector.
!--------------------------------------------------------------------
      
!     Only the part of U which is owned by this processor is included in
!     norm calculation, i.e. U(:,1:cS(tF)%ptr+cS(tF)%n-1)
!     In order to have the correct answer it is needed that COMMU has
!     been done before calling this function (or the ansesters of U
!     are passed through COMMU)
      FUNCTION FSILS_NORMV(dof, nNo, commu, U)
 
      INCLUDE "FSILS_STD.h"
       
      INTEGER, INTENT(IN) :: dof, nNo
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: U(dof,nNo)
      
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_NORMV
            
      FSILS_NORMV = 0D0
      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMV)
         DO i=1, nNo
            FSILS_NORMV = FSILS_NORMV + U(1,i)*U(1,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMV)
         DO i=1, nNo
            FSILS_NORMV = FSILS_NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMV)
         DO i=1, nNo
            FSILS_NORMV = FSILS_NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i) + &
     &         U(3,i)*U(3,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMV)
         DO i=1, nNo
            FSILS_NORMV = FSILS_NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i) + &
     &         U(3,i)*U(3,i) + U(4,i)*U(4,i)
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMV)
         DO i=1, nNo
            FSILS_NORMV = FSILS_NORMV + SUM(U(:,i)*U(:,i))
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      IF (commu%nTasks .NE. 1) THEN
         CALL MPI_ALLREDUCE(FSILS_NORMV, tmp, 1, mpreal, MPI_SUM,       &
     &      commu%comm, ierr)
         FSILS_NORMV = tmp
      END IF
      FSILS_NORMV = SQRT(FSILS_NORMV)

      RETURN
      END FUNCTION FSILS_NORMV

!====================================================================
      
      FUNCTION FSILS_NORMS(nNo, commu, U)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: nNo
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: U(nNo)
      
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_NORMS
      
      FSILS_NORMS = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_NORMS)
      DO i=1, nNo
         FSILS_NORMS = FSILS_NORMS + U(i)*U(i)
      END DO
!$OMP END PARALLEL DO         

      IF (commu%nTasks .NE. 1) THEN
         CALL MPI_ALLREDUCE(FSILS_NORMS, tmp, 1, mpreal, MPI_SUM,       &
     &      commu%comm, ierr)
         FSILS_NORMS = tmp
      END IF
      FSILS_NORMS = SQRT(FSILS_NORMS)

      RETURN
      END FUNCTION FSILS_NORMS



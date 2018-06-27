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
!     The dot product between two scaler-vector container vectors are
!     calculated here.
!--------------------------------------------------------------------

!     Only the part of U and V which are owned by this processor is 
!     included in dot product calculation
!     In order to have the correct answer it is needed that COMMU has
!     been done before calling this function (or the ansesters of U and
!     V are passed through COMMU)
      FUNCTION FSILS_DOTV(dof, nNo, commu, U, V)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: dof, nNo
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: V(dof,nNo), U(dof,nNo)
      
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_DOTV

      FSILS_DOTV = 0D0
      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) +   &
     &         U(3,i)*V(3,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) +   &
     &         U(3,i)*V(3,i) + U(4,i)*V(4,i)
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + SUM(U(:,i)*V(:,i))
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      IF (commu%nTasks .EQ. 1) RETURN
      CALL MPI_ALLREDUCE(FSILS_DOTV, tmp, 1, mpreal, MPI_SUM,           &
     &   commu%comm, ierr)

      FSILS_DOTV = tmp

      RETURN
      END FUNCTION FSILS_DOTV

!====================================================================
      
      FUNCTION FSILS_DOTS(nNo, commu, U, V)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: nNo
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: V(nNo), U(nNo)
 
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_DOTS

      FSILS_DOTS = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTS)
      DO i=1, nNo
         FSILS_DOTS = FSILS_DOTS + U(i)*V(i)
      END DO
!$OMP END PARALLEL DO         

      IF (commu%nTasks .EQ. 1) RETURN
      CALL MPI_ALLREDUCE(FSILS_DOTS, tmp, 1, mpreal, MPI_SUM,           &
     &   commu%comm, ierr)

      FSILS_DOTS = tmp

      RETURN
      END FUNCTION FSILS_DOTS

!####################################################################
      
      FUNCTION FSILS_NCDOTV(dof, nNo, U, V) RESULT(FSILS_DOTV)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: dof, nNo
      REAL(KIND=8), INTENT(IN) :: V(dof,nNo), U(dof,nNo)
      
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_DOTV

      FSILS_DOTV = 0D0
      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) +   &
     &         U(3,i)*V(3,i)
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) +   &
     &         U(3,i)*V(3,i) + U(4,i)*V(4,i)
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTV)
         DO i=1, nNo
            FSILS_DOTV = FSILS_DOTV + SUM(U(:,i)*V(:,i))
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      RETURN
      END FUNCTION FSILS_NCDOTV

!====================================================================
      
      FUNCTION FSILS_NCDOTS(nNo, U, V) RESULT(FSILS_DOTS)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: nNo
      REAL(KIND=8), INTENT(IN) :: V(nNo), U(nNo)
 
      INTEGER i, ierr
      REAL(KIND=8) tmp, FSILS_DOTS

      FSILS_DOTS = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:FSILS_DOTS)
      DO i=1, nNo
         FSILS_DOTS = FSILS_DOTS + U(i)*V(i)
      END DO
!$OMP END PARALLEL DO         

      RETURN
      END FUNCTION FSILS_NCDOTS

!####################################################################
 
      SUBROUTINE FSILS_MULTDOTV(dof, nNo, tnNo, nV, commu, U, V, res)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: dof, nNo, tnNo, nV
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: U(dof,tnNo,nV), V(dof,tnNo)
      REAL(KIND=8), INTENT(OUT) :: res(nV)
      
      INTEGER i, j, ierr, iV
      REAL(KIND=8) tmp(nV)

      res = 0D0
      SELECT CASE(dof)
      CASE(1)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
         DO iV=1, nV
            DO i=1, nNo
               res(iV) = res(iV) + U(1,i,iV)*V(1,i)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(2)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
         DO iV=1, nV
            DO i=1, nNo
               res(iV) = res(iV) + U(1,i,iV)*V(1,i) + U(2,i,iV)*V(2,i)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(3)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
         DO iV=1, nV
            DO i=1, nNo
               res(iV) = res(iV) + U(1,i,iV)*V(1,i) + U(2,i,iV)*V(2,i) +&
     &            U(3,i,iV)*V(3,i)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE(4)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
         DO iV=1, nV
            DO i=1, nNo
               res(iV) = res(iV) + U(1,i,iV)*V(1,i) + U(2,i,iV)*V(2,i) +&
     &            U(3,i,iV)*V(3,i) + U(4,i,iV)*V(4,i)
            END DO
         END DO
!$OMP END PARALLEL DO         
      CASE DEFAULT 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
         DO iV=1, nV
            DO i=1, nNo
               DO j=1, dof
                  res(iV) = res(iV) + U(j,i,iV)*V(j,i)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO         
      END SELECT

      IF (commu%nTasks .EQ. 1) RETURN
      CALL MPI_ALLREDUCE(res, tmp, nV, mpreal, MPI_SUM, commu%comm,ierr)
      
      res = tmp

      RETURN
      END SUBROUTINE FSILS_MULTDOTV

!====================================================================
  
      SUBROUTINE FSILS_MULTDOTS(nNo, tnNo, nV, commu, U, V, res)
 
      INCLUDE "FSILS_STD.h"
      
      INTEGER, INTENT(IN) :: nNo, tnNo, nV
      TYPE(FSILS_commuType), INTENT(IN) :: commu
      REAL(KIND=8), INTENT(IN) :: U(tnNo,nV), V(tnNo)
      REAL(KIND=8), INTENT(OUT) :: res(nV)
      
      INTEGER i, ierr, iV
      REAL(KIND=8) tmp(nV)

      res = 0D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iV) SCHEDULE(GUIDED)
!$OMP&   REDUCTION(+:res)
      DO iV=1, nV
         DO i=1, nNo
            res(iV) = res(iV) + U(i,iV)*V(i)
         END DO
      END DO
!$OMP END PARALLEL DO         

      IF (commu%nTasks .EQ. 1) RETURN
      CALL MPI_ALLREDUCE(res, tmp, nV, mpreal, MPI_SUM, commu%comm,ierr)
      
      res = tmp

      RETURN
      END SUBROUTINE FSILS_MULTDOTS

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
!     To syncronize the data on the boundaries between the processors.
!--------------------------------------------------------------------
      
!     This a both way communication with three main part:
!     1 - rTmp {in master} = R          {from slave}
!     2 - R    {in master} = R + rTmp   {both from master}
!     3 - rTmp {in master} = R          {from master}
!     4 - R    {in slave}  = rTmp       {from master}

      SUBROUTINE FSILS_COMMUV(lhs, dof, R)
 
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER, INTENT(IN) :: dof
      REAL(KIND=8), INTENT(INOUT) :: R(dof,lhs%nNo)
 
      INTEGER nReq, i, j, k, ierr, stat(mpsts*lhs%nReq)
      REAL(KIND=8), ALLOCATABLE :: sB(:,:,:), rB(:,:,:), sReq(:),       &
     &   rReq(:)

      IF (lhs%commu%nTasks .EQ. 1) RETURN
      
      nReq = lhs%nReq
      i = MAXVAL(lhs%cS%n)
      ALLOCATE(sB(dof,i,nReq), rB(dof,i,nReq), rReq(nReq), sReq(nReq))
      
      DO i=1, nReq
         DO j=1, lhs%cS(i)%n
            k = lhs%cS(i)%ptr(j)
            sB(:,j,i) = R(:,k)
         END DO
      END DO

      DO i=1, nReq
         CALL MPI_IRECV(rB(:,:,i), lhs%cS(i)%n*dof, mpreal,             &
     &      lhs%cS(i)%iP-1, 1, lhs%commu%comm, rReq(i), ierr)
         CALL MPI_ISEND(sB(:,:,i), lhs%cS(i)%n*dof, mpreal,             &
     &      lhs%cS(i)%iP-1, 1, lhs%commu%comm, sReq(i), ierr)
      END DO

      DO i=1, nReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO

      DO i=1, nReq
         DO j=1, lhs%cS(i)%n
            k = lhs%cS(i)%ptr(j)
            R(:,k) = R(:,k) + rB(:,j,i)
         END DO
      END DO
      
      DO i=1, nReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO

      RETURN
      END SUBROUTINE FSILS_COMMUV

!====================================================================
      
      SUBROUTINE FSILS_COMMUS(lhs, R)
 
      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      REAL(KIND=8), INTENT(INOUT) :: R(lhs%nNo)
      
      INTEGER nReq, i, j, k, ierr, stat(mpsts*lhs%nReq)
      REAL(KIND=8), ALLOCATABLE :: sB(:,:), rB(:,:), sReq(:), rReq(:)

      IF (lhs%commu%nTasks .EQ. 1) RETURN
 
      nReq = lhs%nReq
      i = MAXVAL(lhs%cS%n)
      ALLOCATE(sB(i,nReq), rB(i,nReq), rReq(nReq), sReq(nReq))
      
      DO i=1, nReq
         DO j=1, lhs%cS(i)%n
            k = lhs%cS(i)%ptr(j)
            sB(j,i) = R(k)
         END DO
      END DO

      DO i=1, nReq
         CALL MPI_IRECV(rB(:,i), lhs%cS(i)%n, mpreal, lhs%cS(i)%iP-1, 1,&
     &      lhs%commu%comm, rReq(i), ierr)
         CALL MPI_ISEND(sB(:,i), lhs%cS(i)%n, mpreal, lhs%cS(i)%iP-1, 1,&
     &      lhs%commu%comm, sReq(i), ierr)
      END DO

      DO i=1, nReq
         CALL MPI_WAIT(rReq(i), stat, ierr)
      END DO

      DO i=1, nReq
         DO j=1, lhs%cS(i)%n
            k = lhs%cS(i)%ptr(j)
            R(k) = R(k) + rB(j,i)
         END DO
      END DO
      
      DO i=1, nReq
         CALL MPI_WAIT(sReq(i), stat, ierr)
      END DO

      RETURN
      END SUBROUTINE FSILS_COMMUS
      

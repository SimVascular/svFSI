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
!     To broadcast a variable to all processors      
!--------------------------------------------------------------------

      SUBROUTINE FSILS_BCAST(u, commu)

      INCLUDE "FSILS_STD.h"

      REAL(KIND=8), INTENT(INOUT) :: u
      TYPE(FSILS_commuType), INTENT(IN) :: commu

      INTEGER ierr
      REAL(KIND=8) uG

      IF (commu%nTasks .GT. 1) THEN
         CALL MPI_ALLREDUCE(u, uG, 1, mpreal, MPI_SUM, commu%comm, ierr)
         u = uG
      END IF

      RETURN
      END SUBROUTINE FSILS_BCAST
      
!====================================================================

      SUBROUTINE FSILS_BCASTV(n, u, commu)

      INCLUDE "FSILS_STD.h"

      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), INTENT(INOUT) :: u(n)
      TYPE(FSILS_commuType), INTENT(IN) :: commu

      INTEGER ierr
      REAL(KIND=8), ALLOCATABLE :: uG(:)

      IF (commu%nTasks .GT. 1) THEN
         ALLOCATE(uG(n))
         CALL MPI_ALLREDUCE(u, uG, n, mpreal, MPI_SUM, commu%comm, ierr)
         u = uG
      END IF

      RETURN
      END SUBROUTINE FSILS_BCASTV

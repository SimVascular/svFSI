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
!     Communication structure is created here.
!--------------------------------------------------------------------

      SUBROUTINE FSILS_COMMU_CREATE(commu, commi)
      
      INCLUDE "FSILS_STD.h"
   
      TYPE(FSILS_commuType), INTENT(INOUT) :: commu
      INTEGER, INTENT(IN) :: commi

      INTEGER ierr

      IF (commu%foC) THEN
         PRINT *, "FSILS: COMMU is not free, you may use",              &
     &      " FSILS_COMMU_FREE to free it"
         STOP "FSILS: FATAL ERROR"
      END IF

!     Some of these parameters are set for sequential version
      commu%foC    = .TRUE.
      commu%comm   = commi
      commu%nTasks = 1
      commu%task   = 0
      commu%master = 0

      CALL MPI_COMM_RANK(commi, commu%task, ierr)
      CALL MPI_COMM_SIZE(commi, commu%nTasks, ierr)
      CALL MPI_ALLREDUCE(commu%task, commu%master, 1, mpint, MPI_MIN,   &
     &   commi, ierr)
      
      IF (commu%master .NE. 0) THEN
         STOP "Master ID is not zero - might cause problems"
      END IF

      commu%masF = .FALSE.
      commu%tF   = commu%task + 1
      IF (commu%task .EQ. commu%master) THEN
         commu%masF = .TRUE.
      END IF

      RETURN
      END SUBROUTINE FSILS_COMMU_CREATE

!====================================================================

      SUBROUTINE FSILS_COMMU_FREE(commu)
      
      INCLUDE "FSILS_STD.h"
   
      TYPE(FSILS_commuType), INTENT(INOUT) :: commu

      IF (.NOT.commu%foC) STOP 'COMMU is not created yet to be freed'
      commu%foC  = .FALSE.

      RETURN
      END SUBROUTINE FSILS_COMMU_FREE

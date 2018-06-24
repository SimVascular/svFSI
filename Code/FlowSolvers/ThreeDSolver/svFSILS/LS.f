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
!     This is to construct the structure for LS.
!--------------------------------------------------------------------

      SUBROUTINE FSILS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr,   &
     &   dimKry, relTolIn, absTolIn, maxItrIn)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lsType), INTENT(INOUT) :: ls
      INTEGER, INTENT(IN) :: LS_type
      REAL(KIND=8), INTENT(IN), OPTIONAL :: relTol, absTol, relTolIn(2),&
     &   absTolIn(2)
      INTEGER, INTENT(IN), OPTIONAL :: maxItr, dimKry, maxItrIn(2)

      IF (ls%foC) THEN
         PRINT *, "FSILS: LS is not free You may use FSILS_LS_FREE to", &
     &      " free this structure"
         STOP "FSILS: FATAL ERROR"
      END IF

      ls%foC     = .TRUE.
      ls%LS_type = LS_type

      SELECT CASE (LS_type)
         CASE (LS_TYPE_NS)
            ls%RI%relTol = 4D-1
            ls%GM%relTol = 1D-2
            ls%CG%relTol = 2D-1
            ls%RI%mItr = 10
            ls%GM%mItr = 2
            ls%CG%mItr = 500
            ls%GM%sD   = 100
         CASE (LS_TYPE_GMRES)
            ls%RI%relTol = 1D-1
            ls%RI%mItr   = 4
            ls%RI%sD     = 250
         CASE (LS_TYPE_CG)
            ls%RI%reltol = 1D-2
            ls%RI%mItr   = 1000
         CASE (LS_TYPE_BICGS)
            ls%RI%reltol = 1D-2
            ls%RI%mItr   = 500
         CASE DEFAULT
            PRINT *, 'FSILS: LS_TYPE is not defined'
            STOP "FSILS: FATAL ERROR"
      END SELECT
      ls%RI%absTol = 1D-10
      ls%GM%absTol = 1D-10
      ls%CG%absTol = 1D-10

      IF (PRESENT(relTol)) ls%RI%relTol = relTol
      IF (PRESENT(absTol)) ls%RI%absTol = absTol
      IF (PRESENT(maxItr)) ls%RI%mItr   = maxItr

      IF (PRESENT(dimKry)) THEN
         ls%RI%sD = dimKry
         ls%GM%sD = dimKry
      END IF
      IF (PRESENT(relTolIn)) THEN
         ls%GM%relTol = relTolIn(1)
         ls%CG%relTol = relTolIn(2)
      END IF
      IF (PRESENT(absTolIn)) THEN
         ls%GM%absTol = absTolIn(1)
         ls%CG%absTol = absTolIn(2)
      END IF
      IF (PRESENT(maxItrIn)) THEN
         ls%GM%mItr = maxItrIn(1)
         ls%CG%mItr = maxItrIn(2)
      END IF

      RETURN
      END SUBROUTINE FSILS_LS_CREATE

!====================================================================

      SUBROUTINE FSILS_LS_FREE (ls)

      INCLUDE "FSILS_STD.h"

      TYPE(FSILS_lsType), INTENT(INOUT) :: ls

      IF (.NOT.ls%foC) THEN
         PRINT *, 'FSILS: LS is not created yet to be freed'
         STOP "FSILS: FATAL ERROR"
      END IF

      ls%foC  = .FALSE.

      RETURN
      END SUBROUTINE FSILS_LS_FREE


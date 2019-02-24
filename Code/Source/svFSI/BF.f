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
!     This routine computes body force on a mesh.
!
!--------------------------------------------------------------------

      SUBROUTINE GETBF()
      USE COMMOD
      IMPLICIT NONE

      LOGICAL flag
      INTEGER iM

      flag = .FALSE.
      DO iM=1, nMsh
         IF (ALLOCATED(msh(iM)%bf)) THEN
            flag = .TRUE.
            EXIT
         END IF
      END DO
      IF (flag) THEN
         IF (ALLOCATED(Bfg)) DEALLOCATE(Bfg)
         ALLOCATE(Bfg(nsd,tnNo))
         Bfg = 0D0
         DO iM=1, nMsh
            IF (ALLOCATED(msh(iM)%bf)) CALL GETBFL(msh(iM), Bfg)
         END DO
      END IF

      RETURN
      END SUBROUTINE GETBF
!--------------------------------------------------------------------
      SUBROUTINE GETBFL(lM, lBf)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(INOUT) :: lBf(nsd,tnNo)

      INTEGER :: a, Ac
      REAL(KIND=8), ALLOCATABLE :: tmpY(:,:), tmpA(:,:)

      ALLOCATE(tmpY(lM%bf%dof,lM%nNo), tmpA(lM%bf%dof,lM%nNo))
      CALL IGBC(lM%bf, tmpY, tmpA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         lBf(:,Ac) = tmpY(:,a)
      END DO

      DEALLOCATE(tmpY, tmpA)

      RETURN
      END SUBROUTINE GETBFL
!####################################################################

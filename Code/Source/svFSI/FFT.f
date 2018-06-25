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
!     Here Fast Fourier Transfer (FFT) and IFFT are calculated as well
!     as values based on the linear interpolation.
!      
!--------------------------------------------------------------------

      SUBROUTINE FFT(fid, np, gt)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: fid, np
      TYPE(fcType), INTENT(INOUT) :: gt

      INTEGER i, n
      REAL(KIND=8) tmp, kn, ko, s

      REAL(KIND=8), ALLOCATABLE :: t(:), q(:)

      ALLOCATE (t(np), q(np))
      DO i=1, np
         READ (fid,*) t(i), q(i)
      END DO

      gt%ti = t(1)
      gt%T  = t(np) - t(1)
      gt%qi = q(1)
      gt%qs = (q(np) - q(1))/gt%T

      DO i=1, np
         t(i) = t(i) - gt%ti
         q(i) = q(i) - gt%qi - gt%qs*t(i)
      END DO

      DO n=1, gt%n
         tmp = REAL(n-1,8)
         gt%r(n) = 0D0
         gt%i(n) = 0D0
         DO i=1, np-1
            ko = 2D0*pi*tmp*t(i)/gt%T
            kn = 2D0*pi*tmp*t(i+1)/gt%T
            s  = (q(i+1) - q(i))/(t(i+1) - t(i))

            IF (n .EQ. 1) THEN
               gt%r(n) = gt%r(n) + 5D-1*(t(i+1)-t(i))*(q(i+1)+q(i))
            ELSE
               gt%r(n) = gt%r(n) + s*(COS(kn) - COS(ko))
               gt%i(n) = gt%i(n) - s*(SIN(kn) - SIN(ko))
            END IF
         END DO

         IF (n .EQ. 1) THEN
            gt%r(n) = gt%r(n)/gt%T
         ELSE
            gt%r(n) = 5D-1*gt%r(n)*gt%T/(pi*pi*tmp*tmp)
            gt%i(n) = 5D-1*gt%i(n)*gt%T/(pi*pi*tmp*tmp)
         END IF
      END DO

      RETURN
      END SUBROUTINE FFT
!-------------------------------------------------------------------- 
!     This is to calculate flow rate and flow acceleration (IFFT)
      PURE SUBROUTINE IFFT(gt, Y, dY)
      
      USE COMMOD
      
      IMPLICIT NONE

      TYPE(fcType), INTENT(IN) :: gt
      REAL(KIND=8), INTENT(OUT) :: Y, dY

      INTEGER i
      REAL(KIND=8) t, tmp, K, kd

      t    = DMOD(time - gt%ti, gt%T)
      tmp  = 2D0*pi/gt%T
      Y    = gt%qi + t*gt%qs
      dY   = gt%qs
      DO i=1, gt%n
         kd = tmp*REAL(i-1,8)
         K  = t*kd
         Y  =  Y +  gt%r(i)*COS(K) - gt%i(i)*SIN(K)
         dY = dY - (gt%r(i)*SIN(K) + gt%i(i)*COS(K))*kd
      END DO

      RETURN 
      END SUBROUTINE IFFT

!--------------------------------------------------------------------
!     This routine is for calculating values by the inverse of general 
!     BC
      PURE SUBROUTINE IGBC(gm, Y, dY)
      
      USE COMMOD
      
      IMPLICIT NONE

      TYPE(MBType), INTENT(IN) :: gm
      REAL(KIND=8), INTENT(OUT) :: Y(gm%dof,SIZE(gm%d,2)),
     2   dY(gm%dof,SIZE(gm%d,2))

      INTEGER a, i
      REAL(KIND=8) t, tmp, delT

      t = DMOD(time,gm%period)
      DO i=1, gm%nTP - 1
         IF (gm%t(i+1) .GE. t) THEN
            Y  = 0D0
            dY = 0D0
            EXIT
         END IF
      END DO
      delT = gm%t(i+1) - gm%t(i)
      tmp  = (t - gm%t(i))/delT
      DO a=1, SIZE(gm%d,2)
         Y (:,a) = tmp*gm%d(:,a,i+1) + gm%d(:,a,i)*(1D0-tmp)
         dY(:,a) =    (gm%d(:,a,i+1) - gm%d(:,a,i))/delT
      END DO

      RETURN
      END SUBROUTINE IGBC

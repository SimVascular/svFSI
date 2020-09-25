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
      INTEGER(KIND=IKIND), INTENT(IN) :: fid, np
      TYPE(fcType), INTENT(INOUT) :: gt

      INTEGER(KIND=IKIND) i, j, n
      REAL(KIND=RKIND) tmp, kn, ko, s

      REAL(KIND=RKIND), ALLOCATABLE :: t(:), q(:,:)

      ALLOCATE (t(np), q(gt%d,np))
      DO i=1, np
         READ (fid,*) t(i), (q(j,i), j=1, gt%d)
      END DO

      gt%ti    = t(1)
      gt%T     = t(np) - t(1)
      gt%qi(:) = q(:,1)
      gt%qs(:) = (q(:,np) - q(:,1))/gt%T

      DO i=1, np
         t(i) = t(i) - gt%ti
         DO j=1, gt%d
            q(j,i) = q(j,i) - gt%qi(j) - gt%qs(j)*t(i)
         END DO
      END DO

      DO n=1, gt%n
         tmp = REAL(n-1, KIND=RKIND)
         gt%r(:,n) = 0._RKIND
         gt%i(:,n) = 0._RKIND
         DO i=1, np-1
            ko = 2._RKIND*pi*tmp*t(i)/gt%T
            kn = 2._RKIND*pi*tmp*t(i+1)/gt%T
            DO j=1, gt%d
               s = (q(j,i+1) - q(j,i))/(t(i+1) - t(i))

               IF (n .EQ. 1) THEN
                  gt%r(j,n) = gt%r(j,n) + 0.5_RKIND*(t(i+1) - t(i))
     2               * (q(j,i+1) + q(j,i))
               ELSE
                  gt%r(j,n) = gt%r(j,n) + s*(COS(kn) - COS(ko))
                  gt%i(j,n) = gt%i(j,n) - s*(SIN(kn) - SIN(ko))
               END IF
            END DO
         END DO

         IF (n .EQ. 1) THEN
            gt%r(:,n) = gt%r(:,n)/gt%T
         ELSE
            gt%r(:,n) = 0.5_RKIND*gt%r(:,n)*gt%T/(pi*pi*tmp*tmp)
            gt%i(:,n) = 0.5_RKIND*gt%i(:,n)*gt%T/(pi*pi*tmp*tmp)
         END IF
      END DO

      DEALLOCATE(t, q)

      RETURN
      END SUBROUTINE FFT
!--------------------------------------------------------------------
!     This is to calculate flow rate and flow acceleration (IFFT)
      PURE SUBROUTINE IFFT(gt, Y, dY)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fcType), INTENT(IN) :: gt
      REAL(KIND=RKIND), INTENT(OUT) :: Y(gt%d), dY(gt%d)

      INTEGER(KIND=IKIND) i, j
      REAL(KIND=RKIND) t, tmp, K, dk

      IF (gt%lrmp) THEN
         t = time - gt%ti
         IF (t .LE. 0._RKIND) THEN
            t = MAX(t, 0._RKIND)
         ELSE
            t = MIN(t, gt%T)
         END IF
         DO i=1, gt%d
            Y(i)  = gt%qi(i) + t*gt%qs(i)
            dY(i) = gt%qs(i)
         END DO
      ELSE
         t   = MOD(time - gt%ti, gt%T)
         tmp = 2._RKIND*pi/gt%T
         DO j=1, gt%d
            Y(j)  = gt%qi(j) + t*gt%qs(j)
            dY(j) = gt%qs(j)
         END DO

         DO i=1, gt%n
            dk = tmp*REAL(i-1, KIND=RKIND)
            K  = t*dk
            DO j=1, gt%d
               Y(j)  =  Y(j) +  gt%r(j,i)*COS(K) - gt%i(j,i)*SIN(K)
               dY(j) = dY(j) - (gt%r(j,i)*SIN(K) + gt%i(j,i)*COS(K))*dk
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE IFFT
!####################################################################
!     This routine is for calculating values by the inverse of general
!     BC
      PURE SUBROUTINE IGBC(gm, Y, dY)
      USE COMMOD
      IMPLICIT NONE
      TYPE(MBType), INTENT(IN) :: gm
      REAL(KIND=RKIND), INTENT(OUT) :: Y(gm%dof,SIZE(gm%d,2)),
     2   dY(gm%dof,SIZE(gm%d,2))

      INTEGER(KIND=IKIND) a, i
      REAL(KIND=RKIND) t, tmp, delT

      t = MOD(time, gm%period)
      DO i=1, gm%nTP - 1
         IF (gm%t(i+1) .GE. t) THEN
            Y  = 0._RKIND
            dY = 0._RKIND
            EXIT
         END IF
      END DO
      delT = gm%t(i+1) - gm%t(i)
      tmp  = (t - gm%t(i))/delT
      DO a=1, SIZE(gm%d,2)
         Y (:,a) = tmp*gm%d(:,a,i+1) + gm%d(:,a,i)*(1._RKIND-tmp)
         dY(:,a) = (gm%d(:,a,i+1) - gm%d(:,a,i))/delT
      END DO

      RETURN
      END SUBROUTINE IGBC
!####################################################################

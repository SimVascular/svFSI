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
!-----------------------------------------------------------------------
!
!     These subroutines process data for non-electrophysiology type
!     (decoupled), uniformly activated, excitation-contraction coupling.
!
!
!-----------------------------------------------------------------------

!     Initialize decoupled excitation-contraction coupling
      SUBROUTINE EC_DCPLD_INIT(ec)
      USE ECMOD
      USE COMMOD, ONLY : eccModelType, IKIND, RKIND, std
      IMPLICIT NONE
      TYPE(eccModelType), INTENT(INOUT) :: ec

      INTEGER(KIND=IKIND) slen

      ec%Ya = 0._RKIND

!     Overwrite parameters with user-provided input file
      slen = LEN(TRIM(ec%fpar_in))
      IF (slen .GT. 0) THEN
         std = " Reading decoupled excitation-contraction parameters"//
     2      " from input file"
         CALL EC_READPARFF(ec%fpar_in)
      END IF

      RETURN
      END SUBROUTINE EC_DCPLD_INIT
!####################################################################
!     Compute the activation parameter (active stress/fiber shortening)
!     for decoupled excitation-contraction. Dual time-stepping is
!     employed depending on user inputs.
      SUBROUTINE EC_DCPLD_GETY(ec)
      USE ECMOD
      USE COMMOD, ONLY : eccModelType, IKIND, RKIND, tIntType_FE,
     2   tIntType_RK4, tIntType_BE, err, dt, time
      IMPLICIT NONE
      TYPE(eccModelType), INTENT(INOUT) :: ec

      INTEGER i, nt
      REAL(KIND=RKIND) :: t, ts

!     Total time steps
      nt = NINT(dt/ec%dt, KIND=IKIND)

!     Start time
      ts = time - dt

!     Excitation-contraction coupling due to active stress
      IF (ec%astress) THEN
         SELECT CASE (ec%odeS%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = ts + (REAL(i-1,KIND=RKIND)*ec%dt)
               CALL EC_ACTVSTRS_FE(t, ec%dt, ec%Ya)
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = ts + (REAL(i-1,KIND=RKIND)*ec%dt)
               CALL EC_ACTVSTRS_RK(t, ec%dt, ec%Ya)
            END DO

         CASE (tIntType_BE)
            DO i=1, nt
               t = ts + (REAL(i-1,KIND=RKIND)*ec%dt)
               CALL EC_ACTVSTRS_BE(t, ec%dt, ec%Ya)
            END DO

         END SELECT

         IF (ISNAN(ec%Ya)) THEN
            err = " NaN occurence (ec%Ya). Aborted!"
         END IF

!     Excitation-contraction coupling due to active strain
      ELSE IF (ec%astrain) THEN
         CALL EC_ACTVSTRN(time, ec%Ya)

      END IF

      RETURN
      END SUBROUTINE EC_DCPLD_GETY

!####################################################################

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
!     This module defines data structures for decoupled and uniformly
!     activated excitation model for excitation-contraction coupling.
!
!     Reference for active stress model:
!        Pfaller, M. R., et al.(2019). The importance of the pericardium
!        for cardiac biomechanics: from physiology to computational
!        modeling. Biomechanics and Modeling in Mechanobiology,
!        18(2), 503–529. https://doi.org/10.1007/s10237-018-1098-4
!
!     Reference for active strain model:
!        Barbarotta, L., et al.(2018). A transmurally heterogeneous
!        orthotropic activation model for ventricular contraction and
!        its numerical validation. International Journal for Numerical
!        Methods in Biomedical Engineering, 34(12), 1–24.
!        https://doi.org/10.1002/cnm.3137
!
!-----------------------------------------------------------------------

      MODULE ECMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL, pi
      IMPLICIT NONE

      PRIVATE

      INCLUDE "PARAMS_EC_DCPLD.f"

      PUBLIC :: EC_READPARFF
      PUBLIC :: EC_ACTVSTRS_FE
      PUBLIC :: EC_ACTVSTRS_RK
      PUBLIC :: EC_ACTVSTRS_BE
      PUBLIC :: EC_ACTVSTRN_SIN2

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE EC_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))

!     Excitaton-contraction coupling: active stress
      CALL GETRVAL(fid, "sigm0", sigm0)
      CALL GETRVAL(fid, "min_alpha", min_alpha)
      CALL GETRVAL(fid, "max_alpha", max_alpha)
      CALL GETRVAL(fid, "t_sys", t_sys)
      CALL GETRVAL(fid, "t_dia", t_dia)
      CALL GETRVAL(fid, "gamma", gamma)

!     Excitaton-contraction coupling: active strain
      CALL GETRVAL(fid, "gf_min", gf_min)
      CALL GETRVAL(fid, "ta_s", ta_s)
      CALL GETRVAL(fid, "ta_e", ta_e)

      CLOSE(fid)

      RETURN
      END SUBROUTINE EC_READPARFF
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using forward Euler integration
      SUBROUTINE EC_ACTVSTRS_FE(t, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) f

      CALL EC_ASTRS_GETF(t, Ta, f)
      Ta = Ta + (dt*f)

      RETURN
      END SUBROUTINE EC_ACTVSTRS_FE
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using RK4 integration
      SUBROUTINE EC_ACTVSTRS_RK(t, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: dt2, dt6, t_rk, Ta_rk, f_rk(4)

      dt2  = dt / 2._RKIND
      dt6  = dt / 6._RKIND

      t_rk  = t
      Ta_rk = Ta
      CALL EC_ASTRS_GETF(t_rk, Ta_rk, f_rk(1))

      t_rk  = t  + dt2
      Ta_rk = Ta + (dt2*f_rk(1))
      CALL EC_ASTRS_GETF(t_rk, Ta_rk, f_rk(2))

      t_rk  = t  + dt2
      Ta_rk = Ta + (dt2*f_rk(2))
      CALL EC_ASTRS_GETF(t_rk, Ta_rk, f_rk(3))

      t_rk  = t  + dt
      Ta_rk = Ta + (dt*f_rk(3))
      CALL EC_ASTRS_GETF(t_rk, Ta_rk, f_rk(4))

      Ta = Ta + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE EC_ACTVSTRS_RK
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using backward Euler integration
      SUBROUTINE EC_ACTVSTRS_BE(t, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: sp, sm, ft, a, ap

      sp = 0.5_RKIND*(1._RKIND + TANH((t-t_sys)/gamma))
      sm = 0.5_RKIND*(1._RKIND - TANH((t-t_dia)/gamma))
      ft = sp * sm

      a  = max_alpha*ft + min_alpha*(1._RKIND - ft)
      ap = MAX(a, 0._RKIND)

      Ta = (Ta + (dt*sigm0*ap)) / (1._RKIND + ABS(a)*dt)

      RETURN
      END SUBROUTINE EC_ACTVSTRS_BE
!-----------------------------------------------------------------------
      SUBROUTINE EC_ASTRS_GETF(t, Ta, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t, Ta
      REAL(KIND=RKIND), INTENT(OUT) :: f

      REAL(KIND=RKIND) :: sp, sm, ft, a, ap

      sp = 0.5_RKIND*(1._RKIND + TANH((t-t_sys)/gamma))
      sm = 0.5_RKIND*(1._RKIND - TANH((t-t_dia)/gamma))
      ft = sp * sm

      a  = max_alpha*ft + min_alpha*(1._RKIND - ft)
      ap = MAX(a, 0._RKIND)

      f  = -ABS(a)*Ta + sigm0*ap

      RETURN
      END SUBROUTINE EC_ASTRS_GETF
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     strain model (prescribed sine-squared function)
      SUBROUTINE EC_ACTVSTRN_SIN2(t, gf)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t
      REAL(KIND=RKIND), INTENT(INOUT) :: gf

      REAL(KIND=RKIND) a

      gf = 0._RKIND
      IF ((t.GE.ta_s) .AND. (t.LE.ta_e)) THEN
         a  = (t - ta_s)/(ta_e - ta_s)
         gf = gf_min * (SIN(pi*a)**2)
      END IF

      RETURN
      END SUBROUTINE EC_ACTVSTRN_SIN2
!-----------------------------------------------------------------------
      SUBROUTINE GETRVAL(fileId, skwrd, rVal)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: fileId
      CHARACTER(LEN=*), INTENT(IN) :: skwrd
      REAL(KIND=RKIND), INTENT(INOUT) :: rVal

      INTEGER(KIND=IKIND) :: slen, i, ios
      CHARACTER(LEN=stdL) :: sline, scmd, sval

      REWIND(fileId)
      slen = LEN(TRIM(skwrd))
      DO
         READ(fileId,'(A)',END=001) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         DO i=1, slen
            IF (sline(i:i) .EQ. ':') EXIT
         END DO

         IF (i .GE. slen) THEN
            STOP "Error: inconsistent input file format"
         END IF

         scmd = sline(1:i-1)
         sval = sline(i+1:slen)
         sval = ADJUSTC(sval)

!        Remove any trailing comments
         slen = LEN(TRIM(sval))
         DO i=1, slen
            IF (sval(i:i) .EQ. '#') EXIT
         END DO
         sval = TRIM(ADJUSTC(sval(1:i-1)))

         IF (TRIM(skwrd) .EQ. TRIM(scmd)) THEN
            READ(sval,*,IOSTAT=ios) rval
            IF (ios .NE. 0) THEN
               WRITE(*,'(A)') " Error: while reading "//TRIM(skwrd)
               STOP
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

! 001  WRITE(*,'(A)') " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"
!      STOP

      END SUBROUTINE GETRVAL
!-----------------------------------------------------------------------
!     Removes any leading spaces or tabs
      PURE FUNCTION ADJUSTC(str)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: str
      CHARACTER(LEN=LEN(str)) ADJUSTC

      INTEGER(KIND=IKIND) i

      DO i=1, LEN(str)
         IF (str(i:i) .NE. " " .AND. str(i:i) .NE. "  ") EXIT
      END DO
      IF (i .GT. LEN(str)) THEN
         ADJUSTC = ""
      ELSE
         ADJUSTC = str(i:)
      END IF

      RETURN
      END FUNCTION ADJUSTC
!-----------------------------------------------------------------------
      END MODULE ECMOD
!#######################################################################


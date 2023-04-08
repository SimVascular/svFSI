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
!     This module defines data structures for Aliev-Panfilov cellular
!     activation model for cardiac electrophysiology. Active stress
!     model is used for excitation-contraction coupling.
!
!     Reference for Aliev-Panfilov electrophysiology model:
!        Goktepe, S., & Kuhl, E. (2009). Computational modeling of
!        cardiac electrophysiology: A novel finite element approach.
!        Int. J. Numer. Meth. Engng, 79, 156–178.
!        https://doi.org/10.1002/nme
!
!     Reference for active stress model:
!        Goktepe, S., & Kuhl, E. (2010). Electromechanics of the heart:
!        A unified approach to the strongly coupled excitation-
!        contraction problem. Computational Mechanics, 45(2–3), 227–243.
!        https://doi.org/10.1007/s00466-009-0434-z
!
!-----------------------------------------------------------------------

      MODULE APMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      INCLUDE "PARAMS_AP.f"

      PUBLIC :: AP_INIT
      PUBLIC :: AP_READPARFF
      PUBLIC :: AP_INTEGFE
      PUBLIC :: AP_INTEGRK
      PUBLIC :: AP_INTEGCN2
      PUBLIC :: AP_ACTVSTRS_FE
      PUBLIC :: AP_ACTVSTRS_RK
      PUBLIC :: AP_ACTVSTRS_BE

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE AP_INIT(nX, X)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)

      X(1) = Voffset
      X(2) = 1.E-3_RKIND

      RETURN
      END SUBROUTINE AP_INIT
!-----------------------------------------------------------------------
      SUBROUTINE AP_READPARFF(fname)
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: fname

      INTEGER fid

      fid = 1528

      OPEN(fid, FILE=TRIM(fname))

!     Scaling factors
      CALL GETRVAL(fid, "Vscale", Vscale)
      CALL GETRVAL(fid, "Tscale", Tscale)
      CALL GETRVAL(fid, "Voffset", Voffset)

!     Cellular activation model parameters
      CALL GETRVAL(fid, "alpha", alpha)
      CALL GETRVAL(fid, "a", a)
      CALL GETRVAL(fid, "b", b)
      CALL GETRVAL(fid, "c", c)
      CALL GETRVAL(fid, "mu1", mu1)
      CALL GETRVAL(fid, "mu2", mu2)

!     Electromechanics coupling parameters
      CALL GETRVAL(fid, "Vrest", Vrest)
      CALL GETRVAL(fid, "Vcrit", Vcrit)
      CALL GETRVAL(fid, "K_T", K_T)
      CALL GETRVAL(fid, "eps_0", eps_0)
      CALL GETRVAL(fid, "eps_i", eps_i)
      CALL GETRVAL(fid, "xi_T", xi_T)

!     Electrophysiology model parameters
      CALL GETRVAL(fid, "Cm", Cm)
      CALL GETRVAL(fid, "sV", sV)
      CALL GETRVAL(fid, "rho", rho)

      CLOSE(fid)

      RETURN
      END SUBROUTINE AP_READPARFF
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE AP_INTEGFE(nX, X, Ts, Ti, Istim, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      REAL(KIND=RKIND) :: t, dt, f(nX), Isac, fext

      t    = Ts / Tscale
      dt   = Ti / Tscale
      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
      CALL AP_GETF(nX, X, f, fext)
      X = X + dt*f
      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE AP_INTEGRK(nX, X, Ts, Ti, Istim, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      REAL(KIND=RKIND) :: t, dt, dt6, Isac, fext, Xrk(nX), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      dt6  = dt / 6._RKIND

      Isac = Ksac * (Vrest - X(1))
      fext = (Istim + Isac) * Tscale / Vscale
      X(1) = (X(1) - Voffset)/Vscale

!     RK4: 1st pass
      Xrk  = X
      CALL AP_GETF(nX, Xrk, frk(:,1), fext)

!     RK4: 2nd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,1)
      CALL AP_GETF(nX, Xrk, frk(:,2), fext)

!     RK4: 3rd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,2)
      CALL AP_GETF(nX, Xrk, frk(:,3), fext)

!     RK4: 4th pass
      Xrk  = X + dt*frk(:,3)
      CALL AP_GETF(nX, Xrk, frk(:,4), fext)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE AP_INTEGRK
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE AP_INTEGCN2(nX, Xn, Ts, Ti, Istim, Ksac, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      INTEGER(KIND=IKIND), INTENT(INOUT) :: IPAR(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xn(nX), RPAR(2)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim, Ksac

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INTEGER(KIND=IKIND) :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: t, dt, Isac, fext, atol, rtol, Xk(nX), fn(nX),
     2   fk(nX), rK(nX), Im(nX,nX), JAC(nX,nX), rmsA, rmsR

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      t     = Ts / Tscale
      dt    = Ti / Tscale
      Isac  = Ksac * (Vrest - Xn(1))
      fext  = (Istim + Isac) * Tscale / Vscale
      Xn(1) = (Xn(1) - Voffset)/Vscale
      Im    = MAT_ID(nX)

      CALL AP_GETF(nX, Xn, fn, fext)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL AP_GETF(nX, Xk, fk, fext)
         rK(:) = Xk(:) - Xn(:) - 0.5_RKIND*dt*(fk(:) + fn(:))

         rmsA = 0._RKIND
         rmsR = 0._RKIND
         DO i=1, nX
            rmsA = rmsA + rK(i)**2._RKIND
            rmsR = rmsR + ( rK(i) / (Xk(i)+eps) )**2._RKIND
         END DO
         rmsA = SQRT(rmsA/REAL(nX, KIND=RKIND))
         rmsR = SQRT(rmsR/REAL(nX, KIND=RKIND))

         l1   = k .GT. itMax
         l2   = rmsA .LE. atol
         l3   = rmsR .LE. rtol
         IF (l1 .OR. l2 .OR. l3) EXIT

         CALL AP_GETJ(nX, Xk, JAC, Ksac*Tscale)
         JAC   = Im - 0.5_RKIND*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL AP_GETF(nX, Xn, fn, fext)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE AP_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETF(n, X, f, fext)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), fext
      REAL(KIND=RKIND), INTENT(OUT) :: f(n)

      f(1) = X(1)*(c*(X(1)-alpha)*(1._RKIND-X(1)) - X(2)) + fext

      f(2) = (a + mu1*X(2)/(mu2 + X(1))) *
     2       (-X(2) - c*X(1)*(X(1) - b - 1._RKIND))

      RETURN
      END SUBROUTINE AP_GETF
!-----------------------------------------------------------------------
      SUBROUTINE AP_GETJ(n, X, JAC, Ksac)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), Ksac
      REAL(KIND=RKIND), INTENT(OUT) :: JAC(n,n)

      REAL(KIND=RKIND) :: n1, n2, n3

      JAC(:,:) = 0._RKIND

      n1 = X(1) - alpha
      n2 = 1._RKIND - X(1)
      JAC(1,1) = c * (n1*n2 + X(1)*(n2 - n1)) - X(2) - Ksac

      JAC(1,2) = -X(1)

      n1 = mu1*X(2)/(mu2 + X(1))
      n2 = n1 / (mu2 + X(1))
      n3 = X(2) + c*X(1)*(X(1) - b - 1._RKIND)

      JAC(2,1) = n2*n3 - c*(a+n1)*(2._RKIND*X(1) - b - 1._RKIND)

      n1 = mu1/(mu2 + X(1))
      n2 = a + n1*X(2)
      n3 = -n3
      JAC(2,2) = n1*n3 - n2

      RETURN
      END SUBROUTINE AP_GETJ
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using forward Euler integration
      SUBROUTINE AP_ACTVSTRS_FE(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: epsX, f

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      CALL AP_ASTRS_GETF(X, epsX, Ta, f)
      Ta = Ta + (dt*f)

      RETURN
      END SUBROUTINE AP_ACTVSTRS_FE
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using RK4 integration
      SUBROUTINE AP_ACTVSTRS_RK(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: dt6, epsX, Ta_rk, f_rk(4)

      dt6  = dt / 6._RKIND
      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      Ta_rk = Ta
      CALL AP_ASTRS_GETF(X, epsX, Ta_rk, f_rk(1))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(1))
      CALL AP_ASTRS_GETF(X, epsX, Ta_rk, f_rk(2))

      Ta_rk = Ta + (0.5_RKIND*dt*f_rk(2))
      CALL AP_ASTRS_GETF(X, epsX, Ta_rk, f_rk(3))

      Ta_rk = Ta + (dt*f_rk(3))
      CALL AP_ASTRS_GETF(X, epsX, Ta_rk, f_rk(4))

      Ta = Ta + dt6*(f_rk(1) + 2._RKIND*(f_rk(2) + f_rk(3)) + f_rk(4))

      RETURN
      END SUBROUTINE AP_ACTVSTRS_RK
!-----------------------------------------------------------------------
!     Compute activation force for electromechanics based on active
!     stress model using backward Euler integration
      SUBROUTINE AP_ACTVSTRS_BE(X, dt, Ta)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, dt
      REAL(KIND=RKIND), INTENT(INOUT) :: Ta

      REAL(KIND=RKIND) :: epsX

      epsX = EXP(-EXP(-xi_T*(X - Vcrit)))
      epsX = eps_0 + (eps_i - eps_0)*epsX

      Ta = (Ta + (dt*epsX*K_T*(X-Vrest))) / (1._RKIND + epsX*dt)

      RETURN
      END SUBROUTINE AP_ACTVSTRS_BE
!-----------------------------------------------------------------------
      SUBROUTINE AP_ASTRS_GETF(X, eX, Ta, f)
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: X, eX, Ta
      REAL(KIND=RKIND), INTENT(OUT) :: f

      f = eX*(K_T*(X-Vrest) - Ta)

      RETURN
      END SUBROUTINE AP_ASTRS_GETF
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
      END MODULE APMOD
!#######################################################################


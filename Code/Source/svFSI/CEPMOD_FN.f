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
!     This module defines data structures for Fitzhugh-Nagumo cellular
!     activation model for cardiac electrophysiology.
!
!     Reference for Aliev-Panfilov electrophysiology model:
!        Goktepe, S., & Kuhl, E. (2009). Computational modeling of
!        cardiac electrophysiology: A novel finite element approach.
!        Int. J. Numer. Meth. Engng, 79, 156–178.
!        https://doi.org/10.1002/nme
!
!-----------------------------------------------------------------------

      MODULE FNMOD
      USE TYPEMOD
      USE UTILMOD, ONLY : stdL
      IMPLICIT NONE

      PRIVATE

      INCLUDE "PARAMS_FN.f"

      PUBLIC :: FN_INIT
      PUBLIC :: FN_READPARFF
      PUBLIC :: FN_INTEGFE
      PUBLIC :: FN_INTEGRK
      PUBLIC :: FN_INTEGCN2

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE FN_INIT(nX, X)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)

      X(:) = 1.E-3_RKIND

      RETURN
      END SUBROUTINE FN_INIT
!-----------------------------------------------------------------------
      SUBROUTINE FN_READPARFF(fname)
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

      CLOSE(fid)

      RETURN
      END SUBROUTINE FN_READPARFF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Time integration performed using Forward Euler method
      SUBROUTINE FN_INTEGFE(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim

      REAL(KIND=RKIND) :: t, dt, f(nX), fext

      t    = Ts / Tscale
      dt   = Ti / Tscale
      fext = Istim * Tscale / Vscale

      X(1) = (X(1) - Voffset)/Vscale
      CALL FN_GETF(nX, X, f, fext)
      X(:) = X(:) + dt*f(:)
      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE FN_INTEGFE
!-----------------------------------------------------------------------
!     Time integration performed using 4th order Runge-Kutta method
      SUBROUTINE FN_INTEGRK(nX, X, Ts, Ti, Istim)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim

      REAL(KIND=RKIND) :: t, dt, dt6, fext, Xrk(nX), frk(nX,4)

      t    = Ts / Tscale
      dt   = Ti / Tscale
      dt6  = dt / 6._RKIND

      fext = Istim * Tscale / Vscale
      X(1) = (X(1) - Voffset)/Vscale

!     RK4: 1st pass
      Xrk  = X
      CALL FN_GETF(nX, Xrk, frk(:,1), fext)

!     RK4: 2nd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,1)
      CALL FN_GETF(nX, Xrk, frk(:,2), fext)

!     RK4: 3rd pass
      Xrk  = X + 0.5_RKIND*dt*frk(:,2)
      CALL FN_GETF(nX, Xrk, frk(:,3), fext)

!     RK4: 4th pass
      Xrk  = X + dt*frk(:,3)
      CALL FN_GETF(nX, Xrk, frk(:,4), fext)

      X = X + dt6*(frk(:,1) + 2._RKIND*(frk(:,2) + frk(:,3)) + frk(:,4))

      X(1) = X(1)*Vscale + Voffset

      RETURN
      END SUBROUTINE FN_INTEGRK
!-----------------------------------------------------------------------
!     Time integration performed using Crank-Nicholson method
      SUBROUTINE FN_INTEGCN2(nX, Xn, Ts, Ti, Istim, IPAR, RPAR)
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nX
      INTEGER(KIND=IKIND), INTENT(INOUT) :: IPAR(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: Xn(nX), RPAR(2)
      REAL(KIND=RKIND), INTENT(IN) :: Ts, Ti, Istim

      REAL(KIND=RKIND), PARAMETER :: eps = EPSILON(eps)

      INTEGER(KIND=IKIND) :: i, k, itMax
      LOGICAL :: l1, l2, l3
      REAL(KIND=RKIND) :: t, dt, fext, atol, rtol, Xk(nX), fn(nX),
     2   fK(nX), rK(nX), Im(nX,nX), JAC(nX,nX), rmsA, rmsR

      itMax = IPAR(1)
      atol  = RPAR(1)
      rtol  = RPAR(2)

      t     = Ts / Tscale
      dt    = Ti / Tscale
      fext  = Istim * Tscale / Vscale
      Xn(1) = (Xn(1) - Voffset)/Vscale
      Im    = MAT_ID(nX)

      CALL FN_GETF(nX, Xn, fn, fext)

      k  = 0
      Xk = Xn
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      t  = Ts + dt
      DO
         k = k + 1
         CALL FN_GETF(nX, Xk, fk, fext)
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

         CALL FN_GETJ(nX, Xk, JAC)
         JAC   = Im - 0.5_RKIND*dt*JAC
         JAC   = MAT_INV(JAC, nX)
         rK(:) = MATMUL(JAC, rK)
         Xk(:) = Xk(:) - rK(:)
      END DO
      Xn(:) = Xk(:)
      CALL FN_GETF(nX, Xn, fn, fext)
      Xn(1) = Xn(1)*Vscale + Voffset

      IF (.NOT.l2 .AND. .NOT.l3) IPAR(2) = IPAR(2) + 1

      RETURN
      END SUBROUTINE FN_INTEGCN2
!-----------------------------------------------------------------------
      SUBROUTINE FN_GETF(n, X, f, fext)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n), fext
      REAL(KIND=RKIND), INTENT(OUT) :: f(n)

      f(1) = c * ( X(1)*(X(1)-alpha)*(1._RKIND-X(1)) - X(2) ) + fext

      f(2) = X(1) - b*X(2) + a

      RETURN
      END SUBROUTINE FN_GETF
!-----------------------------------------------------------------------
      SUBROUTINE FN_GETJ(n, X, JAC)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: n
      REAL(KIND=RKIND), INTENT(IN) :: X(n)
      REAL(KIND=RKIND), INTENT(OUT) :: JAC(n,n)

      REAL(KIND=RKIND) :: n1, n2

      n1 = -3._RKIND*X(1)**2._RKIND
      n2 = 2._RKIND*(1._RKIND+alpha)*X(1)
      JAC(1,1) = c * (n1 + n2 - 1._RKIND)
      JAC(1,2) = -c
      JAC(2,1) = 1._RKIND
      JAC(2,2) = -b

      RETURN
      END SUBROUTINE FN_GETJ
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
               STOP " Error: while reading "//TRIM(skwrd)
            END IF
            EXIT
         END IF
      END DO

 001  RETURN

!  001  STOP " Error: EOF reached while finding command <"//
!     2   TRIM(skwrd)//">"

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
      END MODULE FNMOD
!#######################################################################


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
!     This routine embodies formulation for solving electrophysiology
!     cellular activation model at the nodal level.
!
!-----------------------------------------------------------------------

!     Initialize Cardiac Electrophysiology Model
      SUBROUTINE CEPINIT()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER a, iEq, iDmn, cPhys, dID, nX

      REAL(KIND=8), ALLOCATABLE :: Xl(:), sA(:), sF(:,:)

      DO iEq=1, nEq
         IF (eq(iEq)%phys .NE. phys_CEP) CYCLE

         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sA(tnNo), sF(nXion,tnNo))
            sA = 0D0
            sF = 0D0
            DO a=1, tnNo
               IF (.NOT.ISDOMAIN(iEq, a, phys_CEP)) CYCLE
               DO iDmn=1, eq(iEq)%nDmn
                  cPhys = eq(iEq)%dmn(iDmn)%phys
                  dID   = eq(iEq)%dmn(iDmn)%Id
                  IF (cPhys.NE.phys_CEP .OR. .NOT.BTEST(dmnId(a),dID))
     2                CYCLE
                  nX = eq(iEq)%dmn(iDmn)%cep%nX
                  ALLOCATE(Xl(nX))
                  CALL CEPINITL(eq(iEq)%dmn(iDmn)%cep, nX, Xl)
                  sA(a) = sA(a) + 1.0D0
                  sF(1:nX,a) = sF(1:nX,a) + Xl(:)
                  DEALLOCATE(Xl)
               END DO
            END DO
            CALL COMMU(sA)
            CALL COMMU(sF)
            DO a=1, tnNo
               IF (.NOT.ISZERO(sA(a)))
     2            Xion(:,a) = sF(:,a)/sA(a)
            END DO
            DEALLOCATE(sA, sF)
         ELSE
            DO a=1, tnNo
               IF (.NOT.ISDOMAIN(iEq, a, phys_CEP)) CYCLE
               nX = eq(iEq)%dmn(1)%cep%nX
               ALLOCATE(Xl(nX))
               CALL CEPINITL(eq(iEq)%dmn(1)%cep, nX, Xl)
               Xion(1:nX,a) = Xl(:)
               DEALLOCATE(Xl)
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE CEPINIT
!-----------------------------------------------------------------------
      SUBROUTINE CEPINITL(cep, nX, X)
      USE CEPMOD
      IMPLICIT NONE
      TYPE(cepModelType), INTENT(IN) :: cep
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(OUT) :: X(nX)

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         CALL AP_INIT(nX, X)

      CASE (cepModel_FN)
         CALL FN_INIT(nX, X)

      CASE (cepModel_TTP)
         CALL TTP_INIT(cep%imyo, nX, X)

      CASE (cepModel_BO)
         CALL BO_INIT(nX, X)

      END SELECT

      RETURN
      END SUBROUTINE CEPINITL
!#######################################################################
!     State variable integration
      SUBROUTINE CEPINTEG(iEq, iDof, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iEq, iDof
      REAL(KIND=8), INTENT(IN) :: Dg(tDof,tnNo)

      LOGICAL :: IPASS = .TRUE.
      INTEGER :: a, Ac, iM, iDmn, cPhys, dID, nX
      REAL(KIND=8) :: yl

      REAL(KIND=8), ALLOCATABLE :: Xl(:), I4f(:), sA(:), sY(:), sF(:,:)
      SAVE IPASS

      ALLOCATE(I4f(tnNo))
      I4f = 0D0

!     Electromechanics: get fiber stretch for stretch activated currents
      IF (cem%cpld) THEN
         DO iM=1, nMsh
            IF (msh(iM)%nFn .NE. 0) THEN
               ALLOCATE(sA(msh(iM)%nNo))
               sA = 0D0
               CALL FIBSTRETCH(iEq, msh(iM), Dg, sA)
               DO a=1, msh(iM)%nNo
                  Ac = msh(iM)%gN(a)
                  I4f(Ac) = sA(a)
               END DO
               DEALLOCATE(sA)
            END IF
         END DO
      END IF

!     Ignore first pass as Xion is already initialized
      IF (IPASS) THEN
         IPASS = .FALSE.
      ELSE
!        Copy action potential after diffusion as first state variable
         DO Ac=1, tnNo
            Xion(1,Ac) = Yo(iDof,Ac)
         END DO
      END IF

!     Integrate electric potential based on cellular activation model
      IF (ALLOCATED(dmnId)) THEN
         ALLOCATE(sA(tnNo), sF(nXion,tnNo), sY(tnNo))
         sA = 0D0
         sF = 0D0
         sY = 0D0
         DO Ac=1, tnNo
            IF (.NOT.ISDOMAIN(iEq, Ac, phys_CEP)) CYCLE
            DO iDmn=1, eq(iEq)%nDmn
               cPhys = eq(iEq)%dmn(iDmn)%phys
               dID   = eq(iEq)%dmn(iDmn)%Id
               IF (cPhys.NE.phys_CEP .OR. .NOT.BTEST(dmnId(Ac),dID))
     2             CYCLE
               nX = eq(iEq)%dmn(iDmn)%cep%nX
               ALLOCATE(Xl(nX))
               Xl = Xion(1:nX,Ac)
               yl = 0D0
               IF (cem%cpld) yl = cem%Ya(Ac)
               CALL CEPINTEGL(eq(iEq)%dmn(iDmn)%cep, nX, Xl, time-dt,
     2            yl, I4f(Ac))
               sA(Ac) = sA(Ac) + 1.0D0
               sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
               IF (cem%cpld) sY(Ac) = sY(Ac) + yl
               DEALLOCATE(Xl)
            END DO
         END DO
         CALL COMMU(sA)
         CALL COMMU(sF)
         CALL COMMU(sY)
         DO Ac=1, tnNo
            IF (.NOT.ISZERO(sA(Ac))) THEN
               Xion(:,Ac) = sF(:,Ac)/sA(Ac)
               IF (cem%cpld) cem%Ya(Ac) = sY(Ac) / sA(Ac)
            END IF
         END DO
         DEALLOCATE(sA, sF, sY)
      ELSE
         DO Ac=1, tnNo
            IF (.NOT.ISDOMAIN(iEq, Ac, phys_CEP)) CYCLE
            nX = eq(iEq)%dmn(1)%cep%nX
            ALLOCATE(Xl(nX))
            Xl = Xion(1:nX,Ac)
            yl = 0D0
            IF (cem%cpld) yl = cem%Ya(Ac)
            CALL CEPINTEGL(eq(iEq)%dmn(1)%cep, nX, Xl, time-dt, yl,
     2         I4f(Ac))
            Xion(1:nX,Ac) = Xl(:)
            IF (cem%cpld) cem%Ya(Ac) = yl
            DEALLOCATE(Xl)
         END DO
      END IF

      DO Ac=1, tnNo
         Yo(iDof,Ac) = Xion(1,Ac)
      END DO

      DEALLOCATE(I4f)

      RETURN
      END SUBROUTINE CEPINTEG
!-----------------------------------------------------------------------
!     Integrate local electrophysiology variables from t1 to t1+dt. Also
!     integrate excitation-activation variables form coupled electro-
!     mechanics. The equations are integrated at domain nodes.
      SUBROUTINE CEPINTEGL(cep, nX, X, t1, yl, I4f)
      USE CEPMOD
      USE UTILMOD, ONLY : eps
      USE COMMOD, ONLY : dt
      IMPLICIT NONE
      TYPE(cepModelType), INTENT(IN) :: cep
      INTEGER, INTENT(IN) :: nX
      REAL(KIND=8), INTENT(IN) :: t1, I4f
      REAL(KIND=8), INTENT(INOUT) :: X(nX), yl

      INTEGER i, icl, nt
      REAL(KIND=8) :: t, Ts, Te, Istim, Ksac

      INTEGER, ALLOCATABLE :: IPAR(:)
      REAL(KIND=8), ALLOCATABLE :: RPAR(:)

!     Feedback coefficient for stretch-activated-currents
      IF (I4f .GT. 1.0D0) THEN
         Ksac = cep%Ksac * (SQRT(I4f) - 1.0D0)
      ELSE
         Ksac = 0D0
      END IF

!     Total time steps
      nt = NINT(dt/cep%dt)

!     External stimulus duration
      icl = MAX(FLOOR(t1/cep%Istim%CL), 0)
      Ts  = cep%Istim%Ts + REAL(icl, KIND=8)*cep%Istim%CL
      Te  = Ts + cep%Istim%Td

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL AP_INTEGFE(nX, X, t, cep%dt, Istim, Ksac)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL AP_ACTVSTRS(X(1), cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL AP_INTEGRK(nX, X, t, cep%dt, Istim, Ksac)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL AP_ACTVSTRS(X(1), cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL AP_INTEGCN2(nX, X, t, cep%dt, Istim, Ksac, IPAR,
     2            RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL AP_ACTVSTRS(X(1), cep%dt, yl)
               END IF
            END DO
         END SELECT

      CASE (cepModel_FN)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL FN_INTEGFE(nX, X, t, cep%dt, Istim)
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL FN_INTEGRK(nX, X, t, cep%dt, Istim)
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL FN_INTEGCN2(nX, X, t, cep%dt, Istim, IPAR, RPAR)
            END DO
         END SELECT

      CASE (cepModel_TTP)
         ALLOCATE(IPAR(2), RPAR(18))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL TTP_INTEGFE(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL TTP_ACTVSTRS(X(4), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL TTP_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL TTP_INTEGRK(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL TTP_ACTVSTRS(X(4), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL TTP_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL TTP_INTEGCN2(cep%imyo, nX, X, t, cep%dt, Istim,
     2            Ksac, IPAR, RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL TTP_ACTVSTRS(X(4), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL TTP_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO
         END SELECT

      CASE (cepModel_BO)
         ALLOCATE(IPAR(2), RPAR(5))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0D0
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL BO_INTEGFE(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL BO_ACTVSTRS(X(1), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL BO_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL BO_INTEGRK(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL BO_ACTVSTRS(X(1), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL BO_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1,KIND=8) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0D0
               END IF
               CALL BO_INTEGCN2(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            IPAR, RPAR)

!              Electromechanics excitation-activation
               IF (cem%aStress) THEN
                  CALL BO_ACTVSTRS(X(1), cep%dt, yl)
               ELSE IF (cem%aStrain) THEN
                  CALL BO_ACTVSTRN(X(4), I4f, cep%dt, yl)
               END IF
            END DO
         END SELECT

      END SELECT

      IF (ISNAN(X(1)) .OR. ISNAN(yl)) THEN
         WRITE(*,'(A)') " NaN occurence. Aborted!"
         CALL STOPSIM()
      END IF

      DEALLOCATE(IPAR, RPAR)

      RETURN
      END SUBROUTINE CEPINTEGL
!####################################################################

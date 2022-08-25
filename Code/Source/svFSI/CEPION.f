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

      INTEGER(KIND=IKIND) a, iEq, iDmn, cPhys, dID, nX, nG

      REAL(KIND=RKIND), ALLOCATABLE :: Xl(:), Xgl(:), sA(:), sF(:,:)

      DO iEq=1, nEq
         IF (eq(iEq)%phys .NE. phys_CEP) CYCLE

!        Average cellular activation state variables at domain
!        interfaces for multi-domain problem
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sA(tnNo), sF(nXion,tnNo))
            sA = 0._RKIND
            sF = 0._RKIND
            DO a=1, tnNo
               IF (.NOT.ISDOMAIN(iEq, a, phys_CEP)) CYCLE
               DO iDmn=1, eq(iEq)%nDmn
                  cPhys = eq(iEq)%dmn(iDmn)%phys
                  dID   = eq(iEq)%dmn(iDmn)%Id
                  IF (cPhys.NE.phys_CEP .OR. .NOT.BTEST(dmnId(a),dID))
     2                CYCLE
                  nX = eq(iEq)%dmn(iDmn)%cep%nX
                  nG = eq(iEq)%dmn(iDmn)%cep%nG
                  ALLOCATE(Xl(nX), Xgl(nG))
                  CALL CEPINITL(eq(iEq)%dmn(iDmn)%cep, nX, nG, Xl, Xgl)
                  sA(a) = sA(a) + 1._RKIND
                  sF(1:nX,a)  = sF(1:nX,a) + Xl(:)
                  sF(nX+1:nX+nG,a) = sF(nX+1:nX+nG,a) + Xgl(:)
                  DEALLOCATE(Xl, Xgl)
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
               nG = eq(iEq)%dmn(1)%cep%nG
               ALLOCATE(Xl(nX), Xgl(nG))
               CALL CEPINITL(eq(iEq)%dmn(1)%cep, nX, nG, Xl, Xgl)
               Xion(1:nX,a) = Xl(:)
               Xion(nX+1:nX+nG,a) = Xgl(:)
               DEALLOCATE(Xl, Xgl)
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE CEPINIT
!-----------------------------------------------------------------------
      SUBROUTINE CEPINITL(cep, nX, nG, X, Xg)
      USE CEPMOD
      USE COMMOD, ONLY : std
      IMPLICIT NONE
      TYPE(cepModelType), INTENT(IN) :: cep
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      REAL(KIND=RKIND), INTENT(OUT) :: X(nX), Xg(nG)

      INTEGER(KIND=IKIND) slen

      slen = LEN(TRIM(cep%fpar_in))
      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         CALL AP_INIT(nX, X)
         IF (slen .GT. 0) THEN
            CALL AP_READPARFF(cep%fpar_in)
            std = " Reading Aliev-Panfilov model parameters"//
     2         " from input file"
         END IF

      CASE (cepModel_BO)
         CALL BO_INIT(nX, X)
         IF (slen .GT. 0) THEN
            CALL BO_READPARFF(cep%fpar_in)
            std = " Reading Bueno-Orovio model parameters"//
     2         " from input file"
         END IF

      CASE (cepModel_FN)
         CALL FN_INIT(nX, X)
         IF (slen .GT. 0) THEN
            CALL FN_READPARFF(cep%fpar_in)
            std = " Reading Fitzhugh-Nagumo model parameters"//
     2         " from input file"
         END IF

      CASE (cepModel_TTP)
         CALL TTP_INIT(cep%imyo, nX, nG, X, Xg)
         IF (slen .GT. 0) THEN
            CALL TTP_READPARFF(cep%fpar_in)
            std = " Reading tenTusscher-Panfilov model parameters"//
     2         " from input file"
         END IF

      END SELECT

      RETURN
      END SUBROUTINE CEPINITL
!#######################################################################
!     Time integration of the state variables in the cellular activation
!     model for cardiac electrophysiology. This function is called
!     during the predictor step in PIC.f (subroutine PICP) and only
!     when the equation solved is `cep'.
      SUBROUTINE CEPINTEG(iEq, iDof, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, iDof
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      LOGICAL :: flag, IPASS = .TRUE.
      INTEGER(KIND=IKIND) :: a, Ac, iM, iDmn, cPhys, dID, nX, nG
      REAL(KIND=RKIND) :: yl

      REAL(KIND=RKIND), ALLOCATABLE :: I4f(:), Xl(:), Xgl(:), sA(:),
     2   sY(:), sF(:,:)
      SAVE IPASS

      ALLOCATE(I4f(tnNo))
      I4f = 0._RKIND

!     Excitation-contraction coupling: get fiber stretch (I4f) for
!     stretch activated currents and also for active strain model
      IF (ecCpld) THEN
         DO iM=1, nMsh
            IF (msh(iM)%nFn .NE. 0) THEN
               ALLOCATE(sA(msh(iM)%nNo))
               sA = 0._RKIND
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
!     Copy action potential after diffusion as the first state variable
         DO Ac=1, tnNo
            Xion(1,Ac) = Yo(iDof,Ac)
         END DO
      END IF

!     Integrate activation potential based on cellular activation model
!     In a multi-domain problem, integrate action potential at each node
!     and average them across domain interfaces
      IF (ALLOCATED(dmnId)) THEN
         ALLOCATE(sA(tnNo), sF(nXion,tnNo), sY(tnNo))
         sA = 0._RKIND
         sF = 0._RKIND
         sY = 0._RKIND
         DO Ac=1, tnNo
            IF (.NOT.ISDOMAIN(iEq, Ac, phys_CEP)) CYCLE
            DO iDmn=1, eq(iEq)%nDmn
               cPhys = eq(iEq)%dmn(iDmn)%phys
               dID   = eq(iEq)%dmn(iDmn)%Id

               flag = (cPhys .NE. phys_CEP)
     2           .OR. (.NOT.BTEST(dmnId(Ac),dID))
     3           .OR. (.NOT.eq(iEq)%dmn(iDmn)%ec%caCpld)
               IF (flag) CYCLE

               nX = eq(iEq)%dmn(iDmn)%cep%nX
               nG = eq(iEq)%dmn(iDmn)%cep%nG

!              Local copies
               ALLOCATE(Xl(nX), Xgl(nG))
               Xl  = Xion(1:nX,Ac)
               Xgl = Xion(nX+1:nX+nG,Ac)

!              Excitation-contraction
               yl  = 0._RKIND
               IF (ecCpld) yl = ec_Ya(Ac)

               CALL CEPINTEGL(time-dt, nX, nG, eq(iEq)%dmn(iDmn)%cep,
     2            eq(iEq)%dmn(iDmn)%ec, I4f(Ac), Xl, Xgl, yl)

               sA(Ac) = sA(Ac) + 1._RKIND
               sF(1:nX,Ac) = sF(1:nX,Ac) + Xl(:)
               sF(nX+1:nX+nG,Ac) = sF(nX+1:nX+nG,Ac) + Xgl(:)
               IF (ecCpld) sY(Ac) = sY(Ac) + yl

               DEALLOCATE(Xl, Xgl)
            END DO
         END DO

         CALL COMMU(sA)
         CALL COMMU(sF)
         IF (ecCpld) CALL COMMU(sY)

         DO Ac=1, tnNo
            IF (.NOT.ISZERO(sA(Ac))) THEN
               Xion(:,Ac) = sF(:,Ac)/sA(Ac)
               IF (ecCpld) ec_Ya(Ac) = sY(Ac) / sA(Ac)
            END IF
         END DO

         DEALLOCATE(sA, sF, sY)
      ELSE
         DO Ac=1, tnNo
            flag = (.NOT.ISDOMAIN(iEq, Ac, phys_CEP))
     2        .OR. (.NOT.eq(iEq)%dmn(1)%ec%caCpld)
            IF (flag) CYCLE

            nX = eq(iEq)%dmn(1)%cep%nX
            nG = eq(iEq)%dmn(1)%cep%nG

!           Local copies
            ALLOCATE(Xl(nX), Xgl(nG))
            Xl  = Xion(1:nX,Ac)
            Xgl = Xion(nX+1:nX+nG,Ac)

!           Excitation-contraction coupling
            yl = 0._RKIND
            IF (ecCpld) yl = ec_Ya(Ac)

            CALL CEPINTEGL(time-dt, nX, nG, eq(iEq)%dmn(1)%cep,
     2         eq(iEq)%dmn(1)%ec, I4f(Ac), Xl, Xgl, yl)

            Xion(1:nX,Ac) = Xl(:)
            Xion(nX+1:nX+nG,Ac) = Xgl(:)
            IF (ecCpld) ec_Ya(Ac) = yl

            DEALLOCATE(Xl, Xgl)
         END DO
      END IF

      DO Ac=1, tnNo
         Yo(iDof,Ac) = Xion(1,Ac)
      END DO

      DEALLOCATE(I4f)

      RETURN
      END SUBROUTINE CEPINTEG
!-----------------------------------------------------------------------
!     Integrate local electrophysiology state variables from t1 to t1+dt
!     Also integrate excitation-contraction variables for electro-
!     mechanics modeling. The equations are integrated at domain nodes.
      SUBROUTINE CEPINTEGL(t1, nX, nG, cep, ec, I4f, X, Xg, yl)
      USE CEPMOD
      USE UTILMOD, ONLY : eps
      USE COMMOD, ONLY : dt, err, eccModelType
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: t1, I4f
      INTEGER(KIND=IKIND), INTENT(IN) :: nX, nG
      TYPE(cepModelType), INTENT(IN) :: cep
      TYPE(eccModelType), INTENT(IN) :: ec
      REAL(KIND=RKIND), INTENT(INOUT) :: X(nX), Xg(nG), yl

      INTEGER(KIND=IKIND) i, icl, nt
      REAL(KIND=RKIND) :: t, Ts, Te, Istim, Ksac, X0

      INTEGER(KIND=IKIND), ALLOCATABLE :: IPAR(:)
      REAL(KIND=RKIND), ALLOCATABLE :: RPAR(:)

!     Feedback coefficient for stretch-activated-currents
      IF (I4f .GT. 1._RKIND) THEN
         Ksac = cep%Ksac * (SQRT(I4f) - 1._RKIND)
      ELSE
         Ksac = 0._RKIND
      END IF

!     Total time steps
      nt = NINT(dt/cep%dt, KIND=IKIND)

!     External stimulus duration
      icl = MAX( FLOOR(t1/cep%Istim%CL, KIND=IKIND), 0)
      Ts  = cep%Istim%Ts + REAL(icl, KIND=RKIND)*cep%Istim%CL
      Te  = Ts + cep%Istim%Td

      SELECT CASE (cep%cepType)
      CASE (cepModel_AP)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0._RKIND
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

               X0 = X(1)
               CALL AP_INTEGFE(nX, X, t, cep%dt, Istim, Ksac)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

               X0 = X(1)
               CALL AP_INTEGRK(nX, X, t, cep%dt, Istim, Ksac)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

               X0 = X(1)
               CALL AP_INTEGCN2(nX, X, t, cep%dt, Istim, Ksac, IPAR,
     2            RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL AP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL AP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL AP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF
            END DO
         END SELECT

      CASE (cepModel_BO)
         ALLOCATE(IPAR(2), RPAR(5))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0._RKIND
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (ec%astress) THEN
                  X0 = X(1)
               ELSE IF (ec%astrain) THEN
                  X0 = X(4)
               END IF

               CALL BO_INTEGFE(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (ec%astress) THEN
                  X0 = X(1)
               ELSE IF (ec%astrain) THEN
                  X0 = X(4)
               END IF

               CALL BO_INTEGRK(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Copy old state variable for explicit coupling with
!              excitation-contraction model
               IF (ec%astress) THEN
                  X0 = X(1)
               ELSE IF (ec%astrain) THEN
                  X0 = X(4)
               END IF

               CALL BO_INTEGCN2(cep%imyo, nX, X, t, cep%dt, Istim, Ksac,
     2            IPAR, RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(1)
                     CALL BO_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL BO_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL BO_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL BO_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO

         END SELECT

      CASE (cepModel_FN)
         ALLOCATE(IPAR(2), RPAR(2))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0._RKIND
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF
               CALL FN_INTEGFE(nX, X, t, cep%dt, Istim)
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF
               CALL FN_INTEGRK(nX, X, t, cep%dt, Istim)
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF
               CALL FN_INTEGCN2(nX, X, t, cep%dt, Istim, IPAR, RPAR)
            END DO
         END SELECT

      CASE (cepModel_TTP)
         ALLOCATE(IPAR(2), RPAR(18))
         IPAR(1) = cep%odes%maxItr
         IPAR(2) = 0
         RPAR(:) = 0._RKIND
         RPAR(1) = cep%odes%absTol
         RPAR(2) = cep%odes%relTol

         SELECT CASE (cep%odes%tIntType)
         CASE (tIntType_FE)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGFE(cep%imyo, nX, nG, X, Xg, cep%dt, Istim,
     2            Ksac, RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO

         CASE (tIntType_RK4)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGRK(cep%imyo, nX, nG, X, Xg, cep%dt, Istim,
     2            Ksac, RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO

         CASE (tIntType_CN2)
            DO i=1, nt
               t = t1 + REAL(i-1, KIND=RKIND) * cep%dt
               IF (t.GE.Ts-eps .AND. t.LE.Te+eps) THEN
                  Istim = cep%Istim%A
               ELSE
                  Istim = 0._RKIND
               END IF

!              Integrate local state variables
               X0 = X(4)
               CALL TTP_INTEGCN2(cep%imyo, nX, nG, X, Xg, cep%dt, Istim,
     2            Ksac, IPAR, RPAR)

!              Excitation-contraction coupling due to active stress
               IF (ec%astress) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRS_FE(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRS_RK(X0, cep%dt, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRS_BE(X0, cep%dt, yl)
                  END IF
               END IF

!              Excitation-contraction coupling due to active strain
               IF (ec%astrain) THEN
                  IF (ec%odeS%tIntType .EQ. tIntType_FE) THEN
                     CALL TTP_ACTVSTRN_FE(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_RK4) THEN
                     CALL TTP_ACTVSTRN_RK(X0, cep%dt, I4f, yl)

                  ELSE IF (ec%odeS%tIntType .EQ. tIntType_BE) THEN
                     X0 = X(4)
                     CALL TTP_ACTVSTRN_BE(X0, cep%dt, I4f, yl,
     2                  ec%odeS%maxItr, ec%odeS%absTol, ec%odeS%relTol)
                  END IF
               END IF
            END DO
         END SELECT

      END SELECT

      IF (ISNAN(X(1)) .OR. ISNAN(yl)) THEN
         err = " NaN occurence. Aborted!"
      END IF

      DEALLOCATE(IPAR, RPAR)

      RETURN
      END SUBROUTINE CEPINTEGL
!####################################################################

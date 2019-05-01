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
!     This routine contains predictor/initiator/corrector routines, as
!     part of time integration scheme.
!
!--------------------------------------------------------------------

!     This is the predictor
      SUBROUTINE PICP
      USE COMMOD
      IMPLICIT NONE

      INTEGER iEq, s, e
      REAL(KIND=8) coef, ctime

!     Prestress initialization
      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_preSt) THEN
            pS0 = pS0 + pSn
            eq(iEq)%iNorm = 0D0
            eq(iEq)%pNorm = HUGE(coef)
            Ao = 0D0
            Yo = 0D0
            Do = 0D0
         END IF
      END DO

!     Immersed body treatment (explicit)
      IF (ibFlag) THEN
         ib%callD(1) = CPUT()

!        Set IB forces to zero, except for feedback force
         ib%R = 0D0

!        Compute FSI forcing (ib%R) for immersed bodies (IFEM)
         CALL IB_CALCFFSI(Do)

!        Treat IB dirichlet boundaries using penalty forces
c         CALL IB_SETBCPEN()

!        Update IB location and tracers
         ib%Ao = ib%An
         ib%Yo = ib%Yn
         ib%Uo = ib%Un
         ib%callD(2) = CPUT()
         CALL IB_UPDATE(Do)
         ctime = CPUT()
         ib%callD(2) = ctime - ib%callD(2)
         ib%callD(1) = ctime - ib%callD(1)
      END IF

      DO iEq=1, nEq
         s = eq(iEq)%s
         e = eq(iEq)%e
         coef = (eq(iEq)%gam - 1D0)/eq(iEq)%gam
         An(s:e,:) = Ao(s:e,:)*coef

!        electrophysiology
         IF (eq(iEq)%phys .EQ. phys_CEP) THEN
            CALL CEPION(iEq, e)
         END IF

         Yn(s:e,:) = Yo(s:e,:)

         IF (dFlag) THEN
            IF ( .NOT. sstEq) THEN
!              struct, lElas, FSI (struct, mesh)
               coef = dt*dt*(5D-1*eq(iEq)%gam - eq(iEq)%beta)
     2            /(eq(iEq)%gam - 1D0)
               Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
            ELSE
!              vms_struct, FSI
               IF (eq(iEq)%phys .EQ. phys_vms_struct .OR.
     2             eq(iEq)%phys .EQ. phys_FSI) THEN
                  coef = (eq(iEq)%gam - 1D0)/eq(iEq)%gam
                  Ad(:,:)   = Ad(:,:)*coef
                  Dn(s:e,:) = Do(s:e,:)
               ELSE IF (eq(iEq)%phys .EQ. phys_mesh) THEN
!              mesh
                  coef = dt*dt*(5D-1*eq(iEq)%gam - eq(iEq)%beta)
     2               /(eq(iEq)%gam - 1D0)
                  Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
               ELSE
                  err = "Undefined physics (PICP)"
               END IF
            END IF
         ELSE
            Dn(s:e,:) = Do(s:e,:)
         END IF
      END DO

      RETURN
      END SUBROUTINE PICP
!====================================================================
!     This is the initiator
      SUBROUTINE PICI(Ag, Yg, Dg)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER s, e, i, a
      REAL(KIND=8) coef(4)

      dof         = eq(cEq)%dof
      eq(cEq)%itr = eq(cEq)%itr + 1

      DO i=1, nEq
         s       = eq(i)%s
         e       = eq(i)%e
         coef(1) = 1D0 - eq(i)%am
         coef(2) = eq(i)%am
         coef(3) = 1D0 - eq(i)%af
         coef(4) = eq(i)%af

         DO a=1, tnNo
            Ag(s:e,a) = Ao(s:e,a)*coef(1) + An(s:e,a)*coef(2)
            Yg(s:e,a) = Yo(s:e,a)*coef(3) + Yn(s:e,a)*coef(4)
            Dg(s:e,a) = Do(s:e,a)*coef(3) + Dn(s:e,a)*coef(4)
         END DO

!        Reset pressure variable initiator
         IF ( (eq(i)%phys .EQ. phys_fluid .OR.
     2         eq(i)%phys .EQ. phys_CMM   .OR.
     3         eq(i)%phys .EQ. phys_FSI) .AND. .NOT.sstEq ) THEN
            DO a=1, tnNo
               Yg(e,a) = Yn(e,a)
            END DO
         END IF

         IF (eq(i)%phys .EQ. phys_preSt) pSn(:,:) = 0D0
      END DO

      RETURN
      END SUBROUTINE PICI
!====================================================================
!     This is the corrector. Decision for next eqn is also made here
      SUBROUTINE PICC
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: l1, l2, l3, l4, l5
      INTEGER :: s, e, a, Ac
      REAL(KIND=8) :: coef(5), r1, r2, dUl(nsd)
      CHARACTER(LEN=stdL) :: sCmd

      s       = eq(cEq)%s
      e       = eq(cEq)%e
      coef(1) = eq(cEq)%gam*dt
      coef(2) = eq(cEq)%af*coef(1)
      coef(3) = eq(cEq)%beta*dt*dt
      coef(4) = 1D0 / eq(cEq)%am
      coef(5) = coef(2)*coef(4)

      IF ( (eq(cEq)%phys .EQ. phys_fluid .OR.
     2      eq(cEq)%phys .EQ. phys_CMM   .OR.
     3      eq(cEq)%phys .EQ. phys_FSI) .AND. .NOT.sstEq ) THEN
         DO a=1, tnNo
            An(s:e,a)   = An(s:e,a)   - R(:,a)
            Yn(s:e-1,a) = Yn(s:e-1,a) - R(1:dof-1,a)*coef(1)
            Yn(e,a)     = Yn(e,a)     - R(dof,a)*coef(2)
            Dn(s:e,a)   = Dn(s:e,a)   - R(:,a)*coef(3)
         END DO

      ELSE IF (sstEq) THEN
!        vms_struct, FSI (vms_struct)
         IF (eq(cEq)%phys .EQ. phys_vms_struct .OR.
     2       eq(cEq)%phys .EQ. phys_FSI) THEN
            DO a=1, tnNo
               An(s:e,a)   = An(s:e,a)   - R(:,a)
               Yn(s:e,a)   = Yn(s:e,a)   - R(:,a)*coef(1)
               dUl(:)      = Rd(:,a)*coef(4) + R(1:dof-1,a)*coef(5)
               Ad(:,a)     = Ad(:,a)     - dUl(:)
               Dn(s:e-1,a) = Dn(s:e-1,a) - dUl(:)*coef(1)
            END DO

         ELSE IF (eq(cEq)%phys .EQ. phys_mesh) THEN
            DO a=1, tnNo
               An(s:e,a)   = An(s:e,a) - R(:,a)
               Yn(s:e,a)   = Yn(s:e,a) - R(:,a)*coef(1)
               Dn(s:e,a)   = Dn(s:e,a) - R(:,a)*coef(3)
            END DO

         ELSE
            err = "Undefined phys (PICC)"
         END IF

      ELSE
         DO a=1, tnNo
            An(s:e,a) = An(s:e,a) - R(:,a)
            Yn(s:e,a) = Yn(s:e,a) - R(:,a)*coef(1)
            Dn(s:e,a) = Dn(s:e,a) - R(:,a)*coef(3)
         END DO
      END IF

      IF (eq(cEq)%phys .EQ. phys_FSI) THEN
         s = eq(2)%s
         e = eq(2)%e
         DO Ac=1, tnNo
            IF (ISDOMAIN(cEq, Ac, phys_struct) .OR.
     2          ISDOMAIN(cEq, Ac, phys_vms_struct)) THEN
               An(s:e,Ac) = An(1:nsd,Ac)
               Yn(s:e,Ac) = Yn(1:nsd,Ac)
               Dn(s:e,Ac) = Dn(1:nsd,Ac)
            END IF
         END DO
      END IF

!     Update prestress at the nodes and re-initialize
      IF (eq(cEq)%phys .EQ. phys_preSt) THEN
         CALL COMMU(pSn)
         CALL COMMU(pSa)
         DO a=1, tnNo
            IF (.NOT.ISZERO(pSa(a))) THEN
               pSn(:,a) = pSn(:,a) / pSa(a)
            END IF
         END DO
         pSa = 0D0
!     Stop the simulation if the norm of displacements is less than a
!     given tolerance
         s = eq(cEq)%s
         e = eq(cEq)%e
         r1 = 0D0
         DO a=1, tnNo
            r1 = r1 + NORM(Dn(s:e,a))
         END DO
         r1 = SQRT(cm%reduce(r1))
         IF (cm%mas()) THEN
            WRITE(1000,'(A)') STR(cTS)//" "//STR(eq(cEq)%itr)//" "//
     2         STR(r1)
            CALL FLUSH(1000)
            IF(r1 .LT. 1E-6) THEN
               sCmd = "touch  "//TRIM(stopTrigName)
               CALL SYSTEM(TRIM(sCmd))
            END IF
         END IF
      END IF

!     Filter out the non-wall displacements if solving the CMM equation
      IF (eq(cEQ)%phys .EQ. phys_CMM) CALL CMM_DISPF(eq(cEq), Dn)

      IF (ISZERO(eq(cEq)%FSILS%RI%iNorm)) eq(cEq)%FSILS%RI%iNorm = eps
      IF (ISZERO(eq(cEq)%iNorm)) eq(cEq)%iNorm = eq(cEq)%FSILS%RI%iNorm
      IF (eq(cEq)%itr .EQ. 1) THEN
         eq(cEq)%pNorm = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm
      END IF
      r1 = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm
      r2 = 2D1*LOG10(r1/eq(cEq)%pNorm)

      l1 = eq(cEq)%itr .GE. eq(cEq)%maxItr
      l2 = r1 .LE. eq(cEq)%tol
      l3 = r1 .LE. eq(cEq)%tol*eq(cEq)%pNorm
      l4 = eq(cEq)%itr .GE. eq(cEq)%minItr
      l5 = r2 .LE. eq(cEq)%dBr
      IF (l1 .OR. ((l2.OR.l3).AND.l4.AND.l5)) eq(cEq)%ok = .TRUE.
      IF (ALL(eq%ok)) RETURN

      IF (eq(cEq)%coupled) THEN
         cEq = cEq + 1
         IF (ALL(.NOT.eq%coupled .OR. eq%ok)) THEN
            DO WHILE (cEq .LE. nEq)
               IF (.NOT.eq(cEq)%coupled) EXIT
               cEq = cEq + 1
            END DO
         ELSE
            IF (cEq .GT. nEq) cEq = 1
            DO WHILE (.NOT.eq(cEq)%coupled)
               cEq = cEq + 1
               IF (cEq .GT. nEq) cEq = 1
            END DO
         END IF
      ELSE
         IF (eq(cEq)%ok) cEq = cEq + 1
      END IF

      RETURN
      END SUBROUTINE PICC
!====================================================================


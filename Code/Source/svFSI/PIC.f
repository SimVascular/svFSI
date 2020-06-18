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

      INTEGER(KIND=IKIND) iEq, s, e
      REAL(KIND=RKIND) coef, ctime

!     Prestress initialization
      IF (pstEq) THEN
         pS0 = pS0 + pSn
         Ao = 0._RKIND
         Yo = 0._RKIND
         Do = 0._RKIND
      END IF

!     Immersed body treatment (explicit)
      IF (ibFlag) THEN
         ib%callD(1) = CPUT()

!        Set IB forces to zero, except for feedback force
         ib%R = 0._RKIND

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
         coef = (eq(iEq)%gam - 1._RKIND)/eq(iEq)%gam
         An(s:e,:) = Ao(s:e,:)*coef

!        electrophysiology
         IF (eq(iEq)%phys .EQ. phys_CEP) THEN
            CALL CEPINTEG(iEq, e, Do)
         END IF

         Yn(s:e,:) = Yo(s:e,:)

         IF (dFlag) THEN
            IF (.NOT.sstEq) THEN
!              struct, lElas, FSI (struct, mesh)
               coef = dt*dt*(0.5_RKIND*eq(iEq)%gam - eq(iEq)%beta)
     2            /(eq(iEq)%gam - 1._RKIND)
               Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
            ELSE
!              ustruct, FSI
               IF (eq(iEq)%phys .EQ. phys_ustruct .OR.
     2             eq(iEq)%phys .EQ. phys_FSI) THEN
                  coef = (eq(iEq)%gam - 1._RKIND)/eq(iEq)%gam
                  Ad(:,:)   = Ad(:,:)*coef
                  Dn(s:e,:) = Do(s:e,:)
               ELSE IF (eq(iEq)%phys .EQ. phys_mesh) THEN
!              mesh
                  coef = dt*dt*(0.5_RKIND*eq(iEq)%gam - eq(iEq)%beta)
     2               /(eq(iEq)%gam - 1._RKIND)
                  Dn(s:e,:) = Do(s:e,:) + Yn(s:e,:)*dt + An(s:e,:)*coef
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
      REAL(KIND=RKIND), INTENT(INOUT) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) s, e, i, a
      REAL(KIND=RKIND) coef(4)

      dof         = eq(cEq)%dof
      eq(cEq)%itr = eq(cEq)%itr + 1

      DO i=1, nEq
         s       = eq(i)%s
         e       = eq(i)%e
         coef(1) = 1._RKIND - eq(i)%am
         coef(2) = eq(i)%am
         coef(3) = 1._RKIND - eq(i)%af
         coef(4) = eq(i)%af

         DO a=1, tnNo
            Ag(s:e,a) = Ao(s:e,a)*coef(1) + An(s:e,a)*coef(2)
            Yg(s:e,a) = Yo(s:e,a)*coef(3) + Yn(s:e,a)*coef(4)
            Dg(s:e,a) = Do(s:e,a)*coef(3) + Dn(s:e,a)*coef(4)
         END DO

!        Reset pressure variable initiator
         IF ( (eq(i)%phys .EQ. phys_fluid .OR.
     2         eq(i)%phys .EQ. phys_CMM   .OR.
     3         eq(i)%phys .EQ. phys_FSI) .AND.
     4        .NOT.sstEq .AND. .NOT.cmmInit ) THEN
            DO a=1, tnNo
               Yg(e,a) = Yn(e,a)
            END DO
         END IF
      END DO

      IF (pstEq) THEN
         pSn(:,:) = 0._RKIND
         pSa(:)   = 0._RKIND
      END IF

      RETURN
      END SUBROUTINE PICI
!====================================================================
!     This is the corrector. Decision for next eqn is also made here
      SUBROUTINE PICC
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: l1, l2, l3, l4
      INTEGER(KIND=IKIND) :: s, e, a, Ac
      REAL(KIND=RKIND) :: coef(5), r1, dUl(nsd)

      s       = eq(cEq)%s
      e       = eq(cEq)%e
      coef(1) = eq(cEq)%gam*dt
      coef(2) = eq(cEq)%af*coef(1)
      coef(3) = eq(cEq)%beta*dt*dt
      coef(4) = 1._RKIND / eq(cEq)%am
      coef(5) = coef(2)*coef(4)

      IF ( (eq(cEq)%phys .EQ. phys_fluid .OR.
     2      eq(cEq)%phys .EQ. phys_CMM   .OR.
     3      eq(cEq)%phys .EQ. phys_FSI) .AND.
     4     .NOT.sstEq .AND. .NOT.cmmInit ) THEN
         DO a=1, tnNo
            An(s:e,a)   = An(s:e,a)   - R(:,a)
            Yn(s:e-1,a) = Yn(s:e-1,a) - R(1:dof-1,a)*coef(1)
            Yn(e,a)     = Yn(e,a)     - R(dof,a)*coef(2)
            Dn(s:e,a)   = Dn(s:e,a)   - R(:,a)*coef(3)
         END DO
      ELSE IF (sstEq) THEN
!        ustruct, FSI (ustruct)
         IF (eq(cEq)%phys .EQ. phys_ustruct .OR.
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
         END IF
      ELSE
         DO a=1, tnNo
            An(s:e,a) = An(s:e,a) - R(:,a)
            Yn(s:e,a) = Yn(s:e,a) - R(:,a)*coef(1)
            Dn(s:e,a) = Dn(s:e,a) - R(:,a)*coef(3)
         END DO
      END IF

      IF ((eq(cEq)%phys .EQ. phys_ustruct) .OR.
     2    (eq(cEq)%phys .EQ. phys_stokes)) THEN
         CALL PICETH()
      END IF

      IF (eq(cEq)%phys .EQ. phys_FSI) THEN
         s = eq(2)%s
         e = eq(2)%e
         DO Ac=1, tnNo
            IF (ISDOMAIN(cEq, Ac, phys_struct) .OR.
     2          ISDOMAIN(cEq, Ac, phys_ustruct) .OR.
     3          ISDOMAIN(cEq, Ac, phys_lElas)) THEN
               An(s:e,Ac) = An(1:nsd,Ac)
               Yn(s:e,Ac) = Yn(1:nsd,Ac)
               Dn(s:e,Ac) = Dn(1:nsd,Ac)
            END IF
         END DO
      END IF

!     Update Xion for cardiac electrophysiology
      IF (eq(cEq)%phys .EQ. phys_CEP) THEN
         s = eq(cEq)%s
         DO a=1, tnNo
            Xion(1,a) = Yn(s,a)
         END DO
      END IF

!     Update prestress at the nodes and re-initialize
      IF (pstEq) THEN
         CALL COMMU(pSn)
         CALL COMMU(pSa)
         DO a=1, tnNo
            IF (.NOT.ISZERO(pSa(a))) THEN
               pSn(:,a) = pSn(:,a) / pSa(a)
            END IF
         END DO
         pSa = 0._RKIND
      END IF

!     Filter out the non-wall displacements for CMM equation
      IF (eq(cEq)%phys.EQ.phys_CMM .AND. .NOT.cmmInit) THEN
         DO a=1, tnNo
            r1 = REAL(cmmBdry(a), KIND=RKIND)
            Dn(s:e-1,a) = Dn(s:e-1,a)*r1
         END DO
      END IF

      IF (ISZERO(eq(cEq)%FSILS%RI%iNorm)) eq(cEq)%FSILS%RI%iNorm = eps
      IF (ISZERO(eq(cEq)%iNorm)) eq(cEq)%iNorm = eq(cEq)%FSILS%RI%iNorm
      IF (eq(cEq)%itr .EQ. 1) THEN
         eq(cEq)%pNorm = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm
      END IF
      r1 = eq(cEq)%FSILS%RI%iNorm/eq(cEq)%iNorm

      l1 = eq(cEq)%itr .GE. eq(cEq)%maxItr
      l2 = r1 .LE. eq(cEq)%tol
      l3 = r1 .LE. eq(cEq)%tol*eq(cEq)%pNorm
      l4 = eq(cEq)%itr .GE. eq(cEq)%minItr
      IF (l1 .OR. ((l2.OR.l3).AND.l4)) eq(cEq)%ok = .TRUE.
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
!     Pressure correction at edge nodes for Taylor-Hood type element
!     via interpolation
      SUBROUTINE PICETH()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL THflag
      INTEGER(KIND=IKIND) a, b, e, g, s, iM, Ac, eType, eNoN, eNoNq
      REAL(KIND=RKIND) Jac, eVol, p, xp(nsd), xi0(nsd), xi(nsd),
     2   ksix(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), xql(:,:), pl(:), Nq(:),
     2   Nqx(:,:), sA(:), sF(:)

      THflag = .FALSE.
      DO iM=1, nMsh
         IF (msh(iM)%nFs .EQ. 2) THEN
            THflag = .TRUE.
            EXIT
         END IF
      END DO
      IF (.NOT.THflag) RETURN

      ALLOCATE(sA(tnNo), sF(tnNo))
      sF(:) = 0._RKIND
      sA(:) = 0._RKIND

      s = eq(cEq)%s
      DO iM=1, nMsh
         IF (msh(iM)%nFs .EQ. 1) CYCLE

         eType = msh(iM)%fs(2)%eType

         eNoN  = msh(iM)%fs(1)%eNoN
         eNoNq = msh(iM)%fs(2)%eNoN
         ALLOCATE(xl(nsd,eNoN), xql(nsd,eNoNq), pl(eNoNq), Nq(eNoNq),
     2      Nqx(nsd,eNoNq))

         xi0 = 0._RKIND
         DO g=1, msh(iM)%fs(2)%nG
            xi0 = xi0 + msh(iM)%fs(2)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(msh(iM)%fs(2)%nG, KIND=RKIND)

         DO e=1, msh(iM)%nEl
            cDmn = DOMAIN(msh(iM), cEq, e)
            IF ((eq(cEq)%dmn(cDmn)%phys .NE. phys_ustruct) .AND.
     2          (eq(cEq)%dmn(cDmn)%phys .NE. phys_stokes)) CYCLE

            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,e)
               xl(:,a) = x(:,Ac)
            END DO

            DO a=1, eNoNq
               Ac = msh(iM)%IEN(a,e)
               pl(a)    = Yn(s+nsd,Ac)
               xql(:,a) = xl(:,a)
            END DO

            eVol = 0._RKIND
            DO g=1, msh(iM)%fs(2)%nG
               IF (g.EQ.1 .OR. .NOT.msh(iM)%fs(2)%lShpF) THEN
                  CALL GNN(eNoNq, nsd, msh(iM)%fs(2)%Nx(:,:,g), xql,
     2               Nqx, Jac, ksix)
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               END IF
               eVol = eVol + msh(iM)%fs(2)%w(g)*Jac
            END DO

            DO a=eNoNq+1, eNoN
               Ac = msh(iM)%IEN(a,e)
               xp = xl(:,a)

               xi = xi0
               CALL GETNNX(eType, eNoNq, xql, msh(iM)%fs(2)%xib,
     2            msh(iM)%fs(2)%Nb, xp, xi, Nq, Nqx)

               p = 0._RKIND
               DO b=1, eNoNq
                  p = p + pl(b)*Nq(b)
               END DO

               sF(Ac) = sF(Ac) + p*eVol
               sA(Ac) = sA(Ac) + eVol
            END DO

         END DO ! e-loop
         DEALLOCATE(xl, xql, pl, Nq, Nqx)
      END DO ! iM-loop

      CALL COMMU(sA)
      CALL COMMU(sF)

      DO a=1, tnNo
         IF (.NOT.ISZERO(sA(a))) Yn(s+nsd,a) = sF(a)/sA(a)
      END DO

      DEALLOCATE(sA, sF)

      RETURN
      END SUBROUTINE PICETH
!====================================================================


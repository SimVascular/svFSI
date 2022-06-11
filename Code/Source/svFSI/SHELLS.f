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
!     This routines is for solving nonlinear shell mechanics problem
!     using finite elements and IGA.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_SHELL(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, b, e, g, Ac, eNoN, cPhys

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN
      IF (lM%eType .EQ. eType_TRI3) eNoN = 2*eNoN

!     Initialize tensor operations
      CALL TEN_INIT(2)

!     SHELLS: dof = nsd
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), lR(dof,eNoN),
     3   lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_shell) CYCLE

!        Create local copies
         xl  = 0._RKIND
         al  = 0._RKIND
         yl  = 0._RKIND
         dl  = 0._RKIND
         bfl = 0._RKIND
         DO a=1, eNoN
            IF (a .LE. lM%eNoN) THEN
               Ac = lM%IEN(a,e)
               ptr(a) = Ac
            ELSE
               b  = a - lM%eNoN
               Ac = lM%eIEN(b,e)
               ptr(a) = Ac
               IF (Ac .EQ. 0) CYCLE
            END IF
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
         END DO

         IF (lM%eType .EQ. eType_TRI3) THEN
!           Constant strain triangles, no numerical integration
            CALL SHELLCST(lM, e, eNoN, al, yl, dl, xl, bfl, ptr)

         ELSE
            lR = 0._RKIND
            lK = 0._RKIND

!           Update shape functions for NURBS elements
            IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!           Gauss integration
            DO g=1, lM%nG
               CALL SHELL3D(lM, g, eNoN, al, yl, dl, xl, bfl, lR, lK)
            END DO

!           Assembly
#ifdef WITH_TRILINOS
            IF (eq(cEq)%assmTLS) THEN
               CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
            ELSE
#endif
               CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
            END IF
#endif
         END IF
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, dl, bfl, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_SHELL
!####################################################################
!     Construct shell mechanics for constant strain triangle elements
      SUBROUTINE SHELLCST (lM, e, eNoN, al, yl, dl, xl, bfl, ptr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e, eNoN, ptr(eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), xl(3,eNoN), bfl(3,eNoN)

      LOGICAL :: setIt(3)
      INTEGER(KIND=IKIND) :: i, j, k, a, b, g
      REAL(KIND=RKIND) :: amd, afl, rho, dmp, ht, w, Jac0, Jac, fb(3),
     2   ud(3), nV0(3), nV(3), x0(3,eNoN), xc(3,eNoN), aCov(3,2),
     3   aCov0(3,2), aCnv(3,2), aCnv0(3,2), aa_0(2,2), aa_x(2,2),
     4   bb_0(2,2), bb_x(2,2), Sm(3,2), Dm(3,3,3), Bm(3,3,lM%eNoN),
     5   Bb(3,3,eNoN), D0Bm(3,3,lM%eNoN), D1Bm(3,3,lM%eNoN),
     6   D1Bb(3,3,eNoN), D2Bb(3,3,eNoN), BtS, NxSNx, BtDB

      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:),
     2   lK(:,:,:), tmpX(:,:)

!     Note that for triangular elements, eNoN=6 and lM%eNoN=3
      ALLOCATE(N(lM%eNoN), Nx(2,lM%eNoN), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN))

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp   = eq(cEq)%dmn(cDmn)%prop(damping)
      ht    = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      amd  = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl  = eq(cEq)%af*eq(cEq)%beta*dt*dt

      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

!---------------------------------------------------------------------
!     Get the reference configuration
      x0(:,:) = xl(:,:)

!     Get the current configuration
      DO a=1, eNoN
         xc(1,a) = x0(1,a) + dl(i,a)
         xc(2,a) = x0(2,a) + dl(j,a)
         xc(3,a) = x0(3,a) + dl(k,a)
      END DO
      Nx(:,:) = lM%Nx(:,:,1)

!---------------------------------------------------------------------
!     Covariant and contravariant bases in reference config
      ALLOCATE(tmpX(nsd,lM%eNoN))
      tmpX = x0(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, Nx, tmpX, nV0, aCov0, aCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      tmpX = xc(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, Nx, tmpX, nV, aCov, aCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac
      DEALLOCATE(tmpX)

!---------------------------------------------------------------------
!     Compute metric tensor in reference and current config
      aa_0 = 0._RKIND
      aa_x = 0._RKIND
      DO g=1, nsd
         aa_0(1,1) = aa_0(1,1) + aCov0(g,1)*aCov0(g,1)
         aa_0(1,2) = aa_0(1,2) + aCov0(g,1)*aCov0(g,2)
         aa_0(2,1) = aa_0(2,1) + aCov0(g,2)*aCov0(g,1)
         aa_0(2,2) = aa_0(2,2) + aCov0(g,2)*aCov0(g,2)

         aa_x(1,1) = aa_x(1,1) + aCov(g,1)*aCov(g,1)
         aa_x(1,2) = aa_x(1,2) + aCov(g,1)*aCov(g,2)
         aa_x(2,1) = aa_x(2,1) + aCov(g,2)*aCov(g,1)
         aa_x(2,2) = aa_x(2,2) + aCov(g,2)*aCov(g,2)
      END DO

!---------------------------------------------------------------------
!     Define variation in membrane strain only for the main element
      Bm = 0._RKIND
      DO a=1, lM%eNoN
         Bm(1,1,a) = Nx(1,a)*aCov(1,1)
         Bm(1,2,a) = Nx(1,a)*aCov(2,1)
         Bm(1,3,a) = Nx(1,a)*aCov(3,1)

         Bm(2,1,a) = Nx(2,a)*aCov(1,2)
         Bm(2,2,a) = Nx(2,a)*aCov(2,2)
         Bm(2,3,a) = Nx(2,a)*aCov(3,2)

         Bm(3,1,a) = Nx(2,a)*aCov(1,1) + Nx(1,a)*aCov(1,2)
         Bm(3,2,a) = Nx(2,a)*aCov(2,1) + Nx(1,a)*aCov(2,2)
         Bm(3,3,a) = Nx(2,a)*aCov(3,1) + Nx(1,a)*aCov(3,2)
      END DO

!     For the boundary elements, zero-out Bm for fixed/clamped BC.
      setIt = .FALSE.
      a = lM%eNoN + 1
      DO WHILE (a .LE. eNoN)
         IF (ptr(a) .EQ. 0) THEN
            b = a - lM%eNoN
            IF (BTEST(lM%sbc(b,e),bType_fix)) setIt(b) = .TRUE.
         END IF
         a = a + 1
      END DO

      DO a=1, lM%eNoN
         IF (setIt(a)) THEN
            DO b=1, lM%eNoN
               IF (a .EQ. b) CYCLE
               Bm(:,:,b) = 0._RKIND
            END DO
         END IF
      END DO

!---------------------------------------------------------------------
!     Compute curvature coefficients for bending strain and its
!     variation for CST elements
      CALL SHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x, Bb)

!---------------------------------------------------------------------
!     Compute stress resultants by integrating 2nd Piola Kirchhoff
!     stress and elasticity tensors through the shell thickness. These
!     resultants are computed in Voigt notation.
      CALL SHL_STRS_RES(eq(cEq)%dmn(cDmn), aa_0, aa_x, bb_0, bb_x, Sm,
     2   Dm)

!---------------------------------------------------------------------
!     Contribution to tangent matrices: Dm * Bm, Dm*Bb
      DO a=1, lM%eNoN
         D0Bm(:,:,a) = MATMUL(Dm(:,:,1), Bm(:,:,a))
         D1Bm(:,:,a) = MATMUL(Dm(:,:,2), Bm(:,:,a))
      END DO

      DO a=1, eNoN
         D1Bb(:,:,a) = MATMUL(Dm(:,:,2), Bb(:,:,a))
         D2Bb(:,:,a) = MATMUL(Dm(:,:,3), Bb(:,:,a))
      END DO

!---------------------------------------------------------------------
!     Contribution to residue and stiffness matrices due to inertia and
!     body forces
      lR = 0._RKIND
      lK = 0._RKIND
      DO g=1, lM%nG
         N = lM%N(:,g)
         w = lM%w(g)*Jac0*ht

!        Acceleration and mass damping at the integration point
         ud = -fb
         DO a=1, lM%eNoN
            ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
            ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))
            ud(3) = ud(3) + N(a)*(rho*(al(k,a)-bfl(3,a)) + dmp*yl(k,a))
         END DO

!        Local residue
         DO a=1, lM%eNoN
            lR(1,a) = lR(1,a) + N(a)*w*ud(1)
            lR(2,a) = lR(2,a) + N(a)*w*ud(2)
            lR(3,a) = lR(3,a) + N(a)*w*ud(3)
         END DO

!        Local stiffness contribution from mass matrix
         DO b=1, lM%eNoN
            DO a=1, lM%eNoN
               BtS = w*amd*N(a)*N(b)
               lK(1,a,b) = lK(1,a,b) + BtS
               lK(dof+2,a,b) = lK(dof+2,a,b) + BtS
               lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + BtS
            END DO
         END DO
      END DO

!---------------------------------------------------------------------
!     Contribution to residue from membrane strain
      w = Jac0 * 0.5_RKIND
      DO a=1, lM%eNoN
         BtS = Bm(1,1,a)*Sm(1,1) + Bm(2,1,a)*Sm(2,1) + Bm(3,1,a)*Sm(3,1)
         lR(1,a) = lR(1,a) + w*BtS

         BtS = Bm(1,2,a)*Sm(1,1) + Bm(2,2,a)*Sm(2,1) + Bm(3,2,a)*Sm(3,1)
         lR(2,a) = lR(2,a) + w*BtS

         BtS = Bm(1,3,a)*Sm(1,1) + Bm(2,3,a)*Sm(2,1) + Bm(3,3,a)*Sm(3,1)
         lR(3,a) = lR(3,a) + w*BtS
      END DO

!     Contribution to residue from bending strain
      DO a=1, eNoN
         BtS = Bb(1,1,a)*Sm(1,2) + Bb(2,1,a)*Sm(2,2) + Bb(3,1,a)*Sm(3,2)
         lR(1,a) = lR(1,a) + w*BtS

         BtS = Bb(1,2,a)*Sm(1,2) + Bb(2,2,a)*Sm(2,2) + Bb(3,2,a)*Sm(3,2)
         lR(2,a) = lR(2,a) + w*BtS

         BtS = Bb(1,3,a)*Sm(1,2) + Bb(2,3,a)*Sm(2,2) + Bb(3,3,a)*Sm(3,2)
         lR(3,a) = lR(3,a) + w*BtS
      END DO

!---------------------------------------------------------------------
!     Contribution to stiffness from membrane-membrane interactions
      w = afl*Jac0*0.5_RKIND
      DO b=1, lM%eNoN
         DO a=1, lM%eNoN
            NxSNx = Nx(1,a)*Nx(1,b)*Sm(1,1) + Nx(2,a)*Nx(2,b)*Sm(2,1)
     2            + Nx(1,a)*Nx(2,b)*Sm(3,1) + Nx(2,a)*Nx(1,b)*Sm(3,1)

            BtDB = Bm(1,1,a)*D0Bm(1,1,b) + Bm(2,1,a)*D0Bm(2,1,b)
     2           + Bm(3,1,a)*D0Bm(3,1,b)
            lK(1,a,b) = lK(1,a,b) + w*(BtDB + NxSNx)

            BtDB = Bm(1,1,a)*D0Bm(1,2,b) + Bm(2,1,a)*D0Bm(2,2,b)
     2           + Bm(3,1,a)*D0Bm(3,2,b)
            lK(2,a,b) = lK(2,a,b) + w*BtDB

            BtDB = Bm(1,1,a)*D0Bm(1,3,b) + Bm(2,1,a)*D0Bm(2,3,b)
     2           + Bm(3,1,a)*D0Bm(3,3,b)
            lK(3,a,b) = lK(3,a,b) + w*BtDB

            BtDB = Bm(1,2,a)*D0Bm(1,1,b) + Bm(2,2,a)*D0Bm(2,1,b)
     2           + Bm(3,2,a)*D0Bm(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB

            BtDB = Bm(1,2,a)*D0Bm(1,2,b) + Bm(2,2,a)*D0Bm(2,2,b)
     2           + Bm(3,2,a)*D0Bm(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*(BtDB + NxSNx)

            BtDB = Bm(1,2,a)*D0Bm(1,3,b) + Bm(2,2,a)*D0Bm(2,3,b)
     2           + Bm(3,2,a)*D0Bm(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D0Bm(1,1,b) + Bm(2,3,a)*D0Bm(2,1,b)
     2           + Bm(3,3,a)*D0Bm(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D0Bm(1,2,b) + Bm(2,3,a)*D0Bm(2,2,b)
     2           + Bm(3,3,a)*D0Bm(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D0Bm(1,3,b) + Bm(2,3,a)*D0Bm(2,3,b)
     2           + Bm(3,3,a)*D0Bm(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + w*(BtDB + NxSNx)
         END DO
      END DO

!     Contribution to stiffness from membrane-bending interactions
      DO b=1, eNoN
         DO a=1, lM%eNoN
            BtDB = Bm(1,1,a)*D1Bb(1,1,b) + Bm(2,1,a)*D1Bb(2,1,b)
     2           + Bm(3,1,a)*D1Bb(3,1,b)
            lK(1,a,b) = lK(1,a,b) + w*BtDB

            BtDB = Bm(1,1,a)*D1Bb(1,2,b) + Bm(2,1,a)*D1Bb(2,2,b)
     2           + Bm(3,1,a)*D1Bb(3,2,b)
            lK(2,a,b) = lK(2,a,b) + w*BtDB

            BtDB = Bm(1,1,a)*D1Bb(1,3,b) + Bm(2,1,a)*D1Bb(2,3,b)
     2           + Bm(3,1,a)*D1Bb(3,3,b)
            lK(3,a,b) = lK(3,a,b) + w*BtDB

            BtDB = Bm(1,2,a)*D1Bb(1,1,b) + Bm(2,2,a)*D1Bb(2,1,b)
     2           + Bm(3,2,a)*D1Bb(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB

            BtDB = Bm(1,2,a)*D1Bb(1,2,b) + Bm(2,2,a)*D1Bb(2,2,b)
     2           + Bm(3,2,a)*D1Bb(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB

            BtDB = Bm(1,2,a)*D1Bb(1,3,b) + Bm(2,2,a)*D1Bb(2,3,b)
     2           + Bm(3,2,a)*D1Bb(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D1Bb(1,1,b) + Bm(2,3,a)*D1Bb(2,1,b)
     2           + Bm(3,3,a)*D1Bb(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D1Bb(1,2,b) + Bm(2,3,a)*D1Bb(2,2,b)
     2           + Bm(3,3,a)*D1Bb(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB

            BtDB = Bm(1,3,a)*D1Bb(1,3,b) + Bm(2,3,a)*D1Bb(2,3,b)
     2           + Bm(3,3,a)*D1Bb(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + w*BtDB
         END DO
      END DO

!     Contribution to stiffness from bending-membrane interactions
      DO b=1, lM%eNoN
         DO a=1, eNoN
            BtDB = Bb(1,1,a)*D1Bm(1,1,b) + Bb(2,1,a)*D1Bm(2,1,b)
     2           + Bb(3,1,a)*D1Bm(3,1,b)
            lK(1,a,b) = lK(1,a,b) + w*BtDB

            BtDB = Bb(1,1,a)*D1Bm(1,2,b) + Bb(2,1,a)*D1Bm(2,2,b)
     2           + Bb(3,1,a)*D1Bm(3,2,b)
            lK(2,a,b) = lK(2,a,b) + w*BtDB

            BtDB = Bb(1,1,a)*D1Bm(1,3,b) + Bb(2,1,a)*D1Bm(2,3,b)
     2           + Bb(3,1,a)*D1Bm(3,3,b)
            lK(3,a,b) = lK(3,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D1Bm(1,1,b) + Bb(2,2,a)*D1Bm(2,1,b)
     2           + Bb(3,2,a)*D1Bm(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D1Bm(1,2,b) + Bb(2,2,a)*D1Bm(2,2,b)
     2           + Bb(3,2,a)*D1Bm(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D1Bm(1,3,b) + Bb(2,2,a)*D1Bm(2,3,b)
     2           + Bb(3,2,a)*D1Bm(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D1Bm(1,1,b) + Bb(2,3,a)*D1Bm(2,1,b)
     2           + Bb(3,3,a)*D1Bm(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D1Bm(1,2,b) + Bb(2,3,a)*D1Bm(2,2,b)
     2           + Bb(3,3,a)*D1Bm(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D1Bm(1,3,b) + Bb(2,3,a)*D1Bm(2,3,b)
     2           + Bb(3,3,a)*D1Bm(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + w*BtDB
         END DO
      END DO

!     Contribution to stiffness from bending-bending interactions
      DO b=1, eNoN
         DO a=1, eNoN
            BtDB = Bb(1,1,a)*D2Bb(1,1,b) + Bb(2,1,a)*D2Bb(2,1,b)
     2           + Bb(3,1,a)*D2Bb(3,1,b)
            lK(1,a,b) = lK(1,a,b) + w*BtDB

            BtDB = Bb(1,1,a)*D2Bb(1,2,b) + Bb(2,1,a)*D2Bb(2,2,b)
     2           + Bb(3,1,a)*D2Bb(3,2,b)
            lK(2,a,b) = lK(2,a,b) + w*BtDB

            BtDB = Bb(1,1,a)*D2Bb(1,3,b) + Bb(2,1,a)*D2Bb(2,3,b)
     2           + Bb(3,1,a)*D2Bb(3,3,b)
            lK(3,a,b) = lK(3,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D2Bb(1,1,b) + Bb(2,2,a)*D2Bb(2,1,b)
     2           + Bb(3,2,a)*D2Bb(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D2Bb(1,2,b) + Bb(2,2,a)*D2Bb(2,2,b)
     2           + Bb(3,2,a)*D2Bb(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*BtDB

            BtDB = Bb(1,2,a)*D2Bb(1,3,b) + Bb(2,2,a)*D2Bb(2,3,b)
     2           + Bb(3,2,a)*D2Bb(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D2Bb(1,1,b) + Bb(2,3,a)*D2Bb(2,1,b)
     2           + Bb(3,3,a)*D2Bb(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D2Bb(1,2,b) + Bb(2,3,a)*D2Bb(2,2,b)
     2           + Bb(3,3,a)*D2Bb(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*BtDB

            BtDB = Bb(1,3,a)*D2Bb(1,3,b) + Bb(2,3,a)*D2Bb(2,3,b)
     2           + Bb(3,3,a)*D2Bb(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + w*BtDB
         END DO
      END DO

!---------------------------------------------------------------------
!     Global assembly
#ifdef WITH_TRILINOS
      IF (eq(cEq)%assmTLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#endif
         CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
      END IF
#endif

      DEALLOCATE(N, Nx, lR, lK)

      RETURN
      END SUBROUTINE SHELLCST
!--------------------------------------------------------------------
!     This subroutine computes bending strain, Eb, and its variational
!     derivative, Bb, for CST elements
      SUBROUTINE SHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x, Bb)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e, ptr(6)
      REAL(KIND=RKIND), INTENT(INOUT) ::  x0(3,6), xc(3,6)
      REAL(KIND=RKIND), INTENT(OUT) :: bb_0(2,2), bb_x(2,2), Bb(3,3,6)

      LOGICAL :: bFlag, lFix(3)
      INTEGER(KIND=IKIND) :: i, j, p, f, eNoN

      REAL(KIND=RKIND) :: Jac0, Jac, cI, aIi, nV0(3), nV(3), eI(3),
     2   nI(3), eIeI(3,3), nInI(3,3), eIaP(3,3), Im(3,3), aCov(3,2),
     3   aCnv(3,2), aCov0(3,2), aCnv0(3,2)

      REAL(KIND=RKIND) :: a0(3,6), a(3,6), adg0(3,3), adg(3,3),
     2   xi0(3,3), xi(3,3), Tm0(3,3), Tm(3,3), v0(3), v(3), B1b(3,6),
     3   Nm(3,3), Mm(3,3,2), H1(6,18), H2(18,18), H3(3,18), H1H2(6,18),
     4   H3H2(3,18), Bb1(3,18), tmpA(3,3)

      eNoN = 2*lM%eNoN
!     Boundary element check
      bFlag = .FALSE.
      DO j=lM%eNoN+1, eNoN
         IF (ptr(j) .EQ. 0) THEN
            bFlag = .TRUE.
            EXIT
         END IF
      END DO

!     Edge vectors of the main element (reference config)
      a0(:,1) = x0(:,3) - x0(:,2)
      a0(:,2) = x0(:,1) - x0(:,3)
      a0(:,3) = x0(:,2) - x0(:,1)

!     Edge vectors of the main element (current config)
      a(:,1)  = xc(:,3) - xc(:,2)
      a(:,2)  = xc(:,1) - xc(:,3)
      a(:,3)  = xc(:,2) - xc(:,1)

!     Covariant and contravariant bases in reference config
      tmpA = x0(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV0, aCov0, aCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      tmpA = xc(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV, aCov, aCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac

!     Update the position vector of the `artificial' or `ghost' nodes
!     depending on the boundary condition.
      IF (bFlag) THEN
         DO j=lM%eNoN+1, eNoN
            IF (ptr(j) .NE. 0) CYCLE
            i = j - lM%eNoN
            p = i - 1
            IF (i .EQ. 1) p = 3
!           Reference config
!           eI = eI0 = aI0/|aI0| (reference config)
            aIi   = 1._RKIND/SQRT(NORM(a0(:,i)))
            eI(:) = a0(:,i) * aIi
!           nI = nI0 = eI0 x n0 (reference config)
            nI(1) = eI(2)*nV0(3) - eI(3)*nV0(2)
            nI(2) = eI(3)*nV0(1) - eI(1)*nV0(3)
            nI(3) = eI(1)*nV0(2) - eI(2)*nV0(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            x0(1,j) = 2._RKIND*(nInI(1,1)*a0(1,p) + nInI(1,2)*a0(2,p) +
     2         nInI(1,3)*a0(3,p)) + x0(1,i)
            x0(2,j) = 2._RKIND*(nInI(2,1)*a0(1,p) + nInI(2,2)*a0(2,p) +
     2         nInI(2,3)*a0(3,p)) + x0(2,i)
            x0(3,j) = 2._RKIND*(nInI(3,1)*a0(1,p) + nInI(3,2)*a0(2,p) +
     2         nInI(3,3)*a0(3,p)) + x0(3,i)

!           Current config
!           eI = aI/|aI| (current config)
            aIi   = 1._RKIND/SQRT(NORM(a(:,i)))
            eI(:) = a(:,i)*aIi
!           nI = eI x n (currnt config)
            nI(1) = eI(2)*nV(3) - eI(3)*nV(2)
            nI(2) = eI(3)*nV(1) - eI(1)*nV(3)
            nI(3) = eI(1)*nV(2) - eI(2)*nV(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            xc(1,j) = 2._RKIND*(nInI(1,1)*a(1,p) + nInI(1,2)*a(2,p) +
     2         nInI(1,3)*a(3,p)) + xc(1,i)
            xc(2,j) = 2._RKIND*(nInI(2,1)*a(1,p) + nInI(2,2)*a(2,p) +
     2         nInI(2,3)*a(3,p)) + xc(2,i)
            xc(3,j) = 2._RKIND*(nInI(3,1)*a(1,p) + nInI(3,2)*a(2,p) +
     2         nInI(3,3)*a(3,p)) + xc(3,i)

            IF (BTEST(lM%sbc(i,e),bType_fix)) xc(:,j) = x0(:,j)
         END DO
      END IF

!     Edge vector of surrounding nodes (reference config)
      a0(:,4)   = x0(:,4) - x0(:,2)
      a0(:,5)   = x0(:,5) - x0(:,3)
      a0(:,6)   = x0(:,6) - x0(:,1)

!     Edge vector of surrounding nodes (current config)
      a(:,4)    = xc(:,4) - xc(:,2)
      a(:,5)    = xc(:,5) - xc(:,3)
      a(:,6)    = xc(:,6) - xc(:,1)

!     a.gCnv (reference config)
      adg0(1,1) = NORM(a0(:,4), aCnv0(:,1)) ! xi_4
      adg0(1,2) = NORM(a0(:,5), aCnv0(:,1)) ! xi_5
      adg0(1,3) = NORM(a0(:,6), aCnv0(:,1)) ! xi_6

      adg0(2,1) = NORM(a0(:,4), aCnv0(:,2)) ! eta_4
      adg0(2,2) = NORM(a0(:,5), aCnv0(:,2)) ! eta_5
      adg0(2,3) = NORM(a0(:,6), aCnv0(:,2)) ! eta_6

      adg0(3,1) = NORM(a0(:,4), nV0)        ! z_4
      adg0(3,2) = NORM(a0(:,5), nV0)        ! z_5
      adg0(3,3) = NORM(a0(:,6), nV0)        ! z_6

!     a.gCnv (current config)
      adg(1,1)  = NORM(a(:,4), aCnv(:,1)) ! xi_4
      adg(1,2)  = NORM(a(:,5), aCnv(:,1)) ! xi_5
      adg(1,3)  = NORM(a(:,6), aCnv(:,1)) ! xi_6

      adg(2,1)  = NORM(a(:,4), aCnv(:,2)) ! eta_4
      adg(2,2)  = NORM(a(:,5), aCnv(:,2)) ! eta_5
      adg(2,3)  = NORM(a(:,6), aCnv(:,2)) ! eta_6

      adg(3,1)  = NORM(a(:,4), nV)        ! z_4
      adg(3,2)  = NORM(a(:,5), nV)        ! z_5
      adg(3,3)  = NORM(a(:,6), nV)        ! z_6

!     Xi matrix (reference config)
      xi0(:,:)  = adg0(:,:)
      xi0(1,3)  = adg0(1,3) + 1._RKIND    ! xi_6
      xi0(2,1)  = adg0(2,1) + 1._RKIND    ! eta_4

!     Xi matrix (current config)
      xi(:,:)   = adg(:,:)
      xi(1,3)   = adg(1,3) + 1._RKIND     ! xi_6
      xi(2,1)   = adg(2,1) + 1._RKIND     ! eta_4

!     Tmat and inverse (reference config)
      DO i=1, 3
         Tm0(i,1) = xi0(1,i)*(xi0(1,i) - 1._RKIND) ! xi**2 - xi
         Tm0(i,2) = xi0(2,i)*(xi0(2,i) - 1._RKIND) ! eta**2 - eta
         Tm0(i,3) = xi0(1,i)*xi0(2,i)              ! xi * eta
      END DO
      Tm0 = MAT_INV(Tm0, 3)

!     Tmat and inverse (current config)
      DO i=1, 3
         Tm(i,1) = xi(1,i)*(xi(1,i) - 1._RKIND)    ! xi**2 - xi
         Tm(i,2) = xi(2,i)*(xi(2,i) - 1._RKIND)    ! eta**2 - eta
         Tm(i,3) = xi(1,i)*xi(2,i)                 ! xi * eta
      END DO
      Tm = MAT_INV(Tm, 3)

!     v = Inv(T) * z (reference config)
      v0(1) = Tm0(1,1)*xi0(3,1) + Tm0(1,2)*xi0(3,2) + Tm0(1,3)*xi0(3,3)
      v0(2) = Tm0(2,1)*xi0(3,1) + Tm0(2,2)*xi0(3,2) + Tm0(2,3)*xi0(3,3)
      v0(3) = Tm0(3,1)*xi0(3,1) + Tm0(3,2)*xi0(3,2) + Tm0(3,3)*xi0(3,3)

!     Curvature coefficients (ref. config)
      bb_0(1,1) = 2._RKIND*v0(1)
      bb_0(2,2) = 2._RKIND*v0(2)
      bb_0(1,2) = v0(3)
      bb_0(2,1) = bb_0(1,2)

!     v = Inv(T) * z (current config)
      v(1) = Tm(1,1)*xi(3,1) + Tm(1,2)*xi(3,2) + Tm(1,3)*xi(3,3)
      v(2) = Tm(2,1)*xi(3,1) + Tm(2,2)*xi(3,2) + Tm(2,3)*xi(3,3)
      v(3) = Tm(3,1)*xi(3,1) + Tm(3,2)*xi(3,2) + Tm(3,3)*xi(3,3)

!     Curvature coefficients (current config)
      bb_x(1,1) = 2._RKIND*v(1)
      bb_x(2,2) = 2._RKIND*v(2)
      bb_x(1,2) = v(3)
      bb_x(2,1) = bb_x(1,2)

!     Now compute variation in bending strain
!     B1 bar
      B1b(:,1) = -Tm(:,1) * ((2._RKIND*xi(1,1)-1._RKIND)*v(1) +
     2   xi(2,1)*v(3))
      B1b(:,2) = -Tm(:,1) * ((2._RKIND*xi(2,1)-1._RKIND)*v(2) +
     2   xi(1,1)*v(3))

      B1b(:,3) = -Tm(:,2) * ((2._RKIND*xi(1,2)-1._RKIND)*v(1) +
     2   xi(2,2)*v(3))
      B1b(:,4) = -Tm(:,2) * ((2._RKIND*xi(2,2)-1._RKIND)*v(2) +
     2   xi(1,2)*v(3))

      B1b(:,5) = -Tm(:,3) * ((2._RKIND*xi(1,3)-1._RKIND)*v(1) +
     2   xi(2,3)*v(3))
      B1b(:,6) = -Tm(:,3) * ((2._RKIND*xi(2,3)-1._RKIND)*v(2) +
     2   xi(1,3)*v(3))

!     H1
      H1 = 0._RKIND
      H1(1, 1: 3) =  aCnv(:,1)*adg(2,1)
      H1(2, 1: 3) =  aCnv(:,2)*adg(2,1)
      H1(3, 1: 3) =  aCnv(:,1)*adg(2,2)
      H1(4, 1: 3) =  aCnv(:,2)*adg(2,2)
      H1(5, 1: 3) =  aCnv(:,1)*adg(2,3)
      H1(6, 1: 3) =  aCnv(:,2)*adg(2,3)

      H1(1, 4: 6) = -aCnv(:,1)*adg(1,1)
      H1(2, 4: 6) = -aCnv(:,2)*adg(1,1)
      H1(3, 4: 6) = -aCnv(:,1)*adg(1,2)
      H1(4, 4: 6) = -aCnv(:,2)*adg(1,2)
      H1(5, 4: 6) = -aCnv(:,1)*adg(1,3)
      H1(6, 4: 6) = -aCnv(:,2)*adg(1,3)

      H1(1,10:12) =  aCnv(:,1)
      H1(2,10:12) =  aCnv(:,2)
      H1(3,13:15) =  aCnv(:,1)
      H1(4,13:15) =  aCnv(:,2)
      H1(5,16:18) =  aCnv(:,1)
      H1(6,16:18) =  aCnv(:,2)

!     H2
      H2 = 0._RKIND
      H2( 1, 4) = -1._RKIND
      H2( 2, 5) = -1._RKIND
      H2( 3, 6) = -1._RKIND
      H2( 4, 7) = -1._RKIND
      H2( 5, 8) = -1._RKIND
      H2( 6, 9) = -1._RKIND
      H2( 7, 1) = -1._RKIND
      H2( 8, 2) = -1._RKIND
      H2( 9, 3) = -1._RKIND

      H2( 1, 7) =  1._RKIND
      H2( 2, 8) =  1._RKIND
      H2( 3, 9) =  1._RKIND
      H2( 4, 1) =  1._RKIND
      H2( 5, 2) =  1._RKIND
      H2( 6, 3) =  1._RKIND
      H2( 7, 4) =  1._RKIND
      H2( 8, 5) =  1._RKIND
      H2( 9, 6) =  1._RKIND

      H2(10, 4) = -1._RKIND
      H2(11, 5) = -1._RKIND
      H2(12, 6) = -1._RKIND
      H2(13, 7) = -1._RKIND
      H2(14, 8) = -1._RKIND
      H2(15, 9) = -1._RKIND
      H2(16, 1) = -1._RKIND
      H2(17, 2) = -1._RKIND
      H2(18, 3) = -1._RKIND
      do i=10, 18
         H2(i,i) = 1._RKIND
      end do

!     N matrix
      Nm = MAT_ID(3) - MAT_DYADPROD(nV, nV, 3)
!     M1, M2 matrices
      Mm(:,:,:) = 0._RKIND
      DO i=1, 2
         Mm(1,2,i) = -aCov(3,i)
         Mm(1,3,i) =  aCov(2,i)
         Mm(2,3,i) = -aCov(1,i)

!        Skew-symmetric
         Mm(2,1,i) = -Mm(1,2,i)
         Mm(3,1,i) = -Mm(1,3,i)
         Mm(3,2,i) = -Mm(2,3,i)
      END DO

!     H3 matrix
      H3 = 0._RKIND
      tmpA = MATMUL(Nm, Mm(:,:,1))
      tmpA = -tmpA / Jac
      H3(1,1) = a(1,4)*tmpA(1,1) + a(2,4)*tmpA(2,1) + a(3,4)*tmpA(3,1)
      H3(1,2) = a(1,4)*tmpA(1,2) + a(2,4)*tmpA(2,2) + a(3,4)*tmpA(3,2)
      H3(1,3) = a(1,4)*tmpA(1,3) + a(2,4)*tmpA(2,3) + a(3,4)*tmpA(3,3)

      H3(2,1) = a(1,5)*tmpA(1,1) + a(2,5)*tmpA(2,1) + a(3,5)*tmpA(3,1)
      H3(2,2) = a(1,5)*tmpA(1,2) + a(2,5)*tmpA(2,2) + a(3,5)*tmpA(3,2)
      H3(2,3) = a(1,5)*tmpA(1,3) + a(2,5)*tmpA(2,3) + a(3,5)*tmpA(3,3)

      H3(3,1) = a(1,6)*tmpA(1,1) + a(2,6)*tmpA(2,1) + a(3,6)*tmpA(3,1)
      H3(3,2) = a(1,6)*tmpA(1,2) + a(2,6)*tmpA(2,2) + a(3,6)*tmpA(3,2)
      H3(3,3) = a(1,6)*tmpA(1,3) + a(2,6)*tmpA(2,3) + a(3,6)*tmpA(3,3)

      tmpA = MATMUL(Nm, Mm(:,:,2))
      tmpA = -tmpA / Jac
      H3(1,4) = a(1,4)*tmpA(1,1) + a(2,4)*tmpA(2,1) + a(3,4)*tmpA(3,1)
      H3(1,5) = a(1,4)*tmpA(1,2) + a(2,4)*tmpA(2,2) + a(3,4)*tmpA(3,2)
      H3(1,6) = a(1,4)*tmpA(1,3) + a(2,4)*tmpA(2,3) + a(3,4)*tmpA(3,3)

      H3(2,4) = a(1,5)*tmpA(1,1) + a(2,5)*tmpA(2,1) + a(3,5)*tmpA(3,1)
      H3(2,5) = a(1,5)*tmpA(1,2) + a(2,5)*tmpA(2,2) + a(3,5)*tmpA(3,2)
      H3(2,6) = a(1,5)*tmpA(1,3) + a(2,5)*tmpA(2,3) + a(3,5)*tmpA(3,3)

      H3(3,4) = a(1,6)*tmpA(1,1) + a(2,6)*tmpA(2,1) + a(3,6)*tmpA(3,1)
      H3(3,5) = a(1,6)*tmpA(1,2) + a(2,6)*tmpA(2,2) + a(3,6)*tmpA(3,2)
      H3(3,6) = a(1,6)*tmpA(1,3) + a(2,6)*tmpA(2,3) + a(3,6)*tmpA(3,3)

      H3(1,10:12) = nV(:)
      H3(2,13:15) = nV(:)
      H3(3,16:18) = nV(:)

!     Variation in bending strain (Bb = -2*(B1b*H1*H2 + Tinv*H3*H2))
      H1H2 = MATMUL(H1, H2)
      Bb1  = MATMUL(B1b, H1H2)
      H3H2 = MATMUL(H3, H2)
      Bb1  = -2._RKIND*(Bb1 + MATMUL(Tm, H3H2))
      Bb   = RESHAPE(Bb1, SHAPE(Bb))

!     Update Bb for boundary elements
      IF (bFlag) THEN
         lFix = .FALSE.
         Im   = MAT_ID(3)
         DO j=lM%eNoN+1, eNoN
            IF (ptr(j) .NE. 0) CYCLE
            i = j - lM%eNoN
            p = i - 1
            f = i + 1
            IF (i .EQ. 1) p = 3
            IF (i .EQ. 3) f = 1
            IF (BTEST(lM%sbc(i,e),bType_fix)) THEN
!              eI = eI0 = aI0/|aI0| (reference config)
               aIi   = 1._RKIND/SQRT(NORM(a0(:,i)))
               eI(:) = a0(:,i) * aIi
!              nI = nI0 = eI0 x n0 (reference config)
               nI(1) = eI(2)*nV0(3) - eI(3)*nV0(2)
               nI(2) = eI(3)*nV0(1) - eI(1)*nV0(3)
               nI(3) = eI(1)*nV0(2) - eI(2)*nV0(1)
               nInI  = MAT_DYADPROD(nI, nI, 3)
            ELSE
!              eI = aI/|aI| (current config)
               aIi   = 1._RKIND/SQRT(NORM(a(:,i)))
               eI(:) = a(:,i)*aIi
!              nI = eI x n (currnt config)
               nI(1) = eI(2)*nV(3) - eI(3)*nV(2)
               nI(2) = eI(3)*nV(1) - eI(1)*nV(3)
               nI(3) = eI(1)*nV(2) - eI(2)*nV(1)
               cI    = NORM(a(:,i),a(:,p))*aIi*aIi
               nInI  = MAT_DYADPROD(nI, nI, 3)
               eIeI  = MAT_DYADPROD(eI, eI, 3)
               eIaP  = MAT_DYADPROD(eI, a(:,p), 3)
            END IF

!           Update Bb now
!           Free boundary conditions: assumed that the `artificial'
!           triangle is always located in the plane of the main element
!           of the patch
            IF (BTEST(lM%sbc(i,e),bType_free)) THEN
!              E_I
               IF (.NOT.lFix(i)) THEN
                  tmpA = -Im + 2._RKIND*eIeI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
!              E_P
               IF (.NOT.lFix(p)) THEN
                  tmpA = -2._RKIND*(cI*Im - 2._RKIND*cI*eIeI + aIi*eIaP)
                  Bb(:,:,p) = Bb(:,:,p) + MATMUL(Bb(:,:,j), tmpA)
               END IF
!              E_F
               IF (.NOT.lFix(f)) THEN
                  tmpA = 2._RKIND*((1._RKIND-cI)*Im -
     2               (1._RKIND-2._RKIND*cI)*eIeI -aIi*eIaP)
                  Bb(:,:,f) = Bb(:,:,f) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               Bb(:,:,j) = 0._RKIND

!           Hinged boundary conditions: a special case of simple support
!           in which no translation displacements are allowed.
            ELSE IF (BTEST(lM%sbc(i,e),bType_hing)) THEN
!              E_I
               IF (.NOT.lFix(i)) THEN
                  tmpA = -Im + 2._RKIND*eIeI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               lFix(p)   = .TRUE.
               lFix(f)   = .TRUE.
               Bb(:,:,p) = 0._RKIND
               Bb(:,:,f) = 0._RKIND
               Bb(:,:,j) = 0._RKIND

!           Fixed boundary condition: no displacements and no rotations
!           are allowed.
            ELSE IF (BTEST(lM%sbc(i,e),bType_fix)) THEN
               IF (.NOT.lFix(i)) THEN
                  tmpA = Im - 2._RKIND*nInI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               lFix(p)   = .TRUE.
               lFix(f)   = .TRUE.
               Bb(:,:,f) = 0._RKIND
               Bb(:,:,p) = 0._RKIND

!           Symmetric BCs (need to be verified)
            ELSE IF (BTEST(lM%sbc(i,e),bType_symm)) THEN
               IF (.NOT.lFix(i)) THEN
                  tmpA = Im - 2._RKIND*nInI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               Bb(:,:,j) = 0._RKIND
               tmpA(:,1) = eI(:)
               tmpA(:,2) = nV(:)
               tmpA(:,3) = nI(:)
               Bb(:,:,f) = MATMUL(Bb(:,:,f), tmpA)
               Bb(:,3,f) = 0._RKIND
               Bb(:,:,p) = MATMUL(Bb(:,:,p), tmpA)
               Bb(:,3,p) = 0._RKIND
            END IF
         END DO
      END IF

      RETURN
      END SUBROUTINE SHELLBENDCST
!####################################################################
!     Construct shell mechanics for higher order elements/NURBS
      SUBROUTINE SHELL3D (lM, g, eNoN, al, yl, dl, xl, bfl, lR, lK)
      USE COMMOD
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: g, eNoN
      REAL(KIND=RKIND), INTENT(IN) :: al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), xl(3,eNoN), bfl(3,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, k, l, a, b
      REAL(KIND=RKIND) :: amd, afl, rho, dmp, ht, w, wh, Jac0, Jac,
     2   fb(3), ud(3), nV0(3), nV(3), N(eNoN), x0(3,eNoN), xc(3,eNoN),
     3   Nx(2,eNoN), Nxx(3,eNoN), aCov0(3,2), aCov(3,2), aCnv0(3,2),
     4   aCnv(3,2), r0_xx(2,2,3), r_xx(2,2,3), aa_0(2,2), aa_x(2,2),
     5   bb_0(2,2), bb_x(2,2), Sm(3,2), Dm(3,3,3), Kc(3,3), Nm(3,3),
     6   Mm(3,3), KNmMm(3,3,2), Bm(3,3,eNoN), Bb(3,3,eNoN),
     7   D0Bm(3,3,eNoN), D1Bm(3,3,eNoN), D1Bb(3,3,eNoN), D2Bb(3,3,eNoN),
     8   T1, BmS, BbS, NxSNx, BmDBm, BmDBb, BbDBm, BbDBb

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp   = eq(cEq)%dmn(cDmn)%prop(damping)
      ht    = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      amd  = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl  = eq(cEq)%af*eq(cEq)%beta*dt*dt

      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

!     Get the reference configuration
      x0(:,:) = xl(:,:)

!     Get the current configuration
      DO a=1, eNoN
         xc(1,a) = x0(1,a) + dl(i,a)
         xc(2,a) = x0(2,a) + dl(j,a)
         xc(3,a) = x0(3,a) + dl(k,a)
      END DO

!     Define shape functions and their derivatives at Gauss point
      IF (lM%eType .EQ. eType_NRB) THEN
         N   = lM%N(:,g)
         Nx  = lM%Nx(:,:,g)
         Nxx = lM%Nxx(:,:,g)
      ELSE
         N   = lM%fs(1)%N(:,g)
         Nx  = lM%fs(1)%Nx(:,:,g)
         Nxx = lM%fs(1)%Nxx(:,:,g)
      END IF

c=====================================================================
c    TODO: Might have to call GNNxx for Jacobian transformation. Check
c    formulation again.
c
c     CALL GNNxx(3, eNoN, 2, lM%fs(1)%Nx(:,:,g), lM%fs(1)%Nxx(:,:,g),
c        xl, Nx, Nxx)
c
c=====================================================================

!---------------------------------------------------------------------
!     Compute preliminaries on the reference configuration
!     Covariant and contravariant bases (reference config)
      CALL GNNS(eNoN, Nx, x0, nV0, aCov0, aCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Second derivatives for computing curvature coeffs. (ref. config)
      r0_xx(:,:,:) = 0._RKIND
      DO a=1, eNoN
         r0_xx(1,1,:) = r0_xx(1,1,:) + Nxx(1,a)*x0(:,a)
         r0_xx(2,2,:) = r0_xx(2,2,:) + Nxx(2,a)*x0(:,a)
         r0_xx(1,2,:) = r0_xx(1,2,:) + Nxx(3,a)*x0(:,a)
      END DO
      r0_xx(2,1,:) = r0_xx(1,2,:)

!     Compute metric tensor and curvature coefficients (ref. config)
      aa_0 = 0._RKIND
      bb_0 = 0._RKIND
      DO l=1, nsd
         aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)
         aa_0(1,2) = aa_0(1,2) + aCov0(l,1)*aCov0(l,2)
         aa_0(2,1) = aa_0(2,1) + aCov0(l,2)*aCov0(l,1)
         aa_0(2,2) = aa_0(2,2) + aCov0(l,2)*aCov0(l,2)

         bb_0(1,1) = bb_0(1,1) + r0_xx(1,1,l)*nV0(l)
         bb_0(1,2) = bb_0(1,2) + r0_xx(1,2,l)*nV0(l)
         bb_0(2,1) = bb_0(2,1) + r0_xx(2,1,l)*nV0(l)
         bb_0(2,2) = bb_0(2,2) + r0_xx(2,2,l)*nV0(l)
      END DO

!---------------------------------------------------------------------
!     Now compute preliminaries on the current configuration
!     Covariant and contravariant bases (current/spatial config)
      CALL GNNS(eNoN, Nx, xc, nV, aCov, aCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac

!     Second derivatives for computing curvature coeffs. (cur. config)
      r_xx(:,:,:) = 0._RKIND
      DO a=1, eNoN
         r_xx(1,1,:) = r_xx(1,1,:) + Nxx(1,a)*xc(:,a)
         r_xx(2,2,:) = r_xx(2,2,:) + Nxx(2,a)*xc(:,a)
         r_xx(1,2,:) = r_xx(1,2,:) + Nxx(3,a)*xc(:,a)
      END DO
      r_xx(2,1,:) = r_xx(1,2,:)

!     Compute metric tensor and curvature coefficients (cur. config)
      aa_x = 0._RKIND
      bb_x = 0._RKIND
      DO l=1, nsd
         aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)
         aa_x(1,2) = aa_x(1,2) + aCov(l,1)*aCov(l,2)
         aa_x(2,1) = aa_x(2,1) + aCov(l,2)*aCov(l,1)
         aa_x(2,2) = aa_x(2,2) + aCov(l,2)*aCov(l,2)

         bb_x(1,1) = bb_x(1,1) + r_xx(1,1,l)*nV(l)
         bb_x(1,2) = bb_x(1,2) + r_xx(1,2,l)*nV(l)
         bb_x(2,1) = bb_x(2,1) + r_xx(2,1,l)*nV(l)
         bb_x(2,2) = bb_x(2,2) + r_xx(2,2,l)*nV(l)
      END DO

!---------------------------------------------------------------------
!     Compute stress resultants by integrating 2nd Piola Kirchhoff
!     stress and elasticity tensors through the shell thickness. These
!     resultants are computed in Voigt notation.
      CALL SHL_STRS_RES(eq(cEq)%dmn(cDmn), aa_0, aa_x, bb_0, bb_x, Sm,
     2   Dm)

!---------------------------------------------------------------------
!     Variation in the membrane strain
      Bm = 0._RKIND
      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*aCov(1,1)
         Bm(1,2,a) = Nx(1,a)*aCov(2,1)
         Bm(1,3,a) = Nx(1,a)*aCov(3,1)

         Bm(2,1,a) = Nx(2,a)*aCov(1,2)
         Bm(2,2,a) = Nx(2,a)*aCov(2,2)
         Bm(2,3,a) = Nx(2,a)*aCov(3,2)

         Bm(3,1,a) = Nx(2,a)*aCov(1,1) + Nx(1,a)*aCov(1,2)
         Bm(3,2,a) = Nx(2,a)*aCov(2,1) + Nx(1,a)*aCov(2,2)
         Bm(3,3,a) = Nx(2,a)*aCov(3,1) + Nx(1,a)*aCov(3,2)
      END DO

!---------------------------------------------------------------------
!     Variation in the bending strain
!     dB = -(B1 + B2) du; B1 = N_xx * n;
!     B2 = (r_xx Nm M1 Nx - r_xx N M2 Nx)

!     Second derivatives of the position vector (current)
      Kc(1,:) = r_xx(1,1,:)
      Kc(2,:) = r_xx(2,2,:)
      Kc(3,:) = r_xx(1,2,:)  + r_xx(2,1,:)

!     N matrix
      Nm = MAT_ID(3) - MAT_DYADPROD(nV, nV, 3)
      Nm = Nm / Jac

!     M1, M2 matrices
      DO l=1, 2
         Mm(:,:) = 0._RKIND
         Mm(1,2) = -aCov(3,l)
         Mm(1,3) =  aCov(2,l)
         Mm(2,3) = -aCov(1,l)

!        Skew-symmetric
         Mm(2,1) = -Mm(1,2)
         Mm(3,1) = -Mm(1,3)
         Mm(3,2) = -Mm(2,3)

         KNmMm(:,:,l) = MATMUL(Kc, MATMUL(Nm, Mm))
      END DO

!     Define variation in bending strain tensor (Bb), Voigt notation
      Bb = 0._RKIND
      DO a=1, eNoN
         Bb(1,:,a) = -Nxx(1,a)*nV(:)
         Bb(2,:,a) = -Nxx(2,a)*nV(:)
         Bb(3,:,a) = -Nxx(3,a)*nV(:)*2._RKIND
      END DO

      DO a=1, eNoN
         Bb(:,:,a) = Bb(:,:,a) + Nx(1,a)*KNmMm(:,:,2) -
     2      Nx(2,a)*KNmMm(:,:,1)
      END DO

!---------------------------------------------------------------------
!     Contribution to tangent matrices: Dm * Bm, Dm*Bb
      DO a=1, eNoN
         D0Bm(:,:,a) = MATMUL(Dm(:,:,1), Bm(:,:,a))
         D1Bm(:,:,a) = MATMUL(Dm(:,:,2), Bm(:,:,a))
         D1Bb(:,:,a) = MATMUL(Dm(:,:,2), Bb(:,:,a))
         D2Bb(:,:,a) = MATMUL(Dm(:,:,3), Bb(:,:,a))
      END DO

!---------------------------------------------------------------------
!     Acceleration and mass damping at the integration point
      ud = -fb
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))
         ud(3) = ud(3) + N(a)*(rho*(al(k,a)-bfl(3,a)) + dmp*yl(k,a))
      END DO

!     Local residue
      w  = lM%w(g)*Jac0
      wh = w*ht
      DO a=1, eNoN
         BmS = Bm(1,1,a)*Sm(1,1) + Bm(2,1,a)*Sm(2,1) + Bm(3,1,a)*Sm(3,1)
         BbS = Bb(1,1,a)*Sm(1,2) + Bb(2,1,a)*Sm(2,2) + Bb(3,1,a)*Sm(3,2)
         lR(1,a) = lR(1,a) + wh*N(a)*ud(1) + w*(BmS + BbS)

         BmS = Bm(1,2,a)*Sm(1,1) + Bm(2,2,a)*Sm(2,1) + Bm(3,2,a)*Sm(3,1)
         BbS = Bb(1,2,a)*Sm(1,2) + Bb(2,2,a)*Sm(2,2) + Bb(3,2,a)*Sm(3,2)
         lR(2,a) = lR(2,a) + wh*N(a)*ud(2) + w*(BmS + BbS)

         BmS = Bm(1,3,a)*Sm(1,1) + Bm(2,3,a)*Sm(2,1) + Bm(3,3,a)*Sm(3,1)
         BbS = Bb(1,3,a)*Sm(1,2) + Bb(2,3,a)*Sm(2,2) + Bb(3,3,a)*Sm(3,2)
         lR(3,a) = lR(3,a) + wh*N(a)*ud(3) + w*(BmS + BbS)
      END DO

!     Local stiffness
      amd = wh*amd
      afl = w*afl
      DO b=1, eNoN
         DO a=1, eNoN
!           Contribution from inertia and geometric stiffness
            NxSNx = Nx(1,a)*Nx(1,b)*Sm(1,1) + Nx(2,a)*Nx(2,b)*Sm(2,1)
     2            + Nx(1,a)*Nx(2,b)*Sm(3,1) + Nx(2,a)*Nx(1,b)*Sm(3,1)
            T1 = amd*N(a)*N(b) + afl*NxSNx

            lK(1,a,b) = lK(1,a,b) + T1
            lK(dof+2,a,b) = lK(dof+2,a,b) + T1
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + T1

!           Contribution from material stiffness
            BmDBm = Bm(1,1,a)*D0Bm(1,1,b) + Bm(2,1,a)*D0Bm(2,1,b)
     2            + Bm(3,1,a)*D0Bm(3,1,b)
            BmDBb = Bm(1,1,a)*D1Bb(1,1,b) + Bm(2,1,a)*D1Bb(2,1,b)
     2            + Bm(3,1,a)*D1Bb(3,1,b)
            BbDBm = Bb(1,1,a)*D1Bm(1,1,b) + Bb(2,1,a)*D1Bm(2,1,b)
     2            + Bb(3,1,a)*D1Bm(3,1,b)
            BbDBb = Bb(1,1,a)*D2Bb(1,1,b) + Bb(2,1,a)*D2Bb(2,1,b)
     2            + Bb(3,1,a)*D2Bb(3,1,b)
            lK(1,a,b) = lK(1,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb)

            BmDBm = Bm(1,1,a)*D0Bm(1,2,b) + Bm(2,1,a)*D0Bm(2,2,b)
     2            + Bm(3,1,a)*D0Bm(3,2,b)
            BmDBb = Bm(1,1,a)*D1Bb(1,2,b) + Bm(2,1,a)*D1Bb(2,2,b)
     2            + Bm(3,1,a)*D1Bb(3,2,b)
            BbDBm = Bb(1,1,a)*D1Bm(1,2,b) + Bb(2,1,a)*D1Bm(2,2,b)
     2            + Bb(3,1,a)*D1Bm(3,2,b)
            BbDBb = Bb(1,1,a)*D2Bb(1,2,b) + Bb(2,1,a)*D2Bb(2,2,b)
     2            + Bb(3,1,a)*D2Bb(3,2,b)
            lK(2,a,b) = lK(2,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb)

            BmDBm = Bm(1,1,a)*D0Bm(1,3,b) + Bm(2,1,a)*D0Bm(2,3,b)
     2            + Bm(3,1,a)*D0Bm(3,3,b)
            BmDBb = Bm(1,1,a)*D1Bb(1,3,b) + Bm(2,1,a)*D1Bb(2,3,b)
     2            + Bm(3,1,a)*D1Bb(3,3,b)
            BbDBm = Bb(1,1,a)*D1Bm(1,3,b) + Bb(2,1,a)*D1Bm(2,3,b)
     2            + Bb(3,1,a)*D1Bm(3,3,b)
            BbDBb = Bb(1,1,a)*D2Bb(1,3,b) + Bb(2,1,a)*D2Bb(2,3,b)
     2            + Bb(3,1,a)*D2Bb(3,3,b)
            lK(3,a,b) = lK(3,a,b) + afl*(BmDBm + BmDBb + BbDBm + BbDBb)

            BmDBm = Bm(1,2,a)*D0Bm(1,1,b) + Bm(2,2,a)*D0Bm(2,1,b)
     2            + Bm(3,2,a)*D0Bm(3,1,b)
            BmDBb = Bm(1,2,a)*D1Bb(1,1,b) + Bm(2,2,a)*D1Bb(2,1,b)
     2            + Bm(3,2,a)*D1Bb(3,1,b)
            BbDBm = Bb(1,2,a)*D1Bm(1,1,b) + Bb(2,2,a)*D1Bm(2,1,b)
     2            + Bb(3,2,a)*D1Bm(3,1,b)
            BbDBb = Bb(1,2,a)*D2Bb(1,1,b) + Bb(2,2,a)*D2Bb(2,1,b)
     2            + Bb(3,2,a)*D2Bb(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + afl*(BmDBm + BmDBb + BbDBm
     2         + BbDBb)

            BmDBm = Bm(1,2,a)*D0Bm(1,2,b) + Bm(2,2,a)*D0Bm(2,2,b)
     2            + Bm(3,2,a)*D0Bm(3,2,b)
            BmDBb = Bm(1,2,a)*D1Bb(1,2,b) + Bm(2,2,a)*D1Bb(2,2,b)
     2            + Bm(3,2,a)*D1Bb(3,2,b)
            BbDBm = Bb(1,2,a)*D1Bm(1,2,b) + Bb(2,2,a)*D1Bm(2,2,b)
     2            + Bb(3,2,a)*D1Bm(3,2,b)
            BbDBb = Bb(1,2,a)*D2Bb(1,2,b) + Bb(2,2,a)*D2Bb(2,2,b)
     2            + Bb(3,2,a)*D2Bb(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + afl*(BmDBm + BmDBb + BbDBm
     2         + BbDBb)

            BmDBm = Bm(1,2,a)*D0Bm(1,3,b) + Bm(2,2,a)*D0Bm(2,3,b)
     2            + Bm(3,2,a)*D0Bm(3,3,b)
            BmDBb = Bm(1,2,a)*D1Bb(1,3,b) + Bm(2,2,a)*D1Bb(2,3,b)
     2            + Bm(3,2,a)*D1Bb(3,3,b)
            BbDBm = Bb(1,2,a)*D1Bm(1,3,b) + Bb(2,2,a)*D1Bm(2,3,b)
     2            + Bb(3,2,a)*D1Bm(3,3,b)
            BbDBb = Bb(1,2,a)*D2Bb(1,3,b) + Bb(2,2,a)*D2Bb(2,3,b)
     2            + Bb(3,2,a)*D2Bb(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + afl*(BmDBm + BmDBb + BbDBm
     2         + BbDBb)

            BmDBm = Bm(1,3,a)*D0Bm(1,1,b) + Bm(2,3,a)*D0Bm(2,1,b)
     2            + Bm(3,3,a)*D0Bm(3,1,b)
            BmDBb = Bm(1,3,a)*D1Bb(1,1,b) + Bm(2,3,a)*D1Bb(2,1,b)
     2            + Bm(3,3,a)*D1Bb(3,1,b)
            BbDBm = Bb(1,3,a)*D1Bm(1,1,b) + Bb(2,3,a)*D1Bm(2,1,b)
     2            + Bb(3,3,a)*D1Bm(3,1,b)
            BbDBb = Bb(1,3,a)*D2Bb(1,1,b) + Bb(2,3,a)*D2Bb(2,1,b)
     2            + Bb(3,3,a)*D2Bb(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + afl*(BmDBm + BmDBb
     2         + BbDBm + BbDBb)

            BmDBm = Bm(1,3,a)*D0Bm(1,2,b) + Bm(2,3,a)*D0Bm(2,2,b)
     2            + Bm(3,3,a)*D0Bm(3,2,b)
            BmDBb = Bm(1,3,a)*D1Bb(1,2,b) + Bm(2,3,a)*D1Bb(2,2,b)
     2            + Bm(3,3,a)*D1Bb(3,2,b)
            BbDBm = Bb(1,3,a)*D1Bm(1,2,b) + Bb(2,3,a)*D1Bm(2,2,b)
     2            + Bb(3,3,a)*D1Bm(3,2,b)
            BbDBb = Bb(1,3,a)*D2Bb(1,2,b) + Bb(2,3,a)*D2Bb(2,2,b)
     2            + Bb(3,3,a)*D2Bb(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + afl*(BmDBm + BmDBb
     2         + BbDBm + BbDBb)

            BmDBm = Bm(1,3,a)*D0Bm(1,3,b) + Bm(2,3,a)*D0Bm(2,3,b)
     2            + Bm(3,3,a)*D0Bm(3,3,b)
            BmDBb = Bm(1,3,a)*D1Bb(1,3,b) + Bm(2,3,a)*D1Bb(2,3,b)
     2            + Bm(3,3,a)*D1Bb(3,3,b)
            BbDBm = Bb(1,3,a)*D1Bm(1,3,b) + Bb(2,3,a)*D1Bm(2,3,b)
     2            + Bb(3,3,a)*D1Bm(3,3,b)
            BbDBb = Bb(1,3,a)*D2Bb(1,3,b) + Bb(2,3,a)*D2Bb(2,3,b)
     2            + Bb(3,3,a)*D2Bb(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + afl*(BmDBm + BmDBb
     2         + BbDBm + BbDBb)
         END DO
      END DO

      RETURN
      END SUBROUTINE SHELL3D
!####################################################################
!     Compute stress resultants for shell elements
      SUBROUTINE SHL_STRS_RES(lDmn, aa_0, aa_x, bb_0, bb_x, Sm, Dm)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(IN) :: aa_0(2,2), aa_x(2,2), bb_0(2,2),
     2   bb_x(2,2)
      REAL(KIND=RKIND), INTENT(OUT) :: Sm(3,2), Dm(3,3,3)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: g
      REAL(KIND=RKIND) :: ht, nu, xis, xi(3), wh(3), wl(3), gg_0(2,2),
     2   gg_x(2,2), Sml(3), Dml(3,3)

!     Set shell thickness
      ht = lDmn%prop(shell_thickness)
      nu = lDmn%prop(poisson_ratio)

!     Check for incompressibility
      flag = .FALSE.
      IF (ISZERO(nu-0.5_RKIND)) flag = .TRUE.

!     Set integration parameters (Gauss coordinates and weights)
      wh(1) = 5._RKIND/9._RKIND
      wh(2) = 5._RKIND/9._RKIND
      wh(3) = 8._RKIND/9._RKIND

      xi(1) = -SQRT(0.6_RKIND)
      xi(2) =  SQRT(0.6_RKIND)
      xi(3) =  0._RKIND

!     Initialize stress-resultants
      Sm = 0._RKIND
      Dm = 0._RKIND

      Sml = 0._RKIND
      Dml = 0._RKIND

!     Gauss integration through shell thickness
      DO g=1, 3
!        Local shell thickness coordinate
         xis = .5_RKIND*ht*xi(g)

!        Metric coefficients in shell continuum (ref/cur)
         gg_0(:,:) = aa_0(:,:) - 2._RKIND*xis*bb_0(:,:)
         gg_x(:,:) = aa_x(:,:) - 2._RKIND*xis*bb_x(:,:)

!        Get 2nd Piola-Kirchhoff and elasticity tensors
         IF (flag) THEN
!           For incompressible materials
            CALL GETPK2CC_SHLi(lDmn, gg_0, gg_x, Sml, Dml)
         ELSE
!           For compressible materials
            CALL GETPK2CC_SHLc(lDmn, gg_0, gg_x, Sml, Dml)
         END IF

         wl(1) = .5_RKIND*wh(g)*ht
         wl(2) = wl(1) * xis
         wl(3) = wl(2) * xis

         Sm(:,1) = Sm(:,1) + wl(1)*Sml(:)
         Sm(:,2) = Sm(:,2) + wl(2)*Sml(:)

         Dm(:,:,1) = Dm(:,:,1) + wl(1)*Dml(:,:)
         Dm(:,:,2) = Dm(:,:,2) + wl(2)*Dml(:,:)
         Dm(:,:,3) = Dm(:,:,3) + wl(3)*Dml(:,:)
      END DO

      RETURN
      END SUBROUTINE SHL_STRS_RES
!####################################################################
!     Set follower pressure load/net traction on shells. The traction
!     on shells is treated as body force and the subroutine is called
!     from BF.f
      SUBROUTINE SHELLBF (eNoN, w, N, Nx, dl, xl, tfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   dl(tDof,eNoN), xl(3,eNoN), tfl(eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: T1, afl, wl, tfn, nV(3), gCov(3,2), gCnv(3,2),
     2   xc(3,eNoN), lKP(3)

      afl  = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

!     Get the current configuration and traction vector
      tfn = 0._RKIND
      DO a=1, eNoN
         xc(1,a) = xl(1,a) + dl(i,a)
         xc(2,a) = xl(2,a) + dl(j,a)
         xc(3,a) = xl(3,a) + dl(k,a)

         tfn = tfn + N(a)*tfl(a)
      END DO
      wl = w * tfn

!     Covariant and contravariant bases in current config
      CALL GNNS(eNoN, Nx, xc, nV, gCov, gCnv)

!     Local residue
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nV(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nV(2)
         lR(3,a) = lR(3,a) - wl*N(a)*nV(3)
      END DO

!     Local stiffness: mass matrix and stiffness contribution due to
!     follower traction load
      T1 = afl*wl*0.5_RKIND
      DO b=1, eNoN
         DO a=1, eNoN
            lKp(:) = gCov(:,1)*(N(b)*Nx(2,a) - N(a)*Nx(2,b))
     2             - gCov(:,2)*(N(b)*Nx(1,a) - N(a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) - T1*lKp(3)
            lK(3,a,b) = lK(3,a,b) + T1*lKp(2)

            lK(dof+1,a,b) = lK(dof+1,a,b) + T1*lKp(3)
            lK(dof+3,a,b) = lK(dof+3,a,b) - T1*lKp(1)

            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) - T1*lKp(2)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1*lKp(1)
         END DO
      END DO

      RETURN
      END SUBROUTINE SHELLBF
!####################################################################


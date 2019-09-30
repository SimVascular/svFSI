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
!     using linear triangle finite elements and IGA.
!
!--------------------------------------------------------------------

      SUBROUTINE SHELLTRI (lM, e, eNoN, al, yl, dl, xl, bfl, ptr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: e, eNoN, ptr(eNoN)
      REAL(KIND=8), INTENT(IN) :: al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), xl(3,eNoN), bfl(3,eNoN)

      LOGICAL :: bFlag, setIt(3)
      INTEGER :: i, j, k, a, b, g
      REAL(KIND=8) :: rho, dmp, elM, nu, ht, T1, amd, afl, w, Jac0, Jac,
     2   ud(3), fb(3), nV0(3), nV(3), gCov0(3,2), gCnv0(3,2), gCov(3,2),
     3   gCnv(3,2), x0(3,eNoN), xc(3,eNoN), eLoc(3,3), Jm(2,2), Qm(3,3),
     4   Dm(3,3), Em(3), Eb(3), DEm(3), DEb(3), Bm(3,3,eNoN),
     5   Bb(3,3,eNoN), DBm(3,3,eNoN), DBb(3,3,eNoN), BtDE, BtDB, NxSNx

      REAL(KIND=8), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:),
     2   tmpX(:,:)

!     Note that for triangular elements, eNoN=6 and lM%eNoN=3
      ALLOCATE(N(lM%eNoN), Nx(2,lM%eNoN), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN))

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp   = eq(cEq)%dmn(cDmn)%prop(damping)
      elM   = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu    = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      ht    = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      amd  = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl  = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

!     Determine if current element is a boundary element
      bFlag = .FALSE.
      DO a=lM%eNoN+1, eNoN
         IF (ptr(a) .EQ. 0) THEN
            bFlag = .TRUE.
            EXIT
         END IF
      END DO

!     Get the reference configuration
      x0(:,:) = xl(:,:)

!     Get the current configuration
      DO a=1, eNoN
         xc(1,a) = x0(1,a) + dl(i,a)
         xc(2,a) = x0(2,a) + dl(j,a)
         xc(3,a) = x0(3,a) + dl(k,a)
      END DO
      Nx(:,:) = lM%Nx(:,:,1)

!     Covariant and contravariant bases in reference config
      ALLOCATE(tmpX(nsd,lM%eNoN))
      tmpX = x0(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, Nx, tmpX, nV0, gCov0, gCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      tmpX = xc(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, Nx, tmpX, nV, gCov, gCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac
      DEALLOCATE(tmpX)

!     Define local coordinates, Jacobian (J) and its inverse tensor
!     defined in Voigt notation such that Q*E = J^{-T} * E * J^{-1})
      eLoc(:,1) = gCov0(:,1)/SQRT(NORM(gCov0(:,1)))
      eLoc(:,3) = nV0(:)

      eLoc(1,2) = eLoc(2,3)*eLoc(3,1) - eLoc(3,3)*eLoc(2,1)
      eLoc(2,2) = eLoc(3,3)*eLoc(1,1) - eLoc(1,3)*eLoc(3,1)
      eLoc(3,2) = eLoc(1,3)*eLoc(2,1) - eLoc(2,3)*eLoc(1,1)

      Jm(1,1) = NORM(gCov0(:,1), eLoc(:,1))
      Jm(1,2) = NORM(gCov0(:,2), eLoc(:,1))
      Jm(2,1) = NORM(gCov0(:,1), eLoc(:,2))
      Jm(2,2) = NORM(gCov0(:,2), eLoc(:,2))

      Qm = 0D0
      Qm(1,1) = 1D0/Jm(1,1)**2
      Qm(2,1) = (-Jm(1,2)/(Jm(1,1)*Jm(2,2)))**2
      Qm(2,2) = 1D0/Jm(2,2)**2
      Qm(2,3) = -Jm(1,2)/(Jm(1,1)*Jm(2,2)**2)
      Qm(3,1) = -2D0*Jm(1,2)/(Jm(2,2)*Jm(1,1)**2)
      Qm(3,3) = 1D0/(Jm(1,1)*Jm(2,2))

!     Material tensor D, isotropic St Venant Kirchhoff model
      Dm = 0D0
      Dm(1,1) = 1D0
      Dm(1,2) = nu
      Dm(2,1) = nu
      Dm(2,2) = 1D0
      Dm(3,3) = 5D-1*(1D0-nu)
      Dm = elM/(1D0-nu**2) * Dm

!     Compute Qt * D * Q
      Dm = MATMUL(Dm, Qm)
      Dm = MATMUL(TRANSPOSE(Qm), Dm)

!     Membrane strain and its variation
!     Define membrane strain tensor (Em), Voigt notation
      Em = 0D0
      DO g=1, nsd
         Em(1) = Em(1) + gCov(g,1)*gCov(g,1) - gCov0(g,1)*gCov0(g,1)
         Em(2) = Em(2) + gCov(g,2)*gCov(g,2) - gCov0(g,2)*gCov0(g,2)
         Em(3) = Em(3) + gCov(g,1)*gCov(g,2) - gCov0(g,1)*gCov0(g,2)
      END DO
      Em(1) = 5D-1 * Em(1)
      Em(2) = 5D-1 * Em(2)

!     Define variation in membrane strain only for the main element
      Bm = 0D0
      DO a=1, lM%eNoN
         Bm(1,1,a) = Nx(1,a)*gCov(1,1)
         Bm(1,2,a) = Nx(1,a)*gCov(2,1)
         Bm(1,3,a) = Nx(1,a)*gCov(3,1)

         Bm(2,1,a) = Nx(2,a)*gCov(1,2)
         Bm(2,2,a) = Nx(2,a)*gCov(2,2)
         Bm(2,3,a) = Nx(2,a)*gCov(3,2)

         Bm(3,1,a) = Nx(2,a)*gCov(1,1) + Nx(1,a)*gCov(1,2)
         Bm(3,2,a) = Nx(2,a)*gCov(2,1) + Nx(1,a)*gCov(2,2)
         Bm(3,3,a) = Nx(2,a)*gCov(3,1) + Nx(1,a)*gCov(3,2)
      END DO

!     Zero-out Bm for fixed BC on boundary elements
      IF (bFlag) THEN
         setIt = .FALSE.
         DO a=lM%eNoN+1, eNoN
            IF (ptr(a) .EQ. 0) THEN
               b = a - lM%eNoN
               IF (BTEST(lM%sbc(b,e),bType_fix)) setIt(b) = .TRUE.
            END IF
         END DO

         DO a=1, lM%eNoN
            IF (setIt(a)) THEN
               DO b=1, lM%eNoN
                  IF (a .EQ. b) CYCLE
                  Bm(:,:,b) = 0D0
               END DO
            END IF
         END DO
      END IF

!     Bending strain and its variation for triangular elements
      CALL SHELLBENDTRI(lM, e, ptr, x0, xc, Eb, Bb)

!     Contribution to residue and tangent matrices from membrane strain
!     D * Em
      DEm(1) = Dm(1,1)*Em(1) + Dm(1,2)*Em(2) + Dm(1,3)*Em(3)
      DEm(2) = Dm(2,1)*Em(1) + Dm(2,2)*Em(2) + Dm(2,3)*Em(3)
      DEm(3) = Dm(3,1)*Em(1) + Dm(3,2)*Em(2) + Dm(3,3)*Em(3)

!     D * Bm
      DO a=1, lM%eNoN
         DBm(1,1,a) = Dm(1,1)*Bm(1,1,a) + Dm(1,2)*Bm(2,1,a) +
     2      Dm(1,3)*Bm(3,1,a)
         DBm(1,2,a) = Dm(1,1)*Bm(1,2,a) + Dm(1,2)*Bm(2,2,a) +
     2      Dm(1,3)*Bm(3,2,a)
         DBm(1,3,a) = Dm(1,1)*Bm(1,3,a) + Dm(1,2)*Bm(2,3,a) +
     2      Dm(1,3)*Bm(3,3,a)

         DBm(2,1,a) = Dm(2,1)*Bm(1,1,a) + Dm(2,2)*Bm(2,1,a) +
     2      Dm(2,3)*Bm(3,1,a)
         DBm(2,2,a) = Dm(2,1)*Bm(1,2,a) + Dm(2,2)*Bm(2,2,a) +
     2      Dm(2,3)*Bm(3,2,a)
         DBm(2,3,a) = Dm(2,1)*Bm(1,3,a) + Dm(2,2)*Bm(2,3,a) +
     2      Dm(2,3)*Bm(3,3,a)

         DBm(3,1,a) = Dm(3,1)*Bm(1,1,a) + Dm(3,2)*Bm(2,1,a) +
     2      Dm(3,3)*Bm(3,1,a)
         DBm(3,2,a) = Dm(3,1)*Bm(1,2,a) + Dm(3,2)*Bm(2,2,a) +
     2      Dm(3,3)*Bm(3,2,a)
         DBm(3,3,a) = Dm(3,1)*Bm(1,3,a) + Dm(3,2)*Bm(2,3,a) +
     2      Dm(3,3)*Bm(3,3,a)
      END DO

!     Contribution to residue and tangent matrices from bending strain
!     D * Eb
      DEb(1) = Dm(1,1)*Eb(1) + Dm(1,2)*Eb(2) + Dm(1,3)*Eb(3)
      DEb(2) = Dm(2,1)*Eb(1) + Dm(2,2)*Eb(2) + Dm(2,3)*Eb(3)
      DEb(3) = Dm(3,1)*Eb(1) + Dm(3,2)*Eb(2) + Dm(3,3)*Eb(3)

!     D * Bb
      DO a=1, eNoN
         DBb(1,1,a) = Dm(1,1)*Bb(1,1,a) + Dm(1,2)*Bb(2,1,a) +
     2      Dm(1,3)*Bb(3,1,a)
         DBb(1,2,a) = Dm(1,1)*Bb(1,2,a) + Dm(1,2)*Bb(2,2,a) +
     2      Dm(1,3)*Bb(3,2,a)
         DBb(1,3,a) = Dm(1,1)*Bb(1,3,a) + Dm(1,2)*Bb(2,3,a) +
     2      Dm(1,3)*Bb(3,3,a)

         DBb(2,1,a) = Dm(2,1)*Bb(1,1,a) + Dm(2,2)*Bb(2,1,a) +
     2      Dm(2,3)*Bb(3,1,a)
         DBb(2,2,a) = Dm(2,1)*Bb(1,2,a) + Dm(2,2)*Bb(2,2,a) +
     2      Dm(2,3)*Bb(3,2,a)
         DBb(2,3,a) = Dm(2,1)*Bb(1,3,a) + Dm(2,2)*Bb(2,3,a) +
     2      Dm(2,3)*Bb(3,3,a)

         DBb(3,1,a) = Dm(3,1)*Bb(1,1,a) + Dm(3,2)*Bb(2,1,a) +
     2      Dm(3,3)*Bb(3,1,a)
         DBb(3,2,a) = Dm(3,1)*Bb(1,2,a) + Dm(3,2)*Bb(2,2,a) +
     2      Dm(3,3)*Bb(3,2,a)
         DBb(3,3,a) = Dm(3,1)*Bb(1,3,a) + Dm(3,2)*Bb(2,3,a) +
     2      Dm(3,3)*Bb(3,3,a)
      END DO

!     Contribution to residue and stiffness matrices due to inertia and
!     body forces
      lR = 0D0
      lK = 0D0
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
               T1 = w*amd*N(a)*N(b)
               lK(1,a,b) = lK(1,a,b) + T1
               lK(dof+2,a,b) = lK(dof+2,a,b) + T1
               lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + T1
            END DO
         END DO
      END DO

!     Contribution from membrane strain related terms only on the main
!     triangular element
      T1 = Jac0 * ht/2D0
      DO a=1, lM%eNoN
         BtDE = Bm(1,1,a)*DEm(1) + Bm(2,1,a)*DEm(2) + Bm(3,1,a)*DEm(3)
         lR(1,a) = lR(1,a) + T1*BtDE

         BtDE = Bm(1,2,a)*DEm(1) + Bm(2,2,a)*DEm(2) + Bm(3,2,a)*DEm(3)
         lR(2,a) = lR(2,a) + T1*BtDE

         BtDE = Bm(1,3,a)*DEm(1) + Bm(2,3,a)*DEm(2) + Bm(3,3,a)*DEm(3)
         lR(3,a) = lR(3,a) + T1*BtDE

         DO b=1, lM%eNoN
!           Geometric stiffness
            NxSNx = ( Nx(1,a)*Nx(1,b)*DEm(1) + Nx(2,a)*Nx(2,b)*DEm(2) +
     2                Nx(1,a)*Nx(2,b)*DEm(3) + Nx(2,a)*Nx(1,b)*DEm(3) )

!           Material stiffness
            BtDB = Bm(1,1,a)*DBm(1,1,b) + Bm(2,1,a)*DBm(2,1,b) +
     2             Bm(3,1,a)*DBm(3,1,b)
            lK(1,a,b) = lK(1,a,b) + afl*T1*(BtDB + NxSNx)

            BtDB = Bm(1,1,a)*DBm(1,2,b) + Bm(2,1,a)*DBm(2,2,b) +
     2             Bm(3,1,a)*DBm(3,2,b)
            lK(2,a,b) = lK(2,a,b) + afl*T1*BtDB

            BtDB = Bm(1,1,a)*DBm(1,3,b) + Bm(2,1,a)*DBm(2,3,b) +
     2             Bm(3,1,a)*DBm(3,3,b)
            lK(3,a,b) = lK(3,a,b) + afl*T1*BtDB

            BtDB = Bm(1,2,a)*DBm(1,1,b) + Bm(2,2,a)*DBm(2,1,b) +
     2             Bm(3,2,a)*DBm(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + afl*T1*BtDB

            BtDB = Bm(1,2,a)*DBm(1,2,b) + Bm(2,2,a)*DBm(2,2,b) +
     2             Bm(3,2,a)*DBm(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + afl*T1*(BtDB + NxSNx)

            BtDB = Bm(1,2,a)*DBm(1,3,b) + Bm(2,2,a)*DBm(2,3,b) +
     2             Bm(3,2,a)*DBm(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + afl*T1*BtDB

            BtDB = Bm(1,3,a)*DBm(1,1,b) + Bm(2,3,a)*DBm(2,1,b) +
     2             Bm(3,3,a)*DBm(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + afl*T1*BtDB

            BtDB = Bm(1,3,a)*DBm(1,2,b) + Bm(2,3,a)*DBm(2,2,b) +
     2             Bm(3,3,a)*DBm(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + afl*T1*BtDB

            BtDB = Bm(1,3,a)*DBm(1,3,b) + Bm(2,3,a)*DBm(2,3,b) +
     2             Bm(3,3,a)*DBm(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + afl*T1*(BtDB + NxSNx)
         END DO
      END DO

!     Contribution from bending strain related terms for the main
!     triangle and its neighbors. Geometric stiffness is ignored
      T1 = Jac0 * ht**3 /24D0
      DO a=1, eNoN
         BtDE = Bb(1,1,a)*DEb(1) + Bb(2,1,a)*DEb(2) + Bb(3,1,a)*DEb(3)
         lR(1,a) = lR(1,a) + T1*BtDE

         BtDE = Bb(1,2,a)*DEb(1) + Bb(2,2,a)*DEb(2) + Bb(3,2,a)*DEb(3)
         lR(2,a) = lR(2,a) + T1*BtDE

         BtDE = Bb(1,3,a)*DEb(1) + Bb(2,3,a)*DEb(2) + Bb(3,3,a)*DEb(3)
         lR(3,a) = lR(3,a) + T1*BtDE

         DO b=1, eNoN
            BtDB = Bb(1,1,a)*DBb(1,1,b) + Bb(2,1,a)*DBb(2,1,b) +
     2             Bb(3,1,a)*DBb(3,1,b)
            lK(1,a,b) = lK(1,a,b) + afl*T1*BtDB

            BtDB = Bb(1,1,a)*DBb(1,2,b) + Bb(2,1,a)*DBb(2,2,b) +
     2             Bb(3,1,a)*DBb(3,2,b)
            lK(2,a,b) = lK(2,a,b) + afl*T1*BtDB

            BtDB = Bb(1,1,a)*DBb(1,3,b) + Bb(2,1,a)*DBb(2,3,b) +
     2             Bb(3,1,a)*DBb(3,3,b)
            lK(3,a,b) = lK(3,a,b) + afl*T1*BtDB

            BtDB = Bb(1,2,a)*DBb(1,1,b) + Bb(2,2,a)*DBb(2,1,b) +
     2             Bb(3,2,a)*DBb(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + afl*T1*BtDB

            BtDB = Bb(1,2,a)*DBb(1,2,b) + Bb(2,2,a)*DBb(2,2,b) +
     2             Bb(3,2,a)*DBb(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + afl*T1*BtDB

            BtDB = Bb(1,2,a)*DBb(1,3,b) + Bb(2,2,a)*DBb(2,3,b) +
     2             Bb(3,2,a)*DBb(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + afl*T1*BtDB

            BtDB = Bb(1,3,a)*DBb(1,1,b) + Bb(2,3,a)*DBb(2,1,b) +
     2             Bb(3,3,a)*DBb(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + afl*T1*BtDB

            BtDB = Bb(1,3,a)*DBb(1,2,b) + Bb(2,3,a)*DBb(2,2,b) +
     2             Bb(3,3,a)*DBb(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + afl*T1*BtDB

            BtDB = Bb(1,3,a)*DBb(1,3,b) + Bb(2,3,a)*DBb(2,3,b) +
     2             Bb(3,3,a)*DBb(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + afl*T1*BtDB
         END DO
      END DO

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
      END SUBROUTINE SHELLTRI
!--------------------------------------------------------------------
!     This subroutine computes bending strain, Eb, and its variational
!     derivative, Bb, for triangular shell elements
      SUBROUTINE SHELLBENDTRI(lM, e, ptr, x0, xc, Eb, Bb)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: e, ptr(6)
      REAL(KIND=8), INTENT(INOUT) ::  x0(3,6), xc(3,6), Eb(3), Bb(3,3,6)

      LOGICAL :: bFlag, lFix(3)
      INTEGER :: i, j, p, f, eNoN

      REAL(KIND=8) :: Jac0, Jac, cI, aIi, nV0(3), nV(3), eI(3), nI(3),
     2   eIeI(3,3), nInI(3,3), eIaP(3,3), Im(3,3), gCov(3,2), gCnv(3,2),
     3   gCov0(3,2), gCnv0(3,2)

      REAL(KIND=8) :: a0(3,6), a(3,6), adg0(3,3), adg(3,3), xi0(3,3),
     2   xi(3,3), Tm0(3,3), Tm(3,3), v0(3), v(3), B1b(3,6), Nm(3,3),
     2   Mm(3,3,2), H1(6,18), H2(18,18), H3(3,18), H1H2(6,18),
     3   H3H2(3,18), Bb1(3,18), tmpA(3,3)

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
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV0, gCov0, gCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      tmpA = xc(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV, gCov, gCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac

!     Update position vector of surrounding nodes, if boundary elements
      IF (bFlag) THEN
         DO j=lM%eNoN+1, eNoN
            IF (ptr(j) .NE. 0) CYCLE
            i = j - lM%eNoN
            p = i - 1
            IF (i .EQ. 1) p = 3
!           Reference config
!           eI = eI0 = aI0/|aI0| (reference config)
            aIi   = 1D0/SQRT(NORM(a0(:,i)))
            eI(:) = a0(:,i) * aIi
!           nI = nI0 = eI0 x n0 (reference config)
            nI(1) = eI(2)*nV0(3) - eI(3)*nV0(2)
            nI(2) = eI(3)*nV0(1) - eI(1)*nV0(3)
            nI(3) = eI(1)*nV0(2) - eI(2)*nV0(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            x0(1,j) = 2D0*(nInI(1,1)*a0(1,p) + nInI(1,2)*a0(2,p) +
     2         nInI(1,3)*a0(3,p)) + x0(1,i)
            x0(2,j) = 2D0*(nInI(2,1)*a0(1,p) + nInI(2,2)*a0(2,p) +
     2         nInI(2,3)*a0(3,p)) + x0(2,i)
            x0(3,j) = 2D0*(nInI(3,1)*a0(1,p) + nInI(3,2)*a0(2,p) +
     2         nInI(3,3)*a0(3,p)) + x0(3,i)

!           Current config
!           eI = aI/|aI| (current config)
            aIi   = 1D0/SQRT(NORM(a(:,i)))
            eI(:) = a(:,i)*aIi
!           nI = eI x n (currnt config)
            nI(1) = eI(2)*nV(3) - eI(3)*nV(2)
            nI(2) = eI(3)*nV(1) - eI(1)*nV(3)
            nI(3) = eI(1)*nV(2) - eI(2)*nV(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            xc(1,j) = 2D0*(nInI(1,1)*a(1,p) + nInI(1,2)*a(2,p) +
     2         nInI(1,3)*a(3,p)) + xc(1,i)
            xc(2,j) = 2D0*(nInI(2,1)*a(1,p) + nInI(2,2)*a(2,p) +
     2         nInI(2,3)*a(3,p)) + xc(2,i)
            xc(3,j) = 2D0*(nInI(3,1)*a(1,p) + nInI(3,2)*a(2,p) +
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
      adg0(1,1) = NORM(a0(:,4), gCnv0(:,1)) ! xi_4
      adg0(1,2) = NORM(a0(:,5), gCnv0(:,1)) ! xi_5
      adg0(1,3) = NORM(a0(:,6), gCnv0(:,1)) ! xi_6

      adg0(2,1) = NORM(a0(:,4), gCnv0(:,2)) ! eta_4
      adg0(2,2) = NORM(a0(:,5), gCnv0(:,2)) ! eta_5
      adg0(2,3) = NORM(a0(:,6), gCnv0(:,2)) ! eta_6

      adg0(3,1) = NORM(a0(:,4), nV0)        ! z_4
      adg0(3,2) = NORM(a0(:,5), nV0)        ! z_5
      adg0(3,3) = NORM(a0(:,6), nV0)        ! z_6

!     a.gCnv (current config)
      adg(1,1)  = NORM(a(:,4), gCnv(:,1)) ! xi_4
      adg(1,2)  = NORM(a(:,5), gCnv(:,1)) ! xi_5
      adg(1,3)  = NORM(a(:,6), gCnv(:,1)) ! xi_6

      adg(2,1)  = NORM(a(:,4), gCnv(:,2)) ! eta_4
      adg(2,2)  = NORM(a(:,5), gCnv(:,2)) ! eta_5
      adg(2,3)  = NORM(a(:,6), gCnv(:,2)) ! eta_6

      adg(3,1)  = NORM(a(:,4), nV)        ! z_4
      adg(3,2)  = NORM(a(:,5), nV)        ! z_5
      adg(3,3)  = NORM(a(:,6), nV)        ! z_6

!     Xi matrix (reference config)
      xi0(:,:)  = adg0(:,:)
      xi0(1,1)  = adg0(1,1) + 1D0
      xi0(2,2)  = adg0(2,2) + 1D0

!     Xi matrix (current config)
      xi(:,:)   = adg(:,:)
      xi(1,1)   = adg(1,1) + 1D0
      xi(2,2)   = adg(2,2) + 1D0

!     Tmat and inverse (reference config)
      DO i=1, 3
         Tm0(i,1) = xi0(1,i)**2 - xi0(1,i) ! xi**2 - xi
         Tm0(i,2) = xi0(2,i)**2 - xi0(2,i) ! eta**2 - eta
         Tm0(i,3) = xi0(1,i)*xi0(2,i)      ! xi * eta
      END DO
      Tm0 = MAT_INV(Tm0, 3)

!     Tmat and inverse (current config)
      DO i=1, 3
         Tm(i,1) = xi(1,i)**2 - xi(1,i) ! xi**2 - xi
         Tm(i,2) = xi(2,i)**2 - xi(2,i) ! eta**2 - eta
         Tm(i,3) = xi(1,i)*xi(2,i)      ! xi * eta
      END DO
      Tm = MAT_INV(Tm, 3)

!     v = Inv(T) * z (reference config)
      v0(1) = Tm0(1,1)*xi0(3,1) + Tm0(1,2)*xi0(3,2) + Tm0(1,3)*xi0(3,3)
      v0(2) = Tm0(2,1)*xi0(3,1) + Tm0(2,2)*xi0(3,2) + Tm0(2,3)*xi0(3,3)
      v0(3) = Tm0(3,1)*xi0(3,1) + Tm0(3,2)*xi0(3,2) + Tm0(3,3)*xi0(3,3)

!     v = Inv(T) * z (current config)
      v(1) = Tm(1,1)*xi(3,1) + Tm(1,2)*xi(3,2) + Tm(1,3)*xi(3,3)
      v(2) = Tm(2,1)*xi(3,1) + Tm(2,2)*xi(3,2) + Tm(2,3)*xi(3,3)
      v(3) = Tm(3,1)*xi(3,1) + Tm(3,2)*xi(3,2) + Tm(3,3)*xi(3,3)

!     Bending strain, Eb = 2*(v0-v) = 2*(Tinv0*z0 - Tinv*z)
      Eb(:) = 2D0 * (v0(:) - v(:))

!     Now compute variation in bending strain
!     B1 bar
      B1b(:,1) = -Tm(:,1) * ((2D0*xi(1,1)-1D0)*v(1) + xi(2,1)*v(3))
      B1b(:,2) = -Tm(:,1) * ((2D0*xi(2,1)-1D0)*v(2) + xi(1,1)*v(3))

      B1b(:,3) = -Tm(:,2) * ((2D0*xi(1,2)-1D0)*v(1) + xi(2,2)*v(3))
      B1b(:,4) = -Tm(:,2) * ((2D0*xi(2,2)-1D0)*v(2) + xi(1,2)*v(3))

      B1b(:,5) = -Tm(:,3) * ((2D0*xi(1,3)-1D0)*v(1) + xi(2,3)*v(3))
      B1b(:,6) = -Tm(:,3) * ((2D0*xi(2,3)-1D0)*v(2) + xi(1,3)*v(3))

!     H1
      H1 = 0D0
      H1(1, 4: 6) =  gCnv(:,1)*adg(2,1)
      H1(2, 4: 6) =  gCnv(:,2)*adg(2,1)
      H1(3, 4: 6) =  gCnv(:,1)*adg(2,2)
      H1(4, 4: 6) =  gCnv(:,2)*adg(2,2)
      H1(5, 4: 6) =  gCnv(:,1)*adg(2,3)
      H1(6, 4: 6) =  gCnv(:,2)*adg(2,3)

      H1(1, 7: 9) = -gCnv(:,1)*adg(1,1)
      H1(2, 7: 9) = -gCnv(:,2)*adg(1,1)
      H1(3, 7: 9) = -gCnv(:,1)*adg(1,2)
      H1(4, 7: 9) = -gCnv(:,2)*adg(1,2)
      H1(5, 7: 9) = -gCnv(:,1)*adg(1,3)
      H1(6, 7: 9) = -gCnv(:,2)*adg(1,3)

      H1(1,10:12) =  gCnv(:,1)
      H1(2,10:12) =  gCnv(:,2)
      H1(3,13:15) =  gCnv(:,1)
      H1(4,13:15) =  gCnv(:,2)
      H1(5,16:18) =  gCnv(:,1)
      H1(6,16:18) =  gCnv(:,2)

!     H2
      H2 = 0D0
      H2( 1, 4) = -1D0
      H2( 2, 5) = -1D0
      H2( 3, 6) = -1D0
      H2( 4, 7) = -1D0
      H2( 5, 8) = -1D0
      H2( 6, 9) = -1D0
      H2( 7, 1) = -1D0
      H2( 8, 2) = -1D0
      H2( 9, 3) = -1D0

      H2( 1, 7) =  1D0
      H2( 2, 8) =  1D0
      H2( 3, 9) =  1D0
      H2( 4, 1) =  1D0
      H2( 5, 2) =  1D0
      H2( 6, 3) =  1D0
      H2( 7, 4) =  1D0
      H2( 8, 5) =  1D0
      H2( 9, 6) =  1D0

      H2(10, 4) = -1D0
      H2(11, 5) = -1D0
      H2(12, 6) = -1D0
      H2(13, 7) = -1D0
      H2(14, 8) = -1D0
      H2(15, 9) = -1D0
      H2(16, 1) = -1D0
      H2(17, 2) = -1D0
      H2(18, 3) = -1D0
      do i=10, 18
         H2(i,i) = 1D0
      end do

!     N matrix
      Nm = MAT_ID(3) - MAT_DYADPROD(nV, nV, 3)
!     M1, M2 matrices
      Mm(:,:,:) = 0D0
      DO i=1, 2
         Mm(1,2,i) = -gCov(3,i)
         Mm(1,3,i) =  gCov(2,i)
         Mm(2,3,i) = -gCov(1,i)

!        Skew-symmetric
         Mm(2,1,i) = -Mm(1,2,i)
         Mm(3,1,i) = -Mm(1,3,i)
         Mm(3,2,i) = -Mm(2,3,i)
      END DO

!     H3 matrix
      H3 = 0D0
      tmpA = MATMUL(Nm, Mm(:,:,1))
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

      tmpA = MATMUL(Nm, Mm(:,:,2))
      tmpA = -tmpA / Jac
      H3(1,7) = a(1,4)*tmpA(1,1) + a(2,4)*tmpA(2,1) + a(3,4)*tmpA(3,1)
      H3(1,8) = a(1,4)*tmpA(1,2) + a(2,4)*tmpA(2,2) + a(3,4)*tmpA(3,2)
      H3(1,9) = a(1,4)*tmpA(1,3) + a(2,4)*tmpA(2,3) + a(3,4)*tmpA(3,3)

      H3(2,7) = a(1,5)*tmpA(1,1) + a(2,5)*tmpA(2,1) + a(3,5)*tmpA(3,1)
      H3(2,8) = a(1,5)*tmpA(1,2) + a(2,5)*tmpA(2,2) + a(3,5)*tmpA(3,2)
      H3(2,9) = a(1,5)*tmpA(1,3) + a(2,5)*tmpA(2,3) + a(3,5)*tmpA(3,3)

      H3(3,7) = a(1,6)*tmpA(1,1) + a(2,6)*tmpA(2,1) + a(3,6)*tmpA(3,1)
      H3(3,8) = a(1,6)*tmpA(1,2) + a(2,6)*tmpA(2,2) + a(3,6)*tmpA(3,2)
      H3(3,9) = a(1,6)*tmpA(1,3) + a(2,6)*tmpA(2,3) + a(3,6)*tmpA(3,3)

      H3(1,10:12) = nV(:)
      H3(2,13:15) = nV(:)
      H3(3,16:18) = nV(:)

!     Variation in bending strain (Bb = -2*(B1b*H1*H2 + Tinv*H3*H2))
      H1H2 = MATMUL(H1, H2)
      Bb1  = MATMUL(B1b, H1H2)
      H3H2 = MATMUL(H3, H2)
      Bb1  = -2D0*(Bb1 + MATMUL(Tm, H3H2))

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
               aIi   = 1D0/SQRT(NORM(a0(:,i)))
               eI(:) = a0(:,i) * aIi
!              nI = nI0 = eI0 x n0 (reference config)
               nI(1) = eI(2)*nV0(3) - eI(3)*nV0(2)
               nI(2) = eI(3)*nV0(1) - eI(1)*nV0(3)
               nI(3) = eI(1)*nV0(2) - eI(2)*nV0(1)
               nInI  = MAT_DYADPROD(nI, nI, 3)
            ELSE
!              eI = aI/|aI| (current config)
               aIi   = 1D0/SQRT(NORM(a(:,i)))
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
            IF (BTEST(lM%sbc(i,e),bType_free)) THEN
!              E_I
               IF (.NOT.lFix(i)) THEN
                  tmpA = -Im + 2D0*eIeI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
!              E_P
               IF (.NOT.lFix(p)) THEN
                  tmpA = -2D0*(cI*Im - 2D0*cI*eIeI + aIi*eIaP)
                  Bb(:,:,p) = Bb(:,:,p) + MATMUL(Bb(:,:,j), tmpA)
               END IF
!              E_F
               IF (.NOT.lFix(f)) THEN
                  tmpA = 2D0*((1D0-cI)*Im - (1D0-2D0*cI)*eIeI -aIi*eIaP)
                  Bb(:,:,f) = Bb(:,:,f) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               Bb(:,:,j) = 0D0
            ELSE IF (BTEST(lM%sbc(i,e),bType_hing)) THEN
!              E_I
               IF (.NOT.lFix(i)) THEN
                  tmpA = -Im + 2D0*eIeI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               lFix(p)   = .TRUE.
               lFix(f)   = .TRUE.
               Bb(:,:,p) = 0D0
               Bb(:,:,f) = 0D0
               Bb(:,:,j) = 0D0
            ELSE IF (BTEST(lM%sbc(i,e),bType_fix)) THEN
               IF (.NOT.lFix(i)) THEN
                  tmpA = Im - 2D0*nInI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               lFix(p)   = .TRUE.
               lFix(f)   = .TRUE.
               Bb(:,:,f) = 0D0
               Bb(:,:,p) = 0D0
            ELSE IF (BTEST(lM%sbc(i,e),bType_symm)) THEN
               IF (.NOT.lFix(i)) THEN
                  tmpA = Im - 2D0*nInI
                  Bb(:,:,i) = Bb(:,:,i) + MATMUL(Bb(:,:,j), tmpA)
               END IF
               Bb(:,:,j) = 0D0
               tmpA(:,1) = eI(:)
               tmpA(:,2) = nV(:)
               tmpA(:,3) = nI(:)
               Bb(:,:,f) = MATMUL(Bb(:,:,f), tmpA)
               Bb(:,3,f) = 0D0
               Bb(:,:,p) = MATMUL(Bb(:,:,p), tmpA)
               Bb(:,3,p) = 0D0
            END IF
         END DO
      END IF

      RETURN
      END SUBROUTINE SHELLBENDTRI
!####################################################################
!     Set follower pressure load/net traction on shells
      SUBROUTINE SHELLFP (eNoN, w, N, Nx, dl, xl, tfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), dl(tDof,eNoN),
     2   xl(3,eNoN), tfl(eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, k, a, b
      REAL(KIND=8) :: T1, afl, wl, tfn, nV(3), gCov(3,2), gCnv(3,2),
     2   xc(3,eNoN), lKP(3)

      afl  = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

!     Get the current configuration and traction vector
      tfn = 0D0
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
      T1 = afl*wl*5D-1
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
      END SUBROUTINE SHELLFP
!####################################################################
      SUBROUTINE SHELLNRB (lM, g, eNoN, al, yl, dl, xl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: g, eNoN
      REAL(KIND=8), INTENT(IN) :: al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), xl(3,eNoN), bfl(3,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, k, l, a, b
      REAL(KIND=8) :: rho, dmp, elM, nu, ht, amd, afl, w, wb, Jac0, Jac,
     2   ud(3), fb(3), nV0(3), nV(3), gCov0(3,2), gCnv0(3,2), gCov(3,2),
     3   gCnv(3,2), N(eNoN), Nx(2,eNoN), Nxx(3,eNoN), x0(3,eNoN),
     4   xc(3,eNoN),Jm(2,2), eLoc(3,3), Qm(3,3), Dm(3,3), Em(3), Eb(3),
     5   DEm(3), DEb(3), Bm(3,3,eNoN), Bb(3,3,eNoN), K0(3,3), Kc(3,3),
     6   Nm(3,3), Mm(3,3,2), NmMm(3,3,2), KNmMm(3,3,2), DBm(3,3,eNoN),
     7   DBb(3,3,eNoN), BtDEm, BtDEb, NxSNx, BtDBm, BtDBb, T1

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp   = eq(cEq)%dmn(cDmn)%prop(damping)
      elM   = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu    = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
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
      N   = lM%N(:,g)
      Nx  = lM%Nx(:,:,g)
      Nxx = lM%Nxx(:,:,g)

!     Covariant and contravariant bases in reference config
      CALL GNNS(eNoN, Nx, x0, nV0, gCov0, gCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      CALL GNNS(eNoN, Nx, xc, nV, gCov, gCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac

!     Define local coordinates, Jacobian (J) and its inverse tensor
!     defined in Voigt notation such that Q*E = J^{-T} * E * J^{-1})
      eLoc(:,1) = gCov0(:,1)/SQRT(NORM(gCov0(:,1)))
      eLoc(:,3) = nV0(:)

      eLoc(1,2) = eLoc(2,3)*eLoc(3,1) - eLoc(3,3)*eLoc(2,1)
      eLoc(2,2) = eLoc(3,3)*eLoc(1,1) - eLoc(1,3)*eLoc(3,1)
      eLoc(3,2) = eLoc(1,3)*eLoc(2,1) - eLoc(2,3)*eLoc(1,1)

      Jm(1,1) = NORM(gCov0(:,1), eLoc(:,1))
      Jm(1,2) = NORM(gCov0(:,2), eLoc(:,1))
      Jm(2,1) = NORM(gCov0(:,1), eLoc(:,2))
      Jm(2,2) = NORM(gCov0(:,2), eLoc(:,2))

      Qm = 0D0
      Qm(1,1) = 1D0/Jm(1,1)**2
      Qm(2,1) = (-Jm(1,2)/(Jm(1,1)*Jm(2,2)))**2
      Qm(2,2) = 1D0/Jm(2,2)**2
      Qm(2,3) = -Jm(1,2)/(Jm(1,1)*Jm(2,2)**2)
      Qm(3,1) = -2D0*Jm(1,2)/(Jm(2,2)*Jm(1,1)**2)
      Qm(3,3) = 1D0/(Jm(1,1)*Jm(2,2))

!     Material tensor D, isotropic St Venant Kirchhoff model
      Dm = 0D0
      Dm(1,1) = 1D0
      Dm(1,2) = nu
      Dm(2,1) = nu
      Dm(2,2) = 1D0
      Dm(3,3) = 5D-1*(1D0-nu)
      Dm = elM/(1D0-nu**2) * Dm

!     Compute Qt * D * Q
      Dm = MATMUL(Dm, Qm)
      Dm = MATMUL(TRANSPOSE(Qm), Dm)

!     Membrane strain and its variation
!     Define membrane strain tensor (Em), Voigt notation
      Em = 0D0
      DO l=1, nsd
         Em(1) = Em(1) + gCov(l,1)*gCov(l,1) - gCov0(l,1)*gCov0(l,1)
         Em(2) = Em(2) + gCov(l,2)*gCov(l,2) - gCov0(l,2)*gCov0(l,2)
         Em(3) = Em(3) + gCov(l,1)*gCov(l,2) - gCov0(l,1)*gCov0(l,2)
      END DO
      Em(1) = 5D-1 * Em(1)
      Em(2) = 5D-1 * Em(2)

!     Define variation in membrane strain
      Bm = 0D0
      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*gCov(1,1)
         Bm(1,2,a) = Nx(1,a)*gCov(2,1)
         Bm(1,3,a) = Nx(1,a)*gCov(3,1)

         Bm(2,1,a) = Nx(2,a)*gCov(1,2)
         Bm(2,2,a) = Nx(2,a)*gCov(2,2)
         Bm(2,3,a) = Nx(2,a)*gCov(3,2)

         Bm(3,1,a) = Nx(2,a)*gCov(1,1) + Nx(1,a)*gCov(1,2)
         Bm(3,2,a) = Nx(2,a)*gCov(2,1) + Nx(1,a)*gCov(2,2)
         Bm(3,3,a) = Nx(2,a)*gCov(3,1) + Nx(1,a)*gCov(3,2)
      END DO

!     Bending strain and its variation
!     Second derivative of x, reference and current configurations
!     Compute second derivatives of x in reference and current configs
      K0 = 0D0
      Kc = 0D0
      DO a=1, eNoN
         K0(1,:) = K0(1,:) + Nxx(1,a)*x0(:,a)
         K0(2,:) = K0(2,:) + Nxx(2,a)*x0(:,a)
         K0(3,:) = K0(3,:) + Nxx(3,a)*x0(:,a)*2D0

         Kc(1,:) = Kc(1,:) + Nxx(1,a)*xc(:,a)
         Kc(2,:) = Kc(2,:) + Nxx(2,a)*xc(:,a)
         Kc(3,:) = Kc(3,:) + Nxx(3,a)*xc(:,a)*2D0
      END DO

!     Define bending strain tensor (Eb), Voigt notation
      Eb = 0D0
      DO l=1, nsd
         Eb(1) = Eb(1) + K0(l,1)*nV0(l) - Kc(l,1)*nV(l)
         Eb(2) = Eb(2) + K0(l,2)*nV0(l) - Kc(l,2)*nV(l)
         Eb(3) = Eb(3) + K0(l,3)*nV0(l) - Kc(l,3)*nV(l)
      END DO

!     Preliminaries to define variation in Eb
!     N matrix
      Nm = MAT_ID(3) - MAT_DYADPROD(nV, nV, 3)
      Nm = Nm / Jac
!     M1, M2 matrices
      Mm(:,:,:) = 0D0
      DO l=1, 2
         Mm(1,2,l) = -gCov(3,l)
         Mm(1,3,l) =  gCov(2,l)
         Mm(2,3,l) = -gCov(1,l)

!        Skew-symmetric
         Mm(2,1,l) = -Mm(1,2,l)
         Mm(3,1,l) = -Mm(1,3,l)
         Mm(3,2,l) = -Mm(2,3,l)

         NmMm(:,:,l)  = MATMUL(Nm, Mm(:,:,l))
         KNmMm(:,:,l) = MATMUL(Kc, NmMm(:,:,l))
      END DO

!     Define variation in bending strain tensor (Bb), Voigt notation
      Bb = 0D0
      DO a=1, eNoN
         Bb(1,:,a) = -Nxx(1,a)*nV(:)
         Bb(2,:,a) = -Nxx(2,a)*nV(:)
         Bb(3,:,a) = -Nxx(3,a)*nV(:)*2D0
      END DO

      DO a=1, eNoN
         Bb(:,:,a) = Bb(:,:,a) + Nx(1,a)*KNmMm(:,:,2) -
     2      Nx(2,a)*KNmMm(:,:,1)
      END DO

!     Contribution to residue and tangent matrices from membrane strain
!     D * Em
      DEm(1) = Dm(1,1)*Em(1) + Dm(1,2)*Em(2) + Dm(1,3)*Em(3)
      DEm(2) = Dm(2,1)*Em(1) + Dm(2,2)*Em(2) + Dm(2,3)*Em(3)
      DEm(3) = Dm(3,1)*Em(1) + Dm(3,2)*Em(2) + Dm(3,3)*Em(3)

!     D * Bm
      DO a=1, eNoN
         DBm(1,1,a) = Dm(1,1)*Bm(1,1,a) + Dm(1,2)*Bm(2,1,a) +
     2      Dm(1,3)*Bm(3,1,a)
         DBm(1,2,a) = Dm(1,1)*Bm(1,2,a) + Dm(1,2)*Bm(2,2,a) +
     2      Dm(1,3)*Bm(3,2,a)
         DBm(1,3,a) = Dm(1,1)*Bm(1,3,a) + Dm(1,2)*Bm(2,3,a) +
     2      Dm(1,3)*Bm(3,3,a)

         DBm(2,1,a) = Dm(2,1)*Bm(1,1,a) + Dm(2,2)*Bm(2,1,a) +
     2      Dm(2,3)*Bm(3,1,a)
         DBm(2,2,a) = Dm(2,1)*Bm(1,2,a) + Dm(2,2)*Bm(2,2,a) +
     2      Dm(2,3)*Bm(3,2,a)
         DBm(2,3,a) = Dm(2,1)*Bm(1,3,a) + Dm(2,2)*Bm(2,3,a) +
     2      Dm(2,3)*Bm(3,3,a)

         DBm(3,1,a) = Dm(3,1)*Bm(1,1,a) + Dm(3,2)*Bm(2,1,a) +
     2      Dm(3,3)*Bm(3,1,a)
         DBm(3,2,a) = Dm(3,1)*Bm(1,2,a) + Dm(3,2)*Bm(2,2,a) +
     2      Dm(3,3)*Bm(3,2,a)
         DBm(3,3,a) = Dm(3,1)*Bm(1,3,a) + Dm(3,2)*Bm(2,3,a) +
     2      Dm(3,3)*Bm(3,3,a)
      END DO

!     Contribution to residue and tangent matrices from bending strain
!     D * Eb
      DEb(1) = Dm(1,1)*Eb(1) + Dm(1,2)*Eb(2) + Dm(1,3)*Eb(3)
      DEb(2) = Dm(2,1)*Eb(1) + Dm(2,2)*Eb(2) + Dm(2,3)*Eb(3)
      DEb(3) = Dm(3,1)*Eb(1) + Dm(3,2)*Eb(2) + Dm(3,3)*Eb(3)

!     D * Bb
      DO a=1, eNoN
         DBb(1,1,a) = Dm(1,1)*Bb(1,1,a) + Dm(1,2)*Bb(2,1,a) +
     2      Dm(1,3)*Bb(3,1,a)
         DBb(1,2,a) = Dm(1,1)*Bb(1,2,a) + Dm(1,2)*Bb(2,2,a) +
     2      Dm(1,3)*Bb(3,2,a)
         DBb(1,3,a) = Dm(1,1)*Bb(1,3,a) + Dm(1,2)*Bb(2,3,a) +
     2      Dm(1,3)*Bb(3,3,a)

         DBb(2,1,a) = Dm(2,1)*Bb(1,1,a) + Dm(2,2)*Bb(2,1,a) +
     2      Dm(2,3)*Bb(3,1,a)
         DBb(2,2,a) = Dm(2,1)*Bb(1,2,a) + Dm(2,2)*Bb(2,2,a) +
     2      Dm(2,3)*Bb(3,2,a)
         DBb(2,3,a) = Dm(2,1)*Bb(1,3,a) + Dm(2,2)*Bb(2,3,a) +
     2      Dm(2,3)*Bb(3,3,a)

         DBb(3,1,a) = Dm(3,1)*Bb(1,1,a) + Dm(3,2)*Bb(2,1,a) +
     2      Dm(3,3)*Bb(3,1,a)
         DBb(3,2,a) = Dm(3,1)*Bb(1,2,a) + Dm(3,2)*Bb(2,2,a) +
     2      Dm(3,3)*Bb(3,2,a)
         DBb(3,3,a) = Dm(3,1)*Bb(1,3,a) + Dm(3,2)*Bb(2,3,a) +
     2      Dm(3,3)*Bb(3,3,a)
      END DO

!     Acceleration and mass damping at the integration point
      ud = -fb
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))
         ud(3) = ud(3) + N(a)*(rho*(al(k,a)-bfl(3,a)) + dmp*yl(k,a))
      END DO

!     Local residue
      w  = lM%w(g)*Jac0*ht
      wb = w*ht**2 / 12D0
      DO a=1, eNoN
         BtDEm = Bm(1,1,a)*DEm(1) + Bm(2,1,a)*DEm(2) + Bm(3,1,a)*DEm(3)
         BtDEb = Bb(1,1,a)*DEb(1) + Bb(2,1,a)*DEb(2) + Bb(3,1,a)*DEb(3)
         lR(1,a) = lR(1,a) + w*N(a)*ud(1) + w*BtDEm + wb*BtDEb

         BtDEm = Bm(1,2,a)*DEm(1) + Bm(2,2,a)*DEm(2) + Bm(3,2,a)*DEm(3)
         BtDEb = Bb(1,2,a)*DEb(1) + Bb(2,2,a)*DEb(2) + Bb(3,2,a)*DEb(3)
         lR(2,a) = lR(2,a) + w*N(a)*ud(2) + w*BtDEm + wb*BtDEb

         BtDEm = Bm(1,3,a)*DEm(1) + Bm(2,3,a)*DEm(2) + Bm(3,3,a)*DEm(3)
         BtDEb = Bb(1,3,a)*DEb(1) + Bb(2,3,a)*DEb(2) + Bb(3,3,a)*DEb(3)
         lR(3,a) = lR(3,a) + w*N(a)*ud(3) + w*BtDEm + wb*BtDEb
      END DO

!     Local stiffness
      DO b=1, eNoN
         DO a=1, eNoN
            NxSNx = ( Nx(1,a)*Nx(1,b)*DEm(1) + Nx(2,a)*Nx(2,b)*DEm(2) +
     2                Nx(1,a)*Nx(2,b)*DEm(3) + Nx(2,a)*Nx(1,b)*DEm(3) )
            T1 = w*( amd*N(a)*N(b) + afl*NxSNx )

            lK(1,a,b) = lK(1,a,b) + T1
            lK(dof+2,a,b) = lK(dof+2,a,b) + T1
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + T1

            BtDBm = Bm(1,1,a)*DBm(1,1,b) + Bm(2,1,a)*DBm(2,1,b) +
     2              Bm(3,1,a)*DBm(3,1,b)
            BtDBb = Bb(1,1,a)*DBb(1,1,b) + Bb(2,1,a)*DBb(2,1,b) +
     2              Bb(3,1,a)*DBb(3,1,b)
            lK(1,a,b) = lK(1,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,1,a)*DBm(1,2,b) + Bm(2,1,a)*DBm(2,2,b) +
     2              Bm(3,1,a)*DBm(3,2,b)
            BtDBb = Bb(1,1,a)*DBb(1,2,b) + Bb(2,1,a)*DBb(2,2,b) +
     2              Bb(3,1,a)*DBb(3,2,b)
            lK(2,a,b) = lK(2,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,1,a)*DBm(1,3,b) + Bm(2,1,a)*DBm(2,3,b) +
     2              Bm(3,1,a)*DBm(3,3,b)
            BtDBb = Bb(1,1,a)*DBb(1,3,b) + Bb(2,1,a)*DBb(2,3,b) +
     2              Bb(3,1,a)*DBb(3,3,b)
            lK(3,a,b) = lK(3,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,2,a)*DBm(1,1,b) + Bm(2,2,a)*DBm(2,1,b) +
     2              Bm(3,2,a)*DBm(3,1,b)
            BtDBb = Bb(1,2,a)*DBb(1,1,b) + Bb(2,2,a)*DBb(2,1,b) +
     2              Bb(3,2,a)*DBb(3,1,b)
            lK(dof+1,a,b) = lK(dof+1,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,2,a)*DBm(1,2,b) + Bm(2,2,a)*DBm(2,2,b) +
     2              Bm(3,2,a)*DBm(3,2,b)
            BtDBb = Bb(1,2,a)*DBb(1,2,b) + Bb(2,2,a)*DBb(2,2,b) +
     2              Bb(3,2,a)*DBb(3,2,b)
            lK(dof+2,a,b) = lK(dof+2,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,2,a)*DBm(1,3,b) + Bm(2,2,a)*DBm(2,3,b) +
     2              Bm(3,2,a)*DBm(3,3,b)
            BtDBb = Bb(1,2,a)*DBb(1,3,b) + Bb(2,2,a)*DBb(2,3,b) +
     2              Bb(3,2,a)*DBb(3,3,b)
            lK(dof+3,a,b) = lK(dof+3,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,3,a)*DBm(1,1,b) + Bm(2,3,a)*DBm(2,1,b) +
     2              Bm(3,3,a)*DBm(3,1,b)
            BtDBb = Bb(1,3,a)*DBb(1,1,b) + Bb(2,3,a)*DBb(2,1,b) +
     2              Bb(3,3,a)*DBb(3,1,b)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,3,a)*DBm(1,2,b) + Bm(2,3,a)*DBm(2,2,b) +
     2              Bm(3,3,a)*DBm(3,2,b)
            BtDBb = Bb(1,3,a)*DBb(1,2,b) + Bb(2,3,a)*DBb(2,2,b) +
     2              Bb(3,3,a)*DBb(3,2,b)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + afl*(w*BtDBm + wb*BtDBb)

            BtDBm = Bm(1,3,a)*DBm(1,3,b) + Bm(2,3,a)*DBm(2,3,b) +
     2              Bm(3,3,a)*DBm(3,3,b)
            BtDBb = Bb(1,3,a)*DBb(1,3,b) + Bb(2,3,a)*DBb(2,3,b) +
     2              Bb(3,3,a)*DBb(3,3,b)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + afl*(w*BtDBm + wb*BtDBb)
         END DO
      END DO

      RETURN
      END SUBROUTINE SHELLNRB
!####################################################################


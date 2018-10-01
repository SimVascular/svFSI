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
!     This routines is for solving nonlinear structural problem using
!     velocity based formulation.
!
!--------------------------------------------------------------------

      SUBROUTINE USTRUCT3D(eNoN, w, Je, N, Nx, al, yl, dl, fNl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, Je, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), fNl(3*nFn,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, k, l, a, b, iFn
      REAL(KIND=8) :: bf(3), am, af, aff, v(3), vd(3), vx(3,3), p, pd,
     2   px(3), fl(3,nFn), F(3,3), Jac, Fi(3,3), rho, beta, drho, dbeta,
     3   Siso(3,3), CCiso(3,3,3,3), tauM, tauC, rC, rCl, rM(3), Dm(6,6),
     4   Pdev(3,3), Bm(6,3,eNoN), DBm(6,3), NxFi(3,eNoN), VxFi(3,3), T1,
     5   T2, T3, PxFi(3), VxNx(3,eNoN), BtDB, NxNx, NxSNx, rMNx(eNoN),
     6   Ku

      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      bf(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      bf(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      aff     = af*af/am

!     {i,j,k} := velocity dofs; {l} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1
      l       = k + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0D0
      vd = -bf
      vx = 0D0
      p  = 0D0
      pd = 0D0
      px = 0D0
      fl = 0D0
      F  = MAT_ID(3)
      DO a=1, eNoN
         v(1)    = v(1)  + N(a)*yl(i,a)
         v(2)    = v(2)  + N(a)*yl(j,a)
         v(3)    = v(3)  + N(a)*yl(k,a)

         vd(1)   = vd(1) + N(a)*al(i,a)
         vd(2)   = vd(2) + N(a)*al(j,a)
         vd(3)   = vd(3) + N(a)*al(k,a)

         vx(1,1) = vx(1,1) + Nx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nx(2,a)*yl(i,a)
         vx(1,3) = vx(1,3) + Nx(3,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)
         vx(2,3) = vx(2,3) + Nx(3,a)*yl(j,a)
         vx(3,1) = vx(3,1) + Nx(1,a)*yl(k,a)
         vx(3,2) = vx(3,2) + Nx(2,a)*yl(k,a)
         vx(3,3) = vx(3,3) + Nx(3,a)*yl(k,a)

         p       = p     + N(a)*yl(l,a)
         pd      = pd    + N(a)*al(l,a)
         px(1)   = px(1) + Nx(1,a)*yl(l,a)
         px(2)   = px(2) + Nx(2,a)*yl(l,a)
         px(3)   = px(3) + Nx(3,a)*yl(l,a)

         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(1,3)  = F(1,3) + Nx(3,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)
         F(2,3)  = F(2,3) + Nx(3,a)*dl(j,a)
         F(3,1)  = F(3,1) + Nx(1,a)*dl(k,a)
         F(3,2)  = F(3,2) + Nx(2,a)*dl(k,a)
         F(3,3)  = F(3,3) + Nx(3,a)*dl(k,a)

         DO iFn=1, nFn
            b = (iFn-1)*3
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
            fl(3,iFn) = fl(3,iFn) + N(a)*fNl(b+3,a)
         END DO
      END DO
      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(stModel, p, rho, beta, drho, dbeta)

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     elasticity tensor (CCiso)
      CALL GETPK2CCdev(stModel, F, nFn, fl, Siso, CCiso)

!     Compute stabilization parameters
      CALL GETTAU(stModel, Je, tauM, tauC)

!     Deviatoric 1st Piola-Kirchhoff tensor (P)
      Pdev = MATMUL(F, Siso)

!     Convert to Voigt Notation
      Dm(1,1) = CCiso(1,1,1,1)
      Dm(1,2) = CCiso(1,1,2,2)
      Dm(1,3) = CCiso(1,1,3,3)
      Dm(1,4) = CCiso(1,1,1,2)
      Dm(1,5) = CCiso(1,1,2,3)
      Dm(1,6) = CCiso(1,1,3,1)

      Dm(2,2) = CCiso(2,2,2,2)
      Dm(2,3) = CCiso(2,2,3,3)
      Dm(2,4) = CCiso(2,2,1,2)
      Dm(2,5) = CCiso(2,2,2,3)
      Dm(2,6) = CCiso(2,2,3,1)

      Dm(3,3) = CCiso(3,3,3,3)
      Dm(3,4) = CCiso(3,3,1,2)
      Dm(3,5) = CCiso(3,3,2,3)
      Dm(3,6) = CCiso(3,3,3,1)

      Dm(4,4) = CCiso(1,2,1,2)
      Dm(4,5) = CCiso(1,2,2,3)
      Dm(4,6) = CCiso(1,2,3,1)

      Dm(5,5) = CCiso(2,3,2,3)
      Dm(5,6) = CCiso(2,3,3,1)

      Dm(6,6) = CCiso(3,1,3,1)

      DO a=2, 6
         DO b=1, a-1
            Dm(a,b) = Dm(b,a)
         END DO
      END DO

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)
         Bm(1,3,a) = Nx(1,a)*F(3,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)
         Bm(2,3,a) = Nx(2,a)*F(3,2)

         Bm(3,1,a) = Nx(3,a)*F(1,3)
         Bm(3,2,a) = Nx(3,a)*F(2,3)
         Bm(3,3,a) = Nx(3,a)*F(3,3)

         Bm(4,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(4,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
         Bm(4,3,a) = (Nx(1,a)*F(3,2) + F(3,1)*Nx(2,a))

         Bm(5,1,a) = (Nx(2,a)*F(1,3) + F(1,2)*Nx(3,a))
         Bm(5,2,a) = (Nx(2,a)*F(2,3) + F(2,2)*Nx(3,a))
         Bm(5,3,a) = (Nx(2,a)*F(3,3) + F(3,2)*Nx(3,a))

         Bm(6,1,a) = (Nx(3,a)*F(1,1) + F(1,3)*Nx(1,a))
         Bm(6,2,a) = (Nx(3,a)*F(2,1) + F(2,3)*Nx(1,a))
         Bm(6,3,a) = (Nx(3,a)*F(3,1) + F(3,3)*Nx(1,a))
      END DO

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1) + Nx(3,a)*Fi(3,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2) + Nx(3,a)*Fi(3,2)
         NxFi(3,a) = Nx(1,a)*Fi(1,3) + Nx(2,a)*Fi(2,3) + Nx(3,a)*Fi(3,3)
      END DO

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1) + vx(1,3)*Fi(3,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2) + vx(1,3)*Fi(3,2)
      VxFi(1,3) = vx(1,1)*Fi(1,3) + vx(1,2)*Fi(2,3) + vx(1,3)*Fi(3,3)

      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1) + vx(2,3)*Fi(3,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2) + vx(2,3)*Fi(3,2)
      VxFi(2,3) = vx(2,1)*Fi(1,3) + vx(2,2)*Fi(2,3) + vx(2,3)*Fi(3,3)

      VxFi(3,1) = vx(3,1)*Fi(1,1) + vx(3,2)*Fi(2,1) + vx(3,3)*Fi(3,1)
      VxFi(3,2) = vx(3,1)*Fi(1,2) + vx(3,2)*Fi(2,2) + vx(3,3)*Fi(3,2)
      VxFi(3,3) = vx(3,1)*Fi(1,3) + vx(3,2)*Fi(2,3) + vx(3,3)*Fi(3,3)

      PxFi(1)   = px(1)*Fi(1,1) + px(2)*Fi(2,1) + px(3)*Fi(3,1)
      PxFi(2)   = px(1)*Fi(1,2) + px(2)*Fi(2,2) + px(3)*Fi(3,2)
      PxFi(3)   = px(1)*Fi(1,3) + px(2)*Fi(2,3) + px(3)*Fi(3,3)

      rC    = beta*pd + VxFi(1,1) + VxFi(2,2) + VxFi(3,3)
      rCl   = -p + tauC*rC
      rM(1) = rho*vd(1) + PxFi(1)
      rM(2) = rho*vd(2) + PxFi(2)
      rM(3) = rho*vd(3) + PxFi(3)

!     Local residue
      DO a=1, eNoN
         T1 = Jac*rho*vd(1)*N(a)
         T2 = Pdev(1,1)*Nx(1,a) + Pdev(1,2)*Nx(2,a) + Pdev(1,3)*Nx(3,a)
         T3 = Jac*rCl*NxFi(1,a)
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(2)*N(a)
         T2 = Pdev(2,1)*Nx(1,a) + Pdev(2,2)*Nx(2,a) + Pdev(2,3)*Nx(3,a)
         T3 = Jac*rCl*NxFi(2,a)
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(3)*N(a)
         T2 = Pdev(3,1)*Nx(1,a) + Pdev(3,2)*Nx(2,a) + Pdev(3,3)*Nx(3,a)
         T3 = Jac*rCl*NxFi(3,a)
         lR(3,a) = lR(3,a) + w*(T1 + T2 + T3)

         rMNx(a) = rM(1)*NxFi(1,a) + rM(2)*NxFi(2,a) + rM(3)*NxFi(3,a)
         lR(4,a) = lR(4,a) + w*Jac*(N(a)*rC + tauM*rMNx(a))
      END DO

      DO a=1, eNoN
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a) +
     2               VxFi(3,1)*NxFi(3,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a) +
     2               VxFi(3,2)*NxFi(3,a)
         VxNx(3,a) = VxFi(1,3)*NxFi(1,a) + VxFi(2,3)*NxFi(2,a) +
     2               VxFi(3,3)*NxFi(3,a)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx  = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b) +
     2              NxFi(3,a)*NxFi(3,b)
            NxSNx = Nx(1,a)*Siso(1,1)*Nx(1,b) +
     2         Nx(1,a)*Siso(1,2)*Nx(2,b) + Nx(1,a)*Siso(1,3)*Nx(3,b) +
     3         Nx(2,a)*Siso(2,1)*Nx(1,b) + Nx(2,a)*Siso(2,2)*Nx(2,b) +
     4         Nx(2,a)*Siso(2,3)*Nx(3,b) + Nx(3,a)*Siso(3,1)*Nx(1,b) +
     5         Nx(3,a)*Siso(3,2)*Nx(2,b) + Nx(3,a)*Siso(3,3)*Nx(3,b)
            DBm   = MATMUL(Dm, Bm(:,:,b))

!           dM1_dV1 + af/am *dM_1/dU_1
            BtDB = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2             Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     3             Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            T1   = Jac*rho*vd(1)*N(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b)
            Ku   = aff*(T1 + T2 + BtDB + NxSNx)

            T1   = am*Jac*rho*N(a)*N(b)
            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b) + T1
            lK(1,a,b) = lK(1,a,b) + w*(T2 + Ku)

!           dM_1/dV_2 + af/am *dM_1/dU_2
            BtDB = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2             Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     3             Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            T1   = Jac*rho*vd(1)*N(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(2,b)
            lK(2,a,b) = lK(2,a,b) + w*(T2 + Ku)

!           dM_1/dV_3 + af/am *dM_1/dU_3
            BtDB = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2             Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     3             Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            T1   = Jac*rho*vd(1)*N(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(3,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(3,b) - NxFi(3,a)*NxFi(1,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(3,b)
            lK(3,a,b) = lK(3,a,b) + w*(T2 + Ku)

!           dM_2/dV_1 + af/am *dM_2/dU_1
            BtDB = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2             Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     3             Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            T1   = Jac*rho*vd(2)*N(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(1,b)
            lK(5,a,b) = lK(5,a,b) + w*(T2 + Ku)

!           dM_2/dV_2 + af/am *dM_2/dU_2
            BtDB = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2             Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     3             Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            T1   = Jac*rho*vd(2)*N(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(2,b)
            Ku   = aff*(T1 + T2 + BtDB + NxSNx)

            T1   = am*Jac*rho*N(a)*N(b)
            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(2,b) + T1
            lK(6,a,b) = lK(6,a,b) + w*(T2 + Ku)

!           dM_2/dV_3 + af/am *dM_2/dU_3
            BtDB = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2             Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     3             Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            T1   = Jac*rho*vd(2)*N(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(3,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(3,b) - NxFi(3,a)*NxFi(2,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(3,b)
            lK(7,a,b) = lK(7,a,b) + w*(T2 + Ku)

!           dM_3/dV_1 + af/am *dM_3/dU_1
            BtDB = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2             Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     3             Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            T1   = Jac*rho*vd(3)*N(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(3,a)*NxFi(1,b) - NxFi(1,a)*NxFi(3,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(3,a)*NxFi(1,b)
            lK(9,a,b) = lK(9,a,b) + w*(T2 + Ku)

!           dM_3/dV_2 + af/am *dM_3/dU_2
            BtDB = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2             Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     3             Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            T1   = Jac*rho*vd(3)*N(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(3,a)*NxFi(2,b) - NxFi(2,a)*NxFi(3,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(3,a)*NxFi(2,b)
            lK(10,a,b) = lK(10,a,b) + w*(T2 + Ku)

!           dM_3/dV_3 + af/am *dM_3/dU_3
            BtDB = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2             Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     3             Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            T1   = Jac*rho*vd(3)*N(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(3,b)
            Ku   = aff*(T1 + T2 + BtDB + NxSNx)

            T1   = am*Jac*rho*N(a)*N(b)
            T2   = af*Jac*tauC*rho*NxFi(3,a)*NxFi(3,b) + T1
            lK(11,a,b) = lK(11,a,b) + w*(T2 + Ku)

!           dC/dV_1 + af/am *dC/dU_1
            T1   = N(a)*(rC*NxFi(1,b) - VxNx(1,b))
            T2   = tauM*(rMNx(a)*NxFi(1,b) - rMNx(b)*NxFi(1,a))
            T3   = -tauM*NxNx*PxFi(1)
            Ku   = aff*(T1 + T2 + T3)

            T2   = (am*tauM*rho)*NxFi(1,a)*N(b) + af*N(a)*NxFi(1,b)
            lK(13,a,b) = lK(13,a,b) + w*Jac*(T2 + Ku)

!           dC/dV_2 + af/am *dC/dU_2
            T1   = N(a)*(rC*NxFi(2,b) - VxNx(2,b))
            T2   = tauM*(rMNx(a)*NxFi(2,b) - rMNx(b)*NxFi(2,a))
            T3   = -tauM*NxNx*PxFi(2)
            Ku   = aff*(T1 + T2 + T3)

            T2 = (am*tauM*rho)*NxFi(2,a)*N(b) + af*N(a)*NxFi(2,b)
            lK(14,a,b) = lK(14,a,b) + w*Jac*(T2 + Ku)

!           dC/dV_3 + af/am *dC/dU_3
            T1   = N(a)*(rC*NxFi(3,b) - VxNx(3,b))
            T2   = tauM*(rMNx(a)*NxFi(3,b) - rMNx(b)*NxFi(3,a))
            T3   = -tauM*NxNx*PxFi(3)
            Ku   = aff*(T1 + T2 + T3)

            T2 = (am*tauM*rho)*NxFi(3,a)*N(b) + af*N(a)*NxFi(3,b)
            lK(15,a,b) = lK(15,a,b) + w*Jac*(T2 + Ku)

!           dM_1/dP
            T1 = am*tauC*beta + af*(tauC*dbeta*pd - 1.0D0)
            T2 = T1*NxFi(1,a)*N(b) + af*drho*vd(1)*N(a)*N(b)
            lK(4,a,b) = lK(4,a,b) + w*Jac*T2

!           dM_2/dP
            T2 = T1*NxFi(2,a)*N(b) + af*drho*vd(2)*N(a)*N(b)
            lK(8,a,b) = lK(8,a,b) + w*Jac*T2

!           dM_3/dP
            T2 = T1*NxFi(3,a)*N(b) + af*drho*vd(3)*N(a)*N(b)
            lK(12,a,b) = lK(12,a,b) + w*Jac*T2

!           dC/dP
            T1 = (am*beta + af*dbeta*pd)*N(a)*N(b)
            T2 = NxFi(1,a)*vd(1) + NxFi(2,a)*vd(2) + NxFi(3,a)*vd(3)
            T3 = T1 + af*tauM*(NxNx + drho*T2*N(b))
            lK(16,a,b) = lK(16,a,b) + w*Jac*T3
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE USTRUCT2D(eNoN, w, Je, N, Nx, al, yl, dl, fNl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, Je, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), fNl(2*nFn,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, k, a, b, iFn
      REAL(KIND=8) :: bf(2), am, af, aff, v(2), vd(2), vx(2,2), p, pd,
     2   px(2), fl(2,nFn), F(2,2), Jac, Fi(2,2), rho, beta, drho, dbeta,
     3   Siso(2,2), CCiso(2,2,2,2), tauM, tauC, rC, rCl, rM(2), Dm(3,3),
     4   Pdev(2,2), Bm(3,2,eNoN), DBm(3,2), NxFi(2,eNoN), VxFi(2,2), T1,
     5   T2, T3, PxFi(2), VxNx(2,eNoN), BtDB, NxNx, NxSNx, rMNx(eNoN),
     6   Ku

      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      bf(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      aff     = af*af/am

!     {i,j} := velocity dofs; {k} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0D0
      vd = -bf
      vx = 0D0
      p  = 0D0
      pd = 0D0
      px = 0D0
      fl = 0D0
      F  = MAT_ID(2)
      DO a=1, eNoN
         v(1)    = v(1)  + N(a)*yl(i,a)
         v(2)    = v(2)  + N(a)*yl(j,a)

         vd(1)   = vd(1) + N(a)*al(i,a)
         vd(2)   = vd(2) + N(a)*al(j,a)

         vx(1,1) = vx(1,1) + Nx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nx(2,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)

         p       = p     + N(a)*yl(k,a)
         pd      = pd    + N(a)*al(k,a)
         px(1)   = px(1) + Nx(1,a)*yl(k,a)
         px(2)   = px(2) + Nx(2,a)*yl(k,a)

         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)

         DO iFn=1, nFn
            b = (iFn-1)*2
            fl(1,iFn) = fl(1,iFn) + N(a)*fNl(b+1,a)
            fl(2,iFn) = fl(2,iFn) + N(a)*fNl(b+2,a)
         END DO
      END DO
      Jac = F(1,1)*F(2,2) - F(1,2)*F(2,1)
      Fi  = MAT_INV(F, 2)

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(stModel, p, rho, beta, drho, dbeta)

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     elasticity tensor (CCiso)
      CALL GETPK2CCdev(stModel, F, nFn, fl, Siso, CCiso)

!     Compute stabilization parameters
      CALL GETTAU(stModel, Je, tauM, tauC)

!     Deviatoric 1st Piola-Kirchhoff tensor (P)
      Pdev = MATMUL(F, Siso)

!     Convert to Voigt Notation
      Dm(1,1) = CCiso(1,1,1,1)
      Dm(1,2) = CCiso(1,1,2,2)
      Dm(1,3) = CCiso(1,1,1,2)

      Dm(2,2) = CCiso(2,2,2,2)
      Dm(2,3) = CCiso(2,2,1,2)

      Dm(3,3) = CCiso(1,2,1,2)

      Dm(2,1) = Dm(1,2)
      Dm(3,1) = Dm(1,3)
      Dm(3,2) = Dm(2,3)

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)

         Bm(3,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(3,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
      END DO

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2)
      END DO

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2)

      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2)

      PxFi(1)   = px(1)*Fi(1,1) + px(2)*Fi(2,1)
      PxFi(2)   = px(1)*Fi(1,2) + px(2)*Fi(2,2)

      rC    = beta*pd + VxFi(1,1) + VxFi(2,2)
      rCl   = -p + tauC*rC
      rM(1) = rho*vd(1) + PxFi(1)
      rM(2) = rho*vd(2) + PxFi(2)

!     Local residue
      DO a=1, eNoN
         T1 = Jac*rho*vd(1)*N(a)
         T2 = Pdev(1,1)*Nx(1,a) + Pdev(1,2)*Nx(2,a)
         T3 = Jac*rCl*NxFi(1,a)
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(2)*N(a)
         T2 = Pdev(2,1)*Nx(1,a) + Pdev(2,2)*Nx(2,a)
         T3 = Jac*rCl*NxFi(2,a)
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)

         rMNx(a) = rM(1)*NxFi(1,a) + rM(2)*NxFi(2,a)
         lR(3,a) = lR(3,a) + w*Jac*(N(a)*rC + tauM*rMNx(a))
      END DO

      DO a=1, eNoN
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx  = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b)
            NxSNx = Nx(1,a)*Siso(1,1)*Nx(1,b)
     2            + Nx(1,a)*Siso(1,2)*Nx(2,b)
     3            + Nx(2,a)*Siso(2,1)*Nx(1,b)
     4            + Nx(2,a)*Siso(2,2)*Nx(2,b)
            DBm(1,1) = Dm(1,1)*Bm(1,1,b) + Dm(1,2)*Bm(2,1,b) +
     2         Dm(1,3)*Bm(3,1,b)
            DBm(1,2) = Dm(1,1)*Bm(1,2,b) + Dm(1,2)*Bm(2,2,b) +
     2         Dm(1,3)*Bm(3,2,b)

            DBm(2,1) = Dm(2,1)*Bm(1,1,b) + Dm(2,2)*Bm(2,1,b) +
     2         Dm(2,3)*Bm(3,1,b)
            DBm(2,2) = Dm(2,1)*Bm(1,2,b) + Dm(2,2)*Bm(2,2,b) +
     2         Dm(2,3)*Bm(3,2,b)

            DBm(3,1) = Dm(3,1)*Bm(1,1,b) + Dm(3,2)*Bm(2,1,b) +
     2         Dm(3,3)*Bm(3,1,b)
            DBm(3,2) = Dm(3,1)*Bm(1,2,b) + Dm(3,2)*Bm(2,2,b) +
     2         Dm(3,3)*Bm(3,2,b)

!           dM_1/dV_1 + af/am *dM_1/dU_1
            BtDB = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1)
            T1   = Jac*rho*vd(1)*N(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b)
            Ku   = aff*(T1 + T2 + BtDB + NxSNx)

            T1   = am*Jac*rho*N(a)*N(b)
            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b) + T1
            lK(1,a,b) = lK(1,a,b) + w*(T2 + Ku)

!           dM_1/dV_2 + af/am *dM_1/dU_2
            BtDB = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2             Bm(3,1,a)*DBm(3,2)
            T1   = Jac*rho*vd(1)*N(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(2,b)
            lK(2,a,b) = lK(2,a,b) + w*(T2 + Ku)

!           dM_2/dV_1 + af/am *dM_2/dU_1
            BtDB = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2             Bm(3,2,a)*DBm(3,1)
            T1   = Jac*rho*vd(2)*N(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))
            Ku   = aff*(T1 + T2 + T3 + BtDB)

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(1,b)
            lK(4,a,b) = lK(4,a,b) + w*(T2 + Ku)

!           dM_2/dV_2 + af/am *dM_2/dU_2
            BtDB = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2             Bm(3,2,a)*DBm(3,2)
            T1   = Jac*rho*vd(2)*N(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(2,b)
            Ku   = aff*(T1 + T2 + BtDB + NxSNx)

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(2,b) + T1
            lK(5,a,b) = lK(5,a,b) + w*(T2 + Ku)

!           dC/dV_1 + af/am *dC/dU_1
            T1   = N(a)*(rC*NxFi(1,b) - VxNx(1,b))
            T2   = tauM*(rMNx(a)*NxFi(1,b) - rMNx(b)*NxFi(1,a))
            T3   = -tauM*NxNx*PxFi(1)
            Ku   = aff*(T1 + T2 + T3)

            T2   = (am*tauM*rho)*NxFi(1,a)*N(b) + af*N(a)*NxFi(1,b)
            lK(7,a,b) = lK(7,a,b) + w*Jac*(T2 + Ku)

!           dC/dV_2 + af/am *dC/dU_2
            T1   = N(a)*(rC*NxFi(2,b) - VxNx(2,b))
            T2   = tauM*(rMNx(a)*NxFi(2,b) - rMNx(b)*NxFi(2,a))
            T3   = -tauM*NxNx*PxFi(2)
            Ku   = aff*(T1 + T2 + T3)

            T2 = (am*tauM*rho)*NxFi(2,a)*N(b) + af*N(a)*NxFi(2,b)
            lK(8,a,b) = lK(8,a,b) + w*Jac*(T2 + Ku)

!           dM_1/dP
            T1 = am*tauC*beta + af*(tauC*dbeta*pd - 1.0D0)
            T2 = T1*NxFi(1,a)*N(b) + af*drho*vd(1)*N(a)*N(b)
            lK(3,a,b) = lK(3,a,b) + w*Jac*T2

!           dM_2/dP
            T2 = T1*NxFi(2,a)*N(b) + af*drho*vd(2)*N(a)*N(b)
            lK(6,a,b) = lK(6,a,b) + w*Jac*T2

!           dC/dP
            T1 = (am*beta + af*dbeta*pd)*N(a)*N(b)
            T2 = NxFi(1,a)*vd(1) + NxFi(2,a)*vd(2)
            T3 = T1 + af*tauM*(NxNx + drho*T2*N(b))
            lK(9,a,b) = lK(9,a,b) + w*Jac*T3
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT2D
!####################################################################
      SUBROUTINE BUSTRUCTNEU(lFa, hg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: hg(tnNo), Dg(tDof,tnNo)

      INTEGER :: a, b, e, g, iM, Ac, Ec, eNoN, eNoNb, cPhys
      REAL(KIND=8) :: w, Jac, xp(nsd), xi(nsd), nV(nsd), ksix(nsd,nsd),
     2   rt
      LOGICAL :: flag

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:), Nx(:,:),
     2   dl(:,:), hl(:), lR(:,:), lK(:,:,:)

      iM    = lFa%iM
      eNoNb = lFa%eNoN
      eNoN  = msh(iM)%eNoN

      ALLOCATE(xl(nsd,eNoN), N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN),
     2   ptr(eNoN), dl(tDof,eNoN), hl(eNoN),lR(dof,eNoN),
     3   lK(dof*dof,eNoN,eNoN))

      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_ustruct) CYCLE

         lR    = 0D0
         lK    = 0D0
         DO a=1, eNoN
            Ac      = msh(iM)%IEN(a,Ec)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            dl(:,a) = Dg(:,Ac)
            hl(a)   = hg(Ac)
         END DO

         DO g=1, lFa%nG
            xp = 0D0
            DO a=1, eNoNb
               Ac = lFa%IEN(a,e)
               xp = xp + x(:,Ac)*lFa%N(a,g)
            END DO

            xi = msh(iM)%xi(:,1)
            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, flag)
!           If Newton's method fails
            IF (.NOT.flag) err = "Newton method failed to converge. "//
     2         "(BUSTRUCTNEU)"

            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)
            a = 0
            DO b=1, eNoN
               IF (N(b).GT.-1E-4 .AND. N(b).LT.1.0001D0) a = a + 1
            END DO
            IF (a .NE. eNoN) err =
     2         "Error in computing face derivatives (BUSTRUCTNEU)"

            IF (g.EQ.1 .OR. .NOT.msh(iM)%lShpF)
     2         CALL GNN(eNoN, Nxi, xl, Nx, rt, ksix)

            CALL GNNB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
            nV  = nV / Jac
            w   = lFa%w(g)*Jac

            IF (nsd .EQ. 3) THEN
               CALL BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
            ELSE
               CALL BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
            END IF
         END DO

#ifdef WITH_TRILINOS
         IF (useTrilinosAssemAndLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO

      RETURN
      END SUBROUTINE BUSTRUCTNEU
!--------------------------------------------------------------------
      SUBROUTINE BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(3)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, k, a, b
      REAL(KIND=8) :: af, Jac, h, wl, afl, T1, F(3,3), Fi(3,3), nFi(3),
     2   NxFi(3,eNoN)

      af = eq(cEq)%af*eq(cEq)%gam*dt
      af = af * af / eq(cEq)%am
      i  = eq(cEq)%s
      j  = i + 1
      k  = j + 1

      h  = 0D0
      F  = MAT_ID(3)
      DO a=1, eNoN
         h       = h      + N(a)*hl(a)

         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(1,3)  = F(1,3) + Nx(3,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)
         F(2,3)  = F(2,3) + Nx(3,a)*dl(j,a)
         F(3,1)  = F(3,1) + Nx(1,a)*dl(k,a)
         F(3,2)  = F(3,2) + Nx(2,a)*dl(k,a)
         F(3,3)  = F(3,3) + Nx(3,a)*dl(k,a)
      END DO
      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1) + Nx(3,a)*Fi(3,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2) + Nx(3,a)*Fi(3,2)
         NxFi(3,a) = Nx(1,a)*Fi(1,3) + Nx(2,a)*Fi(2,3) + Nx(3,a)*Fi(3,3)
      END DO

      nFi(1) = nV(1)*Fi(1,1) + nV(2)*Fi(2,1) + nV(3)*Fi(3,1)
      nFi(2) = nV(1)*Fi(1,2) + nV(2)*Fi(2,2) + nV(3)*Fi(3,2)
      nFi(3) = nV(1)*Fi(1,3) + nV(2)*Fi(2,3) + nV(3)*Fi(3,3)

      wl  = w*Jac*h
      afl = af*wl
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)
         lR(3,a) = lR(3,a) - wl*N(a)*nFi(3)

         DO b=1, eNoN
            T1 = nFi(1)*NxFi(2,b) - nFi(2)*NxFi(1,b)
            lK(2,a,b) = lK(2,a,b) - afl*N(a)*T1

            T1 = nFi(1)*NxFi(3,b) - nFi(3)*NxFi(1,b)
            lK(3,a,b) = lK(3,a,b) - afl*N(a)*T1

            T1 = nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b)
            lK(5,a,b) = lK(5,a,b) - afl*N(a)*T1

            T1 = nFi(2)*NxFi(3,b) - nFi(3)*NxFi(2,b)
            lK(7,a,b) = lK(7,a,b) - afl*N(a)*T1

            T1 = nFi(3)*NxFi(1,b) - nFi(1)*NxFi(3,b)
            lK(9,a,b) = lK(9,a,b) - afl*N(a)*T1

            T1 = nFi(3)*NxFi(2,b) - nFi(2)*NxFi(3,b)
            lK(10,a,b) = lK(10,a,b) - afl*N(a)*T1
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(2)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: i, j, a, b
      REAL(KIND=8) :: af, Jac, h, wl, afl, T1, F(2,2), Fi(2,2), nFi(2),
     2   NxFi(2,eNoN)

      af = eq(cEq)%af*eq(cEq)%gam*dt
      af = af * af / eq(cEq)%am
      i  = eq(cEq)%s
      j  = i + 1

      h  = 0D0
      F  = MAT_ID(2)
      DO a=1, eNoN
         h       = h      + N(a)*hl(a)

         F(1,1)  = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nx(2,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nx(2,a)*dl(j,a)
      END DO
      Jac = F(1,1)*F(2,2) - F(1,2)*F(2,1)
      Fi  = MAT_INV(F, 2)

      DO a=1, eNoN
         NxFi(1,a) = Nx(1,a)*Fi(1,1) + Nx(2,a)*Fi(2,1)
         NxFi(2,a) = Nx(1,a)*Fi(1,2) + Nx(2,a)*Fi(2,2)
      END DO

      nFi(1) = nV(1)*Fi(1,1) + nV(2)*Fi(2,1)
      nFi(2) = nV(1)*Fi(1,2) + nV(2)*Fi(2,2)

      wl  = w*Jac*h
      afl = af*wl
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)

         DO b=1, eNoN
            T1 = nFi(1)*NxFi(2,b) - nFi(2)*NxFi(1,b)
            lK(2,a,b) = lK(2,a,b) - afl*N(a)*T1

            T1 = nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b)
            lK(4,a,b) = lK(4,a,b) - afl*N(a)*T1
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT2D
!####################################################################

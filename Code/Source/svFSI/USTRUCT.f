!--------------------------------------------------------------------
!
!     This routines is for solving nonlinear structural problem using
!     velocity based formulation
!
!--------------------------------------------------------------------
!####################################################################
      SUBROUTINE USTRUCT3D(eNoN, w, Je, N, Nx, al, yl, dl, adl, fNl,
     2   lRd, lKd, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, Je, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), adl(3,eNoN),
     3   fNl(nFn*3,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lRd(3,eNoN), lR(dof,eNoN),
     2   lKd(3*dof,eNoN,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: a, b, i, j, k, l, iFn
      REAL(KIND=8) :: rho, drho, bf(3), am, af, v(3), vd(3), vx(3,3), p,
     2   pd, px(3), ud(3), fl(3,nFn), F(3,3), Jac, Fi(3,3), beta, dbeta,
     3   Siso(3,3), CCiso(3,3,3,3), tauM, tauC, rC, rCl, rM(3), Dm(6,6),
     4   Pdev(3,3), Bm(6,3,eNoN), DBm(6,3), NxFi(3,eNoN), VxFi(3,3), T1,
     5   T2, T3, PxFi(3), VxNx(3,eNoN), NxNx, NxSNx, NxPx, BtDB, drb,
     6   rMNx(eNoN)

      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      bf(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      bf(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt

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
      ud = 0D0
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
         vx(2,1) = vx(2,1) + Nx(2,a)*yl(i,a)
         vx(3,1) = vx(3,1) + Nx(3,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)
         vx(3,2) = vx(3,2) + Nx(3,a)*yl(j,a)
         vx(1,3) = vx(1,3) + Nx(1,a)*yl(k,a)
         vx(2,3) = vx(2,3) + Nx(2,a)*yl(k,a)
         vx(3,3) = vx(3,3) + Nx(3,a)*yl(k,a)

         p       = p     + N(a)*yl(l,a)
         pd      = pd    + N(a)*al(l,a)
         px(1)   = px(1) + Nx(1,a)*yl(l,a)
         px(2)   = px(2) + Nx(2,a)*yl(l,a)
         px(3)   = px(3) + Nx(3,a)*yl(l,a)

         ud(1)   = ud(1) + N(a)*adl(1,a)
         ud(2)   = ud(2) + N(a)*adl(2,a)
         ud(3)   = ud(3) + N(a)*adl(3,a)

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
      drb = drho*beta + rho*dbeta

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     elasticity tensor (CCiso)
      CALL GETPK2CCdev(stModel, F, fl, Siso, CCiso)

!     Compute stabilization parameters
      CALL GETTAU(stModel, Je, F, fl, tauM, tauC)

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

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(2,1)*Fi(2,1) + vx(3,1)*Fi(3,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(2,1)*Fi(2,2) + vx(3,1)*Fi(3,2)
      VxFi(1,3) = vx(1,1)*Fi(1,3) + vx(2,1)*Fi(2,3) + vx(3,1)*Fi(3,3)

      VxFi(2,1) = vx(1,2)*Fi(1,1) + vx(2,2)*Fi(2,1) + vx(3,2)*Fi(3,1)
      VxFi(2,2) = vx(1,2)*Fi(1,2) + vx(2,2)*Fi(2,2) + vx(3,2)*Fi(3,2)
      VxFi(2,3) = vx(1,2)*Fi(1,3) + vx(2,2)*Fi(2,3) + vx(3,2)*Fi(3,3)

      VxFi(3,1) = vx(1,3)*Fi(1,1) + vx(2,3)*Fi(2,1) + vx(3,3)*Fi(3,1)
      VxFi(3,2) = vx(1,3)*Fi(1,2) + vx(2,3)*Fi(2,2) + vx(3,3)*Fi(3,2)
      VxFi(3,3) = vx(1,3)*Fi(1,3) + vx(2,3)*Fi(2,3) + vx(3,3)*Fi(3,3)

      PxFi(1)   = px(1)*Fi(1,1) + px(2)*Fi(2,1) + px(3)*Fi(3,1)
      PxFi(2)   = px(1)*Fi(1,2) + px(2)*Fi(2,2) + px(3)*Fi(3,2)
      PxFi(3)   = px(1)*Fi(1,3) + px(2)*Fi(2,3) + px(3)*Fi(3,3)

      DO a=1, eNoN
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a) +
     2               VxFi(3,1)*NxFi(3,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a) +
     2               VxFi(3,2)*NxFi(3,a)
         VxNx(3,a) = VxFi(1,3)*NxFi(1,a) + VxFi(2,3)*NxFi(2,a) +
     2               VxFi(3,3)*NxFi(3,a)
      END DO

      rC    = beta*pd + VxFi(1,1) + VxFi(2,2) + VxFi(3,3)
      rCl   = -p + tauC*rho*rC
      rM(1) = vd(1) + PxFi(1)/rho
      rM(2) = vd(2) + PxFi(2)/rho
      rM(3) = vd(3) + PxFi(3)/rho

!     Local residue
      DO a=1, eNoN
         lRd(1,a) = lRd(1,a) + w*N(a)*Jac*(ud(1) - v(1))
         lRd(2,a) = lRd(2,a) + w*N(a)*Jac*(ud(2) - v(2))
         lRd(3,a) = lRd(3,a) + w*N(a)*Jac*(ud(3) - v(3))

         T1 = N(a)*rho*vd(1)*Jac
         T2 = Nx(1,a)*Pdev(1,1) + Nx(2,a)*Pdev(1,2) + Nx(3,a)*Pdev(1,3)
         T3 = NxFi(1,a)*rCl*Jac
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = N(a)*rho*vd(2)*Jac
         T2 = Nx(1,a)*Pdev(2,1) + Nx(2,a)*Pdev(2,2) + Nx(3,a)*Pdev(2,3)
         T3 = NxFi(2,a)*rCl*Jac
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)

         T1 = N(a)*rho*vd(3)*Jac
         T2 = Nx(1,a)*Pdev(3,1) + Nx(2,a)*Pdev(3,2) + Nx(3,a)*Pdev(3,3)
         T3 = NxFi(3,a)*rCl*Jac
         lR(3,a) = lR(3,a) + w*(T1 + T2 + T3)

         rMNx(a) = rM(1)*NxFi(1,a) + rM(2)*NxFi(2,a) + rM(3)*NxFi(3,a)
         lR(4,a) = lR(4,a) + w*Jac*(N(a)*rC + tauM*rMNx(a))
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx  = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b) +
     2              NxFi(3,a)*NxFi(3,b)
            NxPx  = NxFi(1,a)*PxFi(1) + NxFi(2,a)*PxFi(2) +
     2              NxFi(3,a)*PxFi(3)
            NxSNx = Nx(1,a)*Siso(1,1)*Nx(1,b) +
     2         Nx(1,a)*Siso(1,2)*Nx(2,b) + Nx(1,a)*Siso(1,3)*Nx(3,b) +
     3         Nx(2,a)*Siso(2,1)*Nx(1,b) + Nx(2,a)*Siso(2,2)*Nx(2,b) +
     4         Nx(2,a)*Siso(2,3)*Nx(3,b) + Nx(3,a)*Siso(3,1)*Nx(1,b) +
     5         Nx(3,a)*Siso(3,2)*Nx(2,b) + Nx(3,a)*Siso(3,3)*Nx(3,b)
            DBm   = MATMUL(Dm, Bm(:,:,b))

!           dM_1/dU_1
            BtDB  = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     3              Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            T1    = N(a)*rho*vd(1)*NxFi(1,b)*Jac
            T2    = -NxFi(1,a)*tauC*rho*VxNx(1,b)*Jac
            lKd(1,a,b) = lKd(1,a,b) + w*af*(BtDB + T1 + T2 + NxSNx)

!           dM_1/dU_2
            BtDB  = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2              Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     3              Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            T1    = N(a)*rho*vd(1)*NxFi(2,b)*Jac
            T2    = -NxFi(1,a)*tauC*rho*VxNx(2,b)*Jac
            T3    = (NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))*rCl*Jac
            lKd(2,a,b) = lKd(2,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_1/dU_3
            BtDB  = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2              Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     3              Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            T1    = N(a)*rho*vd(1)*NxFi(3,b)*Jac
            T2    = -NxFi(1,a)*tauC*rho*VxNx(3,b)*Jac
            T3    = (NxFi(1,a)*NxFi(3,b) - NxFi(3,a)*NxFi(1,b))*rCl*Jac
            lKd(3,a,b) = lKd(3,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_2/dU_1
            BtDB  = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2              Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     3              Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            T1    = N(a)*rho*vd(2)*NxFi(1,b)*Jac
            T2    = -NxFi(2,a)*tauC*rho*VxNx(1,b)*Jac
            T3    = (NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))*rCl*Jac
            lKd(4,a,b) = lKd(4,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_2/dU_2
            BtDB  = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2              Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     3              Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            T1    = N(a)*rho*vd(2)*NxFi(2,b)*Jac
            T2    = -NxFi(2,a)*tauC*rho*VxNx(2,b)*Jac
            lKd(5,a,b) = lKd(5,a,b) + w*af*(BtDB + T1 + T2 + NxSNx)

!           dM_2/dU_3
            BtDB  = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2              Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     3              Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            T1    = N(a)*rho*vd(2)*NxFi(3,b)*Jac
            T2    = -NxFi(2,a)*tauC*rho*VxNx(3,b)*Jac
            T3    = (NxFi(2,a)*NxFi(3,b) - NxFi(3,a)*NxFi(2,b))*rCl*Jac
            lKd(6,a,b) = lKd(6,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_3/dU_1
            BtDB  = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2              Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     3              Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            T1    = N(a)*rho*vd(3)*NxFi(1,b)*Jac
            T2    = -NxFi(3,a)*tauC*rho*VxNx(1,b)*Jac
            T3    = (NxFi(3,a)*NxFi(1,b) - NxFi(1,a)*NxFi(3,b))*rCl*Jac
            lKd(7,a,b) = lKd(7,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_3/dU_2
            BtDB  = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2              Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     3              Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            T1    = N(a)*rho*vd(3)*NxFi(2,b)*Jac
            T2    = -NxFi(3,a)*tauC*rho*VxNx(2,b)*Jac
            T3    = (NxFi(3,a)*NxFi(2,b) - NxFi(2,a)*NxFi(3,b))*rCl*Jac
            lKd(8,a,b) = lKd(8,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_3/dU_3
            BtDB  = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2              Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     3              Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            T1    = N(a)*rho*vd(3)*NxFi(3,b)*Jac
            T2    = -NxFi(3,a)*tauC*rho*VxNx(3,b)*Jac
            lKd(9,a,b) = lKd(9,a,b) + w*af*(BtDB + T1 + T2 + NxSNx)

!           dC/dU_1
            T1    = N(a)*(rC*NxFi(1,b) - VxNx(1,b))
            T2    = tauM*(rMNx(a)*NxFi(1,b) - NxFi(1,a)*rMNx(b))
            T3    = -tauM*NxNx*PxFi(1)/rho
            lKd(10,a,b) = lKd(10,a,b) + w*af*Jac*(T1 + T2 + T3)

!           dC/dU_2
            T1    = N(a)*(rC*NxFi(2,b) - VxNx(2,b))
            T2    = tauM*(rMNx(a)*NxFi(2,b) - NxFi(2,a)*rMNx(b))
            T3    = -tauM*NxNx*PxFi(2)/rho
            lKd(11,a,b) = lKd(11,a,b) + w*af*Jac*(T1 + T2 + T3)

!           dC/dU_3
            T1    = N(a)*(rC*NxFi(3,b) - VxNx(3,b))
            T2    = tauM*(rMNx(a)*NxFi(3,b) - NxFi(3,a)*rMNx(b))
            T3    = -tauM*NxNx*PxFi(3)/rho
            lKd(12,a,b) = lKd(12,a,b) + w*af*Jac*(T1 + T2 + T3)

!           dM_1/dV_1
            T1 = am*rho*N(a)*N(b)
            T2 = af*tauC*rho*NxFi(1,a)*NxFi(1,b) + T1
            lK(1,a,b)   = lK(1,a,b)  + w*Jac*T2

!           dM_1/dV_2
            T2 = af*tauC*rho*NxFi(1,a)*NxFi(2,b)
            lK(2,a,b)   = lK(2,a,b)  + w*Jac*T2

!           dM_1/dV_3
            T2 = af*tauC*rho*NxFi(1,a)*NxFi(3,b)
            lK(3,a,b)   = lK(3,a,b)  + w*Jac*T2

!           dM_2/dV_1
            T2 = af*tauC*rho*NxFi(2,a)*NxFi(1,b)
            lK(5,a,b)   = lK(5,a,b)  + w*Jac*T2

!           dM_2/dV_2
            T2 = af*tauC*rho*NxFi(2,a)*NxFi(2,b) + T1
            lK(6,a,b)   = lK(6,a,b)  + w*Jac*T2

!           dM_2/dV_7
            T2 = af*tauC*rho*NxFi(2,a)*NxFi(3,b)
            lK(7,a,b)   = lK(7,a,b)  + w*Jac*T2

!           dM_3/dV_1
            T2 = af*tauC*rho*NxFi(3,a)*NxFi(1,b)
            lK(9,a,b)   = lK(9,a,b)  + w*Jac*T2

!           dM_3/dV_2
            T2 = af*tauC*rho*NxFi(3,a)*NxFi(2,b)
            lK(10,a,b)  = lK(10,a,b) + w*Jac*T2

!           dM_3/dV_3
            T2 = af*tauC*rho*NxFi(3,a)*NxFi(3,b) + T1
            lK(11,a,b)  = lK(11,a,b) + w*Jac*T2

!           dM_1/dP
            T1 = am*tauC*rho*beta + af*(tauC*pd*drb - 1D0)
            T2 = NxFi(1,a)*N(b)*T1
            lK(4,a,b)   = lK(4,a,b)  + w*Jac*T2

!           dM_2/dP
            T2 = NxFi(2,a)*N(b)*T1
            lK(8,a,b)   = lK(8,a,b)  + w*Jac*T2

!           dM_3/dP
            T2 = NxFi(3,a)*N(b)*T1
            lK(12,a,b)  = lK(12,a,b) + w*Jac*T2

!           dC/dV_1
            T2 = (am*tauM/rho)*NxFi(1,a)*N(b) + af*N(a)*NxFi(1,b)
            lK(13,a,b)  = lK(13,a,b) + w*Jac*T2

!           dC/dV_2
            T2 = (am*tauM/rho)*NxFi(2,a)*N(b) + af*N(a)*NxFi(2,b)
            lK(14,a,b)  = lK(14,a,b) + w*Jac*T2

!           dC/dV_3
            T2 = (am*tauM/rho)*NxFi(3,a)*N(b) + af*N(a)*NxFi(3,b)
            lK(15,a,b)  = lK(15,a,b) + w*Jac*T2

!           dC/dP
            T1 = (am*beta + af*dbeta*pd)*N(a)*N(b)
            T2 = T1 + af*tauM/rho *(NxNx - NxPx*N(b)*drho/rho)
            lK(16,a,b)  = lK(16,a,b) + w*Jac*T2
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE USTRUCT2D(eNoN, w, Je, N, Nx, al, yl, dl, adl, fNl,
     2   lRd, lKd, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, Je, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), adl(2,eNoN),
     3   fNl(nFn*2,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lRd(2,eNoN), lR(dof,eNoN),
     2   lKd(2*dof,eNoN,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: a, b, i, j, k, iFn
      REAL(KIND=8) :: rho, drho, bf(2), am, af, v(2), vd(2), vx(2,2), p,
     2   pd, px(2), ud(2), fl(2,nFn), F(2,2), Jac, Fi(2,2), beta, dbeta,
     3   Siso(2,2), CCiso(2,2,2,2), tauM, tauC, rC, rCl, rM(2), Dm(3,3),
     4   Pdev(2,2), Bm(3,2,eNoN), DBm(3,2), NxFi(2,eNoN), VxFi(2,2), T1,
     5   T2, T3, PxFi(2), VxNx(2,eNoN), NxNx, NxSNx, NxPx, BtDB, drb,
     6   rMNx(eNoN)

      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      bf(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt

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
      ud = 0D0
      fl = 0D0
      F  = MAT_ID(2)
      DO a=1, eNoN
         v(1)    = v(1)  + N(a)*yl(i,a)
         v(2)    = v(2)  + N(a)*yl(j,a)

         vd(1)   = vd(1) + N(a)*al(i,a)
         vd(2)   = vd(2) + N(a)*al(j,a)

         vx(1,1) = vx(1,1) + Nx(1,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nx(2,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)

         p       = p     + N(a)*yl(k,a)
         pd      = pd    + N(a)*al(k,a)
         px(1)   = px(1) + Nx(1,a)*yl(k,a)
         px(2)   = px(2) + Nx(2,a)*yl(k,a)

         ud(1)   = ud(1) + N(a)*adl(1,a)
         ud(2)   = ud(2) + N(a)*adl(2,a)

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
      drb = drho*beta + rho*dbeta

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     elasticity tensor (CCiso)
      CALL GETPK2CCdev(stModel, F, fl, Siso, CCiso)

!     Compute stabilization parameters
      CALL GETTAU(stModel, Je, F, fl, tauM, tauC)

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

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(2,1)*Fi(2,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(2,1)*Fi(2,2)

      VxFi(2,1) = vx(1,2)*Fi(1,1) + vx(2,2)*Fi(2,1)
      VxFi(2,2) = vx(1,2)*Fi(1,2) + vx(2,2)*Fi(2,2)

      PxFi(1)   = px(1)*Fi(1,1) + px(2)*Fi(2,1)
      PxFi(2)   = px(1)*Fi(1,2) + px(2)*Fi(2,2)

      DO a=1, eNoN
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a)
      END DO

      rC    = beta*pd + VxFi(1,1) + VxFi(2,2)
      rCl   = -p + tauC*rho*rC
      rM(1) = vd(1) + PxFi(1)/rho
      rM(2) = vd(2) + PxFi(2)/rho

!     Local residue
      DO a=1, eNoN
         lRd(1,a) = lRd(1,a) + w*N(a)*Jac*(ud(1) - v(1))
         lRd(2,a) = lRd(2,a) + w*N(a)*Jac*(ud(2) - v(2))

         T1 = N(a)*rho*vd(1)*Jac
         T2 = Nx(1,a)*Pdev(1,1) + Nx(2,a)*Pdev(1,2)
         T3 = NxFi(1,a)*rCl*Jac
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = N(a)*rho*vd(2)*Jac
         T2 = Nx(1,a)*Pdev(2,1) + Nx(2,a)*Pdev(2,2)
         T3 = NxFi(2,a)*rCl*Jac
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)

         rMNx(a) = rM(1)*NxFi(1,a) + rM(2)*NxFi(2,a)
         lR(3,a) = lR(3,a) + w*Jac*(N(a)*rC + tauM*rMNx(a))
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx  = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b)
            NxPx  = NxFi(1,a)*PxFi(1) + NxFi(2,a)*PxFi(2)
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

!           dM_1/dU_1
            BtDB  = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1)
            T1    = N(a)*rho*vd(1)*NxFi(1,b)*Jac
            T2    = -NxFi(1,a)*tauC*rho*VxNx(1,b)*Jac
            lKd(1,a,b) = lKd(1,a,b) + w*af*(BtDB + T1 + T2 + NxSNx)

!           dM_1/dU_2
            BtDB  = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2              Bm(3,1,a)*DBm(3,2)
            T1    = N(a)*rho*vd(1)*NxFi(2,b)*Jac
            T2    = -NxFi(1,a)*tauC*rho*VxNx(2,b)*Jac
            T3    = (NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))*rCl*Jac
            lKd(2,a,b) = lKd(2,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_2/dU_1
            BtDB  = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2              Bm(3,2,a)*DBm(3,1)
            T1    = N(a)*rho*vd(2)*NxFi(1,b)*Jac
            T2    = -NxFi(2,a)*tauC*rho*VxNx(1,b)*Jac
            T3    = (NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))*rCl*Jac
            lKd(3,a,b) = lKd(3,a,b) + w*af*(BtDB + T1 + T2 + T3)

!           dM_2/dU_2
            BtDB  = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2              Bm(3,2,a)*DBm(3,2)
            T1    = N(a)*rho*vd(2)*NxFi(2,b)*Jac
            T2    = -NxFi(2,a)*tauC*rho*VxNx(2,b)*Jac
            lKd(4,a,b) = lKd(4,a,b) + w*af*(BtDB + T1 + T2 + NxSNx)

!           dC/dU_1
            T1    = N(a)*(rC*NxFi(1,b) - VxNx(1,b))
            T2    = tauM*(rMNx(a)*NxFi(1,b) - NxFi(1,a)*rMNx(b))
            T3    = -tauM*NxNx*PxFi(1)/rho
            lKd(5,a,b) = lKd(5,a,b) + w*af*Jac*(T1 + T2 + T3)

!           dC/dU_2
            T1    = N(a)*(rC*NxFi(2,b) - VxNx(2,b))
            T2    = tauM*(rMNx(a)*NxFi(2,b) - NxFi(2,a)*rMNx(b))
            T3    = -tauM*NxNx*PxFi(2)/rho
            lKd(6,a,b) = lKd(6,a,b) + w*af*Jac*(T1 + T2 + T3)

!           dM_1/dV_1
            T1 = am*rho*N(a)*N(b)
            T2 = af*tauC*rho*NxFi(1,a)*NxFi(1,b) + T1
            lK(1,a,b)  = lK(1,a,b) + w*Jac*T2

!           dM_1/dV_2
            T2 = af*tauC*rho*NxFi(1,a)*NxFi(2,b)
            lK(2,a,b)  = lK(2,a,b) + w*Jac*T2

!           dM_2/dV_1
            T2 = af*tauC*rho*NxFi(2,a)*NxFi(1,b)
            lK(4,a,b)  = lK(4,a,b) + w*Jac*T2

!           dM_2/dV_2
            T2 = af*tauC*rho*NxFi(2,a)*NxFi(2,b) + T1
            lK(5,a,b)  = lK(5,a,b) + w*Jac*T2

!           dM_1/dP
            T1 = am*tauC*rho*beta + af*(tauC*pd*drb - 1D0)
            T2 = NxFi(1,a)*N(b)*T1
            lK(3,a,b)  = lK(3,a,b) + w*Jac*T2

!           dM_2/dP
            T2 = NxFi(2,a)*N(b)*T1
            lK(6,a,b)  = lK(6,a,b) + w*Jac*T2

!           dC/dV_1
            T2 = (am*tauM/rho)*NxFi(1,a)*N(b) + af*N(a)*NxFi(1,b)
            lK(7,a,b)  = lK(7,a,b) + w*Jac*T2

!           dC/dV_2
            T2 = (am*tauM/rho)*NxFi(2,a)*N(b) + af*N(a)*NxFi(2,b)
            lK(8,a,b)  = lK(8,a,b) + w*Jac*T2

!           dC/dP
            T1 = (am*beta + af*dbeta*pd)*N(a)*N(b)
            T2 = T1 + af*tauM/rho *(NxNx - NxPx*N(b)*drho/rho)
            lK(9,a,b)  = lK(9,a,b) + w*Jac*T2
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
      LOGICAL :: l1, l2

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:), Nx(:,:),
     2   dl(:,:), hl(:), lRd(:,:), lR(:,:), lKd(:,:,:), lK(:,:,:)

      iM    = lFa%iM
      eNoNb = lFa%eNoN
      eNoN  = msh(iM)%eNoN

      ALLOCATE(xl(nsd,eNoN), N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN),
     2   ptr(eNoN), dl(tDof,eNoN), hl(eNoN), lRd(nsd,eNoN),lR(dof,eNoN),
     3   lKd(nsd*dof,eNoN,eNoN), lK(dof*dof,eNoN,eNoN))

      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_ustruct) CYCLE

         lRd   = 0D0
         lR    = 0D0
         lKd   = 0D0
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
            CALL GETXI(msh(iM)%eType, eNoN, xl, xp, xi, l1)
            CALL GETGNN(nsd, msh(iM)%eType, eNoN, xi, N, Nxi)

            a  = 0
            rt = 0D0
            DO b=1, eNoN
               rt = rt + N(b)
               IF (N(b).GT.-1E-4 .AND. N(b).LT.1.0001D0) a = a + 1
            END DO
            l2 = rt.GE.0.9999D0 .AND. rt.LE.1.0001D0
            IF (.NOT.l1 .OR. .NOT.l2 .OR. a.NE.eNoN) err =
     2         " Error in computing face derivatives (BUSTRUCTNEU)"

            IF (g.EQ.1 .OR. .NOT.msh(iM)%lShpF)
     2         CALL GNN(eNoN, Nxi, xl, Nx, rt, ksix)

            CALL GNNB(lFa, e, g, nV)
            Jac = SQRT(NORM(nV))
            nV  = nV / Jac
            w   = lFa%w(g)*Jac

            IF (nsd .EQ. 3) THEN
               CALL BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lKd)
            ELSE
               CALL BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lKd)
            END IF
         END DO
!$OMP CRITICAL
         CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lRd)
         CALL DOASSEM(eNoN, ptr, lK, lR)
!$OMP END CRITICAL
      END DO

      RETURN
      END SUBROUTINE BUSTRUCTNEU
!--------------------------------------------------------------------
      SUBROUTINE BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(3)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lKd(3*dof,eNoN,eNoN)

      INTEGER :: a, b, i, j, k
      REAL(KIND=8) :: af, Jac, h, F(3,3), Fi(3,3), nFi(3), NxFi(3,eNoN),
     2   wl, T1, T2

      af = eq(cEq)%af*eq(cEq)%gam*dt
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

      wl = w*h*Jac
      DO a=1, eNoN
         T1 = wl*N(a)
         lR(1,a) = lR(1,a) - T1*nFi(1)
         lR(2,a) = lR(2,a) - T1*nFi(2)
         lR(3,a) = lR(3,a) - T1*nFi(3)

         T2 = T1*af
         DO b=1, eNoN
            T1 = T2*(nFi(1)*NxFi(2,b) - nFi(2)*NxFi(1,b))
            lKd(2,a,b) = lKd(2,a,b) - T1

            T1 = T2*(nFi(1)*NxFi(3,b) - nFi(3)*NxFi(1,b))
            lKd(3,a,b) = lKd(3,a,b) - T1

            T1 = T2*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lKd(4,a,b) = lKd(4,a,b) - T1

            T1 = T2*(nFi(2)*NxFi(3,b) - nFi(3)*NxFi(2,b))
            lKd(6,a,b) = lKd(6,a,b) - T1

            T1 = T2*(nFi(3)*NxFi(1,b) - nFi(1)*NxFi(3,b))
            lKd(7,a,b) = lKd(7,a,b) - T1

            T1 = T2*(nFi(3)*NxFi(2,b) - nFi(2)*NxFi(3,b))
            lKd(8,a,b) = lKd(8,a,b) - T1
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(2)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lKd(2*dof,eNoN,eNoN)

      INTEGER :: a, b, i, j
      REAL(KIND=8) :: af, Jac, h, F(2,2), Fi(2,2), nFi(2), NxFi(2,eNoN),
     2   wl, T1, T2

      af = eq(cEq)%af*eq(cEq)%gam*dt
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

      wl = w*h*Jac
      DO a=1, eNoN
         T1 = wl*N(a)
         lR(1,a) = lR(1,a) - T1*nFi(1)
         lR(2,a) = lR(2,a) - T1*nFi(2)

         T2 = T1*af
         DO b=1, eNoN
            T1 = T2*(nFi(1)*NxFi(2,b) - nFi(2)*NxFi(1,b))
            lKd(2,a,b) = lKd(2,a,b) - T1

            T1 = T2*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lKd(4,a,b) = lKd(4,a,b) - T1
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT2D
!####################################################################
      SUBROUTINE USTRUCT_DOASSEM (d, eqN, lKd, lRd)
      USE COMMOD, ONLY: dof, nsd, rowPtr, colPtr, Rd, Kd
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: d, eqN(d)
      REAL(KIND=8), INTENT(IN) :: lKd(nsd*dof,d,d), lRd(nsd,d)

      INTEGER a, b, ptr, rowN, colN, left, right

      DO a=1, d
         rowN = eqN(a)
         Rd(:,rowN) = Rd(:,rowN) + lRd(:,a)
         DO b=1, d
            colN  = eqN(b)
            left  = rowPtr(rowN)
            right = rowPtr(rowN+1)
            ptr   = (right + left)/2
            DO WHILE (colN .NE. colPtr(ptr))
               IF (colN .GT. colPtr(ptr)) THEN
                  left  = ptr
               ELSE
                  right = ptr
               END IF
               ptr = (right + left)/2
            END DO
            Kd(:,ptr) = Kd(:,ptr) + lKd(:,a,b)
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT_DOASSEM
!--------------------------------------------------------------------
      SUBROUTINE USTRUCT_updateValR
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER :: i, j, c
      REAL(KIND=8) :: ami, afl, KdRd(dof)

      ami = 1D0/eq(cEq)%am
      afl = eq(cEq)%af*eq(cEq)%gam*dt * ami

      IF (nsd .EQ. 3) THEN
!        Update stiffness matrices (dM/dU, dC/dU only)
         DO i=1, lhs%nnz
            Val(1,i)  = Val(1,i)  + afl*Kd(1,i)
            Val(2,i)  = Val(2,i)  + afl*Kd(2,i)
            Val(3,i)  = Val(3,i)  + afl*Kd(3,i)

            Val(5,i)  = Val(5,i)  + afl*Kd(4,i)
            Val(6,i)  = Val(6,i)  + afl*Kd(5,i)
            Val(7,i)  = Val(7,i)  + afl*Kd(6,i)

            Val(9,i)  = Val(9,i)  + afl*Kd(7,i)
            Val(10,i) = Val(10,i) + afl*Kd(8,i)
            Val(11,i) = Val(11,i) + afl*Kd(9,i)

            Val(13,i) = Val(13,i) + afl*Kd(10,i)
            Val(14,i) = Val(14,i) + afl*Kd(11,i)
            Val(15,i) = Val(15,i) + afl*Kd(12,i)
         END DO

!        Update RHS
         DO i=1, tnNo
            KdRd = 0D0
            DO j=rowPtr(i), rowPtr(i+1)-1
               c = colPtr(j)
               KdRd(1) = KdRd(1) + Kd(1 ,j)*Rd(1,c) + Kd(2 ,j)*Rd(2,c)
     2                           + Kd(3 ,j)*Rd(3,c)
               KdRd(2) = KdRd(2) + Kd(4 ,j)*Rd(1,c) + Kd(5 ,j)*Rd(2,c)
     2                           + Kd(6 ,j)*Rd(3,c)
               KdRd(3) = KdRd(3) + Kd(7 ,j)*Rd(1,c) + Kd(8 ,j)*Rd(2,c)
     2                           + Kd(9 ,j)*Rd(3,c)
               KdRd(4) = KdRd(4) + Kd(10,j)*Rd(1,c) + Kd(11,j)*Rd(2,c)
     2                           + Kd(12,j)*Rd(3,c)
            END DO
            R(1,i) = R(1,i) - ami*KdRd(1)
            R(2,i) = R(2,i) - ami*KdRd(2)
            R(3,i) = R(3,i) - ami*KdRd(3)
            R(4,i) = R(4,i) - ami*KdRd(4)
         END DO
      ELSE
!        Update stiffness matrices (dM/dU, dC/dU only)
         DO i=1, lhs%nnz
            Val(1,i)  = Val(1,i)  + afl*Kd(1,i)
            Val(2,i)  = Val(2,i)  + afl*Kd(2,i)

            Val(4,i)  = Val(4,i)  + afl*Kd(3,i)
            Val(5,i)  = Val(5,i)  + afl*Kd(4,i)

            Val(7,i)  = Val(7,i)  + afl*Kd(5,i)
            Val(8,i)  = Val(8,i)  + afl*Kd(6,i)
         END DO

!        Update RHS
         DO i=1, tnNo
            KdRd = 0D0
            DO j=rowPtr(i), rowPtr(i+1)-1
               c = colPtr(j)
               KdRd(1) = KdRd(1) + Kd(1 ,j)*Rd(1,c) + Kd(2 ,j)*Rd(2,c)
               KdRd(2) = KdRd(2) + Kd(3 ,j)*Rd(1,c) + Kd(4 ,j)*Rd(2,c)
               KdRd(3) = KdRd(3) + Kd(5 ,j)*Rd(1,c) + Kd(6 ,j)*Rd(2,c)
            END DO
            R(1,i) = R(1,i) - ami*KdRd(1)
            R(2,i) = R(2,i) - ami*KdRd(2)
            R(3,i) = R(3,i) - ami*KdRd(3)
         END DO
      END IF

      RETURN
      END SUBROUTINE USTRUCT_updateValR
!####################################################################

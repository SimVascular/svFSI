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
!     These routines are for solving the nonlinear structural dynamics
!     problem (velocity-pressure based formulation).
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_uSOLID(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      LOGICAL vmsStab
      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys, iFn, nFn
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), fN(:,:), ya_l(:), lR(:,:), lK(:,:,:),
     3   lKd(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nqx(:,:)

      eNoN = lM%eNoN
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
      END IF

!     USTRUCT: dof = nsd+1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(nsd,eNoN), fN(nsd,nFn), ya_l(eNoN),
     3   lR(dof,eNoN), lK(dof*dof,eNoN,eNoN), lKd(dof*nsd,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_ustruct) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         fN   = 0._RKIND
         ya_l = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
            IF (ALLOCATED(lM%fN)) THEN
               DO iFn=1, nFn
                  fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
               END DO
            END IF
            IF (cem%cpld) ya_l(a) = cem%Ya(Ac)
         END DO

!        Initialize residue and tangents
         lR  = 0._RKIND
         lK  = 0._RKIND
         lKd = 0._RKIND

!        Set function spaces for velocity/displacement and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND
         Nqx      = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(1)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL USTRUCT3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn, w,
     2            Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl, dl, bfl,
     3            fN, ya_l, lR, lK, lKd)

            ELSE IF (nsd .EQ. 2) THEN
               CALL USTRUCT2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, nFn, w,
     2            Jac, fs(1)%N(:,g), fs(2)%N(:,g), Nwx, al, yl, dl, bfl,
     3            fN, ya_l, lR, lK, lKd)

            END IF
         END DO ! g: loop

!        Set function spaces for velocity/displacement and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 2)

!        Gauss integration 2
         DO g=1, fs(2)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(2)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL USTRUCT3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, Jac,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, dl, bfl,
     3            lR, lK, lKd)

            ELSE IF (nsd .EQ. 2) THEN
               CALL USTRUCT2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, Jac,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, dl, bfl,
     3            lR, lK, lKd)

            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nqx)

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) err = "Cannot assemble USTRUCT using "//
     2      "Trilinos"
#endif
         CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lK, lR)
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

      RETURN
      END SUBROUTINE CONSTRUCT_uSOLID
!####################################################################
      SUBROUTINE USTRUCT3D_M(vmsFlag, eNoNw, eNoNq, nFn, w, Je, Nw, Nq,
     2   Nwx, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, Je, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), al(tDof,eNoNw), yl(tDof,eNoNw), dl(tDof,eNoNw),
     3   bfl(3,eNoNw), fN(3,nFn), ya_l(eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lKd(dof*3,eNoNw,eNoNw), lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, l, a, b
      REAL(KIND=RKIND) :: fb(3), am, af, afm, v(3), vd(3), vx(3,3), p,
     2   pd, F(3,3), Jac, Fi(3,3), mu, rho, beta, drho, dbeta, ya_g, Ja,
     3   Siso(3,3), Dm(6,6), tauM, tauC, rC, rCl, Pdev(3,3), DBm(6,3),
     4   Bm(6,3,eNoNw), NxFi(3,eNoNw), VxFi(3,3), VxNx(3,eNoNw), BtDB,
     5   NxSNx, NxNx, Ku, T1, T2, T3, T4, r13, r23, ddev(3,3),
     6   Pvis(3,3), PvNx(3,eNoNw)

!     Define parameters
      mu      = eq(cEq)%dmn(cDmn)%prop(solid_viscosity)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      afm     = af / am

!     {i,j,k} := velocity dofs; {l} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1
      l       = k + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0._RKIND
      vd = -fb
      vx = 0._RKIND
      F  = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      ya_g   = 0._RKIND
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)
         v(3)    = v(3)  + Nw(a)*yl(k,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))
         vd(3)   = vd(3) + Nw(a)*(al(k,a)-bfl(3,a))

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(1,3) = vx(1,3) + Nwx(3,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)
         vx(2,3) = vx(2,3) + Nwx(3,a)*yl(j,a)
         vx(3,1) = vx(3,1) + Nwx(1,a)*yl(k,a)
         vx(3,2) = vx(3,2) + Nwx(2,a)*yl(k,a)
         vx(3,3) = vx(3,3) + Nwx(3,a)*yl(k,a)

         F(1,1)  = F(1,1) + Nwx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nwx(2,a)*dl(i,a)
         F(1,3)  = F(1,3) + Nwx(3,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nwx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nwx(2,a)*dl(j,a)
         F(2,3)  = F(2,3) + Nwx(3,a)*dl(j,a)
         F(3,1)  = F(3,1) + Nwx(1,a)*dl(k,a)
         F(3,2)  = F(3,2) + Nwx(2,a)*dl(k,a)
         F(3,3)  = F(3,3) + Nwx(3,a)*dl(k,a)

         ya_g    = ya_g + Nw(a)*ya_l(a)
      END DO
      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

!     Pressure and its time derivative
      p  = 0._RKIND
      pd = 0._RKIND
      DO a=1, eNoNq
         p  = p  + Nq(a)*yl(l,a)
         pd = pd + Nq(a)*al(l,a)
      END DO

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     isochoric elasticity tensor in Voigt notation (Dm)
      CALL GETPK2CCdev(eq(cEq)%dmn(cDmn), F, nFn, fN, ya_g, Siso, Dm,Ja)

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(eq(cEq)%dmn(cDmn), p, rho, beta, drho, dbeta, Ja)

!     Compute stabilization parameters
      IF (vmsFlag) THEN
         CALL GETTAU(eq(cEq)%dmn(cDmn), Jac, Je, tauM, tauC)
      ELSE
         tauM = 0._RKIND
         tauC = 0._RKIND
      END IF

!     Deviatoric 1st Piola-Kirchhoff tensor (P)
      Pdev = MATMUL(F, Siso)

!     Viscous contribution
      ddev = 2._RKIND*mu*Jac*MAT_DEV(MAT_SYMM(vx,3), 3)
      Pvis = MATMUL(ddev, TRANSPOSE(Fi))

!     Total 1st Piola-Kirchhoff stress
      Pdev = Pdev + Pvis

      DO a=1, eNoNw
         Bm(1,1,a) = Nwx(1,a)*F(1,1)
         Bm(1,2,a) = Nwx(1,a)*F(2,1)
         Bm(1,3,a) = Nwx(1,a)*F(3,1)

         Bm(2,1,a) = Nwx(2,a)*F(1,2)
         Bm(2,2,a) = Nwx(2,a)*F(2,2)
         Bm(2,3,a) = Nwx(2,a)*F(3,2)

         Bm(3,1,a) = Nwx(3,a)*F(1,3)
         Bm(3,2,a) = Nwx(3,a)*F(2,3)
         Bm(3,3,a) = Nwx(3,a)*F(3,3)

         Bm(4,1,a) = (Nwx(1,a)*F(1,2) + F(1,1)*Nwx(2,a))
         Bm(4,2,a) = (Nwx(1,a)*F(2,2) + F(2,1)*Nwx(2,a))
         Bm(4,3,a) = (Nwx(1,a)*F(3,2) + F(3,1)*Nwx(2,a))

         Bm(5,1,a) = (Nwx(2,a)*F(1,3) + F(1,2)*Nwx(3,a))
         Bm(5,2,a) = (Nwx(2,a)*F(2,3) + F(2,2)*Nwx(3,a))
         Bm(5,3,a) = (Nwx(2,a)*F(3,3) + F(3,2)*Nwx(3,a))

         Bm(6,1,a) = (Nwx(3,a)*F(1,1) + F(1,3)*Nwx(1,a))
         Bm(6,2,a) = (Nwx(3,a)*F(2,1) + F(2,3)*Nwx(1,a))
         Bm(6,3,a) = (Nwx(3,a)*F(3,1) + F(3,3)*Nwx(1,a))
      END DO

      DO a=1, eNoNw
         NxFi(1,a) = Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1) +
     2      Nwx(3,a)*Fi(3,1)
         NxFi(2,a) = Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2) +
     2      Nwx(3,a)*Fi(3,2)
         NxFi(3,a) = Nwx(1,a)*Fi(1,3) + Nwx(2,a)*Fi(2,3) +
     2      Nwx(3,a)*Fi(3,3)

         PvNx(1,a) = ddev(1,1)*NxFi(1,a) + ddev(1,2)*NxFi(2,a) +
     2      ddev(1,3)*NxFi(3,a)
         PvNx(2,a) = ddev(2,1)*NxFi(1,a) + ddev(2,2)*NxFi(2,a) +
     2      ddev(2,3)*NxFi(3,a)
         PvNx(3,a) = ddev(3,1)*NxFi(1,a) + ddev(3,2)*NxFi(2,a) +
     2      ddev(3,3)*NxFi(3,a)
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

      rC  = beta*pd + VxFi(1,1) + VxFi(2,2) + VxFi(3,3)
      rCl = -p + tauC*rC

!     Local residue
      DO a=1, eNoNw
         T1 = Jac*rho*vd(1)*Nw(a)
         T2 = Pdev(1,1)*Nwx(1,a) + Pdev(1,2)*Nwx(2,a) +
     2      Pdev(1,3)*Nwx(3,a)
         T3 = Jac*rCl*NxFi(1,a)
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(2)*Nw(a)
         T2 = Pdev(2,1)*Nwx(1,a) + Pdev(2,2)*Nwx(2,a) +
     2      Pdev(2,3)*Nwx(3,a)
         T3 = Jac*rCl*NxFi(2,a)
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(3)*Nw(a)
         T2 = Pdev(3,1)*Nwx(1,a) + Pdev(3,2)*Nwx(2,a) +
     2      Pdev(3,3)*Nwx(3,a)
         T3 = Jac*rCl*NxFi(3,a)
         lR(3,a) = lR(3,a) + w*(T1 + T2 + T3)
      END DO

      DO a=1, eNoNw
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a) +
     2               VxFi(3,1)*NxFi(3,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a) +
     2               VxFi(3,2)*NxFi(3,a)
         VxNx(3,a) = VxFi(1,3)*NxFi(1,a) + VxFi(2,3)*NxFi(2,a) +
     2               VxFi(3,3)*NxFi(3,a)
      END DO

!     Tangent (stiffness) matrices
      r13 = 1._RKIND / 3._RKIND
      r23 = 2._RKIND / 3._RKIND
      DO b=1, eNoNw
         DO a=1, eNoNw
            NxSNx = Nwx(1,a)*Siso(1,1)*Nwx(1,b)
     2       + Nwx(1,a)*Siso(1,2)*Nwx(2,b) + Nwx(1,a)*Siso(1,3)*Nwx(3,b)
     3       + Nwx(2,a)*Siso(2,1)*Nwx(1,b) + Nwx(2,a)*Siso(2,2)*Nwx(2,b)
     4       + Nwx(2,a)*Siso(2,3)*Nwx(3,b) + Nwx(3,a)*Siso(3,1)*Nwx(1,b)
     5       + Nwx(3,a)*Siso(3,2)*Nwx(2,b) + Nwx(3,a)*Siso(3,3)*Nwx(3,b)
            DBm   = MATMUL(Dm, Bm(:,:,b))

            NxNx = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b)
     2           + NxFi(3,a)*NxFi(3,b)

!           dM1_dV1 + af/am *dM_1/dU_1
            BtDB = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2             Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     3             Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            T1   = Jac*rho*vd(1)*Nw(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b)
            T3   = PvNx(1,a)*NxFi(1,b) - NxFi(1,a)*PvNx(1,b)
            Ku   = w*af*(T1 + T2 + T3 + BtDB + NxSNx)
            lKd(1,a,b) = lKd(1,a,b) + Ku

            T1   = am*Jac*rho*Nw(a)*Nw(b)
            T2   = T1 + af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b)
            T3   = T2 + af*mu*Jac*(r13*NxFi(1,a)*NxFi(1,b) + NxNx)
            lK(1,a,b)  = lK(1,a,b) + w*T3 + afm*Ku

!           dM_1/dV_2 + af/am *dM_1/dU_2
            BtDB = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2             Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     3             Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            T1   = Jac*rho*vd(1)*Nw(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))
            T4   = PvNx(1,a)*NxFi(2,b) - NxFi(2,a)*PvNx(1,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(2,a,b) = lKd(2,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(2,b)
            T3   = T2 + af*mu*Jac*(NxFi(2,a)*NxFi(1,b)
     2           - r23*NxFi(1,a)*NxFi(2,b))
            lK(2,a,b) = lK(2,a,b) + w*T3 + afm*Ku

!           dM_1/dV_3 + af/am *dM_1/dU_3
            BtDB = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2             Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     3             Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            T1   = Jac*rho*vd(1)*Nw(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(3,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(3,b) - NxFi(3,a)*NxFi(1,b))
            T4   = PvNx(1,a)*NxFi(3,b) - NxFi(3,a)*PvNx(1,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(3,a,b) = lKd(3,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(3,b)
            T3   = T2 + af*mu*Jac*(NxFi(3,a)*NxFi(1,b)
     2           - r23*NxFi(1,a)*NxFi(3,b))
            lK(3,a,b) = lK(3,a,b) + w*T3 + afm*Ku

!           dM_2/dV_1 + af/am *dM_2/dU_1
            BtDB = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2             Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     3             Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            T1   = Jac*rho*vd(2)*Nw(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))
            T4   = PvNx(2,a)*NxFi(1,b) - NxFi(1,a)*PvNx(2,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(4,a,b) = lKd(4,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(1,b)
            T3   = T2 + af*mu*Jac*(NxFi(1,a)*NxFi(2,b)
     2           - r23*NxFi(2,a)*NxFi(1,b))
            lK(5,a,b) = lK(5,a,b) + w*T3 + afm*Ku

!           dM_2/dV_2 + af/am *dM_2/dU_2
            BtDB = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2             Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     3             Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            T1   = Jac*rho*vd(2)*Nw(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(2,b)
            T3   = PvNx(2,a)*NxFi(2,b) - NxFi(2,a)*PvNx(2,b)
            Ku   = w*af*(T1 + T2 + T3 + BtDB + NxSNx)
            lKd(5,a,b) = lKd(5,a,b) + Ku

            T1   = am*Jac*rho*Nw(a)*Nw(b)
            T2   = T1 + af*Jac*tauC*rho*NxFi(2,a)*NxFi(2,b)
            T3   = T2 + af*mu*Jac*(r13*NxFi(2,a)*NxFi(2,b) + NxNx)
            lK(6,a,b) = lK(6,a,b) + w*T3 + afm*Ku

!           dM_2/dV_3 + af/am *dM_2/dU_3
            BtDB = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2             Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     3             Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            T1   = Jac*rho*vd(2)*Nw(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(3,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(3,b) - NxFi(3,a)*NxFi(2,b))
            T4   = PvNx(2,a)*NxFi(3,b) - NxFi(3,a)*PvNx(2,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(6,a,b) = lKd(6,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(3,b)
            T3   = T2 + af*mu*Jac*(NxFi(3,a)*NxFi(2,b)
     2           - r23*NxFi(2,a)*NxFi(3,b))
            lK(7,a,b) = lK(7,a,b) + w*T3 + afm*Ku

!           dM_3/dV_1 + af/am *dM_3/dU_1
            BtDB = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2             Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     3             Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            T1   = Jac*rho*vd(3)*Nw(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(3,a)*NxFi(1,b) - NxFi(1,a)*NxFi(3,b))
            T4   = PvNx(3,a)*NxFi(1,b) - NxFi(1,a)*PvNx(3,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(7,a,b) = lKd(7,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(3,a)*NxFi(1,b)
            T3   = T2 + af*mu*Jac*(NxFi(1,a)*NxFi(3,b)
     2           - r23*NxFi(3,a)*NxFi(1,b))
            lK(9,a,b) = lK(9,a,b) + w*T3 + afm*Ku

!           dM_3/dV_2 + af/am *dM_3/dU_2
            BtDB = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2             Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     3             Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            T1   = Jac*rho*vd(3)*Nw(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(3,a)*NxFi(2,b) - NxFi(2,a)*NxFi(3,b))
            T4   = PvNx(3,a)*NxFi(2,b) - NxFi(2,a)*PvNx(3,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(8,a,b) = lKd(8,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(3,a)*NxFi(2,b)
            T3   = T2 + af*mu*Jac*(NxFi(2,a)*NxFi(3,b)
     2           - r23*NxFi(3,a)*NxFi(2,b))
            lK(10,a,b) = lK(10,a,b) + w*T3 + afm*Ku

!           dM_3/dV_3 + af/am *dM_3/dU_3
            BtDB = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2             Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     3             Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            T1   = Jac*rho*vd(3)*Nw(a)*NxFi(3,b)
            T2   = -tauC*Jac*NxFi(3,a)*VxNx(3,b)
            T3   = PvNx(3,a)*NxFi(3,b) - NxFi(3,a)*PvNx(3,b)
            Ku   = w*af*(T1 + T2 + T3 + BtDB + NxSNx)
            lKd(9,a,b) = lKd(9,a,b) + Ku

            T1   = am*Jac*rho*Nw(a)*Nw(b)
            T2   = T1 + af*Jac*tauC*rho*NxFi(3,a)*NxFi(3,b)
            T3   = T2 + af*mu*Jac*(r13*NxFi(3,a)*NxFi(3,b) + NxNx)
            lK(11,a,b) = lK(11,a,b) + w*T3 + afm*Ku
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
!           dM_1/dP
            T1 = am*tauC*beta + af*(tauC*dbeta*pd - 1._RKIND)
            T2 = T1*NxFi(1,a)*Nq(b) + af*drho*vd(1)*Nw(a)*Nq(b)
            lK(4,a,b) = lK(4,a,b) + w*Jac*T2

!           dM_2/dP
            T2 = T1*NxFi(2,a)*Nq(b) + af*drho*vd(2)*Nw(a)*Nq(b)
            lK(8,a,b) = lK(8,a,b) + w*Jac*T2

!           dM_3/dP
            T2 = T1*NxFi(3,a)*Nq(b) + af*drho*vd(3)*Nw(a)*Nq(b)
            lK(12,a,b) = lK(12,a,b) + w*Jac*T2
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT3D_M
!--------------------------------------------------------------------
      SUBROUTINE USTRUCT2D_M(vmsFlag, eNoNw, eNoNq, nFn, w, Je, Nw, Nq,
     2   Nwx, al, yl, dl, bfl, fN, ya_l, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, Je, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), al(tDof,eNoNw), yl(tDof,eNoNw), dl(tDof,eNoNw),
     3   bfl(2,eNoNw), fN(2,nFn), ya_l(eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lKd(dof*2,eNoNw,eNoNw), lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: fb(2), am, af, afm, v(2), vd(2), vx(2,2), p,
     2   pd, F(2,2), Jac, Fi(2,2), mu, rho, beta, drho, dbeta, ya_g, Ja,
     3   Siso(2,2), Dm(3,3), tauM, tauC, rC, rCl, Pdev(2,2), DBm(3,2),
     4   Bm(3,2,eNoNw), NxFi(2,eNoNw), VxFi(2,2), VxNx(2,eNoNw), BtDB,
     5   NxSNx, NxNx, Ku, T1, T2, T3, T4, ddev(2,2), Pvis(2,2),
     6   PvNx(2,eNoNw)

!     Define parameters
      mu      = eq(cEq)%dmn(cDmn)%prop(solid_viscosity)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      afm     = af / am

!     {i,j} := velocity dofs; {k} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0._RKIND
      vd = -fb
      vx = 0._RKIND
      F  = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      ya_g   = 0._RKIND
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)

         F(1,1)  = F(1,1) + Nwx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nwx(2,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nwx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nwx(2,a)*dl(j,a)

         ya_g    = ya_g + Nw(a)*ya_l(a)
      END DO
      Jac = MAT_DET(F, 2)
      Fi  = MAT_INV(F, 2)

!     Pressure and its time derivative
      p  = 0._RKIND
      pd = 0._RKIND
      DO a=1, eNoNq
         p  = p  + Nq(a)*yl(k,a)
         pd = pd + Nq(a)*al(k,a)
      END DO

!     Compute deviatoric 2nd Piola-Kirchhoff stress tensor (Siso) and
!     isochoric elasticity tensor in Voigt notation (Dm)
      CALL GETPK2CCdev(eq(cEq)%dmn(cDmn), F, nFn, fN, ya_g, Siso, Dm,Ja)

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(eq(cEq)%dmn(cDmn), p, rho, beta, drho, dbeta, Ja)

!     Compute stabilization parameters
      IF (vmsFlag) THEN
         CALL GETTAU(eq(cEq)%dmn(cDmn), Jac, Je, tauM, tauC)
      ELSE
         tauM = 0._RKIND
         tauC = 0._RKIND
      END IF

!     Deviatoric 1st Piola-Kirchhoff tensor (P)
      Pdev = MATMUL(F, Siso)

!     Viscous contribution
      ddev = 2._RKIND*mu*Jac*MAT_DEV(MAT_SYMM(vx,2), 2)
      Pvis = MATMUL(ddev, TRANSPOSE(Fi))

!     Total 1st Piola-Kirchhoff stress
      Pdev = Pdev + Pvis

      DO a=1, eNoNw
         Bm(1,1,a) = Nwx(1,a)*F(1,1)
         Bm(1,2,a) = Nwx(1,a)*F(2,1)

         Bm(2,1,a) = Nwx(2,a)*F(1,2)
         Bm(2,2,a) = Nwx(2,a)*F(2,2)

         Bm(3,1,a) = (Nwx(1,a)*F(1,2) + F(1,1)*Nwx(2,a))
         Bm(3,2,a) = (Nwx(1,a)*F(2,2) + F(2,1)*Nwx(2,a))
      END DO

      DO a=1, eNoNw
         NxFi(1,a) = Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1)
         NxFi(2,a) = Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2)

         PvNx(1,a) = ddev(1,1)*NxFi(1,a) + ddev(1,2)*NxFi(2,a)
         PvNx(2,a) = ddev(2,1)*NxFi(1,a) + ddev(2,2)*NxFi(2,a)
      END DO

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2)

      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2)

      rC  = beta*pd + VxFi(1,1) + VxFi(2,2)
      rCl = -p + tauC*rC

!     Local residue
      DO a=1, eNoNw
         T1 = Jac*rho*vd(1)*Nw(a)
         T2 = Pdev(1,1)*Nwx(1,a) + Pdev(1,2)*Nwx(2,a)
         T3 = Jac*rCl*NxFi(1,a)
         lR(1,a) = lR(1,a) + w*(T1 + T2 + T3)

         T1 = Jac*rho*vd(2)*Nw(a)
         T2 = Pdev(2,1)*Nwx(1,a) + Pdev(2,2)*Nwx(2,a)
         T3 = Jac*rCl*NxFi(2,a)
         lR(2,a) = lR(2,a) + w*(T1 + T2 + T3)
      END DO

      DO a=1, eNoNw
         VxNx(1,a) = VxFi(1,1)*NxFi(1,a) + VxFi(2,1)*NxFi(2,a)
         VxNx(2,a) = VxFi(1,2)*NxFi(1,a) + VxFi(2,2)*NxFi(2,a)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            NxSNx = Nwx(1,a)*Siso(1,1)*Nwx(1,b)
     2            + Nwx(1,a)*Siso(1,2)*Nwx(2,b)
     3            + Nwx(2,a)*Siso(2,1)*Nwx(1,b)
     4            + Nwx(2,a)*Siso(2,2)*Nwx(2,b)
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

            NxNx = NxFi(1,a)*NxFi(1,b) + NxFi(2,a)*NxFi(2,b)

!           dM1_dV1 + af/am *dM_1/dU_1
            BtDB = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1)
            T1   = Jac*rho*vd(1)*Nw(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(1,b)
            T3   = PvNx(1,a)*NxFi(1,b) - NxFi(1,a)*PvNx(1,b)
            Ku   = w*af*(T1 + T2 + T3 + BtDB + NxSNx)
            lKd(1,a,b) = lKd(1,a,b) + Ku

            T1   = am*Jac*rho*Nw(a)*Nw(b)
            T2   = T1 + af*Jac*tauC*rho*NxFi(1,a)*NxFi(1,b)
            T3   = T2 + af*mu*Jac*NxNx
            lK(1,a,b)  = lK(1,a,b) + w*T3 + afm*Ku

!           dM_1/dV_2 + af/am *dM_1/dU_2
            BtDB = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2             Bm(3,1,a)*DBm(3,2)
            T1   = Jac*rho*vd(1)*Nw(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(1,a)*VxNx(2,b)
            T3   = Jac*rCl*(NxFi(1,a)*NxFi(2,b) - NxFi(2,a)*NxFi(1,b))
            T4   = PvNx(1,a)*NxFi(2,b) - NxFi(2,a)*PvNx(1,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(2,a,b) = lKd(2,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(1,a)*NxFi(2,b)
            T3   = T2 + af*mu*Jac*(NxFi(2,a)*NxFi(1,b)
     2           - NxFi(1,a)*NxFi(2,b))
            lK(2,a,b) = lK(2,a,b) + w*T3 + afm*Ku

!           dM_2/dV_1 + af/am *dM_2/dU_1
            BtDB = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2             Bm(3,2,a)*DBm(3,1)
            T1   = Jac*rho*vd(2)*Nw(a)*NxFi(1,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(1,b)
            T3   = Jac*rCl*(NxFi(2,a)*NxFi(1,b) - NxFi(1,a)*NxFi(2,b))
            T4   = PvNx(2,a)*NxFi(1,b) - NxFi(1,a)*PvNx(2,b)
            Ku   = w*af*(T1 + T2 + T3 + T4 + BtDB)
            lKd(3,a,b) = lKd(3,a,b) + Ku

            T2   = af*Jac*tauC*rho*NxFi(2,a)*NxFi(1,b)
            T3   = T2 + af*mu*Jac*(NxFi(1,a)*NxFi(2,b)
     2           - NxFi(2,a)*NxFi(1,b))
            lK(4,a,b) = lK(4,a,b) + w*T3 + afm*Ku

!           dM_2/dV_2 + af/am *dM_2/dU_2
            BtDB = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2             Bm(3,2,a)*DBm(3,2)
            T1   = Jac*rho*vd(2)*Nw(a)*NxFi(2,b)
            T2   = -tauC*Jac*NxFi(2,a)*VxNx(2,b)
            T3   = PvNx(2,a)*NxFi(2,b) - NxFi(2,a)*PvNx(2,b)
            Ku   = w*af*(T1 + T2 + T3 + BtDB + NxSNx)
            lKd(4,a,b) = lKd(4,a,b) + Ku

            T1   = am*Jac*rho*Nw(a)*Nw(b)
            T2   = T1 + af*Jac*tauC*rho*NxFi(2,a)*NxFi(2,b)
            T3   = T2 + af*mu*Jac*NxNx
            lK(5,a,b) = lK(5,a,b) + w*T3 + afm*Ku
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
!           dM_1/dP
            T1 = am*tauC*beta + af*(tauC*dbeta*pd - 1._RKIND)
            T2 = T1*NxFi(1,a)*Nq(b) + af*drho*vd(1)*Nw(a)*Nq(b)
            lK(3,a,b) = lK(3,a,b) + w*Jac*T2

!           dM_2/dP
            T2 = T1*NxFi(2,a)*Nq(b) + af*drho*vd(2)*Nw(a)*Nq(b)
            lK(6,a,b) = lK(6,a,b) + w*Jac*T2
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT2D_M
!####################################################################
      SUBROUTINE USTRUCT3D_C(vmsFlag, eNoNw, eNoNq, w, Je, Nw, Nq, Nwx,
     2   Nqx, al, yl, dl, bfl, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Je, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), Nqx(3,eNoNq), al(tDof,eNoNw), yl(tDof,eNoNw),
     3   dl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lKd(dof*3,eNoNw,eNoNw), lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, l, a, b
      REAL(KIND=RKIND) :: fb(3), am, af, afm, v(3), vd(3), vx(3,3), p,
     2   pd, px(3), F(3,3), Jac, Fi(3,3), rho, beta, drho, dbeta, tauM,
     3   tauC, rC, rM(3), NwxFi(3,eNoNw), NqxFi(3,eNoNq),VxFi(3,3),
     4   PxFi(3), rMNqx(eNoNq), rMNwx(eNoNw), VxNwx(3,eNoNw), NxNx, T1,
     5   T2, T3, Ku

!     Define parameters
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      afm     = af / am

!     {i,j,k} := velocity dofs; {l} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1
      l       = k + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0._RKIND
      vd = -fb
      vx = 0._RKIND
      F  = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)
         v(3)    = v(3)  + Nw(a)*yl(k,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))
         vd(3)   = vd(3) + Nw(a)*(al(k,a)-bfl(3,a))

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(1,3) = vx(1,3) + Nwx(3,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)
         vx(2,3) = vx(2,3) + Nwx(3,a)*yl(j,a)
         vx(3,1) = vx(3,1) + Nwx(1,a)*yl(k,a)
         vx(3,2) = vx(3,2) + Nwx(2,a)*yl(k,a)
         vx(3,3) = vx(3,3) + Nwx(3,a)*yl(k,a)

         F(1,1)  = F(1,1) + Nwx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nwx(2,a)*dl(i,a)
         F(1,3)  = F(1,3) + Nwx(3,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nwx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nwx(2,a)*dl(j,a)
         F(2,3)  = F(2,3) + Nwx(3,a)*dl(j,a)
         F(3,1)  = F(3,1) + Nwx(1,a)*dl(k,a)
         F(3,2)  = F(3,2) + Nwx(2,a)*dl(k,a)
         F(3,3)  = F(3,3) + Nwx(3,a)*dl(k,a)
      END DO
      Jac = MAT_DET(F, 3)
      Fi  = MAT_INV(F, 3)

!     Pressure and its gradients
      p  = 0._RKIND
      pd = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p       = p     + Nq(a)*yl(l,a)
         pd      = pd    + Nq(a)*al(l,a)
         px(1)   = px(1) + Nqx(1,a)*yl(l,a)
         px(2)   = px(2) + Nqx(2,a)*yl(l,a)
         px(3)   = px(3) + Nqx(3,a)*yl(l,a)
      END DO

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(eq(cEq)%dmn(cDmn), p, rho, beta, drho, dbeta,
     2   1._RKIND)

!     Compute stabilization parameters
      IF (vmsFlag) THEN
         CALL GETTAU(eq(cEq)%dmn(cDmn), Jac, Je, tauM, tauC)
      ELSE
         tauM = 0._RKIND
         tauC = 0._RKIND
      END IF

      DO a=1, eNoNw
         NwxFi(1,a) = Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1) +
     2      Nwx(3,a)*Fi(3,1)
         NwxFi(2,a) = Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2) +
     2      Nwx(3,a)*Fi(3,2)
         NwxFi(3,a) = Nwx(1,a)*Fi(1,3) + Nwx(2,a)*Fi(2,3) +
     2      Nwx(3,a)*Fi(3,3)
      END DO

      DO a=1, eNoNq
         NqxFi(1,a) = Nqx(1,a)*Fi(1,1) + Nqx(2,a)*Fi(2,1) +
     2      Nqx(3,a)*Fi(3,1)
         NqxFi(2,a) = Nqx(1,a)*Fi(1,2) + Nqx(2,a)*Fi(2,2) +
     2      Nqx(3,a)*Fi(3,2)
         NqxFi(3,a) = Nqx(1,a)*Fi(1,3) + Nqx(2,a)*Fi(2,3) +
     2      Nqx(3,a)*Fi(3,3)
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
      rM(1) = rho*vd(1) + PxFi(1)
      rM(2) = rho*vd(2) + PxFi(2)
      rM(3) = rho*vd(3) + PxFi(3)

!     Local residue
      DO a=1, eNoNq
         rMNqx(a) = rM(1)*NqxFi(1,a) + rM(2)*NqxFi(2,a) +
     2      rM(3)*NqxFi(3,a)
         lR(4,a) = lR(4,a) + w*Jac*(Nq(a)*rC + tauM*rMNqx(a))
      END DO

      DO a=1, eNoNw
         rMNwx(a) = rM(1)*NwxFi(1,a) + rM(2)*NwxFi(2,a) +
     2      rM(3)*NwxFi(3,a)

         VxNwx(1,a) = VxFi(1,1)*NwxFi(1,a) + VxFi(2,1)*NwxFi(2,a) +
     2      VxFi(3,1)*NwxFi(3,a)
         VxNwx(2,a) = VxFi(1,2)*NwxFi(1,a) + VxFi(2,2)*NwxFi(2,a) +
     2      VxFi(3,2)*NwxFi(3,a)
         VxNwx(3,a) = VxFi(1,3)*NwxFi(1,a) + VxFi(2,3)*NwxFi(2,a) +
     2      VxFi(3,3)*NwxFi(3,a)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNq
            NxNx  = NqxFi(1,a)*NwxFi(1,b) + NqxFi(2,a)*NwxFi(2,b) +
     2         NqxFi(3,a)*NwxFi(3,b)

!           dC/dV_1 + af/am *dC/dU_1
            T1   = Nq(a)*(rC*NwxFi(1,b) - VxNwx(1,b))
            T2   = tauM*(rMNqx(a)*NwxFi(1,b) - rMNwx(b)*NqxFi(1,a))
            T3   = -tauM*NxNx*PxFi(1)
            Ku   = w*af*Jac*(T1 + T2 + T3)
            lKd(10,a,b) = lKd(10,a,b) + Ku

            T2   = (am*tauM*rho)*NqxFi(1,a)*Nw(b) + af*Nq(a)*NwxFi(1,b)
            lK(13,a,b) = lK(13,a,b) + w*Jac*T2 + afm*Ku

!           dC/dV_2 + af/am *dC/dU_2
            T1   = Nq(a)*(rC*NwxFi(2,b) - VxNwx(2,b))
            T2   = tauM*(rMNqx(a)*NwxFi(2,b) - rMNwx(b)*NqxFi(2,a))
            T3   = -tauM*NxNx*PxFi(2)
            Ku   = w*af*Jac*(T1 + T2 + T3)
            lKd(11,a,b) = lKd(11,a,b) + Ku

            T2 = (am*tauM*rho)*NqxFi(2,a)*Nw(b) + af*Nq(a)*NwxFi(2,b)
            lK(14,a,b) = lK(14,a,b) + w*Jac*T2 + afm*Ku

!           dC/dV_3 + af/am *dC/dU_3
            T1   = Nq(a)*(rC*NwxFi(3,b) - VxNwx(3,b))
            T2   = tauM*(rMNqx(a)*NwxFi(3,b) - rMNwx(b)*NqxFi(3,a))
            T3   = -tauM*NxNx*PxFi(3)
            Ku   = w*af*Jac*(T1 + T2 + T3)
            lKd(12,a,b) = lKd(12,a,b) + Ku

            T2 = (am*tauM*rho)*NqxFi(3,a)*Nw(b) + af*Nq(a)*NwxFi(3,b)
            lK(15,a,b) = lK(15,a,b) + w*Jac*T2 + afm*Ku
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNq
!           dC/dP
            NxNx = NqxFi(1,a)*NqxFi(1,b) + NqxFi(2,a)*NqxFi(2,b)
     2           + NqxFi(3,a)*NqxFi(3,b)

            T1 = (am*beta + af*dbeta*pd)*Nq(a)*Nq(b)
            T2 = NqxFi(1,a)*vd(1) + NqxFi(2,a)*vd(2) + NqxFi(3,a)*vd(3)
            T3 = T1 + af*tauM*(NxNx + drho*T2*Nq(b))
            lK(16,a,b) = lK(16,a,b) + w*Jac*T3
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT3D_C
!--------------------------------------------------------------------
      SUBROUTINE USTRUCT2D_C(vmsFlag, eNoNw, eNoNq, w, Je, Nw, Nq, Nwx,
     2   Nqx, al, yl, dl, bfl, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Je, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), Nqx(2,eNoNq), al(tDof,eNoNw), yl(tDof,eNoNw),
     3   dl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lKd(dof*2,eNoNw,eNoNw), lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: fb(2), am, af, afm, v(2), vd(2), vx(2,2), p,
     2   pd, px(2), F(2,2), Jac, Fi(2,2), rho, beta, drho, dbeta, tauM,
     3   tauC, rC, rM(2), NwxFi(2,eNoNw), NqxFi(2,eNoNq),VxFi(2,2),
     4   PxFi(2), rMNqx(eNoNq), rMNwx(eNoNw), VxNwx(2,eNoNw), NxNx, T1,
     5   T2, T3, Ku

!     Define parameters
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      am      = eq(cEq)%am
      af      = eq(cEq)%af*eq(cEq)%gam*dt
      afm     = af / am

!     {i,j} := velocity dofs; {k} := pressure dof
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     Inertia (velocity and acceleration), body force, fiber directions,
!     and deformation tensor (F) at integration point
      v  = 0._RKIND
      vd = -fb
      vx = 0._RKIND
      F  = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)

         F(1,1)  = F(1,1) + Nwx(1,a)*dl(i,a)
         F(1,2)  = F(1,2) + Nwx(2,a)*dl(i,a)
         F(2,1)  = F(2,1) + Nwx(1,a)*dl(j,a)
         F(2,2)  = F(2,2) + Nwx(2,a)*dl(j,a)
      END DO
      Jac = MAT_DET(F, 2)
      Fi  = MAT_INV(F, 2)

!     Pressure and its gradients
      p  = 0._RKIND
      pd = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p       = p     + Nq(a)*yl(k,a)
         pd      = pd    + Nq(a)*al(k,a)
         px(1)   = px(1) + Nqx(1,a)*yl(k,a)
         px(2)   = px(2) + Nqx(2,a)*yl(k,a)
      END DO

!     Compute rho and beta depending on the volumetric penalty model
      CALL GVOLPEN(eq(cEq)%dmn(cDmn), p, rho, beta, drho, dbeta,
     2   1._RKIND)

!     Compute stabilization parameters
      IF (vmsFlag) THEN
         CALL GETTAU(eq(cEq)%dmn(cDmn), Jac, Je, tauM, tauC)
      ELSE
         tauM = 0._RKIND
         tauC = 0._RKIND
      END IF

      DO a=1, eNoNw
         NwxFi(1,a) = Nwx(1,a)*Fi(1,1) + Nwx(2,a)*Fi(2,1)
         NwxFi(2,a) = Nwx(1,a)*Fi(1,2) + Nwx(2,a)*Fi(2,2)
      END DO

      DO a=1, eNoNq
         NqxFi(1,a) = Nqx(1,a)*Fi(1,1) + Nqx(2,a)*Fi(2,1)
         NqxFi(2,a) = Nqx(1,a)*Fi(1,2) + Nqx(2,a)*Fi(2,2)
      END DO

      VxFi(1,1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1)
      VxFi(1,2) = vx(1,1)*Fi(1,2) + vx(1,2)*Fi(2,2)

      VxFi(2,1) = vx(2,1)*Fi(1,1) + vx(2,2)*Fi(2,1)
      VxFi(2,2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2)

      PxFi(1)   = px(1)*Fi(1,1) + px(2)*Fi(2,1)
      PxFi(2)   = px(1)*Fi(1,2) + px(2)*Fi(2,2)

      rC    = beta*pd + VxFi(1,1) + VxFi(2,2)
      rM(1) = rho*vd(1) + PxFi(1)
      rM(2) = rho*vd(2) + PxFi(2)

!     Local residue
      DO a=1, eNoNq
         rMNqx(a) = rM(1)*NqxFi(1,a) + rM(2)*NqxFi(2,a)
         lR(3,a) = lR(3,a) + w*Jac*(Nq(a)*rC + tauM*rMNqx(a))
      END DO

      DO a=1, eNoNw
         rMNwx(a) = rM(1)*NwxFi(1,a) + rM(2)*NwxFi(2,a)

         VxNwx(1,a) = VxFi(1,1)*NwxFi(1,a) + VxFi(2,1)*NwxFi(2,a)
         VxNwx(2,a) = VxFi(1,2)*NwxFi(1,a) + VxFi(2,2)*NwxFi(2,a)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNq
            NxNx  = NqxFi(1,a)*NwxFi(1,b) + NqxFi(2,a)*NwxFi(2,b)

!           dC/dV_1 + af/am *dC/dU_1
            T1   = Nq(a)*(rC*NwxFi(1,b) - VxNwx(1,b))
            T2   = tauM*(rMNqx(a)*NwxFi(1,b) - rMNwx(b)*NqxFi(1,a))
            T3   = -tauM*NxNx*PxFi(1)
            Ku   = w*af*Jac*(T1 + T2 + T3)
            lKd(5,a,b) = lKd(5,a,b) + Ku

            T2   = (am*tauM*rho)*NqxFi(1,a)*Nw(b) + af*Nq(a)*NwxFi(1,b)
            lK(7,a,b) = lK(7,a,b) + w*Jac*T2 + afm*Ku

!           dC/dV_2 + af/am *dC/dU_2
            T1   = Nq(a)*(rC*NwxFi(2,b) - VxNwx(2,b))
            T2   = tauM*(rMNqx(a)*NwxFi(2,b) - rMNwx(b)*NqxFi(2,a))
            T3   = -tauM*NxNx*PxFi(2)
            Ku   = w*af*Jac*(T1 + T2 + T3)
            lKd(6,a,b) = lKd(6,a,b) + Ku

            T2 = (am*tauM*rho)*NqxFi(2,a)*Nw(b) + af*Nq(a)*NwxFi(2,b)
            lK(8,a,b) = lK(8,a,b) + w*Jac*T2 + afm*Ku
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNq
!           dC/dP
            NxNx = NqxFi(1,a)*NqxFi(1,b) + NqxFi(2,a)*NqxFi(2,b)

            T1 = (am*beta + af*dbeta*pd)*Nq(a)*Nq(b)
            T2 = NqxFi(1,a)*vd(1) + NqxFi(2,a)*vd(2)
            T3 = T1 + af*tauM*(NxNx + drho*T2*Nq(b))
            lK(9,a,b) = lK(9,a,b) + w*Jac*T3
         END DO
      END DO

      RETURN
      END SUBROUTINE USTRUCT2D_C
!####################################################################
      SUBROUTINE BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(3)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lKd(dof*3,eNoN,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: af, afm, h, Jac, wl, Ku, F(3,3), Fi(3,3),
     2   nFi(3), NxFi(3,eNoN)

      af     = eq(cEq)%af*eq(cEq)%gam*dt
      afm    = af / eq(cEq)%am
      i      = eq(cEq)%s
      j      = i + 1
      k      = j + 1

      h      = 0._RKIND
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
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

      wl = w*Jac*h
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)
         lR(3,a) = lR(3,a) - wl*N(a)*nFi(3)

         DO b=1, eNoN
            Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lKd(2,a,b) = lKd(2,a,b) + Ku
            lK(2,a,b)  = lK(2,a,b)  + afm*Ku

            lKd(4,a,b) = lKd(4,a,b) - Ku
            lK(5,a,b)  = lK(5,a,b)  - afm*Ku

            Ku = wl*af*N(a)*(nFi(3)*NxFi(1,b) - nFi(1)*NxFi(3,b))
            lKd(3,a,b) = lKd(3,a,b) + Ku
            lK(3,a,b)  = lK(3,a,b)  + afm*Ku

            lKd(7,a,b) = lKd(7,a,b) - Ku
            lK(9,a,b)  = lK(9,a,b)  - afm*Ku

            Ku = wl*af*N(a)*(nFi(3)*NxFi(2,b) - nFi(2)*NxFi(3,b))
            lKd(6,a,b) = lKd(6,a,b) + Ku
            lK(7,a,b)  = lK(7,a,b)  + afm*Ku

            lKd(8,a,b) = lKd(8,a,b) - Ku
            lK(10,a,b) = lK(10,a,b) - afm*Ku
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT3D
!--------------------------------------------------------------------
      SUBROUTINE BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), hl(eNoN),
     2   dl(tDof,eNoN), nV(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lKd(dof*2,eNoN,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: i, j, a, b
      REAL(KIND=RKIND) :: af, afm, h, Jac, wl, Ku, F(2,2), Fi(2,2),
     2   nFi(2), NxFi(2,eNoN)

      af     = eq(cEq)%af*eq(cEq)%gam*dt
      afm    = af / eq(cEq)%am
      i      = eq(cEq)%s
      j      = i + 1

      h      = 0._RKIND
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
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

      wl = w*Jac*h
      DO a=1, eNoN
         lR(1,a) = lR(1,a) - wl*N(a)*nFi(1)
         lR(2,a) = lR(2,a) - wl*N(a)*nFi(2)

         DO b=1, eNoN
            Ku = wl*af*N(a)*(nFi(2)*NxFi(1,b) - nFi(1)*NxFi(2,b))
            lKd(2,a,b) = lKd(2,a,b) + Ku
            lK(2,a,b)  = lK(2,a,b)  + afm*Ku

            lKd(3,a,b) = lKd(3,a,b) - Ku
            lK(4,a,b)  = lK(4,a,b)  - afm*Ku
         END DO
      END DO

      RETURN
      END SUBROUTINE BUSTRUCT2D
!####################################################################
      SUBROUTINE USTRUCT_DOASSEM(d, eqN, lKd, lK, lR)
      USE TYPEMOD
      USE COMMOD, ONLY: dof, nsd, rowPtr, colPtr, idMap, Kd, Val, R
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: d, eqN(d)
      REAL(KIND=RKIND), INTENT(IN) :: lKd(dof*nsd,d,d), lK(dof*dof,d,d),
     2   lR(dof,d)

      INTEGER(KIND=IKIND) a, b, ptr, rowN, colN

      DO a=1, d
!        Momentum equation residue is assembled at mapped rows
         rowN = idMap(eqn(a))
         R(1:nsd,rowN) = R(1:nsd,rowN) + lR(1:nsd,a)

!        Continuity equation residue is assembled at unmapped rows
         rowN = eqn(a)
         R(nsd+1,rowN) = R(nsd+1,rowN) + lR(nsd+1,a)
      END DO

      IF (nsd .EQ. 3) THEN
!        Stiffness matrix (A) is assembled using mapped rows and columns
!        Gradient matrix (B) is assembled using mapped row but unmapped
!        column
         DO a=1, d
            rowN = idMap(eqN(a))
!           A - matrix
            DO b=1, d
               colN = idMap(eqN(b))
               CALL GETCOLPTR()
               Kd (1 :9 ,ptr) = Kd (1 :9 ,ptr) + lKd(1 :9 ,a,b)
               Val(1 :3 ,ptr) = Val(1 :3 ,ptr) + lK (1 :3 ,a,b)
               Val(5 :7 ,ptr) = Val(5 :7 ,ptr) + lK (5 :7 ,a,b)
               Val(9 :11,ptr) = Val(9 :11,ptr) + lK (9 :11,a,b)
            END DO

!           B - matrix
            DO b=1, d
               colN = eqN(b)
               CALL GETCOLPTR()
               Val(4 ,ptr) = Val(4 ,ptr) + lK(4 ,a,b)
               Val(8 ,ptr) = Val(8 ,ptr) + lK(8 ,a,b)
               Val(12,ptr) = Val(12,ptr) + lK(12,a,b)
            END DO
         END DO

!        Divergence matrix (C) is assembled using unmapped rows but
!        mapped columns. Mass matrix (D) is assembled using unmapped
!        rows and columns
         DO a=1, d
            rowN = eqN(a)
!           C - matrix
            DO b=1, d
               colN = idMap(eqN(b))
               CALL GETCOLPTR()
               Kd (10:12,ptr) = Kd (10:12,ptr) + lKd(10:12,a,b)
               Val(13:15,ptr) = Val(13:15,ptr) + lK (13:15,a,b)
            END DO

!           D - matrix
            DO b=1, d
               colN = eqN(b)
               CALL GETCOLPTR()
               Val(16,ptr) = Val(16,ptr) + lK(16,a,b)
            END DO
         END DO

      ELSE
!        Stiffness matrix (A) is assembled using mapped rows and columns
!        Gradient matrix (B) is assembled using mapped row but unmapped
!        column
         DO a=1, d
            rowN = idMap(eqN(a))
!           A - matrix
            DO b=1, d
               colN = idMap(eqN(b))
               CALL GETCOLPTR()
               Kd (1 :4 ,ptr) = Kd (1 :4 ,ptr) + lKd(1 :4 ,a,b)
               Val(1 :2 ,ptr) = Val(1 :2 ,ptr) + lK (1 :2 ,a,b)
               Val(4 :5 ,ptr) = Val(4 :5 ,ptr) + lK (4 :5 ,a,b)
            END DO

!           B - matrix
            DO b=1, d
               colN = eqN(b)
               CALL GETCOLPTR()
               Val(3,ptr) = Val(3,ptr) + lK(3,a,b)
               Val(6,ptr) = Val(6,ptr) + lK(6,a,b)
            END DO
         END DO

!        Divergence matrix (C) is assembled using unmapped rows but
!        mapped columns. Mass matrix (D) is assembled using unmapped
!        rows and columns
         DO a=1, d
            rowN = eqN(a)
!           C - matrix
            DO b=1, d
               colN = idMap(eqN(b))
               CALL GETCOLPTR()
               Kd (5:6,ptr) = Kd (5:6,ptr) + lKd(5:6,a,b)
               Val(7:8,ptr) = Val(7:8,ptr) + lK (7:8,a,b)
            END DO

!           D - matrix
            DO b=1, d
               colN = eqN(b)
               CALL GETCOLPTR()
               Val(9,ptr) = Val(9,ptr) + lK(9,a,b)
            END DO
         END DO
      END IF

      RETURN
      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE GETCOLPTR()
         IMPLICIT NONE

         INTEGER(KIND=IKIND) left, right

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

         RETURN
         END SUBROUTINE GETCOLPTR
!--------------------------------------------------------------------
      END SUBROUTINE USTRUCT_DOASSEM
!####################################################################
      SUBROUTINE USTRUCTR(Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, i, c, s
      REAL(KIND=RKIND) :: amg, ami

      REAL(KIND=RKIND), ALLOCATABLE :: KU(:,:)

      IF (eq(cEq)%phys .NE. phys_ustruct .AND.
     2    eq(cEq)%phys .NE. phys_FSI) RETURN

      s   = eq(cEq)%s
      amg = (eq(cEq)%gam-eq(cEq)%am) / (eq(cEq)%gam-1._RKIND)
      ami = 1._RKIND/eq(cEq)%am

      IF (eq(cEq)%itr .GT. 1) THEN
         Rd = 0._RKIND
      ELSE
         DO a=1, tnNo
            IF (.NOT.ISDOMAIN(cEq, a, phys_ustruct)) CYCLE
            DO i=1, nsd
               Rd(i,a) = amg*Ad(i,a) - Yg(s+i-1,a)
            END DO
         END DO

         IF (nsd .EQ. 3) THEN
            ALLOCATE(KU(4,tnNo))
            KU = 0._RKIND
            DO a=1, tnNo
               IF (.NOT.ISDOMAIN(cEq, a, phys_ustruct)) CYCLE
               DO i=rowPtr(a), rowPtr(a+1)-1
                  c = colPtr(i)
                  KU(1,a) = KU(1,a) + Kd(1 ,i)*Rd(1,c) +
     2               Kd(2 ,i)*Rd(2,c) + Kd(3 ,i)*Rd(3,c)
                  KU(2,a) = KU(2,a) + Kd(4 ,i)*Rd(1,c) +
     2               Kd(5 ,i)*Rd(2,c) + Kd(6 ,i)*Rd(3,c)
                  KU(3,a) = KU(3,a) + Kd(7 ,i)*Rd(1,c) +
     2               Kd(8 ,i)*Rd(2,c) + Kd(9 ,i)*Rd(3,c)
                  KU(4,a) = KU(4,a) + Kd(10,i)*Rd(1,c) +
     2               Kd(11,i)*Rd(2,c) + Kd(12,i)*Rd(3,c)
               END DO
            END DO

            CALL COMMU(KU)

            DO a=1, tnNo
               R(1,a) = R(1,a) - ami*KU(1,a)
               R(2,a) = R(2,a) - ami*KU(2,a)
               R(3,a) = R(3,a) - ami*KU(3,a)
               R(4,a) = R(4,a) - ami*KU(4,a)
            END DO
            DEALLOCATE(KU)
         ELSE
            ALLOCATE(KU(3,tnNo))
            KU = 0._RKIND
            DO a=1, tnNo
               IF (.NOT.ISDOMAIN(cEq, a, phys_ustruct)) CYCLE
               DO i=rowPtr(a), rowPtr(a+1)-1
                  c = colPtr(i)
                  KU(1,a) = KU(1,a) + Kd(1,i)*Rd(1,c) + Kd(2,i)*Rd(2,c)
                  KU(2,a) = KU(2,a) + Kd(3,i)*Rd(1,c) + Kd(4,i)*Rd(2,c)
                  KU(3,a) = KU(3,a) + Kd(5,i)*Rd(1,c) + Kd(6,i)*Rd(2,c)
               END DO
            END DO

            CALL COMMU(KU)

            DO a=1, tnNo
               R(1,a) = R(1,a) - ami*KU(1,a)
               R(2,a) = R(2,a) - ami*KU(2,a)
               R(3,a) = R(3,a) - ami*KU(3,a)
            END DO
            DEALLOCATE(KU)
         END IF
      END  IF

      RETURN
      END SUBROUTINE USTRUCTR
!####################################################################

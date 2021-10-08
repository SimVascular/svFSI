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
!     This is for solving fluid transport equation solving Navier-Stokes
!     equations. Dirichlet boundary conditions are either treated
!     strongly or weakly.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_FLUID(lM, Ag, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)

      LOGICAL :: vmsStab
      INTEGER(KIND=IKIND) a, e, g, l, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:)

      eNoN = lM%eNoN

      IF (lM%nFs .EQ. 1) THEN
         vmsStab = .TRUE.
      ELSE
         vmsStab = .FALSE.
      END IF

!     l = 3, if nsd==2 ; else 6;
      l = nsymd

!     FLUID: dof = nsd+1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   bfl(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_fluid) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
         END DO

!        Initialize residue and tangents
         lR = 0._RKIND
         lK = 0._RKIND

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN),
     2      Nwxx(l,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND
         Nqx      = 0._RKIND
         Nwxx     = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g),
     2            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
            END IF
            w = fs(1)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)

             ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
            END IF
         END DO ! g: loop

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, vmsStab, 2)

!        Gauss integration 2
         DO g=1, fs(2)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(2)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_C(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nwxx, Nqx)

!        Assembly
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, bfl, lR, lK)

      CALL DESTROY(fs(1))
      CALL DESTROY(fs(2))

      RETURN
      END SUBROUTINE CONSTRUCT_FLUID
!####################################################################
      SUBROUTINE FLUID3D_M(vmsFlag, eNoNw, eNoNq, w, Kxi, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Kxi(3,3), Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), Nqx(3,eNoNq), Nwxx(6,eNoNw), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) a, b, k
      REAL(KIND=RKIND) ctM, ctC, amd, wl, wr, rho, tauM, tauC, tauB, kT,
     2   kS, kU, divU, gam, mu, mu_s, mu_g, p, pa, u(3), ud(3), px(3),
     3   f(3), up(3), ua(3), ux(3,3), uxx(3,3,3), es(3,3), es_x(3,3,3),
     4   esNx(3,eNoNw), mu_x(3), rV(3), rS(3), rM(3,3), updu(3,3,eNoNw),
     5   uNx(eNoNw), upNx(eNoNw), uaNx(eNoNw), d2u2(3), NxNx, T1, T2

      ctM  = 1._RKIND
      ctC  = 36._RKIND

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1
      wr   = w*rho

!     Note that indices are not selected based on the equation because
!     fluid equation always come first
!     Velocity and its gradients, inertia (acceleration & body force)
      ud  = -f
      u   = 0._RKIND
      ux  = 0._RKIND
      uxx = 0._RKIND
      DO a=1, eNoNw
         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
         ud(3)   = ud(3) + Nw(a)*(al(3,a)-bfl(3,a))

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)
         u(3)    = u(3) + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nwx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nwx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nwx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nwx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nwx(3,a)*yl(3,a)

         uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a)
         uxx(2,1,2) = uxx(2,1,2) + Nwxx(2,a)*yl(1,a)
         uxx(3,1,3) = uxx(3,1,3) + Nwxx(3,a)*yl(1,a)
         uxx(2,1,1) = uxx(2,1,1) + Nwxx(4,a)*yl(1,a)
         uxx(3,1,2) = uxx(3,1,2) + Nwxx(5,a)*yl(1,a)
         uxx(1,1,3) = uxx(1,1,3) + Nwxx(6,a)*yl(1,a)

         uxx(1,2,1) = uxx(1,2,1) + Nwxx(1,a)*yl(2,a)
         uxx(2,2,2) = uxx(2,2,2) + Nwxx(2,a)*yl(2,a)
         uxx(3,2,3) = uxx(3,2,3) + Nwxx(3,a)*yl(2,a)
         uxx(2,2,1) = uxx(2,2,1) + Nwxx(4,a)*yl(2,a)
         uxx(3,2,2) = uxx(3,2,2) + Nwxx(5,a)*yl(2,a)
         uxx(1,2,3) = uxx(1,2,3) + Nwxx(6,a)*yl(2,a)

         uxx(1,3,1) = uxx(1,3,1) + Nwxx(1,a)*yl(3,a)
         uxx(2,3,2) = uxx(2,3,2) + Nwxx(2,a)*yl(3,a)
         uxx(3,3,3) = uxx(3,3,3) + Nwxx(3,a)*yl(3,a)
         uxx(2,3,1) = uxx(2,3,1) + Nwxx(4,a)*yl(3,a)
         uxx(3,3,2) = uxx(3,3,2) + Nwxx(5,a)*yl(3,a)
         uxx(1,3,3) = uxx(1,3,3) + Nwxx(6,a)*yl(3,a)
      END DO
      divU = ux(1,1) + ux(2,2) + ux(3,3)

      uxx(1,1,2) = uxx(2,1,1)
      uxx(2,1,3) = uxx(3,1,2)
      uxx(3,1,1) = uxx(1,1,3)

      uxx(1,2,2) = uxx(2,2,1)
      uxx(2,2,3) = uxx(3,2,2)
      uxx(3,2,1) = uxx(1,2,3)

      uxx(1,3,2) = uxx(2,3,1)
      uxx(2,3,3) = uxx(3,3,2)
      uxx(3,3,1) = uxx(1,3,3)

      d2u2(1) = uxx(1,1,1) + uxx(2,1,2) + uxx(3,1,3)
      d2u2(2) = uxx(1,2,1) + uxx(2,2,2) + uxx(3,2,3)
      d2u2(3) = uxx(1,3,1) + uxx(2,3,2) + uxx(3,3,3)

!     Pressure and its gradient
      p  = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p  = p + Nq(a)*yl(4,a)

         px(1) = px(1) + Nqx(1,a)*yl(4,a)
         px(2) = px(2) + Nqx(2,a)*yl(4,a)
         px(3) = px(3) + Nqx(3,a)*yl(4,a)
      END DO

!     Update convection velocity relative to mesh velocity
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(5,a)
            u(2) = u(2) - Nw(a)*yl(6,a)
            u(3) = u(3) - Nw(a)*yl(7,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)
      es(2,1) = ux(2,1) + ux(1,2)
      es(3,2) = ux(3,2) + ux(2,3)
      es(1,3) = ux(1,3) + ux(3,1)
      es(1,2) = es(2,1)
      es(2,3) = es(3,2)
      es(3,1) = es(1,3)

      DO a=1, eNoNw
        esNx(1,a) = es(1,1)*Nwx(1,a) + es(2,1)*Nwx(2,a)
     2            + es(3,1)*Nwx(3,a)
        esNx(2,a) = es(1,2)*Nwx(1,a) + es(2,2)*Nwx(2,a)
     2            + es(3,2)*Nwx(3,a)
        esNx(3,a) = es(1,3)*Nwx(1,a) + es(2,3)*Nwx(2,a)
     2            + es(3,3)*Nwx(3,a)
      END DO

      DO k=1, 3
         es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k)
         es_x(2,2,k) = uxx(2,2,k) + uxx(2,2,k)
         es_x(3,3,k) = uxx(3,3,k) + uxx(3,3,k)
         es_x(2,1,k) = uxx(2,1,k) + uxx(1,2,k)
         es_x(3,2,k) = uxx(3,2,k) + uxx(2,3,k)
         es_x(1,3,k) = uxx(1,3,k) + uxx(3,1,k)

         es_x(1,2,k) = es_x(2,1,k)
         es_x(2,3,k) = es_x(3,2,k)
         es_x(3,1,k) = es_x(1,3,k)
      END DO

      mu_x(1) = (es_x(1,1,1)*es(1,1) + es_x(2,2,1)*es(2,2)
     2        +  es_x(3,3,1)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,1)*es(2,1) + es_x(3,2,1)*es(3,2)
     4        +  es_x(1,3,1)*es(1,3)

      mu_x(2) = (es_x(1,1,2)*es(1,1) + es_x(2,2,2)*es(2,2)
     2        +  es_x(3,3,2)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,2)*es(2,1) + es_x(3,2,2)*es(3,2)
     4        +  es_x(1,3,2)*es(1,3)

      mu_x(3) = (es_x(1,1,3)*es(1,1) + es_x(2,2,3)*es(2,2)
     2        +  es_x(3,3,3)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,3)*es(2,1) + es_x(3,2,3)*es(3,2)
     4        +  es_x(1,3,3)*es(1,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1) + es(3,1)*es(3,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2) + es(3,2)*es(3,2)
     3    + es(1,3)*es(1,3) + es(2,3)*es(2,3) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_g := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

      IF (ISZERO(gam)) THEN
         mu_g = 0._RKIND
      ELSE
         mu_g = mu_g/gam
      END IF
      mu_x(:) = mu_g * mu_x(:)

!     Stabilization parameters
      kT = 4._RKIND*(ctM/dt)**2._RKIND

      kU = u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1) + u(3)*u(1)*Kxi(3,1)
     2   + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2) + u(3)*u(2)*Kxi(3,2)
     3   + u(1)*u(3)*Kxi(1,3) + u(2)*u(3)*Kxi(2,3) + u(3)*u(3)*Kxi(3,3)

      kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1) + Kxi(3,1)*Kxi(3,1)
     2   + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2) + Kxi(3,2)*Kxi(3,2)
     3   + Kxi(1,3)*Kxi(1,3) + Kxi(2,3)*Kxi(2,3) + Kxi(3,3)*Kxi(3,3)
      kS = ctC * kS * (mu/rho)**2._RKIND

      tauM = 1._RKIND / (rho * SQRT( kT + kU + kS ))

      rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + u(3)*ux(3,1)
      rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + u(3)*ux(3,2)
      rV(3) = ud(3) + u(1)*ux(1,3) + u(2)*ux(2,3) + u(3)*ux(3,3)

      rS(1) = mu_x(1)*es(1,1) + mu_x(2)*es(2,1) + mu_x(3)*es(3,1)
     2      + mu*d2u2(1)
      rS(2) = mu_x(1)*es(1,2) + mu_x(2)*es(2,2) + mu_x(3)*es(3,2)
     2      + mu*d2u2(2)
      rS(3) = mu_x(1)*es(1,3) + mu_x(2)*es(2,3) + mu_x(3)*es(3,3)
     2      + mu*d2u2(3)

      up(1) = -tauM*(rho*rV(1) + px(1) - rS(1))
      up(2) = -tauM*(rho*rV(2) + px(2) - rS(2))
      up(3) = -tauM*(rho*rV(3) + px(3) - rS(3))

      IF (vmsFlag) THEN
         tauC = 1._RKIND / (tauM * (Kxi(1,1) + Kxi(2,2) + Kxi(3,3)))
         tauB = up(1)*up(1)*Kxi(1,1) + up(2)*up(1)*Kxi(2,1)
     2        + up(3)*up(1)*Kxi(3,1) + up(1)*up(2)*Kxi(1,2)
     3        + up(2)*up(2)*Kxi(2,2) + up(3)*up(2)*Kxi(3,2)
     4        + up(1)*up(3)*Kxi(1,3) + up(2)*up(3)*Kxi(2,3)
     5        + up(3)*up(3)*Kxi(3,3)
         IF (ISZERO(tauB)) tauB = eps
         tauB = rho/SQRT(tauB)

         ua(1) = u(1) + up(1)
         ua(2) = u(2) + up(2)
         ua(3) = u(3) + up(3)
         pa    = p - tauC*divU
      ELSE
         tauC  = 0._RKIND
         tauB  = 0._RKIND
         ua(:) = u(:)
         pa    = p
      END IF

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      rV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))

      rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa
      rM(2,1) = mu*es(2,1) - rho*up(1)*ua(2) + rV(1)*up(2)
      rM(3,1) = mu*es(3,1) - rho*up(1)*ua(3) + rV(1)*up(3)

      rM(1,2) = mu*es(1,2) - rho*up(2)*ua(1) + rV(2)*up(1)
      rM(2,2) = mu*es(2,2) - rho*up(2)*ua(2) + rV(2)*up(2) - pa
      rM(3,2) = mu*es(3,2) - rho*up(2)*ua(3) + rV(2)*up(3)

      rM(1,3) = mu*es(1,3) - rho*up(3)*ua(1) + rV(3)*up(1)
      rM(2,3) = mu*es(2,3) - rho*up(3)*ua(2) + rV(3)*up(2)
      rM(3,3) = mu*es(3,3) - rho*up(3)*ua(3) + rV(3)*up(3) - pa

      rV(1) = ud(1) + ua(1)*ux(1,1) + ua(2)*ux(2,1) + ua(3)*ux(3,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*ux(2,2) + ua(3)*ux(3,2)
      rV(3) = ud(3) + ua(1)*ux(1,3) + ua(2)*ux(2,3) + ua(3)*ux(3,3)

!     Local residue
      DO a=1, eNoNw
         lR(1,a) = lR(1,a) + wr*Nw(a)*rV(1) + w*(Nwx(1,a)*rM(1,1)
     2      + Nwx(2,a)*rM(2,1) + Nwx(3,a)*rM(3,1))

         lR(2,a) = lR(2,a) + wr*Nw(a)*rV(2) + w*(Nwx(1,a)*rM(1,2)
     2      + Nwx(2,a)*rM(2,2) + Nwx(3,a)*rM(3,2))

         lR(3,a) = lR(3,a) + wr*Nw(a)*rV(3) + w*(Nwx(1,a)*rM(1,3)
     2      + Nwx(2,a)*rM(2,3) + Nwx(3,a)*rM(3,3))

!        Quantities used for stiffness matrix
         uNx(a)  = u(1)*Nwx(1,a)  + u(2)*Nwx(2,a)  + u(3)*Nwx(3,a)
         upNx(a) = up(1)*Nwx(1,a) + up(2)*Nwx(2,a) + up(3)*Nwx(3,a)
         IF (vmsFlag) THEN
            uaNx(a) = uNx(a) + upNx(a)
         ELSE
            uaNx(a) = uNx(a)
         END IF

         T1 = -rho*uNx(a) + mu*(Nwxx(1,a) + Nwxx(2,a) + Nwxx(3,a))
     2      + mu_x(1)*Nwx(1,a) + mu_x(2)*Nwx(2,a) + mu_x(3)*Nwx(3,a)

         updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1
         updu(2,1,a) = mu_x(2)*Nwx(1,a) + d2u2(1)*mu_g*esNx(2,a)
         updu(3,1,a) = mu_x(3)*Nwx(1,a) + d2u2(1)*mu_g*esNx(3,a)

         updu(1,2,a) = mu_x(1)*Nwx(2,a) + d2u2(2)*mu_g*esNx(1,a)
         updu(2,2,a) = mu_x(2)*Nwx(2,a) + d2u2(2)*mu_g*esNx(2,a) + T1
         updu(3,2,a) = mu_x(3)*Nwx(2,a) + d2u2(2)*mu_g*esNx(3,a)

         updu(1,3,a) = mu_x(1)*Nwx(3,a) + d2u2(3)*mu_g*esNx(1,a)
         updu(2,3,a) = mu_x(2)*Nwx(3,a) + d2u2(3)*mu_g*esNx(2,a)
         updu(3,3,a) = mu_x(3)*Nwx(3,a) + d2u2(3)*mu_g*esNx(3,a) + T1
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            rM(1,1) = Nwx(1,a)*Nwx(1,b)
            rM(2,1) = Nwx(2,a)*Nwx(1,b)
            rM(3,1) = Nwx(3,a)*Nwx(1,b)
            rM(1,2) = Nwx(1,a)*Nwx(2,b)
            rM(2,2) = Nwx(2,a)*Nwx(2,b)
            rM(3,2) = Nwx(3,a)*Nwx(2,b)
            rM(1,3) = Nwx(1,a)*Nwx(3,b)
            rM(2,3) = Nwx(2,a)*Nwx(3,b)
            rM(3,3) = Nwx(3,a)*Nwx(3,b)

            NxNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b)
     2           + Nwx(3,a)*Nwx(3,b)

            T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx(a))
     2         + rho*Nw(a)*(uNx(b)+upNx(b)) + tauB*upNx(a)*upNx(b)

!           dRm_a1/du_b1
            T2 = (mu + tauC)*rM(1,1) + esNx(1,a)*mu_g*esNx(1,b)
     2         - rho*tauM*uaNx(a)*updu(1,1,b)
            lK(1,a,b)  = lK(1,a,b)  + wl*(T2 + T1)

!           dRm_a1/du_b2
            T2 = mu*rM(2,1) + tauC*rM(1,2) + esNx(1,a)*mu_g*esNx(2,b)
     2         - rho*tauM*uaNx(a)*updu(2,1,b)
            lK(2,a,b)  = lK(2,a,b)  + wl*(T2)

!           dRm_a1/du_b3
            T2 = mu*rM(3,1) + tauC*rM(1,3) + esNx(1,a)*mu_g*esNx(3,b)
     2         - rho*tauM*uaNx(a)*updu(3,1,b)
            lK(3,a,b)  = lK(3,a,b)  + wl*(T2)

!           dRm_a2/du_b1
            T2 = mu*rM(1,2) + tauC*rM(2,1) + esNx(2,a)*mu_g*esNx(1,b)
     2         - rho*tauM*uaNx(a)*updu(1,2,b)
            lK(5,a,b)  = lK(5,a,b)  + wl*(T2)

!           dRm_a2/du_b2
            T2 = (mu + tauC)*rM(2,2) + esNx(2,a)*mu_g*esNx(2,b)
     2         - rho*tauM*uaNx(a)*updu(2,2,b)
            lK(6,a,b)  = lK(6,a,b)  + wl*(T2 + T1)

!           dRm_a2/du_b3
            T2 = mu*rM(3,2) + tauC*rM(2,3) + esNx(2,a)*mu_g*esNx(3,b)
     2         - rho*tauM*uaNx(a)*updu(3,2,b)
            lK(7,a,b)  = lK(7,a,b)  + wl*(T2)

!           dRm_a3/du_b1
            T2 = mu*rM(1,3) + tauC*rM(3,1) + esNx(3,a)*mu_g*esNx(1,b)
     2         - rho*tauM*uaNx(a)*updu(1,3,b)
            lK(9,a,b)  = lK(9,a,b)  + wl*(T2)

!           dRm_a3/du_b2
            T2 = mu*rM(2,3) + tauC*rM(3,2) + esNx(3,a)*mu_g*esNx(2,b)
     2         - rho*tauM*uaNx(a)*updu(2,3,b)
            lK(10,a,b) = lK(10,a,b) + wl*(T2)

!           dRm_a3/du_b3
            T2 = (mu + tauC)*rM(3,3) + esNx(3,a)*mu_g*esNx(3,b)
     2         - rho*tauM*uaNx(a)*updu(3,3,b)
            lK(11,a,b) = lK(11,a,b) + wl*(T2 + T1)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T1 = rho*tauM*uaNx(a)
!           dRm_a1/dp_b
            lK(4,a,b)  = lK(4,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1)
!           dRm_a2/dp_b
            lK(8,a,b)  = lK(8,a,b)  - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T1)
!           dRm_a3/dp_b
            lK(12,a,b) = lK(12,a,b) - wl*(Nwx(3,a)*Nq(b) - Nqx(3,b)*T1)
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID3D_M
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D_M(vmsFlag, eNoNw, eNoNq, w, Kxi, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Kxi(2,2), Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), Nqx(2,eNoNq), Nwxx(3,eNoNw), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) a, b, k
      REAL(KIND=RKIND) ctM, ctC, amd, wl, wr, rho, tauM, tauC, tauB, kT,
     2   kS, kU, divU, gam, mu, mu_s, mu_g, p, pa, u(2), ud(2), px(2),
     3   f(2), up(2), ua(2), ux(2,2), uxx(2,2,2), es(2,2), es_x(2,2,2),
     4   esNx(2,eNoNw), mu_x(2), rV(2), rS(2), rM(2,2), updu(2,2,eNoNw),
     5   uNx(eNoNw), upNx(eNoNw), uaNx(eNoNw), NxNx, d2u2(2), T1, T2

      ctM  = 1._RKIND
      ctC  = 36._RKIND

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1
      wr   = w*rho

!     Note that indices are not selected based on the equation because
!     fluid equation always come first
!     Velocity and its gradients, inertia (acceleration & body force)
      ud  = -f
      u   = 0._RKIND
      ux  = 0._RKIND
      uxx = 0._RKIND
      DO a=1, eNoNw
         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

         uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a)
         uxx(2,1,2) = uxx(2,1,2) + Nwxx(2,a)*yl(1,a)
         uxx(2,1,1) = uxx(2,1,1) + Nwxx(3,a)*yl(1,a)

         uxx(1,2,1) = uxx(1,2,1) + Nwxx(1,a)*yl(2,a)
         uxx(2,2,2) = uxx(2,2,2) + Nwxx(2,a)*yl(2,a)
         uxx(2,2,1) = uxx(2,2,1) + Nwxx(3,a)*yl(2,a)
      END DO
      divU = ux(1,1) + ux(2,2)

      uxx(1,1,2) = uxx(2,1,1)
      uxx(1,2,2) = uxx(2,2,1)

      d2u2(1) = uxx(1,1,1) + uxx(2,1,2)
      d2u2(2) = uxx(1,2,1) + uxx(2,2,2)

!     Pressure and its gradient
      p  = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p  = p + Nq(a)*yl(3,a)

         px(1) = px(1) + Nqx(1,a)*yl(3,a)
         px(2) = px(2) + Nqx(2,a)*yl(3,a)
      END DO

!     Update convection velocity relative to mesh velocity
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(4,a)
            u(2) = u(2) - Nw(a)*yl(5,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(2,1) = ux(2,1) + ux(1,2)
      es(1,2) = es(2,1)

      DO a=1, eNoNw
        esNx(1,a) = es(1,1)*Nwx(1,a) + es(2,1)*Nwx(2,a)
        esNx(2,a) = es(1,2)*Nwx(1,a) + es(2,2)*Nwx(2,a)
      END DO

      DO k=1, 2
         es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k)
         es_x(2,2,k) = uxx(2,2,k) + uxx(2,2,k)
         es_x(2,1,k) = uxx(2,1,k) + uxx(1,2,k)

         es_x(1,2,k) = es_x(2,1,k)
      END DO

      mu_x(1) = (es_x(1,1,1)*es(1,1) + es_x(2,2,1)*es(2,2))*0.5_RKIND
     2        +  es_x(2,1,1)*es(2,1)

      mu_x(2) = (es_x(1,1,2)*es(1,1) + es_x(2,2,2)*es(2,2))*0.5_RKIND
     2        +  es_x(2,1,2)*es(2,1)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_g := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

      IF (ISZERO(gam)) THEN
         mu_g = 0._RKIND
      ELSE
         mu_g = mu_g/gam
      END IF
      mu_x(:) = mu_g * mu_x(:)

!     Stabilization parameters
      kT = 4._RKIND*(ctM/dt)**2._RKIND

      kU = u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1)
     2   + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2)

      kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
     2   + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2)
      kS = ctC * kS * (mu/rho)**2._RKIND

      tauM = 1._RKIND / (rho * SQRT( kT + kU + kS ))

      rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1)
      rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2)

      rS(1) = mu_x(1)*es(1,1) + mu_x(2)*es(2,1) + mu*d2u2(1)
      rS(2) = mu_x(1)*es(1,2) + mu_x(2)*es(2,2) + mu*d2u2(2)

      up(1) = -tauM*(rho*rV(1) + px(1) - rS(1))
      up(2) = -tauM*(rho*rV(2) + px(2) - rS(2))

      IF (vmsFlag) THEN
         tauC = 1._RKIND / (tauM * (Kxi(1,1) + Kxi(2,2)))
         tauB = up(1)*up(1)*Kxi(1,1) + up(2)*up(1)*Kxi(2,1)
     2        + up(1)*up(2)*Kxi(1,2) + up(2)*up(2)*Kxi(2,2)
         IF (ISZERO(tauB)) tauB = eps
         tauB = rho/SQRT(tauB)

         ua(1) = u(1) + up(1)
         ua(2) = u(2) + up(2)
         pa    = p - tauC*divU
      ELSE
         tauC  = 0._RKIND
         tauB  = 0._RKIND
         ua(:) = u(:)
         pa    = p
      END IF

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))

      rM(1,1) = mu*es(1,1) - rho*up(1)*ua(1) + rV(1)*up(1) - pa
      rM(2,1) = mu*es(2,1) - rho*up(1)*ua(2) + rV(1)*up(2)

      rM(1,2) = mu*es(1,2) - rho*up(2)*ua(1) + rV(2)*up(1)
      rM(2,2) = mu*es(2,2) - rho*up(2)*ua(2) + rV(2)*up(2) - pa

      rV(1) = ud(1) + ua(1)*ux(1,1) + ua(2)*ux(2,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*ux(2,2)

!     Local residue
      DO a=1, eNoNw
         lR(1,a) = lR(1,a) + wr*Nw(a)*rV(1) + w*(Nwx(1,a)*rM(1,1)
     2      + Nwx(2,a)*rM(2,1))

         lR(2,a) = lR(2,a) + wr*Nw(a)*rV(2) + w*(Nwx(1,a)*rM(1,2)
     2      + Nwx(2,a)*rM(2,2))

!        Quantities used for stiffness matrix
         uNx(a)  = u(1)*Nwx(1,a)  + u(2)*Nwx(2,a)
         upNx(a) = up(1)*Nwx(1,a) + up(2)*Nwx(2,a)
         IF (vmsFlag) THEN
            uaNx(a) = uNx(a) + upNx(a)
         ELSE
            uaNx(a) = uNx(a)
         END IF

         T1 = -rho*uNx(a) + mu*(Nwxx(1,a) + Nwxx(2,a))
     2      + mu_x(1)*Nwx(1,a) + mu_x(2)*Nwx(2,a)

         updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1
         updu(2,1,a) = mu_x(2)*Nwx(1,a) + d2u2(1)*mu_g*esNx(2,a)

         updu(1,2,a) = mu_x(1)*Nwx(2,a) + d2u2(2)*mu_g*esNx(1,a)
         updu(2,2,a) = mu_x(2)*Nwx(2,a) + d2u2(2)*mu_g*esNx(2,a) + T1
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            rM(1,1) = Nwx(1,a)*Nwx(1,b)
            rM(2,1) = Nwx(2,a)*Nwx(1,b)
            rM(1,2) = Nwx(1,a)*Nwx(2,b)
            rM(2,2) = Nwx(2,a)*Nwx(2,b)

            NxNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b)

            T1 = mu*NxNx + rho*amd*Nw(b)*(Nw(a) + rho*tauM*uaNx(a))
     2         + rho*Nw(a)*(uNx(b)+upNx(b)) + tauB*upNx(a)*upNx(b)

!           dRm_a1/du_b1
            T2 = (mu + tauC)*rM(1,1) + esNx(1,a)*mu_g*esNx(1,b)
     2         - rho*tauM*uaNx(a)*updu(1,1,b)
            lK(1,a,b)  = lK(1,a,b)  + wl*(T2 + T1)

!           dRm_a1/du_b2
            T2 = mu*rM(2,1) + tauC*rM(1,2) + esNx(1,a)*mu_g*esNx(2,b)
     2         - rho*tauM*uaNx(a)*updu(2,1,b)
            lK(2,a,b)  = lK(2,a,b)  + wl*(T2)

!           dRm_a2/du_b1
            T2 = mu*rM(1,2) + tauC*rM(2,1) + esNx(2,a)*mu_g*esNx(1,b)
     2         - rho*tauM*uaNx(a)*updu(1,2,b)
            lK(4,a,b)  = lK(4,a,b)  + wl*(T2)

!           dRm_a2/du_b2
            T2 = (mu + tauC)*rM(2,2) + esNx(2,a)*mu_g*esNx(2,b)
     2         - rho*tauM*uaNx(a)*updu(2,2,b)
            lK(5,a,b)  = lK(5,a,b)  + wl*(T2 + T1)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T1 = rho*tauM*uaNx(a)
!           dRm_a1/dp_b
            lK(3,a,b)  = lK(3,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T1)
!           dRm_a2/dp_b
            lK(6,a,b)  = lK(6,a,b)  - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T1)
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID2D_M
!####################################################################
      SUBROUTINE FLUID3D_C(vmsFlag, eNoNw, eNoNq, w, Kxi, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Kxi(3,3), Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), Nqx(3,eNoNq), Nwxx(6,eNoNw), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) a, b, k
      REAL(KIND=RKIND) ctM, ctC, amd, wl, rho, tauM, kT, kS, kU, divU,
     2   gam, mu, mu_s, mu_g, u(3), ud(3), px(3), f(3), up(3), ux(3,3),
     3   uxx(3,3,3), es(3,3), es_x(3,3,3), esNx(3,eNoNw), mu_x(3), uNx,
     4   rV(3), rS(3), updu(3,3,eNoNw), d2u2(3), upNx, NxNx, T1, T2

      ctM  = 1._RKIND
      ctC  = 36._RKIND

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

!     Note that indices are not selected based on the equation because
!     fluid equation always come first
!     Velocity and its gradients, inertia (acceleration & body force)
      ud  = -f
      u   = 0._RKIND
      ux  = 0._RKIND
      uxx = 0._RKIND
      DO a=1, eNoNw
         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
         ud(3)   = ud(3) + Nw(a)*(al(3,a)-bfl(3,a))

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)
         u(3)    = u(3) + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nwx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nwx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nwx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nwx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nwx(3,a)*yl(3,a)

         uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a)
         uxx(2,1,2) = uxx(2,1,2) + Nwxx(2,a)*yl(1,a)
         uxx(3,1,3) = uxx(3,1,3) + Nwxx(3,a)*yl(1,a)
         uxx(2,1,1) = uxx(2,1,1) + Nwxx(4,a)*yl(1,a)
         uxx(3,1,2) = uxx(3,1,2) + Nwxx(5,a)*yl(1,a)
         uxx(1,1,3) = uxx(1,1,3) + Nwxx(6,a)*yl(1,a)

         uxx(1,2,1) = uxx(1,2,1) + Nwxx(1,a)*yl(2,a)
         uxx(2,2,2) = uxx(2,2,2) + Nwxx(2,a)*yl(2,a)
         uxx(3,2,3) = uxx(3,2,3) + Nwxx(3,a)*yl(2,a)
         uxx(2,2,1) = uxx(2,2,1) + Nwxx(4,a)*yl(2,a)
         uxx(3,2,2) = uxx(3,2,2) + Nwxx(5,a)*yl(2,a)
         uxx(1,2,3) = uxx(1,2,3) + Nwxx(6,a)*yl(2,a)

         uxx(1,3,1) = uxx(1,3,1) + Nwxx(1,a)*yl(3,a)
         uxx(2,3,2) = uxx(2,3,2) + Nwxx(2,a)*yl(3,a)
         uxx(3,3,3) = uxx(3,3,3) + Nwxx(3,a)*yl(3,a)
         uxx(2,3,1) = uxx(2,3,1) + Nwxx(4,a)*yl(3,a)
         uxx(3,3,2) = uxx(3,3,2) + Nwxx(5,a)*yl(3,a)
         uxx(1,3,3) = uxx(1,3,3) + Nwxx(6,a)*yl(3,a)
      END DO
      divU = ux(1,1) + ux(2,2) + ux(3,3)

      uxx(1,1,2) = uxx(2,1,1)
      uxx(2,1,3) = uxx(3,1,2)
      uxx(3,1,1) = uxx(1,1,3)

      uxx(1,2,2) = uxx(2,2,1)
      uxx(2,2,3) = uxx(3,2,2)
      uxx(3,2,1) = uxx(1,2,3)

      uxx(1,3,2) = uxx(2,3,1)
      uxx(2,3,3) = uxx(3,3,2)
      uxx(3,3,1) = uxx(1,3,3)

      d2u2(1) = uxx(1,1,1) + uxx(2,1,2) + uxx(3,1,3)
      d2u2(2) = uxx(1,2,1) + uxx(2,2,2) + uxx(3,2,3)
      d2u2(3) = uxx(1,3,1) + uxx(2,3,2) + uxx(3,3,3)

!     Pressure and its gradient
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(4,a)
         px(2) = px(2) + Nqx(2,a)*yl(4,a)
         px(3) = px(3) + Nqx(3,a)*yl(4,a)
      END DO

!     Update convection velocity relative to mesh velocity
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(5,a)
            u(2) = u(2) - Nw(a)*yl(6,a)
            u(3) = u(3) - Nw(a)*yl(7,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)
      es(2,1) = ux(2,1) + ux(1,2)
      es(3,2) = ux(3,2) + ux(2,3)
      es(1,3) = ux(1,3) + ux(3,1)
      es(1,2) = es(2,1)
      es(2,3) = es(3,2)
      es(3,1) = es(1,3)

      DO a=1, eNoNw
        esNx(1,a) = es(1,1)*Nwx(1,a) + es(2,1)*Nwx(2,a)
     2            + es(3,1)*Nwx(3,a)
        esNx(2,a) = es(1,2)*Nwx(1,a) + es(2,2)*Nwx(2,a)
     2            + es(3,2)*Nwx(3,a)
        esNx(3,a) = es(1,3)*Nwx(1,a) + es(2,3)*Nwx(2,a)
     2            + es(3,3)*Nwx(3,a)
      END DO

      DO k=1, 3
         es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k)
         es_x(2,2,k) = uxx(2,2,k) + uxx(2,2,k)
         es_x(3,3,k) = uxx(3,3,k) + uxx(3,3,k)
         es_x(2,1,k) = uxx(2,1,k) + uxx(1,2,k)
         es_x(3,2,k) = uxx(3,2,k) + uxx(2,3,k)
         es_x(1,3,k) = uxx(1,3,k) + uxx(3,1,k)

         es_x(1,2,k) = es_x(2,1,k)
         es_x(2,3,k) = es_x(3,2,k)
         es_x(3,1,k) = es_x(1,3,k)
      END DO

      mu_x(1) = (es_x(1,1,1)*es(1,1) + es_x(2,2,1)*es(2,2)
     2        +  es_x(3,3,1)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,1)*es(2,1) + es_x(3,2,1)*es(3,2)
     4        +  es_x(1,3,1)*es(1,3)

      mu_x(2) = (es_x(1,1,2)*es(1,1) + es_x(2,2,2)*es(2,2)
     2        +  es_x(3,3,2)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,2)*es(2,1) + es_x(3,2,2)*es(3,2)
     4        +  es_x(1,3,2)*es(1,3)

      mu_x(3) = (es_x(1,1,3)*es(1,1) + es_x(2,2,3)*es(2,2)
     2        +  es_x(3,3,3)*es(3,3))*0.5_RKIND
     3        +  es_x(2,1,3)*es(2,1) + es_x(3,2,3)*es(3,2)
     4        +  es_x(1,3,3)*es(1,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1) + es(3,1)*es(3,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2) + es(3,2)*es(3,2)
     3    + es(1,3)*es(1,3) + es(2,3)*es(2,3) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_x := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

      IF (ISZERO(gam)) THEN
         mu_g = 0._RKIND
      ELSE
         mu_g = mu_g/gam
      END IF
      mu_x(:) = mu_g * mu_x(:)

      IF (vmsFlag) THEN
!        Stabilization parameters
         kT = 4._RKIND*(ctM/dt)**2._RKIND

         kU = u(1)*u(1)*Kxi(1,1) +u(2)*u(1)*Kxi(2,1) +u(3)*u(1)*Kxi(3,1)
     2      + u(1)*u(2)*Kxi(1,2) +u(2)*u(2)*Kxi(2,2) +u(3)*u(2)*Kxi(3,2)
     3      + u(1)*u(3)*Kxi(1,3) +u(2)*u(3)*Kxi(2,3) +u(3)*u(3)*Kxi(3,3)

         kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1) + Kxi(3,1)*Kxi(3,1)
     2      + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2) + Kxi(3,2)*Kxi(3,2)
     3      + Kxi(1,3)*Kxi(1,3) + Kxi(2,3)*Kxi(2,3) + Kxi(3,3)*Kxi(3,3)
         kS = ctC * kS * (mu/rho)**2._RKIND

         tauM = 1._RKIND / (rho * SQRT( kT + kU + kS ))

         rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + u(3)*ux(3,1)
         rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + u(3)*ux(3,2)
         rV(3) = ud(3) + u(1)*ux(1,3) + u(2)*ux(2,3) + u(3)*ux(3,3)

         rS(1) = mu_x(1)*es(1,1) + mu_x(2)*es(2,1) + mu_x(3)*es(3,1)
     2         + mu*d2u2(1)
         rS(2) = mu_x(1)*es(1,2) + mu_x(2)*es(2,2) + mu_x(3)*es(3,2)
     2         + mu*d2u2(2)
         rS(3) = mu_x(1)*es(1,3) + mu_x(2)*es(2,3) + mu_x(3)*es(3,3)
     2         + mu*d2u2(3)

         up(1) = -tauM*(rho*rV(1) + px(1) - rS(1))
         up(2) = -tauM*(rho*rV(2) + px(2) - rS(2))
         up(3) = -tauM*(rho*rV(3) + px(3) - rS(3))

         DO a=1, eNoNw
            uNx = u(1)*Nwx(1,a) + u(2)*Nwx(2,a) + u(3)*Nwx(3,a)

            T1  = -rho*uNx + mu*(Nwxx(1,a) + Nwxx(2,a) + Nwxx(3,a))
     2          + mu_x(1)*Nwx(1,a) + mu_x(2)*Nwx(2,a) + mu_x(3)*Nwx(3,a)

            updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1
            updu(2,1,a) = mu_x(2)*Nwx(1,a) + d2u2(1)*mu_g*esNx(2,a)
            updu(3,1,a) = mu_x(3)*Nwx(1,a) + d2u2(1)*mu_g*esNx(3,a)

            updu(1,2,a) = mu_x(1)*Nwx(2,a) + d2u2(2)*mu_g*esNx(1,a)
            updu(2,2,a) = mu_x(2)*Nwx(2,a) + d2u2(2)*mu_g*esNx(2,a) + T1
            updu(3,2,a) = mu_x(3)*Nwx(2,a) + d2u2(2)*mu_g*esNx(3,a)

            updu(1,3,a) = mu_x(1)*Nwx(3,a) + d2u2(3)*mu_g*esNx(1,a)
            updu(2,3,a) = mu_x(2)*Nwx(3,a) + d2u2(3)*mu_g*esNx(2,a)
            updu(3,3,a) = mu_x(3)*Nwx(3,a) + d2u2(3)*mu_g*esNx(3,a) + T1
         END DO
      ELSE
         tauM  = 0._RKIND
         up(:) = 0._RKIND
         updu  = 0._RKIND
      END IF

!     Local residue
      DO a=1, eNoNq
         upNx    = up(1)*Nqx(1,a) + up(2)*Nqx(2,a) + up(3)*Nqx(3,a)
         lR(4,a) = lR(4,a) + w*(Nq(a)*divU - upNx)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         T1 = rho*amd*Nw(b)
         DO a=1, eNoNq
!           dRc_a/dU_b1
            T2 = Nqx(1,a)*(updu(1,1,b) - T1) + Nqx(2,a)*updu(1,2,b)
     2         + Nqx(3,a)*updu(1,3,b)
            lK(13,a,b) = lK(13,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2)

!           dRc_a/dU_b2
            T2 = Nqx(1,a)*updu(2,1,b) + Nqx(2,a)*(updu(2,2,b) - T1)
     2         + Nqx(3,a)*updu(2,3,b)
            lK(14,a,b) = lK(14,a,b) + wl*(Nq(a)*Nwx(2,b) - tauM*T2)

!           dRc_a/dU_b3
            T2 = Nqx(1,a)*updu(3,1,b) + Nqx(2,a)*updu(3,2,b)
     2         + Nqx(3,a)*(updu(3,3,b) - T1)
            lK(15,a,b) = lK(15,a,b) + wl*(Nq(a)*Nwx(3,b) - tauM*T2)
         END DO
      END DO

      IF (vmsFlag) THEN
         DO b=1, eNoNq
            DO a=1, eNoNq
!              dC/dP
               NxNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b)
     2              + Nqx(3,a)*Nqx(3,b)
               lK(16,a,b) = lK(16,a,b) + wl*tauM*NxNx
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE FLUID3D_C
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D_C(vmsFlag, eNoNw, eNoNq, w, Kxi, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: vmsFlag
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Kxi(2,2), Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), Nqx(2,eNoNq), Nwxx(3,eNoNw), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) a, b, k
      REAL(KIND=RKIND) ctM, ctC, amd, wl, rho, tauM, kT, kS, kU, divU,
     2   gam, mu, mu_s, mu_g, u(2), ud(2), px(2), f(2), up(2), ux(2,2),
     3   uxx(2,2,2), es(2,2), es_x(2,2,2), esNx(2,eNoNw), mu_x(2), uNx,
     4   rV(2), rS(2), updu(2,2,eNoNw), d2u2(2), upNx, NxNx, T1, T2

      ctM  = 1._RKIND
      ctC  = 36._RKIND

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

!     Note that indices are not selected based on the equation because
!     fluid equation always come first
!     Velocity and its gradients, inertia (acceleration & body force)
      ud  = -f
      u   = 0._RKIND
      ux  = 0._RKIND
      uxx = 0._RKIND
      DO a=1, eNoNw
         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

         uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a)
         uxx(2,1,2) = uxx(2,1,2) + Nwxx(2,a)*yl(1,a)
         uxx(2,1,1) = uxx(2,1,1) + Nwxx(3,a)*yl(1,a)

         uxx(1,2,1) = uxx(1,2,1) + Nwxx(1,a)*yl(2,a)
         uxx(2,2,2) = uxx(2,2,2) + Nwxx(2,a)*yl(2,a)
         uxx(2,2,1) = uxx(2,2,1) + Nwxx(3,a)*yl(2,a)
      END DO
      divU = ux(1,1) + ux(2,2)

      uxx(1,1,2) = uxx(2,1,1)
      uxx(1,2,2) = uxx(2,2,1)

      d2u2(1) = uxx(1,1,1) + uxx(2,1,2)
      d2u2(2) = uxx(1,2,1) + uxx(2,2,2)

!     Pressure and its gradient
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(3,a)
         px(2) = px(2) + Nqx(2,a)*yl(3,a)
      END DO

!     Update convection velocity relative to mesh velocity
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(4,a)
            u(2) = u(2) - Nw(a)*yl(5,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,1) = ux(2,1) + ux(1,2)
      es(2,2) = ux(2,2) + ux(2,2)
      es(1,2) = es(2,1)

      DO a=1, eNoNw
        esNx(1,a) = es(1,1)*Nwx(1,a) + es(2,1)*Nwx(2,a)
        esNx(2,a) = es(1,2)*Nwx(1,a) + es(2,2)*Nwx(2,a)
      END DO

      DO k=1, 2
         es_x(1,1,k) = uxx(1,1,k) + uxx(1,1,k)
         es_x(2,2,k) = uxx(2,2,k) + uxx(2,2,k)
         es_x(2,1,k) = uxx(2,1,k) + uxx(1,2,k)

         es_x(1,2,k) = es_x(2,1,k)
      END DO

      mu_x(1) = (es_x(1,1,1)*es(1,1) + es_x(2,2,1)*es(2,2))*0.5_RKIND
     3        +  es_x(2,1,1)*es(2,1)

      mu_x(2) = (es_x(1,1,2)*es(1,1) + es_x(2,2,2)*es(2,2))*0.5_RKIND
     3        +  es_x(2,1,2)*es(2,1)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_x := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

      IF (ISZERO(gam)) THEN
         mu_g = 0._RKIND
      ELSE
         mu_g = mu_g/gam
      END IF
      mu_x(:) = mu_g * mu_x(:)

      IF (vmsFlag) THEN
!        Stabilization parameters
         kT = 4._RKIND*(ctM/dt)**2_RKIND

         kU = u(1)*u(1)*Kxi(1,1) + u(2)*u(1)*Kxi(2,1)
     2      + u(1)*u(2)*Kxi(1,2) + u(2)*u(2)*Kxi(2,2)

         kS = Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
     2      + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2)
         kS = ctC * kS * (mu/rho)**2._RKIND

         tauM = 1._RKIND / (rho * SQRT( kT + kU + kS ))

         rV(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1)
         rV(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2)

         rS(1) = mu_x(1)*es(1,1) + mu_x(2)*es(2,1) + mu*d2u2(1)
         rS(2) = mu_x(1)*es(1,2) + mu_x(2)*es(2,2) + mu*d2u2(2)

         up(1) = -tauM*(rho*rV(1) + px(1) - rS(1))
         up(2) = -tauM*(rho*rV(2) + px(2) - rS(2))

         DO a=1, eNoNw
            uNx = u(1)*Nwx(1,a) + u(2)*Nwx(2,a)

            T1  = -rho*uNx + mu*(Nwxx(1,a) + Nwxx(2,a))
     2          + mu_x(1)*Nwx(1,a) + mu_x(2)*Nwx(2,a)

            updu(1,1,a) = mu_x(1)*Nwx(1,a) + d2u2(1)*mu_g*esNx(1,a) + T1
            updu(2,1,a) = mu_x(2)*Nwx(1,a) + d2u2(1)*mu_g*esNx(2,a)

            updu(1,2,a) = mu_x(1)*Nwx(2,a) + d2u2(2)*mu_g*esNx(1,a)
            updu(2,2,a) = mu_x(2)*Nwx(2,a) + d2u2(2)*mu_g*esNx(2,a) + T1
         END DO
      ELSE
         tauM  = 0._RKIND
         up(:) = 0._RKIND
         updu  = 0._RKIND
      END IF

!     Local residue
      DO a=1, eNoNq
         upNx    = up(1)*Nqx(1,a) + up(2)*Nqx(2,a)
         lR(3,a) = lR(3,a) + w*(Nq(a)*divU - upNx)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         T1 = rho*amd*Nw(b)
         DO a=1, eNoNq
!           dRc_a/dU_b1
            T2 = Nqx(1,a)*(updu(1,1,b) - T1) + Nqx(2,a)*updu(1,2,b)
            lK(7,a,b) = lK(7,a,b) + wl*(Nq(a)*Nwx(1,b) - tauM*T2)

!           dRc_a/dU_b2
            T2 = Nqx(1,a)*updu(2,1,b) + Nqx(2,a)*(updu(2,2,b) - T1)
            lK(8,a,b) = lK(8,a,b) + wl*(Nq(a)*Nwx(2,b) - tauM*T2)
         END DO
      END DO

      IF (vmsFlag) THEN
         DO b=1, eNoNq
            DO a=1, eNoNq
!              dC/dP
               NxNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b)
               lK(9,a,b) = lK(9,a,b) + wl*tauM*NxNx
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE FLUID2D_C
!####################################################################
      PURE SUBROUTINE BFLUID(eNoN, w, N, y, h, nV, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), y(tDof), h, nV(nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i, j
      REAL(KIND=RKIND) T1, wl, hc(nsd), udn, u(nsd)

      wl  = w*eq(cEq)%af*eq(cEq)%gam*dt
      udn = 0._RKIND
      IF (mvMsh) THEN
         DO i=1, nsd
            j    = i + nsd + 1
            u(i) = y(i) - y(j)
            udn  = udn + u(i)*nV(i)
         END DO
      ELSE
         DO i=1, nsd
            u(i) = y(i)
            udn  = udn + u(i)*nV(i)
         END DO
      END IF

      udn = 0.5_RKIND*eq(cEq)%dmn(cDmn)%prop(backflow_stab)*
     2   eq(cEq)%dmn(cDmn)%prop(fluid_density)*(udn - ABS(udn))
      hc  = h*nV + udn*u

!     Here the loop is started for constructing left and right hand side
      IF (nsd .EQ. 2) THEN
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            DO b=1, eNoN
               T1        = wl*N(a)*N(b)*udn
               lK(1,a,b) = lK(1,a,b) - T1
               lK(5,a,b) = lK(5,a,b) - T1
            END DO
         END DO
      ELSE
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - w*N(a)*hc(1)
            lR(2,a) = lR(2,a) - w*N(a)*hc(2)
            lR(3,a) = lR(3,a) - w*N(a)*hc(3)
            DO b=1, eNoN
               T1 = wl*N(a)*N(b)*udn
               lK(1,a,b)  = lK(1,a,b)  - T1
               lK(6,a,b)  = lK(6,a,b)  - T1
               lK(11,a,b) = lK(11,a,b) - T1
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE BFLUID
!####################################################################
      SUBROUTINE BWFLUID3D(eNoNw, eNoNq, w, Nw, Nq, Nwx, yl, ub, nV,
     2   tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), yl(tDof,eNoNw), ub(3), nV(3), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, tauT, tauN, wl, p, un, ubn, gam, mu,
     2   mu_s, mu_g, u(3), uh(3), ux(3,3), es(3,3), sgmn(3), NxN(eNoNw),
     3   T1, T2, T3, rV(3), rM(3,3)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      tauT = tauB(1)
      tauN = tauB(2)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wl   = w * T1

      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoNw
         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)
         u(3)    = u(3) + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nwx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nwx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nwx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nwx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nwx(3,a)*yl(3,a)

         Nxn(a)  = Nwx(1,a)*nV(1) + Nwx(2,a)*nV(2) + Nwx(3,a)*nV(3)
      END DO

      p = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(4,a)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoNw
            uh(1) = uh(1) + Nw(a)*yl(5,a)
            uh(2) = uh(2) + Nw(a)*yl(6,a)
            uh(3) = uh(3) + Nw(a)*yl(7,a)
         END DO
      END IF
      un = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2) + (u(3)-uh(3))*nV(3)
      un = (ABS(un) - un) * 0.5_RKIND

      u(:) = u(:) - ub(:)
      ubn  = u(1)*nV(1) + u(2)*nV(2) + u(3)*nV(3)

!     Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)
      es(2,1) = ux(2,1) + ux(1,2)
      es(3,2) = ux(3,2) + ux(2,3)
      es(1,3) = ux(1,3) + ux(3,1)
      es(1,2) = es(2,1)
      es(2,3) = es(3,2)
      es(3,1) = es(1,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1) + es(3,1)*es(3,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2) + es(3,2)*es(3,2)
     3    + es(1,3)*es(1,3) + es(2,3)*es(2,3) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_g := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

!     sigma.n (deviatoric)
      sgmn(1) = mu*(es(1,1)*nV(1) + es(2,1)*nV(2) + es(3,1)*nV(3))
      sgmn(2) = mu*(es(1,2)*nV(1) + es(2,2)*nV(2) + es(3,2)*nV(3))
      sgmn(3) = mu*(es(1,3)*nV(1) + es(2,3)*nV(2) + es(3,3)*nV(3))

      rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1)
     2      + (tauN-tauT)*ubn*nV(1)
      rV(2) = p*nV(2) - sgmn(2) + (tauT + rho*un)*u(2)
     2      + (tauN-tauT)*ubn*nV(2)
      rV(3) = p*nV(3) - sgmn(3) + (tauT + rho*un)*u(3)
     2      + (tauN-tauT)*ubn*nV(3)

      rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) )
      rM(2,1) = -mu*( u(1)*nV(2) + u(2)*nV(1) )
      rM(3,1) = -mu*( u(1)*nV(3) + u(3)*nV(1) )

      rM(1,2) = -mu*( u(2)*nV(1) + u(1)*nV(2) )
      rM(2,2) = -mu*( u(2)*nV(2) + u(2)*nV(2) )
      rM(3,2) = -mu*( u(2)*nV(3) + u(3)*nV(2) )

      rM(1,3) = -mu*( u(3)*nV(1) + u(1)*nV(3) )
      rM(2,3) = -mu*( u(3)*nV(2) + u(2)*nV(3) )
      rM(3,3) = -mu*( u(3)*nV(3) + u(3)*nV(3) )

!     Local residue (momentum)
      DO a=1, eNoNw
         lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(1,a)*rM(1,1)
     2      + Nwx(2,a)*rM(2,1) + Nwx(3,a)*rM(3,1))

         lR(2,a) = lR(2,a) + w*(Nw(a)*rV(2) + Nwx(1,a)*rM(1,2)
     2      + Nwx(2,a)*rM(2,2) + Nwx(3,a)*rM(3,2))

         lR(3,a) = lR(3,a) + w*(Nw(a)*rV(3) + Nwx(1,a)*rM(1,3)
     2      + Nwx(2,a)*rM(2,3) + Nwx(3,a)*rM(3,3))
      END DO

!     Local residue (continuity)
      DO a=1, eNoNq
         lR(4,a) = lR(4,a) - w*Nq(a)*ubn
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            T3 = Nw(a)*Nw(b)
            T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b))

!           dRm_a1/du_b1
            T2 = (tauN - tauT)*T3*nV(1)*nV(1)
     2         - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1))
            lK(1,a,b)  = lK(1,a,b)  + wl*(T2 + T1)

!           dRm_a1/du_b2
            T2 = (tauN - tauT)*T3*nV(1)*nV(2)
     2         - mu*(Nw(a)*Nwx(1,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(1))
            lK(2,a,b)  = lK(2,a,b)  + wl*T2

!           dRm_a1/du_b3
            T2 = (tauN - tauT)*T3*nV(1)*nV(3)
     2         - mu*(Nw(a)*Nwx(1,b)*nV(3) + Nw(b)*Nwx(3,a)*nV(1))
            lK(3,a,b)  = lK(3,a,b)  + wl*T2

!           dRm_a2/du_b1
            T2 = (tauN - tauT)*T3*nV(2)*nV(1)
     2         - mu*(Nw(a)*Nwx(2,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(2))
            lK(5,a,b)  = lK(5,a,b)  + wl*T2

!           dRm_a2/du_b2
            T2 = (tauN - tauT)*T3*nV(2)*nV(2)
     2         - mu*(Nw(a)*Nwx(2,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(2))
            lK(6,a,b)  = lK(6,a,b)  + wl*(T2 + T1)

!           dRm_a2/du_b3
            T2 = (tauN - tauT)*T3*nV(2)*nV(3)
     2         - mu*(Nw(a)*Nwx(2,b)*nV(3) + Nw(b)*Nwx(3,a)*nV(2))
            lK(7,a,b)  = lK(7,a,b)  + wl*T2

!           dRm_a3/du_b1
            T2 = (tauN - tauT)*T3*nV(3)*nV(1)
     2         - mu*(Nw(a)*Nwx(3,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(3))
            lK(9,a,b)  = lK(9,a,b)  + wl*T2

!           dRm_a3/du_b2
            T2 = (tauN - tauT)*T3*nV(3)*nV(2)
     2         - mu*(Nw(a)*Nwx(3,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(3))
            lK(10,a,b) = lK(10,a,b) + wl*T2

!           dRm_a3/du_b3
            T2 = (tauN - tauT)*T3*nV(3)*nV(3)
     2         - mu*(Nw(a)*Nwx(3,b)*nV(3) + Nw(b)*Nwx(3,a)*nV(3))
            lK(11,a,b) = lK(11,a,b) + wl*(T2 + T1)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T1 = wl*Nw(a)*Nq(b)
!           dRm_a1/dp_b
            lK(4,a,b)  = lK(4,a,b)  + T1*nV(1)
!           dRm_a2/dp_b
            lK(8,a,b)  = lK(8,a,b)  + T1*nV(2)
!           dRm_a3/dp_b
            lK(12,a,b) = lK(12,a,b) + T1*nV(3)
         END DO
      END DO

      DO b=1, eNoNw
         DO a=1, eNoNq
            T1 = wl*Nq(a)*Nw(b)
!           dRc_a/du_b1
            lK(13,a,b) = lK(13,a,b) - T1*nV(1)
!           dRc_a/du_b2
            lK(14,a,b) = lK(14,a,b) - T1*nV(2)
!           dRc_a/du_b3
            lK(15,a,b) = lK(15,a,b) - T1*nV(3)
         END DO
      END DO

      RETURN
      END SUBROUTINE BWFLUID3D
!--------------------------------------------------------------------
      SUBROUTINE BWFLUID2D(eNoNw, eNoNq, w, Nw, Nq, Nwx, yl, ub, nV,
     2   tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), yl(tDof,eNoNw), ub(2), nV(2), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, tauT, tauN, wl, p, un, ubn, gam, mu,
     2   mu_s, mu_g, u(2), uh(2), ux(2,2), es(2,2), sgmn(2), NxN(eNoNw),
     3   T1, T2, T3, rV(2), rM(2,2)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      tauT = tauB(1)
      tauN = tauB(2)

      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wl   = w * T1

      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoNw
         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

         Nxn(a)  = Nwx(1,a)*nV(1) + Nwx(2,a)*nV(2)
      END DO

      p = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(3,a)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoNw
            uh(1) = uh(1) + Nw(a)*yl(4,a)
            uh(2) = uh(2) + Nw(a)*yl(5,a)
         END DO
      END IF
      un = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2)
      un = (ABS(un) - un) * 0.5_RKIND

      u(:) = u(:) - ub(:)
      ubn  = u(1)*nV(1) + u(2)*nV(2)

!     Strain rate tensor 2*e_ij := (u_i,j + u_j,i)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(2,1) = ux(2,1) + ux(1,2)
      es(1,2) = es(2,1)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(2,1)*es(2,1)
     2    + es(1,2)*es(1,2) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
!     The returned mu_g := (d\mu / d\gamma)
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, mu, mu_s, mu_g)

!     sigma.n (deviatoric)
      sgmn(1) = mu*(es(1,1)*nV(1) + es(2,1)*nV(2))
      sgmn(2) = mu*(es(1,2)*nV(1) + es(2,2)*nV(2))

      rV(1) = p*nV(1) - sgmn(1) + (tauT + rho*un)*u(1)
     2      + (tauN-tauT)*ubn*nV(1)
      rV(2) = p*nV(2) - sgmn(2) + (tauT + rho*un)*u(2)
     2      + (tauN-tauT)*ubn*nV(2)

      rM(1,1) = -mu*( u(1)*nV(1) + u(1)*nV(1) )
      rM(2,1) = -mu*( u(1)*nV(2) + u(2)*nV(1) )
      rM(1,2) = -mu*( u(2)*nV(1) + u(1)*nV(2) )
      rM(2,2) = -mu*( u(2)*nV(2) + u(2)*nV(2) )

!     Local residue (momentum)
      DO a=1, eNoNw
         lR(1,a) = lR(1,a) + w*(Nw(a)*rV(1) + Nwx(1,a)*rM(1,1)
     2      + Nwx(2,a)*rM(2,1))

         lR(2,a) = lR(2,a) + w*(Nw(a)*rV(2) + Nwx(1,a)*rM(1,2)
     2      + Nwx(2,a)*rM(2,2))
      END DO

!     Local residue (continuity)
      DO a=1, eNoNq
         lR(3,a) = lR(3,a) - w*Nq(a)*ubn
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            T3 = Nw(a)*Nw(b)
            T1 = (tauT + rho*un)*T3 - mu*(Nw(a)*Nxn(b) + Nxn(a)*Nw(b))

!           dRm_a1/du_b1
            T2 = (tauN - tauT)*T3*nV(1)*nV(1)
     2         - mu*(Nw(a)*Nwx(1,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(1))
            lK(1,a,b) = lK(1,a,b) + wl*(T2 + T1)

!           dRm_a1/du_b2
            T2 = (tauN - tauT)*T3*nV(1)*nV(2)
     2         - mu*(Nw(a)*Nwx(1,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(1))
            lK(2,a,b) = lK(2,a,b) + wl*T2

!           dRm_a2/du_b1
            T2 = (tauN - tauT)*T3*nV(2)*nV(1)
     2         - mu*(Nw(a)*Nwx(2,b)*nV(1) + Nw(b)*Nwx(1,a)*nV(2))
            lK(4,a,b) = lK(4,a,b) + wl*T2

!           dRm_a2/du_b2
            T2 = (tauN - tauT)*T3*nV(2)*nV(2)
     2         - mu*(Nw(a)*Nwx(2,b)*nV(2) + Nw(b)*Nwx(2,a)*nV(2))
            lK(5,a,b) = lK(5,a,b)  + wl*(T2 + T1)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T1 = wl*Nw(a)*Nq(b)
!           dRm_a1/dp_b
            lK(3,a,b) = lK(3,a,b) + T1*nV(1)
!           dRm_a2/dp_b
            lK(6,a,b) = lK(6,a,b) + T1*nV(2)
         END DO
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNq
            T1 = wl*Nq(a)*Nw(b)
!           dRc_a/du_b1
            lK(7,a,b) = lK(7,a,b) - T1*nV(1)
!           dRc_a/du_b2
            lK(8,a,b) = lK(8,a,b) - T1*nV(2)
         END DO
      END DO

      RETURN
      END SUBROUTINE BWFLUID2D
!####################################################################
      SUBROUTINE GETVISCOSITY(lDmn, gamma, mu, mu_s, mu_x)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(INOUT)  :: gamma
      REAL(KIND=RKIND), INTENT(OUT) :: mu, mu_s, mu_x

      REAL(KIND=RKIND) :: mu_i, mu_o, lam, a, n, T1, T2

      SELECT CASE (lDmn%visc%viscType)
      CASE (viscType_Const)
         mu   = lDmn%visc%mu_i
         mu_s = mu
         mu_x = 0._RKIND

      CASE (viscType_CY)
         mu_i = lDmn%visc%mu_i
         mu_o = lDmn%visc%mu_o
         lam  = lDmn%visc%lam
         a    = lDmn%visc%a
         n    = lDmn%visc%n

         T1   = 1._RKIND + (lam*gamma)**a
         T2   = T1**((n-1._RKIND)/a)
         mu   = mu_i + (mu_o-mu_i)*T2
         mu_s = mu_i

         T1   = T2/T1
         T2   = lam**a * gamma**(a-1._RKIND) * T1
         mu_x = (mu_o-mu_i)*(n-1._RKIND)*T2

      CASE (viscType_Cass)
         mu_i = lDmn%visc%mu_i
         mu_o = lDmn%visc%mu_o
         lam  = lDmn%visc%lam

         IF (gamma .LT. lam) THEN
            mu_o  = mu_o/SQRT(lam)
            gamma = lam
         ELSE
            mu_o = mu_o/SQRT(gamma)
         END IF
         mu   = (mu_i + mu_o) * (mu_i + mu_o)
         mu_s = mu_i*mu_i
         mu_x = 2._RKIND*mu_o*(mu_o + mu_i)/gamma

      END SELECT

      RETURN
      END SUBROUTINE GETVISCOSITY
!####################################################################

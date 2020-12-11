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

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      REAL(KIND=RKIND) :: lStab
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nwxx(:,:), Nqx(:,:)
      INTEGER(KIND=IKIND) :: l

      eNoN = lM%eNoN

      IF (lM%nFs .EQ. 1) THEN
         lStab = 1._RKIND
      ELSE
         lStab = 0._RKIND
      END IF

!     l = 3 if nsd == 2 else 6
      l = 3*(nsd-1)

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
         CALL GETTHOODFS(fs, lM, lStab, 1)

!        Define element coordinates appropriate for function spaces
         ALLOCATE(xwl(nsd,fs(1)%eNoN), Nwx(nsd,fs(1)%eNoN))
         ALLOCATE(Nwxx(l,fs(1)%eNoN))
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)
         Nwx      = 0._RKIND 
         Nqx      = 0._RKIND 
         Nwxx     = 0._RKIND

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               IF (lStab .LT. 1._RKIND-eps) 
     2            CALL GNNxx(l, fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), 
     3            fs(1)%Nxx(:,:,g), xwl, Nwx, Nwxx)
            END IF
            w = fs(1)%w(g) * Jac

            IF (g.EQ.1 .OR. .NOT.fs(2)%lShpF) THEN
               CALL GNN(fs(2)%eNoN, nsd, fs(2)%Nx(:,:,g), xql, Nqx, Jac,
     2            ksix)
            END IF

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D_M(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix, 
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl, 
     3            bfl, lR, lK)
            ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_M(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix, 
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl, 
     3            bfl, lR, lK)
            END IF
         END DO ! g: loop

!        Set function spaces for velocity and pressure.
         CALL GETTHOODFS(fs, lM, lStab, 2)

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
               CALL FLUID3D_C(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, bfl, lR,
     3            lK)
            ELSE IF (nsd .EQ. 2) THEN
               CALL FLUID2D_C(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, bfl, lR,
     3            lK)
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

      RETURN
      END SUBROUTINE CONSTRUCT_FLUID
!####################################################################
      SUBROUTINE FLUID3D_M(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(3,eNoNw), Nqx(3,eNoNq), Nwxx(6,eNoNw), 
     3   al(tDof,eNoNw), yl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      REAL(KIND=RKIND), PARAMETER :: ct(2) = (/1._RKIND, 36._RKIND/)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, fb(3), wr, wrl, wl, amd, ud(3), u(3),
     2   uxx(3), up(3), ux(3,3), p, px(3), es(3,3), gam, es_x(3,eNoNw),
     3   s, nu, nu_s, nu_x, kU, kS, tauM, tauC, tauB, divU, ua(3),
     4   rV(3), rM(3,3), T1, T2, NxdNx, udNx(eNoNw), updNx(eNoNw),
     5   uadNx(eNoNw), sumxx(eNoNw)

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      wr   = w*rho
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wrl  = wr*T1
      wl   = w*T1

!     Velocity and its gradients, body force
      u   =  0._RKIND
      ux  =  0._RKIND
      uxx = 0._RKIND
      ud  = -fb
      DO a=1, eNoNw
         u(1)    = u(1)  + Nw(a)*yl(1,a)
         u(2)    = u(2)  + Nw(a)*yl(2,a)
         u(3)    = u(3)  + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nwx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nwx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nwx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nwx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nwx(3,a)*yl(3,a)

         sumxx(a)= Nwxx(1,a) + Nwxx(2,a) + Nwxx(3,a)
         uxx(1)  = uxx(1)  + sumxx(a)*yl(1,a)
         uxx(2)  = uxx(2)  + sumxx(a)*yl(2,a)
         uxx(3)  = uxx(3)  + sumxx(a)*yl(3,a)

         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
         ud(3)   = ud(3) + Nw(a)*(al(3,a)-bfl(3,a))
      END DO
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(5,a)
            u(2) = u(2) - Nw(a)*yl(6,a)
            u(3) = u(3) - Nw(a)*yl(7,a)
         END DO
      END IF
!     Pressure
      p  = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(4,a)
         px(1) = px(1) + Nqx(1,a)*yl(4,a)
         px(2) = px(2) + Nqx(2,a)*yl(4,a)
         px(3) = px(3) + Nqx(3,a)*yl(4,a)
      END DO

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)

      es(1,2) = ux(1,2) + ux(2,1)
      es(1,3) = ux(1,3) + ux(3,1)
      es(2,3) = ux(2,3) + ux(3,2)

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)
    
      DO a=1, eNoNw
         es_x(1,a) = es(1,1)*Nwx(1,a)+es(1,2)*Nwx(2,a)+es(1,3)*Nwx(3,a)
         es_x(2,a) = es(2,1)*Nwx(1,a)+es(2,2)*Nwx(2,a)+es(2,3)*Nwx(3,a)
         es_x(3,a) = es(3,1)*Nwx(1,a)+es(3,2)*Nwx(2,a)+es(3,3)*Nwx(3,a)
      END DO
 
!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(1,3)*es(1,3) +
     2      es(2,1)*es(2,1) + es(2,2)*es(2,2) + es(2,3)*es(2,3) +
     3      es(3,1)*es(3,1) + es(3,2)*es(3,2) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0._RKIND
      ELSE
         nu_x = nu_x/rho/gam
      END IF
      s  = nu/eq(cEq)%dmn(cDmn)%prop(permeability)

!     Stabilization coefficients   
!     SUPG for both P2P1 and P1P1   
      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(3)*u(1)*ksix(3,1) + u(1)*u(2)*ksix(1,2)
     3   + u(2)*u(2)*ksix(2,2) + u(3)*u(2)*ksix(3,2)
     4   + u(1)*u(3)*ksix(1,3) + u(2)*u(3)*ksix(2,3)
     5   + u(3)*u(3)*ksix(3,3)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(3,1)*ksix(3,1) + ksix(1,2)*ksix(1,2)
     3   + ksix(2,2)*ksix(2,2) + ksix(3,2)*ksix(3,2)
     4   + ksix(1,3)*ksix(1,3) + ksix(2,3)*ksix(2,3)
     5   + ksix(3,3)*ksix(3,3)

      tauM = 1._RKIND / SQRT( (2._RKIND*ct(1)/dt)**2._RKIND + kU +
     2   ct(2)*nu_s*nu_s*kS + s*s ) !* lStab

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + u(3)*ux(3,1) - nu*uxx(1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + u(3)*ux(3,2) - nu*uxx(2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3)
     2      + u(3)*ux(3,3) - nu*uxx(3) + s*u(3))

      tauC = ksix(1,1) + ksix(2,2) + ksix(3,3)
      tauC = 1._RKIND/(tauM+eps)/tauC * lStab

      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1)
     2     + up(3)*up(1)*ksix(3,1) + up(1)*up(2)*ksix(1,2)
     3     + up(2)*up(2)*ksix(2,2) + up(3)*up(2)*ksix(3,2)
     4     + up(1)*up(3)*ksix(1,3) + up(2)*up(3)*ksix(2,3)
     5     + up(3)*up(3)*ksix(3,3)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1._RKIND/SQRT(tauB) * lStab
      
      ua(1) = u(1) + up(1) * lStab
      ua(2) = u(2) + up(2) * lStab
      ua(3) = u(3) + up(3) * lStab

      divU = ux(1,1) + ux(2,2) + ux(3,3)

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      rV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))

      rM(1,1) = nu*es(1,1) - up(1)*ua(1) + rV(1)*up(1) + divU*tauC
     2        - p/rho
      rM(2,1) = nu*es(2,1) - up(2)*ua(1) + rV(2)*up(1)
      rM(3,1) = nu*es(3,1) - up(3)*ua(1) + rV(3)*up(1)

      rM(1,2) = nu*es(1,2) - up(1)*ua(2) + rV(1)*up(2)
      rM(2,2) = nu*es(2,2) - up(2)*ua(2) + rV(2)*up(2) + divU*tauC
     2        - p/rho
      rM(3,2) = nu*es(3,2) - up(3)*ua(2) + rV(3)*up(2)

      rM(1,3) = nu*es(1,3) - up(1)*ua(3) + rV(1)*up(3)
      rM(2,3) = nu*es(2,3) - up(2)*ua(3) + rV(2)*up(3)
      rM(3,3) = nu*es(3,3) - up(3)*ua(3) + rV(3)*up(3) + divU*tauC
     2        - p/rho

      rV(1) = ud(1) + ua(1)*(s+ux(1,1)) + ua(2)*ux(2,1) + ua(3)*ux(3,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*(s+ux(2,2)) + ua(3)*ux(3,2)
      rV(3) = ud(3) + ua(1)*ux(1,3) + ua(2)*ux(2,3) + ua(3)*(s+ux(3,3))

!     Local residue
      DO a=1, eNoNw
         udNx(a)  = u(1)*Nwx(1,a)  + u(2)*Nwx(2,a)  + u(3)*Nwx(3,a)
         updNx(a) = up(1)*Nwx(1,a) + up(2)*Nwx(2,a) + up(3)*Nwx(3,a)
         uadNx(a) = updNx(a) * lStab + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*Nw(a) + rM(1,1)*Nwx(1,a)
     2      + rM(1,2)*Nwx(2,a) + rM(1,3)*Nwx(3,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*Nw(a) + rM(2,1)*Nwx(1,a)
     2      + rM(2,2)*Nwx(2,a) + rM(2,3)*Nwx(3,a))

         lR(3,a) = lR(3,a) + wr*(rV(3)*Nw(a) + rM(3,1)*Nwx(1,a)
     2      + rM(3,2)*Nwx(2,a) + rM(3,3)*Nwx(3,a))
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

            NxdNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b) + 
     2              Nwx(3,a)*Nwx(3,b)

            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + Nw(a)*((s+amd)*Nw(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*Nw(b) - nu*sumxx(b))

!           dM/dU
            lK(1,a,b)  = lK(1,a,b)  + wrl*((nu + tauC)*rM(1,1) + T1 +
     2         nu_x*es_x(1,a)*es_x(1,b))
            lK(2,a,b)  = lK(2,a,b)  + wrl*(nu*rM(2,1) + tauC*rM(1,2) +
     2         nu_x*es_x(1,a)*es_x(2,b))
            lK(3,a,b)  = lK(3,a,b)  + wrl*(nu*rM(3,1) + tauC*rM(1,3) +
     2         nu_x*es_x(1,a)*es_x(3,b))

            lK(5,a,b)  = lK(5,a,b)  + wrl*(nu*rM(1,2) + tauC*rM(2,1) +
     2         nu_x*es_x(2,a)*es_x(1,b))
            lK(6,a,b)  = lK(6,a,b)  + wrl*((nu + tauC)*rM(2,2) + T1 +
     2         nu_x*es_x(2,a)*es_x(2,b))
            lK(7,a,b)  = lK(7,a,b)  + wrl*(nu*rM(3,2) + tauC*rM(2,3) +
     2         nu_x*es_x(2,a)*es_x(3,b))

            lK(9,a,b)  = lK(9,a,b)  + wrl*(nu*rM(1,3) + tauC*rM(3,1) +
     2         nu_x*es_x(3,a)*es_x(1,b))
            lK(10,a,b) = lK(10,a,b) + wrl*(nu*rM(2,3) + tauC*rM(3,2) +
     2         nu_x*es_x(3,a)*es_x(2,b))
            lK(11,a,b) = lK(11,a,b) + wrl*((nu + tauC)*rM(3,3) + T1 +
     2         nu_x*es_x(3,a)*es_x(3,b))

         END DO 
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T2 = tauM*uadNx(a)
!           dM/dP
            lK(4,a,b)  = lK(4,a,b)  - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T2)
            lK(8,a,b)  = lK(8,a,b)  - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T2)
            lK(12,a,b) = lK(12,a,b) - wl*(Nwx(3,a)*Nq(b) - Nqx(3,b)*T2)
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID3D_M
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D_M(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, Nwxx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(2,eNoNw), Nqx(2,eNoNq), Nwxx(3,eNoNw), 
     3   al(tDof,eNoNw), yl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      REAL(KIND=RKIND), PARAMETER :: ct(2) = (/1._RKIND, 36._RKIND/)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, fb(2), wr, wrl, wl, amd, ud(2), u(2),
     2   uxx(2), up(2), ux(2,2), p, px(2), es(2,2), gam, es_x(2,eNoNw),
     3   s, nu, nu_s, nu_x, kU, kS, tauM, tauC, tauB, divU, ua(2),
     4   rV(2), rM(2,2), T1, T2, NxdNx, udNx(eNoNw), updNx(eNoNw),
     5   uadNx(eNoNw), sumxx(eNoNw)

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)

      wr   = w*rho
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wrl  = wr*T1
      wl   = w*T1

!     Velocity and its gradients, body force
      u   =  0._RKIND
      ux  =  0._RKIND
      uxx =  0._RKIND
      ud  = -fb
      DO a=1, eNoNw

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

         sumxx(a)= Nwxx(1,a) + Nwxx(2,a)
         uxx(1)  = uxx(1)  + sumxx(a)*yl(1,a)
         uxx(2)  = uxx(2)  + sumxx(a)*yl(2,a)         

         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
      END DO
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(4,a)
            u(2) = u(2) - Nw(a)*yl(5,a)
         END DO
      END IF

!     Pressure
      p  = 0._RKIND
      px = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(3,a)

         px(1) = px(1) + Nqx(1,a)*yl(3,a)
         px(2) = px(2) + Nqx(2,a)*yl(3,a)
      END DO

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

      DO a=1, eNoNw
        es_x(1,a) = es(1,1)*Nwx(1,a) + es(1,2)*Nwx(2,a)
        es_x(2,a) = es(2,1)*Nwx(1,a) + es(2,2)*Nwx(2,a)
      END DO

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0._RKIND
      ELSE
         nu_x = nu_x/rho/gam
      END IF
      s  = nu/eq(cEq)%dmn(cDmn)%prop(permeability)

!     Stabilization coefficients
!     SUPG for both P2P1 and P1P1         
      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)

      tauM = 1._RKIND / SQRT( (2._RKIND*ct(1)/dt)**2._RKIND + kU +
     2   ct(2)*nu_s*nu_s*kS + s*s ) !* lStab

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      - nu*uxx(1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      - nu*uxx(2) + s*u(2))

      tauC = ksix(1,1) + ksix(2,2)
      tauC = 1._RKIND/(tauM+eps)/tauC * lStab

      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1)
     2     + up(1)*up(2)*ksix(1,2) + up(2)*up(2)*ksix(2,2)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1._RKIND/SQRT(tauB) * lStab

      ua(1) = u(1) + up(1) * lStab
      ua(2) = u(2) + up(2) * lStab

      divU = ux(1,1) + ux(2,2)

      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))

      rM(1,1) = nu*es(1,1) - up(1)*ua(1) + rV(1)*up(1) + divU*tauC
     2        - p/rho
      rM(2,1) = nu*es(2,1) - up(2)*ua(1) + rV(2)*up(1)

      rM(1,2) = nu*es(1,2) - up(1)*ua(2) + rV(1)*up(2)
      rM(2,2) = nu*es(2,2) - up(2)*ua(2) + rV(2)*up(2) + divU*tauC
     2        - p/rho

      rV(1) = ud(1) + ua(1)*(s + ux(1,1)) + ua(2)*ux(2,1)
      rV(2) = ud(2) + ua(1)*ux(1,2) + ua(2)*(s + ux(2,2))

!     Local residue
      DO a=1, eNoNw
         udNx(a)  = u(1)*Nwx(1,a)  + u(2)*Nwx(2,a)
         updNx(a) = up(1)*Nwx(1,a) + up(2)*Nwx(2,a)
         uadNx(a) = updNx(a) * lStab + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*Nw(a) + rM(1,1)*Nwx(1,a)
     2      + rM(1,2)*Nwx(2,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*Nw(a) + rM(2,1)*Nwx(1,a)
     2      + rM(2,2)*Nwx(2,a))

      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            rM(1,1) = Nwx(1,a)*Nwx(1,b)
            rM(2,1) = Nwx(2,a)*Nwx(1,b)
            rM(1,2) = Nwx(1,a)*Nwx(2,b)
            rM(2,2) = Nwx(2,a)*Nwx(2,b)

            NxdNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b)

            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + Nw(a)*((s+amd)*Nw(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*Nw(b) - nu*sumxx(b))

!           dM/dU
            lK(1,a,b)  = lK(1,a,b) + wrl*((nu + tauC)*rM(1,1) + T1 +
     2         nu_x*es_x(1,a)*es_x(1,b))
            lK(2,a,b)  = lK(2,a,b) + wrl*(nu*rM(2,1) + tauC*rM(1,2) +
     2         nu_x*es_x(1,a)*es_x(2,b))

            lK(4,a,b)  = lK(4,a,b) + wrl*(nu*rM(1,2) + tauC*rM(2,1) +
     2         nu_x*es_x(2,a)*es_x(1,b))
            lK(5,a,b)  = lK(5,a,b) + wrl*((nu + tauC)*rM(2,2) + T1 +
     2         nu_x*es_x(2,a)*es_x(2,b))
         END DO 
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            T2 = tauM*uadNx(a)
!           dM/dP
            lK(3,a,b)  = lK(3,a,b) - wl*(Nwx(1,a)*Nq(b) - Nqx(1,b)*T2)
            lK(6,a,b)  = lK(6,a,b) - wl*(Nwx(2,a)*Nq(b) - Nqx(2,b)*T2)
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID2D_M
!####################################################################
      SUBROUTINE FLUID3D_C(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(3,eNoNw), Nqx(3,eNoNq), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      REAL(KIND=RKIND), PARAMETER :: ct(2) = (/1._RKIND, 36._RKIND/)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, fb(3), wl, amd, ud(3), u(3), 
     2   up(3), ux(3,3), px(3), es(3,3), gam, s, nu, nu_s,
     3   nu_x, kU, kS, tauM, divU, T1, NxdNx, updNx

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

!     Velocity and its gradients, body force
      u  =  0._RKIND
      ux =  0._RKIND
      ud = -fb
      DO a=1, eNoNw
         u(1)    = u(1)  + Nw(a)*yl(1,a)
         u(2)    = u(2)  + Nw(a)*yl(2,a)
         u(3)    = u(3)  + Nw(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nwx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nwx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nwx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nwx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nwx(3,a)*yl(3,a)

         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
         ud(3)   = ud(3) + Nw(a)*(al(3,a)-bfl(3,a))
      END DO
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(5,a)
            u(2) = u(2) - Nw(a)*yl(6,a)
            u(3) = u(3) - Nw(a)*yl(7,a)
         END DO
      END IF

!     Pressure
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(4,a)
         px(2) = px(2) + Nqx(2,a)*yl(4,a)
         px(3) = px(3) + Nqx(3,a)*yl(4,a)
      END DO

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)

      es(1,2) = ux(1,2) + ux(2,1)
      es(1,3) = ux(1,3) + ux(3,1)
      es(2,3) = ux(2,3) + ux(3,2)

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(1,3)*es(1,3) +
     2      es(2,1)*es(2,1) + es(2,2)*es(2,2) + es(2,3)*es(2,3) +
     3      es(3,1)*es(3,1) + es(3,2)*es(3,2) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0._RKIND
      ELSE
         nu_x = nu_x/rho/gam
      END IF
      s  = nu/eq(cEq)%dmn(cDmn)%prop(permeability)

!     Stabilization coefficients      
      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(3)*u(1)*ksix(3,1) + u(1)*u(2)*ksix(1,2)
     3   + u(2)*u(2)*ksix(2,2) + u(3)*u(2)*ksix(3,2)
     4   + u(1)*u(3)*ksix(1,3) + u(2)*u(3)*ksix(2,3)
     5   + u(3)*u(3)*ksix(3,3)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(3,1)*ksix(3,1) + ksix(1,2)*ksix(1,2)
     3   + ksix(2,2)*ksix(2,2) + ksix(3,2)*ksix(3,2)
     4   + ksix(1,3)*ksix(1,3) + ksix(2,3)*ksix(2,3)
     5   + ksix(3,3)*ksix(3,3)

      tauM = 1._RKIND / SQRT( (2._RKIND*ct(1)/dt)**2._RKIND + kU +
     2       ct(2)*nu_s*nu_s*kS + s*s ) * lStab

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + u(3)*ux(3,1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + u(3)*ux(3,2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3)
     2      + u(3)*ux(3,3) + s*u(3))

      divU = ux(1,1) + ux(2,2) + ux(3,3)

!     Local residue
      DO a=1, eNoNq
         updNx   = up(1)*Nqx(1,a) + up(2)*Nqx(2,a) + up(3)*Nqx(3,a)
         lR(4,a) = lR(4,a) + w*(Nq(a)*divU - updNx)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         T1 = tauM*(amd*Nw(b) + u(1)*Nwx(1,b) + u(2)*Nwx(2,b) + 
     2      u(3)*Nwx(3,b))
         DO a=1, eNoNq
!           dC/dU
            lK(13,a,b) = lK(13,a,b) + wl*(Nwx(1,b)*Nq(a) + Nqx(1,a)*T1)
            lK(14,a,b) = lK(14,a,b) + wl*(Nwx(2,b)*Nq(a) + Nqx(2,a)*T1)
            lK(15,a,b) = lK(15,a,b) + wl*(Nwx(3,b)*Nq(a) + Nqx(3,a)*T1)
         END DO 
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNq
            NxdNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b) + 
     2              Nqx(3,a)*Nqx(3,b)
!           dC/dP
            lK(16,a,b) = lK(16,a,b) + wl*tauM*NxdNx/rho
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID3D_C
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D_C(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(2,eNoNw), Nqx(2,eNoNq), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      REAL(KIND=RKIND), PARAMETER :: ct(2) = (/1._RKIND, 36._RKIND/)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, fb(2), wl, amd, ud(2), u(2), 
     2   up(2), ux(2,2), px(2), es(2,2), gam, s, nu, nu_s, 
     3   nu_x, kU, kS, tauM, divU, T1, NxdNx, updNx

!     Define parameters
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

!     Velocity and its gradients, body force
      u  =  0._RKIND
      ux =  0._RKIND
      ud = -fb
      DO a=1, eNoNw

         u(1)    = u(1) + Nw(a)*yl(1,a)
         u(2)    = u(2) + Nw(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nwx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nwx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nwx(2,a)*yl(2,a)

         ud(1)   = ud(1) + Nw(a)*(al(1,a)-bfl(1,a))
         ud(2)   = ud(2) + Nw(a)*(al(2,a)-bfl(2,a))
      END DO
      IF (mvMsh) THEN
         DO a=1, eNoNw
            u(1) = u(1) - Nw(a)*yl(4,a)
            u(2) = u(2) - Nw(a)*yl(5,a)
         END DO
      END IF

!     Pressure
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(3,a)
         px(2) = px(2) + Nqx(2,a)*yl(3,a)
      END DO

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0._RKIND
      ELSE
         nu_x = nu_x/rho/gam
      END IF
      s  = nu/eq(cEq)%dmn(cDmn)%prop(permeability)

!     Stabilization coefficients      
      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)

      tauM = 1._RKIND / SQRT( (2._RKIND*ct(1)/dt)**2._RKIND + kU +
     2       ct(2)*nu_s*nu_s*kS + s*s ) * lStab

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + s*u(2))
      
      divU = ux(1,1) + ux(2,2)

!     Local residue
      DO a=1, eNoNq
         updNx   = up(1)*Nqx(1,a) + up(2)*Nqx(2,a)
         lR(3,a) = lR(3,a) + w*(Nq(a)*divU - updNx)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         T1 = tauM*(amd*Nw(b) + u(1)*Nwx(1,b) + u(2)*Nwx(2,b))
         DO a = 1, eNoNq
!           dC/dU
            lK(7,a,b) = lK(7,a,b) + wl*(Nwx(1,b)*Nq(a) + Nqx(1,a)*T1)
            lK(8,a,b) = lK(8,a,b) + wl*(Nwx(2,b)*Nq(a) + Nqx(2,a)*T1)
         END DO 
      END DO 

      DO b = 1, eNoNq
         DO a = 1, eNoNq
            NxdNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b)
!           dC/dP
            lK(9,a,b) = lK(9,a,b) + wl*tauM*NxdNx/rho
         END DO 
      END DO 

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
      SUBROUTINE BWFLUID3D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   yl(tDof,eNoN), ub(3), nV(3), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn,
     2   ubn, u(3), uh(3), ux(3,3), sigman(3), Nxn(eNoN), rV(3),
     3   rM(3,3), nu_s, es(3,3), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w * rho
      wl   = w * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0._RKIND
      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoN
         p = p + N(a)*yl(4,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*yl(3,a)

         Nxn(a)  = Nx(1,a)*nV(1) + Nx(2,a)*nV(2) + Nx(3,a)*nV(3)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(5,a)
            uh(2) = uh(2) + N(a)*yl(6,a)
            uh(3) = uh(3) + N(a)*yl(7,a)
         END DO
      END IF
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2) + (u(3)-ub(3))*nV(3)
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2) + (u(3)-uh(3))*nV(3)
      uhn = (ABS(uhn) - uhn) * 0.5_RKIND

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)
      es(3,3) = ux(3,3) + ux(3,3)

      es(1,2) = ux(1,2) + ux(2,1)
      es(1,3) = ux(1,3) + ux(3,1)
      es(2,3) = ux(2,3) + ux(3,2)

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(1,3)*es(1,3) +
     2      es(2,1)*es(2,1) + es(2,2)*es(2,2) + es(2,3)*es(2,3) +
     3      es(3,1)*es(3,1) + es(3,2)*es(3,2) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu = nu/rho

      sigman(1) = nu*(es(1,1)*nV(1) + es(1,2)*nV(2) + es(1,3)*nV(3))
     2   - (p/rho)*nV(1)
      sigman(2) = nu*(es(2,1)*nV(1) + es(2,2)*nV(2) + es(2,3)*nV(3))
     2   - (p/rho)*nV(2)
      sigman(3) = nu*(es(3,1)*nV(1) + es(3,2)*nV(2) + es(3,3)*nV(3))
     2   - (p/rho)*nV(3)

      rV(1) = -sigman(1) + (tauT + uhn)*(u(1)-ub(1))
      rV(2) = -sigman(2) + (tauT + uhn)*(u(2)-ub(2))
      rV(3) = -sigman(3) + (tauT + uhn)*(u(3)-ub(3))

      DO a=1, eNoN
c         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
c     2      nu*(Nxn(a)*(u(1)-ub(1)) + Nx(1,a)*ubn) )
         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
     2      nu*(u(1)-ub(1))*(Nxn(a) + Nx(1,a)*nV(1)) )

c         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
c     2      nu*(Nxn(a)*(u(2)-ub(2)) + Nx(2,a)*ubn) )
         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
     2      nu*(u(2)-ub(2))*(Nxn(a) + Nx(2,a)*nV(2)) )

c         lR(3,a) = lR(3,a) + wr*( N(a)*(rV(3) + (tauT-tauN)*ubn*nV(3)) -
c     2      nu*(Nxn(a)*(u(3)-ub(3)) + Nx(3,a)*ubn) )
         lR(3,a) = lR(3,a) + wr*( N(a)*(rV(3) + (tauT-tauN)*ubn*nV(3)) -
     2      nu*(u(3)-ub(3))*(Nxn(a) + Nx(3,a)*nV(3)) )

         lR(4,a) = lR(4,a) - w*N(a)*ubn
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            T1 = (tauT + uhn)*N(a)*N(b) - nu*(N(a)*Nxn(b) + N(b)*Nxn(a))

            rM(1,1) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(1) -
     2         nu*(N(a)*Nx(1,b)*nV(1) + Nx(1,a)*N(b)*nV(1))
c            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
c     2         nu*(N(a)*Nx(1,b)*nV(2) + Nx(1,a)*N(b)*nV(2))
c            rM(1,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
c     2         nu*(N(a)*Nx(1,b)*nV(3) + Nx(1,a)*N(b)*nV(3))
            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
     2         nu*N(a)*Nx(1,b)*nV(2)
            rM(1,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
     2         nu*N(a)*Nx(1,b)*nV(3)

c            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
c     2         nu*(N(a)*Nx(2,b)*nV(1) + Nx(2,a)*N(b)*nV(1))
            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
     2         nu*N(a)*Nx(2,b)*nV(1)
            rM(2,2) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(2) -
     2         nu*(N(a)*Nx(2,b)*nV(2) + Nx(2,a)*N(b)*nV(2))
c            rM(2,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
c     2         nu*(N(a)*Nx(2,b)*nV(3) + Nx(2,a)*N(b)*nV(3))
            rM(2,3) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(3) -
     2         nu*N(a)*Nx(2,b)*nV(3)

c            rM(3,1) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(1) -
c     2         nu*(N(a)*Nx(3,b)*nV(1) + Nx(3,a)*N(b)*nV(1))
c            rM(3,2) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(2) -
c     2         nu*(N(a)*Nx(3,b)*nV(2) + Nx(3,a)*N(b)*nV(2))
            rM(3,1) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(1) -
     2         nu*N(a)*Nx(3,b)*nV(1)
            rM(3,2) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(2) -
     2         nu*N(a)*Nx(3,b)*nV(2)
            rM(3,3) = (tauT-tauN)*N(a)*N(b)*nV(3)*nV(3) -
     2         nu*(N(a)*Nx(3,b)*nV(3) + Nx(3,a)*N(b)*nV(3))

!           dM/dU
            lK(1,a,b) = lK(1,a,b) + wrl*(T1 + rM(1,1))
            lK(2,a,b) = lK(2,a,b) + wrl*rM(1,2)
            lK(3,a,b) = lK(3,a,b) + wrl*rM(1,3)

            lK(5,a,b) = lK(5,a,b) + wrl*rM(2,1)
            lK(6,a,b) = lK(6,a,b) + wrl*(T1 + rM(2,2))
            lK(7,a,b) = lK(7,a,b) + wrl*rM(2,3)

            lK(9,a,b)  = lK(9,a,b)  + wrl*rM(3,1)
            lK(10,a,b) = lK(10,a,b) + wrl*rM(3,2)
            lK(11,a,b) = lK(11,a,b) + wrl*(T1 + rM(3,3))

!           dM/dP
            lK(4,a,b)  = lK(4,a,b)  + wl*N(a)*N(b)*nV(1)
            lK(8,a,b)  = lK(8,a,b)  + wl*N(a)*N(b)*nV(2)
            lK(12,a,b) = lK(12,a,b) + wl*N(a)*N(b)*nV(3)

!           dC/dU
            lK(13,a,b) = lK(13,a,b) - wl*N(a)*N(b)*nV(1)
            lK(14,a,b) = lK(14,a,b) - wl*N(a)*N(b)*nV(2)
            lK(15,a,b) = lK(15,a,b) - wl*N(a)*N(b)*nV(3)
         END DO
      END DO

      RETURN
      END SUBROUTINE BWFLUID3D
!--------------------------------------------------------------------
      SUBROUTINE BWFLUID2D(eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   yl(tDof,eNoN), ub(2), nV(2), tauB(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b
      REAL(KIND=RKIND) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn,
     2   ubn, u(2), uh(2), ux(2,2), sigman(2), Nxn(eNoN), rV(2),
     3   rM(2,2), nu_s, es(2,2), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w*rho
      wl   = w  * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0._RKIND
      u    = 0._RKIND
      ux   = 0._RKIND
      DO a=1, eNoN
         p = p + N(a)*yl(3,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)

         Nxn(a)  = Nx(1,a)*nV(1) + Nx(2,a)*nV(2)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(4,a)
            uh(2) = uh(2) + N(a)*yl(5,a)
         END DO
      END IF
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2)
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2)
      uhn = (ABS(uhn) - uhn) * 0.5_RKIND

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu = nu/rho

      sigman(1) = -(p/rho)*nV(1) + nu*(es(1,1)*nV(1) + es(1,2)*nV(2))
      sigman(2) = -(p/rho)*nV(2) + nu*(es(2,1)*nV(1) + es(2,2)*nV(2))

      rV(1) = -sigman(1) + (tauT + uhn)*(u(1)-ub(1))
      rV(2) = -sigman(2) + (tauT + uhn)*(u(2)-ub(2))

      DO a=1, eNoN
c         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
c     2      nu*(Nxn(a)*(u(1)-ub(1)) + Nx(1,a)*ubn) )
         lR(1,a) = lR(1,a) + wr*( N(a)*(rV(1) + (tauT-tauN)*ubn*nV(1)) -
     2      nu*(u(1)-ub(1))*(Nxn(a) + Nx(1,a)*nV(1)) )

c         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
c     2      nu*(Nxn(a)*(u(2)-ub(2)) + Nx(2,a)*ubn) )
         lR(2,a) = lR(2,a) + wr*( N(a)*(rV(2) + (tauT-tauN)*ubn*nV(2)) -
     2      nu*(u(2)-ub(2))*(Nxn(a) + Nx(2,a)*nV(2)) )

         lR(3,a) = lR(3,a) - w*N(a)*ubn
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            T1 = (tauT + uhn)*N(a)*N(b) - nu*(N(a)*Nxn(b) + N(b)*Nxn(a))

            rM(1,1) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(1) -
     2         nu*(N(a)*Nx(1,b)*nV(1) + Nx(1,a)*N(b)*nV(1))
c            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
c     2         nu*(N(a)*Nx(1,b)*nV(2) + Nx(1,a)*N(b)*nV(2))
            rM(1,2) = (tauT-tauN)*N(a)*N(b)*nV(1)*nV(2) -
     2         nu*N(a)*Nx(1,b)*nV(2)

c            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
c     2         nu*(N(a)*Nx(2,b)*nV(1) + Nx(2,a)*N(b)*nV(1))
            rM(2,1) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(1) -
     2         nu*N(a)*Nx(2,b)*nV(1)
            rM(2,2) = (tauT-tauN)*N(a)*N(b)*nV(2)*nV(2) -
     2         nu*(N(a)*Nx(2,b)*nV(2) + Nx(2,a)*N(b)*nV(2))

!           dM/dU
            lK(1,a,b) = lK(1,a,b) + wrl*(T1 + rM(1,1))
            lK(2,a,b) = lK(2,a,b) + wrl*rM(1,2)

            lK(4,a,b) = lK(4,a,b) + wrl*rM(2,1)
            lK(5,a,b) = lK(5,a,b) + wrl*(T1 + rM(2,2))

!           dM/dP
            lK(3,a,b) = lK(4,a,b) + wl*N(a)*N(b)*nV(1)
            lK(6,a,b) = lK(8,a,b) + wl*N(a)*N(b)*nV(2)

!           dC/dU
            lK(7,a,b) = lK(7,a,b) - wl*N(a)*N(b)*nV(1)
            lK(8,a,b) = lK(8,a,b) - wl*N(a)*N(b)*nV(2)
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

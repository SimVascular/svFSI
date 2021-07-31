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
!     These routines are for solving the Stokes equations.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_STOKES(lM, Ag, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)
      TYPE(fsType) :: fs(2)

      LOGICAL lStab
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xwl(:,:), xql(:,:), Nwx(:,:),
     2   Nqx(:,:)

      eNoN = lM%eNoN

      IF (lM%nFs .EQ. 1) THEN
         lStab = .TRUE.
      ELSE
         lStab = .FALSE.
      END IF

!     STOKES: dof = nsd+1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   bfl(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_stokes) CYCLE

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
         ALLOCATE(xql(nsd,fs(2)%eNoN), Nqx(nsd,fs(2)%eNoN))
         xwl(:,:) = xl(:,:)
         xql(:,:) = xl(:,1:fs(2)%eNoN)

!        Gauss integration 1
         DO g=1, fs(1)%nG
            IF (g.EQ.1 .OR. .NOT.fs(1)%lShpF) THEN
               CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac,
     2            ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = fs(1)%w(g) * Jac

            IF (nsd .EQ. 3) THEN
               CALL STOKES3D_M(fs(1)%eNoN, fs(2)%eNoN, w, fs(1)%N(:,g),
     2            fs(2)%N(:,g), Nwx, al, yl, bfl, lR, lK)
            ELSE IF (nsd .EQ. 2) THEN
               CALL STOKES2D_M(fs(1)%eNoN, fs(2)%eNoN, w, fs(1)%N(:,g),
     2            fs(2)%N(:,g), Nwx, al, yl, bfl, lR, lK)
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
               CALL STOKES3D_C(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, bfl, lR,
     3            lK)
            ELSE IF (nsd .EQ. 2) THEN
               CALL STOKES2D_C(lStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, al, yl, bfl, lR,
     3            lK)
            END IF
         END DO ! g: loop

         DEALLOCATE(xwl, xql, Nwx, Nqx)

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
      END SUBROUTINE CONSTRUCT_STOKES
!####################################################################
      SUBROUTINE STOKES3D_M(eNoNw, eNoNq, w, Nw, Nq, Nwx, al, yl, bfl,
     2   lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(3,eNoNw), al(tDof,eNoNw), yl(tDof,eNoNw), bfl(nsd,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, l, a, b
      REAL(KIND=RKIND) :: mu, af, wf, wm, fb(3), vd(3), v(3), vx(3,3),
     2   p, es(3,3), rM, NxNx

!     Define parameters
      mu    = eq(cEq)%dmn(cDmn)%visc%mu_i
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)
      af    = eq(cEq)%af*eq(cEq)%gam*dt
      wf    = w * af
      wm    = w * eq(cEq)%am

!     {i,j,k} := velocity dofs; {l} := pressure dof
      i = eq(cEq)%s
      j = i + 1
      k = j + 1
      l = k + 1

!     Velocity and its gradients, body force
      v  =  0._RKIND
      vx =  0._RKIND
      vd = -fb
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)
         v(3)    = v(3)  + Nw(a)*yl(k,a)

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(1,3) = vx(1,3) + Nwx(3,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)
         vx(2,3) = vx(2,3) + Nwx(3,a)*yl(j,a)
         vx(3,1) = vx(3,1) + Nwx(1,a)*yl(k,a)
         vx(3,2) = vx(3,2) + Nwx(2,a)*yl(k,a)
         vx(3,3) = vx(3,3) + Nwx(3,a)*yl(k,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))
         vd(3)   = vd(3) + Nw(a)*(al(k,a)-bfl(3,a))
      END DO

!     Pressure
      p = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(l,a)
      END DO

!     Viscous stress tensor 2*mu*e_ij := mu*(u_ij + u_ji)
      es(1,1) = mu*(vx(1,1) + vx(1,1))
      es(2,2) = mu*(vx(2,2) + vx(2,2))
      es(3,3) = mu*(vx(3,3) + vx(3,3))

      es(1,2) = mu*(vx(1,2) + vx(2,1))
      es(1,3) = mu*(vx(1,3) + vx(3,1))
      es(2,3) = mu*(vx(2,3) + vx(3,2))

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)

!     Local residue
      DO a=1, eNoNw
         rM = Nwx(1,a)*(-p + es(1,1)) + Nwx(2,a)*es(1,2) +
     2      Nwx(3,a)*es(1,3)
         lR(1,a) = lR(1,a) + w*(Nw(a)*vd(1) + rM)

         rM = Nwx(1,a)*es(2,1) + Nwx(2,a)*(-p + es(2,2)) +
     2      Nwx(3,a)*es(2,3)
         lR(2,a) = lR(2,a) + w*(Nw(a)*vd(2) + rM)

         rM = Nwx(1,a)*es(3,1) + Nwx(2,a)*es(3,2) + Nwx(3,a)*
     3      (-p + es(3,3))
         lR(3,a) = lR(3,a) + w*(Nw(a)*vd(3) + rM)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            NxNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b) +
     2         Nwx(3,a)*Nwx(3,b)

            lK(1,a,b)  = lK(1,a,b)  + wf*mu*(Nwx(1,a)*Nwx(1,b) + NxNx)
     2                              + wm*Nw(a)*Nw(b)
            lK(2,a,b)  = lK(2,a,b)  + wf*mu*(Nwx(2,a)*Nwx(1,b))
            lK(3,a,b)  = lK(3,a,b)  + wf*mu*(Nwx(3,a)*Nwx(1,b))

            lK(5,a,b)  = lK(5,a,b)  + wf*mu*(Nwx(1,a)*Nwx(2,b))
            lK(6,a,b)  = lK(6,a,b)  + wf*mu*(Nwx(2,a)*Nwx(2,b) + NxNx)
     2                              + wm*Nw(a)*Nw(b)
            lK(7,a,b)  = lK(7,a,b)  + wf*mu*(Nwx(3,a)*Nwx(2,b))

            lK(9,a,b)  = lK(9,a,b)  + wf*mu*(Nwx(1,a)*Nwx(3,b))
            lK(10,a,b) = lK(10,a,b) + wf*mu*(Nwx(2,a)*Nwx(3,b))
            lK(11,a,b) = lK(11,a,b) + wf*mu*(Nwx(3,a)*Nwx(3,b) + NxNx)
     2                              + wm*Nw(a)*Nw(b)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            lK(4,a,b)  = lK(4,a,b)  - wf*Nwx(1,a)*Nq(b)
            lK(8,a,b)  = lK(8,a,b)  - wf*Nwx(2,a)*Nq(b)
            lK(12,a,b) = lK(12,a,b) - wf*Nwx(3,a)*Nq(b)
         END DO
      END DO

      RETURN
      END SUBROUTINE STOKES3D_M
!--------------------------------------------------------------------
      SUBROUTINE STOKES2D_M(eNoNw, eNoNq, w, Nw, Nq, Nwx, al, yl, bfl,
     2   lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, Nw(eNoNw), Nq(eNoNq),
     2   Nwx(2,eNoNw), al(tDof,eNoNw), yl(tDof,eNoNw), bfl(nsd,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: mu, af, wf, wm, fb(2), vd(2), v(2), vx(2,2),
     2   p, es(2,2), rM, NxNx

!     Define parameters
      mu    = eq(cEq)%dmn(cDmn)%visc%mu_i
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      af    = eq(cEq)%af*eq(cEq)%gam*dt
      wf    = w * af
      wm    = w * eq(cEq)%am

!     {i,j} := velocity dofs; {k} := pressure dof
      i = eq(cEq)%s
      j = i + 1
      k = j + 1

!     Velocity and its gradients, body force
      v  =  0._RKIND
      vx =  0._RKIND
      vd = -fb
      DO a=1, eNoNw
         v(1)    = v(1)  + Nw(a)*yl(i,a)
         v(2)    = v(2)  + Nw(a)*yl(j,a)

         vx(1,1) = vx(1,1) + Nwx(1,a)*yl(i,a)
         vx(1,2) = vx(1,2) + Nwx(2,a)*yl(i,a)
         vx(2,1) = vx(2,1) + Nwx(1,a)*yl(j,a)
         vx(2,2) = vx(2,2) + Nwx(2,a)*yl(j,a)

         vd(1)   = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2)   = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))
      END DO

!     Pressure
      p = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(k,a)
      END DO

!     Viscous stress tensor 2*mu*e_ij := mu*(u_ij + u_ji)
      es(1,1) = mu*(vx(1,1) + vx(1,1))
      es(2,2) = mu*(vx(2,2) + vx(2,2))
      es(1,2) = mu*(vx(1,2) + vx(2,1))
      es(2,1) = es(1,2)

!     Local residue
      DO a=1, eNoNw
         rM = Nwx(1,a)*(-p + es(1,1)) + Nwx(2,a)*es(1,2)
         lR(1,a) = lR(1,a) + w*(Nw(a)*vd(1) + rM)

         rM = Nwx(1,a)*es(2,1) + Nwx(2,a)*(-p + es(2,2))
         lR(2,a) = lR(2,a) + w*(Nw(a)*vd(2) + rM)
      END DO

!     Tangent (stiffness) matrices
      DO b=1, eNoNw
         DO a=1, eNoNw
            NxNx = Nwx(1,a)*Nwx(1,b) + Nwx(2,a)*Nwx(2,b)

            lK(1,a,b) = lK(1,a,b) + wf*mu*(Nwx(1,a)*Nwx(1,b) + NxNx)
     2                            + wm*Nw(a)*Nw(b)
            lK(2,a,b) = lK(2,a,b) + wf*mu*(Nwx(2,a)*Nwx(1,b))

            lK(4,a,b) = lK(4,a,b) + wf*mu*(Nwx(1,a)*Nwx(2,b))
            lK(5,a,b) = lK(5,a,b) + wf*mu*(Nwx(2,a)*Nwx(2,b) + NxNx)
     2                            + wm*Nw(a)*Nw(b)
         END DO
      END DO

      DO b=1, eNoNq
         DO a=1, eNoNw
            lK(3,a,b) = lK(3,a,b) - wf*Nwx(1,a)*Nq(b)
            lK(6,a,b) = lK(6,a,b) - wf*Nwx(2,a)*Nq(b)
         END DO
      END DO

      RETURN
      END SUBROUTINE STOKES2D_M
!####################################################################
      SUBROUTINE STOKES3D_C(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(3,eNoNw), Nqx(3,eNoNq), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(3,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, l, a, b
      REAL(KIND=RKIND) :: mu, wm, wf, div, vd(3), fb(3), px(3), kS, ctM,
     2   tauM, rMv(3), rM, NxNx

      mu    = eq(cEq)%dmn(cDmn)%visc%mu_i
      ctM   = eq(cEq)%dmn(cDmn)%prop(ctau_M)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3) = eq(cEq)%dmn(cDmn)%prop(f_z)
      wm    = w * eq(cEq)%am
      wf    = w * eq(cEq)%af*eq(cEq)%gam*dt

!     {i,j,k} := velocity dofs; {l} := pressure dof
      i = eq(cEq)%s
      j = i + 1
      k = j + 1
      l = k + 1

!     Divergence of velocity
      div =  0._RKIND
      vd  = -fb
      DO a=1, eNoNw
         vd(1) = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2) = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))
         vd(3) = vd(3) + Nw(a)*(al(k,a)-bfl(3,a))

         div   = div + Nwx(1,a)*yl(i,a) + Nwx(2,a)*yl(j,a) +
     2      Nwx(3,a)*yl(k,a)
      END DO

!     Pressure gradient
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(l,a)
         px(2) = px(2) + Nqx(2,a)*yl(l,a)
         px(3) = px(3) + Nqx(3,a)*yl(l,a)
      END DO

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(3,1)*ksix(3,1) + ksix(1,2)*ksix(1,2)
     3   + ksix(2,2)*ksix(2,2) + ksix(3,2)*ksix(3,2)
     4   + ksix(1,3)*ksix(1,3) + ksix(2,3)*ksix(2,3)
     5   + ksix(3,3)*ksix(3,3)
      kS = SQRT(kS)

!     kS is an estimate of (he**-2)
      IF (lStab) THEN
         tauM = ctM / (2._RKIND*mu*kS)
      ELSE
         tauM = 0._RKIND
      END IF

!     Local residue
      rMv(:) = vd(:) + px(:)
      DO a=1, eNoNq
         rM = rMv(1)*Nqx(1,a) + rMv(2)*Nqx(2,a) + rMv(3)*Nqx(3,a)
         lR(4,a) = lR(4,a) + w*(Nq(a)*div + tauM*rM)
      END DO

!     Tangent (stiffness) matrices
      wm = wm * tauM
      DO b=1, eNoNw
         DO a=1, eNoNq
            lK(13,a,b) = lK(13,a,b) + wm*Nqx(1,a)*Nw(b)
     2                              + wf*Nq(a)*Nwx(1,b)
            lK(14,a,b) = lK(14,a,b) + wm*Nqx(2,a)*Nw(b)
     2                              + wf*Nq(a)*Nwx(2,b)
            lK(15,a,b) = lK(15,a,b) + wm*Nqx(3,a)*Nw(b)
     2                              + wf*Nq(a)*Nwx(3,b)
         END DO
      END DO

      IF (lStab) THEN
         wf = wf * tauM
         DO b=1, eNoNq
            DO a=1, eNoNq
               NxNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b)
     2              + Nqx(3,a)*Nqx(3,b)
               lK(16,a,b) = lK(16,a,b) + wf*NxNx
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE STOKES3D_C
!--------------------------------------------------------------------
      SUBROUTINE STOKES2D_C(lStab, eNoNw, eNoNq, w, ksix, Nw, Nq, Nwx,
     2   Nqx, al, yl, bfl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: lStab
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNw, eNoNq
      REAL(KIND=RKIND), INTENT(IN) :: w, ksix(nsd,nsd), Nw(eNoNw),
     2   Nq(eNoNq), Nwx(2,eNoNw), Nqx(2,eNoNq), al(tDof,eNoNw),
     3   yl(tDof,eNoNw), bfl(2,eNoNw)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoNw),
     2   lK(dof*dof,eNoNw,eNoNw)

      INTEGER(KIND=IKIND) :: i, j, k, a, b
      REAL(KIND=RKIND) :: mu, wm, wf, div, vd(2), fb(2), px(2), kS, ctM,
     2   tauM, rMv(2), rM, NxNx

      mu    = eq(cEq)%dmn(cDmn)%visc%mu_i
      ctM   = eq(cEq)%dmn(cDmn)%prop(ctau_M)
      fb(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      wm    = w * eq(cEq)%am
      wf    = w * eq(cEq)%af*eq(cEq)%gam*dt

!     {i,j} := velocity dofs; {k} := pressure dof
      i = eq(cEq)%s
      j = i + 1
      k = j + 1

!     Divergence of velocity
      div = 0._RKIND
      vd  = -fb
      DO a=1, eNoNw
         vd(1) = vd(1) + Nw(a)*(al(i,a)-bfl(1,a))
         vd(2) = vd(2) + Nw(a)*(al(j,a)-bfl(2,a))

         div   = div + Nwx(1,a)*yl(i,a) + Nwx(2,a)*yl(j,a)
      END DO

!     Pressure gradient
      px = 0._RKIND
      DO a=1, eNoNq
         px(1) = px(1) + Nqx(1,a)*yl(k,a)
         px(2) = px(2) + Nqx(2,a)*yl(k,a)
      END DO

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)
      kS = SQRT(kS)

!     kS is an estimate of (he**-2)
      IF (lStab) THEN
         tauM = ctM / (2._RKIND*mu*kS)
      ELSE
         tauM = 0._RKIND
      END IF

!     Local residue
      rMv(:) = vd(:) + px(:)
      DO a=1, eNoNq
         rM = rMv(1)*Nqx(1,a) + rMv(2)*Nqx(2,a)
         lR(3,a) = lR(3,a) + w*(Nq(a)*div + tauM*rM)
      END DO

!     Tangent (stiffness) matrices
      wm = wm * tauM
      DO b=1, eNoNw
         DO a=1, eNoNq
            lK(7,a,b) = lK(7,a,b) + wm*Nqx(1,a)*Nw(b)
     2                            + wf*Nq(a)*Nwx(1,b)
            lK(8,a,b) = lK(8,a,b) + wm*Nqx(2,a)*Nw(b)
     2                            + wf*Nq(a)*Nwx(2,b)
         END DO
      END DO

      IF (lStab) THEN
         wf = wf * tauM
         DO b=1, eNoNq
            DO a=1, eNoNq
               NxNx = Nqx(1,a)*Nqx(1,b) + Nqx(2,a)*Nqx(2,b)
               lK(9,a,b) = lK(9,a,b) + wf*NxNx
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE STOKES2D_C
!####################################################################

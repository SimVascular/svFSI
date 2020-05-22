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
!     These subroutines implement the Coupled Momentum Method (CMM).
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_CMM(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, vwp(2), ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), bfl(:,:), pS0l(:,:), pSl(:), vwpl(:,:), N(:), Nx(:,:),
     3   lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN
      IF (lM%eType.NE.eType_TET .OR. nsd.NE.3) err =
     2   " CMM equation is constructed for 3D tetrahedral meshes only"

!     CMM: dof = nsd+1
!     CMM init: dof = nsd
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(nsd,eNoN), bfl(nsd,eNoN), pS0l(nstd,eNoN), pSl(nstd),
     3   vwpl(2,eNoN), N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN),
     4   lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_cmm) CYCLE

!        Create local copies
         pS0l = 0._RKIND
         vwpl = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            yl(:,a)  = Yg(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
            IF (ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
            IF (cmmVarWall) vwpl(:,a) = varWallProps(:,Ac)
         END DO

         IF (cmmInit) THEN
            pSl(:) = 0._RKIND
            vwp(:) = 0._RKIND
            DO a=1, eNoN
               pSl(:) = pSl(:) + pS0l(:,a)
               vwp(:) = vwp(:) + vwpl(:,a)
            END DO
            pSl(:) = pSl(:)/REAL(eNoN, KIND=RKIND)
            vwp(:) = vwp(:)/REAL(eNoN, KIND=RKIND)
            CALL CMMi(lM, eNoN, al, dl, xl, bfl, pSl, vwp, ptr)
         ELSE
!           Gauss integration
            lR = 0._RKIND
            lK = 0._RKIND
            DO g=1, lM%nG
               IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
                  CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               END IF
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
               w = lM%w(g) * Jac
               N = lM%N(:,g)

               CALL CMM3D(eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)
            END DO ! g: loop

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
         END IF
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, yl, dl, bfl, pS0l, pSl, vwpl, N, Nx, lR,
     2   lK)

      RETURN
      END SUBROUTINE CONSTRUCT_CMM
!####################################################################
!     CMM initialization (interior)
      SUBROUTINE CMMi(lM, eNoN, al, dl, xl, bfl, pS0l, vwp, ptr)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, ptr(eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: al(tDof,eNoN), dl(tDof,eNoN),
     2   xl(3,eNoN), bfl(3,eNoN), pS0l(6), vwp(2)

      INTEGER(KIND=IKIND) a, g, Ac
      REAL(KIND=RKIND) :: w, Jac, nV(3), pSl(6), gC(3,2)

      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      ALLOCATE(N(eNoN), Nx(2,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
      lK  = 0._RKIND
      lR  = 0._RKIND

      Nx  = lM%Nx(:,:,1)
      CALL GNNS(eNoN, Nx, xl, nV, gC, gC)
      Jac = SQRT(NORM(nV))
      nV  = nV / Jac

!     Internal stresses (stiffness) contribution
      pSl = 0._RKIND
      CALL CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, vwp, pSl, lR, lK)

!     Inertia and body forces (mass) contribution
      DO g=1, lM%nG
         N(:) = lM%N(:,g)
         w = lM%w(g)*Jac
         CALL CMM_MASS(eNoN, w, N, al, bfl, vwp, lR, lK)

!        Prestress
         IF (pstEq) THEN
            DO a=1, eNoN
               Ac = ptr(a)
               pSn(:,Ac) = pSn(:,Ac) + w*N(a)*pSl(:)
               pSa(Ac)   = pSa(Ac)   + w*N(a)
            END DO
         END IF
      END DO

!     Now doing the assembly part
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
      END SUBROUTINE CMMi
!--------------------------------------------------------------------
!     CMM initialization (boundary elements)
      SUBROUTINE CMMb(lFa, e, eNoN, al, dl, xl, bfl, pS0l, vwp, ptr)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER(KIND=IKIND), INTENT(IN) :: e, eNoN, ptr(eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: al(tDof,eNoN), dl(tDof,eNoN),
     2   xl(3,eNoN), bfl(3,eNoN), pS0l(6), vwp(2)

      INTEGER(KIND=IKIND) g
      REAL(KIND=RKIND) :: w, Jac, nV(3), pSl(6)

      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      ALLOCATE(N(eNoN), Nx(2,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
      lK  = 0._RKIND
      lR  = 0._RKIND

!     Internal stresses (stiffness) contribution
      pSl = 0._RKIND
      Nx  = lFa%Nx(:,:,1)
      CALL CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, vwp, pSl, lR, lK)

!     Inertia and body forces (mass) contribution
      DO g=1, lFa%nG
         CALL GNNB(lFa, e, g, nV)
         Jac = SQRT(NORM(nV))
         nV  = nV / Jac
         w   = lFa%w(g)*Jac
         N   = lFa%N(:,g)
         CALL CMM_MASS(eNoN, w, N, al, bfl, vwp, lR, lK)
      END DO

!     Now doing the assembly part
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
      END SUBROUTINE CMMb
!####################################################################
      SUBROUTINE CMM3D (eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), bfl(3,eNoN), ksix(3,3)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      REAL(KIND=RKIND), PARAMETER :: ct(2) = (/1._RKIND, 36._RKIND/)

      INTEGER(KIND=IKIND) a, b
      REAL(KIND=RKIND) tauM, tauC, tauB, kS, kU, nu, rho, T1, T2, T3,
     2  divU, amd, wl, wr, wrl, s, p, u(3), ud(3), px(3), f(3), up(3),
     3  ua(3), ux(3,3), udNx(eNoN), updNx(eNoN), uadNx(eNoN), rV(3),
     4  rM(3,3), NxdNx, nu_s, es(3,3), gam, nu_x, es_x(3,eNoN),
     5  bdry(eNoN)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      wr   = w*rho
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wrl  = wr*T1
      wl   = w*T1

      bdry = 0._RKIND
c      DO a=1, eNoN
c         bdry(a) = REAL(cmmBdry(ptr(a)), KIND=RKIND)
c      END DO

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0._RKIND
      u  = 0._RKIND
      ud = -f
      px = 0._RKIND
      ux = 0._RKIND
      DO a=1, eNoN
         p  = p + N(a)*yl(4,a)

         ud(1) = ud(1) + N(a)*(al(1,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(2,a)-bfl(2,a))
         ud(3) = ud(3) + N(a)*(al(3,a)-bfl(3,a))

         px(1) = px(1) + Nx(1,a)*yl(4,a)
         px(2) = px(2) + Nx(2,a)*yl(4,a)
         px(3) = px(3) + Nx(3,a)*yl(4,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)

c         u(1) = u(1) + N(a)*yl(1,a)*(1._RKIND-bdry(a))
c         u(2) = u(2) + N(a)*yl(2,a)*(1._RKIND-bdry(a))
c         u(3) = u(3) + N(a)*yl(3,a)*(1._RKIND-bdry(a))

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*yl(2,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*yl(3,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*yl(3,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*yl(3,a)
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

      DO a=1, eNoN
        es_x(1,a) = es(1,1)*Nx(1,a) + es(1,2)*Nx(2,a) + es(1,3)*Nx(3,a)
        es_x(2,a) = es(2,1)*Nx(1,a) + es(2,2)*Nx(2,a) + es(2,3)*Nx(3,a)
        es_x(3,a) = es(3,1)*Nx(1,a) + es(3,2)*Nx(2,a) + es(3,3)*Nx(3,a)
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
     2   ct(2)*nu_s*nu_s*kS + s*s)

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + u(3)*ux(3,1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + u(3)*ux(3,2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3)
     2      + u(3)*ux(3,3) + s*u(3))

      tauC = ksix(1,1) + ksix(2,2) + ksix(3,3)
      tauC = 1._RKIND/tauM/tauC

      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1)
     2     + up(3)*up(1)*ksix(3,1) + up(1)*up(2)*ksix(1,2)
     3     + up(2)*up(2)*ksix(2,2) + up(3)*up(2)*ksix(3,2)
     4     + up(1)*up(3)*ksix(1,3) + up(2)*up(3)*ksix(2,3)
     5     + up(3)*up(3)*ksix(3,3)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1._RKIND/SQRT(tauB)

      divU = ux(1,1) + ux(2,2) + ux(3,3)

      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      ua(3) = u(3) + up(3)

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

      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)  + u(3)*Nx(3,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a) + up(3)*Nx(3,a)
         uadNx(a) = updNx(a) + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*N(a) + rM(1,1)*Nx(1,a)
     2      + rM(1,2)*Nx(2,a) + rM(1,3)*Nx(3,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*N(a) + rM(2,1)*Nx(1,a)
     2      + rM(2,2)*Nx(2,a) + rM(2,3)*Nx(3,a))

         lR(3,a) = lR(3,a) + wr*(rV(3)*N(a) + rM(3,1)*Nx(1,a)
     2      + rM(3,2)*Nx(2,a) + rM(3,3)*Nx(3,a))

         lR(4,a) = lR(4,a) + w*(N(a)*divU - updNx(a))
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(3,1) = Nx(3,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)
            rM(3,2) = Nx(3,a)*Nx(2,b)
            rM(1,3) = Nx(1,a)*Nx(3,b)
            rM(2,3) = Nx(2,a)*Nx(3,b)
            rM(3,3) = Nx(3,a)*Nx(3,b)

            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + N(a)*((s+amd)*N(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*N(b))

            T2 = tauM*udNx(a)
            T3 = tauM*(amd*N(b) + udNx(b))
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

!           dM/dP
            lK(4,a,b)  = lK(4,a,b)  - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2)
            lK(8,a,b)  = lK(8,a,b)  - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)
            lK(12,a,b) = lK(12,a,b) - wl*(Nx(3,a)*N(b) - Nx(3,b)*T2)
!           dC/dU
            lK(13,a,b) = lK(13,a,b) + wl*(Nx(1,b)*N(a) + Nx(1,a)*T3)
            lK(14,a,b) = lK(14,a,b) + wl*(Nx(2,b)*N(a) + Nx(2,a)*T3)
            lK(15,a,b) = lK(15,a,b) + wl*(Nx(3,b)*N(a) + Nx(3,a)*T3)
!           dC/dP
            lK(16,a,b) = lK(16,a,b) + wl*(tauM*NxdNx)/rho
         END DO
      END DO

      RETURN
      END SUBROUTINE CMM3D
!####################################################################
!     Contribution from internal stresses (stiffness)
      SUBROUTINE CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, vwp, pSl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: Nx(2,eNoN), xl(3,eNoN),
     2   dl(tDof,eNoN), pS0l(6), vwp(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: pSl(6), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      REAL(KIND=RKIND), PARAMETER :: kT = 5._RKIND/6._RKIND

      INTEGER(KIND=IKIND) :: a, b, i, j, k
      REAL(KIND=RKIND) :: elM, nu, ht, lam, mu, afl, Jac, Jaci, T1,
     2   nV(3), gC(3,2), gCi(3,2), thet(3,3), thetT(3,3), phi(9,9),
     3   xloc(3,3), ul(9), Bm(5,9), BmT(9,5), Dm(5,5), DBm(5,9),
     4   pSm(3,3), S0l(5), Sl(5), BtS(9), Ke(9,9)

      nu  = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      IF (cmmVarWall) THEN
         ht  = vwp(1)
         elM = vwp(2)
      ELSE
         ht  = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
         elM = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      END IF

      lam = elM/(1._RKIND - nu*nu)
      mu  = 0.5_RKIND*elM/(1._RKIND + nu)

      afl = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i   = eq(cEq)%s
      j   = i + 1
      k   = j + 1

      CALL GNNS(eNoN, Nx, xl, nV, gC, gCi)
      Jac   = SQRT(NORM(nV))
      Jaci  = 1._RKIND/Jac
      nV(:) = nV(:) * Jaci

!     Rotation matrix
      thet(1,:) = gC(:,1)/SQRT(NORM(gC(:,1)))
      thet(3,:) = nV(:)
      thet(2,1) = thet(3,2)*thet(1,3) - thet(3,3)*thet(1,2)
      thet(2,2) = thet(3,3)*thet(1,1) - thet(3,1)*thet(1,3)
      thet(2,3) = thet(3,1)*thet(1,2) - thet(3,2)*thet(1,1)
      thetT = TRANSPOSE(thet)

!     Define phi matrix
      phi = 0._RKIND
      phi(1:3,1:3) = thet
      phi(4:6,4:6) = thet
      phi(7:9,7:9) = thet

!     Transform the global coordinates into local frame. Copy
!     displacements into a vector.
      DO a=1, eNoN
         xloc(1,a) = thet(1,1)*xl(1,a) + thet(1,2)*xl(2,a) +
     2      thet(1,3)*xl(3,a)
         xloc(2,a) = thet(2,1)*xl(1,a) + thet(2,2)*xl(2,a) +
     2      thet(2,3)*xl(3,a)
         xloc(3,a) = thet(3,1)*xl(1,a) + thet(3,2)*xl(2,a) +
     2      thet(3,3)*xl(3,a)

         b = (a-1)*eNoN
         ul(b+1) = dl(i,a)
         ul(b+2) = dl(j,a)
         ul(b+3) = dl(k,a)
      END DO

!     If prestress is present, convert into full matrix, and transform
!     into local coordinates, convert back into voigt notation
      pSm = 0._RKIND
      pSm(1,1) = pS0l(1)
      pSm(2,2) = pS0l(2)
      pSm(3,3) = pS0l(3)
      pSm(1,2) = pS0l(4)
      pSm(1,3) = pS0l(5)
      pSm(2,3) = pS0l(6)

      pSm(2,1) = pSm(1,2)
      pSm(3,1) = pSm(1,3)
      pSm(3,2) = pSm(2,3)

      pSm = MATMUL(pSm, thetT)
      pSm = MATMUL(thet, pSm)

      S0l(1) = pSm(1,1)
      S0l(2) = pSm(2,2)
      S0l(3) = pSm(1,2)
      S0l(4) = pSm(1,3)
      S0l(5) = pSm(2,3)

!     B matrix
      Bm = 0._RKIND
      Bm(1,1) = (xloc(2,2) - xloc(2,3))*Jaci  ! y23 = y2 - y3
      Bm(1,4) = (xloc(2,3) - xloc(2,1))*Jaci  ! y31 = y3 - y1
      Bm(1,7) = (xloc(2,1) - xloc(2,2))*Jaci  ! y12 = y1 - y2

      Bm(2,2) = (xloc(1,3) - xloc(1,2))*Jaci  ! x32 = x3 - x2
      Bm(2,5) = (xloc(1,1) - xloc(1,3))*Jaci  ! x13 = x1 - x3
      Bm(2,8) = (xloc(1,2) - xloc(1,1))*Jaci  ! x21 = x2 - x1

      Bm(3,1) = Bm(2,2)
      Bm(3,2) = Bm(1,1)
      Bm(3,4) = Bm(2,5)
      Bm(3,5) = Bm(1,4)
      Bm(3,7) = Bm(2,8)
      Bm(3,8) = Bm(1,7)

      Bm(4,3) = Bm(1,1)
      Bm(4,6) = Bm(1,4)
      Bm(4,9) = Bm(1,7)

      Bm(5,3) = Bm(2,2)
      Bm(5,6) = Bm(2,5)
      Bm(5,9) = Bm(2,8)

!     Transform B using phi and compute its transpose
      Bm  = MATMUL(Bm, phi)
      BmT = TRANSPOSE(Bm)

!     Material tensor, D
      Dm  = 0._RKIND
      Dm(1,1) = lam
      Dm(1,2) = lam*nu
      Dm(3,3) = mu
      Dm(4,4) = mu*kT

      Dm(2,1) = Dm(1,2)
      Dm(2,2) = Dm(1,1)
      Dm(5,5) = Dm(4,4)

!     D*Bm
      DBm = MATMUL(Dm, Bm)

!     Stress tensor in local frame
      Sl = 0._RKIND
      DO a=1, 9
         Sl(1) = Sl(1) + DBm(1,a)*ul(a)
         Sl(2) = Sl(2) + DBm(2,a)*ul(a)
         Sl(3) = Sl(3) + DBm(3,a)*ul(a)
         Sl(4) = Sl(4) + DBm(4,a)*ul(a)
         Sl(5) = Sl(5) + DBm(5,a)*ul(a)
      END DO

!     Prestress is updated with new stress and pulled to global frame
      pSl = 0._RKIND
      IF (pstEq) THEN
         pSm(:,:) = 0._RKIND
         pSm(1,1) = Sl(1)
         pSm(2,2) = Sl(2)
         pSm(1,2) = Sl(3)
         pSm(1,3) = Sl(4)
         pSm(2,3) = Sl(5)

         pSm(2,1) = pSm(1,2)
         pSm(3,1) = pSm(1,3)
         pSm(3,2) = pSm(2,3)

         pSm = MATMUL(pSm,  thet)
         pSm = MATMUL(thetT, pSm)

         pSl(1) = pSm(1,1)
         pSl(2) = pSm(2,2)
         pSl(3) = pSm(3,3)
         pSl(4) = pSm(1,2)
         pSl(5) = pSm(1,3)
         pSl(6) = pSm(2,3)
      END IF

!     Internal stress contribution to residue together with prestress
      Sl(:) = Sl(:) + S0l(:)
      BtS = 0._RKIND
      DO a=1, 5
         BtS(1) = BtS(1) + BmT(1,a)*Sl(a)
         BtS(2) = BtS(2) + BmT(2,a)*Sl(a)
         BtS(3) = BtS(3) + BmT(3,a)*Sl(a)
         BtS(4) = BtS(4) + BmT(4,a)*Sl(a)
         BtS(5) = BtS(5) + BmT(5,a)*Sl(a)
         BtS(6) = BtS(6) + BmT(6,a)*Sl(a)
         BtS(7) = BtS(7) + BmT(7,a)*Sl(a)
         BtS(8) = BtS(8) + BmT(8,a)*Sl(a)
         BtS(9) = BtS(9) + BmT(9,a)*Sl(a)
      END DO

!     Now all the required tensors are defined, contributions to
!     residue and stiffness can be computed
      T1 = ht*Jac*0.5_RKIND
      DO a=1, eNoN
         b = (a-1)*eNoN
         lR(1,a) = lR(1,a) + T1*BtS(b+1)
         lR(2,a) = lR(2,a) + T1*BtS(b+2)
         lR(3,a) = lR(3,a) + T1*BtS(b+3)
      END DO

!     Compute element level global stiffness matrix and remapping
      Ke = MATMUL(BmT, DBm)
      T1 = T1*afl
      DO a=1, eNoN
         i = (a-1)*eNoN
         DO b=1, eNoN
            j = (b-1)*eNoN

            lK(1,a,b) = lK(1,a,b) + T1*Ke(i+1,j+1)
            lK(2,a,b) = lK(2,a,b) + T1*Ke(i+1,j+2)
            lK(3,a,b) = lK(3,a,b) + T1*Ke(i+1,j+3)

            lK(dof+1,a,b) = lK(dof+1,a,b) + T1*Ke(i+2,j+1)
            lK(dof+2,a,b) = lK(dof+2,a,b) + T1*Ke(i+2,j+2)
            lK(dof+3,a,b) = lK(dof+3,a,b) + T1*Ke(i+2,j+3)

            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + T1*Ke(i+3,j+1)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + T1*Ke(i+3,j+2)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + T1*Ke(i+3,j+3)
         END DO
      END DO

      RETURN
      END SUBROUTINE CMM_STIFFNESS
!--------------------------------------------------------------------
!     Contribution from inertia and body forces
      SUBROUTINE CMM_MASS(eNoN, w, N, al, bfl, vwp, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), al(tDof,eNoN),
     2   bfl(3,eNoN), vwp(2)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i, j, k
      REAL(KIND=RKIND) rho, ht, wl, am, T1, f(3), ud(3)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      IF (cmmVarWall) THEN
         ht = vwp(1)
      ELSE
         ht = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      END IF

      wl = w * ht * rho
      am = eq(cEq)%am
      i  = eq(cEq)%s
      j  = i + 1
      k  = j + 1

      ud = -f
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(al(i,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(j,a)-bfl(2,a))
         ud(3) = ud(3) + N(a)*(al(k,a)-bfl(3,a))
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + wl*N(a)*ud(1)
         lR(2,a) = lR(2,a) + wl*N(a)*ud(2)
         lR(3,a) = lR(3,a) + wl*N(a)*ud(3)
      END DO

      DO b=1, eNoN
         DO a=1, eNoN
            T1 = wl*am*N(a)*N(b)
            lK(1,a,b) = lK(1,a,b) + T1
            lK(dof+2,a,b) = lK(dof+2,a,b) + T1
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + T1
         END DO
      END DO

      RETURN
      END SUBROUTINE CMM_MASS
!####################################################################
!     Set traction on CMM wall during initialization. The load
!     applied is not a usual follower pressure load, aka, load does not
!     follow deformation (acceptable for small strains - CMM)
      SUBROUTINE BCMMi(eNoN, idof, w, N, Nx, xl, tfl, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, idof
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   xl(3,eNoN), tfl(idof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER(KIND=IKIND) a
      REAL(KIND=RKIND) :: wl, nV(3), gC(3,2)

      REAL(KIND=RKIND), ALLOCATABLE :: tfn(:)

!     Get traction vector
      ALLOCATE(tfn(idof))
      tfn = 0._RKIND
      DO a=1, eNoN
         tfn = tfn + N(a)*tfl(:,a)
      END DO

!     Get surface normal vector (reference configuration)
      CALL GNNS(eNoN, Nx, xl, nV, gC, gC)

!     Local residue
      IF (idof .EQ. 1) THEN
         wl = w * tfn(1)
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - wl*N(a)*nV(1)
            lR(2,a) = lR(2,a) - wl*N(a)*nV(2)
            lR(3,a) = lR(3,a) - wl*N(a)*nV(3)
         END DO
      ELSE
         wl = w * SQRT(NORM(nV(:)))
         DO a=1, eNoN
            lR(1,a) = lR(1,a) - wl*N(a)*tfn(1)
            lR(2,a) = lR(2,a) - wl*N(a)*tfn(2)
            lR(3,a) = lR(3,a) - wl*N(a)*tfn(3)
         END DO
      END IF

      DEALLOCATE(tfn)

      RETURN
      END SUBROUTINE BCMMi
!####################################################################

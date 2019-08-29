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
!     strongly, weakly or using penalty methods for immersed boundaries.
!
!--------------------------------------------------------------------

      SUBROUTINE FLUID3D (eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), al(tDof,eNoN),
     2   yl(tDof,eNoN), bfl(3,eNoN), ksix(3,3)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,36D0/)
      INTEGER a, b
      REAL(KIND=8) tauM, tauC, tauB, kS, kU, nu, rho, T1, T2, T3, divU,
     2  amd, wl, wr, wrl, s, p, u(3), ud(3), px(3), f(3), up(3), ua(3),
     3  ux(3,3), udNx(eNoN), updNx(eNoN), uadNx(eNoN), rV(3), rM(3,3),
     4  NxdNx, nu_s, es(3,3), gam, nu_x, es_x(3,eNoN)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      wr   = w*rho
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wrl  = wr*T1
      wl   = w*T1

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
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
      IF (mvMsh) THEN
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(5,a)
            u(2) = u(2) - N(a)*yl(6,a)
            u(3) = u(3) - N(a)*yl(7,a)
         END DO
      END IF

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
      gam = SQRT(0.5D0*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0D0
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

      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU +
     2   ct(2)*nu_s*nu_s*kS + s*s)

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + u(3)*ux(3,1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + u(3)*ux(3,2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3)
     2      + u(3)*ux(3,3) + s*u(3))

      tauC = ksix(1,1) + ksix(2,2) + ksix(3,3)
      tauC = 1D0/tauM/tauC

      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1)
     2     + up(3)*up(1)*ksix(3,1) + up(1)*up(2)*ksix(1,2)
     3     + up(2)*up(2)*ksix(2,2) + up(3)*up(2)*ksix(3,2)
     4     + up(1)*up(3)*ksix(1,3) + up(2)*up(3)*ksix(2,3)
     5     + up(3)*up(3)*ksix(3,3)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1D0/SQRT(tauB)

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
      END SUBROUTINE FLUID3D
!--------------------------------------------------------------------
      SUBROUTINE FLUID2D (eNoN, w, N, Nx, al, yl, bfl, ksix, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), al(tDof,eNoN),
     2   yl(tDof,eNoN), bfl(2,eNoN), ksix(2,2)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,36D0/)
      INTEGER a, b
      REAL(KIND=8) tauM, tauC, tauB, kS, kU, nu, rho, T1, T2, T3, divU,
     2  amd, wl, wr, wrl, s, p, u(2), ud(2), px(2), f(2), up(2), ua(2),
     3  ux(2,2), udNx(eNoN), updNx(eNoN), uadNx(eNoN), rV(2), rM(2,2),
     4  NxdNx, nu_s, es(2,2), gam, nu_x, es_x(2,eNoN)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)

      wr   = w*rho
      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wrl  = wr*T1
      wl   = w*T1

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         p = p + N(a)*yl(3,a)

         ud(1) = ud(1) + N(a)*(al(1,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(2,a)-bfl(2,a))

         px(1) = px(1) + Nx(1,a)*yl(3,a)
         px(2) = px(2) + Nx(2,a)*yl(3,a)

         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
      END DO
      IF (mvMsh) THEN
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(4,a)
            u(2) = u(2) - N(a)*yl(5,a)
         END DO
      END IF

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

      DO a=1, eNoN
        es_x(1,a) = es(1,1)*Nx(1,a) + es(1,2)*Nx(2,a)
        es_x(2,a) = es(2,1)*Nx(1,a) + es(2,2)*Nx(2,a)
      END DO

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5D0*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(cEq)%dmn(cDmn), gam, nu, nu_s, nu_x)
      nu   = nu/rho
      nu_s = nu_s/rho
      IF (ISZERO(gam)) THEN
         nu_x = 0D0
      ELSE
         nu_x = nu_x/rho/gam
      END IF
      s  = nu/eq(cEq)%dmn(cDmn)%prop(permeability)

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)

      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU +
     2   ct(2)*nu_s*nu_s*kS + s*s)

      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + s*u(2))

      tauC = ksix(1,1) + ksix(2,2)
      tauC = 1D0/tauM/tauC

      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1)
     2     + up(1)*up(2)*ksix(1,2) + up(2)*up(2)*ksix(2,2)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1D0/SQRT(tauB)

      divU = ux(1,1) + ux(2,2)

      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)

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

      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a)
         uadNx(a) = updNx(a) + udNx(a)

         lR(1,a) = lR(1,a) + wr*(rV(1)*N(a) + rM(1,1)*Nx(1,a)
     2      + rM(1,2)*Nx(2,a))

         lR(2,a) = lR(2,a) + wr*(rV(2)*N(a) + rM(2,1)*Nx(1,a)
     2      + rM(2,2)*Nx(2,a))

         lR(3,a) = lR(3,a) + w*(N(a)*divU - updNx(a))
      END DO

      DO a=1, eNoN
         DO b=1, eNoN
            rM(1,1) = Nx(1,a)*Nx(1,b)
            rM(2,1) = Nx(2,a)*Nx(1,b)
            rM(1,2) = Nx(1,a)*Nx(2,b)
            rM(2,2) = Nx(2,a)*Nx(2,b)

            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T1 = nu*NxdNx + tauB*updNx(a)*updNx(b)
     2         + N(a)*((s+amd)*N(b) + uadNx(b))
     3         + tauM*uadNx(a)*(udNx(b) + (s+amd)*N(b))

            T2 = tauM*udNx(a)
            T3 = tauM*(amd*N(b) + udNx(b))

!           dM/dU
            lK(1,a,b)  = lK(1,a,b) + wrl*((nu + tauC)*rM(1,1) + T1 +
     2         nu_x*es_x(1,a)*es_x(1,b))
            lK(2,a,b)  = lK(2,a,b) + wrl*(nu*rM(2,1) + tauC*rM(1,2) +
     2         nu_x*es_x(1,a)*es_x(2,b))

            lK(4,a,b)  = lK(4,a,b) + wrl*(nu*rM(1,2) + tauC*rM(2,1) +
     2         nu_x*es_x(2,a)*es_x(1,b))
            lK(5,a,b)  = lK(5,a,b) + wrl*((nu + tauC)*rM(2,2) + T1 +
     2         nu_x*es_x(2,a)*es_x(2,b))

!           dM/dP
            lK(3,a,b)  = lK(3,a,b) - wl*(Nx(1,a)*N(b) - Nx(1,b)*T2)
            lK(6,a,b)  = lK(6,a,b) - wl*(Nx(2,a)*N(b) - Nx(2,b)*T2)

!           dC/dU
            lK(7,a,b) = lK(7,a,b)  + wl*(Nx(1,b)*N(a) + Nx(1,a)*T3)
            lK(8,a,b) = lK(8,a,b)  + wl*(Nx(2,b)*N(a) + Nx(2,a)*T3)

!           dC/dP
            lK(9,a,b) = lK(9,a,b)  + wl*tauM*NxdNx/rho
         END DO
      END DO

      RETURN
      END SUBROUTINE FLUID2D
!####################################################################
      PURE SUBROUTINE BFLUID (eNoN, w, N, y, h, nV, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), y(tDof), h, nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j
      REAL(KIND=8) T1, wl, hc(nsd), udn, u(nsd)

      wl  = w*eq(cEq)%af*eq(cEq)%gam*dt
      udn = 0D0
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

      udn = eq(cEq)%dmn(cDmn)%prop(backflow_stab)/2D0*
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
      SUBROUTINE BWFLUID3D (eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), yl(tDof,eNoN),
     2   ub(3), nV(3), tauB(2)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: a, b
      REAL(KIND=8) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn, ubn,
     2   u(3), uh(3), ux(3,3), sigman(3), Nxn(eNoN), rV(3), rM(3,3),
     3   nu_s, es(3,3), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w * rho
      wl   = w * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0D0
      u    = 0D0
      ux   = 0D0
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

      uh = 0D0
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(5,a)
            uh(2) = uh(2) + N(a)*yl(6,a)
            uh(3) = uh(3) + N(a)*yl(7,a)
         END DO
      END IF
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2) + (u(3)-ub(3))*nV(3)
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2) + (u(3)-uh(3))*nV(3)
      uhn = (ABS(uhn) - uhn) / 2D0

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
      gam = SQRT(0.5D0*gam)

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
      SUBROUTINE BWFLUID2D (eNoN, w, N, Nx, yl, ub, nV, tauB, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), yl(tDof,eNoN),
     2   ub(2), nV(2), tauB(2)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER :: a, b
      REAL(KIND=8) :: rho, nu, T1, wr, wl, wrl, tauT, tauN, p, uhn, ubn,
     2   u(2), uh(2), ux(2,2), sigman(2), Nxn(eNoN), rV(2), rM(2,2),
     3   nu_s, es(2,2), gam, nu_x

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      T1   = eq(cEq)%af * eq(cEq)%gam * dt
      wr   = w*rho
      wl   = w  * T1
      wrl  = wr * T1
      tauT = tauB(1) / rho
      tauN = tauB(2) / rho

      p    = 0D0
      u    = 0D0
      ux   = 0D0
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

      uh = 0D0
      IF (mvMsh) THEN
         DO a=1, eNoN
            uh(1) = uh(1) + N(a)*yl(4,a)
            uh(2) = uh(2) + N(a)*yl(5,a)
         END DO
      END IF
      uhn = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2)
      ubn = (u(1)-ub(1))*nV(1) + (u(2)-ub(2))*nV(2)
      uhn = (ABS(uhn) - uhn) / 2D0

!     Strain rate tensor 2*e_ij := (u_ij + u_ji)
      es(1,1) = ux(1,1) + ux(1,1)
      es(2,2) = ux(2,2) + ux(2,2)

      es(1,2) = ux(1,2) + ux(2,1)
      es(2,1) = es(1,2)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(2,1)*es(2,1)
     2    + es(2,2)*es(2,2)
      gam = SQRT(0.5D0*gam)

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
      REAL(KIND=8), INTENT(INOUT)  :: gamma
      REAL(KIND=8), INTENT(OUT) :: mu, mu_s, mu_x

      REAL(KIND=8) :: mu_i, mu_o, lam, a, n, T1, T2

      SELECT CASE (lDmn%visc%viscType)
      CASE (viscType_Const)
         mu   = lDmn%visc%mu_i
         mu_s = mu
         mu_x = 0D0

      CASE (viscType_CY)
         mu_i = lDmn%visc%mu_i
         mu_o = lDmn%visc%mu_o
         lam  = lDmn%visc%lam
         a    = lDmn%visc%a
         n    = lDmn%visc%n

         T1   = 1.0D0 + (lam*gamma)**a
         T2   = T1**((n-1.0D0)/a)
         mu   = mu_i + (mu_o-mu_i)*T2
         mu_s = mu_i

         T1   = T2/T1
         T2   = lam**a * gamma**(a-1.0D0) * T1
         mu_x = (mu_o-mu_i)*(n-1.0D0)*T2

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
         mu_x = 2.0D0*mu_o*(mu_o + mu_i)/gamma

      END SELECT

      RETURN
      END SUBROUTINE GETVISCOSITY
!####################################################################

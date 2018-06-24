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
!      This routine embodies formulation for solving
!      Basset–Boussinesq–Oseen equation, which models the particle
!      motion with gravity and inertial effects.
!      
!--------------------------------------------------------------------

      PURE SUBROUTINE BBO3D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)
      
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,1D0,4D0/)

      INTEGER a, b, i, j, k
      REAL(KIND=8) wr, T1, T2, amd, wl, udNx(eNoN), res(nsd), tauM, kU, 
     2   nu, tauB, up(nsd), updNx(eNoN), tmpV(nsd), NxdNx, rM(nsd,nsd), 
     3   tauMs, wlnu, ud(nsd), u(nsd), ux(nsd,nsd), bf(nsd), s, px(nsd),
     4   rho, pRho, pDia

      pDia  = eq(cEq)%dmn(cDmn)%prop(particle_diameter)
      pRho  = eq(cEq)%dmn(cDmn)%prop(particle_density)
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      nu    = eq(cEq)%dmn(cDmn)%prop(viscosity)
      bf(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      bf(3) = eq(cEq)%dmn(cDmn)%prop(f_z) 
      T1    = eq(cEq)%af*eq(cEq)%gam*dt
      amd   = eq(cEq)%am/T1
      i     = eq(cEq)%s
      j     = i + 1
      k     = j + 1
      s     = 18D0*nu/pDia/pDia/(pRho + 5D-1*rho)
      wr    = w*(pRho + 5D-1*rho)
      wl    = wr*T1

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      ud = 0D0
      u  = 0D0
      up = 0D0
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         up(1) = up(1) + N(a)*yl(i,a)
         up(2) = up(2) + N(a)*yl(j,a)
         up(3) = up(3) + N(a)*yl(k,a)
               
         ud(1) = ud(1) + N(a)*al(1,a)
         ud(2) = ud(2) + N(a)*al(2,a)
         ud(3) = ud(3) + N(a)*al(3,a)
               
         u(1)  = u(1)  + N(a)*yl(1,a)
         u(2)  = u(2)  + N(a)*yl(2,a)
         u(3)  = u(3)  + N(a)*yl(3,a)
               
         px(1) = px(1) + Nx(1,a)*yl(3,a)
         px(2) = px(2) + Nx(2,a)*yl(3,a)
         px(3) = px(3) + Nx(3,a)*yl(3,a)

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
 
      bf(1) = s*u(1) + (-rho*px(1) + 5D-1*rho*(ud(1)
     2   + up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
     3   + (pRho - rho)*bf(1))/(pRho + 5D-1*rho)
      bf(2) = s*u(2) + (-rho*px(2) + 5D-1*rho*(ud(2)
     2   + up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
     3   + (pRho - rho)*bf(2))/(pRho + 5D-1*rho)
      bf(3) = s*u(3) + (-rho*px(3) + 5D-1*rho*(ud(3)
     2   + up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))
     3   + (pRho - rho)*bf(3))/(pRho + 5D-1*rho)
       
      u  = up
      ud = 0D0
      ux = 0D0
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*al(i,a)
         ud(2) = ud(2) + N(a)*al(j,a)
         ud(3) = ud(3) + N(a)*al(k,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(i,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(i,a)
         ux(3,1) = ux(3,1) + Nx(3,a)*yl(i,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(j,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(j,a)
         ux(3,2) = ux(3,2) + Nx(3,a)*yl(j,a)
         ux(1,3) = ux(1,3) + Nx(1,a)*yl(k,a)
         ux(2,3) = ux(2,3) + Nx(2,a)*yl(k,a)
         ux(3,3) = ux(3,3) + Nx(3,a)*yl(k,a)
      END DO

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1) 
     2   + u(3)*u(1)*ksix(3,1) + u(1)*u(2)*ksix(1,2) 
     3   + u(2)*u(2)*ksix(2,2) + u(3)*u(2)*ksix(3,2) 
     4   + u(1)*u(3)*ksix(1,3) + u(2)*u(3)*ksix(2,3) 
     5   + u(3)*u(3)*ksix(3,3)

      tauM = ct(3)/SQRT(ct(1)/dt/dt + ct(2)*kU + s*s)
      tauMs = ct(3)/SQRT(ct(2)*kU + s*s)
      
      res(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + u(3)*ux(3,1)
     2   + s*u(1) - bf(1)
      res(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + u(3)*ux(3,2) 
     2   + s*u(2) - bf(2)
      res(3) = ud(3) + u(1)*ux(1,3) + u(2)*ux(2,3) + u(3)*ux(3,3) 
     2   + s*u(3) - bf(3)

!     Discontinuity capturing nu (wx , ux)
      tmpV(1) = (u(1)*ux(1,1) + u(2)*ux(1,2) + u(3)*ux(1,3))
      tmpV(2) = (u(1)*ux(2,1) + u(2)*ux(2,2) + u(3)*ux(2,3))
      tmpV(2) = (u(1)*ux(3,1) + u(2)*ux(3,2) + u(3)*ux(3,3))
      
      nu = tmpV(1)*tmpV(1)*ksix(1,1) + tmpV(2)*tmpV(1)*ksix(2,1) 
     2   + tmpV(3)*tmpV(1)*ksix(3,1) + tmpV(1)*tmpV(2)*ksix(1,2) 
     3   + tmpV(2)*tmpV(2)*ksix(2,2) + tmpV(3)*tmpV(2)*ksix(3,2) 
     4   + tmpV(1)*tmpV(3)*ksix(1,3) + tmpV(2)*tmpV(3)*ksix(2,3) 
     5   + tmpV(3)*tmpV(3)*ksix(3,3)
      IF (ISZERO(nu)) nu = eps
      nu = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))/nu
      nu = SQRT(nu*(res(1)*res(1) + res(2)*res(2) + res(3)*res(3)))/2D0
      IF (ISZERO(kU)) kU = eps
      T1 = (u(1)*u(1) + u(2)*u(2) + u(3)*u(3))/SQRT(kU)/2D0
      IF (nu .GT. T1) nu = T1

!     tauB terms tauB (up.wx , up.ux)
      up = -tauM*res
      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1) 
     2     + up(3)*up(1)*ksix(3,1) + up(1)*up(2)*ksix(1,2) 
     3     + up(2)*up(2)*ksix(2,2) + up(3)*up(2)*ksix(3,2) 
     4     + up(1)*up(3)*ksix(1,3) + up(2)*up(3)*ksix(2,3) 
     5     + up(3)*up(3)*ksix(3,3)
      IF (ISZERO(tauB)) tauB = eps
      tauB = ct(3)/SQRT(tauB)
      
      rM(1,1) = nu*(ux(1,1) + ux(1,1)) 
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) 
      rM(3,1) = nu*(ux(3,1) + ux(1,3)) 
      rM(1,2) = nu*(ux(1,2) + ux(2,1))
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) 
      rM(3,2) = nu*(ux(3,2) + ux(2,3)) 
      rM(1,3) = nu*(ux(1,3) + ux(3,1))
      rM(2,3) = nu*(ux(2,3) + ux(3,2)) 
      rM(3,3) = nu*(ux(3,3) + ux(3,3)) 

      tmpV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1) + up(3)*ux(3,1))
      tmpV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2) + up(3)*ux(3,2))
      tmpV(3) = tauB*(up(1)*ux(1,3) + up(2)*ux(2,3) + up(3)*ux(3,3))
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a)  + u(2)*Nx(2,a)  + u(3)*Nx(3,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a) + up(3)*Nx(3,a)
         T1       = N(a) + udNx(a)*tauM

         lR(1,a) = lR(1,a) + wr*(res(1)*T1 + updNx(a)*tmpV(1) 
     2           + rM(1,1)*Nx(1,a) + rM(1,2)*Nx(2,a) + rM(1,3)*Nx(3,a))
         lR(2,a) = lR(2,a) + wr*(res(2)*T1 + updNx(a)*tmpV(2)
     2           + rM(2,1)*Nx(1,a) + rM(2,2)*Nx(2,a) + rM(2,3)*Nx(3,a))
         lR(3,a) = lR(3,a) + wr*(res(3)*T1 + updNx(a)*tmpV(3)
     2           + rM(3,1)*Nx(1,a) + rM(3,2)*Nx(2,a) + rM(3,3)*Nx(3,a))
      END DO

      wlnu = wl*nu
      DO a=1, eNoN
         T1 = N(a) + udNx(a)*tauM
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T2 = wl*(T1*(udNx(b) + (amd + s)*N(b)) 
     2         + tauB*updNx(a)*updNx(b) + nu*NxdNx)
             
            lK(1,a,b) = lK(1,a,b) + wlnu*Nx(1,a)*Nx(1,b) + T2
            lK(2,a,b) = lK(2,a,b) + wlnu*Nx(2,a)*Nx(1,b)
            lK(3,a,b) = lK(3,a,b) + wlnu*Nx(3,a)*Nx(1,b)
            lK(4,a,b) = lK(4,a,b) + wlnu*Nx(1,a)*Nx(2,b)
            lK(5,a,b) = lK(5,a,b) + wlnu*Nx(2,a)*Nx(2,b) + T2
            lK(6,a,b) = lK(6,a,b) + wlnu*Nx(3,a)*Nx(2,b)
            lK(7,a,b) = lK(7,a,b) + wlnu*Nx(1,a)*Nx(3,b)
            lK(8,a,b) = lK(8,a,b) + wlnu*Nx(2,a)*Nx(3,b)
            lK(9,a,b) = lK(9,a,b) + wlnu*Nx(3,a)*Nx(3,b) + T2
         END DO
      END DO

      RETURN
      END SUBROUTINE BBO3D

!====================================================================
      
      PURE SUBROUTINE BBO2D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)
      
      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,1D0,4D0/)
      INTEGER a, b, i, j
      REAL(KIND=8) wr, T1, T2, amd, wl, udNx(eNoN), res(nsd), tauM, kU, 
     2   nu, tauB, up(nsd), updNx(eNoN), tmpV(nsd), NxdNx, rM(nsd,nsd), 
     3   tauMs, wlnu, ud(nsd), u(nsd), ux(nsd,nsd), bf(nsd), s, rho,
     4   pRho, pDia, px(nsd)

      pDia  = eq(cEq)%dmn(cDmn)%prop(particle_diameter)
      pRho  = eq(cEq)%dmn(cDmn)%prop(particle_density)
      rho   = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      nu    = eq(cEq)%dmn(cDmn)%prop(viscosity)
      bf(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      bf(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      T1    = eq(cEq)%af*eq(cEq)%gam*dt
      amd   = eq(cEq)%am/T1
      i     = eq(cEq)%s
      j     = i + 1
      wr    = w*(pRho + 5D-1*rho)
      wl    = wr*T1
      s     = 18D0*nu/pDia/pDia/(pRho + 5D-1*rho)
 
!     Indices are not selected based on the equation only
!     because fluid equation always come first
      ud = 0D0
      u  = 0D0
      up = 0D0
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         up(1) = up(1) + N(a)*yl(i,a)
         up(2) = up(2) + N(a)*yl(j,a)
               
         ud(1) = ud(1) + N(a)*al(1,a)
         ud(2) = ud(2) + N(a)*al(2,a)
               
         u(1)  = u(1)  + N(a)*yl(1,a)
         u(2)  = u(2)  + N(a)*yl(2,a)
               
         px(1) = px(1) + Nx(1,a)*yl(3,a)
         px(2) = px(2) + Nx(2,a)*yl(3,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(1,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(1,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(2,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(2,a)
      END DO

      bf(1) = s*u(1) + (-rho*px(1) + 5D-1*rho*(ud(1) 
     2   + up(1)*ux(1,1) + up(2)*ux(2,1))
     3   + (pRho - rho)*bf(1))/(pRho + 5D-1*rho)
      bf(2) = s*u(2) + (-rho*px(2) + 5D-1*rho*(ud(2) 
     2   + up(1)*ux(1,2) + up(2)*ux(2,2))
     3   + (pRho - rho)*bf(2))/(pRho + 5D-1*rho)

      u  = up
      ud = 0D0
      ux = 0D0
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*al(i,a)
         ud(2) = ud(2) + N(a)*al(j,a)

         ux(1,1) = ux(1,1) + Nx(1,a)*yl(i,a)
         ux(2,1) = ux(2,1) + Nx(2,a)*yl(i,a)
         ux(1,2) = ux(1,2) + Nx(1,a)*yl(j,a)
         ux(2,2) = ux(2,2) + Nx(2,a)*yl(j,a)
      END DO

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1) 
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      tauM = ct(3)/SQRT(ct(1)/dt/dt + ct(2)*kU + s*s)
      tauMs = ct(3)/SQRT(ct(2)*kU + s*s)
      
      res(1) = ud(1) + u(1)*ux(1,1) + u(2)*ux(2,1) + s*u(1) - bf(1)
      res(2) = ud(2) + u(1)*ux(1,2) + u(2)*ux(2,2) + s*u(2) - bf(2)

!     Discontinuity capturing nu (wx , ux)
      tmpV(1) = (u(1)*ux(1,1) + u(2)*ux(1,2))
      tmpV(2) = (u(1)*ux(2,1) + u(2)*ux(2,2))
      nu = tmpV(1)*tmpV(1)*ksix(1,1) + tmpV(2)*tmpV(1)*ksix(2,1) 
     2   + tmpV(1)*tmpV(2)*ksix(1,2) + tmpV(2)*tmpV(2)*ksix(2,2)
      IF (ISZERO(nu)) nu = eps
      nu = (u(1)*u(1) + u(2)*u(2))/nu
      nu = SQRT(nu*(res(1)*res(1) + res(2)*res(2)))/2D0
      IF (ISZERO(kU)) kU = eps
      T1 = (u(1)*u(1) + u(2)*u(2))/SQRT(kU)/2D0
      IF (nu .GT. T1) nu = T1

!     tauB terms tauB (up.wx , up.ux)
      up = -tauM*res
      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1) 
     2     + up(1)*up(2)*ksix(1,2) + up(2)*up(2)*ksix(2,2)
      IF (ISZERO(tauB)) tauB = eps
      tauB = ct(3)/SQRT(tauB)
      
      rM(1,1) = nu*(ux(1,1) + ux(1,1)) 
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) 
      rM(1,2) = nu*(ux(1,2) + ux(2,1))
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) 

      tmpV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      tmpV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))
      DO a=1, eNoN
         udNx(a)  = u(1)*Nx(1,a) + u(2)*Nx(2,a)
         updNx(a) = up(1)*Nx(1,a) + up(2)*Nx(2,a)
         T1       = N(a) + udNx(a)*tauM

         lR(1,a) = lR(1,a) + wr*(res(1)*T1 + updNx(a)*tmpV(1) 
     2           + rM(1,1)*Nx(1,a) + rM(1,2)*Nx(2,a))
         lR(2,a) = lR(2,a) + wr*(res(2)*T1 + updNx(a)*tmpV(2)
     2           + rM(2,1)*Nx(1,a) + rM(2,2)*Nx(2,a))
      END DO

      wlnu = wl*nu
      DO a=1, eNoN
         T1 = N(a) + udNx(a)*tauM
         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T2 = wl*(T1*(udNx(b) + (amd + s)*N(b)) 
     2         + tauB*updNx(a)*updNx(b) + nu*NxdNx)
             
            lK(1,a,b) = lK(1,a,b) + wlnu*Nx(1,a)*Nx(1,b) + T2
            lK(2,a,b) = lK(2,a,b) + wlnu*Nx(2,a)*Nx(1,b)
            lK(3,a,b) = lK(3,a,b) + wlnu*Nx(1,a)*Nx(2,b)
            lK(4,a,b) = lK(4,a,b) + wlnu*Nx(2,a)*Nx(2,b) + T2
         END DO
      END DO

      RETURN
      END SUBROUTINE BBO2D

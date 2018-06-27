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
!     Neccessary contribution to LHS/RHS are calculated for fluid NS
!     problem  in this routine.
!      
!--------------------------------------------------------------------

      PURE SUBROUTINE FLUID3D (eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), 
     2   lK(dof*dof,eNoN,eNoN) 

      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,3D0/)
      INTEGER a, b
      REAL(KIND=8) tauM, tauC, tauB, kS, kU, nu, rho, T1, T2, T3, divU, 
     2   up(nsd), ua(nsd), udNx(eNoN), updNx(eNoN), uadNx(eNoN), 
     3   rV(nsd), rM(nsd,nsd), NxdNx, amd, wl, wr, wrl, s, ud(nsd), 
     4   u(nsd), p, ux(nsd,nsd), px(nsd), f(nsd)

      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      nu   = eq(cEq)%dmn(cDmn)%prop(viscosity)/rho
      s    = nu/eq(cEq)%dmn(cDmn)%prop(permeability)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z) 
      
      wr  = w*rho
      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      wrl = wr*T1
      wl  = w*T1

!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         p  = p + N(a)*yl(4,a)
 
         ud(1) = ud(1) + N(a)*al(1,a)
         ud(2) = ud(2) + N(a)*al(2,a)
         ud(3) = ud(3) + N(a)*al(3,a)
      
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

      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU + ct(2)*nu*nu*kS + s*s)
      
      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1) 
     2      + u(3)*ux(3,1) + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2) 
     2      + u(3)*ux(3,2) + s*u(2))
      up(3) = -tauM*(ud(3) + px(3)/rho + u(1)*ux(1,3) + u(2)*ux(2,3) 
     2      + u(3)*ux(3,3) + s*u(3))
            
      tauC = ksix(1,1) + ksix(2,2) + ksix(3,3)
      tauC = 1D0/tauM/tauC/16D0
            
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

      rM(1,1) = nu*(ux(1,1) + ux(1,1)) - up(1)*ua(1) + rV(1)*up(1) 
     2        + divU*tauC - p/rho
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) - up(2)*ua(1) + rV(2)*up(1) 
      rM(3,1) = nu*(ux(3,1) + ux(1,3)) - up(3)*ua(1) + rV(3)*up(1) 
      
      rM(1,2) = nu*(ux(1,2) + ux(2,1)) - up(1)*ua(2) + rV(1)*up(2) 
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) - up(2)*ua(2) + rV(2)*up(2) 
     2        + divU*tauC - p/rho
      rM(3,2) = nu*(ux(3,2) + ux(2,3)) - up(3)*ua(2) + rV(3)*up(2) 
      
      rM(1,3) = nu*(ux(1,3) + ux(3,1)) - up(1)*ua(3) + rV(1)*up(3) 
      rM(2,3) = nu*(ux(2,3) + ux(3,2)) - up(2)*ua(3) + rV(2)*up(3) 
      rM(3,3) = nu*(ux(3,3) + ux(3,3)) - up(3)*ua(3) + rV(3)*up(3) 
     2        + divU*tauC - p/rho

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
            lK(1,a,b)  = lK(1,a,b)  + wrl*((nu + tauC)*rM(1,1) + T1)
            lK(2,a,b)  = lK(2,a,b)  + wrl*(nu*rM(2,1) + tauC*rM(1,2))
            lK(3,a,b)  = lK(3,a,b)  + wrl*(nu*rM(3,1) + tauC*rM(1,3))
            
            lK(5,a,b)  = lK(5,a,b)  + wrl*(nu*rM(1,2) + tauC*rM(2,1))
            lK(6,a,b)  = lK(6,a,b)  + wrl*((nu + tauC)*rM(2,2) + T1)
            lK(7,a,b)  = lK(7,a,b)  + wrl*(nu*rM(3,2) + tauC*rM(2,3))
            
            lK(9,a,b)  = lK(9,a,b)  + wrl*(nu*rM(1,3) + tauC*rM(3,1))
            lK(10,a,b) = lK(10,a,b) + wrl*(nu*rM(2,3) + tauC*rM(3,2))
            lK(11,a,b) = lK(11,a,b) + wrl*((nu + tauC)*rM(3,3) + T1)
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

!====================================================================
      
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

!====================================================================
      
      PURE SUBROUTINE FLUID2D (eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), 
     2   lK(dof*dof,eNoN,eNoN) 

      REAL(KIND=8), PARAMETER :: ct(2) = (/1D0,3D0/)
      INTEGER a, b
      REAL(KIND=8) tauM, tauC, tauB, kS, kU, nu, rho, T1, T2, T3, divU, 
     2   up(nsd), ua(nsd), udNx(eNoN), updNx(eNoN), uadNx(eNoN), 
     3   rV(nsd), rM(nsd,nsd), NxdNx, amd, wl, wr, wrl, s, ud(nsd), 
     4   u(nsd), p, ux(nsd,nsd), px(nsd), f(nsd)
  
      rho  = eq(cEq)%dmn(cDmn)%prop(fluid_density)
      nu   = eq(cEq)%dmn(cDmn)%prop(viscosity)/rho
      s    = nu/eq(cEq)%dmn(cDmn)%prop(permeability)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      
      wr  = w*rho
      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      wrl = wr*T1
      wl  = w*T1
 
!     Indices are not selected based on the equation only
!     because fluid equation always come first
      p  = 0D0
      u  = 0D0
      ud = -f
      px = 0D0
      ux = 0D0
      DO a=1, eNoN
         p = p + N(a)*yl(3,a)
      
         ud(1) = ud(1) + N(a)*al(1,a)
         ud(2) = ud(2) + N(a)*al(2,a)
               
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

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1)
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1)
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)
         
      tauM = 1D0/SQRT((2D0*ct(1)/dt)**2D0 + kU + ct(2)*nu*nu*kS + s*s)
            
      up(1) = -tauM*(ud(1) + px(1)/rho + u(1)*ux(1,1) + u(2)*ux(2,1)
     2      + s*u(1))
      up(2) = -tauM*(ud(2) + px(2)/rho + u(1)*ux(1,2) + u(2)*ux(2,2)
     2      + s*u(2))
            
      tauC = ksix(1,1) + ksix(2,2)
      tauC = 1D0/tauM/tauC/16D0
            
      tauB = up(1)*up(1)*ksix(1,1) + up(2)*up(1)*ksix(2,1) 
     2     + up(1)*up(2)*ksix(1,2) + up(2)*up(2)*ksix(2,2)

      IF (ISZERO(tauB)) tauB = eps
      tauB = 1D0/SQRT(tauB)

      divU = ux(1,1) + ux(2,2)
   
      ua(1) = u(1) + up(1)
      ua(2) = u(2) + up(2)
      
      rV(1) = tauB*(up(1)*ux(1,1) + up(2)*ux(2,1))
      rV(2) = tauB*(up(1)*ux(1,2) + up(2)*ux(2,2))

      rM(1,1) = nu*(ux(1,1) + ux(1,1)) - up(1)*ua(1) + rV(1)*up(1) 
     2        + divU*tauC - p/rho
      rM(2,1) = nu*(ux(2,1) + ux(1,2)) - up(2)*ua(1) + rV(2)*up(1) 
      
      rM(1,2) = nu*(ux(1,2) + ux(2,1)) - up(1)*ua(2) + rV(1)*up(2) 
      rM(2,2) = nu*(ux(2,2) + ux(2,2)) - up(2)*ua(2) + rV(2)*up(2) 
     2        + divU*tauC - p/rho
      
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
            lK(1,a,b)  = lK(1,a,b) + wrl*((nu + tauC)*rM(1,1) + T1)
            lK(2,a,b)  = lK(2,a,b) + wrl*(nu*rM(2,1) + tauC*rM(1,2))
            
            lK(4,a,b)  = lK(4,a,b) + wrl*(nu*rM(1,2) + tauC*rM(2,1))
            lK(5,a,b)  = lK(5,a,b) + wrl*((nu + tauC)*rM(2,2) + T1)
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


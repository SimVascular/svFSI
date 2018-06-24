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
!     This is for solving heat equation in a fluid/solid.
!      
!--------------------------------------------------------------------

!     This is for solving heat equation in a fluid.
      PURE SUBROUTINE HEATF3D (eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,3D0,1D0/)
      INTEGER i, a, b
      REAL(KIND=8) nu, tauM, kU, kS, nTx, Tp, udTx, udNx(eNoN), T1, 
     2   amd, wl, Td, Tx(nsd), u(nsd), s
 
      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s

      wl  = w*T1

      u  = 0D0
      Td = -s
      Tx = 0D0
      DO a=1, eNoN
         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
         u(3) = u(3) + N(a)*yl(3,a)
               
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)
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
         
      nTx = ksix(1,1)*Tx(1)*Tx(1) 
     2    + ksix(2,2)*Tx(2)*Tx(2)
     3    + ksix(3,3)*Tx(3)*Tx(3) 
     4    + (ksix(1,2) + ksix(2,1))*Tx(1)*Tx(2)
     5    + (ksix(1,3) + ksix(3,1))*Tx(1)*Tx(3)
     6    + (ksix(2,3) + ksix(3,2))*Tx(2)*Tx(3)
            
      IF (ISZERO(nTx)) nTx = eps
      
      udTx = u(1)*Tx(1) + u(2)*Tx(2) + u(3)*Tx(3)
      
      Tp = ABS(Td + udTx)
      
!     nu = nu + kDC
      nu = nu + Tp/SQRT(nTx)/2D0 

      tauM = ct(4)/SQRT(ct(1)/dt/dt + ct(2)*kU + ct(3)*nu*nu*kS)

      Tp = -tauM*(Td + udTx)
      
      DO a=1,eNoN
         udNx(a) = u(1)*Nx(1,a) + u(2)*Nx(2,a) + u(3)*Nx(3,a)
      END DO

      DO a=1,eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*(Td + udTx) + (Nx(1,a)*Tx(1) 
     2      + Nx(2,a)*Tx(2) + Nx(3,a)*Tx(3))*nu - udNx(a)*Tp)
         
         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(nu*(Nx(1,a)*Nx(1,b) 
     2         + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b))
     3         + (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATF3D
!====================================================================
      PURE SUBROUTINE BHEATF (eNoN, w, N, y, h, nV, lR, lK)

      USE COMMOD

      IMPLICIT NONE
       
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), y(tDof), h, nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN), lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i
      REAL(KIND=8) T1, wl, T, udn
 
      wl = w*eq(cEq)%af*eq(cEq)%gam*dt
      T  = y(eq(cEq)%s)
      
      udn = 0D0
      DO i=1, nsd
         udn = udn + y(i)*nV(i)
      END DO
      udn = 5D-1*(udn - ABS(udn))
      T1  = h - udn*T

!     Here the loop is started for constructing left and right hand side
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*T1
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) - wl*N(a)*N(b)*udn
         END DO
      END DO

      RETURN
      END SUBROUTINE BHEATF

!####################################################################
!     This is for solving the heat equation in a solid      
      PURE SUBROUTINE HEATS3D (eNoN, w, N, Nx, al, yl, lR, lK)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)
      
      INTEGER i, a, b
      REAL(KIND=8) nu, T1, amd, wl, Td, Tx(nsd), s

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s
 
      wl = w*T1
      
      Td = -s
      Tx = 0D0
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)
      
         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)
      END DO

      DO a=1,eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td
     2      + (Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2) + Nx(3,a)*Tx(3))*nu)
         
         DO b=1,eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd 
     2         + nu*(Nx(1,a)*Nx(1,b) +Nx(2,a)*Nx(2,b) +Nx(3,a)*Nx(3,b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATS3D
!====================================================================
      PURE SUBROUTINE BHEATS (eNoN, w, N, h, lR)

      USE COMMOD

      IMPLICIT NONE
       
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER a

!     Here the loop is started for constructing left and right hand side
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*h
      END DO

      RETURN
      END SUBROUTINE BHEATS

!####################################################################

!     2D versions

!####################################################################
!     This is for solving heat equation in a fluid
      PURE SUBROUTINE HEATF2D (eNoN, w, N, Nx, al, yl, ksix, lR, lK)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN), ksix(nsd,nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: ct(4) = (/4D0,1D0,3D0,1D0/)
      INTEGER i, a, b
      REAL(KIND=8) nu, tauM, kU, kS, nTx, Tp, udTx, udNx(eNoN), T1, 
     2   amd, wl, Td, Tx(nsd), u(nsd), s
 
      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s
      
      wl  = w*T1

      u  = 0D0
      Td = -s
      Tx = 0D0
      DO a=1, eNoN
         u(1) = u(1) + N(a)*yl(1,a)
         u(2) = u(2) + N(a)*yl(2,a)
               
         Td = Td + N(a)*al(i,a)
      
         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
      END DO

      IF (mvMsh) THEN
         DO a=1, eNoN
            u(1) = u(1) - N(a)*yl(5,a)
            u(2) = u(2) - N(a)*yl(6,a)
            u(3) = u(3) - N(a)*yl(7,a)
         END DO
      END IF

      kU = u(1)*u(1)*ksix(1,1) + u(2)*u(1)*ksix(2,1) 
     2   + u(1)*u(2)*ksix(1,2) + u(2)*u(2)*ksix(2,2)

      kS = ksix(1,1)*ksix(1,1) + ksix(2,1)*ksix(2,1) 
     2   + ksix(1,2)*ksix(1,2) + ksix(2,2)*ksix(2,2)
         
      nTx = ksix(1,1)*Tx(1)*Tx(1) + ksix(2,2)*Tx(2)*Tx(2)
     4    + (ksix(1,2) + ksix(2,1))*Tx(1)*Tx(2)
      
      IF (ISZERO(nTx)) nTx = eps
      
      udTx = u(1)*Tx(1) + u(2)*Tx(2)
      
      Tp = ABS(Td + udTx)

!     nu = nu + kDC
      nu = nu + Tp/SQRT(nTx)/2D0 

      tauM = ct(4)/SQRT(ct(1)/dt/dt + ct(2)*kU + ct(3)*nu*nu*kS)

      Tp = -tauM*(Td + udTx)
      
      DO a=1, eNoN
         udNx(a) = u(1)*Nx(1,a) + u(2)*Nx(2,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*(Td + udTx)
     2      + (Nx(1,a)*Tx(1) + Nx(2,a)*Tx(2))*nu - udNx(a)*Tp)
         
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(nu*(Nx(1,a)*Nx(1,b) 
     2         + Nx(2,a)*Nx(2,b)) 
     3         + (N(a) + tauM*udNx(a))*(N(b)*amd + udNx(b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATF2D

!####################################################################
!     This is for solving the heat equation in a solid      
      PURE SUBROUTINE HEATS2D (eNoN, w, N, Nx, al, yl, lR, lK)

      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN), 
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)
 
      INTEGER i, a, b
      REAL(KIND=8) nu, T1, amd, wl, Td, Tx(nsd), s

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am/T1
      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      i   = eq(cEq)%s
 
      wl = w*T1
 
      Td = -s
      Tx = 0D0
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)
      
         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td + (Nx(1,a)*Tx(1) 
     2      + Nx(2,a)*Tx(2))*nu)
         
         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd
     2         + nu*(Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)))
         END DO
      END DO

      RETURN
      END SUBROUTINE HEATS2D

!####################################################################

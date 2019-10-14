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
!     This is for solving linear elasticity.
!
!--------------------------------------------------------------------

      PURE SUBROUTINE LELAS3D (eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl,
     2   lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN), al(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(3,eNoN), pS0l(6,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: pSl(6), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j, k
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl,
     2   lDm, ed(6), ud(3), f(3), S0(6), S(6)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      ed = 0D0
      ud = -f
      S0 = 0D0
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(al(i,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(j,a)-bfl(2,a))
         ud(3) = ud(3) + N(a)*(al(k,a)-bfl(3,a))

         ed(1) = ed(1) + Nx(1,a)*dl(i,a)
         ed(2) = ed(2) + Nx(2,a)*dl(j,a)
         ed(3) = ed(3) + Nx(3,a)*dl(k,a)
         ed(4) = ed(4) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
         ed(5) = ed(5) + Nx(3,a)*dl(i,a) + Nx(1,a)*dl(k,a)
         ed(6) = ed(6) + Nx(3,a)*dl(j,a) + Nx(2,a)*dl(k,a)

         S0(1) = S0(1) + N(a)*pS0l(1,a)
         S0(2) = S0(2) + N(a)*pS0l(2,a)
         S0(3) = S0(3) + N(a)*pS0l(3,a)
         S0(4) = S0(4) + N(a)*pS0l(4,a)
         S0(5) = S0(5) + N(a)*pS0l(5,a)
         S0(6) = S0(6) + N(a)*pS0l(6,a)
      END DO
      divD = lambda*(ed(1) + ed(2) + ed(3))

!     Stress in Voigt notation
      S(1) = divD + 2D0*mu*ed(1)
      S(2) = divD + 2D0*mu*ed(2)
      S(3) = divD + 2D0*mu*ed(3)
      S(4) = mu*ed(4)  ! 2*eps_12
      S(5) = mu*ed(5)  ! 2*eps_13
      S(6) = mu*ed(6)  ! 2*eps_23
      pSl  = S

!     Add prestress contribution
      S = pSl + S0

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(1,a)*S(1) +
     2      Nx(2,a)*S(4) + Nx(3,a)*S(5))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(1,a)*S(4) +
     2      Nx(2,a)*S(2) + Nx(3,a)*S(6))

         lR(3,a) = lR(3,a) + w*(rho*N(a)*ud(3) + Nx(1,a)*S(5) +
     2      Nx(2,a)*S(6) + Nx(3,a)*S(3))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = amd*N(a)*N(b)/mu + NxdNx

            lK(1,a,b) = lK(1,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(1,b))

            lK(3,a,b) = lK(3,a,b) + wl*(lDm*Nx(1,a)*Nx(3,b)
     2         + Nx(3,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(2,a)*Nx(2,b))

            lK(dof+3,a,b) = lK(dof+3,a,b) + wl*(lDm*Nx(2,a)*Nx(3,b)
     2         + Nx(3,a)*Nx(2,b))

            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + wl*(lDm*Nx(3,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(3,b))

            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + wl*(lDm*Nx(3,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(3,b))

            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(3,a)*Nx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS3D
!####################################################################
!     This is for boundary elasticity equation
      PURE SUBROUTINE BLELAS (eNoN, w, N, h, nV, lR)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), h, nV(nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER a, i
      REAL(KIND=8) hc(nsd)

      hc = h*nV
      DO a=1,eNoN
         DO i=1, nsd
            lR(i,a) = lR(i,a) - w*N(a)*hc(i)
         END DO
      END DO

      RETURN
      END SUBROUTINE BLELAS
!####################################################################
      PURE SUBROUTINE LELAS2D (eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl,
     2   lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), al(tDof,eNoN),
     2   dl(tDof,eNoN), bfl(2,eNoN), pS0l(3,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: pSl(3), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd, wl,
     2   lDm, ed(3), ud(2), f(2), S0(3), S(3)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      i    = eq(cEq)%s
      j    = i + 1

      lambda = elM*nu/(1D0 + nu)/(1D0 - 2D0*nu)
      mu     = elM/2D0/(1D0 + nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      ed = 0D0
      ud = -f
      S0 = 0D0
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(al(i,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(j,a)-bfl(2,a))

         ed(1) = ed(1) + Nx(1,a)*dl(i,a)
         ed(2) = ed(2) + Nx(2,a)*dl(j,a)
         ed(3) = ed(3) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)

         S0(1) = S0(1) + N(a)*pS0l(1,a)
         S0(2) = S0(2) + N(a)*pS0l(2,a)
         S0(3) = S0(3) + N(a)*pS0l(3,a)
      END DO
      divD = lambda*(ed(1) + ed(2))

!     Stress in Voigt notation
      S(1) = divD + 2D0*mu*ed(1)
      S(2) = divD + 2D0*mu*ed(2)
      S(3) = mu*ed(3)  ! 2*eps_12
      pSl  = S

!     Add prestress contribution
      S = pSl + S0

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(1,a)*S(1) +
     2      Nx(2,a)*S(3))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(1,a)*S(3) +
     2      Nx(2,a)*S(2))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T1 = amd*N(a)*N(b)/mu + NxdNx

            lK(1,a,b) = lK(1,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T1
     2         + (1D0 + lDm)*Nx(2,a)*Nx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS2D
!####################################################################

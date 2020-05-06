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
!     This routines is for solving nonlinear structural mechanics
!     problem using pure displacement-based formulation.
!
!--------------------------------------------------------------------

      SUBROUTINE STRUCT3D (eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2   pS0l, pSl, ya_l, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(3,eNoN),
     3   fN(3,nFn), pS0l(6,eNoN), ya_l(eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: pSl(6)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j, k
      REAL(KIND=RKIND) :: rho, dmp, T1, amd, afl, ya_g, fb(3), ud(3),
     2   NxSNx, BmDBm, F(3,3), S(3,3), P(3,3), Dm(6,6), DBm(6,3),
     3   Bm(6,3,eNoN), CC(3,3,3,3), S0(3,3)
      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp     = eq(cEq)%dmn(cDmn)%prop(damping)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      fb(3)   = eq(cEq)%dmn(cDmn)%prop(f_z)
      amd     = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i       = eq(cEq)%s
      j       = i + 1
      k       = j + 1

!     Inertia, body force and deformation tensor (F)
      ud     = -rho*fb
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      S0     = 0._RKIND
      ya_g   = 0._RKIND
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))
         ud(3) = ud(3) + N(a)*(rho*(al(k,a)-bfl(3,a)) + dmp*yl(k,a))

         F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
         F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
         F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
         F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
         F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
         F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
         F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)

         S0(1,1) = S0(1,1) + N(a)*pS0l(1,a)
         S0(2,2) = S0(2,2) + N(a)*pS0l(2,a)
         S0(3,3) = S0(3,3) + N(a)*pS0l(3,a)
         S0(1,2) = S0(1,2) + N(a)*pS0l(4,a)
         S0(1,3) = S0(1,3) + N(a)*pS0l(5,a)
         S0(2,3) = S0(2,3) + N(a)*pS0l(6,a)

         ya_g    = ya_g + N(a)*ya_l(a)
      END DO
      S0(2,1) = S0(1,2)
      S0(3,1) = S0(1,3)
      S0(3,2) = S0(2,3)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
      CALL GETPK2CC(stModel, F, nFn, fN, ya_g, S, CC)

!     Prestress
      pSl(1) = S(1,1)
      pSl(2) = S(2,2)
      pSl(3) = S(3,3)
      pSl(4) = S(1,2)
      pSl(5) = S(1,3)
      pSl(6) = S(2,3)
      S = S + S0

!     Active stress - electromechanics
      IF (cem%aStress) CALL ACTVSTRESS(ya_g, nFn, fN, S)

!     1st Piola-Kirchhoff tensor (P)
      P = MATMUL(F, S)

!     Convert to Voigt Notation
      Dm(1,1) = CC(1,1,1,1)
      Dm(1,2) = CC(1,1,2,2)
      Dm(1,3) = CC(1,1,3,3)
      Dm(1,4) = CC(1,1,1,2)
      Dm(1,5) = CC(1,1,2,3)
      Dm(1,6) = CC(1,1,3,1)

      Dm(2,2) = CC(2,2,2,2)
      Dm(2,3) = CC(2,2,3,3)
      Dm(2,4) = CC(2,2,1,2)
      Dm(2,5) = CC(2,2,2,3)
      Dm(2,6) = CC(2,2,3,1)

      Dm(3,3) = CC(3,3,3,3)
      Dm(3,4) = CC(3,3,1,2)
      Dm(3,5) = CC(3,3,2,3)
      Dm(3,6) = CC(3,3,3,1)

      Dm(4,4) = CC(1,2,1,2)
      Dm(4,5) = CC(1,2,2,3)
      Dm(4,6) = CC(1,2,3,1)

      Dm(5,5) = CC(2,3,2,3)
      Dm(5,6) = CC(2,3,3,1)

      Dm(6,6) = CC(3,1,3,1)

      DO a=2, 6
         DO b=1, a-1
            Dm(a,b) = Dm(b,a)
         END DO
      END DO

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)
         Bm(1,3,a) = Nx(1,a)*F(3,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)
         Bm(2,3,a) = Nx(2,a)*F(3,2)

         Bm(3,1,a) = Nx(3,a)*F(1,3)
         Bm(3,2,a) = Nx(3,a)*F(2,3)
         Bm(3,3,a) = Nx(3,a)*F(3,3)

         Bm(4,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(4,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
         Bm(4,3,a) = (Nx(1,a)*F(3,2) + F(3,1)*Nx(2,a))

         Bm(5,1,a) = (Nx(2,a)*F(1,3) + F(1,2)*Nx(3,a))
         Bm(5,2,a) = (Nx(2,a)*F(2,3) + F(2,2)*Nx(3,a))
         Bm(5,3,a) = (Nx(2,a)*F(3,3) + F(3,2)*Nx(3,a))

         Bm(6,1,a) = (Nx(3,a)*F(1,1) + F(1,3)*Nx(1,a))
         Bm(6,2,a) = (Nx(3,a)*F(2,1) + F(2,3)*Nx(1,a))
         Bm(6,3,a) = (Nx(3,a)*F(3,1) + F(3,3)*Nx(1,a))
      END DO

!     Local residue and tangent matrices
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(1,a)*P(1,1) +
     2      Nx(2,a)*P(1,2) + Nx(3,a)*P(1,3))
         lR(2,a) = lR(2,a) + w*(N(a)*ud(2) + Nx(1,a)*P(2,1) +
     2      Nx(2,a)*P(2,2) + Nx(3,a)*P(2,3))
         lR(3,a) = lR(3,a) + w*(N(a)*ud(3) + Nx(1,a)*P(3,1) +
     2      Nx(2,a)*P(3,2) + Nx(3,a)*P(3,3))

         DO b=1, eNoN
!           Geometric stiffness
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(2,a)*S(2,1)*Nx(1,b) +
     2              Nx(3,a)*S(3,1)*Nx(1,b) + Nx(1,a)*S(1,2)*Nx(2,b) +
     3              Nx(2,a)*S(2,2)*Nx(2,b) + Nx(3,a)*S(3,2)*Nx(2,b) +
     4              Nx(1,a)*S(1,3)*Nx(3,b) + Nx(2,a)*S(2,3)*Nx(3,b) +
     5              Nx(3,a)*S(3,3)*Nx(3,b)
            T1 = amd*N(a)*N(b) + afl*NxSNx

!           Material Stiffness (Bt*D*B)
            DBm = MATMUL(Dm, Bm(:,:,b))

            BmDBm = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2              Bm(3,1,a)*DBm(3,1) + Bm(4,1,a)*DBm(4,1) +
     2              Bm(5,1,a)*DBm(5,1) + Bm(6,1,a)*DBm(6,1)
            lK(1,a,b) = lK(1,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2              Bm(3,1,a)*DBm(3,2) + Bm(4,1,a)*DBm(4,2) +
     2              Bm(5,1,a)*DBm(5,2) + Bm(6,1,a)*DBm(6,2)
            lK(2,a,b) = lK(2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,1,a)*DBm(1,3) + Bm(2,1,a)*DBm(2,3) +
     2              Bm(3,1,a)*DBm(3,3) + Bm(4,1,a)*DBm(4,3) +
     2              Bm(5,1,a)*DBm(5,3) + Bm(6,1,a)*DBm(6,3)
            lK(3,a,b) = lK(3,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2              Bm(3,2,a)*DBm(3,1) + Bm(4,2,a)*DBm(4,1) +
     2              Bm(5,2,a)*DBm(5,1) + Bm(6,2,a)*DBm(6,1)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2              Bm(3,2,a)*DBm(3,2) + Bm(4,2,a)*DBm(4,2) +
     2              Bm(5,2,a)*DBm(5,2) + Bm(6,2,a)*DBm(6,2)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,2,a)*DBm(1,3) + Bm(2,2,a)*DBm(2,3) +
     2              Bm(3,2,a)*DBm(3,3) + Bm(4,2,a)*DBm(4,3) +
     2              Bm(5,2,a)*DBm(5,3) + Bm(6,2,a)*DBm(6,3)
            lK(dof+3,a,b) = lK(dof+3,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,1) + Bm(2,3,a)*DBm(2,1) +
     2              Bm(3,3,a)*DBm(3,1) + Bm(4,3,a)*DBm(4,1) +
     2              Bm(5,3,a)*DBm(5,1) + Bm(6,3,a)*DBm(6,1)
            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,2) + Bm(2,3,a)*DBm(2,2) +
     2              Bm(3,3,a)*DBm(3,2) + Bm(4,3,a)*DBm(4,2) +
     2              Bm(5,3,a)*DBm(5,2) + Bm(6,3,a)*DBm(6,2)
            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,3,a)*DBm(1,3) + Bm(2,3,a)*DBm(2,3) +
     2              Bm(3,3,a)*DBm(3,3) + Bm(4,3,a)*DBm(4,3) +
     2              Bm(5,3,a)*DBm(5,3) + Bm(6,3,a)*DBm(6,3)
            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + w*(T1 + afl*BmDBm)
         END DO
      END DO

      RETURN
      END SUBROUTINE STRUCT3D
!####################################################################
      SUBROUTINE STRUCT2D (eNoN, nFn, w, N, Nx, al, yl, dl, bfl, fN,
     2   pS0l, pSl, ya_l, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), bfl(2,eNoN),
     3   fN(2,nFn), pS0l(3,eNoN), ya_l(eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: pSl(3)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j
      REAL(KIND=RKIND) :: rho, dmp, T1, amd, afl, ya_g, fb(2), ud(2),
     2   NxSNx, BmDBm, F(2,2), S(2,2), P(2,2), Dm(3,3), DBm(3,2),
     3   Bm(3,2,eNoN), CC(2,2,2,2), S0(2,2)
      TYPE (stModelType) :: stModel

!     Define parameters
      stModel = eq(cEq)%dmn(cDmn)%stM
      rho     = eq(cEq)%dmn(cDmn)%prop(solid_density)
      dmp     = eq(cEq)%dmn(cDmn)%prop(damping)
      fb(1)   = eq(cEq)%dmn(cDmn)%prop(f_x)
      fb(2)   = eq(cEq)%dmn(cDmn)%prop(f_y)
      amd     = eq(cEq)%am*rho + eq(cEq)%af*eq(cEq)%gam*dt*dmp
      afl     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i       = eq(cEq)%s
      j       = i + 1

!     Inertia, body force and deformation tensor (F)
      ud     = -rho*fb
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      S0     = 0._RKIND
      ya_g   = 0._RKIND
      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(rho*(al(i,a)-bfl(1,a)) + dmp*yl(i,a))
         ud(2) = ud(2) + N(a)*(rho*(al(j,a)-bfl(2,a)) + dmp*yl(j,a))

         F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
         F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
         F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
         F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)

         S0(1,1) = S0(1,1) + N(a)*pS0l(1,a)
         S0(2,2) = S0(2,2) + N(a)*pS0l(2,a)
         S0(1,2) = S0(1,2) + N(a)*pS0l(3,a)

         ya_g    = ya_g + N(a)*ya_l(a)
      END DO
      S0(2,1) = S0(1,2)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
      CALL GETPK2CC(stModel, F, nFn, fN, ya_g, S, CC)

!     Prestress
      pSl(1) = S(1,1)
      pSl(2) = S(2,2)
      pSl(3) = S(1,2)
      S = S + S0

!     Active stress - electromechanics
      IF (cem%aStress) CALL ACTVSTRESS(ya_g, nFn, fN, S)

!     1st Piola-Kirchhoff tensor (P)
      P = MATMUL(F, S)

!     Convert to Voigt Notation
      Dm(1,1) = CC(1,1,1,1)
      Dm(1,2) = CC(1,1,2,2)
      Dm(1,3) = CC(1,1,1,2)

      Dm(2,2) = CC(2,2,2,2)
      Dm(2,3) = CC(2,2,1,2)

      Dm(3,3) = CC(1,2,1,2)

      Dm(2,1) = Dm(1,2)
      Dm(3,1) = Dm(1,3)
      Dm(3,2) = Dm(2,3)

      DO a=1, eNoN
         Bm(1,1,a) = Nx(1,a)*F(1,1)
         Bm(1,2,a) = Nx(1,a)*F(2,1)

         Bm(2,1,a) = Nx(2,a)*F(1,2)
         Bm(2,2,a) = Nx(2,a)*F(2,2)

         Bm(3,1,a) = (Nx(1,a)*F(1,2) + F(1,1)*Nx(2,a))
         Bm(3,2,a) = (Nx(1,a)*F(2,2) + F(2,1)*Nx(2,a))
      END DO

!     Local residue and tangent matrices
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*ud(1) + Nx(1,a)*P(1,1) +
     2      Nx(2,a)*P(1,2))
         lR(2,a) = lR(2,a) + w*(N(a)*ud(2) + Nx(1,a)*P(2,1) +
     2      Nx(2,a)*P(2,2))

         DO b=1, eNoN
!           Geometric stiffness
            NxSNx = Nx(1,a)*S(1,1)*Nx(1,b) + Nx(2,a)*S(2,1)*Nx(1,b) +
     2              Nx(1,a)*S(1,2)*Nx(2,b) + Nx(2,a)*S(2,2)*Nx(2,b)
            T1 = amd*N(a)*N(b) + afl*NxSNx

!           Material stiffness (Bt*D*B)
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

            BmDBm = Bm(1,1,a)*DBm(1,1) + Bm(2,1,a)*DBm(2,1) +
     2         Bm(3,1,a)*DBm(3,1)
            lK(1,a,b) = lK(1,a,b) + w*(T1 + afl*BmDBm)

            BmDBm = Bm(1,1,a)*DBm(1,2) + Bm(2,1,a)*DBm(2,2) +
     2         Bm(3,1,a)*DBm(3,2)
            lK(2,a,b) = lK(2,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,1) + Bm(2,2,a)*DBm(2,1) +
     2         Bm(3,2,a)*DBm(3,1)
            lK(dof+1,a,b) = lK(dof+1,a,b) + w*afl*BmDBm

            BmDBm = Bm(1,2,a)*DBm(1,2) + Bm(2,2,a)*DBm(2,2) +
     2         Bm(3,2,a)*DBm(3,2)
            lK(dof+2,a,b) = lK(dof+2,a,b) + w*(T1 + afl*BmDBm)
         END DO
      END DO

      RETURN
      END SUBROUTINE STRUCT2D
!####################################################################

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

      SUBROUTINE CONSTRUCT_LELAS(lM, Ag, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), dl(:,:),
     2   bfl(:,:), pS0l(:,:), pSl(:), N(:), Nx(:,:), lR(:,:), lK(:,:,:),
     3   lVWP(:,:)

      eNoN = lM%eNoN

!     LELAS: dof = nsd
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), dl(tDof,eNoN),
     2   bfl(nsd,eNoN), pS0l(nsymd,eNoN), pSl(nsymd), N(eNoN),
     3   Nx(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     4   lVWP(nvwp,eNoN))


!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_lElas) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         pS0l = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            bfl(:,a) = Bf(:,Ac)
            IF (ALLOCATED(pS0)) pS0l(:,a) = pS0(:,Ac)
!           Varwall properties------------------------------------------
!           Calculate local wall property
            IF (useVarWall) THEN
               lVWP(:,a) = vWP0(:,Ac)
            END IF
!        ---------------------------------------------------------------
         END DO


!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g) * Jac
            N = lM%N(:,g)

            pSl = 0._RKIND
            IF (nsd .EQ. 3) THEN
               CALL LELAS3D(eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR,
     2            lK, lVWP)

            ELSE IF (nsd .EQ. 2) THEN
               CALL LELAS2D(eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR,
     2            lK)

            END IF

!           Prestress
            IF (pstEq) THEN
               DO a=1, eNoN
                  Ac = ptr(a)
                  pSn(:,Ac) = pSn(:,Ac) + w*N(a)*pSl(:)
                  pSa(Ac)   = pSa(Ac)   + w*N(a)
               END DO
            END IF
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
      END DO ! e: loop

      DEALLOCATE(ptr, xl, al, dl, bfl, pS0l, pSl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_LELAS
!####################################################################
      PURE SUBROUTINE LELAS3D (eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl,
     2   lR, lK, lVWP)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), dl(tDof,eNoN), bfl(3,eNoN), pS0l(6,eNoN),
     3   lVWP(nvwp,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: pSl(6), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i, j, k
      REAL(KIND=RKIND) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd,
     2   wl, lDm, ed(6), ud(3), f(3), S0(6), S(6), eVWP(nvwp), Cst(6,6)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

      ed = 0._RKIND
      ud = -f
      S0 = 0._RKIND
      eVWP = 0._RKIND

      DO a=1, eNoN
         ud(1) = ud(1) + N(a)*(al(i,a)-bfl(1,a))
         ud(2) = ud(2) + N(a)*(al(j,a)-bfl(2,a))
         ud(3) = ud(3) + N(a)*(al(k,a)-bfl(3,a))

         ed(1) = ed(1) + Nx(1,a)*dl(i,a)
         ed(2) = ed(2) + Nx(2,a)*dl(j,a)
         ed(3) = ed(3) + Nx(3,a)*dl(k,a)
         ed(4) = ed(4) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
         ed(5) = ed(5) + Nx(3,a)*dl(j,a) + Nx(2,a)*dl(k,a)
         ed(6) = ed(6) + Nx(1,a)*dl(k,a) + Nx(3,a)*dl(i,a)

!     ------------------------------------------------------------------
!     Calculate local wall property
!     Don't use if calculating mesh
         IF (useVarWall .AND. (phys_mesh .NE. eq(cEq)%phys)) THEN
            eVWP(:) = eVWP(:) + N(a)*lVWP(:,a)
         END IF
!     ------------------------------------------------------------------

         S0(1) = S0(1) + N(a)*pS0l(1,a)
         S0(2) = S0(2) + N(a)*pS0l(2,a)
         S0(3) = S0(3) + N(a)*pS0l(3,a)
         S0(4) = S0(4) + N(a)*pS0l(4,a)
         S0(5) = S0(5) + N(a)*pS0l(5,a)
         S0(6) = S0(6) + N(a)*pS0l(6,a)
      END DO


      IF (useVarWall .AND. (phys_mesh .NE. eq(cEq)%phys)) THEN
         elM  = eVWP(1)
         nu   = evWP(2)
      ELSE
         elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
         nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      END IF
      
      lambda = elM*nu / (1._RKIND+nu) / (1._RKIND-2._RKIND*nu)
      mu     = elM * 0.5_RKIND / (1._RKIND+nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      divD = lambda*(ed(1) + ed(2) + ed(3))


!     Stress in Voigt notation
      S(1) = divD + 2._RKIND*mu*ed(1)
      S(2) = divD + 2._RKIND*mu*ed(2)
      S(3) = divD + 2._RKIND*mu*ed(3)
      S(4) = mu*ed(4)  ! 2*eps_12
      S(5) = mu*ed(5)  ! 2*eps_23
      S(6) = mu*ed(6)  ! 2*eps_13

      pSl  = S

!     Add prestress contribution
      S = pSl + S0

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(1,a)*S(1) +
     2      Nx(2,a)*S(4) + Nx(3,a)*S(6))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(1,a)*S(4) +
     2      Nx(2,a)*S(2) + Nx(3,a)*S(5))

         lR(3,a) = lR(3,a) + w*(rho*N(a)*ud(3) + Nx(1,a)*S(6) +
     2      Nx(2,a)*S(5) + Nx(3,a)*S(3))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b) + Nx(3,a)*Nx(3,b)

            T1 = amd*N(a)*N(b)/mu + NxdNx

            lK(1,a,b) = lK(1,a,b) + wl*(T1
     2         + (1._RKIND + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(1,b))

            lK(3,a,b) = lK(3,a,b) + wl*(lDm*Nx(1,a)*Nx(3,b)
     2         + Nx(3,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T1
     2         + (1._RKIND + lDm)*Nx(2,a)*Nx(2,b))

            lK(dof+3,a,b) = lK(dof+3,a,b) + wl*(lDm*Nx(2,a)*Nx(3,b)
     2         + Nx(3,a)*Nx(2,b))

            lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + wl*(lDm*Nx(3,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(3,b))

            lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + wl*(lDm*Nx(3,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(3,b))

            lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + wl*(T1
     2         + (1._RKIND + lDm)*Nx(3,a)*Nx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS3D
!--------------------------------------------------------------------
      PURE SUBROUTINE LELAS2D (eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl,
     2   lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), dl(tDof,eNoN), bfl(2,eNoN), pS0l(3,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: pSl(3), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i, j
      REAL(KIND=RKIND) NxdNx, rho, elM, nu, lambda, mu, divD, T1, amd,
     2   T2, wl, lDm, ed(3), ud(2), f(2), S0(3), S(3)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM  = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu   = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      i    = eq(cEq)%s
      j    = i + 1

      lambda = elM*nu / (1._RKIND+nu) / (1._RKIND-2._RKIND*nu)
      mu     = elM * 0.5_RKIND / (1._RKIND+nu)
      lDm    = lambda/mu
      T1     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd    = eq(cEq)%am/T1*rho
      wl     = w*T1*mu

      ed = 0._RKIND
      ud = -f
      S0 = 0._RKIND
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
      S(1) = divD + 2._RKIND*mu*ed(1)
      S(2) = divD + 2._RKIND*mu*ed(2)
      S(3) = mu*ed(3)  ! 2*eps_12
      pSl  = S

!     Add prestress contribution
      S = pSl + S0

!     Need to add variable wall tangent matrix

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(rho*N(a)*ud(1) + Nx(1,a)*S(1) +
     2      Nx(2,a)*S(3))

         lR(2,a) = lR(2,a) + w*(rho*N(a)*ud(2) + Nx(1,a)*S(3) +
     2      Nx(2,a)*S(2))

         DO b=1, eNoN
            NxdNx = Nx(1,a)*Nx(1,b) + Nx(2,a)*Nx(2,b)

            T2 = amd*N(a)*N(b)/mu + NxdNx

            lK(1,a,b) = lK(1,a,b) + wl*(T2
     2         + (1._RKIND + lDm)*Nx(1,a)*Nx(1,b))

            lK(2,a,b) = lK(2,a,b) + wl*(lDm*Nx(1,a)*Nx(2,b)
     2         + Nx(2,a)*Nx(1,b))

            lK(dof+1,a,b) = lK(dof+1,a,b) + wl*(lDm*Nx(2,a)*Nx(1,b)
     2         + Nx(1,a)*Nx(2,b))

            lK(dof+2,a,b) = lK(dof+2,a,b) + wl*(T2
     2         + (1._RKIND + lDm)*Nx(2,a)*Nx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE LELAS2D
!####################################################################
!     This is for boundary elasticity equation
      PURE SUBROUTINE BLELAS (eNoN, w, N, h, nV, lR)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), h, nV(nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER(KIND=IKIND) a, i
      REAL(KIND=RKIND) hc(nsd)

      hc = h*nV
      DO a=1,eNoN
         DO i=1, nsd
            lR(i,a) = lR(i,a) - w*N(a)*hc(i)
         END DO
      END DO

      RETURN
      END SUBROUTINE BLELAS
!####################################################################

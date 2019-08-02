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
!     These subroutines implement the Coupled Momentum Method (CMM)
!     as an extension of the tools inside MUPFES already. I have made
!     modifications to add CMM as an equation type, and a new boundary
!     type. The boundary treatment will be considered in these
!     subroutines. I chose to implement CMM as a new top level BC
!     that is in the same "category" as Dirichlet or Neumann BC's
!     so that I can have more freedom in what data I can pass to
!     the subroutines, and to avoid breaking the existing code.
!
!--------------------------------------------------------------------

      SUBROUTINE CMMi(lM, eNoN, al, dl, xl, bfl, pS0l, ptr)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: eNoN, ptr(eNoN)
      REAL(KIND=8), INTENT(IN) :: al(tDof,eNoN), dl(tDof,eNoN),
     2   xl(3,eNoN), bfl(3,eNoN), pS0l(6)

      INTEGER a, g, Ac
      REAL(KIND=8) :: w, Jac, nV(3), pSl(6), gC(3,2)

      REAL(KIND=8), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      ALLOCATE(N(eNoN), Nx(2,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
      lK  = 0D0
      lR  = 0D0

      Nx  = lM%Nx(:,:,1)
      CALL GNNS(eNoN, Nx, xl, nV, gC, gC)
      Jac = SQRT(NORM(nV))
      nV  = nV / Jac

!     Internal stresses (stiffness) contribution
      pSl = 0D0
      CALL CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, pSl, lR, lK)

!     Inertia and body forces (mass) contribution
      DO g=1, lM%nG
         N(:) = lM%N(:,g)
         w = lM%w(g)*Jac
         CALL CMM_MASS(eNoN, w, N, al, bfl, lR, lK)

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
      IF (useTrilinosAssemAndLS) THEN
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
      SUBROUTINE CMMb(lFa, e, eNoN, al, dl, xl, bfl, pS0l, ptr)
      USE COMMOD
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER, INTENT(IN) :: e, eNoN, ptr(eNoN)
      REAL(KIND=8), INTENT(IN) :: al(tDof,eNoN), dl(tDof,eNoN),
     2   xl(3,eNoN), bfl(3,eNoN), pS0l(6)

      INTEGER a, g, Ac
      REAL(KIND=8) :: w, Jac, nV(3), pSl(6)

      REAL(KIND=8), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      ALLOCATE(N(eNoN), Nx(2,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
      lK  = 0D0
      lR  = 0D0

!     Internal stresses (stiffness) contribution
      pSl = 0D0
      Nx  = lFa%Nx(:,:,1)
      CALL CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, pSl, lR, lK)

!     Inertia and body forces (mass) contribution
      DO g=1, lFa%nG
         CALL GNNB(lFa, e, g, nV)
         Jac = SQRT(NORM(nV))
         nV  = nV / Jac
         w   = lFa%w(g)*Jac
         N   = lFa%N(:,g)
         CALL CMM_MASS(eNoN, w, N, al, bfl, lR, lK)
      END DO

!     Now doing the assembly part
#ifdef WITH_TRILINOS
      IF (useTrilinosAssemAndLS) THEN
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
!     Contribution from internal stresses to residue and tangents
      SUBROUTINE CMM_STIFFNESS(eNoN, Nx, xl, dl, pS0l, pSl, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: Nx(2,eNoN), xl(3,eNoN), dl(tDof,eNoN),
     2   pS0l(6)
      REAL(KIND=8), INTENT(INOUT) :: pSl(6), lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      REAL(KIND=8), PARAMETER :: kT = 0.833333333333D0

      INTEGER :: a, b, i, j, k
      REAL(KIND=8) :: elM, nu, ht, lam, mu, afl, Jac, Jaci, T1, nV(3),
     2   gC(3,2), gCi(3,2), thet(3,3), thetT(3,3), phi(9,9), xloc(3,3),
     3   ul(9), Bm(5,9), BmT(9,5), Dm(5,5), DBm(5,9), pSm(3,3), S0l(5),
     4   Sl(5), BtS(9), Ke(9,9)

      elM = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu  = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      ht  = eq(cEq)%dmn(cDmn)%prop(shell_thickness)

      lam = elM/(1D0 - nu*nu)
      mu  = 0.5D0*elM/(1D0 + nu)

      afl = eq(cEq)%af*eq(cEq)%beta*dt*dt
      i   = eq(cEq)%s
      j   = i + 1
      k   = j + 1

      CALL GNNS(eNoN, Nx, xl, nV, gC, gCi)
      Jac   = SQRT(NORM(nV))
      Jaci  = 1.0D0/Jac
      nV(:) = nV(:) * Jaci

!     Rotation matrix
      thet(1,:) = gC(:,1)/SQRT(NORM(gC(:,1)))
      thet(3,:) = nV(:)
      thet(2,1) = thet(3,2)*thet(1,3) - thet(3,3)*thet(1,2)
      thet(2,2) = thet(3,3)*thet(1,1) - thet(3,1)*thet(1,3)
      thet(2,3) = thet(3,1)*thet(1,2) - thet(3,2)*thet(1,1)
      thetT = TRANSPOSE(thet)

!     Define phi matrix
      phi = 0D0
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
      pSm = 0D0
      pSm(1,1) = pS0l(1)
      pSm(2,2) = pS0l(2)
      pSm(3,3) = pS0l(3)
      pSm(1,2) = pS0l(4)
      pSm(1,3) = pS0l(5)
      pSm(2,3) = pS0l(6)

      pSm(2,1) = pSm(1,2)
      pSm(3,1) = pSm(1,3)
      pSm(3,2) = pSm(2,3)

      pSm = MATMUL(pSm,  thetT)
      pSm = MATMUL(thet,   pSm)

      S0l(1) = pSm(1,1)
      S0l(2) = pSm(2,2)
      S0l(3) = pSm(1,2)
      S0l(4) = pSm(1,3)
      S0l(5) = pSm(2,3)

!     B matrix
      Bm = 0D0
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
      Dm  = 0D0
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
      Sl = 0D0
      DO a=1, 9
         Sl(1) = Sl(1) + DBm(1,a)*ul(a)
         Sl(2) = Sl(2) + DBm(2,a)*ul(a)
         Sl(3) = Sl(3) + DBm(3,a)*ul(a)
         Sl(4) = Sl(4) + DBm(4,a)*ul(a)
         Sl(5) = Sl(5) + DBm(5,a)*ul(a)
      END DO

!     Prestress is updated with new stress and pulled to global frame
      pSl = 0D0
      IF (pstEq) THEN
         pSm(:,:) = 0D0
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
      BtS = 0D0
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
      T1 = ht*Jac/2.0D0
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
!     Contribution from inertia and body forces to residue and tangent
!     matrix
      SUBROUTINE CMM_MASS(eNoN, w, N, al, bfl, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), al(tDof,eNoN), bfl(3,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN),
     2   lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j, k
      REAL(KIND=8) rho, ht, wl, am, T1, f(3), ud(3)

      rho  = eq(cEq)%dmn(cDmn)%prop(solid_density)
      ht   = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      f(1) = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2) = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3) = eq(cEq)%dmn(cDmn)%prop(f_z)

      wl   = w * ht * rho
      am   = eq(cEq)%am
      i    = eq(cEq)%s
      j    = i + 1
      k    = j + 1

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
!     applied is not a usual follower pressure load aka load does not
!     follow deformation (acceptable for small strains - CMM)
      SUBROUTINE BCMMi(eNoN, idof, w, N, Nx, xl, tfl, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN, idof
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN), xl(3,eNoN),
     2   tfl(idof,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER a
      REAL(KIND=8) :: wl, nV(3), gC(3,2)

      REAL(KIND=8), ALLOCATABLE :: tfn(:)

!     Get traction vector
      ALLOCATE(tfn(idof))
      tfn = 0D0
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

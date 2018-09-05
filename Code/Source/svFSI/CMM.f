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
!     This routine constructs the wall stiffness contributions to the
!     LHS matrix and RHS residual using details from Alberto's thesis
      SUBROUTINE CMM_STIFFNESS(lFa, e, xl, dl, S0, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER, INTENT(IN) :: e
      REAL(KIND=8), INTENT(IN) :: xl(nsd,lFa%eNoN), dl(tDof,lFa%eNoN),
     2 S0(nstd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,lFa%eNoN),
     2   lK(dof*dof,lFa%eNoN,lFa%eNoN)

      INTEGER a, i, j, k, b, k1, k2
      REAL(KIND=8) x_1(nsd), x_2(nsd), nV(nsd), x_cont(nsd,2),
     2   x_local(nsd, lFa%eNoN), thet(nsd,nsd), phiT(nsd*nsd, nsd*nsd),
     3   area, elM, nu, lambda, mu, Ts, thick, trans_stab,
     4   B_mat(5, 9), D_mat(5, 5), BT_mat(9, 5), K_temp(9,9), R_temp(9),
     5   phi(nsd*nsd, nsd*nsd), Jac, d_local(nsd, lFa%eNoN),
     6   S_global(nsd, nsd), S_local(nsd, nsd), s_res(9),
     7   thetT(nsd,nsd)

      elM   = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu    = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      thick = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      trans_stab = (5D0/6D0)
      i     = eq(cEq)%s
      j     = i + 1
      k     = j + 1

      lambda = elM/(1D0 + nu)/(1D0 - nu)
      mu     = elM/2D0/(1D0 + nu)
      Ts     = eq(cEq)%af*eq(cEq)%beta*dt*dt

      ! First, calculate the area using the relation area = 0.5 * mag(x_1 x x_2)
      ! where x_1 and x_2 are vectors originating from a vertex of the triangle
      x_1 = xl(:,2) - xl(:,1)
      x_2 = xl(:,3) - xl(:,1)
      x_cont(:,1) = x_1
      x_cont(:,2) = x_2
      nV = CROSS(x_cont)

      Jac = SQRT(NORM(x_1))
      x_1 = x_1/Jac
      Jac = SQRT(NORM(nV))
      nV = nV/Jac

      x_cont(:,1) = nV
      x_cont(:,2) = x_1

      x_2 = CROSS(x_cont)
      Jac = SQRT(NORM(x_2))
      x_2 = x_2/Jac

      thet(1,:) = x_1
      thet(2,:) = x_2
      thet(3,:) = nV

      DO a=1,lFa%eNoN
        x_local(1,a) = DOT_PRODUCT(xl(:,a), thet(1,:))
        x_local(2,a) = DOT_PRODUCT(xl(:,a), thet(2,:))
        x_local(3,a) = DOT_PRODUCT(xl(:,a), thet(3,:))

        d_local(1,a) = DOT_PRODUCT(dl(:,a), thet(1,:))
        d_local(2,a) = DOT_PRODUCT(dl(:,a), thet(2,:))
        d_local(3,a) = DOT_PRODUCT(dl(:,a), thet(3,:))
      END DO

      Jac = (x_local(1,1)-x_local(1,3))*(x_local(2,2)-x_local(2,3)) -
     2      (x_local(2,1)-x_local(2,3))*(x_local(1,2)-x_local(1,3))
      area = 5D-1 * Jac

      phi = 0D0
      phi(1:3,1:3) = thet
      phi(4:6,4:6) = thet
      phi(7:9,7:9) = thet
      phiT = TRANSPOSE(phi)

      ! Forming the B matrices
      B_mat = 0D0
      B_mat(1,1) = x_local(2,2) - x_local(2,3)
      B_mat(2,2) = x_local(1,3) - x_local(1,2)
      B_mat(3,1) = B_mat(2,2)
      B_mat(3,2) = B_mat(1,1)
      B_mat(4,3) = B_mat(1,1)
      B_mat(5,3) = B_mat(2,2)
      B_mat(1,4) = x_local(2,3) - x_local(2,1)
      B_mat(2,5) = x_local(1,1) - x_local(1,3)
      B_mat(3,4) = B_mat(2,5)
      B_mat(3,5) = B_mat(1,4)
      B_mat(4,6) = B_mat(1,4)
      B_mat(5,6) = B_mat(2,5)
      B_mat(1,7) = x_local(2,1) - x_local(2,2)
      B_mat(2,8) = x_local(1,2) - x_local(1,1)
      B_mat(3,7) = B_mat(2,8)
      B_mat(3,8) = B_mat(1,7)
      B_mat(4,9) = B_mat(1,7)
      B_mat(5,9) = B_mat(2,8)

      B_mat = 5D-1*(1/area)*B_mat
      BT_mat = TRANSPOSE(B_mat)

      ! Forming the D Matrix
      D_mat = 0D0
      D_mat(1,1) = 1D0
      D_mat(1,2) = nu
      D_mat(2,1) = nu
      D_mat(2,2) = 1D0
      D_mat(3,3) = 5D-1 * (1D0 - nu)
      D_mat(4,4) = 5D-1 * trans_stab * (1D0 - nu)
      D_mat(5,5) = 5D-1 * trans_stab * (1D0 - nu)

      D_mat = lambda*D_mat

      ! Form the element stiffness
      BT_mat = MATMUL(BT_mat, D_mat)
      K_temp = MATMUL(BT_mat,B_mat)

      ! Rotate back to the global coordinates
      K_temp = MATMUL(phiT,K_temp)
      K_temp = MATMUL(K_temp,phi)

      ! Multiply by area and thickness
      K_temp = thick*area*K_temp

      ! Multiply stiffness by local displacements to get residual
      R_temp = 0D0
      DO a=1, lFa%eNoN
        k1 = (a-1)*lFa%eNoN + 1

        R_temp = R_temp + MATMUL(K_temp(:,k1:k1+2),dl(i:k,a))
      END DO

      ! Multiply K_temp by af * beta * dt^2
      K_temp = K_temp * Ts

      ! Add in the pre-stress terms
      S_global = 0D0
      S_global(1,1) = S0(1)
      S_global(2,2) = S0(2)
      S_global(3,3) = S0(3)
      S_global(1,2) = S0(4)
      S_global(2,1) = S0(4)
      S_global(1,3) = S0(5)
      S_global(3,1) = S0(5)
      S_global(2,3) = S0(6)
      S_global(3,2) = S0(6)

      ! Rotate this into the local coordinate system
      thetT = TRANSPOSE(thet)
      S_local = 0D0
      S_local = MATMUL(thet, S_global)
      S_local = MATMUL(S_local, thetT)

      ! Apply the components of the local pre-stress to the residual
      s_res = 0D0
      s_res(1) = B_mat(1,1)*S_local(1,1) + B_mat(2,2)*S_local(1,2)
      s_res(2) = B_mat(1,1)*S_local(1,2) + B_mat(2,2)*S_local(2,2)
      s_res(3) = B_mat(1,1)*S_local(1,3) + B_mat(2,2)*S_local(3,2)
      s_res(4) = B_mat(1,4)*S_local(1,1) + B_mat(2,5)*S_local(1,2)
      s_res(5) = B_mat(1,4)*S_local(1,2) + B_mat(2,5)*S_local(2,2)
      s_res(6) = B_mat(1,4)*S_local(1,3) + B_mat(2,5)*S_local(3,2)
      s_res(7) = B_mat(1,7)*S_local(1,1) + B_mat(2,8)*S_local(1,2)
      s_res(8) = B_mat(1,7)*S_local(1,2) + B_mat(2,8)*S_local(2,2)
      s_res(9) = B_mat(1,7)*S_local(1,3) + B_mat(2,8)*S_local(3,2)

      ! Rotate these local residuals back to the global coordinate system
      s_res = MATMUL(phiT, s_res)
      s_res = s_res*area*thick

      ! Add these local element contributions to the elements stiffness
      ! and residual
      DO a=1, lFa%eNoN
        k1 = (a-1)*lFa%eNoN + 1
        lR(1,a) = lR(1,a) + R_temp(k1) + s_res(k1)
        lR(2,a) = lR(2,a) + R_temp(k1+1) + s_res(k1+1)
        lR(3,a) = lR(3,a) + R_temp(k1+2) + s_res(k1+2)
        DO b=1, lFa%eNoN
          k2 = (b-1)*lFa%eNoN + 1
          lK(1,a,b) = lK(1,a,b) + K_temp(k1,k2)
          lK(2,a,b) = lK(2,a,b) + K_temp(k1,k2+1)
          lK(3,a,b) = lK(3,a,b) + K_temp(k1,k2+2)
          lK(dof+1,a,b) = lK(dof+1,a,b) + K_temp(k1+1,k2)
          lK(dof+2,a,b) = lK(dof+2,a,b) + K_temp(k1+1,k2+1)
          lK(dof+3,a,b) = lK(dof+3,a,b) + K_temp(k1+1,k2+2)
          lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + K_temp(k1+2,k2)
          lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + K_temp(k1+2,k2+1)
          lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + K_temp(k1+2,k2+2)
        END DO
      END DO

      RETURN
      END SUBROUTINE CMM_STIFFNESS
!-----------------------------------------------------------------------
!     This subroutine computes the mass matrix contributions to the LHS
!     matrix and RHS residual
      PURE SUBROUTINE CMM_MASS(eNoN, w, N, Nx, phi, al, yl, dl, nV, lR,
     2 lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: w, N(eNoN), Nx(nsd-1,eNoN),
     2              al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
     3              nV(nsd), phi(nsd*nsd, nsd*nsd)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,eNoN),
     2 lK(dof*dof,eNoN,eNoN)

      INTEGER a, b, i, j, k, k1, k2
      REAL(KIND=8) NxdNx, rho, elM, nu, lambda, mu, Ts, Tfl, amd_s, ws,
     2             wf, amd_f, thick, trans_stab, ed(6), ud(nsd), f(nsd),
     3             udn, hc(nsd), y(nsd), y_n(nsd), R_temp(nsd*eNoN),
     4             K_temp(nsd*eNoN, nsd*eNoN), phiT(nsd*nsd,nsd*nsd)

      rho   = eq(cEq)%dmn(cDmn)%prop(solid_density)
      elM   = eq(cEq)%dmn(cDmn)%prop(elasticity_modulus)
      nu    = eq(cEq)%dmn(cDmn)%prop(poisson_ratio)
      thick = eq(cEq)%dmn(cDmn)%prop(shell_thickness)
      f(1)  = eq(cEq)%dmn(cDmn)%prop(f_x)
      f(2)  = eq(cEq)%dmn(cDmn)%prop(f_y)
      f(3)  = eq(cEq)%dmn(cDmn)%prop(f_z)
      i     = eq(cEq)%s
      j     = i + 1
      k     = j + 1

      lambda = elM/(1D0 + nu)/(1D0 - nu)
      mu     = elM/2D0/(1D0 + nu)
      Ts     = eq(cEq)%af*eq(cEq)%beta*dt*dt
      amd_s  = eq(cEq)%am/Ts*rho
      ws     = w*Ts

      trans_stab = (5D0/6D0)

      ed = 0D0
      ud = -f
      DO a=1, eNoN
          ! Function time derivatives
          ud(1) = ud(1) + N(a)*al(i,a)
          ud(2) = ud(2) + N(a)*al(j,a)
          ud(3) = ud(3) + N(a)*al(k,a)
      END DO

      K_temp = 0D0
      R_temp = 0D0

      ! Add in Mass Matrix contributions
      DO a=1, eNoN
         k1 = (a-1)*eNoN + 1
         R_temp(k1) = R_temp(k1) + w*thick*(rho*N(a)*ud(1))
         R_temp(k1+1) = R_temp(k1+1) + w*thick*(rho*N(a)*ud(2))
         R_temp(k1+2) = R_temp(k1+2) + w*thick*(rho*N(a)*ud(3))
         DO b=1, eNoN
            k2 = (b-1)*eNoN + 1
            K_temp(k1,k2) = K_temp(k1,k2) +
     2         ws*thick*(amd_s*N(a)*N(b))

            K_temp(k1+1,k2+1) = K_temp(k1+1,k2+1) +
     2         ws*thick*(amd_s*N(a)*N(b))

            K_temp(k1+2,k2+2) = K_temp(k1+2,k2+2) +
     2         ws*thick*(amd_s*N(a)*N(b))
         END DO
      END DO

      ! Add contributions from this Gauss point to the element
      ! Stiffness and residual
      DO a=1, eNoN
          k1 = (a-1)*eNoN + 1
          lR(1,a) = lR(1,a) + R_temp(k1)
          lR(2,a) = lR(2,a) + R_temp(k1+1)
          lR(3,a) = lR(3,a) + R_temp(k1+2)
          DO b=1, eNoN
             k2 = (b-1)*eNoN + 1
             lK(1,a,b) = lK(1,a,b) + K_temp(k1,k2)
             lK(2,a,b) = lK(2,a,b) + K_temp(k1,k2+1)
             lK(3,a,b) = lK(3,a,b) + K_temp(k1,k2+2)
             lK(dof+1,a,b) = lK(dof+1,a,b) + K_temp(k1+1,k2)
             lK(dof+2,a,b) = lK(dof+2,a,b) + K_temp(k1+1,k2+1)
             lK(dof+3,a,b) = lK(dof+3,a,b) + K_temp(k1+1,k2+2)
             lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + K_temp(k1+2,k2)
             lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + K_temp(k1+2,k2+1)
             lK(2*dof+3,a,b) = lK(2*dof+3,a,b) + K_temp(k1+2,k2+2)
          END DO
      END DO

      RETURN
      END SUBROUTINE CMM_MASS
!-----------------------------------------------------------------------
!     Filter CMM displacements on the walls with CMM type BCs
      SUBROUTINE CMM_DISPF(lEq, Dg)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      REAL(KIND=8), INTENT(INOUT) :: Dg(tDof,tnNo)

      INTEGER a, Ac, iBc, iM, iFa
      REAL(KIND=8), ALLOCATABLE :: tmpR(:,:)

      ALLOCATE(tmpR(nsd,tnNo))
      tmpR = 0D0
      DO iBc=1, lEq%nBc
         IF (.NOT.BTEST(lEq%bc(iBc)%bType,bType_CMM)) CYCLE
         iFa = lEq%bc(iBc)%iFa
         iM  = lEq%bc(iBc)%iM
         DO a=1, msh(iM)%fa(iFa)%nNo
            Ac = msh(iM)%fa(iFa)%gN(a)
            tmpR(:,Ac) = Dg(1:nsd,Ac)
         END DO
      END DO

      DO a=1, tnNo
         Dg(1:nsd,a) = tmpR(:,a)
      END DO
      DEALLOCATE(tmpR)

      RETURN
      END SUBROUTINE CMM_DISPF
!####################################################################

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
!-----------------------------------------------------------------------
!
!     This routine embodies formulation for solving electrophysiology
!     model equation using operator-splitting method. This routine
!     particularly solves diffusive propagation of action potential.
!
!-----------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_CEP(lM, Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, insd, eNoN, cPhys, iFn, nFn
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   dl(:,:), fN(:,:), N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN

      insd = nsd
      IF (lM%lFib) insd = 1

      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     CEP: dof = 1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), fN(nsd,nFn), N(eNoN), Nx(insd,eNoN),
     3   lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_CEP) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         fN = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            al(:,a) = Ag(:,Ac)
            yl(:,a) = Yg(:,Ac)
            dl(:,a) = Dg(:,Ac)
         END DO

         IF (ALLOCATED(lM%fN)) THEN
            DO iFn=1, nFn
               fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
            END DO
         END IF

!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, insd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g) * Jac
            N = lM%N(:,g)

            IF (insd .EQ. 3) THEN
               CALL CEP3D(eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK)

            ELSE IF (insd .EQ. 2) THEN
               CALL CEP2D(eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK)

            ELSE IF (insd .EQ. 1) THEN
               CALL CEP1D(eNoN, insd, w, N, Nx, al, yl, lR, lK)

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

      DEALLOCATE(ptr, xl, al, yl, dl, fN, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_CEP
!####################################################################
!     This is for solving 3D electrophysiology diffusion equations
      SUBROUTINE CEP3D (eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(3,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), fN(3,nFn)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i
      REAL(KIND=RKIND) :: T1, amd, wl, Diso, Dani(nFn), Vd, Vx(3),
     2   F(3,3), C(3,3), Ci(3,3), Jac, fl(3,nFn), Lf(nFn), D(3,3),
     3   DVx(3), DNx(3,eNoN)

      IF (nFn .LT. eq(cEq)%dmn(cDmn)%cep%nFn) err =
     2   " No. of anisotropic conductivies exceed mesh fibers"

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      DO i=1, nFn
         IF (i .LE. eq(cEq)%dmn(cDmn)%cep%nFn) THEN
            Dani(i) = eq(cEq)%dmn(cDmn)%cep%Dani(i)
         ELSE
            Dani(i) = Dani(i-1)
         END IF
      END DO

!     Compute the isotropic part of diffusion tensor based on spatial
!     isotropy for electromechanics. This models stretch induced changes
!     in conduction velocities
      Lf(:) = 1._RKIND
      IF (ecCpld) THEN
!        Get the displacement degrees of freedom
         DO a=1, nEq
            IF (eq(a)%phys .EQ. phys_struct .OR.
     2          eq(a)%phys .EQ. phys_ustruct) THEN
               i = eq(a)%s
               EXIT
            END IF
         END DO

!        Compute deformation gradient tensor
         F(:,:) = 0._RKIND
         F(1,1) = 1._RKIND
         F(2,2) = 1._RKIND
         F(3,3) = 1._RKIND
         DO a=1, eNoN
            F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
            F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
            F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
            F(2,1) = F(2,1) + Nx(1,a)*dl(i+1,a)
            F(2,2) = F(2,2) + Nx(2,a)*dl(i+1,a)
            F(2,3) = F(2,3) + Nx(3,a)*dl(i+1,a)
            F(3,1) = F(3,1) + Nx(1,a)*dl(i+2,a)
            F(3,2) = F(3,2) + Nx(2,a)*dl(i+2,a)
            F(3,3) = F(3,3) + Nx(3,a)*dl(i+2,a)
         END DO
!        Jacobian
         Jac = MAT_DET(F, 3)

!        Compute the Cauchy-Green tensor and its inverse
         C  = MATMUL(TRANSPOSE(F), F)
         Ci = MAT_INV(C, 3)

!        Compute fiber stretch
         DO i=1, nFn
            Lf(i)   = SQRT(NORM(fN(:,i), MATMUL(C, fN(:,i))))
            fl(:,i) = fN(:,i) / Lf(i)
         END DO

!        Diffusion tensor - spatial isotropy
         Diso    = Diso * Jac
         Dani(:) = Dani(:) * Jac
         D(:,:)  = Diso * Ci(:,:)
      ELSE
         D(:,:)  = 0._RKIND
         D(1,1)  = Diso
         D(2,2)  = Diso
         D(3,3)  = Diso
         fl(:,:) = fN(:,:)
      END IF

!     Compute anisotropic components of diffusion tensor
      DO i=1, nFn
         D(1,1) = D(1,1) + Dani(i)*fl(1,i)*fl(1,i)
         D(1,2) = D(1,2) + Dani(i)*fl(1,i)*fl(2,i)
         D(1,3) = D(1,3) + Dani(i)*fl(1,i)*fl(3,i)

         D(2,1) = D(2,1) + Dani(i)*fl(2,i)*fl(1,i)
         D(2,2) = D(2,2) + Dani(i)*fl(2,i)*fl(2,i)
         D(2,3) = D(2,3) + Dani(i)*fl(2,i)*fl(3,i)

         D(3,1) = D(3,1) + Dani(i)*fl(3,i)*fl(1,i)
         D(3,2) = D(3,2) + Dani(i)*fl(3,i)*fl(2,i)
         D(3,3) = D(3,3) + Dani(i)*fl(3,i)*fl(3,i)
      END DO

      i  = eq(cEq)%s
      Vd = 0._RKIND
      Vx = 0._RKIND
      DO a=1, eNoN
         Vd = Vd + N(a)*al(i,a)

         Vx(1) = Vx(1) + Nx(1,a)*yl(i,a)
         Vx(2) = Vx(2) + Nx(2,a)*yl(i,a)
         Vx(3) = Vx(3) + Nx(3,a)*yl(i,a)

         DNx(1,a) = D(1,1)*Nx(1,a) + D(1,2)*Nx(2,a) + D(1,3)*Nx(3,a)
         DNx(2,a) = D(2,1)*Nx(1,a) + D(2,2)*Nx(2,a) + D(2,3)*Nx(3,a)
         DNx(3,a) = D(3,1)*Nx(1,a) + D(3,2)*Nx(2,a) + D(3,3)*Nx(3,a)
      END DO

      DVx(1) = D(1,1)*Vx(1) + D(1,2)*Vx(2) + D(1,3)*Vx(3)
      DVx(2) = D(2,1)*Vx(1) + D(2,2)*Vx(2) + D(2,3)*Vx(3)
      DVx(3) = D(3,1)*Vx(1) + D(3,2)*Vx(2) + D(3,3)*Vx(3)

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Vd + Nx(1,a)*DVx(1) +
     2      Nx(2,a)*DVx(2) + Nx(3,a)*DVx(3))

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b) + Nx(3,a)*DNx(3,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP3D
!-----------------------------------------------------------------------
!     This is for solving 2D electrophysiology diffusion equation
      SUBROUTINE CEP2D (eNoN, nFn, w, N, Nx, al, yl, dl, fN, lR, lK)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(2,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN), fN(2,nFn)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i
      REAL(KIND=RKIND) :: T1, amd, wl, Diso, Dani(nFn), Vd, Vx(2),
     2   F(2,2), C(2,2), Ci(2,2), Jac, fl(2,nFn), Lf(nFn), D(2,2),
     3   DVx(2), DNx(2,eNoN)

      IF (nFn .LT. eq(cEq)%dmn(cDmn)%cep%nFn) err =
     2   "No. of anisotropic conductivies exceed mesh fibers"

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      wl   = w*T1

      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      DO i=1, nFn
         IF (i .LE. eq(cEq)%dmn(cDmn)%cep%nFn) THEN
            Dani(i) = eq(cEq)%dmn(cDmn)%cep%Dani(i)
         ELSE
            Dani(i) = Dani(i-1)
         END IF
      END DO

!     Compute the isotropic part of diffusion tensor based on spatial
!     isotropy for electromechanics. This models stretch induced changes
!     in conduction velocities
      Lf(:) = 1._RKIND
      IF (ecCpld) THEN
         DO a=1, nEq
            IF (eq(a)%phys .EQ. phys_struct .OR.
     2          eq(a)%phys .EQ. phys_ustruct) THEN
               i = eq(a)%s
               EXIT
            END IF
         END DO

!        Compute deformation gradient tensor
         F(:,:) = 0._RKIND
         F(1,1) = 1._RKIND
         F(2,2) = 1._RKIND
         DO a=1, eNoN
            F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
            F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
            F(2,1) = F(2,1) + Nx(1,a)*dl(i+1,a)
            F(2,2) = F(2,2) + Nx(2,a)*dl(i+1,a)
         END DO
!        Jacobian
         Jac = MAT_DET(F, 2)

!        Compute Cauchy-Green tensor and its inverse
         C  = MATMUL(TRANSPOSE(F), F)
         Ci = MAT_INV(C, 2)

!        Compute fiber stretch
         DO i=1, nFn
            Lf(i)   = SQRT(NORM(fN(:,i), MATMUL(C, fN(:,i))))
            fl(:,i) = fN(:,i) / Lf(i)
         END DO

!        Diffusion tensor - spatial isotropy
         Diso    = Diso * Jac
         Dani(:) = Dani(:) * Jac
         D(:,:)  = Diso * C(:,:)
      ELSE
         D(:,:)  = 0._RKIND
         D(1,1)  = Diso
         D(2,2)  = Diso
         fl(:,:) = fN(:,:)
      END IF

      DO i=1, nFn
         D(1,1) = D(1,1) + Dani(i)*fl(1,i)*fl(1,i)
         D(1,2) = D(1,2) + Dani(i)*fl(1,i)*fl(2,i)

         D(2,1) = D(2,1) + Dani(i)*fl(2,i)*fl(1,i)
         D(2,2) = D(2,2) + Dani(i)*fl(2,i)*fl(2,i)
      END DO

      i  = eq(cEq)%s
      Vd = 0._RKIND
      Vx = 0._RKIND
      DO a=1, eNoN
         Vd = Vd + N(a)*al(i,a)

         Vx(1) = Vx(1) + Nx(1,a)*yl(i,a)
         Vx(2) = Vx(2) + Nx(2,a)*yl(i,a)

         DNx(1,a) = D(1,1)*Nx(1,a) + D(1,2)*Nx(2,a)
         DNx(2,a) = D(2,1)*Nx(1,a) + D(2,2)*Nx(2,a)
      END DO

      DVx(1) = D(1,1)*Vx(1) + D(1,2)*Vx(2)
      DVx(2) = D(2,1)*Vx(1) + D(2,2)*Vx(2)

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Vd + Nx(1,a)*DVx(1) +
     2      Nx(2,a)*DVx(2))

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd +
     2         Nx(1,a)*DNx(1,b) + Nx(2,a)*DNx(2,b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP2D
!-----------------------------------------------------------------------
!     This is for solving 1D electrophysiology diffusion equation
!     for Purkinje fibers
      PURE SUBROUTINE CEP1D (eNoN, insd, w, N, Nx, al, yl, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, insd
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(insd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) a, b, i
      REAL(KIND=RKIND) :: T1, amd, wl, Td, Tx, Diso, DNx(eNoN)

      T1   = eq(cEq)%af*eq(cEq)%gam*dt
      amd  = eq(cEq)%am/T1
      Diso = eq(cEq)%dmn(cDmn)%cep%Diso
      i    = eq(cEq)%s
      wl   = w*T1

      Td = 0._RKIND
      Tx = 0._RKIND
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)
         Tx = Tx + Nx(1,a)*yl(i,a)
         DNx(a) = Diso*Nx(1,a)
      END DO

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*(N(a)*Td + Nx(1,a)*Diso*Tx)

         DO b=1, eNoN
            lK(1,a,b) = lK(1,a,b) + wl*(N(a)*N(b)*amd + Nx(1,a)*DNx(b))
         END DO
      END DO

      RETURN
      END SUBROUTINE CEP1D
!#######################################################################
      PURE SUBROUTINE BCEP (eNoN, w, N, h, lR)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), h
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER(KIND=IKIND) :: a
      REAL(KIND=RKIND) f

      f = w*h

!     Here the loop is started for constructing left and right hand side
      DO a=1, eNoN
         lR(1,a) = lR(1,a) + N(a)*f
      END DO

      RETURN
      END SUBROUTINE BCEP
!#######################################################################

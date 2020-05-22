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
!     This is for solving heat equation in a solid (simple diffusion/
!     Laplace equation)
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_HEATS(lM, Ag, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      eNoN = lM%eNoN

!     HEATS: dof = 1
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_heatS) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)  = Ac
            xl(:,a) = x(:,Ac)
            al(:,a) = Ag(:,Ac)
            yl(:,a) = Yg(:,Ac)
         END DO

!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
            END IF
            IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            w = lM%w(g) * Jac
            N = lM%N(:,g)

            IF (nsd .EQ. 3) THEN
               CALL HEATS3D(eNoN, w, N, Nx, al, yl, lR, lK)

            ELSE IF (nsd .EQ. 2) THEN
               CALL HEATS2D(eNoN, w, N, Nx, al, yl, lR, lK)

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

      DEALLOCATE(ptr, xl, al, yl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_HEATS
!####################################################################
!     This is for solving the heat equation in a solid
      PURE SUBROUTINE HEATS3D (eNoN, w, N, Nx, al, yl, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, T1, amd, wl, Td, Tx(nsd), s, rho

      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      rho = eq(cEq)%dmn(cDmn)%prop(solid_density)

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am * rho/T1
      i   = eq(cEq)%s

      wl = w*T1

      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
         Tx(3) = Tx(3) + Nx(3,a)*yl(i,a)
      END DO
      Td = Td * rho

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
!--------------------------------------------------------------------
!     This is for solving the heat equation in a solid
      PURE SUBROUTINE HEATS2D (eNoN, w, N, Nx, al, yl, lR, lK)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), Nx(nsd,eNoN),
     2   al(tDof,eNoN), yl(tDof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(1,eNoN), lK(1,eNoN,eNoN)

      INTEGER(KIND=IKIND) i, a, b
      REAL(KIND=RKIND) nu, T1, amd, wl, Td, Tx(nsd), s, rho

      nu  = eq(cEq)%dmn(cDmn)%prop(conductivity)
      s   = eq(cEq)%dmn(cDmn)%prop(source_term)
      rho = eq(cEq)%dmn(cDmn)%prop(solid_density)

      T1  = eq(cEq)%af*eq(cEq)%gam*dt
      amd = eq(cEq)%am * rho/T1
      i   = eq(cEq)%s

      wl = w*T1

      Td = -s
      Tx = 0._RKIND
      DO a=1, eNoN
         Td = Td + N(a)*al(i,a)

         Tx(1) = Tx(1) + Nx(1,a)*yl(i,a)
         Tx(2) = Tx(2) + Nx(2,a)*yl(i,a)
      END DO
      Td = Td * rho

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
!--------------------------------------------------------------------
      PURE SUBROUTINE BHEATS (eNoN, w, N, h, lR)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: w, N(eNoN), h
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN)

      INTEGER(KIND=IKIND) a

      DO a=1, eNoN
         lR(1,a) = lR(1,a) + w*N(a)*h
      END DO

      RETURN
      END SUBROUTINE BHEATS
!####################################################################

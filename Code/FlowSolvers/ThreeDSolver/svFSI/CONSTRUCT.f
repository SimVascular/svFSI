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
!     This routine calculates the parameters that are required for
!     sub-equations and calculates them based on general variables.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT(lM, al, yl, dl, dol, xl, fNl, ptr, e)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: ptr(lM%eNoN), e
      REAL(KIND=8), INTENT(IN) :: al(tDof,lM%eNoN), yl(tDof,lM%eNoN),
     2   dl(tDof,lM%eNoN), dol(nsd,lM%eNoN), fNl(nsd,lM%eNoN)
      REAL(KIND=8), INTENT(INOUT) :: xl(nsd,lM%eNoN)

      INTEGER g, cPhys, eNoN
      REAL(KIND=8) w, Jac, ksix(nsd,nsd)

      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:), Nx(:,:),
     2   dc(:,:)

      eNoN = lM%eNoN
      ALLOCATE(lK(dof*dof,eNoN,eNoN), lR(dof,eNoN), Nx(nsd,eNoN),
     2   N(eNoN), dc(tDof,eNoN))
!     Updating the current domain
      cDmn = DOMAIN(lM, cEq, e)
!     Updating the shape functions, if neccessary
      IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
!     Setting intial values
      lK    = 0D0
      lR    = 0D0
      cPhys = eq(cEq)%dmn(cDmn)%phys
!     If neccessary correct X, if mesh is moving
      IF (mvMsh) THEN
         IF (cPhys .EQ. phys_mesh) THEN
!     For mesh the reference configuration is the one at beginnign of
!     the time step
            xl = xl + dol
         ELSE IF (cPhys .NE. phys_struct) THEN
!     Otherwise we use the most recent configuration
            xl = xl + dl(nsd+2:2*nsd+1,:)
         END IF
      END IF

      DO g=1, lM%nG
         IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
            CALL GNN(eNoN, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
         END IF
         IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
         w = lM%w(g)*Jac
         N = lM%N(:,g)

         SELECT CASE (cPhys)
         CASE (phys_fluid)
            IF (nsd .EQ. 3) THEN
               CALL FLUID3D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            ELSE
               CALL FLUID2D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            END IF

         CASE (phys_heatS)
            IF (nsd .EQ. 3) THEN
               CALL HEATS3D(eNoN, w, N, Nx, al, yl, lR, lK)
            ELSE
               CALL HEATS2D(eNoN, w, N, Nx, al, yl, lR, lK)
            END IF

         CASE (phys_heatF)
            IF (nsd .EQ. 3) THEN
               CALL HEATF3D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            ELSE
               CALL HEATF2D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            END IF

         CASE (phys_lElas)
            IF (nsd .EQ. 3) THEN
               CALL LELAS3D(eNoN, w, N, Nx, al, dl, lR, lK)
            ELSE
               CALL LELAS2D(eNoN, w, N, Nx, al, dl, lR, lK)
            END IF

         CASE (phys_struct)
            IF (nsd .EQ. 3) THEN
               CALL STRUCT3D(eNoN, w, N, Nx, al, yl, dl, fNl, lR, lK)
            ELSE
               CALL STRUCT2D(eNoN, w, N, Nx, al, yl, dl, fNl, lR, lK)
            END IF

         CASE (phys_mesh)
            w = w/Jac
            IF (nsd .EQ. 3) THEN
               dc(5:7,:) = dl(5:7,:) - dol
               CALL LELAS3D(eNoN, w, N, Nx, al, dc, lR, lK)
            ELSE
               dc(4:5,:) = dl(4:5,:) - dol
               CALL LELAS2D(eNoN, w, N, Nx, al, dc, lR, lK)
            END IF

         CASE (phys_BBO)
            IF (nsd .EQ. 3) THEN
               CALL BBO3D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            ELSE
               CALL BBO2D(eNoN, w, N, Nx, al, yl, ksix, lR, lK)
            END IF
         CASE DEFAULT
            err = "Undefined phys in CONSTRUCT"
         END SELECT
      END DO

#ifdef WITH_TRILINOS
!     Now doing the assembly part
      IF (useTrilinosAssemAndLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#endif
!$OMP CRITICAL
         CALL DOASSEM(eNoN, ptr, lK, lR)
!$OMP END CRITICAL
#ifdef WITH_TRILINOS
      END IF
#endif

      RETURN
      END SUBROUTINE CONSTRUCT

!####################################################################

      SUBROUTINE BCONSTRUCT(lFa, yl, dl, hl, ptr, e)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER, INTENT(IN) :: ptr(lFa%eNoN), e
      REAL(KIND=8), INTENT(IN) :: yl(tDof,lFa%eNoN), hl(lFa%eNoN),
     2   dl(tDof,lFa%eNoN)
      INTEGER a, g, cPhys, eNoN, iM
      REAL(KIND=8) w, h, nV(nsd), y(tDof), Jac

      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:)

      eNoN = lFa%eNoN
      iM   = lFa%iM
      ALLOCATE(lK(dof*dof,eNoN,eNoN), lR(dof,eNoN), N(eNoN))
!     Updating the current domain
      cDmn = DOMAIN(msh(iM), cEq, lFa%gE(e))
!     Updating the shape functions, if neccessary
      IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)
!     If neccessary correct X, if mesh is moving
      lK    = 0D0
      lR    = 0D0
      cPhys = eq(cEq)%dmn(cDmn)%phys

      DO g=1, lFa%nG
         CALL GNNB(lFa, e, g, nV)
         Jac = SQRT(NORM(nV))
         nV  = nV/Jac
         w   = lFa%w(g)*Jac
         N   = lFa%N(:,g)

         h = 0D0
         y = 0D0
         DO a=1, eNoN
            h = h + N(a)*hl(a)
            y = y + N(a)*yl(:,a)
         END DO

         SELECT CASE (cPhys)
         CASE (phys_fluid)
            CALL BFLUID(eNoN, w, N, y, h, nV, lR, lK)
         CASE (phys_heatS)
            CALL BHEATS(eNoN, w, N, h, lR)
         CASE (phys_heatF)
            CALL BHEATF(eNoN, w, N, y, h, nV, lR, lK)
         CASE (phys_lElas)
            CALL BLELAS(eNoN, w, N, h, nV, lR)
         CASE (phys_struct)
            CALL BLELAS(eNoN, w, N, h, nV, lR)
         CASE DEFAULT
            err = "Undefined phys in BCONSTRUCT"
         END SELECT
      END DO

#ifdef WITH_TRILINOS
!     Now doing the assembly part
      IF (useTrilinosAssemAndLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#else
!$OMP CRITICAL
         CALL DOASSEM(eNoN, ptr, lK, lR)
!$OMP END CRITICAL
#endif
#ifdef WITH_TRILINOS
      END IF
#endif
      RETURN
      END SUBROUTINE BCONSTRUCT

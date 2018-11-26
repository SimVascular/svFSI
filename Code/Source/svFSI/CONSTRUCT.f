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

      SUBROUTINE CONSTRUCT(lM, e, eNoN, al, yl, dl, dol, xl, fNl, pS0l,
     2   bfl, ptr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN) :: e, eNoN, ptr(eNoN)
      REAL(KIND=8), INTENT(IN) :: al(tDof,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), dol(nsd,eNoN), fNl(nFn*nsd,eNoN),
     3   pS0l(nstd,eNoN), bfl(nsd,eNoN)
      REAL(KIND=8), INTENT(INOUT) :: xl(nsd,eNoN)

      INTEGER a, g, Ac, cPhys, insd
      REAL(KIND=8) w, Jac, ksix(nsd,nsd), ctime

      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:), Nx(:,:),
     2   dc(:,:), pSl(:), lKd(:,:,:)

      TYPE(fCellType) :: fCell

      insd = nsd
      IF (lM%lFib) insd = 1

      ALLOCATE(lK(dof*dof,eNoN,eNoN), lR(dof,eNoN), Nx(insd,eNoN),
     2   N(eNoN), dc(tDof,eNoN), pSl(nstd))
!     Updating the current domain
      cDmn = DOMAIN(lM, cEq, e)
!     Updating the shape functions, if neccessary
      IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
!     Setting intial values
      lK    = 0D0
      lR    = 0D0
      pSl   = 0D0
      cPhys = eq(cEq)%dmn(cDmn)%phys

      IF (cPhys .EQ. phys_ustruct) THEN
         ALLOCATE(lKd(dof*nsd,eNoN,eNoN))
         lKd = 0D0
      END IF

!     If neccessary correct X, if mesh is moving
      IF (mvMsh) THEN
         IF (cPhys .EQ. phys_mesh) THEN
!     For mesh the reference configuration is the one at beginning of
!     the time step
            xl = xl + dol
         ELSE IF (cPhys .NE. phys_struct) THEN
!     Otherwise we use the most recent configuration
            xl = xl + dl(nsd+2:2*nsd+1,:)
         END IF
      END IF

!     Finite cell integration for Nitsche formulation
      IF (lM%iGC(e) .EQ. 1) THEN
         IF (lM%lShl .OR. lM%lFib) RETURN
         IF (ib%fcFlag) THEN
            ctime = CPUT()
            CALL FC_INIT(fCell, lM, xl)
            CALL FC_SET(fCell, lM, xl, ib%Uo)
            CALL FC_CONSTRUCT(fCell, al, yl, xl, lR, lK)
            CALL DOASSEM(eNoN, ptr, lK, lR)
            ib%callD(1) = ib%callD(1) + CPUT() - ctime
            RETURN
         END IF
      END IF

      IF (cPhys .EQ. phys_shell) THEN
         IF (lM%eType .EQ. eType_TRI) THEN
!        No need for numerical integration for constant strain triangles
            CALL SHELLTRI(lM, e, eNoN, al, yl, dl, xl, fNl, ptr)
            DEALLOCATE(lR, lK, N, Nx, dc)
            RETURN
         ELSE IF (lM%eType .EQ. eType_NRB) THEN
!        For NURBS, perform Gauss integration
            DO g=1, lM%nG
               CALL SHELLNRB(lM, g, eNoN, al, yl, dl, xl, fNl, lR,lK)
            END DO
!        Perform global assembly
            CALL DOASSEM(eNoN, ptr, lK, lR)
            DEALLOCATE(lR, lK, N, Nx, dc)
            RETURN
         END IF
      END IF

      DO g=1, lM%nG
         IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
            CALL GNN(eNoN, insd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
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

         CASE (phys_CMM)
            IF(nsd .EQ. 3) THEN
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

         CASE (phys_struct, phys_preSt)
            IF (nsd .EQ. 3) THEN
               CALL STRUCT3D(eNoN, w, N, Nx, al, yl, dl, bfl, fNl, pS0l,
     2            pSl, lR, lK)
            ELSE
               CALL STRUCT2D(eNoN, w, N, Nx, al, yl, dl, bfl, fNl, pS0l,
     2            pSl, lR, lK)
            END IF

!      Map pSl values to global nodal vector
            IF (cPhys .EQ. phys_preSt) THEN
               DO a=1, eNoN
                  Ac = ptr(a)
                  pSn(:,Ac) = pSn(:,Ac) + w*N(a)*pSl(:)
                  pSa(Ac)   = pSa(Ac)   + w*N(a)
               END DO
            END IF

         CASE (phys_ustruct)
            IF (nsd .EQ. 3) THEN
               CALL USTRUCT3D(eNoN, w, Jac, N, Nx, al, yl, dl, bfl, fNl,
     2            lR, lK, lKd)
            ELSE
               CALL USTRUCT2D(eNoN, w, Jac, N, Nx, al, yl, dl, bfl, fNl,
     2            lR, lK, lKd)
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

         CASE (phys_CEP)
            IF (insd .EQ. 3) THEN
               CALL CEP3D(eNoN, w, N, Nx, al, yl, fNl, lR, lK)

            ELSE IF (insd .EQ. 2) THEN
               CALL CEP2D(eNoN, w, N, Nx, al, yl, fNl, lR, lK)

            ELSE IF (insd .EQ. 1) THEN
               CALL CEP1D(eNoN, insd, w, N, Nx, al, yl, lR, lK)
            END IF

         CASE DEFAULT
            err = "Undefined phys in CONSTRUCT"
         END SELECT
      END DO

      IF (cPhys .EQ. phys_ustruct) THEN
         CALL USTRUCT_DOASSEM(eNoN, ptr, lKd)
      END IF

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

      RETURN
      END SUBROUTINE CONSTRUCT
!####################################################################
      SUBROUTINE BCONSTRUCT(lFa, yl, hl, ptr, e)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER, INTENT(IN) :: ptr(lFa%eNoN), e
      REAL(KIND=8), INTENT(IN) :: yl(tDof,lFa%eNoN), hl(lFa%eNoN)

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
         CASE (phys_fluid, phys_CMM)
            CALL BFLUID(eNoN, w, N, y, h, nV, lR, lK)

         CASE (phys_heatS)
            CALL BHEATS(eNoN, w, N, h, lR)

         CASE (phys_heatF)
            CALL BHEATF(eNoN, w, N, y, h, nV, lR, lK)

         CASE (phys_lElas)
            CALL BLELAS(eNoN, w, N, h, nV, lR)

         CASE (phys_struct)
            CALL BLELAS(eNoN, w, N, h, nV, lR)

         CASE (phys_CEP)
            CALL BCEP(eNoN, w, N, h, lR)

         CASE (phys_preSt)
            CALL BLELAS(eNoN, w, N, h, nV, lR)

         CASE (phys_shell)
            CALL BLELAS(eNoN, w, N, h, nV, lR)

         CASE DEFAULT
            err = "Undefined phys in BCONSTRUCT"
         END SELECT
      END DO

!     Now doing the assembly part
#ifdef WITH_TRILINOS
      IF (useTrilinosAssemAndLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#else
         CALL DOASSEM(eNoN, ptr, lK, lR)
#endif
#ifdef WITH_TRILINOS
      END IF
#endif

      RETURN
      END SUBROUTINE BCONSTRUCT
!####################################################################
!     This subroutines actually calculates the LHS and RHS
!     contributions and adds their contributions for the Coupled
!     Momentum Method (CMM)
      SUBROUTINE CMM_CONSTRUCT(lFa, al, yl, dl, xl, pS0l, ptr, e)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa
      INTEGER, INTENT(IN) :: ptr(lFa%eNon), e
      REAL(KIND=8), INTENT(IN) :: al(tDof,lFa%eNoN),yl(tDof,lFa%eNoN),
     2   dl(tDof,lFa%eNoN), xl(nsd,lFa%eNoN), pS0l(nstd,lFa%eNoN)

      INTEGER g, cPhys, eNoN, iM, a, b
      REAL(KIND=8) w, nV(nsd), Jac, phi(nsd*nsd,nsd*nsd), S0(nstd)

      REAL(KIND=8), ALLOCATABLE :: lK(:,:,:), lR(:,:), N(:), Nx(:,:)

      eNoN = lFa%eNoN
      iM = lFa%iM
      ALLOCATE(lK(dof*dof,eNoN,eNoN), lR(dof,eNoN), N(eNoN),
     2   Nx(nsd-1,eNoN))

!     Updating the current domain
      cDmn = DOMAIN(msh(iM), cEq, e)

!     Updating the shape functions, if neccesary
      IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)

      lK = 0D0
      lR = 0D0
      cPhys = eq(cEq)%dmn(cDmn)%phys

!     Compute the averaged pre-stressed tensor, S0
      S0 = 0D0
      DO a=1, eNoN
         DO b=1, nstd
            S0(b) = S0(b) + pS0l(b,a)/REAL(eNoN, KIND=8)
         END DO
      END DO

!     Finding the local basis vectors for the element, and the local
!     versions of the acceleration, velocity, and displacement
      CALL CMM_STIFFNESS(lFa, e, xl, dl, S0, lR, lK)

!     Add in the mass matrix contribution
      DO g=1, lFa%nG
         CALL GNNB(lFa, e, g, nV)
         Jac = SQRT(NORM(nV))
         nV  = nV/Jac
         w   = lFa%w(g)*Jac
         N   = lFa%N(:,g)
         CALL CMM_MASS(eNoN, w, N, Nx, phi, al, yl, dl, nV, lR, lK)
      END DO

!     Now doing the assembly part
#ifdef WITH_TRILINOS
      IF (useTrilinosAssemAndLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#else
         CALL DOASSEM(eNoN, ptr, lK, lR)
#endif
#ifdef WITH_TRILINOS
      END IF
#endif

      RETURN
      END SUBROUTINE CMM_CONSTRUCT
!####################################################################

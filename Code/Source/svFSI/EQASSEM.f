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
!     This routine assembles the equation on a given mesh.
!
!--------------------------------------------------------------------

      SUBROUTINE GLOBALEQASSEM(lM, Ag, Yg, Dg)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      SELECT CASE (eq(cEq)%phys)
      CASE (phys_fluid)
         CALL CONSTRUCT_FLUID(lM, Ag, Yg)

      CASE (phys_heatF)
         CALL CONSTRUCT_HEATF(lM, Ag, Yg)

      CASE (phys_heatS)
         CALL CONSTRUCT_HEATS(lM, Ag, Yg)

      CASE (phys_lElas)
         CALL CONSTRUCT_LELAS(lM, Ag, Dg)

      CASE (phys_struct)
         CALL CONSTRUCT_dSOLID(lM, Ag, Yg, Dg)

      CASE (phys_ustruct)
         CALL CONSTRUCT_uSOLID(lM, Ag, Yg, Dg)

      CASE (phys_CMM)
         CALL CONSTRUCT_CMM(lM, Ag, Yg, Dg)

      CASE (phys_shell)
         CALL CONSTRUCT_SHELL(lM, Ag, Yg, Dg)

      CASE (phys_FSI)
         CALL CONSTRUCT_FSI(lM, Ag, Yg, Dg)

      CASE (phys_mesh)
         CALL CONSTRUCT_MESH(lM, Ag, Dg)

      CASE (phys_CEP)
         CALL CONSTRUCT_CEP(lM, Ag, Yg, Dg)

      CASE (phys_stokes)
         CALL CONSTRUCT_STOKES(lM, Ag, Yg)

      CASE DEFAULT
         err = " Undefined physics selection for assembly"
      END SELECT

      RETURN
      END SUBROUTINE GLOBALEQASSEM
!####################################################################
!     Construct Neumann BCs
      SUBROUTINE BASSEMNEUBC(lFa, hg, Yg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: hg(tnNo), Yg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, iM, cPhys, eNoN
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), hl(:), yl(:,:), lR(:,:),
     2   lK(:,:,:)

      iM   = lFa%iM
      eNoN = lFa%eNoN
      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys

         ALLOCATE(ptr(eNoN), N(eNoN), hl(eNoN), yl(tDof,eNoN),
     2      lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
         lK = 0._RKIND
         lR = 0._RKIND
         DO a=1, eNoN
            Ac      = lFa%IEN(a,e)
            ptr(a)  = Ac
            yl(:,a) = Yg(:,Ac)
            hl(a)   = hg(Ac)
         END DO

!        Updating the shape functions, if neccessary
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM), lFa, e)

         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, nsd-1, eNoN, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV))
            nV  = nV/Jac
            w   = lFa%w(g)*Jac
            N   = lFa%N(:,g)

            h = 0._RKIND
            y = 0._RKIND
            DO a=1, eNoN
               h = h + N(a)*hl(a)
               y = y + N(a)*yl(:,a)
            END DO

            SELECT CASE (cPhys)
            CASE (phys_fluid)
               CALL BFLUID(eNoN, w, N, y, h, nV, lR, lK)

            CASE (phys_CMM)
               CALL BFLUID(eNoN, w, N, y, h, nV, lR, lK)

            CASE (phys_heatS)
               CALL BHEATS(eNoN, w, N, h, lR)

            CASE (phys_heatF)
               CALL BHEATF(eNoN, w, N, y, h, nV, lR, lK)

            CASE (phys_lElas)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_struct)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_ustruct)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_shell)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_mesh)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_stokes)
               CALL BLELAS(eNoN, w, N, h, nV, lR)

            CASE (phys_CEP)
               CALL BCEP(eNoN, w, N, h, lR)

            CASE DEFAULT
               err = "Undefined phys in BCONSTRUCT"
            END SELECT
         END DO

!        Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
         END IF
#endif
         DEALLOCATE(ptr, N, hl, yl, lR, lK)
      END DO

      RETURN
      END SUBROUTINE BASSEMNEUBC
!--------------------------------------------------------------------
!     For struct/ustruct - construct follower pressure load.
!     We use Nanson's formula to take change in normal direction with
!     deformation into account. Additional calculations based on mesh
!     need to be performed.
      SUBROUTINE BNEUFOLWP(lFa, hg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: hg(tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, iM, cPhys, eNoN, eNoNb
      REAL(KIND=RKIND) w, Jac, xp(nsd), xi(nsd), xi0(nsd), nV(nsd),
     2   ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: hl(:), xl(:,:), dl(:,:), N(:),
     2   Nxi(:,:), Nx(:,:), lR(:,:), lK(:,:,:), lKd(:,:,:)

      iM    = lFa%iM
      eNoN  = msh(iM)%eNoN
      eNoNb = lFa%eNoN
      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys

         ALLOCATE(ptr(eNoN), hl(eNoN), xl(nsd,eNoN), dl(tDof,eNoN),
     2      N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN),
     3      lK(dof*dof,eNoN,eNoN))
         lR = 0._RKIND
         lK = 0._RKIND
         IF (cPhys .EQ. phys_ustruct) THEN
            ALLOCATE(lKd(dof*nsd,eNoN,eNoN))
            lKd = 0._RKIND
         END IF

!        Create local copies
         DO a=1, eNoN
            Ac      = msh(iM)%IEN(a,Ec)
            ptr(a)  = Ac
            hl(a)   = hg(Ac)
            xl(:,a) = x(:,Ac) ! reference config position vector at each node
            dl(:,a) = Dg(:,Ac) ! displacement vector at each node
         END DO

!        Initialize parameteric coordinate for Newton's iterations
         xi0 = 0._RKIND
         DO g=1, msh(iM)%nG
            xi0 = xi0 + msh(iM)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(msh(iM)%nG, KIND=RKIND)

         DO g=1, lFa%nG
            xp = 0._RKIND
            DO a=1, eNoNb
               Ac = lFa%IEN(a,e)
               xp = xp + x(:,Ac)*lFa%N(a,g)
            END DO

            xi = xi0
            CALL GETNNX(msh(iM)%eType, eNoN, xl, msh(iM)%xib,
     2         msh(iM)%Nb, xp, xi, N, Nxi)

            IF (g.EQ.1 .OR. .NOT.msh(iM)%lShpF)
     2         CALL GNN(eNoN, nsd, Nxi, xl, Nx, Jac, ksix)

!           Get a vector (nV) at element "e" and Gauss point
!           "g" of face "lFa" that is the normal weigthed by Jac, i.e.
!           Jac = SQRT(NORM(n)). 
            CALL GNNB(lFa, e, g, nsd-1, eNoNb, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV)) ! Extract Jacobian
!           AB 5/11/22: I believe this is the Jacobian of the mapping from parent
!           surface element to ref configuration surface element, so this encodes 
!           the area of the ref configuration surface element
            nV  = nV / Jac ! Normalize element surface normal
            w   = lFa%w(g)*Jac ! Scale Gauss point weights by Jacobian

            IF (cPhys .EQ. phys_ustruct) THEN
               IF (nsd .EQ. 3) THEN
                  CALL BUSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK,
     2               lKd)
               ELSE
                  CALL BUSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK,
     2               lKd)
               END IF
            ELSE IF (cPhys .EQ. phys_struct) THEN
               IF (nsd .EQ. 3) THEN
                  CALL BSTRUCT3D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
               ELSE
                  CALL BSTRUCT2D(eNoN, w, N, Nx, dl, hl, nV, lR, lK)
               END IF
            END IF
         END DO

!        Now doing the assembly part
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
         ELSE
#endif
            IF (cPhys .EQ. phys_ustruct) THEN
               CALL USTRUCT_DOASSEM(eNoN, ptr, lKd, lK, lR)
               DEALLOCATE(lKd)

            ELSE IF (cPhys .EQ. phys_struct) THEN
               CALL DOASSEM(eNoN, ptr, lK, lR)

            END IF
#ifdef WITH_TRILINOS
         END IF
#endif
         DEALLOCATE(ptr, hl, xl, dl, N, Nxi, Nx, lR, lK)
      END DO

      RETURN
      END SUBROUTINE BNEUFOLWP
!####################################################################

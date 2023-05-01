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

      CASE (phys_shell)
         CALL CONSTRUCT_SHELL(lM, Ag, Yg, Dg)

      CASE (phys_CMM)
         CALL CONSTRUCT_CMM(lM, Ag, Yg, Dg)

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
      SUBROUTINE BNEUFOLWP(lBc, lFa, hg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(IN) :: lBc
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
      DO e=1, lFa%nEl   ! Begin loop over elements
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
            xl(:,a) = x(:,Ac)    ! reference config position vector at each node
            dl(:,a) = Dg(:,Ac)   ! displacement vector at each node
         END DO

!        Initialize parameteric coordinate for Newton's iterations
         xi0 = 0._RKIND
         DO g=1, msh(iM)%nG
            xi0 = xi0 + msh(iM)%xi(:,g)
         END DO
         xi0 = xi0 / REAL(msh(iM)%nG, KIND=RKIND)

         DO g=1, lFa%nG ! Begin loop over Gauss point
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

!           Get a vector (nV) at element e and Gauss point g of face lFa that 
!           is the normal weighted by Jac, i.e. Jac = SQRT(NORM(n)), where 
!           NORM(u) gives the SQUARE of the Euclidean norm of u
            CALL GNNB(lFa, e, g, nsd-1, eNoNb, lFa%Nx(:,:,g), nV)
            Jac = SQRT(NORM(nV)) ! Extract Jacobian
!           Note, this is the Jacobian of the mapping from parent surface 
!           element to ref configuration surface element, so this encodes 
!           the area of the ref configuration surface element
            nV  = nV / Jac       ! Normalize element surface normal
            w   = lFa%w(g)*Jac   ! Scale Gauss point weights by Jacobian

!           Calculate residual and tangent contributions due to follower pressure
!           load. These are stored in local arrays lR and lK
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
         END DO   ! End loop over Gauss point

!        Now assemble contributions into global residual and stiffness arrays,
!        R and Val
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
      END DO ! End loop over elements

!     AB 5/16/22:
!     Now update surface integrals involved in coupled/resistance BC contribution to
!     stiffness matrix to reflect deformed geometry. The value of this
!     integral is stored in lhs%face%val. Since we are using
!     the deformed geometry to compute the contribution of the pressure
!     load to the residual vector (i.e. follower pressure), we must also use the 
!     deformed geometry to compute the contribution of the resistance BC to the tangent
!     matrix
      IF (BTEST(lBc%bType, bType_res)) THEN
         CALL FSILSUPD(lBc, lFa, lBc%lsPtr)
      END IF
      

      RETURN
      END SUBROUTINE BNEUFOLWP

! ----------------------------------------------------------------------
!     AB 5/16/22:
!     Update the surface integral involved in the coupled/resistance BC input to the 
!     linear solver to take into account the deformed geometry.
!     This integral is sV = int_Gammat (Na * n_i) (See Moghadam 2013 eq. 27.)
!     This function recomputes this integral and updates the variable 
!     lhs%face%val with the new value, which is eventually used in ADDBCMUL() 
!     to add the resistance BC contribution to the matrix-vector product 
!     of the tangent matrix and an arbitrary vector. 
!     This code was more or less copied from BAFINI.f::FSILSINI(). The
!     major differences is that I call GNNB() with the 'n' flag to get 
!     the weighted normal in the current configuration, rather than the 
!     weighted normal in the reference configuration, and that I call a
!     new function FSILS_BC_UPDATE() to update lhs values rather than
!     reallocate and recompute everything.
      SUBROUTINE FSILSUPD(lBc, lFa, lsPtr)   
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(INOUT) :: lsPtr
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) a, e, Ac, g, iM, i, nNo
      REAL(KIND=RKIND) n(nsd)
      LOGICAL :: eDrn

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:)
      REAL(KIND=RKIND), ALLOCATABLE :: sV(:,:), sVl(:,:)
      CHARACTER cfg

      iM  = lFa%iM
      nNo = lFa%nNo
      ALLOCATE(sVl(nsd,nNo), sV(nsd,tnNo), gNodes(nNo))
      DO a=1, nNo
         gNodes(a) = lFa%gN(a)
      END DO

      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (lBc%weakDir) THEN
            lBc%lsPtr = 0
         ELSE
!            lsPtr     = lsPtr + 1
!            lBc%lsPtr = lsPtr
            sVl = 0._RKIND
            eDrn = .FALSE.
            DO i=1, nsd
               IF (lBc%eDrn(i) .NE. 0) THEN
                  eDrn = .TRUE.
                  EXIT
               END IF
            END DO
            IF (eDrn) THEN
               sVl = 1._RKIND
               DO i=1, nsd
                  IF (lBc%eDrn(i) .NE. 0) sVl(i,:) = 0._RKIND
               END DO
            END IF
!           Probably should be _UPDATE, but I don't think this is called anyway
            CALL FSILS_BC_CREATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Dir,
     2         gNodes, sVl)
         END IF
      ELSE IF (BTEST(lBc%bType,bType_Neu)) THEN
!        AB 5/13/22: This is where integrals in Moghadam et al. 
!        eq. 27 are computed.
         IF (BTEST(lBc%bType,bType_res)) THEN ! If resistance BC (or cpl BC)
            sV = 0._RKIND
            DO e=1, lFa%nEl ! Loop over elements on face
               IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM),lFa,e) ! If NURBS
               DO g=1, lFa%nG ! Loop over Gauss point
!                 Get weighted normal vector in current config
                  cfg = 'n'
                  CALL GNNB(lFa, e, g, nsd-1, lFa%eNoN, lFa%Nx(:,:,g),
     2             n, cfg)

                  DO a=1, lFa%eNoN     ! Loop over nodes  in element
                     Ac = lFa%IEN(a,e) ! Extract global nodal index
                     IF (Ac .NE. 0) THEN 
                        sV(:,Ac) = sV(:,Ac) + lFa%N(a,g)*lFa%w(g)*n ! Integral of shape function times weighted normal
                     END IF
                  END DO
               END DO
            END DO
            DO a=1, lFa%nNo
               Ac       = lFa%gN(a)
               sVl(:,a) = sV(:,Ac)
            END DO
!            lsPtr     = lsPtr + 1
!            lBc%lsPtr = lsPtr

!           Fills lhs%face(i) variables, including val if sVl exists
            CALL FSILS_BC_UPDATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Neu,
     2         gNodes, sVl)
         ELSE
            lBc%lsPtr = 0
         END IF
      ELSE IF (BTEST(lBc%bType,bType_trac)) THEN
         lBc%lsPtr = 0
      ELSE IF (BTEST(lBc%bType,bType_CMM)) THEN
!         lsPtr     = lsPtr + 1
!         lBc%lsPtr = lsPtr

         nNo = 0
         DO a=1, lFa%nNo
            IF (ISZERO(lBc%gx(a))) nNo = nNo + 1
         END DO
         DEALLOCATE(gNodes, sVl)
         ALLOCATE(sVl(nsd,nNo), gNodes(nNo))
         sVl  = 0._RKIND

         eDrn = .FALSE.
         DO i=1, nsd
            IF (lBc%eDrn(i) .NE. 0) THEN
               eDrn = .TRUE.
               EXIT
            END IF
         END DO
         IF (eDrn) THEN
            sVl = 1._RKIND
            DO i=1, nsd
               IF (lBc%eDrn(i) .NE. 0) sVl(i,:) = 0._RKIND
            END DO
         END IF

         nNo = 0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (ISZERO(lBc%gx(a))) THEN
               nNo = nNo + 1
               gNodes(nNo) = Ac
            END IF
         END DO

!        Probably should be _UPDATE, but I don't think this is called anyway
         CALL FSILS_BC_CREATE(lhs, lsPtr, nNo, nsd, BC_TYPE_Dir, gNodes,
     2      sVl)
      ELSE
         err = "Unxpected bType in FSILSINI"
      END IF

      DEALLOCATE(sVl, sV, gNodes)

      RETURN
      END SUBROUTINE FSILSUPD
!####################################################################

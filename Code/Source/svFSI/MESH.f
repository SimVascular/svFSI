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
!     This is for constructing mesh equations on the fluid domain for
!     moving mesh/FSI problems.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_MESH(lM, Ag, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, is, ie, eNoN, cPhys
      REAL(KIND=RKIND) w, Jac, ksix(nsd,nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), dl(:,:),
     2   dol(:,:), bfl(:,:), pS0l(:,:), pSl(:), N(:), Nx(:,:), lR(:,:),
     3   lK(:,:,:), lVWP(:,:)

      eNoN = lM%eNoN
      IF (.NOT.mvMsh) err = "Mesh equation is solved for moving mesh/"//
     2   "FSI problems only"

!     MESH: dof = nsd
      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), dl(tDof,eNoN),
     2   dol(nsd,eNoN), bfl(nsd,eNoN), pS0l(nsymd,eNoN), pSl(nsymd),
     3   N(eNoN), Nx(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     4   lVWP(nvwp,eNoN))

!     Start and end DOF
      is = nsd + 2
      ie = 2*nsd + 1

!     Loop over all elements of mesh
      DO e=1, lM%nEl
!        Update domain and proceed if domain phys and eqn phys match
         cDmn  = DOMAIN(lM, cEq, e)
         cPhys = eq(cEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_mesh) CYCLE

!        Update shape functions for NURBS
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Create local copies
         bfl  = 0._RKIND
         pS0l = 0._RKIND
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            ptr(a)   = Ac
            xl(:,a)  = x(:,Ac)
            al(:,a)  = Ag(:,Ac)
            dl(:,a)  = Dg(:,Ac)
            dol(:,a) = Do(is:ie,Ac)
!           Variable wall - SCHWARZ July 2021---------------------------
!           Calculate local wall property
            IF (useVarWall) lVWP(:,a) = vWP0(:,Ac)
!           ------------------------------------------------------------
         END DO

!        For MESH, the reference configuration is the one at the
!        beginning of the time step. Update displacements accordingly
         xl(:,:) = xl(:,:) + dol(:,:)
         dl(is:ie,:) = dl(is:ie,:) - dol(:,:)

!        Gauss integration
         lR = 0._RKIND
         lK = 0._RKIND
         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e
            END IF
            w = lM%w(g)
            N = lM%N(:,g)

            pS0l = 0._RKIND
            IF (nsd .EQ. 3) THEN
               CALL LELAS3D(eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR,
     2            lK,lVWP)

            ELSE IF (nsd .EQ. 2) THEN
               CALL LELAS2D(eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR,
     2            lK)

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

      DEALLOCATE(ptr, xl, al, dl, dol, bfl, pS0l, pSl, N, Nx, lR, lK)

      RETURN
      END SUBROUTINE CONSTRUCT_MESH
!####################################################################

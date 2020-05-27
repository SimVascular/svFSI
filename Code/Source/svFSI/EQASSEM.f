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
         CALL CONSTRUCT_STOKES(lM, Yg)

      CASE DEFAULT
         err = " Undefined physics selection for assembly"
      END SELECT

      RETURN
      END SUBROUTINE GLOBALEQASSEM
!####################################################################
      SUBROUTINE BCONSTRUCT(lFa, hg, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: hg(tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, Ec, iM, cPhys, eNoN
      REAL(KIND=RKIND) w, h, nV(nsd), y(tDof), Jac

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:), hl(:), xl(:,:), yl(:,:),
     2   dl(:,:), lR(:,:), lK(:,:,:), lKd(:,:,:)

      iM = lFa%iM
      DO e=1, lFa%nEl
         Ec    = lFa%gE(e)
         cDmn  = DOMAIN(msh(iM), cEq, Ec)
         cPhys = eq(cEq)%dmn(cDmn)%phys

         IF (cPhys .EQ. phys_ustruct) THEN
!           We use Nanson's formula to take change in normal direction
!           with deformation into account. Additional calculations based
!           on mesh element need to be performed.
            eNoN = msh(iM)%eNoN

            ALLOCATE(ptr(eNoN), hl(eNoN), xl(nsd,eNoN), dl(tDof,eNoN),
     2         lR(dof,eNoN), lK(dof*dof,eNoN,eNoN),
     3         lKd(dof*nsd,eNoN,eNoN))
            lR  = 0._RKIND
            lK  = 0._RKIND
            lKd = 0._RKIND
            DO a=1, msh(iM)%eNoN
               Ac      = msh(iM)%IEN(a,Ec)
               ptr(a)  = Ac
               xl(:,a) = x(:,Ac)
               dl(:,a) = Dg(:,Ac)
               hl(a)   = hg(Ac)
            END DO

            CALL BUSTRUCT(lFa, eNoN, e, ptr, xl, dl, hl, lR, lK, lKd)

            DEALLOCATE(ptr, hl, xl, dl, lR, lK, lKd)
         ELSE
            eNoN = lFa%eNoN

            ALLOCATE(ptr(eNoN), N(eNoN), hl(eNoN), yl(tDof,eNoN),
     2         lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))
            lK = 0._RKIND
            lR = 0._RKIND
            DO a=1, eNoN
               Ac      = lFa%IEN(a,e)
               ptr(a)  = Ac
               yl(:,a) = Yg(:,Ac)
               hl(a)   = hg(Ac)
            END DO

!           Updating the shape functions, if neccessary
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

               CASE (phys_CEP)
                  CALL BCEP(eNoN, w, N, h, lR)

               CASE (phys_shell)
                  CALL BLELAS(eNoN, w, N, h, nV, lR)

               CASE (phys_stokes)
                  CALL BLELAS(eNoN, w, N, h, nV, lR)

               CASE DEFAULT
                  err = "Undefined phys in BCONSTRUCT"
               END SELECT
            END DO

!           Now doing the assembly part
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
         END IF
      END DO

      RETURN
      END SUBROUTINE BCONSTRUCT
!####################################################################

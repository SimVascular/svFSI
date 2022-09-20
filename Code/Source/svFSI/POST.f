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
!     Any specific post processing before sending data to outputs is
!     done here.
!
!--------------------------------------------------------------------

!     This is a routine that performs POST/BPOST on all meshes
      SUBROUTINE ALLPOST(res, lY, lD, outGrp, iEq)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(OUT) :: res(maxnsd,tnNo)
      REAL(KIND=RKIND), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: outGrp, iEq

      INTEGER(KIND=IKIND) a, Ac, iM

      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:), tmpVe(:)

      DO iM=1, nMsh
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))
         IF (outGrp.EQ.outGrp_WSS .OR. outGrp.EQ.outGrp_trac) THEN
            CALL BPOST(msh(iM), tmpV, lY, lD, outGrp)
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               res(:,Ac) = tmpV(:,a)
            END DO

         ELSE IF ((outGrp .EQ. outGrp_J)  .OR.
     2            (outGrp .EQ. outGrp_fS) .OR.
     3            (outGrp .EQ. outGrp_mises) ) THEN
            IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
            ALLOCATE(tmpV(1,msh(iM)%nNo), tmpVe(msh(iM)%nEl))
            tmpV  = 0._RKIND
            tmpVe = 0._RKIND

            IF (msh(iM)%lShl) THEN
               CALL SHLPOST(msh(iM), 1, tmpV, tmpVe, lD, iEq, outGrp)
            ELSE
               CALL TPOST(msh(iM), 1, tmpV, tmpVe, lD, lY, iEq, outGrp)
            END  IF

            res  = 0._RKIND
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               res(1,Ac) = tmpV(1,a)
            END DO
            DEALLOCATE(tmpV, tmpVe)
            ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

         ELSE IF (outGrp .EQ. outGrp_divV) THEN
            IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
            ALLOCATE(tmpV(1,msh(iM)%nNo))
            tmpV = 0._RKIND
            CALL DIVPOST(msh(iM), tmpV, lY, lD, iEq)
            res  = 0._RKIND
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               res(1,Ac) = tmpV(1,a)
            END DO
            DEALLOCATE(tmpV)
            ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

         ELSE
            CALL POST(msh(iM), tmpV, lY, lD, outGrp, iEq)
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               res(:,Ac) = tmpV(:,a)
            END DO
         END IF
      END DO

      RETURN
      END SUBROUTINE ALLPOST
!--------------------------------------------------------------------
!     General purpose routine for post processing outputs.
      SUBROUTINE POST(lM, res, lY, lD, outGrp, iEq)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(OUT) :: res(maxnsd,lM%nNo)
      REAL(KIND=RKIND), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: outGrp, iEq

      LOGICAL FSIeq
      INTEGER(KIND=IKIND) a, Ac, e, i, j, eNoN, g, insd
      REAL(KIND=RKIND) rho, kappa, w, Jac, ksix(nsd,nsd), lRes(maxnsd),
     2   q(nsd), u(nsd), p, T, ux(nsd,nsd), gam, mu, mu_s
      COMPLEX(KIND=CXKIND) :: eig(nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), xl(:,:), yl(:,:),
     2   Nx(:,:), N(:)

      FSIeq = .FALSE.
      IF (eq(iEq)%phys .EQ. phys_FSI) FSIeq = .TRUE.

!     Since energy flux can be directly calculated from nodal values.
!     Note that if there are more than one domain, we need to rely on
!     element based calculations
      IF (outGrp.EQ.outGrp_eFlx .AND. .NOT.ALLOCATED(dmnId)) THEN
         rho = eq(iEq)%dmn(1)%prop(fluid_density)
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            p  = lY(nsd+1,Ac)
            u  = lY(1:nsd,Ac)
            res(1:nsd,Ac) = (p + 0.5_RKIND*rho*NORM(u))*u
         END DO
         RETURN
      END IF

!     Other outputs require more calculations
      eNoN  = lM%eNoN
      ALLOCATE (sA(tnNo), sF(maxnsd,tnNo), xl(nsd,eNoN), yl(tDof,eNoN),
     2   Nx(nsd,eNoN), N(eNoN))
      sA   = 0._RKIND
      sF   = 0._RKIND
      lRes = 0._RKIND
      eig  = (0._RKIND, 0._RKIND)
      insd = nsd
      IF (lM%lFib) insd = 1
      DO e=1, lM%nEl
         cDmn = DOMAIN(lM, iEq, e)
         IF (cDmn .EQ. 0) CYCLE
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
!     Finding the norm for all the nodes of this element, including
!     those that don't belong to this face, which will be interpolated
!     from the nodes of the face
         DO a=1, eNoN
            Ac      = lM%IEN(a,e)
            xl(:,a) = x(:,Ac)
            yl(:,a) = lY(:,Ac)
            IF (FSIeq) THEN
               xl(:,a)     = xl(:,a)     + lD(nsd+2:2*nsd+1,Ac)
               yl(1:nsd,a) = yl(1:nsd,a) - lY(nsd+2:2*nsd+1,Ac)
            END IF
         END DO

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, insd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)

            lRes(:) = 0._RKIND

!     Vorticity calculation   ---------------------------------------
            IF (outGrp .EQ. outGrp_vort) THEN
               DO a=1, eNoN
                  IF (nsd .EQ. 2) THEN
                     lRes(3) = lRes(3)+ Nx(1,a)*yl(2,a)- Nx(2,a)*yl(1,a)
                  ELSE
                     lRes(1) = lRes(1)+ Nx(2,a)*yl(3,a)- Nx(3,a)*yl(2,a)
                     lRes(2) = lRes(2)+ Nx(3,a)*yl(1,a)- Nx(1,a)*yl(3,a)
                     lRes(3) = lRes(3)+ Nx(1,a)*yl(2,a)- Nx(2,a)*yl(1,a)
                  END IF
               END DO
!     Vortex Identification Criterion (lamda_ci)
            ELSE IF (outGrp .EQ. outGrp_vortex) THEN
               ux   = 0._RKIND
               DO a=1, eNoN
                  DO i=1, nsd
                     DO j=1, nsd
                        ux(i,j) = ux(i,j) + Nx(i,a)*yl(j,a)
                     END DO
                  END DO
               END DO
               eig = MAT_EIG(ux, nsd)
               lRes(1:nsd) = AIMAG(eig)
               lRes(1) = MAXVAL(lRes(1:nsd))
               lRes(2:maxnsd) = 0._RKIND
!     Energy flux calculation   -------------------------------------
            ELSE IF (outGrp .EQ. outGrp_eFlx) THEN
               rho = eq(iEq)%dmn(cDmn)%prop(fluid_density)
               p   = 0._RKIND
               u   = 0._RKIND
               DO a=1, eNoN
                  p = p + N(a)*yl(nsd+1,a)
                  u = u + N(a)*yl(1:nsd,a)
               END DO
               lRes(1:nsd) = (p + 0.5_RKIND*rho*NORM(u))*u
!     Heat flux calculation   ---------------------------------------
            ELSE IF (outGrp .EQ. outGrp_hFlx) THEN
               kappa = eq(iEq)%dmn(cDmn)%prop(conductivity)
               i     = eq(iEq)%s
               DO a=1, eNoN
               END DO
               IF (eq(iEq)%phys .EQ. phys_heatF) THEN
                  u = 0._RKIND
                  T = 0._RKIND
                  q = 0._RKIND
                  DO a=1, eNoN
                     q = q + Nx(:,a)*yl(i,a)
                     u = u + N(a)*yl(1:nsd,a)
                     T = T + N(a)*yl(i,a)
                  END DO
                  lRes(1:nsd) = u*T - kappa*q
               ELSE
                  q = 0._RKIND
                  DO a=1, eNoN
                     q = q + Nx(:,a)*yl(i,a)
                  END DO
                  lRes(1:nsd) = -kappa*q
               END IF
!     Strain tensor invariants calculation   ------------------------
            ELSE IF (outGrp .EQ. outGrp_stInv) THEN
               ksix = 0._RKIND
               DO a=1, eNoN
                  ksix(1,1) = ksix(1,1) + Nx(1,a)*yl(1,a)
                  ksix(2,1) = ksix(2,1) +(Nx(1,a)*yl(2,a)
     2                                  + Nx(2,a)*yl(1,a))*0.5_RKIND
                  ksix(2,2) = ksix(2,2) + Nx(2,a)*yl(2,a)
                  IF (nsd .EQ. 3) THEN
                     ksix(3,1) = ksix(3,1) +(Nx(1,a)*yl(3,a)
     2                                     + Nx(3,a)*yl(1,a))*0.5_RKIND
                     ksix(3,2) = ksix(3,2) +(Nx(2,a)*yl(3,a)
     2                                     + Nx(3,a)*yl(2,a))*0.5_RKIND
                     ksix(3,3) = ksix(3,3) + Nx(3,a)*yl(3,a)
                  END IF
               END DO
               IF (nsd .EQ. 2) THEN
                  lRes(1) = ksix(1,1) + ksix(2,2)
                  lRes(2) = ksix(1,1)*ksix(2,2) - ksix(2,1)*ksix(2,1)
               ELSE
                  lRes(1) = ksix(1,1) + ksix(2,2) + ksix(3,3)
                  lRes(2) = ksix(1,1)*ksix(2,2) + ksix(2,2)*ksix(3,3)
     2                    + ksix(3,3)*ksix(1,1) - ksix(2,1)*ksix(2,1)
     3                    - ksix(3,1)*ksix(3,1) - ksix(3,2)*ksix(3,2)
                  lRes(3) = ksix(1,1)*ksix(2,2)*ksix(3,3)
     2                    + ksix(2,1)*ksix(3,2)*ksix(3,1)*2._RKIND
     3                    - ksix(1,1)*ksix(3,2)*ksix(3,2)
     4                    - ksix(3,1)*ksix(2,2)*ksix(3,1)
     5                    - ksix(2,1)*ksix(2,1)*ksix(3,3)
               END IF
               lRes = ABS(lRes)
!     Viscosity
            ELSE IF (outGrp .EQ. outGrp_Visc) THEN
               ux = 0._RKIND
               DO a=1, eNoN
                  DO i=1, nsd
                     DO j=1, nsd
                        ux(i,j) = ux(i,j) + Nx(i,a)*yl(j,a)
                     END DO
                  END DO
               END DO

!              Shear rate, gam := (2*e_ij*e_ij)^0.5
               gam = 0._RKIND
               DO i=1, nsd
                  DO j=1, nsd
                     gam = gam + (ux(i,j)+ux(j,i))*(ux(i,j)+ux(j,i))
                  END DO
               END DO
               gam = SQRT(0.5_RKIND*gam)
!              Compute viscosity
               CALL GETVISCOSITY(eq(iEq)%dmn(cDmn), gam, mu, mu_s, mu_s)
               lRes(1) = mu
            ELSE
               err = "Correction is required in POST"
            END IF

!     Mapping Tau into the nodes by assembling it into a local vector
            DO a=1, eNoN
               Ac       = lM%IEN(a,e)
               sA(Ac)   = sA(Ac)   + w*N(a)
               sF(:,Ac) = sF(:,Ac) + w*N(a)*lRes
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac       = lM%gN(a)
         res(:,a) = sF(:,Ac)/sA(Ac)
      END DO

      RETURN
      END SUBROUTINE POST
!--------------------------------------------------------------------
!     General purpose routine for post processing outputs at the
!     faces. Currently this calculates WSS, which is t.n - (n.t.n)n
!     Here t is stress tensor: t = \mu (grad(u) + grad(u)^T)
      SUBROUTINE BPOST(lM, res, lY, lD, outGrp)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(OUT) :: res(maxnsd,lM%nNo)
      REAL(KIND=RKIND), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN) :: outGrp

      LOGICAL FSIeq
      INTEGER(KIND=IKIND) a, Ac, e, Ec, i, j, iEq, iFa, eNoN, g
      REAL(KIND=RKIND) Tdn(nsd), ndTdn, taue(nsd), ux(nsd,nsd), mu, w,
     2   nV(nsd), Jac, ks(nsd,nsd), lRes(maxnsd), p, gam, mu_s
      TYPE(fsType) :: fsP

      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), gnV(:,:),
     2   lnV(:,:), xl(:,:), ul(:,:), pl(:), N(:), Nx(:,:)

      IF (outGrp.NE.outGrp_WSS .AND. outGrp.NE.outGrp_trac) err =
     2   "Invalid output group. Correction is required in BPOST"

      iEq   = 1
      eNoN  = lM%eNoN
      FSIeq = .FALSE.
      IF (eq(iEq)%phys .EQ. phys_FSI) FSIeq = .TRUE.

      ALLOCATE (sA(tnNo), sF(maxnsd,tnNo), xl(nsd,eNoN), ul(nsd,eNoN),
     2   gnV(nsd,tnNo), lnV(nsd,eNoN), N(eNoN), Nx(nsd,eNoN))
      sA   = 0._RKIND
      sF   = 0._RKIND
      gnV  = 0._RKIND
      lRes = 0._RKIND

!     First creating the norm field
      DO iFa=1, lM%nFa
         DO a=1, lM%fa(iFa)%nNo
            Ac        = lM%fa(iFa)%gN(a)
            gnV(:,Ac) = lM%fa(iFa)%nV(:,a)
         END DO
      END DO

!     Update pressure function spaces
      IF (lM%nFs .EQ. 1) THEN
         fsP%nG    = lM%fs(1)%nG
         fsP%eType = lM%fs(1)%eType
         fsP%eNoN  = lM%fs(1)%eNoN
         CALL ALLOCFS(fsP, nsd)
         fsP%w  = lM%fs(1)%w
         fsP%xi = lM%fs(1)%xi
         fsP%N  = lM%fs(1)%N
         fsP%Nx = lM%fs(1)%Nx
      ELSE
         fsP%nG    = lM%fs(1)%nG
         fsP%eType = lM%fs(2)%eType
         fsP%eNoN  = lM%fs(2)%eNoN
         CALL ALLOCFS(fsP, nsd)
         fsP%xi = lM%fs(1)%xi
         DO g=1, fsP%nG
            CALL GETGNN(nsd, fsP%eType, fsP%eNoN, fsP%xi(:,g),
     2         fsP%N(:,g), fsP%Nx(:,:,g))
         END DO
      END IF
      ALLOCATE(pl(fsP%eNoN))

      DO iFa=1, lM%nFa
         DO e=1, lM%fa(iFa)%nEl
            Ec = lM%fa(iFa)%gE(e)
            cDmn = DOMAIN(lM, iEq, Ec)
            IF (cDmn .EQ. 0) CYCLE
            IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, Ec)

!     Finding the norm for all the nodes of this element, including
!     those that don't belong to this face, which will be inerpolated
!     from the nodes of the face
            nV = 0._RKIND
            DO a=1, eNoN
               Ac       = lM%IEN(a,Ec)
               lnV(:,a) = gnV(:,Ac)
               nV       = nV + lnV(:,a)
               xl(:,a)  = x(:,Ac)
               IF (FSIeq) THEN
                  xl(:,a) = xl(:,a)      + lD(nsd+2:2*nsd+1,Ac)
                  ul(:,a) = lY(1:nsd,Ac) - lY(nsd+2:2*nsd+1,Ac)
               ELSE
                  ul(:,a) = lY(1:nsd,Ac)
               END IF
            END DO

            DO a=1, fsP%eNoN
               Ac    = lM%IEN(a,Ec)
               pl(a) = lY(nsd+1,Ac)
            END DO

            nV = nV/lM%fa(iFa)%eNoN
            DO a=1, eNoN
               IF (ISZERO(NORM(lnV(:,a)))) lnV(:,a) = nV
            END DO
            DO g=1, lM%nG
               IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
                  CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ks)
               END IF
               w = lM%w(g)*Jac
               N = lM%N(:,g)

!     Calculating ux = grad(u) and nV at a Gauss point
               ux = 0._RKIND
               nV = 0._RKIND
               DO a=1, eNoN
                  nV = nV + N(a)*lnV(:,a)
                  DO i=1, nsd
                     DO j=1, nsd
                        ux(i,j) = ux(i,j) + Nx(i,a)*ul(j,a)
                     END DO
                  END DO
               END DO

               p = 0._RKIND
               DO a=1, fsP%eNoN
                  p = p + fsP%N(a,g)*pl(a)
               END DO

!              Shear rate, gam := (2*e_ij*e_ij)^0.5
               gam = 0._RKIND
               DO i=1, nsd
                  DO j=1, nsd
                     gam = gam + (ux(i,j)+ux(j,i))*(ux(i,j)+ux(j,i))
                  END DO
               END DO
               gam = SQRT(0.5_RKIND*gam)
!              Compute viscosity
               CALL GETVISCOSITY(eq(iEq)%dmn(cDmn), gam, mu, mu_s, mu_s)

!     Now finding grad(u).n and n.grad(u).n
               Tdn   = 0._RKIND
               ndTdn = 0._RKIND
               DO i=1,nsd
                  DO j=1, nsd
                     Tdn(i) = Tdn(i) + mu*(ux(i,j) + ux(j,i))*nV(j)
                  END DO
                  ndTdn = ndTdn + Tdn(i)*nV(i)
               END DO
               taue = Tdn - ndTdn*nV

               IF (outGrp .EQ. outGrp_WSS) THEN
                  lRes(1:nsd) = -taue
               ELSE IF (outGrp .EQ. outGrp_trac) THEN
                  lRes(1:nsd) = p*nV - Tdn
               END IF

!     Mapping Tau into the nodes by assembling it into a local vector
               DO a=1, eNoN
                  Ac       = lM%IEN(a,Ec)
                  sA(Ac)   = sA(Ac)   + w*N(a)
                  sF(:,Ac) = sF(:,Ac) + w*N(a)*lRes
               END DO
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO iFa=1, lM%nFa
         DO a=1, lM%fa(iFa)%nNo
            Ac = lM%fa(iFa)%gN(a)
            IF (.NOT.ISZERO(sA(Ac))) THEN
               sF(:,Ac) = sF(:,Ac)/sA(Ac)
               sA(Ac)   = 1._RKIND
            END IF
         END DO
      END DO

      res = 0._RKIND
      DO a=1, lM%nNo
         Ac       = lM%gN(a)
         res(:,a) = sF(:,Ac)
      END DO

      RETURN
      END SUBROUTINE BPOST
!####################################################################
!     Routine for post processing stress tensor
      SUBROUTINE TPOST(lM, m, res, resE, lD, lY, iEq, outGrp)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: m, iEq, outGrp
      REAL(KIND=RKIND), INTENT(INOUT) :: res(m,lM%nNo), resE(lM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo), lY(tDof,tnNo)

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, g, Ac, i, j, k, l, cPhys, insd, nFn
      REAL(KIND=RKIND) w, Jac, detF, Je, tmX, ya, Ja, elM, nu, lambda,
     2   mu, p, trS, vmises, xi(nsd), xi0(nsd), xp(nsd), ed(nsymd),
     3   Im(nsd,nsd), F(nsd,nsd), C(nsd,nsd), Eg(nsd,nsd), P1(nsd,nsd),
     4   S(nsd,nsd), sigma(nsd,nsd), Dm(nsymd,nsymd)
      TYPE(fsType) :: fs

      INTEGER, ALLOCATABLE :: eNds(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), dl(:,:), yl(:,:),
     2   tmXl(:), ya_l(:), fN(:,:), resl(:), Nx(:,:), N(:), sA(:),
     3   sF(:,:), sE(:)

      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     For higher order elements, we lower the order of shape functions
!     to compute the quantities at corner nodes of elements. We then
!     use these lower order shape functions to interpolate the values
!     at edge nodes and elements centers (if applicable)
      flag = .FALSE.
      IF (lM%eType.EQ.eType_TRI6  .OR. lM%eType.EQ.eType_QUD8  .OR.
     2    lM%eType.EQ.eType_QUD9  .OR. lM%eType.EQ.eType_TET10 .OR.
     3    lM%eType.EQ.eType_HEX20 .OR. lM%eType .EQ. eType_HEX27) THEN
         flag =.TRUE.
         CALL SETTHOODFS(fs, lM%eType)
      ELSE
         fs%eType = lM%eType
         fs%lShpF = lM%lShpF
         fs%eNoN  = lM%eNoN
         fs%nG    = lM%nG
      END IF
      CALL INITFS(fs, nsd)

      ALLOCATE (sA(tnNo), sF(m,tnNo), sE(lM%nEl), xl(nsd,fs%eNoN),
     2   dl(tDof,fs%eNoN), yl(tDof,fs%eNoN), fN(nsd,nFn), tmXl(fs%eNoN),
     3   ya_l(fs%eNoN), resl(m), Nx(nsd,fs%eNoN), N(fs%eNoN))

      sA   = 0._RKIND
      sF   = 0._RKIND
      sE   = 0._RKIND

      insd = nsd
      IF (lM%lShl) insd = 2

!     Initialize tensor operations
      CALL TEN_INIT(insd)

      IF (lM%lFib) insd = 1

      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_struct .AND.
     2       cPhys .NE. phys_ustruct .AND.
     3       cPhys .NE. phys_lElas) CYCLE

         IF (cPhys .EQ. phys_lElas) THEN
            elM = eq(iEq)%dmn(cDmn)%prop(elasticity_modulus)
            nu  = eq(iEq)%dmn(cDmn)%prop(poisson_ratio)
            lambda = elM*nu/(1._RKIND + nu)/(1._RKIND - 2._RKIND*nu)
            mu     = 0.5_RKIND*elM/(1._RKIND + nu)
         END IF

         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         fN = 0._RKIND
         IF (ALLOCATED(lM%fN)) THEN
            DO l=1, nFn
               fN(:,l) = lM%fN((l-1)*nsd+1:l*nsd,e)
            END DO
         END IF

         dl   = 0._RKIND
         yl   = 0._RKIND
         tmXl = 0._RKIND
         ya_l = 0._RKIND
         DO a=1, fs%eNoN
            Ac = lM%IEN(a,e)
            b  = lM%lN(Ac)
            xl(:,a) = x(:,Ac)
            dl(:,a) = lD(:,Ac)
            yl(:,a) = lY(:,Ac)
            IF (ecCpld) THEN
               IF (ALLOCATED(lM%tmX)) THEN
                  tmXl(a) = lM%tmX(b)
               END IF
               IF (ALLOCATED(ec_Ya)) THEN
                  ya_l(a) = ec_Ya(Ac)
               ELSE
                  ya_l(a) = eq(cEq)%dmn(cDmn)%ec%Ya
               END IF
            END IF
         END DO

!        Gauss integration
         Je = 0._RKIND
         DO g=1, fs%nG
            IF (g.EQ.1 .OR. .NOT.fs%lShpF) THEN
               CALL GNN(fs%eNoN, insd, fs%Nx(:,:,g), xl, Nx, Jac, Im)
            END IF
            w  = fs%w(g)*Jac
            N  = fs%N(:,g)
            Je = Je + w

            Im  = MAT_ID(nsd)
            F   = Im
            tmX = 0._RKIND
            ya  = 0._RKIND
!           Interpolate quantities at the integration point
            DO a=1, fs%eNoN
               IF (nsd .EQ. 3) THEN
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                  F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
                  F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
                  F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
                  F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)
               ELSE
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
               END IF

               tmX = tmX + N(a)*tmXl(a)
               ya  = ya  + N(a)*ya_l(a)
            END DO
            detF = MAT_DET(F, nsd)

            ed = 0._RKIND
            IF (cPhys .EQ. phys_lElas) THEN
               DO a=1, fs%eNoN
                  IF (nsd .EQ. 3) THEN
                     ed(1) = ed(1) + Nx(1,a)*dl(i,a)
                     ed(2) = ed(2) + Nx(2,a)*dl(j,a)
                     ed(3) = ed(3) + Nx(3,a)*dl(k,a)
                     ed(4) = ed(4) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
                     ed(5) = ed(5) + Nx(3,a)*dl(j,a) + Nx(2,a)*dl(k,a)
                     ed(6) = ed(6) + Nx(1,a)*dl(k,a) + Nx(3,a)*dl(i,a)
                  ELSE
                     ed(1) = ed(1) + Nx(1,a)*dl(i,a)
                     ed(2) = ed(2) + Nx(2,a)*dl(j,a)
                     ed(3) = ed(3) + Nx(2,a)*dl(i,a) + Nx(1,a)*dl(j,a)
                  END IF
               END DO
            END IF

            SELECT CASE (outGrp)
            CASE (outGrp_J)
!           Jacobian := determinant of deformation gradient tensor
               resl(1) = detF
               sE(e)   = sE(e) + w*detF

            CASE (outGrp_F)
!           Deformation gradient tensor (F)
               IF (nsd .EQ. 3) THEN
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(1,3)
                  resl(4) = F(2,1)
                  resl(5) = F(2,2)
                  resl(6) = F(2,3)
                  resl(7) = F(3,1)
                  resl(8) = F(3,2)
                  resl(9) = F(3,3)
               ELSE
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(2,1)
                  resl(4) = F(2,2)
               END IF

            CASE (outGrp_strain)
!           Green-Lagrange strain tensor
               IF (cPhys .EQ. phys_lElas) THEN
                  resl(:) = ed(:)
               ELSE
                  C  = MATMUL(TRANSPOSE(F), F)
                  Eg = 0.5_RKIND * (C - Im)
                  ! resl is used to remap Eg
                  IF (nsd .EQ. 3) THEN
                     resl(1) = Eg(1,1)
                     resl(2) = Eg(2,2)
                     resl(3) = Eg(3,3)
                     resl(4) = Eg(1,2)
                     resl(5) = Eg(2,3)
                     resl(6) = Eg(3,1)
                  ELSE
                     resl(1) = Eg(1,1)
                     resl(2) = Eg(2,2)
                     resl(3) = Eg(1,2)
                  END IF
               END IF

            CASE (outGrp_stress, outGrp_cauchy, outGrp_mises)
               sigma = 0._RKIND
               IF (cPhys .EQ. phys_lElas) THEN
                  IF (nsd .EQ. 3) THEN
                     detF = lambda*(ed(1) + ed(2) + ed(3))
                     sigma(1,1) = detF + 2._RKIND*mu*ed(1)
                     sigma(2,2) = detF + 2._RKIND*mu*ed(2)
                     sigma(3,3) = detF + 2._RKIND*mu*ed(3)

                     sigma(1,2) = mu*ed(4)
                     sigma(2,3) = mu*ed(5)
                     sigma(3,1) = mu*ed(6)

                     sigma(2,1) = sigma(1,2)
                     sigma(3,2) = sigma(2,3)
                     sigma(1,3) = sigma(3,1)
                  ELSE
                     detF = lambda*(ed(1) + ed(2))
                     sigma(1,1) = detF + 2._RKIND*mu*ed(1)
                     sigma(2,2) = detF + 2._RKIND*mu*ed(2)
                     sigma(1,2) = mu*ed(3)
                     sigma(2,1) = sigma(1,2)
                  END IF

               ELSE IF (cPhys .EQ. phys_ustruct) THEN
                  p = 0._RKIND
                  DO a=1, fs%eNoN
                     p = p + N(a)*yl(k+1,a)
                  END DO
                  p = (-p)*detF

                  CALL GETPK2CCdev(eq(iEq)%dmn(cDmn), F, nFn, fN, tmX,
     2               ya, S, Dm, Ja)

                  C  = MATMUL(TRANSPOSE(F), F)
                  S  = S + p*MAT_INV(C, nsd)

                  P1 = MATMUL(F, S)
                  sigma = MATMUL(P1, TRANSPOSE(F))
                  IF (.NOT.ISZERO(detF)) sigma(:,:) = sigma(:,:) / detF

               ELSE IF (cPhys .EQ. phys_struct) THEN
                  CALL GETPK2CC(eq(iEq)%dmn(cDmn), F, nFn, fN, tmX, ya,
     2               S, Dm)
                  P1 = MATMUL(F, S)
                  sigma = MATMUL(P1, TRANSPOSE(F))
                  IF (.NOT.ISZERO(detF)) sigma(:,:) = sigma(:,:) / detF

               END IF

!              2nd Piola-Kirchhoff stress tensor
               IF (outGrp .EQ. outGrp_stress) THEN
                  IF (nsd .EQ. 3) THEN
                     resl(1) = S(1,1)
                     resl(2) = S(2,2)
                     resl(3) = S(3,3)
                     resl(4) = S(1,2)
                     resl(5) = S(2,3)
                     resl(6) = S(3,1)
                  ELSE
                     resl(1) = S(1,1)
                     resl(2) = S(2,2)
                     resl(3) = S(1,2)
                  END IF

!              Cauchy stress tensor
               ELSE IF (outGrp .EQ. outGrp_cauchy) THEN
                  IF (nsd .EQ. 3) THEN
                     resl(1) = sigma(1,1)
                     resl(2) = sigma(2,2)
                     resl(3) = sigma(3,3)
                     resl(4) = sigma(1,2)
                     resl(5) = sigma(2,3)
                     resl(6) = sigma(3,1)
                  ELSE
                     resl(1) = sigma(1,1)
                     resl(2) = sigma(2,2)
                     resl(3) = sigma(1,2)
                  END IF

!              Von Mises stress
               ELSE IF (outGrp .EQ. outGrp_mises) THEN
                  trS = MAT_TRACE(sigma, nsd) / REAL(nsd,KIND=RKIND)
                  DO l=1, nsd
                     sigma(l,l) = sigma(l,l) - trS
                  END DO
                  vmises  = SQRT(MAT_DDOT(sigma, sigma, nsd))
                  resl(1) = vmises
                  sE(e)   = sE(e) + w*vmises
               END IF

            CASE (outGrp_fS)
!           Fiber shortening (active strain model)
               resl(1) = ya

            END SELECT

            DO a=1, fs%eNoN
               Ac       = lM%IEN(a,e)
               sA(Ac)   = sA(Ac)   + w*N(a)
               sF(:,Ac) = sF(:,Ac) + w*N(a)*resl(:)
            END DO
         END DO
         IF (.NOT.ISZERO(Je)) sE(e) = sE(e)/Je
      END DO
      resE(:) = sE(:)

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT. iszero(sA(Ac))) THEN
           res(:,a) = res(:,a) + sF(:,Ac)/sA(Ac)
         ENDIF
      END DO

!     For higher order elements, values are interpolated at the edge
!     nodes and element centers using values computed at corners and
!     low-order shape functions
      IF (flag) THEN
         sF = 0._RKIND
         sA = 0._RKIND

         xi0 = 0._RKIND
         DO g=1, fs%nG
            xi0 = xi0 + fs%xi(:,g)
         END DO
         xi0 = xi0 / REAL(fs%nG, KIND=RKIND)

         DEALLOCATE(yl)
         ALLOCATE(yl(m,fs%eNoN), eNds(tnNo))
         eNds = 0
         DO e=1, lM%nEl
            cDmn  = DOMAIN(lM, iEq, e)
            cPhys = eq(iEq)%dmn(cDmn)%phys
            IF (cPhys .NE. phys_struct .AND.
     2          cPhys .NE. phys_ustruct .AND.
     3          cPhys .NE. phys_lElas) CYCLE

            yl = 0._RKIND
            DO a=1, fs%eNoN
               Ac      = lM%IEN(a,e)
               xl(:,a) = x(:,Ac)
               yl(:,a) = res(:,lM%lN(Ac))
            END DO

            Je = 0._RKIND
            DO g=1, fs%nG
               IF (g.EQ.1 .OR. .NOT.fs%lShpF) THEN
                  CALL GNN(fs%eNoN, insd, fs%Nx(:,:,g), xl, Nx, Jac, Im)
               END IF
               Je = Je + fs%w(g)*Jac
            END DO

            DO a=fs%eNoN+1, lM%eNoN
               Ac = lM%IEN(a,e)
               xp = x(:,Ac)

               xi = xi0
               CALL GETNNX(fs%eType, fs%eNoN, xl, fs%xib, fs%Nb, xp, xi,
     2            N, Nx)

!           Note that index `i' is reset here
               resl(:) = 0._RKIND
               DO i=1, fs%eNoN
                  resl(:) = resl(:) + N(i)*yl(:,i)
               END DO
               IF (eNds(Ac) .EQ. 0) eNds(Ac) = 1
               sF(:,Ac) = sF(:,Ac) + resl(:)*Je
               sA(Ac)   = sA(Ac)   + Je
            END DO
         END DO
         CALL COMMU(sF)
         CALL COMMU(sA)
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            IF (eNds(Ac).EQ.1 .AND. .NOT.ISZERO(sA(Ac))) THEN
               res(:,a) = res(:,a) + sF(:,Ac)/sA(Ac)
            END IF
         END DO
         DEALLOCATE(eNds)
      END IF
      DEALLOCATE (sA, sF, sE, xl, dl, yl, tmXl, ya_l, fN, resl, N, Nx)
      CALL DESTROY(fs)
      RETURN
      END SUBROUTINE TPOST
!####################################################################
!     Routine for post processing shell-based quantities
      SUBROUTINE SHLPOST(lM, m, res, resE, lD, iEq, outGrp)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: m, iEq, outGrp
      REAL(KIND=RKIND), INTENT(INOUT) :: res(m,lM%nNo), resE(lM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      LOGICAL incompFlag
      INTEGER(KIND=IKIND) a, b, e, g, i, j, k, l, iFn, Ac, nFn, insd,
     2   eNoN, cPhys, nwg
      REAL(KIND=RKIND) :: w, nu, Jac0, Jac, Je, lam3, detF, trS, vmises,
     2   nV0(3), nV(3), aCov0(3,2), aCov(3,2), aCnv0(3,2), aCnv(3,2),
     3   aa_0(2,2), aa_x(2,2), bb_0(2,2), bb_x(2,2), r0_xx(2,2,3),
     4   r_xx(2,2,3), Sm(3,2), Dm(3,3,3), Im(3,3), F(3,3), C(3,3),
     5   Eg(3,3), S(3,3), PK1(3,3), sigma(3,3)

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), sE(:), resl(:),
     2   dl(:,:), x0(:,:), xc(:,:), fN(:,:), fNa0(:,:), N(:), Nx(:,:),
     3   Nxx(:,:), Bb(:,:,:), tmpX(:,:)

      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     Set shell dimension := 2
      insd = nsd-1

!     Initialize tensor operations
      CALL TEN_INIT(insd)

!     Set eNoN (number of nodes per element)
      eNoN = lM%eNoN
      IF (lM%eType .EQ. eType_TRI3) eNoN = 2*eNoN

!     Allocate arrays
      ALLOCATE(sA(tnNo), sF(m,tnNo), sE(lM%nEl), resl(m), dl(tDof,eNoN),
     2   x0(3,eNoN), xc(3,eNoN), fN(3,nFn), fNa0(2,eNoN), ptr(eNoN),
     3   N(lM%eNoN), Nx(2,lM%eNoN), Nxx(3,lM%eNoN), Bb(3,3,6))

!     Initialize arrays
      sA  = 0._RKIND
      sF  = 0._RKIND
      sE  = 0._RKIND
      Bb  = 0._RKIND
      Nxx = 0._RKIND

!     Compute quantities at the element level and project them to nodes
      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_shell) CYCLE
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!        Get shell properties
         nu = eq(iEq)%dmn(cDmn)%prop(poisson_ratio)

!        Check for incompressibility
         incompFlag = .FALSE.
         IF (ISZERO(nu-0.5_RKIND)) incompFlag = .TRUE.

!        Get the reference configuration and displacement field
         x0 = 0._RKIND
         dl = 0._RKIND
         DO a=1, eNoN
            IF (a .LE. lM%eNoN) THEN
               Ac = lM%IEN(a,e)
               ptr(a) = Ac
            ELSE
               b  = a - lM%eNoN
               Ac = lM%eIEN(b,e)
               ptr(a) = Ac
               IF (Ac .EQ. 0) CYCLE
            END IF
            x0(:,a) = x(:,Ac)
            dl(:,a) = lD(:,Ac)
         END DO

!        Get the current configuration
         xc = 0._RKIND
         DO a=1, eNoN
            xc(1,a) = x0(1,a) + dl(i,a)
            xc(2,a) = x0(2,a) + dl(j,a)
            xc(3,a) = x0(3,a) + dl(k,a)
         END DO

!        Get fiber directions
         fN = 0._RKIND
         IF (ALLOCATED(lM%fN)) THEN
            DO iFn=1, nFn
               fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
            END DO
         END IF

!        Set number of integration points.
!        Note: Gauss integration performed for NURBS elements
!        Not required for constant-strain triangle elements
         IF (lM%eType .EQ. eType_TRI3) THEN
            nwg = 1
         ELSE
            nwg = lM%nG
         END IF

!        Update shapefunctions for NURBS elements
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         Je   = 0._RKIND
         resl = 0._RKIND
         DO g=1, nwg
            IF (lM%eType .NE. eType_TRI3) THEN
               ! Set element shape functions and their derivatives
               IF (lM%eType .EQ. eType_NRB) THEN
                  N   = lM%N(:,g)
                  Nx  = lM%Nx(:,:,g)
                  Nxx = lM%Nxx(:,:,g)
               ELSE
                  N   = lM%fs(1)%N(:,g)
                  Nx  = lM%fs(1)%Nx(:,:,g)
                  Nxx = lM%fs(1)%Nxx(:,:,g)
               END IF

!              Covariant and contravariant bases (ref. config.)
               CALL GNNS(eNoN, Nx, x0, nV0, aCov0, aCnv0)
               Jac0 = SQRT(NORM(nV0))
               nV0  = nV0/Jac0

!              Covariant and contravariant bases (spatial config.)
               CALL GNNS(eNoN, Nx, xc, nV, aCov, aCnv)
               Jac  = SQRT(NORM(nV))
               nV   = nV/Jac

!              Second derivatives for curvature coeffs. (ref. config)
               r0_xx(:,:,:) = 0._RKIND
               r_xx(:,:,:)  = 0._RKIND
               DO a=1, eNoN
                  r0_xx(1,1,:) = r0_xx(1,1,:) + Nxx(1,a)*x0(:,a)
                  r0_xx(2,2,:) = r0_xx(2,2,:) + Nxx(2,a)*x0(:,a)
                  r0_xx(1,2,:) = r0_xx(1,2,:) + Nxx(3,a)*x0(:,a)

                  r_xx(1,1,:) = r_xx(1,1,:) + Nxx(1,a)*xc(:,a)
                  r_xx(2,2,:) = r_xx(2,2,:) + Nxx(2,a)*xc(:,a)
                  r_xx(1,2,:) = r_xx(1,2,:) + Nxx(3,a)*xc(:,a)
               END DO
               r0_xx(2,1,:) = r0_xx(1,2,:)
               r_xx(2,1,:)  = r_xx(1,2,:)

!              Compute metric tensor (aa) and curvature coefficients(bb)
               aa_0 = 0._RKIND
               aa_x = 0._RKIND
               bb_0 = 0._RKIND
               bb_x = 0._RKIND
               DO l=1, nsd
                  aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)
                  aa_0(1,2) = aa_0(1,2) + aCov0(l,1)*aCov0(l,2)
                  aa_0(2,1) = aa_0(2,1) + aCov0(l,2)*aCov0(l,1)
                  aa_0(2,2) = aa_0(2,2) + aCov0(l,2)*aCov0(l,2)

                  aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)
                  aa_x(1,2) = aa_x(1,2) + aCov(l,1)*aCov(l,2)
                  aa_x(2,1) = aa_x(2,1) + aCov(l,2)*aCov(l,1)
                  aa_x(2,2) = aa_x(2,2) + aCov(l,2)*aCov(l,2)

                  bb_0(1,1) = bb_0(1,1) + r0_xx(1,1,l)*nV0(l)
                  bb_0(1,2) = bb_0(1,2) + r0_xx(1,2,l)*nV0(l)
                  bb_0(2,1) = bb_0(2,1) + r0_xx(2,1,l)*nV0(l)
                  bb_0(2,2) = bb_0(2,2) + r0_xx(2,2,l)*nV0(l)

                  bb_x(1,1) = bb_x(1,1) + r_xx(1,1,l)*nV(l)
                  bb_x(1,2) = bb_x(1,2) + r_xx(1,2,l)*nV(l)
                  bb_x(2,1) = bb_x(2,1) + r_xx(2,1,l)*nV(l)
                  bb_x(2,2) = bb_x(2,2) + r_xx(2,2,l)*nV(l)
               END DO

!              Set weight of the Gauss point
               w  = lM%w(g)*Jac0

            ELSE    ! for constant strain triangles
!              Set element shape functions and their derivatives
               N = lM%N(:,g)
               Nx(:,:) = lM%Nx(:,:,1)

!              Covariant and contravariant bases (ref. config.)
               ALLOCATE(tmpX(nsd,lM%eNoN))
               tmpX = x0(:,1:lM%eNoN)
               CALL GNNS(lM%eNoN, Nx, tmpX, nV0, aCov0, aCnv0)
               Jac0 = SQRT(NORM(nV0))
               nV0  = nV0/Jac0

!              Covariant and contravariant bases (spatial config.)
               tmpX = xc(:,1:lM%eNoN)
               CALL GNNS(lM%eNoN, Nx, tmpX, nV, aCov, aCnv)
               Jac = SQRT(NORM(nV))
               nV  = nV/Jac
               DEALLOCATE(tmpX)

!              Compute metric tensor (aa)
               aa_0 = 0._RKIND
               aa_x = 0._RKIND
               DO l=1, nsd
                  aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)
                  aa_0(1,2) = aa_0(1,2) + aCov0(l,1)*aCov0(l,2)
                  aa_0(2,1) = aa_0(2,1) + aCov0(l,2)*aCov0(l,1)
                  aa_0(2,2) = aa_0(2,2) + aCov0(l,2)*aCov0(l,2)

                  aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)
                  aa_x(1,2) = aa_x(1,2) + aCov(l,1)*aCov(l,2)
                  aa_x(2,1) = aa_x(2,1) + aCov(l,2)*aCov(l,1)
                  aa_x(2,2) = aa_x(2,2) + aCov(l,2)*aCov(l,2)
               END DO

               CALL SHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x, Bb,
     2            .FALSE.)

!              Set weight of the Gauss point
               w  = Jac0*0.5_RKIND
            END IF

!           Compute fiber direction in curvature coordinates
            fNa0 = 0._RKIND
            DO iFn=1, nFn
               DO l=1, nsd
                  fNa0(1,iFn) = fNa0(1,iFn) + fN(l,iFn)*aCnv0(l,1)
                  fNa0(2,iFn) = fNa0(2,iFn) + fN(l,iFn)*aCnv0(l,2)
               END DO
            END DO

!           Compute stress resultants and lambda3 (integrated through
!           the shell thickness)
            CALL SHL_STRS_RES(eq(cEq)%dmn(cDmn), nFn, fNa0, aa_0, aa_x,
     2         bb_0, bb_x, lam3, Sm, Dm)

!           Compute 3D deformation gradient tensor in shell continuum
            F = MAT_DYADPROD(aCov(:,1), aCnv0(:,1), 3)
     2        + MAT_DYADPROD(aCov(:,2), aCnv0(:,2), 3)
     3        + lam3*MAT_DYADPROD(nV, nV0, 3)

            detF = MAT_DET(F, nsd)

            Je = Je + w
            Im = MAT_ID(nsd)

            SELECT CASE (outGrp)
            CASE (outGrp_J)
!              Jacobian := determinant of deformation gradient tensor
               resl(1) = detF
               sE(e)   = sE(e) + w*detF

            CASE (outGrp_F)
!              Deformation gradient tensor (F)
               IF (nsd .EQ. 3) THEN
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(1,3)
                  resl(4) = F(2,1)
                  resl(5) = F(2,2)
                  resl(6) = F(2,3)
                  resl(7) = F(3,1)
                  resl(8) = F(3,2)
                  resl(9) = F(3,3)
               ELSE
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(2,1)
                  resl(4) = F(2,2)
               END IF

            CASE (outGrp_strain)
!              Cauchy-Green deformation (3D shell continuum)
               C  = MATMUL(TRANSPOSE(F), F)

!              Green-Lagrange strain tensor (3D shell continuum)
               Eg = 0.5_RKIND * (C - Im)

               ! resl is used to remap Eg
               IF (nsd .EQ. 3) THEN
                  resl(1) = Eg(1,1)
                  resl(2) = Eg(2,2)
                  resl(3) = Eg(3,3)
                  resl(4) = Eg(1,2)
                  resl(5) = Eg(2,3)
                  resl(6) = Eg(3,1)
               ELSE
                  resl(1) = Eg(1,1)
                  resl(2) = Eg(2,2)
                  resl(3) = Eg(1,2)
               END IF

            CASE (outGrp_stress, outGrp_cauchy, outGrp_mises)
               sigma  = 0._RKIND
               S(:,:) = 0._RKIND

!              2nd Piola-Kirchhoff stress (membrane + bending)
               S(1,1) = Sm(1,1) + Sm(1,2)
               S(2,2) = Sm(2,1) + Sm(2,2)
               S(1,2) = Sm(3,1) + Sm(3,2)
               S(2,1) = S(1,2)

!              Cauchy stress
               PK1 = MATMUL(F, S)
               sigma = MATMUL(PK1, TRANSPOSE(F))
               IF (.NOT.ISZERO(detF)) sigma(:,:) = sigma(:,:) / detF

!        2nd Piola-Kirchhoff stress tensor
               IF (outGrp .EQ. outGrp_stress) THEN
                  resl(1) = S(1,1)
                  resl(2) = S(2,2)
                  resl(3) = S(3,3)
                  resl(4) = S(1,2)
                  resl(5) = S(2,3)
                  resl(6) = S(3,1)

!        Cauchy stress tensor
               ELSE IF (outGrp .EQ. outGrp_cauchy) THEN
                  resl(1) = sigma(1,1)
                  resl(2) = sigma(2,2)
                  resl(3) = sigma(3,3)
                  resl(4) = sigma(1,2)
                  resl(5) = sigma(2,3)
                  resl(6) = sigma(3,1)

!        Von Mises stress
               ELSE IF (outGrp .EQ. outGrp_mises) THEN
                  trS = MAT_TRACE(sigma, nsd) / REAL(nsd,KIND=RKIND)
                  DO l=1, nsd
                     sigma(l,l) = sigma(l,l) - trS
                  END DO
                  vmises  = SQRT(MAT_DDOT(sigma, sigma, nsd))
                  resl(1) = vmises
                  sE(e)   = sE(e) + w*vmises
               END IF
            END SELECT

            DO a=1, lM%eNoN
               Ac       = lM%IEN(a,e)
               sA(Ac)   = sA(Ac)   + w*N(a)
               sF(:,Ac) = sF(:,Ac) + w*N(a)*resl(:)
            END DO
         END DO
         IF (.NOT.ISZERO(Je)) sE(e) = sE(e)/Je
      END DO
      resE(:) = sE(:)

!     Exchange data at the shared nodes across processes
      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT. iszero(sA(Ac))) THEN
           res(:,a) = res(:,a) + sF(:,Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE(sA, sF, sE, resl, dl, x0, xc, fN, fNa0, ptr, N, Nx)

      RETURN
      END SUBROUTINE SHLPOST
!####################################################################
!     Routine for post processing fiber directions
      SUBROUTINE FIBDIRPOST(lM, nFn, res, lD, iEq)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, nFn
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(INOUT) :: res(nFn*nsd,lM%nNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) a, b, e, g, Ac, eNoN, i, j, k, l, iFn, cPhys
      REAL(KIND=RKIND) w, Jac, F(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), dl(:,:), fN(:,:),
     2   fl(:,:), Nx(:,:), N(:), sA(:), sF(:,:)

      eNoN = lM%eNoN
      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1

      ALLOCATE (sA(tnNo), sF(nFn*nsd,tnNo), xl(nsd,eNoN), dl(tDof,eNoN),
     2   fN(nsd,lM%nFn), fl(nsd,lM%nFn), Nx(nsd,eNoN), N(eNoN))

      sA = 0._RKIND
      sF = 0._RKIND
      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         DO a=1, eNoN
            Ac       = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            dl(:,a)  = lD(:,Ac)
         END DO

         DO iFn=1, lM%nFn
            fN(:,iFn) = lM%fN((iFn-1)*nsd+1:iFn*nsd,e)
         END DO

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, F)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)

            F  = MAT_ID(nsd)
            DO a=1, eNoN
               IF (nsd .EQ. 3) THEN
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                  F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
                  F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
                  F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
                  F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)
               ELSE
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
               END IF
            END DO

            DO iFn=1, lM%nFn
               fl(:,iFn) = MATMUL(F, fN(:,iFn))
               fl(:,iFn) = fl(:,iFn)/SQRT(NORM(fl(:,iFn)))
            END DO

            DO a=1, eNoN
               Ac     = lM%IEN(a,e)
               sA(Ac) = sA(Ac) + w*N(a)
               DO iFn=1, lM%nFn
                  b = (iFn-1)*nsd
                  DO l=1, nsd
                     sF(b+l,Ac) = sF(b+l,Ac) + w*N(a)*fl(l,iFn)
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT.ISZERO(sA(Ac))) THEN
            res(:,a) = res(:,a) + sF(:,Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE (sA, sF, xl, dl, fN, fl, N, Nx)

      RETURN
      END SUBROUTINE FIBDIRPOST
!--------------------------------------------------------------------
!     Routine for post processing fiber alignment
      SUBROUTINE FIBALGNPOST(lM, res, lD, iEq)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(INOUT) :: res(1,lM%nNo)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, i, j, k, iFn, cPhys
      REAL(KIND=RKIND) w, Jac, sHat, F(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), dl(:,:), fN(:,:),
     2   fl(:,:), Nx(:,:), N(:), sA(:), sF(:)

      eNoN = lM%eNoN
      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1

      ALLOCATE (sA(tnNo), sF(tnNo), xl(nsd,eNoN), dl(tDof,eNoN),
     2   fN(nsd,2), fl(nsd,2), Nx(nsd,eNoN), N(eNoN))

      sA = 0._RKIND
      sF = 0._RKIND
      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_struct .AND.
     2       cPhys .NE. phys_ustruct) CYCLE
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         DO a=1, eNoN
            Ac       = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            dl(:,a)  = lD(:,Ac)
         END DO

         fN(:,1) = lM%fN(1:nsd,e)
         fN(:,2) = lM%fN(nsd+1:2*nsd,e)

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, F)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)

            F  = MAT_ID(nsd)
            DO a=1, eNoN
               IF (nsd .EQ. 3) THEN
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                  F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
                  F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
                  F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
                  F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)
               ELSE
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
               END IF
            END DO

            DO iFn=1, 2
               fl(:,iFn) = MATMUL(F, fN(:,iFn))
               fl(:,iFn) = fl(:,iFn)/SQRT(NORM(fl(:,iFn)))
            END DO
            sHat = NORM(fl(:,1), fl(:,2))

            DO a=1, eNoN
               Ac     = lM%IEN(a,e)
               sA(Ac) = sA(Ac) + w*N(a)
               sF(Ac) = sF(Ac) + w*N(a)*sHat
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT.ISZERO(sA(Ac))) THEN
           res(1,a) = res(1,a) + sF(Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE (sA, sF, xl, dl, fN, fl, N, Nx)

      RETURN
      END SUBROUTINE FIBALGNPOST
!####################################################################
!     Compute fiber stretch based on 4th invariant: I_{4,f}
      SUBROUTINE FIBSTRETCH(iEq, lM, lD, res)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: res(lM%nNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, i, j, k, cPhys
      REAL(KIND=RKIND) w, Jac, I4f, fl(nsd), F(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), dl(:,:), Nx(:,:), N(:),
     2   sA(:), sF(:)

      eNoN = lM%eNoN
      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1

      ALLOCATE (sA(tnNo), sF(tnNo), xl(nsd,eNoN), dl(tDof,eNoN),
     2   Nx(nsd,eNoN), N(eNoN))

      sA = 0._RKIND
      sF = 0._RKIND
      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         DO a=1, eNoN
            Ac       = lM%IEN(a,e)
            xl(:,a)  = x(:,Ac)
            dl(:,a)  = lD(:,Ac)
         END DO

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, F)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)

            F  = MAT_ID(nsd)
            DO a=1, eNoN
               IF (nsd .EQ. 3) THEN
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                  F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
                  F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
                  F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
                  F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)
               ELSE
                  F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                  F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                  F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                  F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
               END IF
            END DO

            fl  = MATMUL(F, lM%fN(1:nsd,e))
            I4f = NORM(fl)
            DO a=1, eNoN
               Ac     = lM%IEN(a,e)
               sA(Ac) = sA(Ac) + w*N(a)
               sF(Ac) = sF(Ac) + w*N(a)*I4f
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      res = 0._RKIND
      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT.ISZERO(sA(Ac))) THEN
            res(a) = res(a) + sF(Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE (sA, sF, xl, dl, N, Nx)

      RETURN
      END SUBROUTINE FIBSTRETCH
!####################################################################
!     Routine for post processing divergence of velocity
      SUBROUTINE DIVPOST(lM, res, lY, lD, iEq)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq
      REAL(KIND=RKIND), INTENT(INOUT) :: res(1,lM%nNo)
      REAL(KIND=RKIND), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, g, Ac, eNoN, i, j, k, cPhys
      REAL(KIND=RKIND) w, Jac, divV, vx(nsd,nsd), ksix(nsd,nsd),
     2   F(nsd,nsd), Fi(nsd,nsd), VxFi(nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), yl(:,:), dl(:,:), N(:),
     2   Nx(:,:), sA(:), sF(:)

      eNoN = lM%eNoN
      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1

      ALLOCATE (sA(tnNo), sF(tnNo), xl(nsd,eNoN), yl(tDof,eNoN),
     2   dl(tDof,eNoN), N(eNoN), Nx(nsd,eNoN))

      sA   = 0._RKIND
      sF   = 0._RKIND
      DO e=1, lM%nEl
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys

         DO a=1, eNoN
            Ac      = lM%IEN(a,e)
            xl(:,a) = x(:,Ac)
            yl(:,a) = lY(:,Ac)
            dl(:,a) = lD(:,Ac)
         END DO

         DO g=1, lM%nG
            IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
               CALL GNN(eNoN, nsd, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)

            IF ((cPhys.EQ.phys_fluid) .OR. (cPhys.EQ.phys_CMM)) THEN
               vx = 0._RKIND
               IF (nsd .EQ. 3) THEN
                  DO a=1, eNoN
                     vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a)
                     vx(2,2) = vx(2,2) + Nx(2,a)*yl(2,a)
                     vx(3,3) = vx(3,3) + Nx(3,a)*yl(3,a)
                  END DO
                  divV = vx(1,1) + vx(2,2) + vx(3,3)
               ELSE
                  DO a=1, eNoN
                     vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a)
                     vx(2,2) = vx(2,2) + Nx(2,a)*yl(2,a)
                  END DO
                  divV = vx(1,1) + vx(2,2)
               END IF

            ELSE IF (cPhys .EQ. phys_ustruct) THEN
               vx = 0._RKIND
               F  = MAT_ID(nsd)
               IF (nsd .EQ. 3) THEN
                  DO a=1, eNoN
                     vx(1,1) = vx(1,1) + Nx(1,a)*yl(i,a)
                     vx(1,2) = vx(1,2) + Nx(2,a)*yl(i,a)
                     vx(1,3) = vx(1,3) + Nx(3,a)*yl(i,a)
                     vx(2,1) = vx(2,1) + Nx(1,a)*yl(j,a)
                     vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)
                     vx(2,3) = vx(2,3) + Nx(3,a)*yl(j,a)
                     vx(3,1) = vx(3,1) + Nx(1,a)*yl(k,a)
                     vx(3,2) = vx(3,2) + Nx(2,a)*yl(k,a)
                     vx(3,3) = vx(3,3) + Nx(3,a)*yl(k,a)

                     F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                     F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                     F(1,3) = F(1,3) + Nx(3,a)*dl(i,a)
                     F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                     F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                     F(2,3) = F(2,3) + Nx(3,a)*dl(j,a)
                     F(3,1) = F(3,1) + Nx(1,a)*dl(k,a)
                     F(3,2) = F(3,2) + Nx(2,a)*dl(k,a)
                     F(3,3) = F(3,3) + Nx(3,a)*dl(k,a)
                  END DO
                  Fi = MAT_INV(F,3)

                  VxFi(1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1) +
     2               vx(1,3)*Fi(3,1)
                  VxFi(2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2) +
     2               vx(2,3)*Fi(3,2)
                  VxFi(3) = vx(3,1)*Fi(1,3) + vx(3,2)*Fi(2,3) +
     2               vx(3,3)*Fi(3,3)

                  divV = VxFi(1) + VxFi(2) + VxFi(3)
               ELSE
                  DO a=1, eNoN
                     vx(1,1) = vx(1,1) + Nx(1,a)*yl(i,a)
                     vx(1,2) = vx(1,2) + Nx(2,a)*yl(i,a)
                     vx(2,1) = vx(2,1) + Nx(1,a)*yl(j,a)
                     vx(2,2) = vx(2,2) + Nx(2,a)*yl(j,a)

                     F(1,1) = F(1,1) + Nx(1,a)*dl(i,a)
                     F(1,2) = F(1,2) + Nx(2,a)*dl(i,a)
                     F(2,1) = F(2,1) + Nx(1,a)*dl(j,a)
                     F(2,2) = F(2,2) + Nx(2,a)*dl(j,a)
                  END DO
                  Fi = MAT_INV(F,2)

                  VxFi(1) = vx(1,1)*Fi(1,1) + vx(1,2)*Fi(2,1)
                  VxFi(2) = vx(2,1)*Fi(1,2) + vx(2,2)*Fi(2,2)

                  divV = VxFi(1) + VxFi(2)
               END IF
            END IF

            DO a=1, eNoN
               Ac     = lM%IEN(a,e)
               sA(Ac) = sA(Ac) + w*N(a)
               sF(Ac) = sF(Ac) + w*N(a)*divV
            END DO
         END DO
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT. iszero(sA(Ac))) THEN
           res(1,a) = res(1,a) + sF(Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE (sA, sF, xl, yl, dl, N, Nx)

      RETURN
      END SUBROUTINE DIVPOST
!####################################################################
!     Postprocessing function - used to convert restart bin files
!     into VTK format.
      SUBROUTINE PPBIN2VTK()
      USE COMMOD
      IMPLICIT NONE

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: iTS, ierr
      REAL(KIND=RKIND) :: rtmp(3)
      CHARACTER(LEN=stdL) :: stmp, fName, sepLine

      sepLine = REPEAT("-", 50)
      std  = TRIM(sepLine)
      std  = CLR(" Post-process conversion from bin to vtk files",3)

      DO iTS=1, nTS
         IF (MOD(iTS,stFileIncr) .EQ. 0) THEN
            WRITE(stmp,'(I3.3)') iTS
            IF (iTS .GE. 1000) stmp = STR(iTS)

!           Ignore if vtu file already exists
            IF (cm%mas()) THEN
               fName = TRIM(saveName)//"_"//TRIM(ADJUSTL(stmp))//".vtu"
               INQUIRE(FILE=TRIM(fName), EXIST=flag)
            END IF
            CALL cm%bcast(flag)
            IF (flag) CYCLE

!           Ignore if bin file does not exist
            fName = TRIM(stFileName)//"_"//TRIM(ADJUSTL(stmp))//".bin"
            INQUIRE(FILE=TRIM(fName), EXIST=flag)
            IF (.NOT.flag) CYCLE

            std = TRIM(sepLine)
            std = " "//CLR("<< BIN ---> VTK >>")//"    "//TRIM(fName)
            CALL INITFROMBIN(fName, rtmp)
            CALL WRITEVTUS(Ao, Yo, Do, .FALSE.)
         END IF
      END DO
      std = TRIM(sepLine)

      CALL FINALIZE
      CALL MPI_FINALIZE(ierr)
      STOP

      END SUBROUTINE PPBIN2VTK
!####################################################################
!     Routine for post processing shell stress tensor
      SUBROUTINE SHLTPOST(lM, m, res, resE, lD, lY, iEq, outGrp)
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: m, iEq, outGrp
      REAL(KIND=RKIND), INTENT(INOUT) :: res(m,lM%nNo), resE(lM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(tDof,tnNo), lY(tDof,tnNo)

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, g, Ac, i, j, k, l, cPhys, insd, nFn,
     2   gs, iFn, li, lj, nwg, eNoN 
      REAL(KIND=RKIND) w, Jac, Jac0, detF, Je, tmX, ya, nu, mu, p, trS,
     2   vmises, xi(nsd), xi0(nsd), xp(nsd), ed(nsymd), Im(nsd,nsd),
     3   F(nsd,nsd), C(nsd,nsd), Eg(nsd,nsd), P1(nsd,nsd), S(nsd,nsd),
     4   sigma(nsd,nsd), aa_0(2,2), aa_x(2,2), aCov0(3,2), aCov(3,2),
     5   aCnv0(3,2), aCnv(3,2), nV0(3), nV(3), r0_xx(2,2,3), 
     6   r_xx(2,2,3), Sml(3), Dml(3,3), ht, wh(3), xis, gg_0(2,2),
     7   gg_x(2,2), bb_0(2,2), bb_x(2,2), xXi(3,2), Gmat(2,2), g33
      TYPE(fsType) :: fs

      INTEGER, ALLOCATABLE :: eNds(:), ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), dl(:,:), yl(:,:),
     2   tmXl(:), ya_l(:), fN(:,:), resl(:), Nx(:,:), N(:), sA(:),
     3   sF(:,:), sE(:), x0(:,:), xc(:,:), fNa0(:,:), Nxx(:,:),
     4   tmpX(:,:)

      dof  = eq(iEq)%dof
      i    = eq(iEq)%s
      j    = i + 1
      k    = j + 1
      nFn  = lM%nFn
      IF (nFn .EQ. 0) nFn = 1

!     For higher order elements, we lower the order of shape functions
!     to compute the quantities at corner nodes of elements. We then
!     use these lower order shape functions to interpolate the values
!     at edge nodes and elements centers (if applicable)
      flag = .FALSE.
      IF (lM%eType.EQ.eType_TRI6  .OR. lM%eType.EQ.eType_QUD8  .OR.
     2    lM%eType.EQ.eType_QUD9  .OR. lM%eType.EQ.eType_TET10 .OR.
     3    lM%eType.EQ.eType_HEX20 .OR. lM%eType .EQ. eType_HEX27) THEN
         flag =.TRUE.
         CALL SETTHOODFS(fs, lM%eType)
      ELSE
         fs%eType = lM%eType
         fs%lShpF = lM%lShpF
         fs%eNoN  = lM%eNoN
         fs%nG    = lM%nG
      END IF
      CALL INITFS(fs, 2)

      eNoN = fs%eNoN
      IF (fs%eType .EQ. eType_TRI3) eNoN = 2*eNoN

      ALLOCATE(sA(tnNo), sF(m,tnNo), sE(lM%nEl), xl(nsd,eNoN),
     2   dl(tDof,eNoN), yl(tDof,eNoN), fN(nsd,nFn), tmXl(eNoN),
     3   ya_l(eNoN), resl(m), Nx(2,fs%eNoN), N(fs%eNoN),
     4   ptr(eNoN), x0(3,eNoN),xc(3,eNoN),fNa0(2,nFn),
     5                   Nxx(3,fs%eNoN))

      sA   = 0._RKIND
      sF   = 0._RKIND
      sE   = 0._RKIND

      insd = nsd
      IF (lM%lShl) insd = 2

!     Initialize tensor operations
      CALL TEN_INIT(insd)

      IF (lM%lFib) insd = 1

      DO e=1, lM%nEl
         cDmn  = DOMAIN(lM, iEq, e)
         cPhys = eq(iEq)%dmn(cDmn)%phys
         IF (cPhys .NE. phys_shell) CYCLE
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

         fN = 0._RKIND
         IF (ALLOCATED(lM%fN)) THEN
            DO l=1, nFn
               fN(:,l) = lM%fN((l-1)*nsd+1:l*nsd,e)
            END DO
         END IF
         dl   = 0._RKIND
         yl   = 0._RKIND
         tmXl = 0._RKIND
         ya_l = 0._RKIND
         DO a=1, eNoN
            IF (a .LE. lM%eNoN) THEN
               Ac = lM%IEN(a,e)
               ptr(a) = Ac
            ELSE
               b  = a - lM%eNoN
               Ac = lM%eIEN(b,e)
               ptr(a) = Ac
               IF (Ac .EQ. 0) CYCLE
            END IF
            b  = lM%lN(Ac)
            xl(:,a) = x(:,Ac)
            dl(:,a) = lD(:,Ac)
            yl(:,a) = lY(:,Ac)
            IF (ecCpld) THEN
               IF (ALLOCATED(lM%tmX)) THEN
                  tmXl(a) = lM%tmX(b)
               END IF
               IF (ALLOCATED(ec_Ya)) THEN
                  ya_l(a) = ec_Ya(Ac)
               ELSE
                  ya_l(a) = eq(cEq)%dmn(cDmn)%ec%Ya
               END IF
            END IF
         END DO

         x0(:,:) = xl(:,:)

         DO a=1, eNoN
            xc(1,a) = x0(1,a) + dl(i,a)
            xc(2,a) = x0(2,a) + dl(j,a)
            xc(3,a) = x0(3,a) + dl(k,a)
         END DO

!        Gauss integration
         Je = 0._RKIND

         IF (fs%eType .NE. eType_TRI3) THEN
            nwg =  fs%nG
         ELSE
            nwg =  1
         END IF
         
         DO g=1, nwg
            IF (fs%eType .NE. eType_TRI3) THEN
            ! Post processing SHELL3D
               IF (lM%eType .EQ. eType_NRB) THEN
                  N   = lM%N(:,g)
                  Nx  = lM%Nx(:,:,g)
                  Nxx = lM%Nxx(:,:,g)
               ELSE
                  N   = lM%fs(1)%N(:,g)
                  Nx  = lM%fs(1)%Nx(:,:,g)
                  Nxx = lM%fs(1)%Nxx(:,:,g)
               END IF

               CALL GNNS(fs%eNoN, Nx, x0, nV0, aCov0, aCnv0)
               Jac0 = SQRT(NORM(nV0))
               nV0  = nV0/Jac0
               
               CALL GNNS(fs%eNoN, Nx, xc, nV, aCov, aCnv)
               Jac = SQRT(NORM(nV))
               nV  = nV/Jac
               r0_xx(:,:,:) = 0._RKIND
               r_xx(:,:,:) = 0._RKIND
               DO a=1, fs%eNoN
                  r0_xx(1,1,:) = r0_xx(1,1,:) + Nxx(1,a)*x0(:,a)
                  r0_xx(2,2,:) = r0_xx(2,2,:) + Nxx(2,a)*x0(:,a)
                  r0_xx(1,2,:) = r0_xx(1,2,:) + Nxx(3,a)*x0(:,a)
                        
                  r_xx(1,1,:) = r_xx(1,1,:) + Nxx(1,a)*xc(:,a)
                  r_xx(2,2,:) = r_xx(2,2,:) + Nxx(2,a)*xc(:,a)
                  r_xx(1,2,:) = r_xx(1,2,:) + Nxx(3,a)*xc(:,a)
               END DO
               r0_xx(2,1,:) = r0_xx(1,2,:)
               r_xx(2,1,:) = r_xx(1,2,:)
               aa_0 = 0._RKIND
               aa_x = 0._RKIND
               bb_0 = 0._RKIND
               bb_x = 0._RKIND
               DO l=1, nsd
                  aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)
                  aa_0(1,2) = aa_0(1,2) + aCov0(l,1)*aCov0(l,2)
                  aa_0(2,1) = aa_0(2,1) + aCov0(l,2)*aCov0(l,1)
                  aa_0(2,2) = aa_0(2,2) + aCov0(l,2)*aCov0(l,2)

                  aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)
                  aa_x(1,2) = aa_x(1,2) + aCov(l,1)*aCov(l,2)
                  aa_x(2,1) = aa_x(2,1) + aCov(l,2)*aCov(l,1)
                  aa_x(2,2) = aa_x(2,2) + aCov(l,2)*aCov(l,2)
                  
                  bb_0(1,1) = bb_0(1,1) + r0_xx(1,1,l)*nV0(l)
                  bb_0(1,2) = bb_0(1,2) + r0_xx(1,2,l)*nV0(l)
                  bb_0(2,1) = bb_0(2,1) + r0_xx(2,1,l)*nV0(l)
                  bb_0(2,2) = bb_0(2,2) + r0_xx(2,2,l)*nV0(l)

                  bb_x(1,1) = bb_x(1,1) + r_xx(1,1,l)*nV(l)
                  bb_x(1,2) = bb_x(1,2) + r_xx(1,2,l)*nV(l)
                  bb_x(2,1) = bb_x(2,1) + r_xx(2,1,l)*nV(l)
                  bb_x(2,2) = bb_x(2,2) + r_xx(2,2,l)*nV(l)
               END DO
            ! END Post processing SHELL3D
            ELSE
            ! Post processing SHELLCST
               N = lM%N(:,g)
               Nx(:,:) = lM%Nx(:,:,1)
               ALLOCATE(tmpX(nsd,fs%eNoN))
               tmpX = x0(:,1:fs%eNoN)
               CALL GNNS(fs%eNoN, Nx, tmpX, nV0, aCov0, aCnv0)
               Jac0 = SQRT(NORM(nV0))
               nV0  = nV0/Jac0

               tmpX = xc(:,1:lM%eNoN)
               CALL GNNS(fs%eNoN, Nx, tmpX, nV, aCov, aCnv)
               Jac = SQRT(NORM(nV))
               nV  = nV/Jac
               DEALLOCATE(tmpX)

               aa_0 = 0._RKIND
               aa_x = 0._RKIND
               DO l=1, nsd
                  aa_0(1,1) = aa_0(1,1) + aCov0(l,1)*aCov0(l,1)
                  aa_0(1,2) = aa_0(1,2) + aCov0(l,1)*aCov0(l,2)
                  aa_0(2,1) = aa_0(2,1) + aCov0(l,2)*aCov0(l,1)
                  aa_0(2,2) = aa_0(2,2) + aCov0(l,2)*aCov0(l,2)

                  aa_x(1,1) = aa_x(1,1) + aCov(l,1)*aCov(l,1)
                  aa_x(1,2) = aa_x(1,2) + aCov(l,1)*aCov(l,2)
                  aa_x(2,1) = aa_x(2,1) + aCov(l,2)*aCov(l,1)
                  aa_x(2,2) = aa_x(2,2) + aCov(l,2)*aCov(l,2)
               END DO

               CALL PSHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x)
            ! END Post processing SHELLCST
            END IF

            fNa0 = 0._RKIND
            DO iFn=1, nFn
               DO l=1, 3
                  fNa0(1,iFn) = fNa0(1,iFn) + fN(l,iFn)*aCnv0(l,1)
                  fNa0(2,iFn) = fNa0(2,iFn) + fN(l,iFn)*aCnv0(l,2)
               END DO
            END DO

            ht = eq(iEq)%dmn(cDmn)%prop(shell_thickness)
            nu = eq(iEq)%dmn(cDmn)%prop(poisson_ratio)
            gg_0(:,:) = aa_0(:,:)
            gg_x(:,:) = aa_x(:,:)
            IF (ISZERO(nu-0.5_RKIND)) THEN
               CALL GETPK2CC_SHLi(eq(iEq)%dmn(cDmn), nFn, fNa0, 
     2               gg_0, gg_x, g33, Sml, Dml)
            ELSE
               CALL GETPK2CC_SHLc(eq(iEq)%dmn(cDmn), nFn, fNa0, 
     2               gg_0, gg_x, g33, Sml, Dml)
            END IF

            F = 0._RKIND
            F = F + MAT_DYADPROD(aCov(:,1),aCnv0(:,1),3)
            F = F + MAT_DYADPROD(aCov(:,2),aCnv0(:,2),3)
            F = F + SQRT(g33) * MAT_DYADPROD(nV0,nV,3)
            detF = MAT_DET(F, nsd)
            w  = fs%w(g)*Jac0
            Je = Je + w
            Im = MAT_ID(nsd)

            ed = 0._RKIND
            
            SELECT CASE (outGrp)
            CASE (outGrp_J)
!           Jacobian := determinant of deformation gradient tensor
               resl(1) = detF
               sE(e)   = sE(e) + w*detF

            CASE (outGrp_F)
!           Deformation gradient tensor (F)
               IF (nsd .EQ. 3) THEN
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(1,3)
                  resl(4) = F(2,1)
                  resl(5) = F(2,2)
                  resl(6) = F(2,3)
                  resl(7) = F(3,1)
                  resl(8) = F(3,2)
                  resl(9) = F(3,3)
               ELSE
                  resl(1) = F(1,1)
                  resl(2) = F(1,2)
                  resl(3) = F(2,1)
                  resl(4) = F(2,2)
               END IF

            CASE (outGrp_strain)
!           Green-Lagrange strain tensor
               C  = MATMUL(TRANSPOSE(F), F)
               Eg = 0.5_RKIND * (C - Im)
               ! resl is used to remap Eg
               IF (nsd .EQ. 3) THEN
                  resl(1) = Eg(1,1)
                  resl(2) = Eg(2,2)
                  resl(3) = Eg(3,3)
                  resl(4) = Eg(1,2)
                  resl(5) = Eg(2,3)
                  resl(6) = Eg(3,1)
               ELSE
                  resl(1) = Eg(1,1)
                  resl(2) = Eg(2,2)
                  resl(3) = Eg(1,2)
               END IF

            CASE (outGrp_stress, outGrp_cauchy, outGrp_mises)
               sigma = 0._RKIND

               S(1,1) = Sml(1)
               S(1,2) = Sml(3)
               S(2,1) = S(1,2)
               S(2,2) = Sml(2)
               S(1,3) = 0._RKIND
               S(3,1) = S(1,3)
               S(2,3) = 0._RKIND
               S(3,2) = S(2,3)
               S(3,3) = 0._RKIND
               P1 = MATMUL(F, S)
               sigma = MATMUL(P1, TRANSPOSE(F))
               IF (.NOT.ISZERO(detF)) sigma(:,:) = sigma(:,:) / detF

!        2nd Piola-Kirchhoff stress tensor
               IF (outGrp .EQ. outGrp_stress) THEN
                  IF (nsd .EQ. 3) THEN
                     resl(1) = S(1,1)
                     resl(2) = S(2,2)
                     resl(3) = S(3,3)
                     resl(4) = S(1,2)
                     resl(5) = S(2,3)
                     resl(6) = S(3,1)
                  ELSE
                     resl(1) = S(1,1)
                     resl(2) = S(2,2)
                     resl(3) = S(1,2)
                  END IF

!        Cauchy stress tensor
               ELSE IF (outGrp .EQ. outGrp_cauchy) THEN
                  IF (nsd .EQ. 3) THEN
                     resl(1) = sigma(1,1)
                     resl(2) = sigma(2,2)
                     resl(3) = sigma(3,3)
                     resl(4) = sigma(1,2)
                     resl(5) = sigma(2,3)
                     resl(6) = sigma(3,1)
                  ELSE
                     resl(1) = sigma(1,1)
                     resl(2) = sigma(2,2)
                     resl(3) = sigma(1,2)
                  END IF

!        Von Mises stress
               ELSE IF (outGrp .EQ. outGrp_mises) THEN
                  trS = MAT_TRACE(sigma, nsd) / REAL(nsd,KIND=RKIND)
                  DO l=1, nsd
                     sigma(l,l) = sigma(l,l) - trS
                  END DO
                  vmises  = SQRT(MAT_DDOT(sigma, sigma, nsd))
                  resl(1) = vmises
                  sE(e)   = sE(e) + w*vmises
               END IF

            CASE (outGrp_fS)
!           Fiber shortening (active strain model)
               resl(1) = ya

            END SELECT

            DO a=1, fs%eNoN
               Ac       = lM%IEN(a,e)
               sA(Ac)   = sA(Ac)   + w*N(a)
               sF(:,Ac) = sF(:,Ac) + w*N(a)*resl(:)
            END DO
         END DO
         IF (.NOT.ISZERO(Je)) sE(e) = sE(e)/Je
      END DO
      resE(:) = sE(:)

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO a=1, lM%nNo
         Ac = lM%gN(a)
         IF (.NOT. iszero(sA(Ac))) THEN
           res(:,a) = res(:,a) + sF(:,Ac)/sA(Ac)
         ENDIF
      END DO

      DEALLOCATE (sA, sF, sE, xl, dl, yl, tmXl, ya_l, fN, resl, N, Nx,
     2                  x0, xc, fNa0, Nxx)
      CALL DESTROY(fs)

      RETURN
      END SUBROUTINE SHLTPOST
!####################################################################
!     This subroutine computes bending strain bb_0 and bb_x for CST
      SUBROUTINE PSHELLBENDCST(lM, e, ptr, x0, xc, bb_0, bb_x)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e, ptr(6)
      REAL(KIND=RKIND), INTENT(INOUT) ::  x0(3,6), xc(3,6)
      REAL(KIND=RKIND), INTENT(OUT) :: bb_0(2,2), bb_x(2,2)

      LOGICAL :: bFlag, lFix(3)
      INTEGER(KIND=IKIND) :: i, j, p, f, eNoN

      REAL(KIND=RKIND) :: Jac0, Jac, cI, aIi, nV0(3), nV(3), eI(3),
     2   nI(3), eIeI(3,3), nInI(3,3), eIaP(3,3), Im(3,3), aCov(3,2),
     3   aCnv(3,2), aCov0(3,2), aCnv0(3,2)

      REAL(KIND=RKIND) :: a0(3,6), a(3,6), adg0(3,3), adg(3,3),
     2   xi0(3,3), xi(3,3), Tm0(3,3), Tm(3,3), v0(3), v(3), B1b(3,6),
     3   Nm(3,3), Mm(3,3,2), H1(6,18), H2(18,18), H3(3,18), H1H2(6,18),
     4   H3H2(3,18), Bb1(3,18), tmpA(3,3)

      eNoN = 2*lM%eNoN
!     Boundary element check
      bFlag = .FALSE.
      DO j=lM%eNoN+1, eNoN
         IF (ptr(j) .EQ. 0) THEN
            bFlag = .TRUE.
            EXIT
         END IF
      END DO

!     Edge vectors of the main element (reference config)
      a0(:,1) = x0(:,3) - x0(:,2)
      a0(:,2) = x0(:,1) - x0(:,3)
      a0(:,3) = x0(:,2) - x0(:,1)

!     Edge vectors of the main element (current config)
      a(:,1)  = xc(:,3) - xc(:,2)
      a(:,2)  = xc(:,1) - xc(:,3)
      a(:,3)  = xc(:,2) - xc(:,1)

!     Covariant and contravariant bases in reference config
      tmpA = x0(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV0, aCov0, aCnv0)
      Jac0 = SQRT(NORM(nV0))
      nV0  = nV0/Jac0

!     Covariant and contravariant bases in current config
      tmpA = xc(:,1:lM%eNoN)
      CALL GNNS(lM%eNoN, lM%Nx(:,:,1), tmpA, nV, aCov, aCnv)
      Jac = SQRT(NORM(nV))
      nV  = nV/Jac

!     Update the position vector of the `artificial' or `ghost' nodes
!     depending on the boundary condition.
      IF (bFlag) THEN
         DO j=lM%eNoN+1, eNoN
            IF (ptr(j) .NE. 0) CYCLE
            i = j - lM%eNoN
            p = i - 1
            IF (i .EQ. 1) p = 3
!           Reference config
!           eI = eI0 = aI0/|aI0| (reference config)
            aIi   = 1._RKIND/SQRT(NORM(a0(:,i)))
            eI(:) = a0(:,i) * aIi
!           nI = nI0 = eI0 x n0 (reference config)
            nI(1) = eI(2)*nV0(3) - eI(3)*nV0(2)
            nI(2) = eI(3)*nV0(1) - eI(1)*nV0(3)
            nI(3) = eI(1)*nV0(2) - eI(2)*nV0(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            x0(1,j) = 2._RKIND*(nInI(1,1)*a0(1,p) + nInI(1,2)*a0(2,p) +
     2         nInI(1,3)*a0(3,p)) + x0(1,i)
            x0(2,j) = 2._RKIND*(nInI(2,1)*a0(1,p) + nInI(2,2)*a0(2,p) +
     2         nInI(2,3)*a0(3,p)) + x0(2,i)
            x0(3,j) = 2._RKIND*(nInI(3,1)*a0(1,p) + nInI(3,2)*a0(2,p) +
     2         nInI(3,3)*a0(3,p)) + x0(3,i)

!           Current config
!           eI = aI/|aI| (current config)
            aIi   = 1._RKIND/SQRT(NORM(a(:,i)))
            eI(:) = a(:,i)*aIi
!           nI = eI x n (currnt config)
            nI(1) = eI(2)*nV(3) - eI(3)*nV(2)
            nI(2) = eI(3)*nV(1) - eI(1)*nV(3)
            nI(3) = eI(1)*nV(2) - eI(2)*nV(1)
!           xJ = xI + 2(nI \ctimes nI)aP
            nInI    = MAT_DYADPROD(nI, nI, 3)
            xc(1,j) = 2._RKIND*(nInI(1,1)*a(1,p) + nInI(1,2)*a(2,p) +
     2         nInI(1,3)*a(3,p)) + xc(1,i)
            xc(2,j) = 2._RKIND*(nInI(2,1)*a(1,p) + nInI(2,2)*a(2,p) +
     2         nInI(2,3)*a(3,p)) + xc(2,i)
            xc(3,j) = 2._RKIND*(nInI(3,1)*a(1,p) + nInI(3,2)*a(2,p) +
     2         nInI(3,3)*a(3,p)) + xc(3,i)

            IF (BTEST(lM%sbc(i,e),bType_fix)) xc(:,j) = x0(:,j)
         END DO
      END IF

!     Edge vector of surrounding nodes (reference config)
      a0(:,4)   = x0(:,4) - x0(:,2)
      a0(:,5)   = x0(:,5) - x0(:,3)
      a0(:,6)   = x0(:,6) - x0(:,1)

!     Edge vector of surrounding nodes (current config)
      a(:,4)    = xc(:,4) - xc(:,2)
      a(:,5)    = xc(:,5) - xc(:,3)
      a(:,6)    = xc(:,6) - xc(:,1)

!     a.gCnv (reference config)
      adg0(1,1) = NORM(a0(:,4), aCnv0(:,1)) ! xi_4
      adg0(1,2) = NORM(a0(:,5), aCnv0(:,1)) ! xi_5
      adg0(1,3) = NORM(a0(:,6), aCnv0(:,1)) ! xi_6

      adg0(2,1) = NORM(a0(:,4), aCnv0(:,2)) ! eta_4
      adg0(2,2) = NORM(a0(:,5), aCnv0(:,2)) ! eta_5
      adg0(2,3) = NORM(a0(:,6), aCnv0(:,2)) ! eta_6

      adg0(3,1) = NORM(a0(:,4), nV0)        ! z_4
      adg0(3,2) = NORM(a0(:,5), nV0)        ! z_5
      adg0(3,3) = NORM(a0(:,6), nV0)        ! z_6

!     a.gCnv (current config)
      adg(1,1)  = NORM(a(:,4), aCnv(:,1)) ! xi_4
      adg(1,2)  = NORM(a(:,5), aCnv(:,1)) ! xi_5
      adg(1,3)  = NORM(a(:,6), aCnv(:,1)) ! xi_6

      adg(2,1)  = NORM(a(:,4), aCnv(:,2)) ! eta_4
      adg(2,2)  = NORM(a(:,5), aCnv(:,2)) ! eta_5
      adg(2,3)  = NORM(a(:,6), aCnv(:,2)) ! eta_6

      adg(3,1)  = NORM(a(:,4), nV)        ! z_4
      adg(3,2)  = NORM(a(:,5), nV)        ! z_5
      adg(3,3)  = NORM(a(:,6), nV)        ! z_6

!     Xi matrix (reference config)
      xi0(:,:)  = adg0(:,:)
      xi0(1,3)  = adg0(1,3) + 1._RKIND    ! xi_6
      xi0(2,1)  = adg0(2,1) + 1._RKIND    ! eta_4

!     Xi matrix (current config)
      xi(:,:)   = adg(:,:)
      xi(1,3)   = adg(1,3) + 1._RKIND     ! xi_6
      xi(2,1)   = adg(2,1) + 1._RKIND     ! eta_4

!     Tmat and inverse (reference config)
      DO i=1, 3
         Tm0(i,1) = xi0(1,i)*(xi0(1,i) - 1._RKIND) ! xi**2 - xi
         Tm0(i,2) = xi0(2,i)*(xi0(2,i) - 1._RKIND) ! eta**2 - eta
         Tm0(i,3) = xi0(1,i)*xi0(2,i)              ! xi * eta
      END DO
      Tm0 = MAT_INV(Tm0, 3)

      ! IF (e.EQ.1) THEN
      !    PRINT *, "POST+++++++++++++++"
      !    PRINT *, x0
      ! END IF 
!     Tmat and inverse (current config)
      DO i=1, 3
         Tm(i,1) = xi(1,i)*(xi(1,i) - 1._RKIND)    ! xi**2 - xi
         Tm(i,2) = xi(2,i)*(xi(2,i) - 1._RKIND)    ! eta**2 - eta
         Tm(i,3) = xi(1,i)*xi(2,i)                 ! xi * eta
      END DO
      Tm = MAT_INV(Tm, 3)

!     v = Inv(T) * z (reference config)
      v0(1) = Tm0(1,1)*xi0(3,1) + Tm0(1,2)*xi0(3,2) + Tm0(1,3)*xi0(3,3)
      v0(2) = Tm0(2,1)*xi0(3,1) + Tm0(2,2)*xi0(3,2) + Tm0(2,3)*xi0(3,3)
      v0(3) = Tm0(3,1)*xi0(3,1) + Tm0(3,2)*xi0(3,2) + Tm0(3,3)*xi0(3,3)
 
!     Curvature coefficients (ref. config)
      bb_0(1,1) = 2._RKIND*v0(1)
      bb_0(2,2) = 2._RKIND*v0(2)
      bb_0(1,2) = v0(3)
      bb_0(2,1) = bb_0(1,2)

!     v = Inv(T) * z (current config)
      v(1) = Tm(1,1)*xi(3,1) + Tm(1,2)*xi(3,2) + Tm(1,3)*xi(3,3)
      v(2) = Tm(2,1)*xi(3,1) + Tm(2,2)*xi(3,2) + Tm(2,3)*xi(3,3)
      v(3) = Tm(3,1)*xi(3,1) + Tm(3,2)*xi(3,2) + Tm(3,3)*xi(3,3)

!     Curvature coefficients (current config)
      bb_x(1,1) = 2._RKIND*v(1)
      bb_x(2,2) = 2._RKIND*v(2)
      bb_x(1,2) = v(3)
      bb_x(2,1) = bb_x(1,2)
      END SUBROUTINE PSHELLBENDCST
!####################################################################

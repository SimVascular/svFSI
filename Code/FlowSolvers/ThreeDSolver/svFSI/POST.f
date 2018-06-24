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

!     General purpose routine for post processing outputs. 
      SUBROUTINE POST(lM, res, lY, lD, outGrp, iEq)
      
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      
      IMPLICIT NONE
      
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=8), INTENT(OUT) :: res(maxnsd,lM%nNo)
      REAL(KIND=8), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER, INTENT(IN) :: outGrp, iEq

      LOGICAL FSIeq
      INTEGER a, Ac, e, i, j, eNoN, g
      REAL(KIND=8) rho, kappa, w, Jac, ksix(nsd,nsd), 
     2   lRes(maxnsd), q(nsd), u(nsd), p, T, ux(nsd,nsd)
      REAL(KIND=8), ALLOCATABLE :: sA(:), sF(:,:), xl(:,:), yl(:,:), 
     2   Nx(:,:), N(:)
      COMPLEX*16 :: eig(nsd)
      
      FSIeq = .FALSE.
      IF (eq(iEq)%phys .EQ. phys_FSI) FSIeq = .TRUE.

!     Since energy flux can be directly calculated from nodal values.
!     Note that if there are more than one domain, we need to rely on
!     element based calculations
      IF (outGrp .EQ. outGrp_eFlx .AND. .NOT.ALLOCATED(dmnId)) THEN
         rho = eq(iEq)%dmn(1)%prop(fluid_density)
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            p  = lY(nsd+1,Ac)
            u  = lY(1:nsd,Ac)
            res(1:nsd,Ac) = (p + 5D-1*rho*NORM(u))*u
         END DO
         RETURN
      END IF

!     Other outputs require more calculations
      eNoN  = lM%eNoN
      ALLOCATE (sA(tnNo), sF(maxnsd,tnNo), xl(nsd,eNoN), 
     2   yl(tDof,eNoN), Nx(nsd,eNoN), N(eNoN))
      sA   = 0D0
      sF   = 0D0
      lRes = 0D0
      eig  = (0D0, 0D0)
      DO e=1, lM%nEl
         cDmn = DOMAIN(lM, iEq, e)
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
               CALL GNN(eNoN, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
            END IF
            w = lM%w(g)*Jac
            N = lM%N(:,g)
!     Vorticity calculation   ---------------------------------------
            IF (outGrp .EQ. outGrp_vort) THEN
               lRes = 0D0
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
               ux   = 0D0
               lRes = 0D0
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
               lRes(2:maxnsd) = 0D0
!     Energy flux calculation   -------------------------------------
            ELSE IF (outGrp .EQ. outGrp_eFlx) THEN
               rho = eq(iEq)%dmn(cDmn)%prop(fluid_density)
               p   = 0D0
               u   = 0D0
               DO a=1, eNoN
                  p = p + N(a)*yl(nsd+1,a)
                  u = u + N(a)*yl(1:nsd,a)
               END DO
               lRes(1:nsd) = (p + 5D-1*rho*NORM(u))*u
!     Heat flux calculation   ---------------------------------------
            ELSE IF (outGrp .EQ. outGrp_hFlx) THEN
               kappa = eq(iEq)%dmn(cDmn)%prop(conductivity)
               i     = eq(iEq)%s
               DO a=1, eNoN
               END DO
               IF (eq(iEq)%phys .EQ. phys_heatF) THEN
                  u = 0D0
                  T = 0D0
                  q = 0D0
                  DO a=1, eNoN
                     q = q + Nx(:,a)*yl(i,a)
                     u = u + N(a)*yl(1:nsd,a)
                     T = T + N(a)*yl(i,a)
                  END DO
                  lRes(1:nsd) = u*T - kappa*q
               ELSE
                  q = 0D0
                  DO a=1, eNoN
                     q = q + Nx(:,a)*yl(i,a)
                  END DO
                  lRes(1:nsd) = -kappa*q
               END IF
!     Strain tensor invariants calculation   ------------------------
            ELSE IF (outGrp .EQ. outGrp_stInv) THEN
               ksix = 0D0
               DO a=1, eNoN
                  ksix(1,1) = ksix(1,1) + Nx(1,a)*yl(1,a)
                  ksix(2,1) = ksix(2,1) +(Nx(1,a)*yl(2,a) 
     2                                  + Nx(2,a)*yl(1,a))/2D0
                  ksix(2,2) = ksix(2,2) + Nx(2,a)*yl(2,a)
                  IF (nsd .EQ. 3) THEN
                     ksix(3,1) = ksix(3,1) +(Nx(1,a)*yl(3,a) 
     2                                     + Nx(3,a)*yl(1,a))/2D0
                     ksix(3,2) = ksix(3,2) +(Nx(2,a)*yl(3,a) 
     2                                     + Nx(3,a)*yl(2,a))/2D0
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
     2                    + ksix(2,1)*ksix(3,2)*ksix(3,1)*2D0
     3                    - ksix(1,1)*ksix(3,2)*ksix(3,2)
     4                    - ksix(3,1)*ksix(2,2)*ksix(3,1)
     5                    - ksix(2,1)*ksix(2,1)*ksix(3,3)
               END IF
               lRes = ABS(lRes)
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
      
!####################################################################
!     General purpose routine for post processing outputs at the
!     faces. Currently this calculats WSS, which is t.n - (n.t.n)n
!     Here t is stress tensor: t = \mu (grad(u) + grad(u)^T)
      SUBROUTINE BPOST(lM, res, lY, lD, outGrp)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=8), INTENT(OUT) :: res(maxnsd,lM%nNo)
      REAL(KIND=8), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER, INTENT(IN) :: outGrp

      LOGICAL FSIeq
      INTEGER a, Ac, e, Ec, i, j, iEq, iFa, eNoN, g
      REAL(KIND=8) Tdn(nsd), ndTdn, taue(nsd), ux(nsd,nsd), mu, w,
     2   nV(nsd), Jac, ksix(nsd,nsd), lRes(maxnsd)
      REAL(KIND=8), ALLOCATABLE :: sA(:), sF(:,:), gnV(:,:), lnV(:,:),
     2   xl(:,:), ul(:,:), Nx(:,:), N(:)
      
      iEq   = 1
      eNoN  = lM%eNoN
      FSIeq = .FALSE.
      IF (eq(iEq)%phys .EQ. phys_FSI) FSIeq = .TRUE.

      ALLOCATE (sA(tnNo), sF(maxnsd,tnNo), gnV(nsd,tnNo), 
     2   lnV(nsd,eNoN), xl(nsd,eNoN), ul(nsd,eNoN), Nx(nsd,eNoN),
     3   N(eNoN))
      sA   = 0D0
      sF   = 0D0
      gnV  = 0D0
      lRes = 0D0
!     First creating the norm field
      DO iFa=1, lM%nFa
         DO a=1, lM%fa(iFa)%nNo
            Ac        = lM%fa(iFa)%gN(a)
            gnV(:,Ac) = lM%fa(iFa)%nV(:,a)
         END DO
      END DO

      DO iFa=1, lM%nFa
         DO e=1, lM%fa(iFa)%nEl
            Ec = lM%fa(iFa)%gE(e)
            cDmn = DOMAIN(lM, iEq, Ec)
            IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, Ec)
!     Finding the norm for all the nodes of this element, including
!     those that don't belong to this face, which will be inerpolated
!     from the nodes of the face
            nV = 0D0
            mu = eq(iEq)%dmn(cDmn)%prop(viscosity)
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
            nV = nV/lM%fa(iFa)%eNoN
            DO a=1, eNoN
               IF (ISZERO(NORM(lnV(:,a)))) lnV(:,a) = nV
            END DO
            DO g=1, lM%nG
               IF (g.EQ.1 .OR. .NOT.lM%lShpF) THEN
                  CALL GNN(eNoN, lM%Nx(:,:,g), xl, Nx, Jac, ksix)
               END IF
               w = lM%w(g)*Jac
               N = lM%N(:,g)

               IF (outGrp .EQ. outGrp_WSS) THEN
!     Calculating ux = grad(u) and nV at a Gauss point
                  ux = 0D0
                  nV = 0D0
                  DO a=1, eNoN
                     nV = nV + N(a)*lnV(:,a)
                     DO i=1, nsd
                        DO j=1, nsd
                           ux(i,j) = ux(i,j) + Nx(i,a)*ul(j,a)
                        END DO
                     END DO
                  END DO

!     Now finding grad(u).n and n.gard(u).n
                  Tdn   = 0D0
                  ndTdn = 0D0
                  DO i=1,nsd
                     DO j=1, nsd
                        Tdn(i) = Tdn(i) + mu*(ux(i,j) + ux(j,i))*nV(j)
                     END DO
                     ndTdn = ndTdn + Tdn(i)*nV(i)
                  END DO
                  taue = Tdn - ndTdn*nV
                  lRes(1:nsd) = -taue
               ELSE
                  err = "Correction is required in BPOST"
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
               sA(Ac)   = 1D0
            END IF
         END DO
      END DO

      res = 0D0
      DO a=1, lM%nNo
         Ac       = lM%gN(a)
         res(:,a) = sF(:,Ac)
      END DO

      RETURN
      END SUBROUTINE BPOST

!####################################################################
!     This is a routine that performs POST/BPOST on all meshes
      SUBROUTINE ALLPOST(res, lY, lD, outGrp, iEq)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=8), INTENT(OUT) :: res(maxnsd,tnNo)
      REAL(KIND=8), INTENT(IN) :: lY(tDof,tnNo), lD(tDof,tnNo)
      INTEGER, INTENT(IN) :: outGrp, iEq

      INTEGER a, Ac, iM
      REAL(KIND=8), ALLOCATABLE :: tmpV(:,:)

      DO iM=1, nMsh
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))
         IF (outGrp .EQ. outGrp_WSS) THEN
            CALL BPOST(msh(iM), tmpV, lY, lD, outGrp)
         ELSE
            CALL POST(msh(iM), tmpV, lY, lD, outGrp, iEq)
         END IF
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            res(:,Ac) = tmpV(:,a)
         END DO
      END DO

      RETURN
      END SUBROUTINE ALLPOST

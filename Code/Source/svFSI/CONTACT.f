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
!     This routine determines which contact model would be used in
!     the surface contact calculation.
!
!--------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_CONTACT(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      SELECT CASE (cntctM%cType)
      CASE (cntctM_penalty)
         CALL CONSTRUCT_CONTACTPNLTY(Dg)
      CASE (cntctM_potential)
         CALL CONSTRUCT_CONTACTPTNL(Dg)
      CASE DEFAULT
         err = "Undefined contact model"
      END SELECT

      END SUBROUTINE CONSTRUCT_CONTACT

!-----------------------------------------------------------------------
!
!     This routine applies penalty-based contact model for possible
!     contacting shell surfaces. 2-pass
!
!--------------------------------------------------------------------  

      SUBROUTINE CONSTRUCT_CONTACTPNLTY(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: i, j, k, l, m, iM, jM, e, a, Ac, b, Bc, g,
     2   eNoN, nNb, insd, maxNnb
      REAL(KIND=RKIND) :: kl, hl, w, Jac, al, c, d, pk, nV1(nsd),
     2   nV2(nsd), x1(nsd), x2(nsd), x12(nsd), xmin(nsd), xmax(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), bBox(:,:), lRBI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), N(:), Nx(:,:),
     2   gCov(:,:), gCnv(:,:), xl(:,:), lR(:,:), lRBc(:,:,:)

      IF (eq(cEq)%phys .NE. phys_shell) RETURN

      i  = eq(cEq)%s
      j  = i + 1
      k  = j + 1
      kl = cntctM%k
      hl = cntctM%h

!     Compute normal vectors at each node in the current configuration
      ALLOCATE(sF(nsd,tnNo), sA(tnNo))
      sF = 0._RKIND
      sA = 0._RKIND
      DO iM=1, nMsh
         IF (.NOT.msh(iM)%lShl) CYCLE
         eNoN = msh(iM)%eNoN
         insd = nsd - 1
         ALLOCATE(Nx(insd,eNoN), N(eNoN), xl(nsd,eNoN), gCov(nsd,insd),
     2      gCnv(nsd,insd))
         Nx = msh(iM)%Nx(:,:,1)
         DO e=1, msh(iM)%nEl
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,e)
               xl(1,a) = x(1,Ac) + Dg(i,Ac)
               xl(2,a) = x(2,Ac) + Dg(j,Ac)
               xl(3,a) = x(3,Ac) + Dg(k,Ac)
            END DO
            CALL GNNS(eNoN, Nx, xl, nV1, gCov, gCnv)
            Jac = SQRT(NORM(nV1(:))) 
            nV1(:) = nV1(:)/Jac
            IF (msh(iM)%lSFp) nV1 = -nV1

            DO g=1, msh(iM)%nG
               w = msh(iM)%w(g) * Jac
               N(:) = msh(iM)%N(:,g)
               DO a=1, eNoN
                  Ac = msh(iM)%IEN(a,e)
                  sA(Ac) = sA(Ac) + w*N(a)
                  sF(:,Ac) = sF(:,Ac) + w*N(a)*nV1(:)
               END DO
            END DO
         END DO
         DEALLOCATE(N, Nx, xl, gCov, gCnv)
      END DO

      CALL COMMU(sF)
      CALL COMMU(sA)

      DO Ac=1, tnNo
         IF (.NOT.ISZERO(sA(Ac))) sF(:,Ac) = sF(:,Ac)/sA(Ac)
         Jac = SQRT(SUM(sF(:,Ac)**2))
         IF (.NOT.ISZERO(Jac)) sF(:,Ac) = sF(:,Ac) / Jac
      END DO

!     Create a bounding box around possible region of contact and bin
!     the box with neighboring nodes
      maxNnb = 15
 101  maxNnb = maxNnb + 5
      IF (ALLOCATED(bBox)) DEALLOCATE(bBox)
      ALLOCATE(bBox(maxNnb,tnNo))
      bBox = 0
      DO iM=1, nMsh
         IF (.NOT.msh(iM)%lShl) CYCLE
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            x1(1) = x(1,Ac) + Dg(i,Ac)
            x1(2) = x(2,Ac) + Dg(j,Ac)
            x1(3) = x(3,Ac) + Dg(k,Ac)

!           Box limits for each node
            xmin(:) = x1(:) - cntctM%c
            xmax(:) = x1(:) + cntctM%c

!           Load the box with neighboring nodes lying within it
            DO jM=1, nMsh
               IF (iM.EQ.jM .OR. .NOT.msh(jM)%lShl) CYCLE
               DO b=1, msh(jM)%nNo
                  Bc = msh(jM)%gN(b)
                  x2(1) = x(1,Bc) + Dg(i,Bc)
                  x2(2) = x(2,Bc) + Dg(j,Bc)
                  x2(3) = x(3,Bc) + Dg(k,Bc)
                  IF (x2(1).GE.xmin(1) .AND. x2(1).LE.xmax(1) .AND.
     2                x2(2).GE.xmin(2) .AND. x2(2).LE.xmax(2) .AND.
     3                x2(3).GE.xmin(3) .AND. x2(3).LE.xmax(3)) THEN
                     DO l=1, maxNnb
                        IF (bBox(l,Ac) .EQ. 0) THEN
                           bBox(l,Ac) = Bc
                           EXIT
                        END IF
                        IF (Bc .GT. bBox(l,Ac)) CYCLE
                        IF (Bc .EQ. bBox(l,Ac)) EXIT
                        IF (bBox(maxNnb,Ac) .NE. 0) GOTO 101
                        DO m=maxNnb, l+1, -1
                           bBox(m,Ac) = bBox(m-1,Ac)
                        END DO
                        bBox(l,Ac) = Bc
                        EXIT
                     END DO
                     IF (l .GT. maxNnb) GOTO 101
                  END IF
               END DO ! b
            END DO ! jM
         END DO ! a
      END DO ! iM

!     Check if any node is strictly involved in contact and compute
!     corresponding penalty forces assembled to the residue
      ALLOCATE(lR(dof,tnNo), incNd(tnNo), lRBc(dof,maxNnb,tnNo), 
     2                lRBI(maxNnb,tnNo))
      lR    = 0._RKIND
      incNd = 0
      lRBc  = 0._RKIND
      lRBI  = 0
      DO Ac=1, tnNo
         IF (bBox(1,Ac) .EQ. 0) CYCLE
         x1(1)  = x(1,Ac) + Dg(i,Ac)
         x1(2)  = x(2,Ac) + Dg(j,Ac)
         x1(3)  = x(3,Ac) + Dg(k,Ac)
         nV1(:) = sF(:,Ac)
         nNb    = 0
         DO a=1, maxNnb
            Bc = bBox(a,Ac)
            IF (Bc .EQ. 0) CYCLE
            x2(1)  = x(1,Bc) + Dg(i,Bc)
            x2(2)  = x(2,Bc) + Dg(j,Bc)
            x2(3)  = x(3,Bc) + Dg(k,Bc)
            nV2(:) = sF(:,Bc)

            x12 = x2(:) - x1(:)
            c   = SQRT(NORM(x12))
            al  = SQRT(ABS(NORM(nV1, nV2)))

            IF (c.LE.cntctM%c .AND. al.GE.cntctM%al) THEN
               d = NORM(x12, nV2)
               flag = .FALSE.
               IF (d.GE.-cntctM%h .AND. d.LT.0._RKIND) THEN
                  pk = 0.5_RKIND*kl/hl * (d + hl)**2._RKIND
                  ! pk = 0.5_RKIND*kl/hl * (EXP(d + hl) - 1._RKIND)
                  flag = .TRUE.
               ELSE IF (d .GE. 0._RKIND) THEN
                  pk = 0.5_RKIND*kl * (hl + 2._RKIND*d)
                  ! pk = 0.5_RKIND*kl/hl * (EXP(d + hl) - 1._RKIND)
                  flag = .TRUE.
               ELSE
                  pk = 0._RKIND
               END IF
               IF (flag) THEN
                  incNd(Ac) = 1
                  nNb = nNb + 1
                  lR(1:nsd,Ac) = lR(1:nsd,Ac) - pk*nV2(:)
                  lRBI(nNb,Ac) = Bc
                  lRBc(1:nsd,nNb,Ac) = pk*nV2(:)
               END IF
            END IF
         END DO
         IF (nNb .NE. 0) lR(:,Ac) = lR(:,Ac) / REAL(nNb, KIND=RKIND)
         lR(:,Ac) = lR(:,Ac) / 2._RKIND
         DO b = 1, maxNnb
            IF (lRBI(b,Ac).EQ.0) CYCLE
            lR(1:nsd,lRBI(b,Ac)) = lR(1:nsd,lRBI(b,Ac)) + 
     2          lRBc(1:nsd,b,Ac) / REAL(nNb, KIND=RKIND) / 2._RKIND
         END DO
      END DO
      DEALLOCATE(sA, sF, bBox)

!     Return if no penalty forces are to be added
      IF (SUM(incNd) .EQ. 0) RETURN

!     Global assembly
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 0) CYCLE
         R(:,Ac) = R(:,Ac) + lR(:,Ac)
      END DO
      DEALLOCATE(lR, incNd)

      RETURN
      END SUBROUTINE CONSTRUCT_CONTACTPNLTY

!####################################################################
!     This is the potential-based contact model; referring to  
!     Kamensky et al. 2018, 2-half-pass 
!--------------------------------
      SUBROUTINE CONSTRUCT_CONTACTPTNL(Dg)
      USE ALLFUN
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: i, j, k, l, m, iM, jM, e, a, b, Ac, Bc,
     2   eNoN, insd, maxNnb, eNoN1, eNoN2, e1, e2, gnEl, giEl, gjEl,
     3   g1, g2
      REAL(KIND=RKIND) :: Jac, nV1(nsd), xmin(nsd), Jac0, nV0(nsd),
     2   xmax(nsd), xa(nsd), xb(nsd), xab(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: bBox(:,:), ptr1(:),
     2    ptr2(:), sI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: sJ0(:),sF(:,:), Nx(:,:),
     2   gCov(:,:), gCnv(:,:), xc(:,:), x1(:,:), x2(:,:), x0(:,:), 
     3   lR(:,:), lR2(:,:), lK1(:,:,:), lK2(:,:,:)

      IF (eq(cEq)%phys .NE. phys_shell) RETURN

      i   = eq(cEq)%s
      j   = i + 1
      k   = j + 1

!     Get global information of each element 
!     Using global element id might cause problems in para. comput.
!     Try and consider to change it to use surf. id + elem. id
      gnEl = 0
      DO iM=1, nMsh
         gnEl = gnEl + msh(iM)%nEl
      END DO
      giEl = 0
      ALLOCATE(sF(nsd,gnEl), sI(2,gnEl), sJ0(gnEl))
      sF = 0._RKIND
      sI = 0
      sJ0 = 0._RKIND
      DO iM=1, nMsh
         IF (.NOT.msh(iM)%lShl) CYCLE
         eNoN = msh(iM)%eNoN
         insd = nsd - 1
         ALLOCATE(Nx(insd,eNoN), gCov(nsd,insd), 
     2      gCnv(nsd,insd), x0(nsd,eNoN), xc(nsd,eNoN))
         Nx = msh(iM)%Nx(:,:,1)
         DO e=1, msh(iM)%nEl
            giEl = giEl + 1
            DO a=1, eNoN
               Ac = msh(iM)%IEN(a,e)
               x0(1,a) = x(1,Ac)
               x0(2,a) = x(2,Ac)
               x0(3,a) = x(3,Ac)
               xc(1,a) = x(1,Ac) + Dg(i,Ac)
               xc(2,a) = x(2,Ac) + Dg(j,Ac)
               xc(3,a) = x(3,Ac) + Dg(k,Ac)
            END DO
            CALL GNNS(eNoN, Nx, x0, nV0, gCov, gCnv)
            Jac0 = SQRT(NORM(nV0(:))) 
            nV0(:) = nV0(:)/Jac0
            sJ0(giEl) = Jac0
            CALL GNNS(eNoN, Nx, xc, nV1, gCov, gCnv)
            Jac = SQRT(NORM(nV1(:))) 
            nV1(:) = nV1(:)/Jac  
            IF (msh(iM)%lSFp) nV1 = -nV1
            sF(:,giEl) = nV1(:) ! Element normal
            sI(1,giEl) = iM 
            sI(2,giEl) = e
         END DO
         DEALLOCATE(Nx, x0, xc, gCov, gCnv)
      END DO

!     Create a bounding box around possible region of contact and bin
!     the box with neighboring elements
      maxNnb = 15
 101  maxNnb = maxNnb + 5

      IF (ALLOCATED(bBox)) DEALLOCATE(bBox)
      ALLOCATE(bBox(maxNnb,gnEl))
      bBox = 0   
      
      giEl = 0
      DO iM=1, nMsh 
         IF (.NOT.msh(iM)%lShl) CYCLE
         IF (ALLOCATED(x1)) DEALLOCATE(x1)
         ALLOCATE(x1(nsd,msh(iM)%eNoN))
         ! Calculate the centroid for ek
         DO e1=1, msh(iM)%nEl
            giEl = giEl + 1 
            xa = 0._RKIND   ! Centroid on surf1
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e1)
               x1(1,a) = x(1,Ac) + Dg(i,Ac)
               x1(2,a) = x(2,Ac) + Dg(j,Ac)
               x1(3,a) = x(3,Ac) + Dg(k,Ac)
               
               xa(:) = xa(:) + x1(:,a)
            END DO
            xa = xa / REAL(msh(iM)%eNoN, KIND=RKIND)
!           Box limits for each element
            xmin(:) = xa(:) - cntctM%Rout
            xmax(:) = xa(:) + cntctM%Rout  
            
!           Load the box with neighboring elements lying within it
            gjEl = 0
            DO jM=1, nMsh
               ! No self contact
               IF (.NOT.msh(jM)%lShl) CYCLE 
               IF (iM.EQ.jM) THEN
                  gjEl = gjEl + msh(jM)%nEl
                  CYCLE
               END IF
               IF (ALLOCATED(x2)) DEALLOCATE(x2)
               ALLOCATE(x2(nsd,msh(jM)%eNoN))
               
               ! Calculate the centroid for el
               DO e2=1, msh(jM)%nEl
                  gjEl = gjEl + 1
                  xb = 0._RKIND   ! Centroid on surf2
                  DO b=1, msh(jM)%eNoN
                     Bc = msh(jM)%IEN(b,e2)
                     x2(1,b) = x(1,Bc) + Dg(i,Bc)
                     x2(2,b) = x(2,Bc) + Dg(j,Bc)
                     x2(3,b) = x(3,Bc) + Dg(k,Bc)
                     
                     xb(:) = xb(:) + x2(:,b)
                  END DO
                  xb = xb / REAL(msh(jM)%eNoN, KIND=RKIND)
       
                  xab = xb(:) - xa(:)

                  l = 0
                  IF (xb(1).GE.xmin(1) .AND. xb(1).LE.xmax(1) .AND.
     2                xb(2).GE.xmin(2) .AND. xb(2).LE.xmax(2) .AND.
     3                xb(3).GE.xmin(3) .AND. xb(3).LE.xmax(3)) THEN
                     DO l=1, maxNnb
                        IF (bBox(l,giEl) .EQ. 0) THEN
                           bBox(l,giEl) = gjEl
                           EXIT
                        END IF
                        IF (gjEl .GT. bBox(l,giEl)) CYCLE
                        IF (gjEl .EQ. bBox(l,giEl)) EXIT
                        IF (bBox(maxNnb,giEl) .NE. 0) GOTO 101
                        DO m=maxNnb, l+1, -1
                           bBox(m,giEl) = bBox(m-1,giEl)
                        END DO
                           
                        bBox(l,giEl) = gjEl
                        EXIT
                     END DO
                  END IF
                  IF (l .GT. maxNnb) GOTO 101
               END DO ! gjEl
               DEALLOCATE(x2)
            END DO ! jM   
         END DO ! giEl
         DEALLOCATE(x1)
      END DO ! iM 

!     Compute the corresponding local residual and stiffness for 
!     each elem
      
      DO giEl=1, gnEl
         IF (bBox(1,giEl) .EQ. 0) CYCLE
         iM = sI(1,giEl)
         eNoN1 = msh(iM)%eNoN
         ALLOCATE(x1(nsd,eNoN1), ptr1(eNoN1), lR(dof,eNoN1),
     2           lK1(dof*dof, eNoN1, eNoN1)) 
         lR = 0._RKIND
         lK1= 0._RKIND 
         e1 = sI(2,giEl)
         xa = 0._RKIND   ! Centroid on surf1
         DO a=1, eNoN1
            Ac = msh(iM)%IEN(a,e1)
            ptr1(a) = Ac
            x1(1,a) = x(1,Ac) + Dg(i,Ac)
            x1(2,a) = x(2,Ac) + Dg(j,Ac)
            x1(3,a) = x(3,Ac) + Dg(k,Ac)

            xa(:) = xa(:) + x1(:,a)
         END DO
         xa = xa / REAL(eNoN1, KIND=RKIND)

         DO l=1, maxNnb
            gjEl = bBox(l,giEl)
            IF (gjEl .EQ. 0) CYCLE
            jM = sI(1,gjEl)
            eNoN2 = msh(jM)%eNoN
            ALLOCATE(x2(nsd,eNoN2),ptr2(eNoN2),lR2(dof,eNoN2),
     2          lK2(dof*dof, eNoN1, eNoN2)) 
            lR2  = 0._RKIND     
            lK2  = 0._RKIND
            e2 = sI(2,gjEl)
            xb = 0._RKIND   ! Centroid on surf2
            DO b=1, eNoN2
               Bc = msh(jM)%IEN(b,e2)
               ptr2(b) = Bc
               x2(1,b) = x(1,Bc) + Dg(i,Bc)
               x2(2,b) = x(2,Bc) + Dg(j,Bc)
               x2(3,b) = x(3,Bc) + Dg(k,Bc)
               
               xb(:) = xb(:) + x2(:,b)
            END DO
            xb = xb / REAL(eNoN2, KIND=RKIND)
            DO g1=1, msh(jM)%nG
               DO g2=1, msh(jM)%nG
                  CALL CONTACTPTNL3D(msh(iM), msh(jM), eNoN1, 
     2                  eNoN2, g1, g2, xa, xb, sJ0(giEl), sJ0(gjEl),
     3                  sF(:,giEl), sF(:,gjEl), lR, lK1, lK2, 
     4                  giEl, gjEl)
                  ! To delete giEl gjEl; for debug
               END DO
            END DO
            ! Assembly lK2
            ! Note that here we only consider the same kinds of elems
            ! are used for both surfaces. If different kinds of elems
            ! are used, this needs to be modified.
            lK2 = 0._RKIND
#ifdef WITH_TRILINOS
            IF (eq(cEq)%assmTLS) THEN
               CALL TRILINOS_DOASSEM(eNoN2, ptr2, lK2, lR2)
            ELSE
#endif
               CALL DOASSEM(eNoN2, ptr2, lK2, lR2)
#ifdef WITH_TRILINOS
            END IF
#endif 
            DEALLOCATE(x2, ptr2, lR2, lK2) 
         END DO 
         ! Assembly lR lK1
         lK1 = 0._RKIND
#ifdef WITH_TRILINOS
         IF (eq(cEq)%assmTLS) THEN
            CALL TRILINOS_DOASSEM(eNoN1, ptr1, lK1, lR)
         ELSE
#endif
            CALL DOASSEM(eNoN1, ptr1, lK1, lR)
#ifdef WITH_TRILINOS
         END IF
#endif 
         DEALLOCATE(x1, ptr1, lR, lK1)     
      END DO
      DEALLOCATE(sJ0, sF, sI, bBox)

      RETURN
      END SUBROUTINE CONSTRUCT_CONTACTPTNL
!----------------------------------------------------------------------     
      SUBROUTINE CONTACTPTNL3D (lM1, lM2, eNoN1, eNoN2, g1, g2, xa, xb,
     2                   J10, J20, nV1, nV2, lR, lK1, lK2, giEl, gjEl)  
      USE COMMOD
      USE ALLFUN
      USE MATFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM1, lM2
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN1, eNoN2, g1, g2, giEl, 
     2                gjEl
      REAL(KIND=RKIND), INTENT(IN) ::  xa(nsd), xb(nsd), J10, J20, 
     2                nV1(3), nV2(3)                           
      REAL(KIND=RKIND), INTENT(INOUT) :: lR(dof,eNoN1), 
     2   lK1(dof*dof,eNoN1,eNoN1), lK2(dof*dof,eNoN1,eNoN2)

      INTEGER(KIND=IKIND) :: a, b  
      REAL(KIND=RKIND) :: fs, nk(3), gk, c1, c2, pl, kl, Rin, Rout, w12,
     2   dfs, fsij(3,3), gap
      REAL(KIND=RKIND), ALLOCATABLE :: N1(:), N2(:)

      ALLOCATE(N1(eNoN1), N2(eNoN2))

      pl = cntctM%p
      kl = cntctM%k
      Rin = cntctM%Rin
      Rout = cntctM%Rout

      nk = xa - xb
      gk = SQRT(NORM(nk))
      ! From the paper; to rethink the nk definition
      IF (NORM(nk,nV2).LT.0) gk = -gk
      nk = nk/gk

      ! Contact model
      c1 = 0.5_RKIND*pl*kl/(Rout - Rin)*(Rin**(-pl-1._RKIND))
      c2 = kl*(Rin**(-pl)) - c1*(Rin - Rout)*(Rin - Rout)
      fs = 0._RKIND
      dfs = 0._RKIND
      IF (gk < Rin) THEN
         fs = kl*(gk**(-pl)) - c2 
         dfs = - pl*kl*(gk**(-pl-1._RKIND))
      ELSE IF (gk < Rout) THEN
         fs = c1*(gk - Rout)*(gk - Rout)
         dfs = 2._RKIND*c1*(gk - Rout)
      END IF
      fsij = fs/gk*MAT_ID(3) + (dfs-fs/gk)*MAT_DYADPROD(nk,nk,3)

      w12 = lM1%w(g1)*lM2%w(g2)*J10*J20
      N1 = lM1%N(:,g1)
      N2 = lM2%N(:,g2)
      DO a=1, eNoN1
         lR(1,a) = lR(1,a) - w12*N1(a)*fs*nk(1)
         lR(2,a) = lR(2,a) - w12*N1(a)*fs*nk(2)
         lR(3,a) = lR(3,a) - w12*N1(a)*fs*nk(3)
      END DO
      
      DO a=1, eNoN1
         DO b=1, eNoN1
            lK1(1,a,b) = lK1(1,a,b) + w12*N1(a)*N1(b)*fsij(1,1)
            lK1(2,a,b) = lK1(2,a,b) + w12*N1(a)*N1(b)*fsij(1,2)
            lK1(3,a,b) = lK1(3,a,b) + w12*N1(a)*N1(b)*fsij(1,3)
            lK1(4,a,b) = lK1(4,a,b) + w12*N1(a)*N1(b)*fsij(2,1)
            lK1(5,a,b) = lK1(5,a,b) + w12*N1(a)*N1(b)*fsij(2,2)
            lK1(6,a,b) = lK1(6,a,b) + w12*N1(a)*N1(b)*fsij(2,3)
            lK1(7,a,b) = lK1(7,a,b) + w12*N1(a)*N1(b)*fsij(3,1)
            lK1(8,a,b) = lK1(8,a,b) + w12*N1(a)*N1(b)*fsij(3,2)
            lK1(9,a,b) = lK1(9,a,b) + w12*N1(a)*N1(b)*fsij(3,3)
         END DO
      END DO
      DO a=1, eNoN1
         DO b=1, eNoN2
            lK2(1,a,b) = lK2(1,a,b) - w12*N1(a)*N2(b)*fsij(1,1)
            lK2(2,a,b) = lK2(2,a,b) - w12*N1(a)*N2(b)*fsij(1,2)
            lK2(3,a,b) = lK2(3,a,b) - w12*N1(a)*N2(b)*fsij(1,3)
            lK2(4,a,b) = lK2(4,a,b) - w12*N1(a)*N2(b)*fsij(2,1)
            lK2(5,a,b) = lK2(5,a,b) - w12*N1(a)*N2(b)*fsij(2,2)
            lK2(6,a,b) = lK2(6,a,b) - w12*N1(a)*N2(b)*fsij(2,3)
            lK2(7,a,b) = lK2(7,a,b) - w12*N1(a)*N2(b)*fsij(3,1)
            lK2(8,a,b) = lK2(8,a,b) - w12*N1(a)*N2(b)*fsij(3,2)
            lK2(9,a,b) = lK2(9,a,b) - w12*N1(a)*N2(b)*fsij(3,3)
         END DO
      END DO
      DEALLOCATE(N1, N2)
      
      END SUBROUTINE CONTACTPTNL3D
!####################################################################
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
!     This routine applies penalty-based contact model for possible
!     contacting shell surfaces.
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

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), bBox(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: sA(:), sF(:,:), N(:), Nx(:,:),
     2   gCov(:,:), gCnv(:,:), xl(:,:), lR(:,:)

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
      ALLOCATE(lR(dof,tnNo), incNd(tnNo))
      lR    = 0._RKIND
      incNd = 0
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

            x12 = x1(:) - x2(:)
            c   = SQRT(NORM(x12))
            al  = SQRT(ABS(NORM(nV1, nV2)))

            IF (c.LE.cntctM%c .AND. al.GE.cntctM%al) THEN
               d = NORM(x12, nV2)
               flag = .FALSE.
               IF (d.GE.-cntctM%h .AND. d.LT.0._RKIND) THEN
                  pk = 0.5_RKIND*kl/hl * (d + hl)**2._RKIND
                  flag = .TRUE.
               ELSE IF (d .GE. 0._RKIND) THEN
                  pk = 0.5_RKIND*kl * (hl + d)
                  flag = .TRUE.
               ELSE
                  pk = 0._RKIND
               END IF
               IF (flag) THEN
                  incNd(Ac) = 1
                  nNb = nNb + 1
                  lR(1:nsd,Ac) = lR(1:nsd,Ac) - pk*nV1(:)
               END IF
            END IF
            ! IF (Ac.EQ.40 .AND. Bc.EQ.627) THEN
            !    PRINT *, Ac, Bc, lR(:,Ac)
            ! ELSE IF (Ac.EQ.301 .AND. Bc.EQ.546) THEN
            !    PRINT *, Ac, Bc, d
            ! END IF
            ! IF (NORM(lR(:, Ac)) > 0.00000001_RKIND) THEN
            !    PRINT *, Ac, Bc, lR(:,Ac)
            ! END IF
         END DO
         IF (nNb .NE. 0) lR(:,Ac) = lR(:,Ac) / REAL(nNb, KIND=RKIND)
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

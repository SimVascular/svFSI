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
!     This routine applies finite cell method for immersed boundaries.
!
!--------------------------------------------------------------------

      SUBROUTINE FC_INIT(lFc, lM, xl)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      TYPE(fCellType), INTENT(INOUT) :: lFc
      REAL(KIND=8), INTENT(IN) :: xl(nsd,lM%eNoN)

      lFc%nSub  = 0
      lFc%ilev  = 0
      lFc%eType = lM%eType
      lFc%eNoN  = lM%eNoN
      lFc%nG    = lM%nG
      ALLOCATE(lFc%x(nsd,lFc%eNoN), lFc%xi(nsd,lFc%eNoN),
     2   lFc%incG(lFc%nG), lFc%w(lFc%nG), lFc%xiGP(nsd,lFc%nG))
      lFc%incG(:) = .TRUE.
      lFc%x  = xl
      lFc%w  = lM%w
      lFc%xiGP = lM%xi
      CALL FC_GETXI(nsd, lFc%eType, lFc%eNoN, lFc%xi)

      RETURN
      END SUBROUTINE FC_INIT
!####################################################################
!     Returns reference coordinates
      PURE SUBROUTINE FC_GETXI(insd, eType, eNoN, xi)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: insd, eType, eNoN
      REAL(KIND=8), INTENT(OUT) :: xi(insd,eNoN)

      REAL(KIND=8) s, t

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_BRK)
         s =  1D0
         t = -1D0
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t

!     2D elements
      CASE(eType_BIL)
         s =  1D0
         t = -1D0
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t

!     1D elements
      CASE(eType_LIN)
         s = 1D0
         xi(1,1) = -s
         xi(1,2) =  s

      END SELECT

      RETURN
      END SUBROUTINE FC_GETXI
!####################################################################
      RECURSIVE SUBROUTINE FC_SET(fCell, lM, xe, Ug)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(mshType), INTENT(IN) :: lM
      TYPE(fCellType), INTENT(INOUT) :: fCell
      REAL(KIND=8), INTENT(IN) :: xe(nsd,fCell%eNoN), Ug(nsd,ib%tnNo)

      INTEGER, PARAMETER :: maxFClev = 2
      INTEGER :: a, i, g, nIn, nOut, eNoN, nG
      REAL(KIND=8) :: xp(nsd)
      LOGICAL :: flag

      REAL(KIND=8), ALLOCATABLE :: xGP(:,:)

      eNoN = fCell%eNoN
      nG   = fCell%nG
      nIn  = 0
      nOut = 0
      DO a=1, eNoN
         xp = fCell%x(:,a)
         flag = .FALSE.
         CALL IB_CHECKINOUT(xp, Ug, flag)
         IF (flag) THEN
            nIn = nIn + 1
         ELSE
            nOut = nOut + 1
         END IF
      END DO

      IF (nIn.EQ.0 .OR. nOut.EQ.0 .OR. fCell%ilev.EQ.maxFClev) THEN
         IF (.NOT.ALLOCATED(fCell%xiGP)) ALLOCATE(fCell%xiGP(nsd,nG))
         fCell%nSub = 0
         fCell%xiGP = 0D0

         ALLOCATE(xGP(nsd,nG))
         xGP = 0D0
         DO g=1, nG
            DO a=1, eNoN
               xGP(:,g) = xGP(:,g) + lM%N(a,g)*fCell%x(:,a)
               fCell%xiGP(:,g) = fCell%xiGP(:,g) +
     2            lM%N(a,g)*fCell%xi(:,a)
            END DO
         END DO

         IF (nIn .EQ. 0) THEN
            fCell%incG = .TRUE.
         ELSE IF (nOut .EQ. 0) THEN
            fCell%incG = .FALSE.
         ELSE IF (fCell%ilev .EQ. maxFClev) THEN
            DO g=1, nG
               flag = .FALSE.
               CALL IB_CHECKINOUT(xGP(:,g), Ug, flag)
               IF (flag) fCell%incG(g) = .FALSE.
            END DO
         END IF
         DEALLOCATE(xGP)
      ELSE
         CALL FC_SUBDIVIDE(fCell)
         DO i=1, fCell%nSub
            CALL FC_SET(fCell%sub(i), lM, xe, Ug)
         END DO
      END IF

      RETURN

      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE FC_SUBDIVIDE(fC)
         IMPLICIT NONE

         TYPE(fCellType), INTENT(INOUT) :: fC

         INTEGER :: i, a

         INTEGER, ALLOCATABLE :: iOrd(:,:)
         REAL(KIND=8), ALLOCATABLE :: xl(:,:), xi(:,:)

         IF (ALLOCATED(fC%xiGP)) DEALLOCATE(fC%xiGP)

         SELECT CASE (fC%eType)
         CASE (eType_BIL)
            fC%nSub = 4
            ALLOCATE(fC%sub(fC%nSub), iOrd(fC%eNoN,fC%nSub), xl(nsd,9),
     2         xi(nsd,9))

            iOrd(:,1) = (/1,5,9,8/)
            iOrd(:,2) = (/5,2,6,9/)
            iOrd(:,3) = (/9,6,3,7/)
            iOrd(:,4) = (/8,9,7,4/)

            xl(:,1:4) = fC%x(:,:)
            xl(:,5)   = (xl(:,1) + xl(:,2))/2D0
            xl(:,6)   = (xl(:,2) + xl(:,3))/2D0
            xl(:,7)   = (xl(:,3) + xl(:,4))/2D0
            xl(:,8)   = (xl(:,4) + xl(:,1))/2D0
            xl(:,9)   = SUM(fC%x, DIM=2)/4D0

            xi(:,1:4) = fC%xi(:,:)
            xi(:,5)   = (xi(:,1) + xi(:,2))/2D0
            xi(:,6)   = (xi(:,2) + xi(:,3))/2D0
            xi(:,7)   = (xi(:,3) + xi(:,4))/2D0
            xi(:,8)   = (xi(:,4) + xi(:,1))/2D0
            xi(:,9)   = SUM(fC%xi, DIM=2)/4D0

            DO i=1, fC%nSub
               fC%sub(i)%eType = fC%eType
               fC%sub(i)%ilev  = fC%ilev + 1
               fC%sub(i)%eNoN  = fC%eNoN
               fC%sub(i)%nG    = fC%nG
               ALLOCATE(fC%sub(i)%incG(fC%sub(i)%nG),
     2            fC%sub(i)%x(nsd,fC%sub(i)%eNoN),
     3            fC%sub(i)%xi(nsd,fC%sub(i)%eNoN),
     4            fC%sub(i)%w(fC%sub(i)%nG))
               fC%sub(i)%incG(:) = .TRUE.
               fC%sub(i)%w(:) = fC%w(:)/4D0
               DO a=1, fC%eNoN
                  fC%sub(i)%x(:,a)  = xl(:,iOrd(a,i))
                  fC%sub(i)%xi(:,a) = xi(:,iOrd(a,i))
               END DO
            END DO

         CASE DEFAULT
            err = "Unknown element type in FC_SUBDIVIDE"

         END SELECT

         RETURN
         END SUBROUTINE FC_SUBDIVIDE
!--------------------------------------------------------------------
      END SUBROUTINE FC_SET
!####################################################################
      RECURSIVE SUBROUTINE FC_vInteg(lFc, xe, f, fI)
      USE COMMOD
      IMPLICIT NONE

      TYPE(fCellType), INTENT(IN) :: lFc
      REAL(KIND=8), INTENT(IN) :: xe(nsd,lFc%eNoN), f(lFc%eNoN)
      REAL(KIND=8), INTENT(INOUT) :: fI

      INTEGER :: g, i, a, eNoN
      REAL(KIND=8) :: Jac, fHat, tmp(nsd,nsd)

      REAL(KIND=8), ALLOCATABLE :: N(:), Nxi(:,:), Nx(:,:)

      IF (lFc%nSub .EQ. 0) THEN
         eNoN = lFc%eNoN
         ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN))
         DO g=1, lFc%nG
            IF (.NOT.lFc%incG(g)) CYCLE
            CALL GETGNN(nsd, lFc%eType, eNoN, lFc%xiGP(:,g), N, Nxi)
            CALL GNN(eNoN, nsd, Nxi, xe, Nx, Jac, tmp)
            fHat = 0D0
            DO a=1, eNoN
               fHat = fHat + f(a)*N(a)
            END DO
            fI = fI + lFc%w(g)*Jac*fHat
         END DO
         DEALLOCATE(N, Nxi, Nx)
      ELSE
         DO i=1, lFc%nSub
            CALL FC_vInteg(lFc%sub(i), xe, f, fI)
         END DO
      END IF

      RETURN
      END SUBROUTINE FC_vInteg
!####################################################################
      RECURSIVE SUBROUTINE FC_CONSTRUCT(lFc, ae, ye, xe, lR, lK)
      USE COMMOD
      IMPLICIT NONE

      TYPE(fCellType), INTENT(IN) :: lFc
      REAL(KIND=8), INTENT(IN) :: ae(tDof,lFc%eNoN), ye(tDof,lFc%eNoN),
     2   xe(nsd,lFc%eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lR(dof,lFc%eNoN)
      REAL(KIND=8), INTENT(INOUT) :: lK(dof*dof,lFc%eNoN,lFc%eNoN)

      INTEGER :: i, g, eNoN
      REAL(KIND=8) :: w, Jac, ksix(nsd,nsd)

      REAL(KIND=8), ALLOCATABLE :: N(:), Nxi(:,:), Nx(:,:)

      IF (lFc%nSub .EQ. 0) THEN
         eNoN = lFc%eNoN
         ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN))
         DO g=1, lFc%nG
            IF (.NOT. lFc%incG(g)) CYCLE
            CALL GETGNN(nsd, lFc%eType, eNoN, lFc%xiGP(:,g), N, Nxi)
            CALL GNN(eNoN, nsd, Nxi, xe, Nx, Jac, ksix)
            w = lFc%w(g)*Jac

            IF (nsd .EQ. 3) THEN
               CALL FLUID3D(eNoN, w, N, Nx, ae, ye, ksix, lR, lK)
            ELSE
               CALL FLUID2D(eNoN, w, N, Nx, ae, ye, ksix, lR, lK)
            END IF
         END DO
         DEALLOCATE(N, Nxi, Nx)
      ELSE
         DO i=1, lFc%nSub
            CALL FC_CONSTRUCT(lFc%sub(i), ae, ye, xe, lR, lK)
         END DO
      END IF

      RETURN
      END SUBROUTINE FC_CONSTRUCT
!####################################################################

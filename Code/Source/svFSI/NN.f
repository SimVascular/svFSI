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
!     This subroutine is mainly intended for calculation of shape
!     function integration. Also, derivative of local coordinate with
!     respect to the parent elements is calculated here.
!
!--------------------------------------------------------------------

!     Here parameters related to the elemnt are set
      SUBROUTINE SELECTELE(lM)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER :: insd, g

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      IF (lM%lFib) insd = 1

      IF (lM%eType .EQ. eType_NRB) THEN
         ALLOCATE(lM%w(lM%nG), lM%N(lM%eNoN,lM%nG),
     2      lM%Nx(insd,lM%eNoN,lM%nG))
         IF (insd .EQ. 2) THEN
            ALLOCATE(lM%Nxx(3,lM%eNoN,lM%nG))
         ELSE IF (insd .EQ. 3) THEN
            ALLOCATE(lM%Nxx(6,lM%eNoN,lM%nG))
         END IF
         RETURN
      END IF

      IF (insd .EQ. 3) THEN
         SELECT CASE (lM%eNoN)
         CASE(8)
            lM%eType   = eType_BRK
            lM%nG      = 8
            lM%vtkType = 12
            lM%nEf     = 6
            lM%lShpF   = .FALSE.
         CASE(6)
            lM%eType   = eType_WDG
            lM%nG      = 6
            lM%vtkType = 13
            lM%nEf     = 3
            lM%lShpF   = .FALSE.
         CASE(4)
            lM%eType   = eType_TET
            lM%nG      = 4
            lM%vtkType = 10
            lM%nEf     = 4
            lM%lShpF   = .TRUE.
         CASE DEFAULT
            err = "Unable to identify combination of nsd and eNoN"
         END SELECT

      ELSE IF (insd .EQ. 2) THEN
         SELECT CASE (lM%eNoN)
         CASE(3)
            lM%eType   = eType_TRI
            lM%nG      = 3
            lM%vtkType = 5
            lM%nEf     = 3
            lM%lShpF   = .TRUE.
         CASE(4)
            lM%eType   = eType_BIL
            lM%nG      = 4
            lM%vtkType = 9
            lM%nEf     = 4
            lM%lShpF   = .FALSE.
         CASE(9)
            lM%eType   = eType_BIQ
            lM%nG      = 9
            lM%vtkType = 28
            lM%nEf     = 4
            lM%lShpF   = .FALSE.
         CASE DEFAULT
            err = "Unable to identify combination of nsd and eNoN"
         END SELECT

      ELSE IF (insd .EQ. 1) THEN
         SELECT CASE (lM%eNoN)
         CASE(2)
            lM%eType   = eType_LIN
            lM%nG      = 2
            lM%vtkType = 3
            lM%nEf     = 2
            lM%lShpF   = .TRUE.
         CASE(3)
            lM%eType   = eType_QUD
            lM%nG      = 3
            lM%vtkType = 3
            lM%nEf     = 2
            lM%lShpF   = .FALSE.
         CASE DEFAULT
            err = "Unable to identify combination of nsd and eNoN"
         END SELECT
      END IF

      ALLOCATE(lM%w(lM%nG), lM%xi(insd,lM%nG), lM%N(lM%eNoN,lM%nG),
     2   lM%Nx(insd,lM%eNoN,lM%nG))

      CALL GETGIP(insd, lM%eType, lM%nG, lM%w, lM%xi)

      DO g=1, lM%nG
         CALL GETGNN(insd, lM%eType, lM%eNoN, lM%xi(:,g), lM%N(:,g),
     2      lM%Nx(:,:,g))
      END DO

      RETURN
      END SUBROUTINE SELECTELE
!--------------------------------------------------------------------
!     This routine selects boundary element type
      SUBROUTINE SELECTELEB(lM, lFa)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER :: insd, g

      insd = nsd - 1
      IF (lM%lShl) insd = insd - 1
      IF (lM%lFib) insd = 0

      IF (lM%eType .EQ. eType_NRB) THEN
         lFa%eType = eType_NRB
         lFa%nG    = lM%nG/lM%bs(lFa%d)%nG
         lFa%eNoN  = lM%eNoN/(lM%bs(lFa%d)%p + 1)

         ALLOCATE(lFa%w(lFa%nG), lFa%N(lFa%eNoN,lFa%nG),
     2      lFa%Nx(insd,lFa%eNoN,lFa%nG))
         IF (insd .EQ. 1) THEN
            ALLOCATE(lFa%Nxx(1,lFa%eNoN,lFa%nG))
         ELSE IF (insd .EQ. 2) THEN
            ALLOCATE(lFa%Nxx(3,lFa%eNoN,lFa%nG))
         END IF
         RETURN
      END IF

      IF (insd .EQ. 2) THEN
         SELECT CASE (lFa%eNoN)
         CASE(4)
            lFa%eType = eType_BIL
            lFa%nG    = 4
         CASE(3)
            lFa%eType = eType_TRI
            lFa%nG    = 3
         CASE DEFAULT
            err = "Unable to identify combination of nsd and eNoN"
         END SELECT
      ELSE IF (insd .EQ. 1) THEN
         SELECT CASE (lFa%eNoN)
         CASE(2)
            lFa%eType = eType_LIN
            lFa%nG    = 2
         CASE(3)
            lFa%eType = eType_QUD
            lFa%nG    = 3
         CASE DEFAULT
            err = "Unable to identify combination of nsd and eNoN"
         END SELECT
      ELSE IF (insd .EQ. 0) THEN
         lFa%eType = eType_PNT
         lFa%nG = 1
      END IF

      ALLOCATE(lFa%w(lFa%nG), lFa%xi(insd,lFa%nG),
     2   lFa%N(lFa%eNoN,lFa%nG), lFa%Nx(insd,lFa%eNoN,lFa%nG))

      CALL GETGIP(insd, lFa%eType, lFa%nG, lFa%w, lFa%xi)

      DO g=1, lFa%nG
         CALL GETGNN(insd, lFa%eType, lFa%eNoN, lFa%xi(:,g),
     2      lFa%N(:,g), lFa%Nx(:,:,g))
      END DO

      RETURN
      END SUBROUTINE SELECTELEB
!####################################################################
!     Returns Gauss integration points in local (ref) coordinates
      PURE SUBROUTINE GETGIP(insd, eType, nG, w, xi)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: insd, eType, nG
      REAL(KIND=8), INTENT(OUT) :: w(nG), xi(insd,nG)

      REAL(KIND=8) s, t, lz, uz

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_BRK)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t
      CASE(eType_TET)
         w = 1D0/24D0
         s = (5D0 + 3D0*SQRT(5D0))/2D1
         t = (5D0 -     SQRT(5D0))/2D1
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
      CASE(eType_WDG)
         w  =  1D0/6D0
         s  =  2D0/3D0
         t  =  1D0/6D0
         uz =  1D0/SQRT(3D0)
         lz = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = lz
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = lz
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = lz
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = uz
         xi(1,5) = t; xi(2,5) = s; xi(3,5) = uz
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = uz

!     2D elements
      CASE(eType_TRI)
         w = 1D0/6D0
         s = 2D0/3D0
         t = 1D0/6D0
         xi(1,1) = t; xi(2,1) = t
         xi(1,2) = s; xi(2,2) = t
         xi(1,3) = t; xi(2,3) = s
      CASE(eType_BIL)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
      CASE(eType_BIQ)
         w(1) = 25D0/81D0; w(2) = 25D0/81D0; w(3) = 25D0/81D0
         w(4) = 25D0/81D0; w(5) = 40D0/81D0; w(6) = 40D0/81D0
         w(7) = 40D0/81D0; w(8) = 40D0/81D0; w(9) = 64D0/81D0
         s    = SQRT(6D-1)
         xi(1,1) =  -s; xi(2,1) =  -s
         xi(1,2) =   s; xi(2,2) =  -s
         xi(1,3) =   s; xi(2,3) =   s
         xi(1,4) =  -s; xi(2,4) =   s
         xi(1,5) = 0D0; xi(2,5) =  -s
         xi(1,6) =   s; xi(2,6) = 0D0
         xi(1,7) = 0D0; xi(2,7) =   s
         xi(1,8) =  -s; xi(2,8) = 0D0
         xi(1,9) = 0D0; xi(2,9) = 0D0

!     1D elements
      CASE(eType_LIN)
         w = 1D0
         s = 1D0/SQRT(3D0)
         xi(1,1) = -s
         xi(1,2) =  s
      CASE(eType_QUD)
         w(1) = 5D0/9D0; w(2) = 5D0/9D0; w(3) = 8D0/9D0
         s = SQRT(6D-1)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0D0
!     0D elements
      CASE(eType_PNT)
         w = 1D0
      END SELECT

      END SUBROUTINE GETGIP
!####################################################################
!     Inverse maps {xp} to {$\xi$} in an element with coordinates {xl}
!     using Newton's method
      SUBROUTINE GETXI(eType, eNoN, xl, xp, xi, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eType, eNoN
      REAL(KIND=8), INTENT(IN)  :: xl(nsd,eNoN), xp(nsd)
      REAL(KIND=8), INTENT(INOUT) :: xi(nsd)
      LOGICAL, INTENT(OUT) :: flag

      INTEGER, PARAMETER :: MAXITR = 5
      REAL(KIND=8), PARAMETER :: RTOL = 1D-6, ATOL = 1D-12

      LOGICAL :: l1, l2, l3
      INTEGER :: itr, i, j, a
      REAL(KIND=8) :: rmsA, rmsR, N(eNoN), Nxi(nsd,eNoN), xiK(nsd),
     2   xK(nsd), rK(nsd), Am(nsd,nsd)

      itr = 0
      xiK = xi
c      WRITE(1000+cm%tF(),'(8X,A)') "Newton iterations.."
c      WRITE(1000+cm%tF(),'(10X,A)',ADVANCE='NO') "Initial xi: "
c      DO i=1, nsd
c         WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xi(i))
c      END DO
c      WRITE(1000+cm%tF(),'(A)')
      DO
         itr = itr + 1

         CALL GETGNN(nsd, eType, eNoN, xiK, N, Nxi)
         xK = 0D0
         DO i=1, nsd
            DO a=1, eNoN
               xK(i) = xK(i) + N(a)*xl(i,a)
            END DO
            rK(i) = xK(i) - xp(i)
         END DO

         rmsA = 0.0D0
         rmsR = 0.0D0
         DO i=1, nsd
            rmsA = rmsA + rK(i)**2.0D0
            rmsR = rmsR + (rK(i) / (xK(i)+eps))**2.0D0
         END DO
         rmsA = SQRT(rmsA/REAL(nsd,KIND=8))
         rmsR = SQRT(rmsR/REAL(nsd,KIND=8))

c         WRITE(1000+cm%tF(),'(10X,A)') "Iter: "//STR(itr)//" "//
c     2      STR(rmsA)//" "//STR(rmsR)

         l1 = itr .GT. MAXITR
         l2 = rmsA .LE. ATOL
         l3 = rmsR .LE. RTOL
         IF (l1 .OR. l2 .OR. l3) EXIT

         Am = 0.0D0
         DO i=1, nsd
            DO j=1, nsd
               DO a=1, eNoN
                  Am(i,j) = Am(i,j) + xl(i,a)*Nxi(j,a)
               END DO
            END DO
         END DO
         Am  = MAT_INV(Am, nsd)
         rK  = MATMUL(Am, rK)
         xiK = xiK - rK
      END DO

      IF (l2 .OR. l3) THEN
!     Newton's method converges
         flag = .TRUE.
c         WRITE(1000+cm%tF(),'(10X,A)') "Success.."
      ELSE
!     Newton's method failed to converge
         flag = .FALSE.
c         WRITE(1000+cm%tF(),'(10X,A)') "Fail.."
      END IF

      xi(:) = xiK(:)

      RETURN
      END SUBROUTINE GETXI
!####################################################################
!     Returns shape functions and derivatives at given natural coords
      PURE SUBROUTINE GETGNN(insd, eType, eNoN, xi, N, Nxi)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: insd, eType, eNoN
      REAL(KIND=8), INTENT(IN)  :: xi(insd)
      REAL(KIND=8), INTENT(OUT) :: N(eNoN), Nxi(insd,eNoN)

      REAL(KIND=8) :: s, t, mx, my, ux, uy, uz, lx, ly, lz

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_BRK)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         uz = 1D0 + xi(3); lz = 1D0 - xi(3)
         N(1) = ux*uy*uz/8D0
         N(2) = lx*uy*uz/8D0
         N(3) = lx*uy*lz/8D0
         N(4) = ux*uy*lz/8D0
         N(5) = ux*ly*uz/8D0
         N(6) = lx*ly*uz/8D0
         N(7) = lx*ly*lz/8D0
         N(8) = ux*ly*lz/8D0

         Nxi(1,1) =  uy*uz/8D0
         Nxi(2,1) =  ux*uz/8D0
         Nxi(3,1) =  ux*uy/8D0
         Nxi(1,2) = -uy*uz/8D0
         Nxi(2,2) =  lx*uz/8D0
         Nxi(3,2) =  lx*uy/8D0
         Nxi(1,3) = -uy*lz/8D0
         Nxi(2,3) =  lx*lz/8D0
         Nxi(3,3) = -lx*uy/8D0
         Nxi(1,4) =  uy*lz/8D0
         Nxi(2,4) =  ux*lz/8D0
         Nxi(3,4) = -ux*uy/8D0
         Nxi(1,5) =  ly*uz/8D0
         Nxi(2,5) = -ux*uz/8D0
         Nxi(3,5) =  ux*ly/8D0
         Nxi(1,6) = -ly*uz/8D0
         Nxi(2,6) = -lx*uz/8D0
         Nxi(3,6) =  lx*ly/8D0
         Nxi(1,7) = -ly*lz/8D0
         Nxi(2,7) = -lx*lz/8D0
         Nxi(3,7) = -lx*ly/8D0
         Nxi(1,8) =  ly*lz/8D0
         Nxi(2,8) = -ux*lz/8D0
         Nxi(3,8) = -ux*ly/8D0
      CASE(eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1D0 - xi(1) - xi(2) - xi(3)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(3,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(3,2) =  0D0
         Nxi(1,3) =  0D0
         Nxi(2,3) =  0D0
         Nxi(3,3) =  1D0
         Nxi(1,4) = -1D0
         Nxi(2,4) = -1D0
         Nxi(3,4) = -1D0
      CASE(eType_WDG)
         ux = xi(1) ; uy = xi(2) ; uz = 1D0 - ux - uy
         s = (1D0 + xi(3))/2D0; t = (1D0 - xi(3))/2D0
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

         Nxi(1,1) =  t
         Nxi(2,1) =  0D0
         Nxi(3,1) = -ux/2D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  t
         Nxi(3,2) = -uy/2D0
         Nxi(1,3) = -t
         Nxi(2,3) = -t
         Nxi(3,3) = -uz/2D0
         Nxi(1,4) =  s
         Nxi(2,4) =  0D0
         Nxi(3,4) =  ux/2D0
         Nxi(1,5) =  0D0
         Nxi(2,5) =  s
         Nxi(3,5) =  uy/2D0
         Nxi(1,6) = -s
         Nxi(2,6) = -s
         Nxi(3,6) =  uz/2D0

!     2D elements
      CASE(eType_TRI)
         N(1) = 1D0 - xi(1) - xi(2)
         N(2) = xi(1)
         N(3) = xi(2)

         Nxi(1,1) = -1D0
         Nxi(2,1) = -1D0
         Nxi(1,2) =  1D0
         Nxi(2,2) =  0D0
         Nxi(1,3) =  0D0
         Nxi(2,3) =  1D0
      CASE(eType_BIL)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0

         Nxi(1,1) =  uy/4D0
         Nxi(2,1) =  ux/4D0
         Nxi(1,2) = -uy/4D0
         Nxi(2,2) =  lx/4D0
         Nxi(1,3) = -ly/4D0
         Nxi(2,3) = -lx/4D0
         Nxi(1,4) =  ly/4D0
         Nxi(2,4) = -ux/4D0
      CASE(eType_BIQ)
         ux = 1D0 + xi(1); mx = xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); my = xi(2); ly = 1D0 - xi(2)
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy

         Nxi(1,1) =  (lx - mx)*my*ly/4D0
         Nxi(2,1) =  (ly - my)*mx*lx/4D0
         Nxi(1,2) = -(ux + mx)*my*ly/4D0
         Nxi(2,2) = -(ly - my)*mx*ux/4D0
         Nxi(1,3) =  (ux + mx)*my*uy/4D0
         Nxi(2,3) =  (uy + my)*mx*ux/4D0
         Nxi(1,4) = -(lx - mx)*my*uy/4D0
         Nxi(2,4) = -(uy + my)*mx*lx/4D0
         Nxi(1,5) = -(lx - ux)*my*ly/2D0
         Nxi(2,5) = -(ly - my)*lx*ux/2D0
         Nxi(1,6) =  (ux + mx)*ly*uy/2D0
         Nxi(2,6) =  (ly - uy)*mx*ux/2D0
         Nxi(1,7) =  (lx - ux)*my*uy/2D0
         Nxi(2,7) =  (uy + my)*lx*ux/2D0
         Nxi(1,8) = -(lx - mx)*ly*uy/2D0
         Nxi(2,8) = -(ly - uy)*mx*lx/2D0
         Nxi(1,9) =  (lx - ux)*ly*uy
         Nxi(2,9) =  (ly - uy)*lx*ux

!     1D elements
      CASE(eType_LIN)
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0

         Nxi(1,1) = -5D-1
         Nxi(1,2) =  5D-1
      CASE(eType_QUD)
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))

         Nxi(1,1) = -5D-1 + xi(1)
         Nxi(1,2) =  5D-1 + xi(1)
         Nxi(1,3) = -2D0*xi(1)

!     0D elements
      CASE(eType_PNT)
         N(1) = 1.0D0
      END SELECT

      RETURN
      END SUBROUTINE GETGNN
!####################################################################
      PURE SUBROUTINE GNN(eNoN, insd, Nxi, x, Nx, Jac, ks)
      USE COMMOD, ONLY: nsd
      USE UTILMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN, insd
      REAL(KIND=8), INTENT(IN) :: Nxi(insd,eNoN), x(nsd,eNoN)
      REAL(KIND=8), INTENT(OUT) :: Nx(insd,eNoN), Jac, ks(nsd,nsd)

      INTEGER a
      REAL(KIND=8) xXi(nsd,insd), xiX(insd,nsd)

      xXi = 0D0
      Jac = 0D0
      Nx  = 0D0
      ks  = 0D0
      IF (insd .EQ. 1) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
         END DO

         Jac = SQRT(NORM(xXi)) + 1D3*eps
         DO a=1, eNoN
            Nx(1,a) = Nxi(1,a)/Jac
         END DO
      ELSE IF (insd .EQ. 2) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
         ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
         ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
         ks(2,1) = ks(1,2)

         DO a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         END DO
      ELSE IF (insd .EQ. 3) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         END DO

         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3)
     2       + xXi(1,2)*xXi(2,3)*xXi(3,1)
     3       + xXi(1,3)*xXi(2,1)*xXi(3,2)
     4       - xXi(1,1)*xXi(2,3)*xXi(3,2)
     5       - xXi(1,2)*xXi(2,1)*xXi(3,3)
     6       - xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1)+xiX(3,1)*xiX(3,1)
         ks(1,2) = xiX(1,2)*xiX(1,1)+xiX(2,2)*xiX(2,1)+xiX(3,2)*xiX(3,1)
         ks(1,3) = xiX(1,3)*xiX(1,1)+xiX(2,3)*xiX(2,1)+xiX(3,3)*xiX(3,1)
         ks(2,2) = xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2)+xiX(3,2)*xiX(3,2)
         ks(2,3) = xiX(1,2)*xiX(1,3)+xiX(2,2)*xiX(2,3)+xiX(3,2)*xiX(3,3)
         ks(3,3) = xiX(1,3)*xiX(1,3)+xiX(2,3)*xiX(2,3)+xiX(3,3)*xiX(3,3)
         ks(2,1) = ks(1,2)
         ks(3,1) = ks(1,3)
         ks(3,2) = ks(2,3)

         DO a=1, eNoN
            Nx(1,a) = Nx(1,a) + Nxi(1,a)*xiX(1,1)
     2                        + Nxi(2,a)*xiX(2,1)
     3                        + Nxi(3,a)*xiX(3,1)

            Nx(2,a) = Nx(2,a) + Nxi(1,a)*xiX(1,2)
     2                        + Nxi(2,a)*xiX(2,2)
     3                        + Nxi(3,a)*xiX(3,2)

            Nx(3,a) = Nx(3,a) + Nxi(1,a)*xiX(1,3)
     2                        + Nxi(2,a)*xiX(2,3)
     3                        + Nxi(3,a)*xiX(3,3)
         END DO
      END IF

      RETURN
      END SUBROUTINE GNN
!--------------------------------------------------------------------
!     Compute shell kinematics: normal vector, covariant & contravariant
!     basis vectors
      SUBROUTINE GNNS(eNoN, Nxi, xl, nV, gCov, gCnv)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: Nxi(nsd-1,eNoN), xl(nsd,eNoN)
      REAL(KIND=8), INTENT(OUT) :: nV(nsd), gCov(nsd,nsd-1),
     2   gCnv(nsd,nsd-1)

      INTEGER a, i, j, insd
      REAL(KIND=8), ALLOCATABLE :: xXi(:,:), Gmat(:,:)

      insd = nsd - 1
      ALLOCATE(xXi(nsd,insd), Gmat(insd,insd))

!     Calculating surface deflation
      xXi = 0D0
      DO a=1, eNoN
         DO i=1, insd
            xXi(:,i) = xXi(:,i) + xl(:,a)*Nxi(i,a)
         END DO
      END DO
      nV = CROSS(xXi)

!     Covariant basis
      gCov = xXi

!     Metric tensor g_i . g_j
      Gmat = 0D0
      DO i=1, insd
         DO j=1, insd
            DO a=1, nsd
               Gmat(i,j) = Gmat(i,j) + gCov(a,i)*gCov(a,j)
            END DO
         END DO
      END DO

!     Contravariant basis
      Gmat = MAT_INV(Gmat, insd)
      gCnv = 0D0
      DO i=1, insd
         DO j=1, insd
            gCnv(:,i) = gCnv(:,i) + Gmat(i,j)*gCov(:,j)
         END DO
      END DO

      RETURN
      END SUBROUTINE GNNS
!--------------------------------------------------------------------
!     This routine returns a vector at element "e" and Gauss point
!     "g" of face "lFa" that is the normal weigthed by Jac, i.e.
!     Jac = SQRT(NORM(n)).
      SUBROUTINE GNNIB(lFa, e, g, n)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e, g
      REAL(KIND=8), INTENT(OUT) :: n(nsd)
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER a, Ac, i, iM, Ec, b, Bc, eNoN, insd
      REAL(KIND=8) v(nsd)

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: lX(:,:), xXi(:,:)

      iM   = lFa%iM
      Ec   = lFa%gE(e)
      eNoN = ib%msh(iM)%eNoN
      insd = nsd - 1

      ALLOCATE(lX(nsd,eNoN), ptr(eNoN), setIt(eNoN))
!     Creating a ptr list that contains pointer to the nodes of elements
!     that are at the face at the beginning of the list and the rest at
!     the end
      setIt = .TRUE.
      DO a=1, lFa%eNoN
         Ac = lFa%IEN(a,e)
         DO b=1, eNoN
            IF (setIt(b)) THEN
               Bc = ib%msh(iM)%IEN(b,Ec)
               IF (Bc .EQ. Ac) EXIT
            END IF
         END DO
         ptr(a)   = b
         setIt(b) = .FALSE.
      END DO
      a = lFa%eNoN
      DO b=1, eNoN
         IF (setIt(b)) THEN
            a      = a + 1
            ptr(a) = b
         END IF
      END DO

!     Correct the position vector
      DO a=1, eNoN
         Ac = ib%msh(iM)%IEN(a,Ec)
         lX(:,a) = ib%x(:,Ac) + ib%Uo(:,Ac)
      END DO

!     Calculating surface deflation
      IF (ib%msh(iM)%lShl) THEN
!        Since the face has only one parametric coordinate (edge), find
!        its normal from cross product of mesh normal and interior edge

!        Update shape functions if NURBS
         IF (ib%msh(iM)%eType .EQ. eType_NRB)
     2      CALL NRBNNX(ib%msh(iM), Ec)

!        Compute adjoining mesh element normal
         ALLOCATE(xXi(nsd,insd))
         xXi = 0D0
         DO a=1, eNoN
            DO i=1, insd
               xXi(:,i) = xXi(:,i) + lX(:,a)*ib%msh(iM)%Nx(i,a,g)
            END DO
         END DO
         v(:) = CROSS(xXi)
         v(:) = v(:) / SQRT(NORM(v))
         DEALLOCATE(xXi)

!        Face element surface deflation
         ALLOCATE(xXi(nsd,1))
         xXi = 0D0
         DO a=1, lFa%eNoN
            b = ptr(a)
            xXi(:,1) = xXi(:,1) + lFa%Nx(1,a,g)*lX(:,b)
         END DO

!        Face normal
         n(1) = v(2)*xXi(3,1) - v(3)*xXi(2,1)
         n(2) = v(3)*xXi(1,1) - v(1)*xXi(3,1)
         n(3) = v(1)*xXi(2,1) - v(2)*xXi(1,1)

!        I choose Gauss point of the mesh element for calculating
!        interior edge
         v(:) = 0D0
         DO a=1, eNoN
            v(:) = v(:) + lX(:,a)*ib%msh(iM)%N(a,g)
         END DO
         a = ptr(1)
         v(:) = lX(:,a) - v(:)
         IF (NORM(n,v) .LT. 0D0) n = -n

         DEALLOCATE(xXi)
         RETURN
      ELSE
         ALLOCATE(xXi(nsd,insd))
         xXi = 0D0
         DO a=1, lFa%eNoN
            b = ptr(a)
            DO i=1, insd
               xXi(:,i) = xXi(:,i) + lFa%Nx(i,a,g)*lX(:,b)
            END DO
         END DO
         n = CROSS(xXi)
         DEALLOCATE(xXi)
      END IF

!     Changing the sign if neccessary. a locates on the face and b
!     outside of the face, in the parent element
      a = ptr(1)
      b = ptr(lFa%eNoN+1)
      v = lX(:,a) - lX(:,b)
      IF (NORM(n,v) .LT. 0D0) n = -n

      RETURN
      END SUBROUTINE GNNIB
!--------------------------------------------------------------------
!     This routine returns a vector at element "e" and Gauss point
!     "g" of face "lFa" that is the normal weigthed by Jac, i.e.
!     Jac = SQRT(NORM(n)).
      SUBROUTINE GNNB(lFa, e, g, n)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e, g
      REAL(KIND=8), INTENT(OUT) :: n(nsd)
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER a, Ac, i, iM, Ec, b, Bc, eNoN, insd
      REAL(KIND=8) v(nsd)

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: lX(:,:), xXi(:,:)

      iM   = lFa%iM
      Ec   = lFa%gE(e)
      eNoN = msh(iM)%eNoN
      insd = nsd - 1

      ALLOCATE(lX(nsd,eNoN), ptr(eNoN), setIt(eNoN))

!     Creating a ptr list that contains pointer to the nodes of elements
!     that are at the face at the beginning of the list and the rest at
!     the end
      setIt = .TRUE.
      DO a=1, lFa%eNoN
         Ac = lFa%IEN(a,e)
         DO b=1, eNoN
            IF (setIt(b)) THEN
               Bc = msh(iM)%IEN(b,Ec)
               IF (Bc .EQ. Ac) EXIT
            END IF
         END DO
         ptr(a)   = b
         setIt(b) = .FALSE.
      END DO
      a = lFa%eNoN
      DO b=1, eNoN
         IF (setIt(b)) THEN
            a      = a + 1
            ptr(a) = b
         END IF
      END DO

!     Correcting the position vector if mesh is moving
      DO a=1, eNoN
         Ac = msh(iM)%IEN(a,Ec)
         lX(:,a) = x(:,Ac)
         IF (mvMsh) lX(:,a) = lX(:,a) + Do(nsd+2:2*nsd+1,Ac)
      END DO

!     Calculating surface deflation
      IF (msh(iM)%lShl) THEN
!        Since the face has only one parametric coordinate (edge), find
!        its normal from cross product of mesh normal and interior edge

!        Update shape functions if NURBS
         IF (msh(iM)%eType .EQ. eType_NRB) CALL NRBNNX(msh(iM), Ec)

!        Compute adjoining mesh element normal
         ALLOCATE(xXi(nsd,insd))
         xXi = 0D0
         DO a=1, eNoN
            DO i=1, insd
               xXi(:,i) = xXi(:,i) + lX(:,a)*msh(iM)%Nx(i,a,g)
            END DO
         END DO
         v(:) = CROSS(xXi)
         v(:) = v(:) / SQRT(NORM(v))
         DEALLOCATE(xXi)

!        Face element surface deflation
         ALLOCATE(xXi(nsd,1))
         xXi = 0D0
         DO a=1, lFa%eNoN
            b = ptr(a)
            xXi(:,1) = xXi(:,1) + lFa%Nx(1,a,g)*lX(:,b)
         END DO

!        Face normal
         n(1) = v(2)*xXi(3,1) - v(3)*xXi(2,1)
         n(2) = v(3)*xXi(1,1) - v(1)*xXi(3,1)
         n(3) = v(1)*xXi(2,1) - v(2)*xXi(1,1)

!        I choose Gauss point of the mesh element for calculating
!        interior edge
         v(:) = 0D0
         DO a=1, eNoN
            v(:) = v(:) + lX(:,a)*msh(iM)%N(a,g)
         END DO
         a = ptr(1)
         v(:) = lX(:,a) - v(:)
         IF (NORM(n,v) .LT. 0D0) n = -n

         DEALLOCATE(xXi)
         RETURN
      ELSE
         ALLOCATE(xXi(nsd,insd))
         xXi = 0D0
         DO a=1, lFa%eNoN
            b = ptr(a)
            DO i=1, insd
               xXi(:,i) = xXi(:,i) + lFa%Nx(i,a,g)*lX(:,b)
            END DO
         END DO
         n = CROSS(xXi)
         DEALLOCATE(xXi)
      END IF

!     Changing the sign if neccessary. a locates on the face and b
!     outside of the face, in the parent element
      a = ptr(1)
      b = ptr(lFa%eNoN+1)
      v = lX(:,a) - lX(:,b)
      IF (NORM(n,v) .LT. 0D0) n = -n

      RETURN
      END SUBROUTINE GNNB
!####################################################################


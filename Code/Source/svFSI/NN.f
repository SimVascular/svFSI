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

      INTEGER(KIND=IKIND) :: insd, g

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
         CASE(10)
            lM%eType   = eType_QTE
            lM%nG      = 15
            lM%vtkType = 24
            lM%nEf     = 4
            lM%lShpF   = .FALSE.
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
         CASE(6)
            lM%eType   = eType_QTR
            lM%nG      = 7
            lM%vtkType = 22
            lM%nEf     = 3
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

      ALLOCATE(lM%xib(2,nsd), lM%Nb(2,lM%eNoN))
      CALL GETNNBNDS(lM%eType, lM%eNoN, lM%xib, lM%Nb)

      RETURN
      END SUBROUTINE SELECTELE
!--------------------------------------------------------------------
!     This routine selects boundary element type
      SUBROUTINE SELECTELEB(lM, lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) :: insd, g

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
         CASE(6)
            lFa%eType = eType_QTR
            lFa%nG    = 7
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
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, nG
      REAL(KIND=RKIND), INTENT(OUT) :: w(nG), xi(insd,nG)

      REAL(KIND=RKIND) s, t, lz, uz

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_BRK)
         w = 1._RKIND
         s =  1._RKIND/SQRT(3._RKIND)
         t = -1._RKIND/SQRT(3._RKIND)
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t
      CASE(eType_TET)
         w = 1._RKIND/24._RKIND
         s = (5._RKIND + 3._RKIND*SQRT(5._RKIND))/20._RKIND
         t = (5._RKIND -          SQRT(5._RKIND))/20._RKIND
         xi(1,1) = t; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = s; xi(2,2) = t; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = s
      CASE(eType_WDG)
         w  =  1._RKIND/6._RKIND
         s  =  2._RKIND/3._RKIND
         t  =  1._RKIND/6._RKIND
         uz =  1._RKIND/SQRT(3._RKIND)
         lz = -1._RKIND/SQRT(3._RKIND)
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = lz
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = lz
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = lz
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = uz
         xi(1,5) = t; xi(2,5) = s; xi(3,5) = uz
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = uz
      CASE(eType_QTE)
         w(1)     = 0.030283678097089_RKIND
         w(2:5)   = 0.006026785714286_RKIND
         w(6:9)   = 0.011645249086029_RKIND
         w(10:15) = 0.010949141561386_RKIND

         s = 0.25_RKIND
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s

         s = 0.333333333333333_RKIND
         t = 0._RKIND
         xi(1,2) = s; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = s
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = s
         xi(1,5) = s; xi(2,5) = s; xi(3,5) = t

         s = 0.090909090909091_RKIND
         t = 0.727272727272727_RKIND
         xi(1,6) = s; xi(2,6) = s; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = s; xi(3,7) = s
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = s
         xi(1,9) = s; xi(2,9) = s; xi(3,9) = t

         s = 0.066550153573664D0
         t = 0.433449846426336D0
         xi(1,10) = t; xi(2,10) = s; xi(3,10) = s
         xi(1,11) = s; xi(2,11) = t; xi(3,11) = s
         xi(1,12) = s; xi(2,12) = s; xi(3,12) = t
         xi(1,13) = s; xi(2,13) = t; xi(3,13) = t
         xi(1,14) = t; xi(2,14) = s; xi(3,14) = t
         xi(1,15) = t; xi(2,15) = t; xi(3,15) = s

!     2D elements
      CASE(eType_TRI)
         w = 1._RKIND/6._RKIND
         s = 2._RKIND/3._RKIND
         t = 1._RKIND/6._RKIND
         xi(1,1) = t; xi(2,1) = t
         xi(1,2) = s; xi(2,2) = t
         xi(1,3) = t; xi(2,3) = s
      CASE(eType_BIL)
         w = 1._RKIND
         s =  1._RKIND/SQRT(3._RKIND)
         t = -1._RKIND/SQRT(3._RKIND)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
      CASE(eType_BIQ)
         w(1) = 25._RKIND/81._RKIND
         w(2) = 25._RKIND/81._RKIND
         w(3) = 25._RKIND/81._RKIND
         w(4) = 25._RKIND/81._RKIND
         w(5) = 40._RKIND/81._RKIND
         w(6) = 40._RKIND/81._RKIND
         w(7) = 40._RKIND/81._RKIND
         w(8) = 40._RKIND/81._RKIND
         w(9) = 64._RKIND/81._RKIND
         s    = SQRT(0.6_RKIND)
         xi(1,1) =  -s;      xi(2,1) =  -s
         xi(1,2) =   s;      xi(2,2) =  -s
         xi(1,3) =   s;      xi(2,3) =   s
         xi(1,4) =  -s;      xi(2,4) =   s
         xi(1,5) = 0._RKIND; xi(2,5) =  -s
         xi(1,6) =   s;      xi(2,6) = 0._RKIND
         xi(1,7) = 0._RKIND; xi(2,7) =   s
         xi(1,8) =  -s;      xi(2,8) = 0._RKIND
         xi(1,9) = 0._RKIND; xi(2,9) = 0._RKIND
      CASE(eType_QTR)
         w(1)   = 0.225000000000000_RKIND * 5E-1_RKIND
         w(2:4) = 0.125939180544827_RKIND * 5E-1_RKIND
         w(5:7) = 0.132394152788506_RKIND * 5E-1_RKIND

         s = 0.333333333333333_RKIND
         xi(1,1) = s; xi(2,1) = s

         s = 0.797426985353087_RKIND
         t = 0.101286507323456_RKIND
         xi(1,2) = t; xi(2,2) = t
         xi(1,3) = s; xi(2,3) = t
         xi(1,4) = t; xi(2,4) = s

         s = 0.059715871789770_RKIND
         t = 0.470142064105115_RKIND
         xi(1,5) = t; xi(2,5) = t
         xi(1,6) = s; xi(2,6) = t
         xi(1,7) = t; xi(2,7) = s

!     1D elements
      CASE(eType_LIN)
         w = 1._RKIND
         s = 1._RKIND/SQRT(3._RKIND)
         xi(1,1) = -s
         xi(1,2) =  s
      CASE(eType_QUD)
         w(1) = 5._RKIND/9._RKIND
         w(2) = 5._RKIND/9._RKIND
         w(3) = 8._RKIND/9._RKIND
         s = SQRT(0.6_RKIND)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0._RKIND
!     0D elements
      CASE(eType_PNT)
         w = 1._RKIND
      END SELECT

      END SUBROUTINE GETGIP
!####################################################################
!     Inverse maps {xp} to {$\xi$} in an element with coordinates {xl}
!     using Newton's method
      SUBROUTINE GETXI(eType, eNoN, xl, xp, xi, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eType, eNoN
      REAL(KIND=RKIND), INTENT(IN)  :: xl(nsd,eNoN), xp(nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: xi(nsd)
      LOGICAL, INTENT(OUT) :: flag

      INTEGER(KIND=IKIND), PARAMETER :: MAXITR = 5
      REAL(KIND=RKIND), PARAMETER :: RTOL = 1.E-6_RKIND
      REAL(KIND=RKIND), PARAMETER :: ATOL = 1.E-12_RKIND

      LOGICAL :: l1, l2, l3
      INTEGER(KIND=IKIND) :: itr, i, j, a
      REAL(KIND=RKIND) :: rmsA, rmsR, N(eNoN), Nxi(nsd,eNoN), xiK(nsd),
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
         xK = 0._RKIND
         DO i=1, nsd
            DO a=1, eNoN
               xK(i) = xK(i) + N(a)*xl(i,a)
            END DO
            rK(i) = xK(i) - xp(i)
         END DO

         rmsA = 0._RKIND
         rmsR = 0._RKIND
         DO i=1, nsd
            rmsA = rmsA + rK(i)**2._RKIND
            rmsR = rmsR + (rK(i) / (xK(i)+eps))**2._RKIND
         END DO
         rmsA = SQRT(rmsA/REAL(nsd, KIND=RKIND))
         rmsR = SQRT(rmsR/REAL(nsd, KIND=RKIND))

c         WRITE(1000+cm%tF(),'(10X,A)') "Iter: "//STR(itr)//" "//
c     2      STR(rmsA)//" "//STR(rmsR)

         l1 = itr .GT. MAXITR
         l2 = rmsA .LE. ATOL
         l3 = rmsR .LE. RTOL
         IF (l1 .OR. l2 .OR. l3) EXIT

         Am = 0._RKIND
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
      INTEGER(KIND=IKIND), INTENT(IN) :: insd, eType, eNoN
      REAL(KIND=RKIND), INTENT(IN)  :: xi(insd)
      REAL(KIND=RKIND), INTENT(OUT) :: N(eNoN), Nxi(insd,eNoN)

      REAL(KIND=RKIND) :: s, t, mx, my, ux, uy, uz, lx, ly, lz

      IF (eType .EQ. eType_NRB) RETURN

!     3D elements
      SELECT CASE(eType)
      CASE(eType_BRK)
         ux = 1._RKIND + xi(1); lx = 1._RKIND - xi(1)
         uy = 1._RKIND + xi(2); ly = 1._RKIND - xi(2)
         uz = 1._RKIND + xi(3); lz = 1._RKIND - xi(3)
         N(1) = ux*uy*uz/8._RKIND
         N(2) = lx*uy*uz/8._RKIND
         N(3) = lx*uy*lz/8._RKIND
         N(4) = ux*uy*lz/8._RKIND
         N(5) = ux*ly*uz/8._RKIND
         N(6) = lx*ly*uz/8._RKIND
         N(7) = lx*ly*lz/8._RKIND
         N(8) = ux*ly*lz/8._RKIND

         Nxi(1,1) =  uy*uz/8._RKIND
         Nxi(2,1) =  ux*uz/8._RKIND
         Nxi(3,1) =  ux*uy/8._RKIND
         Nxi(1,2) = -uy*uz/8._RKIND
         Nxi(2,2) =  lx*uz/8._RKIND
         Nxi(3,2) =  lx*uy/8._RKIND
         Nxi(1,3) = -uy*lz/8._RKIND
         Nxi(2,3) =  lx*lz/8._RKIND
         Nxi(3,3) = -lx*uy/8._RKIND
         Nxi(1,4) =  uy*lz/8._RKIND
         Nxi(2,4) =  ux*lz/8._RKIND
         Nxi(3,4) = -ux*uy/8._RKIND
         Nxi(1,5) =  ly*uz/8._RKIND
         Nxi(2,5) = -ux*uz/8._RKIND
         Nxi(3,5) =  ux*ly/8._RKIND
         Nxi(1,6) = -ly*uz/8._RKIND
         Nxi(2,6) = -lx*uz/8._RKIND
         Nxi(3,6) =  lx*ly/8._RKIND
         Nxi(1,7) = -ly*lz/8._RKIND
         Nxi(2,7) = -lx*lz/8._RKIND
         Nxi(3,7) = -lx*ly/8._RKIND
         Nxi(1,8) =  ly*lz/8._RKIND
         Nxi(2,8) = -ux*lz/8._RKIND
         Nxi(3,8) = -ux*ly/8._RKIND
      CASE(eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1._RKIND - xi(1) - xi(2) - xi(3)

         Nxi(1,1) =  1._RKIND
         Nxi(2,1) =  0._RKIND
         Nxi(3,1) =  0._RKIND
         Nxi(1,2) =  0._RKIND
         Nxi(2,2) =  1._RKIND
         Nxi(3,2) =  0._RKIND
         Nxi(1,3) =  0._RKIND
         Nxi(2,3) =  0._RKIND
         Nxi(3,3) =  1._RKIND
         Nxi(1,4) = -1._RKIND
         Nxi(2,4) = -1._RKIND
         Nxi(3,4) = -1._RKIND
      CASE(eType_WDG)
         ux = xi(1)
         uy = xi(2)
         uz = 1._RKIND - ux - uy
         s = (1._RKIND + xi(3))*0.5_RKIND
         t = (1._RKIND - xi(3))*0.5_RKIND
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

         Nxi(1,1) =  t
         Nxi(2,1) =  0._RKIND
         Nxi(3,1) = -ux*0.5_RKIND
         Nxi(1,2) =  0._RKIND
         Nxi(2,2) =  t
         Nxi(3,2) = -uy*0.5_RKIND
         Nxi(1,3) = -t
         Nxi(2,3) = -t
         Nxi(3,3) = -uz*0.5_RKIND
         Nxi(1,4) =  s
         Nxi(2,4) =  0._RKIND
         Nxi(3,4) =  ux*0.5_RKIND
         Nxi(1,5) =  0._RKIND
         Nxi(2,5) =  s
         Nxi(3,5) =  uy*0.5_RKIND
         Nxi(1,6) = -s
         Nxi(2,6) = -s
         Nxi(3,6) =  uz*0.5_RKIND
      CASE(eType_QTE)
         s     = 1._RKIND - xi(1) - xi(2) - xi(3)
         N(1)  = xi(1)*(2._RKIND*xi(1) - 1._RKIND)
         N(2)  = xi(2)*(2._RKIND*xi(2) - 1._RKIND)
         N(3)  = xi(3)*(2._RKIND*xi(3) - 1._RKIND)
         N(4)  = s    *(2._RKIND*s     - 1._RKIND)
         N(5)  = 4._RKIND*xi(1)*xi(2)
         N(6)  = 4._RKIND*xi(2)*xi(3)
         N(7)  = 4._RKIND*xi(1)*xi(3)
         N(8)  = 4._RKIND*xi(1)*s
         N(9)  = 4._RKIND*xi(2)*s
         N(10) = 4._RKIND*xi(3)*s

         Nxi(1,1)  =  4._RKIND*xi(1) - 1._RKIND
         Nxi(2,1)  =  0._RKIND
         Nxi(3,1)  =  0._RKIND
         Nxi(1,2)  =  0._RKIND
         Nxi(2,2)  =  4._RKIND*xi(2) - 1._RKIND
         Nxi(3,2)  =  0._RKIND
         Nxi(1,3)  =  0._RKIND
         Nxi(2,3)  =  0._RKIND
         Nxi(3,3)  =  4._RKIND*xi(3) - 1._RKIND
         Nxi(1,4)  =  1._RKIND - 4._RKIND*s
         Nxi(2,4)  =  1._RKIND - 4._RKIND*s
         Nxi(3,4)  =  1._RKIND - 4._RKIND*s
         Nxi(1,5)  =  4._RKIND*xi(2)
         Nxi(2,5)  =  4._RKIND*xi(1)
         Nxi(3,5)  =  0._RKIND
         Nxi(1,6)  =  0._RKIND
         Nxi(2,6)  =  4._RKIND*xi(3)
         Nxi(3,6)  =  4._RKIND*xi(2)
         Nxi(1,7)  =  4._RKIND*xi(3)
         Nxi(2,7)  =  0._RKIND
         Nxi(3,7)  =  4._RKIND*xi(1)
         Nxi(1,8)  =  4._RKIND*( s - xi(1) )
         Nxi(2,8)  = -4._RKIND*xi(1)
         Nxi(3,8)  = -4._RKIND*xi(1)
         Nxi(1,9)  = -4._RKIND*xi(2)
         Nxi(2,9)  =  4._RKIND*( s - xi(2) )
         Nxi(3,9)  = -4._RKIND*xi(2)
         Nxi(1,10) = -4._RKIND*xi(3)
         Nxi(2,10) = -4._RKIND*xi(3)
         Nxi(3,10) =  4._RKIND*( s - xi(3) )

!     2D elements
      CASE(eType_TRI)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1._RKIND - xi(1) - xi(2)

         Nxi(1,1) =  1._RKIND
         Nxi(2,1) =  0._RKIND
         Nxi(1,2) =  0._RKIND
         Nxi(2,2) =  1._RKIND
         Nxi(1,3) = -1._RKIND
         Nxi(2,3) = -1._RKIND
      CASE(eType_BIL)
         ux = 1._RKIND + xi(1)
         uy = 1._RKIND + xi(2)
         lx = 1._RKIND - xi(1)
         ly = 1._RKIND - xi(2)

         N(1) = ux*uy/4._RKIND
         N(2) = lx*uy/4._RKIND
         N(3) = lx*ly/4._RKIND
         N(4) = ux*ly/4._RKIND

         Nxi(1,1) =  uy/4._RKIND
         Nxi(2,1) =  ux/4._RKIND
         Nxi(1,2) = -uy/4._RKIND
         Nxi(2,2) =  lx/4._RKIND
         Nxi(1,3) = -ly/4._RKIND
         Nxi(2,3) = -lx/4._RKIND
         Nxi(1,4) =  ly/4._RKIND
         Nxi(2,4) = -ux/4._RKIND
      CASE(eType_BIQ)
         ux = 1._RKIND + xi(1)
         uy = 1._RKIND + xi(2)
         lx = 1._RKIND - xi(1)
         ly = 1._RKIND - xi(2)
         mx = xi(1)
         my = xi(2)

         N(1) =  mx*lx*my*ly/4._RKIND
         N(2) = -mx*ux*my*ly/4._RKIND
         N(3) =  mx*ux*my*uy/4._RKIND
         N(4) = -mx*lx*my*uy/4._RKIND
         N(5) = -lx*ux*my*ly*0.5_RKIND
         N(6) =  mx*ux*ly*uy*0.5_RKIND
         N(7) =  lx*ux*my*uy*0.5_RKIND
         N(8) = -mx*lx*ly*uy*0.5_RKIND
         N(9) =  lx*ux*ly*uy

         Nxi(1,1) =  (lx - mx)*my*ly/4._RKIND
         Nxi(2,1) =  (ly - my)*mx*lx/4._RKIND
         Nxi(1,2) = -(ux + mx)*my*ly/4._RKIND
         Nxi(2,2) = -(ly - my)*mx*ux/4._RKIND
         Nxi(1,3) =  (ux + mx)*my*uy/4._RKIND
         Nxi(2,3) =  (uy + my)*mx*ux/4._RKIND
         Nxi(1,4) = -(lx - mx)*my*uy/4._RKIND
         Nxi(2,4) = -(uy + my)*mx*lx/4._RKIND
         Nxi(1,5) = -(lx - ux)*my*ly*0.5_RKIND
         Nxi(2,5) = -(ly - my)*lx*ux*0.5_RKIND
         Nxi(1,6) =  (ux + mx)*ly*uy*0.5_RKIND
         Nxi(2,6) =  (ly - uy)*mx*ux*0.5_RKIND
         Nxi(1,7) =  (lx - ux)*my*uy*0.5_RKIND
         Nxi(2,7) =  (uy + my)*lx*ux*0.5_RKIND
         Nxi(1,8) = -(lx - mx)*ly*uy*0.5_RKIND
         Nxi(2,8) = -(ly - uy)*mx*lx*0.5_RKIND
         Nxi(1,9) =  (lx - ux)*ly*uy
         Nxi(2,9) =  (ly - uy)*lx*ux
      CASE(eType_QTR)
         s    = 1._RKIND - xi(1) - xi(2)
         N(1) = xi(1)*( 2._RKIND*xi(1) - 1._RKIND )
         N(2) = xi(2)*( 2._RKIND*xi(2) - 1._RKIND )
         N(3) = s    *( 2._RKIND*s     - 1._RKIND )
         N(4) = 4._RKIND*xi(1)*xi(2)
         N(5) = 4._RKIND*xi(2)*s
         N(6) = 4._RKIND*xi(1)*s

         Nxi(1,1) =  4._RKIND*xi(1) - 1._RKIND
         Nxi(2,1) =  0._RKIND
         Nxi(1,2) =  0._RKIND
         Nxi(2,2) =  4._RKIND*xi(2) - 1._RKIND
         Nxi(1,3) =  1._RKIND - 4._RKIND*s
         Nxi(2,3) =  1._RKIND - 4._RKIND*s
         Nxi(1,4) =  4._RKIND*xi(2)
         Nxi(2,4) =  4._RKIND*xi(1)
         Nxi(1,5) = -4._RKIND*xi(2)
         Nxi(2,5) =  4._RKIND*( s - xi(2) )
         Nxi(1,6) =  4._RKIND*( s - xi(1) )
         Nxi(2,6) = -4._RKIND*xi(1)

!     1D elements
      CASE(eType_LIN)
         N(1) = (1._RKIND - xi(1))*0.5_RKIND
         N(2) = (1._RKIND + xi(1))*0.5_RKIND

         Nxi(1,1) = -0.5_RKIND
         Nxi(1,2) =  0.5_RKIND
      CASE(eType_QUD)
         N(1) = -xi(1)*(1._RKIND - xi(1))*0.5_RKIND
         N(2) =  xi(1)*(1._RKIND + xi(1))*0.5_RKIND
         N(3) = (1._RKIND - xi(1))*(1._RKIND + xi(1))

         Nxi(1,1) = -0.5_RKIND + xi(1)
         Nxi(1,2) =  0.5_RKIND + xi(1)
         Nxi(1,3) = -2._RKIND*xi(1)

!     0D elements
      CASE(eType_PNT)
         N(1) = 1._RKIND
      END SELECT

      RETURN
      END SUBROUTINE GETGNN
!--------------------------------------------------------------------
!     Returns shape functions bounds
      PURE SUBROUTINE GETNNBNDS(eType, eNoN, xib, Nb)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eType, eNoN
      REAL(KIND=RKIND), INTENT(OUT) :: xib(2,nsd), Nb(2,eNoN)

      xib(1,:) = -1._RKIND
      xib(2,:) =  1._RKIND

      Nb(1,:)  = 0._RKIND
      Nb(2,:)  = 1._RKIND

!     3D elements
      SELECT CASE(eType)
      CASE(eType_TET)
         xib(1,:) =  0._RKIND

      CASE(eType_WDG)
         xib(1,3) = -1._RKIND

      CASE(eType_QTE)
         xib(1,:) =  0._RKIND

         Nb(1,1)  = -0.125_RKIND
         Nb(1,2)  = -0.125_RKIND
         Nb(1,3)  = -0.125_RKIND
         Nb(1,4)  = -0.125_RKIND

         Nb(2,5)  =  4._RKIND
         Nb(2,6)  =  4._RKIND
         Nb(2,7)  =  4._RKIND
         Nb(2,8)  =  4._RKIND
         Nb(2,9)  =  4._RKIND
         Nb(2,10) =  4._RKIND

!     2D elements
      CASE(eType_TRI)
         xib(1,:) =  0._RKIND

      CASE(eType_BIQ)
         Nb(1,1)  = -0.125_RKIND
         Nb(1,2)  = -0.125_RKIND
         Nb(1,3)  = -0.125_RKIND
         Nb(1,4)  = -0.125_RKIND
         Nb(1,5)  = -0.125_RKIND
         Nb(1,6)  = -0.125_RKIND
         Nb(1,7)  = -0.125_RKIND
         Nb(1,8)  = -0.125_RKIND
         Nb(1,9)  =  0._RKIND

      CASE(eType_QTR)
         xib(1,:) =  0._RKIND

         Nb(1,1)  = -0.125_RKIND
         Nb(1,2)  = -0.125_RKIND
         Nb(1,3)  = -0.125_RKIND

         Nb(2,4)  =  4._RKIND
         Nb(2,5)  =  4._RKIND
         Nb(2,6)  =  4._RKIND

!     1D elements
      CASE(eType_QUD)
         Nb(1,1)  = -0.125_RKIND
         Nb(1,2)  = -0.125_RKIND
         Nb(1,3)  =  0._RKIND

      END SELECT

!     Add a small tolerance around the bounds
      xib(1,:) = xib(1,:) - 1.0E-4_RKIND
      xib(2,:) = xib(2,:) + 1.0E-4_RKIND

      Nb(1,:)  = Nb(1,:)  - 1.0E-4_RKIND
      Nb(2,:)  = Nb(2,:)  + 1.0E-4_RKIND

      RETURN
      END SUBROUTINE GETNNBNDS
!####################################################################
      PURE SUBROUTINE GNN(eNoN, insd, Nxi, x, Nx, Jac, ks)
      USE COMMOD, ONLY: nsd
      USE UTILMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN, insd
      REAL(KIND=RKIND), INTENT(IN) :: Nxi(insd,eNoN), x(nsd,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: Nx(insd,eNoN), Jac, ks(nsd,nsd)

      INTEGER(KIND=IKIND) a
      REAL(KIND=RKIND) xXi(nsd,insd), xiX(insd,nsd)

      xXi = 0._RKIND
      Jac = 0._RKIND
      Nx  = 0._RKIND
      ks  = 0._RKIND
      IF (insd .EQ. 1) THEN
         DO a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
         END DO

         Jac = SQRT(NORM(xXi)) + 1.E+3_RKIND*eps
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
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: Nxi(nsd-1,eNoN), xl(nsd,eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: nV(nsd), gCov(nsd,nsd-1),
     2   gCnv(nsd,nsd-1)

      INTEGER(KIND=IKIND) a, i, j, insd
      REAL(KIND=RKIND), ALLOCATABLE :: xXi(:,:), Gmat(:,:)

      insd = nsd - 1
      ALLOCATE(xXi(nsd,insd), Gmat(insd,insd))

!     Calculating surface deflation
      xXi = 0._RKIND
      DO a=1, eNoN
         DO i=1, insd
            xXi(:,i) = xXi(:,i) + xl(:,a)*Nxi(i,a)
         END DO
      END DO
      nV = CROSS(xXi)

!     Covariant basis
      gCov = xXi

!     Metric tensor g_i . g_j
      Gmat = 0._RKIND
      DO i=1, insd
         DO j=1, insd
            DO a=1, nsd
               Gmat(i,j) = Gmat(i,j) + gCov(a,i)*gCov(a,j)
            END DO
         END DO
      END DO

!     Contravariant basis
      Gmat = MAT_INV(Gmat, insd)
      gCnv = 0._RKIND
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
      INTEGER(KIND=IKIND), INTENT(IN) :: e, g
      REAL(KIND=RKIND), INTENT(OUT) :: n(nsd)
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) a, Ac, i, iM, Ec, b, Bc, eNoN, insd
      REAL(KIND=RKIND) v(nsd)

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lX(:,:), xXi(:,:)

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
         xXi = 0._RKIND
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
         xXi = 0._RKIND
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
         v(:) = 0._RKIND
         DO a=1, eNoN
            v(:) = v(:) + lX(:,a)*ib%msh(iM)%N(a,g)
         END DO
         a = ptr(1)
         v(:) = lX(:,a) - v(:)
         IF (NORM(n,v) .LT. 0._RKIND) n = -n

         DEALLOCATE(xXi)
         RETURN
      ELSE
         ALLOCATE(xXi(nsd,insd))
         xXi = 0._RKIND
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
      IF (NORM(n,v) .LT. 0._RKIND) n = -n

      RETURN
      END SUBROUTINE GNNIB
!--------------------------------------------------------------------
!     This routine returns a vector at element "e" and Gauss point
!     "g" of face "lFa" that is the normal weigthed by Jac, i.e.
!     Jac = SQRT(NORM(n)).
      SUBROUTINE GNNB(lFa, e, g, insd, eNoNb, Nx, n)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: e, g, insd, eNoNb
      REAL(KIND=RKIND), INTENT(OUT) :: n(nsd), Nx(insd,eNoNb)
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) a, Ac, i, iM, Ec, b, Bc, eNoN
      REAL(KIND=RKIND) v(nsd)

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lX(:,:), xXi(:,:)

      iM   = lFa%iM
      Ec   = lFa%gE(e)
      eNoN = msh(iM)%eNoN

      ALLOCATE(lX(nsd,eNoN), ptr(eNoN), setIt(eNoN))

!     Creating a ptr list that contains pointer to the nodes of elements
!     that are at the face at the beginning of the list and the rest at
!     the end
      setIt = .TRUE.
      DO a=1, eNoNb
         Ac = lFa%IEN(a,e)
         DO b=1, eNoN
            IF (setIt(b)) THEN
               Bc = msh(iM)%IEN(b,Ec)
               IF (Bc .EQ. Ac) EXIT
            END IF
         END DO
         IF (b .GT. eNoN) THEN
            WRITE(*,'(A)')
            WRITE(*,'(A)') "=========================================="
            WRITE(*,'(A)') " ERROR: could not find matching face nodes"
            WRITE(*,'(A)',ADVANCE='NO') "    Face "//TRIM(lFa%name)//
     2         " e: "//STR(e)
            DO b=1, eNoNb
               WRITE(*,'(A)',ADVANCE='NO') " "//STR(lFa%IEN(b,e))
            END DO
            WRITE(*,'(A)')
            WRITE(*,'(A)',ADVANCE='NO') "    Mesh "//
     2         TRIM(msh(iM)%name)//" Ec: "//STR(Ec)
            DO b=1, eNoN
               WRITE(*,'(A)',ADVANCE='NO') " "//STR(msh(iM)%IEN(b,Ec))
            END DO
            WRITE(*,'(A)')
            WRITE(*,'(A)') "=========================================="
            WRITE(*,'(A)')
            CALL STOPSIM()
         END IF
         ptr(a)   = b
         setIt(b) = .FALSE.
      END DO
      a = eNoNb
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
         xXi = 0._RKIND
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
         xXi = 0._RKIND
         DO a=1, eNoNb
            b = ptr(a)
            xXi(:,1) = xXi(:,1) + lFa%Nx(1,a,g)*lX(:,b)
         END DO

!        Face normal
         n(1) = v(2)*xXi(3,1) - v(3)*xXi(2,1)
         n(2) = v(3)*xXi(1,1) - v(1)*xXi(3,1)
         n(3) = v(1)*xXi(2,1) - v(2)*xXi(1,1)

!        I choose Gauss point of the mesh element for calculating
!        interior edge
         v(:) = 0._RKIND
         DO a=1, eNoN
            v(:) = v(:) + lX(:,a)*msh(iM)%N(a,g)
         END DO
         a = ptr(1)
         v(:) = lX(:,a) - v(:)
         IF (NORM(n,v) .LT. 0._RKIND) n = -n

         DEALLOCATE(xXi)
         RETURN
      ELSE
         ALLOCATE(xXi(nsd,insd))
         xXi = 0._RKIND
         DO a=1, eNoNb
            b = ptr(a)
            DO i=1, insd
               xXi(:,i) = xXi(:,i) + Nx(i,a)*lX(:,b)
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
      IF (NORM(n,v) .LT. 0._RKIND) n = -n

      RETURN
      END SUBROUTINE GNNB
!####################################################################


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
      
      IF (lM%eType .EQ. eType_NRB) THEN
         ALLOCATE(lM%w(lM%nG), lM%N(lM%eNoN,lM%nG), 
     2      lM%Nx(nsd,lM%eNoN,lM%nG))
         RETURN
      END IF

      IF (nsd .EQ. 3) THEN
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
      ELSE
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
      END IF
      
      ALLOCATE(lM%w(lM%nG), lM%N(lM%eNoN,lM%nG), 
     2   lM%Nx(nsd,lM%eNoN,lM%nG))

      CALL XIGAUSS(nsd, lM%eType, lM%nG, lM%eNoN, lM%w, lM%N, lM%Nx)

      RETURN
      END SUBROUTINE SELECTELE

!--------------------------------------------------------------------
!     This routine selects boundary element type
      SUBROUTINE SELECTELEB(lM, lFa)
      
      USE COMMOD

      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      IF (lM%eType .EQ. eType_NRB) THEN
         lFa%eType = eType_NRB
         lFa%nG    = lM%nG/lM%bs(lFa%d)%nG
         lFa%eNoN  = lM%eNoN/(lM%bs(lFa%d)%p + 1)

         ALLOCATE(lFa%w(lFa%nG), lFa%N(lFa%eNoN,lFa%nG),
     2      lFa%Nx(nsd-1,lFa%eNoN,lFa%nG))
         RETURN
      END IF

      IF (nsd .EQ. 3) THEN
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
      ELSE
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
      END IF
      
      ALLOCATE(lFa%w(lFa%nG), lFa%N(lFa%eNoN,lFa%nG),
     2   lFa%Nx(nsd-1,lFa%eNoN,lFa%nG))

      CALL XIGAUSS(nsd-1, lFa%eType, lFa%nG, lFa%eNoN, lFa%w, lFa%N, 
     2   lFa%Nx)

      RETURN
      END SUBROUTINE SELECTELEB

!####################################################################
!     Assigning values to N and Nx
      PURE SUBROUTINE XIGAUSS(insd, eType, nG, eNoN, w, N, Nx)
      
      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: insd, eType, nG, eNoN
      REAL(KIND=8), INTENT(OUT) :: w(nG), N(eNoN,nG), Nx(insd,eNoN,nG)

      INTEGER g
      REAL(KIND=8) s, t, mx, my, ux, uy, uz, lx, ly, lz
      REAL(KIND=8), ALLOCATABLE :: xi(:,:)

      IF (eType .EQ. eType_NRB) RETURN
      ALLOCATE(xi(insd,nG))

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
         DO g=1, nG
            ux = 1D0 + xi(1,g); lx = 1D0 - xi(1,g)
            uy = 1D0 + xi(2,g); ly = 1D0 - xi(2,g)
            uz = 1D0 + xi(3,g); lz = 1D0 - xi(3,g)
            N(1,g) = ux*uy*uz/8D0
            N(2,g) = lx*uy*uz/8D0
            N(3,g) = lx*uy*lz/8D0
            N(4,g) = ux*uy*lz/8D0
            N(5,g) = ux*ly*uz/8D0
            N(6,g) = lx*ly*uz/8D0
            N(7,g) = lx*ly*lz/8D0
            N(8,g) = ux*ly*lz/8D0

            Nx(1,1,g) =  uy*uz/8D0
            Nx(2,1,g) =  ux*uz/8D0
            Nx(3,1,g) =  ux*uy/8D0
            Nx(1,2,g) = -uy*uz/8D0
            Nx(2,2,g) =  lx*uz/8D0
            Nx(3,2,g) =  lx*uy/8D0
            Nx(1,3,g) = -uy*lz/8D0
            Nx(2,3,g) =  lx*lz/8D0
            Nx(3,3,g) = -lx*uy/8D0
            Nx(1,4,g) =  uy*lz/8D0
            Nx(2,4,g) =  ux*lz/8D0
            Nx(3,4,g) = -ux*uy/8D0
            Nx(1,5,g) =  ly*uz/8D0
            Nx(2,5,g) = -ux*uz/8D0
            Nx(3,5,g) =  ux*ly/8D0
            Nx(1,6,g) = -ly*uz/8D0
            Nx(2,6,g) = -lx*uz/8D0
            Nx(3,6,g) =  lx*ly/8D0
            Nx(1,7,g) = -ly*lz/8D0
            Nx(2,7,g) = -lx*lz/8D0
            Nx(3,7,g) = -lx*ly/8D0
            Nx(1,8,g) =  ly*lz/8D0
            Nx(2,8,g) = -ux*lz/8D0
            Nx(3,8,g) = -ux*ly/8D0
         END DO
      CASE(eType_TET)
         w = 1D0/24D0
         s = (5D0 + 3D0*SQRT(5D0))/2D1
         t = (5D0 -     SQRT(5D0))/2D1
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
         DO g=1, nG
            N(1,g) = xi(1,g)
            N(2,g) = xi(2,g)
            N(3,g) = xi(3,g)
            N(4,g) = 1D0 - xi(1,g) - xi(2,g) - xi(3,g)
            
            Nx(1,1,g) =  1D0
            Nx(2,1,g) =  0D0
            Nx(3,1,g) =  0D0
            Nx(1,2,g) =  0D0
            Nx(2,2,g) =  1D0
            Nx(3,2,g) =  0D0
            Nx(1,3,g) =  0D0
            Nx(2,3,g) =  0D0
            Nx(3,3,g) =  1D0
            Nx(1,4,g) = -1D0
            Nx(2,4,g) = -1D0
            Nx(3,4,g) = -1D0
         END DO
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
         DO g=1, nG
            ux = xi(1,g) ; uy = xi(2,g) ; uz = 1D0 - uz - uy
            s = (1D0 + xi(3,g))/2D0; t = (1D0 - xi(3,g))/2D0
            N(1,g) = ux*t
            N(2,g) = uy*t
            N(3,g) = uz*t
            N(4,g) = ux*s
            N(5,g) = uy*s
            N(6,g) = uz*s
            
            Nx(1,1,g) =  t
            Nx(2,1,g) =  0D0
            Nx(3,1,g) = -ux/2D0
            Nx(1,2,g) =  0D0
            Nx(2,2,g) =  t
            Nx(3,2,g) = -uy/2D0
            Nx(1,3,g) = -t
            Nx(2,3,g) = -t
            Nx(3,3,g) = -uz/2D0
            Nx(1,4,g) =  s
            Nx(2,4,g) =  0D0
            Nx(3,4,g) =  ux/2D0
            Nx(1,5,g) =  0D0
            Nx(2,5,g) =  s
            Nx(3,5,g) =  uy/2D0
            Nx(1,6,g) = -s
            Nx(2,6,g) = -s
            Nx(3,6,g) =  uz/2D0
         END DO

!     2D elements         
      CASE(eType_TRI)
         w = 1D0/6D0
         s = 2D0/3D0
         t = 1D0/6D0
         xi(1,1) = s; xi(2,1) = t
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         DO g=1, nG
            N(1,g) = xi(1,g)
            N(2,g) = xi(2,g)
            N(3,g) = 1D0 - xi(1,g) - xi(2,g)
            
            Nx(1,1,g) =  1D0
            Nx(2,1,g) =  0D0
            Nx(1,2,g) =  0D0
            Nx(2,2,g) =  1D0
            Nx(1,3,g) = -1D0
            Nx(2,3,g) = -1D0
         END DO
      CASE(eType_BIL)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
         DO g=1, nG
            ux = 1D0 + xi(1,g); lx = 1D0 - xi(1,g)
            uy = 1D0 + xi(2,g); ly = 1D0 - xi(2,g)
            N(1,g) = ux*uy/4D0
            N(2,g) = lx*uy/4D0
            N(3,g) = lx*ly/4D0
            N(4,g) = ux*ly/4D0

            Nx(1,1,g) =  uy/4D0
            Nx(2,1,g) =  ux/4D0
            Nx(1,2,g) = -uy/4D0
            Nx(2,2,g) =  lx/4D0
            Nx(1,3,g) = -ly/4D0
            Nx(2,3,g) = -lx/4D0
            Nx(1,4,g) =  ly/4D0
            Nx(2,4,g) = -ux/4D0
         END DO
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
         DO g=1, nG
            ux = 1D0 + xi(1,g); mx = xi(1,g); lx = 1D0 - xi(1,g)
            uy = 1D0 + xi(2,g); my = xi(2,g); ly = 1D0 - xi(2,g)
            N(1,g) =  mx*lx*my*ly/4D0
            N(2,g) = -mx*ux*my*ly/4D0
            N(3,g) =  mx*ux*my*uy/4D0
            N(4,g) = -mx*lx*my*uy/4D0
            N(5,g) = -lx*ux*my*ly/2D0
            N(6,g) =  mx*ux*ly*uy/2D0
            N(7,g) =  lx*ux*my*uy/2D0
            N(8,g) = -mx*lx*ly*uy/2D0
            N(9,g) =  lx*ux*ly*uy
            
            Nx(1,1,g) =  (lx - mx)*my*ly/4D0
            Nx(2,1,g) =  (ly - my)*mx*lx/4D0
            Nx(1,2,g) = -(ux + mx)*my*ly/4D0
            Nx(2,2,g) = -(ly - my)*mx*ux/4D0
            Nx(1,3,g) =  (ux + mx)*my*uy/4D0
            Nx(2,3,g) =  (uy + my)*mx*ux/4D0
            Nx(1,4,g) = -(lx - mx)*my*uy/4D0
            Nx(2,4,g) = -(uy + my)*mx*lx/4D0
            Nx(1,5,g) = -(lx - ux)*my*ly/2D0
            Nx(2,5,g) = -(ly - my)*lx*ux/2D0
            Nx(1,6,g) =  (ux + mx)*ly*uy/2D0
            Nx(2,6,g) =  (ly - uy)*mx*ux/2D0
            Nx(1,7,g) =  (lx - ux)*my*uy/2D0
            Nx(2,7,g) =  (uy + my)*lx*ux/2D0
            Nx(1,8,g) = -(lx - mx)*ly*uy/2D0
            Nx(2,8,g) = -(ly - uy)*mx*lx/2D0
            Nx(1,9,g) =  (lx - ux)*ly*uy
            Nx(2,9,g) =  (ly - uy)*lx*ux
         END DO

!     1D elements         
      CASE(eType_LIN)
         w = 1D0
         s = 1D0/SQRT(3D0)
         xi(1,1) = -s
         xi(1,2) =  s
         DO g=1, nG
            N(1,g) = (1D0 - xi(1,g))/2D0
            N(2,g) = (1D0 + xi(1,g))/2D0
            
            Nx(1,1,g) = -5D-1
            Nx(1,2,g) =  5D-1
         END DO
      CASE(eType_QUD)
         w(1) = 5D0/9D0; w(2) = 5D0/9D0; w(3) = 8D0/9D0
         s = SQRT(6D-1)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0D0
         DO g=1, nG
            N(1,g) = -xi(1,g)*(1D0 - xi(1,g))/2D0
            N(2,g) =  xi(1,g)*(1D0 + xi(1,g))/2D0
            N(3,g) = (1D0 - xi(1,g))*(1D0 + xi(1,g))
            
            Nx(1,1,g) = -5D-1 + xi(1,g)
            Nx(1,2,g) =  5D-1 + xi(1,g)
            Nx(1,3,g) = -2D0*xi(1,g)
         END DO
      END SELECT

      RETURN
      END SUBROUTINE XIGAUSS

!####################################################################
      PURE SUBROUTINE GNN(eNoN, Nxi, x, Nx, Jac, ks)

      USE COMMOD, ONLY: nsd

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=8), INTENT(IN) :: Nxi(nsd,eNoN), x(nsd,eNoN)
      REAL(KIND=8), INTENT(OUT) :: Nx(nsd,eNoN), Jac, ks(nsd,nsd)

      INTEGER a
      REAL(KIND=8) xXi(nsd,nsd), xiX(nsd,nsd)

      Nx  = 0D0
      xXi = 0D0
      IF (nsd .EQ. 2) THEN
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
      ELSE
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
!     This routine returns a vector at element "e" and Gauss point 
!     "g" of face "lFa" that is the normal weigthed by Jac, i.e. 
!     Jac = SQRT(NORM(n)). 
      PURE SUBROUTINE GNNB(lFa, e, g, n)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e, g
      REAL(KIND=8), INTENT(OUT) :: n(nsd)
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER a, Ac, i, iM, Ec, b, Bc, eNoN
      REAL(KIND=8) xXi(nsd,nsd-1), v(nsd)

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=8), ALLOCATABLE :: lX(:,:)

      iM   = lFa%iM
      Ec   = lFa%gE(e)
      eNoN = msh(iM)%eNoN
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
      xXi = 0D0
      DO a=1, lFa%eNoN
         b = ptr(a)
         DO i=1, nsd-1
            xXi(:,i) = xXi(:,i) + lFa%Nx(i,a,g)*lX(:,b)
         END DO
      END DO
      n = CROSS(xXi)

!     Changing the sign if neccessary. a locates on the face and b
!     outside of the face, in the parent element
      a = ptr(1)
      b = ptr(lFa%eNoN+1)
      v = lX(:,a) - lX(:,b)
      IF (NORM(n,v) .LT. 0D0) n = -n 

      RETURN
      END SUBROUTINE GNNB

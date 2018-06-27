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
!     This file contains routines for reading, calculating shape
!     function and computing NURBS
!      
!--------------------------------------------------------------------
!     Calculates shape function values "N" and first derivatives "Nxi"
!     of BSpline "bs" at Gauss points, given the knot "ni"
      PURE SUBROUTINE BSPNNX(ni, bs, N, Nxi)
         
      USE COMMOD, ONLY: bsType, gXi

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ni
      TYPE(bsType), INTENT(IN) :: bs
      REAL(KIND=8), INTENT(OUT) :: N(bs%p+1,bs%nG), Nxi(bs%p+1,bs%nG)

      INTEGER p, i, g
      REAL(KIND=8) xi, saved, tmp, c
      REAL(KIND=8), ALLOCATABLE :: l(:), r(:), Nm(:)
      
      ALLOCATE(l(bs%p), r(bs%p), Nm(bs%p))
      DO g=1, bs%nG
         xi = (bs%xi(ni+1) + bs%xi(ni) 
     2      + (bs%xi(ni+1) - bs%xi(ni))*gXi(g,bs%nG))/2D0
      
         N(1,g) = 1D0
         DO p=1, bs%p
            l(p) = xi - bs%xi(ni-p+1)
            r(p) = bs%xi(ni+p) - xi
            saved = 0D0
            DO i=1, p
               Nm(i) = N(i,g)
               tmp = N(i,g)/(r(i) + l(p-i+1))
               N(i,g) = saved + tmp*r(i)
               saved = tmp*l(p-i+1)
            END DO
            N(p+1,g) = saved
         END DO

!     Nm_(i-1,p) is N_(i,p-1). Based on Nm I will calculate derivative.
!     Note that Nxi is shifted by "ni - p - 1" compared to bs%xi.
         p     = bs%p
         saved = 0D0
         c     = REAL(bs%p,8)*(bs%xi(ni+1) - bs%xi(ni))/2D0
         DO i=1, p
            tmp      = c*Nm(i)/(bs%xi(ni+i) - bs%xi(ni-p+i))
            Nxi(i,g) = saved - tmp
            saved    = tmp
         END DO
         Nxi(p+1,g)  = saved
      END DO

      RETURN
      END SUBROUTINE BSPNNX
!--------------------------------------------------------------------
!     Calculates shape function values "NG" and derivative "NxiG" 
!     of NURB "nrb" at the Gauss points, given the element "e"
      SUBROUTINE NRBNNX(lM, e)
         
      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e
      TYPE(mshType), INTENT(INOUT) :: lM

      TYPE dType  
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
         REAL(KIND=8), ALLOCATABLE :: Nxi(:,:)
      END TYPE dType 

      INTEGER i, Ac, a, ax, ay, az, g, gx, gy, gz
      REAL(KIND=8) s, sx, sy, sz
      TYPE(dType) d(nsd)

!     First calculating shape functions of BSplines in each direction
      DO i=1, nsd
         ALLOCATE(d(i)%N(lM%bs(i)%p+1,lM%bs(i)%nG),
     2          d(i)%Nxi(lM%bs(i)%p+1,lM%bs(i)%nG))
         CALL BSPNNX(lM%INN(i,e), lM%bs(i), d(i)%N, d(i)%Nxi)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (nsd .EQ. 2) THEN
         DO gx=1, lM%bs(1)%nG
            DO gy=1, lM%bs(2)%nG
               g  = g + 1
               a  = 0
               s  = 0D0
               sx = 0D0
               sy = 0D0
               DO ax=1, lM%bs(1)%p + 1
                  DO ay=1, lM%bs(2)%p + 1
                     a  = a + 1
                     Ac = lM%IEN(a,e)
                     Ac = lM%lN(Ac)

                     lM%N(a,g)    = d(1)%N(ax,gx)*  d(2)%N(ay,gy)
     2                            * lM%nW(Ac)
                     lM%Nx(1,a,g) = d(1)%Nxi(ax,gx)*d(2)%N(ay,gy)
     2                            * lM%nW(Ac)
                     lM%Nx(2,a,g) = d(1)%N(ax,gx)*  d(2)%Nxi(ay,gy)
     2                            * lM%nW(Ac)
                     
                     s  = s  + lM%N(a,g)
                     sx = sx + lM%Nx(1,a,g)
                     sy = sy + lM%Nx(2,a,g)
                  END DO
               END DO
               DO a=1, lM%eNoN
                  lM%N(a,g)    = lM%N(a,g)/s
                  lM%Nx(1,a,g) = (lM%Nx(1,a,g) - lM%N(a,g)*sx)/s
                  lM%Nx(2,a,g) = (lM%Nx(2,a,g) - lM%N(a,g)*sy)/s
               END DO
               lM%w(g) = gW(gx,lM%bs(1)%nG)*gW(gy,lM%bs(2)%nG)
            END DO
         END DO
      ELSE
         DO gx=1, lM%bs(1)%nG
            DO gy=1, lM%bs(2)%nG
               DO gz=1, lM%bs(3)%nG
                  g  = g + 1
                  a  = 0
                  s  = 0D0
                  sx = 0D0
                  sy = 0D0
                  sz = 0D0
                  DO ax=1, lM%bs(1)%p + 1
                     DO ay=1, lM%bs(2)%p + 1
                        DO az=1, lM%bs(3)%p + 1
                           a  = a + 1
                           Ac = lM%IEN(a,e)
                           Ac = lM%lN(Ac)
                           
                           lM%N(a,g)    = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(1,a,g) = d(1)%Nxi(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(2,a,g) = d(1)%N(ax,gx)*d(2)%Nxi(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(3,a,g) = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%Nxi(az,gz)*lM%nW(Ac)
                           
                           s  = s  + lM%N(a,g)
                           sx = sx + lM%Nx(1,a,g)
                           sy = sy + lM%Nx(2,a,g)
                           sz = sz + lM%Nx(3,a,g)
                        END DO
                     END DO
                  END DO
                  DO a=1, lM%eNoN
                     lM%N(a,g)    = lM%N(a,g)/s
                     lM%Nx(1,a,g) = (lM%Nx(1,a,g) - lM%N(a,g)*sx)/s
                     lM%Nx(2,a,g) = (lM%Nx(2,a,g) - lM%N(a,g)*sy)/s
                     lM%Nx(3,a,g) = (lM%Nx(3,a,g) - lM%N(a,g)*sz)/s
                  END DO
                  lM%w(g) = gW(gx,lM%bs(1)%nG)*gW(gy,lM%bs(2)%nG)
     2                    * gW(gz,lM%bs(3)%nG)
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE NRBNNX
!--------------------------------------------------------------------
!     Calculates shape function values "N" and derivative "Nx" 
!     of NURB "lM" at boundary the Gauss points, given the element "e"
      SUBROUTINE NRBNNXB(lM, lFa, e)
         
      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      TYPE dType  
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
         REAL(KIND=8), ALLOCATABLE :: Nxi(:,:)
      END TYPE dType 

      INTEGER i, j, Ac, a, ax, ay, g, gx, gy, Ec
      REAL(KIND=8) s, sx, sy
      TYPE(dType) d(nsd)

      Ec = lFa%gE(e)
!     First calculating shape functions of BSplines in each direction,
!     exept the direction that face is normal to
      DO i=1, nsd
         ALLOCATE(d(i)%N(lM%bs(i)%p+1,lM%bs(i)%nG),
     2          d(i)%Nxi(lM%bs(i)%p+1,lM%bs(i)%nG))
         CALL BSPNNX(lM%INN(i,Ec), lM%bs(i), d(i)%N, d(i)%Nxi)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (nsd .EQ. 2) THEN
         IF (lFa%d .EQ. 1) i = 2
         IF (lFa%d .EQ. 2) i = 1
         DO g=1, lM%bs(i)%nG
            s  = 0D0
            sx = 0D0
            DO a=1, lM%bs(i)%p + 1
               Ac = lFa%IEN(a,e)
               Ac = lM%lN(Ac)

               lFa%N(a,g)    = d(i)%N(a,g)*lM%nW(Ac)
               lFa%Nx(1,a,g) = d(i)%Nxi(a,g)*lM%nW(Ac)
 
               s  = s  + lFa%N(a,g)
               sx = sx + lFa%Nx(1,a,g)
            END DO
            DO a=1, lFa%eNoN
               lFa%N(a,g)    = lFa%N(a,g)/s
               lFa%Nx(1,a,g) = (lFa%Nx(1,a,g) - lFa%N(a,g)*sx)/s
            END DO
            lFa%w(g) = gW(g,lM%bs(i)%nG)
         END DO
      ELSE
         IF (lFa%d .EQ. 1) THEN
            i = 2; j = 3
         END IF
         IF (lFa%d .EQ. 2) THEN
            i = 1; j = 3
         END IF
         IF (lFa%d .EQ. 3) THEN
            i = 1; j = 2
         END IF
         DO gx=1, lM%bs(i)%nG
            DO gy=1, lM%bs(j)%nG
               g  = g + 1
               a  = 0
               s  = 0D0
               sx = 0D0
               sy = 0D0
               DO ax=1, lM%bs(i)%p + 1
                  DO ay=1, lM%bs(j)%p + 1
                     a  = a + 1
                     Ac = lFa%IEN(a,e)
                     Ac = lM%lN(Ac)

                     lFa%N(a,g)    = d(i)%N(ax,gx)*  d(j)%N(ay,gy)
     2                             * lM%nW(Ac)
                     lFa%Nx(1,a,g) = d(i)%Nxi(ax,gx)*d(j)%N(ay,gy)
     2                             * lM%nW(Ac)
                     lFa%Nx(2,a,g) = d(i)%N(ax,gx)*  d(j)%Nxi(ay,gy)
     2                             * lM%nW(Ac)
                     
                     s  = s  + lFa%N(a,g)
                     sx = sx + lFa%Nx(1,a,g)
                     sy = sy + lFa%Nx(2,a,g)
                  END DO
               END DO
               DO a=1, lFa%eNoN
                  lFa%N(a,g)    = lFa%N(a,g)/s
                  lFa%Nx(1,a,g) = (lFa%Nx(1,a,g) - lFa%N(a,g)*sx)/s
                  lFa%Nx(2,a,g) = (lFa%Nx(2,a,g) - lFa%N(a,g)*sy)/s
               END DO
               lFa%w(g) = gW(gx,lM%bs(i)%nG)*gW(gy,lM%bs(j)%nG)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE NRBNNXB

!####################################################################
!     Calculates shape function values "N" of BSpline "bs" at "nSl" 
!     sample points uniformly distributed between the knot "ni" and
!     "ni+1"
      PURE SUBROUTINE BSPNNS(ni, bs, N)
 
      USE COMMOD, ONLY: bsType

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ni
      TYPE(bsType), INTENT(IN) :: bs
      REAL(KIND=8), INTENT(OUT) :: N(bs%p+1,bs%nSl)

      INTEGER p, i, s
      REAL(KIND=8) xi, saved, tmp, c
      REAL(KIND=8), ALLOCATABLE :: l(:), r(:)
      
      c = (bs%xi(ni+1) - bs%xi(ni))/REAL(bs%nSl - 1,8)
      ALLOCATE(l(bs%p), r(bs%p))
      DO s=1, bs%nSl
         xi = bs%xi(ni) + c*REAL(s-1,8)
      
         N(1,s) = 1D0
         DO p=1, bs%p
            l(p) = xi - bs%xi(ni-p+1)
            r(p) = bs%xi(ni+p) - xi
            saved = 0D0
            DO i=1, p
               tmp = N(i,s)/(r(i) + l(p-i+1))
               N(i,s) = saved + tmp*r(i)
               saved = tmp*l(p-i+1)
            END DO
            N(p+1,s) = saved
         END DO
      END DO

      RETURN
      END SUBROUTINE BSPNNS
!--------------------------------------------------------------------
!     Calculates shape function values "N" of NURB "nrb" at "nSl" 
!     Sample points, given the element "e"
      SUBROUTINE NRBNNS(lM, N, e)
         
      USE COMMOD

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: e
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(OUT) :: N(lM%eNoN,lM%nSl)

      TYPE dType  
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
      END TYPE dType 

      INTEGER i, Ac, a, ax, ay, az, g, gx, gy, gz
      REAL(KIND=8) s
      TYPE(dType) d(nsd)

!     First calculating shape functions of BSplines in each direction
      DO i=1, nsd
         ALLOCATE(d(i)%N(lM%bs(i)%p+1,lM%bs(i)%nSl))
         CALL BSPNNS(lM%INN(i,e), lM%bs(i), d(i)%N)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (nsd .EQ. 2) THEN
         DO gx=1, lM%bs(1)%nSl
            DO gy=1, lM%bs(2)%nSl
               g = g + 1
               a = 0
               s = 0D0
               DO ax=1, lM%bs(1)%p + 1
                  DO ay=1, lM%bs(2)%p + 1
                     a      = a + 1
                     Ac     = lM%IEN(a,e)
                     Ac     = lM%lN(Ac)
                     N(a,g) = d(1)%N(ax,gx)*d(2)%N(ay,gy)*lM%nW(Ac)
                     s      = s + N(a,g)
                  END DO
               END DO
               N(:,g) = N(:,g)/s
            END DO
         END DO
      ELSE
         DO gx=1, lM%bs(1)%nSl
            DO gy=1, lM%bs(2)%nSl
               DO gz=1, lM%bs(3)%nSl
                  g = g + 1
                  a = 0
                  s = 0D0
                  DO ax=1, lM%bs(1)%p + 1
                     DO ay=1, lM%bs(2)%p + 1
                        DO az=1, lM%bs(3)%p + 1
                           a      = a + 1
                           Ac     = lM%IEN(a,e)
                           Ac     = lM%lN(Ac)
                           N(a,g) = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                            * d(3)%N(az,gz)*lM%nW(Ac)
                           s      = s + N(a,g)
                        END DO
                     END DO
                  END DO
                  N(:,g) = N(:,g)/s
               END DO
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE NRBNNS



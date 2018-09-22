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
!     function and computing NURBS.
!
!--------------------------------------------------------------------

!     Calculates shape function values "N" and first and second order
!     derivatives "Nxi, Nxxi" of BSpline "bs" at Gauss points, given
!     the knot "ni".
!     Reference: The NURBS Book, Les Piegl & Wayne Tiller
      PURE SUBROUTINE BSPNNX(ni, bs, N, Nx, Nxx)

      USE COMMOD, ONLY: bsType, gXi

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ni
      TYPE(bsType), INTENT(IN) :: bs
      REAL(KIND=8), INTENT(OUT) :: N(bs%p+1,bs%nG), Nx(bs%p+1,bs%nG),
     2   Nxx(bs%p+1,bs%nG)

      INTEGER p, i, g, j, k, ik, pk, s1, s2, j1, j2
      REAL(KIND=8) xi, dxi, saved, tmp, d
      REAL(KIND=8), ALLOCATABLE :: l(:), r(:), NdX(:,:), a(:,:),
     2   ders(:,:)

      dxi = (bs%xi(ni+1) - bs%xi(ni))/2D0
      DO g=1, bs%nG
         xi = (bs%xi(ni+1) + bs%xi(ni)
     2      + (bs%xi(ni+1) - bs%xi(ni))*gXi(g,bs%nG))/2D0

!        NdXi: stores basis functions and knot differences
!        ders: stores (p+1) derivatives of (p+1) basis functions
         ALLOCATE(l(bs%p), r(bs%p), NdX(0:bs%p,0:bs%p), a(0:bs%p,2),
     2      ders(0:bs%p,0:bs%p))
         NdX(0,0) = 1D0
         DO p=1, bs%p
            l(p)  = xi - bs%xi(ni+1-p)
            r(p)  = bs%xi(ni+p) - xi
            saved = 0D0
            DO i=0, p-1
               NdX(i,p) = (r(i+1) + l(p-i))
               tmp    = NdX(p-1,i)/NdX(i,p)
               NdX(p,i) = saved + tmp*r(i+1)
               saved  = tmp*l(p-i)
            END DO
            NdX(p,p)  = saved
         END DO

!        Load the basis functions
         ders(:,0) = NdX(bs%p,:)
!        Now compute derivatives
         DO i=0, bs%p
            s1 = 1; s2 = 2
            a(0,1) = 1D0
            DO k=1, bs%p
               d = 0D0
               ik = i - k
               pk = bs%p - k
               IF (i .ge. k) THEN
                  a(0,s2) = a(0,s1)/NdX(ik,pk+1)
                  d = a(0,s2) * NdX(pk,ik)
               END IF
               IF (ik .ge. -1) THEN
                  j1 = 1
               ELSE
                  j1 = -ik
               END IF
               IF (i-1 .le. pk) THEN
                  j2 = k-1
               ELSE
                  j2 = bs%p - i
               END IF
               DO j=j1, j2
                  a(j,s2) = (a(j,s1)-a(j-1,s1))/NdX(ik+j,pk+1)
                  d = d + a(j,s2)*NdX(pk,ik+j)
               END DO
               IF (i .le. pk) THEN
                  a(k,s2) = -a(k-1,s1)/NdX(i,pk+1)
                  d = d + a(k,s2) * NdX(pk,i)
               END IF
               ders(i,k) = d
               j  = s1
               s1 = s2
               s2 = j
            END DO
         END DO

!        Multiply with correction factors
         i = bs%p
         DO k=1, bs%p
            DO j=0, bs%p
               ders(j,k) = ders(j,k)*real(i,kind=8)*(dxi**k)
            END DO
            i = i * (bs%p-k)
         END DO

!        Copy basis functions and first two derivatives at each Gauss pt
         Nxx(:,g) = 0D0
         DO p=1, bs%p+1
            N(p,g)   = ders(p-1,0)
            Nx(p,g)  = ders(p-1,1)
            IF (bs%p .GE. 2) Nxx(p,g) = ders(p-1,2)
         END DO

         DEALLOCATE(l, r, NdX, ders, a)
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
         REAL(KIND=8), ALLOCATABLE :: Nx(:,:)
         REAL(KIND=8), ALLOCATABLE :: Nxx(:,:)
      END TYPE dType

      INTEGER insd, i, Ac, a, ax, ay, az, g, gx, gy, gz, p, nG
      REAL(KIND=8) s, sx, sy, sz, sxx, syy, szz, sxy, syz, szx, Nx, Ny,
     2   Nz
      TYPE(dType), ALLOCATABLE :: d(:)

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      ALLOCATE(d(insd))

!     First calculating shape functions of BSplines in each direction
      DO i=1, insd
         p  = lM%bs(i)%p
         nG = lM%bs(i)%nG
         ALLOCATE(d(i)%N(p+1,nG), d(i)%Nx(p+1,nG), d(i)%Nxx(p+1,nG))
         CALL BSPNNX(lM%INN(i,e), lM%bs(i), d(i)%N, d(i)%Nx, d(i)%Nxx)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (insd .EQ. 2) THEN
         DO gx=1, lM%bs(1)%nG
            DO gy=1, lM%bs(2)%nG
               g   = g + 1
               a   = 0
               s   = 0D0
               sx  = 0D0
               sy  = 0D0
               sxx = 0D0
               syy = 0D0
               sxy = 0D0
               DO ax=1, lM%bs(1)%p + 1
                  DO ay=1, lM%bs(2)%p + 1
                     a  = a + 1
                     Ac = lM%IEN(a,e)
                     Ac = lM%lN(Ac)

                     lM%N(a,g)     = d(1)%N(ax,gx)  * d(2)%N(ay,gy)
     2                             * lM%nW(Ac)
                     lM%Nx(1,a,g)  = d(1)%Nx(ax,gx) * d(2)%N(ay,gy)
     2                             * lM%nW(Ac)
                     lM%Nx(2,a,g)  = d(1)%N(ax,gx)  * d(2)%Nx(ay,gy)
     2                             * lM%nW(Ac)
                     lM%Nxx(1,a,g) = d(1)%Nxx(ax,gx)* d(2)%N(ay,gy)
     2                             * lM%nW(Ac)
                     lM%Nxx(2,a,g) = d(1)%N(ax,gx)  * d(2)%Nxx(ay,gy)
     2                             * lM%nW(Ac)
                     lM%Nxx(3,a,g) = d(1)%Nx(ax,gx) * d(2)%Nx(ay,gy)
     2                             * lM%nW(Ac)

                     s   = s   + lM%N(a,g)
                     sx  = sx  + lM%Nx(1,a,g)
                     sy  = sy  + lM%Nx(2,a,g)
                     sxx = sxx + lM%Nxx(1,a,g)
                     syy = syy + lM%Nxx(2,a,g)
                     sxy = sxy + lM%Nxx(3,a,g)
                  END DO
               END DO
               DO a=1, lM%eNoN
                  Nx = lM%Nx(1,a,g)
                  Ny = lM%Nx(2,a,g)
                  lM%N(a,g)     = lM%N(a,g)/s

                  lM%Nx(1,a,g)  = ( lM%Nx(1,a,g) - lM%N(a,g)*sx ) / s
                  lM%Nx(2,a,g)  = ( lM%Nx(2,a,g) - lM%N(a,g)*sy ) / s

                  lM%Nxx(1,a,g) = ( lM%Nxx(1,a,g) - lM%N(a,g)*sxx
     2               - 2D0*lM%Nx(1,a,g)*sx ) / s
                  lM%Nxx(2,a,g) = ( lM%Nxx(2,a,g) - lM%N(a,g)*syy
     2               - 2D0*lM%Nx(2,a,g)*sy ) / s
                  lM%Nxx(3,a,g) = ( lM%Nxx(3,a,g)*s + lM%N(a,g)*
     2               (2D0*sx*sy - sxy*s) - sx*Ny - sy*Nx ) / s**2
               END DO
               lM%w(g) = gW(gx,lM%bs(1)%nG)*gW(gy,lM%bs(2)%nG)
            END DO
         END DO
      ELSE IF (insd .EQ. 3) THEN
         DO gx=1, lM%bs(1)%nG
            DO gy=1, lM%bs(2)%nG
               DO gz=1, lM%bs(3)%nG
                  g   = g + 1
                  a   = 0
                  s   = 0D0
                  sx  = 0D0
                  sy  = 0D0
                  sz  = 0D0
                  sxx = 0D0
                  syy = 0D0
                  szz = 0D0
                  sxy = 0D0
                  syz = 0D0
                  szx = 0D0
                  DO ax=1, lM%bs(1)%p + 1
                     DO ay=1, lM%bs(2)%p + 1
                        DO az=1, lM%bs(3)%p + 1
                           a  = a + 1
                           Ac = lM%IEN(a,e)
                           Ac = lM%lN(Ac)

                           lM%N(a,g)    = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(1,a,g) = d(1)%Nx(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(2,a,g) = d(1)%N(ax,gx)*d(2)%Nx(ay,gy)
     2                                  * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nx(3,a,g) = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                                  * d(3)%Nx(az,gz)*lM%nW(Ac)
                           lM%Nxx(1,a,g) = d(1)%Nxx(ax,gx)*d(2)%N(ay,gy)
     2                                   * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nxx(2,a,g) = d(1)%N(ax,gx)*d(2)%Nxx(ay,gy)
     2                                   * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nxx(3,a,g) = d(1)%N(ax,gx)*d(2)%N(ay,gy)
     2                                   * d(3)%Nxx(az,gz)*lM%nW(Ac)
                           lM%Nxx(4,a,g) = d(1)%Nx(ax,gx)*d(2)%Nx(ay,gy)
     2                                   * d(3)%N(az,gz)*lM%nW(Ac)
                           lM%Nxx(5,a,g) = d(1)%N(ax,gx)*d(2)%Nx(ay,gy)
     2                                   * d(3)%Nx(az,gz)*lM%nW(Ac)
                           lM%Nxx(6,a,g) = d(1)%Nx(ax,gx)*d(2)%N(ay,gy)
     2                                   * d(3)%Nx(az,gz)*lM%nW(Ac)

                           s   = s   + lM%N(a,g)
                           sx  = sx  + lM%Nx(1,a,g)
                           sy  = sy  + lM%Nx(2,a,g)
                           sz  = sz  + lM%Nx(3,a,g)
                           sxx = sxx + lM%Nxx(1,a,g)
                           syy = syy + lM%Nxx(2,a,g)
                           szz = szz + lM%Nxx(3,a,g)
                           sxy = sxy + lM%Nxx(4,a,g)
                           syz = syz + lM%Nxx(5,a,g)
                           szx = szx + lM%Nxx(6,a,g)
                        END DO
                     END DO
                  END DO
                  DO a=1, lM%eNoN
                     Nx = lM%Nx(1,a,g)
                     Ny = lM%Nx(2,a,g)
                     Nz = lM%Nx(3,a,g)
                     lM%N(a,g)     = lM%N(a,g)/s

                     lM%Nx(1,a,g)  = ( lM%Nx(1,a,g) - lM%N(a,g)*sx ) / s
                     lM%Nx(2,a,g)  = ( lM%Nx(2,a,g) - lM%N(a,g)*sy ) / s
                     lM%Nx(3,a,g)  = ( lM%Nx(3,a,g) - lM%N(a,g)*sz ) / s

                     lM%Nxx(1,a,g) = ( lM%Nxx(1,a,g) - lM%N(a,g)*sxx
     2                  - 2D0*lM%Nx(1,a,g)*sx ) / s
                     lM%Nxx(2,a,g) = ( lM%Nxx(2,a,g) - lM%N(a,g)*syy
     2                  - 2D0*lM%Nx(2,a,g)*sy ) / s
                     lM%Nxx(3,a,g) = ( lM%Nxx(3,a,g) - lM%N(a,g)*szz
     2                  - 2D0*lM%Nx(3,a,g)*sz ) / s
                     lM%Nxx(4,a,g) = ( lM%Nxx(4,a,g)*s + lM%N(a,g)*
     2                  (2D0*sx*sy - sxy*s) - sx*Ny - sy*Nx ) / s**2
                     lM%Nxx(5,a,g) = ( lM%Nxx(5,a,g)*s + lM%N(a,g)*
     2                  (2D0*sy*sz - syz*s) - sy*Nz - sz*Ny ) / s**2
                     lM%Nxx(6,a,g) = ( lM%Nxx(6,a,g)*s + lM%N(a,g)*
     2                  (2D0*sz*sx - szx*s) - sz*Nx - sx*Nz ) / s**2
                  END DO
                  lM%w(g) = gW(gx,lM%bs(1)%nG)*gW(gy,lM%bs(2)%nG)
     2                    * gW(gz,lM%bs(3)%nG)
               END DO
            END DO
         END DO
      END IF

      DEALLOCATE(d)
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
         REAL(KIND=8), ALLOCATABLE :: Nx(:,:)
         REAL(KIND=8), ALLOCATABLE :: Nxx(:,:)
      END TYPE dType

      INTEGER insd, i, j, Ac, a, ax, ay, g, gx, gy, Ec, p, nG
      REAL(KIND=8) s, sx, sy, sxx, syy, sxy, Nx, Ny
      TYPE(dType), ALLOCATABLE :: d(:)

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      ALLOCATE(d(insd))

      Ec = lFa%gE(e)
!     First calculating shape functions of BSplines in each direction,
!     exept the direction that face is normal to
      DO i=1, insd
         p  = lM%bs(i)%p
         nG = lM%bs(i)%nG
         ALLOCATE(d(i)%N(p+1,nG), d(i)%Nx(p+1,nG), d(i)%Nxx(p+1,nG))
         CALL BSPNNX(lM%INN(i,Ec), lM%bs(i), d(i)%N, d(i)%Nx, d(i)%Nxx)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (insd .EQ. 2) THEN
         IF (lFa%d .EQ. 1) i = 2
         IF (lFa%d .EQ. 2) i = 1
         DO g=1, lM%bs(i)%nG
            s   = 0D0
            sx  = 0D0
            sxx = 0D0
            DO a=1, lM%bs(i)%p + 1
               Ac = lFa%IEN(a,e)
               Ac = lM%lN(Ac)

               lFa%N(a,g)     = d(i)%N(a,g)*lM%nW(Ac)
               lFa%Nx(1,a,g)  = d(i)%Nx(a,g)*lM%nW(Ac)
               lFa%Nxx(1,a,g) = d(i)%Nxx(a,g)*lM%nW(Ac)

               s   = s   + lFa%N(a,g)
               sx  = sx  + lFa%Nx(1,a,g)
               sxx = sxx + lFa%Nxx(1,a,g)
            END DO
            DO a=1, lFa%eNoN
               lFa%N(a,g)     = lFa%N(a,g)/s
               lFa%Nx(1,a,g)  = ( lFa%Nx(1,a,g) - lFa%N(a,g)*sx ) / s
               lFa%Nxx(1,a,g) = ( lFa%Nxx(1,a,g) - lFa%N(a,g)*sxx
     2            - 2D0*lFa%Nx(1,a,g)*sx ) / s
            END DO
            lFa%w(g) = gW(g,lM%bs(i)%nG)
         END DO
      ELSE IF (insd .EQ. 3) THEN
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
               g   = g + 1
               a   = 0
               s   = 0D0
               sx  = 0D0
               sy  = 0D0
               sxx = 0D0
               syy = 0D0
               sxy = 0D0
               DO ax=1, lM%bs(i)%p + 1
                  DO ay=1, lM%bs(j)%p + 1
                     a  = a + 1
                     Ac = lFa%IEN(a,e)
                     Ac = lM%lN(Ac)

                     lFa%N(a,g)     = d(i)%N(ax,gx)  * d(j)%N(ay,gy)
     2                              * lM%nW(Ac)
                     lFa%Nx(1,a,g)  = d(i)%Nx(ax,gx) * d(j)%N(ay,gy)
     2                              * lM%nW(Ac)
                     lFa%Nx(2,a,g)  = d(i)%N(ax,gx)  * d(j)%Nx(ay,gy)
     2                              * lM%nW(Ac)
                     lFa%Nxx(1,a,g) = d(i)%Nxx(ax,gx)* d(j)%N(ay,gy)
     2                              * lM%nW(Ac)
                     lFa%Nxx(2,a,g) = d(i)%N(ax,gx)  * d(j)%Nxx(ay,gy)
     2                              * lM%nW(Ac)
                     lFa%Nxx(3,a,g) = d(i)%Nx(ax,gx) * d(j)%Nx(ay,gy)
     2                              * lM%nW(Ac)

                     s   = s   + lFa%N(a,g)
                     sx  = sx  + lFa%Nx(1,a,g)
                     sy  = sy  + lFa%Nx(2,a,g)
                     sxx = sxx + lFa%Nxx(1,a,g)
                     syy = syy + lFa%Nxx(2,a,g)
                     sxy = sxy + lFa%Nxx(3,a,g)
                  END DO
               END DO
               DO a=1, lFa%eNoN
                  Nx = lFa%Nx(1,a,g)
                  Ny = lFa%Nx(2,a,g)
                  lFa%N(a,g)    = lFa%N(a,g)/s

                  lFa%Nx(1,a,g) = ( lFa%Nx(1,a,g) - lFa%N(a,g)*sx ) / s
                  lFa%Nx(2,a,g) = ( lFa%Nx(2,a,g) - lFa%N(a,g)*sy ) / s

                  lFa%Nxx(1,a,g) = ( lFa%Nxx(1,a,g) - lFa%N(a,g)*sxx
     2               - 2D0*lFa%Nx(1,a,g)*sx ) / s
                  lFa%Nxx(2,a,g) = ( lFa%Nxx(2,a,g) - lFa%N(a,g)*syy
     2               - 2D0*lFa%Nx(2,a,g)*sy ) / s
                  lFa%Nxx(3,a,g) = ( lFa%Nxx(3,a,g)*s + lFa%N(a,g)*
     2               (2D0*sx*sy - sxy*s) - sx*Ny - sy*Nx ) / s**2
               END DO
               lFa%w(g) = gW(gx,lM%bs(i)%nG)*gW(gy,lM%bs(j)%nG)
            END DO
         END DO
      END IF

      DEALLOCATE(d)
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

      INTEGER insd, i, Ac, a, ax, ay, az, g, gx, gy, gz
      REAL(KIND=8) s
      TYPE(dType), ALLOCATABLE :: d(:)

      insd = nsd
      IF(lM%lShl) insd = nsd - 1
      ALLOCATE(d(insd))

!     First calculating shape functions of BSplines in each direction
      DO i=1, insd
         ALLOCATE(d(i)%N(lM%bs(i)%p+1,lM%bs(i)%nSl))
         CALL BSPNNS(lM%INN(i,e), lM%bs(i), d(i)%N)
      END DO

!     Combining BSpline shape functions
      g = 0
      IF (insd .EQ. 2) THEN
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
      ELSE IF (insd .EQ. 3) THEN
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

      DEALLOCATE(d)
      RETURN
      END SUBROUTINE NRBNNS
!####################################################################

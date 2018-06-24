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
!     These set of routines perform Delaunay Triangulation for any
!     given 2D point cloud
!
!--------------------------------------------------------------------

      SUBROUTINE DELAUTRI2D(nNo, coords, nel, ien, neigh)
      USE COMMOD, ONLY : err
      USE UTILMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNo
      REAL(KIND=8), INTENT(INOUT) :: coords(2,nNo)
      INTEGER, INTENT(OUT) :: nel
      INTEGER, INTENT(OUT), dimension(3,3*nNo) :: ien,neigh
      INTEGER :: ierr,e,i,j,k,l,lr,lrline,m,m1,m2,n
      INTEGER :: ledg,ltri,redg,rtri,t,top
      INTEGER :: indx(nNo),stack(nNo)
      REAL(KIND=8) :: cmax,tol

      tol = 100.0D+00 * EPSILON(tol)
      ierr = 0

      CALL SORT_HEAP_INDEX(nNo, coords, indx)

      CALL PERMUTE(nNo, indx, coords)

      m1 = 1
      DO i = 2, nNo
         m = m1
         m1 = i
         k = 0

         DO j = 1, 2
            cmax = MAX(ABS(coords(j,m)), ABS(coords(j,m1)))
            IF (tol*(cmax + 1D0) .LT.
     2          ABS(coords(j,m) - coords(j,m1))) THEN
               k = j
               EXIT
            END IF
         END DO

         IF (k .EQ. 0) THEN
            err = "DELAU2D Fatal: Fails for point I="//STR(i)
         END IF
      END DO

      m1 = 1
      m2 = 2
      j = 3
      DO
         IF (nNo .LT. j) THEN
            err = "DELAU2D Fatal error!"
         END IF
         m = j
         lr = LRLINE(coords(1,m ), coords(2,m ),
     2               coords(1,m1), coords(2,m1),
     3               coords(1,m2), coords(2,m2), 0D0)
         IF (lr .NE. 0) EXIT
         j = j + 1
      END DO

      nel = j - 2
      IF (lr .EQ. -1) THEN
         ien(1,1) = m1
         ien(2,1) = m2
         ien(3,1) = m
         neigh(3,1) = -3

         DO i = 2, nel
            m1 = m2
            m2 = i+1
            ien(1,i) = m1
            ien(2,i) = m2
            ien(3,i) = m
            neigh(1,i-1) = -3 * i
            neigh(2,i-1) = i
            neigh(3,i) = i - 1
         END DO

         neigh(1,nel) = -3 * nel - 1
         neigh(2,nel) = -5
         ledg = 2
         ltri = nel
      ELSE
         ien(1,1) = m2
         ien(2,1) = m1
         ien(3,1) = m
         neigh(1,1) = -4

         DO i = 2, nel
            m1 = m2
            m2 = i+1
            ien(1,i) = m2
            ien(2,i) = m1
            ien(3,i) = m
            neigh(3,i-1) = i
            neigh(1,i) = -3 * i - 3
            neigh(2,i) = i - 1
         END DO

         neigh(3,nel) = -3 * nel
         neigh(2,1) = -3 * nel - 2
         ledg = 2
         ltri = 1
      END IF

      top = 0
      DO i = j+1, nNo
         m = i
         m1 = ien(ledg,ltri)
         IF (ledg .LE. 2) THEN
            m2 = ien(ledg+1,ltri)
         ELSE
            m2 = ien(1,ltri)
         END IF

         lr = LRLINE(coords(1,m ), coords(2,m ),
     2               coords(1,m1), coords(2,m1),
     3               coords(1,m2), coords(2,m2), 0D0)
         IF (lr .GT. 0) THEN
            rtri = ltri
            redg = ledg
            ltri = 0
         ELSE
            l = -neigh(ledg,ltri)
            rtri = l / 3
            redg = MOD(l,3) + 1
         END IF

         CALL VBEDG(coords(1,m), coords(2,m), nNo, coords, nel, ien,
     2              neigh, ltri, ledg, rtri, redg )
         n = nel + 1
         l = -neigh(ledg,ltri)

         DO
            t = l / 3
            e = MOD(l, 3) + 1
            l = -neigh(e,t)
            m2 = ien(e,t)

            IF (e .LE. 2) THEN
               m1 = ien(e+1,t)
            ELSE
               m1 = ien(1,t)
            END IF

            nel = nel + 1
            neigh(e,t) = nel
            ien(1,nel) = m1
            ien(2,nel) = m2
            ien(3,nel) = m
            neigh(1,nel) = t
            neigh(2,nel) = nel - 1
            neigh(3,nel) = nel + 1
            top = top + 1

            IF (nNo .LT. top) THEN
               err = "DELAU2D Fatal: Stack overflow"
            END IF

            stack(top) = nel
            IF ((t .EQ. rtri) .AND. (e .EQ. redg)) EXIT
         END DO

         neigh(ledg,ltri) = -3 * n - 1
         neigh(2,n) = -3 * nel - 2
         neigh(3,nel) = -l
         ltri = n
         ledg = 2

         CALL SWAPEC(m, top, ltri, ledg, nNo, coords,
     2               nel, ien, neigh, stack, ierr)

         IF (ierr .NE. 0) THEN
            err = "DELAU2D Fatal: Error return from SWAPEC"
         END IF

      END DO

      DO i = 1, 3
         DO j = 1, nel
            ien(i,j) = indx(ien(i,j))
         END DO
      END DO

      CALL PERM_INVERSE(nNo, indx)

      CALL PERMUTE(nNo, indx, coords)

      END SUBROUTINE DELAUTRI2D

!**************************************************

      SUBROUTINE SORT_HEAP_INDEX(n, a, indx)
      IMPLICIT NONE
      INTEGER, PARAMETER :: ndim=2
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(OUT) :: indx(n)
      REAL(KIND=8), INTENT(IN) :: a(ndim,n)
      INTEGER :: i,ir,j,l,indxt
      REAL(KIND=8) :: aval(ndim)

      DO i=1, n
         indx(i) = i
      END DO

      l = n/2 + 1
      ir = n

      DO
         IF (l .GT. 1) THEN
            l = l - 1
            indxt = indx(l)
            aval(1:ndim) = a(1:ndim,indxt)
         ELSE
            indxt = indx(ir)
            aval(1:ndim) = a(1:ndim,indxt)
            indx(ir) = indx(1)
            ir = ir - 1
            IF (ir .EQ. 1) THEN
               indx(1) = indxt
               EXIT
            END IF
         END IF

         i = l
         j = l + l
         DO WHILE (j .LE. ir)
            IF (j .LT. ir) THEN
               IF ( (  a(1,indx(j)) .LT. a(1,indx(j+1))  ) .OR.
     2              ( (a(1,indx(j)) .EQ. a(1,indx(j+1))) .AND.
     3                (a(2,indx(j)) .LT. a(2,indx(j+1))) ) ) THEN
                  j = j + 1
               END IF
            END IF

            IF ( (  aval(1) .LT. a(1,indx(j))  ) .OR.
     2           ( (aval(1) .EQ. a(1,indx(j))) .AND.
     3             (aval(2) .LT. a(2,indx(j))) ) ) THEN
               indx(i) = indx(j)
               i = j
               j = j + j
            ELSE
               j = ir + 1
            END IF
         END DO

         indx(i) = indxt
      END DO

      END SUBROUTINE SORT_HEAP_INDEX

!**************************************************

      SUBROUTINE PERMUTE(n, p, a)
      USE COMMOD, ONLY : err
      USE UTILMOD
      IMPLICIT NONE
      INTEGER, PARAMETER :: ndim=2, base=1
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: p(n)
      REAL(KIND=8), INTENT(INOUT) :: a(ndim,n)

      INTEGER :: ierror, iget, iput, istart
      REAL(KIND=8) :: a_temp(ndim)

      CALL PERM_CHECK(n, p, base, ierror)

      IF (ierror .NE. 0) THEN
         err = "DELAU2D R82VEC_Permute Fatal error"
      END IF

      DO istart = 1, n
         IF (p(istart) .LT. 0) THEN
            CYCLE
         ELSE IF (p(istart) .EQ. istart) THEN
            p(istart) = - p(istart)
            CYCLE
         ELSE
            a_temp(1:ndim) = a(1:ndim,istart)
            iget = istart
            DO
               iput = iget
               iget = p(iget)
               p(iput) = - p(iput)
               IF ((iget .LT. 1) .OR. (n .LT. iget)) THEN
                  err = "DELAU2D R82VEC_Permute Fatal: "//
     2                   "Permutation index out of range"
               END IF

               IF (iget .EQ. istart) THEN
                  a(1:ndim,iput) = a_temp(1:ndim)
                  EXIT
               END IF

               a(1:ndim,iput) = a(1:ndim,iget)
            END DO
         END IF
      END DO
      p(1:n) = -p(1:n)

      END SUBROUTINE PERMUTE

!**************************************************

      SUBROUTINE PERM_CHECK(n, p, base, ierror)
      USE COMMOD, ONLY : err
      USE UTILMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n,base
      INTEGER, INTENT(IN) :: p(n)
      INTEGER, INTENT(OUT) :: ierror
      INTEGER :: ifind,iseek

      ierror = 0
      DO iseek = base, base+n-1
         ierror = 1
         DO ifind = 1, n
            IF (p(ifind) .EQ. iseek) THEN
               ierror = 0
               EXIT
            END IF
         END DO
         IF (ierror .NE. 0) THEN
            err = "DELAU2D PERM_CHECK Fatal"
         END IF
      END DO

      END SUBROUTINE PERM_CHECK

!**************************************************

      SUBROUTINE PERM_INVERSE(n, p)
      USE COMMOD, ONLY : err
      USE UTILMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: p(n)
      INTEGER, parameter :: base = 1
      INTEGER :: i,i0,i1,i2,ierror,is
      INTEGER, EXTERNAL :: I4_SIGN

      IF (n .LE. 0) THEN
         err = "DELAU2D PERM_INVERSE Fatal: input n <= 0"
      END IF

      CALL PERM_CHECK(n, p, base, ierror)

      IF (ierror .NE. 0) THEN
         err = "DELAU2D PERM_INVERSE Fatal"
      END IF

      is = 1
      DO i = 1, n
         i1 = p(i)
         DO WHILE (i .LT. i1)
            i2 = p(i1)
            p(i1) = -i2
            i1 = i2
         END DO
         is = -I4_SIGN(p(i))
         p(i) = is * ABS(p(i))
      END DO

      DO i = 1, n
         i1 = - p(i)
         IF (i1 .GE. 0) THEN
            i0 = i
            DO
               i2 = p(i1)
               p(i1) = i0
               IF (i2 .LT. 0) EXIT
               i0 = i1
               i1 = i2
            END DO
         END IF
      END DO

      END SUBROUTINE PERM_INVERSE

!**************************************************

      INTEGER FUNCTION LRLINE(xu, yu, xv1, yv1, xv2, yv2, dv)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: xu,yu,xv1,yv1,xv2,yv2,dv
      REAL(KIND=8) :: dx,dxu,dy,dyu,t,tol,tolABS

      tol = 1D2 * EPSILON (tol)
      dx = xv2 - xv1
      dy = yv2 - yv1
      dxu = xu - xv1
      dyu = yu - yv1
      tolABS = tol * max(ABS(dx), ABS(dy), ABS(dxu), ABS(dyu), ABS(dv))
      t = dy*dxu - dx*dyu + dv*dsqrt( dx**2 + dy**2 )

      IF (tolABS .LT. t) THEN
         LRLINE = 1
      ELSE IF (-tolABS .LE. t) THEN
         LRLINE = 0
      ELSE
         LRLINE = -1
      END IF

      END FUNCTION LRLINE

!**************************************************

      SUBROUTINE VBEDG(x, y, nNo, coords, nel, ien,
     2                 neigh, ltri, ledg, rtri, redg )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nNo,nel
      INTEGER, INTENT(IN) :: ien(3,nel),neigh(3,nel)
      INTEGER, INTENT(INOUT) :: ltri,ledg,rtri,redg
      REAL(KIND=8), INTENT(IN) :: x,y
      REAL(KIND=8), INTENT(IN) :: coords(2,nNo)

      INTEGER :: a,b,e,l,lr,t
      INTEGER :: LRLINE
      INTEGER, EXTERNAL :: I4_WRAP
      LOGICAL :: ldone

      IF (ltri .EQ. 0) THEN
         ldone = .FALSE.
         ltri = rtri
         ledg = redg
      ELSE
         ldone = .TRUE.
      END IF

      DO
         l = -neigh(redg,rtri)
         t = l / 3
         e = mod(l, 3) + 1
         a = ien(e,t)
         IF (e .LE. 2) THEN
            b = ien(e+1,t)
         ELSE
            b = ien(1,t)
         END IF
         lr = LRLINE(x, y, coords(1,a), coords(2,a),
     2               coords(1,b), coords(2,b), 0D0)
         IF (lr .LE. 0) EXIT
         rtri = t
         redg = e
      END DO

      IF (ldone) RETURN

      t = ltri
      e = ledg

      DO
         b = ien(e,t)
         e = I4_WRAP(e-1, 1, 3)
         DO WHILE (neigh(e,t) .GT. 0)
            t = neigh(e,t)
            IF (ien(1,t) .EQ. b) THEN
               e = 3
            ELSE IF (ien(2,t) .EQ. b) THEN
               e = 1
            ELSE
               e = 2
            END IF
         END DO
         a = ien(e,t)
         lr = LRLINE(x, y, coords(1,a), coords(2,a),
     2               coords(1,b), coords(2,b), 0D0)

         IF (lr .LE. 0) EXIT
      END DO

      ltri = t
      ledg = e

      END SUBROUTINE VBEDG

!**************************************************

      SUBROUTINE SWAPEC(i, top, btri, bedg, nNo, coords, nel, ien,
     2                  neigh, stack, ierr)
      USE COMMOD, ONLY : err
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,nNo,nel
      INTEGER, INTENT(OUT) :: ierr
      INTEGER, INTENT(INOUT) :: top,btri,bedg
      INTEGER, INTENT(INOUT) :: ien(3,nel),neigh(3,nel)
      REAL(KIND=8), INTENT(IN) :: coords(2,nNo)
      INTEGER :: stack(nNo)

      INTEGER :: a,b,c,e,ee,em1,ep1,f,fm1,fp1
      INTEGER :: l,r,s,swap,t,tt,u
      INTEGER :: DIAEDG
      INTEGER, EXTERNAL :: I4_WRAP
      REAL(KIND=8) :: x,y

      ierr = 0
      x = coords(1,i)
      y = coords(2,i)
      DO
         IF (top .LE. 0) EXIT
         t = stack(top)
         top = top - 1
         IF (ien(1,t) .EQ. i) THEN
            e = 2
            b = ien(3,t)
         ELSE IF (ien(2,t) .EQ. i) THEN
            e = 3
            b = ien(1,t)
         ELSE
            e = 1
            b = ien(2,t)
         END IF

         a = ien(e,t)
         u = neigh(e,t)

         IF (neigh(1,u) .EQ. t) THEN
            f = 1
            c = ien(3,u)
         ELSE IF (neigh(2,u) .EQ. t) THEN
            f = 2
            c = ien(1,u)
         ELSE
            f = 3
            c = ien(2,u)
         END IF

         swap = DIAEDG(x, y, coords(1,a), coords(2,a),
     2                       coords(1,c), coords(2,c),
     3                       coords(1,b), coords(2,b) )

         IF (swap .EQ. 1) THEN
            em1 = I4_WRAP(e-1, 1, 3)
            ep1 = I4_WRAP(e+1, 1, 3)
            fm1 = I4_WRAP(f-1, 1, 3)
            fp1 = I4_WRAP(f+1, 1, 3)

            ien(ep1,t) = c
            ien(fp1,u) = i
            r = neigh(ep1,t)
            s = neigh(fp1,u)
            neigh(ep1,t) = u
            neigh(fp1,u) = t
            neigh(e,t) = s
            neigh(f,u) = r

            IF (0 .LT. neigh(fm1,u)) THEN
               top = top + 1
               stack(top) = u
            END IF

            IF (s .GT. 0) THEN
               IF (neigh(1,s) .EQ. u) THEN
                  neigh(1,s) = t
               ELSE IF (neigh(2,s) .EQ. u) THEN
                  neigh(2,s) = t
               ELSE
                  neigh(3,s) = t
               END IF
               top = top + 1

               IF (nNo .LT. top) THEN
                  err = "DELAU2D SWAPEC Fatal"
               END IF
               stack(top) = t
            ELSE
               IF ((u .EQ. btri) .AND. (fp1 .EQ. bedg)) THEN
                  btri = t
                  bedg = e
               END IF
               l = - (3*t + e - 1)
               tt = t
               ee = em1

               DO WHILE (neigh(ee,tt) .GT. 0)
                  tt = neigh(ee,tt)
                  IF (ien(1,tt) .EQ. a) THEN
                     ee = 3
                  ELSE IF (ien(2,tt) .EQ. a) THEN
                     ee = 1
                  ELSE
                     ee = 2
                  END IF
               END DO
               neigh(ee,tt) = l
            END IF

            IF (r .GT. 0) THEN
               IF (neigh(1,r) .EQ. t) THEN
                  neigh(1,r) = u
               ELSE IF (neigh(2,r) .EQ. t) THEN
                  neigh(2,r) = u
               ELSE
                  neigh(3,r) = u
               END IF
            ELSE
               IF ((t .EQ. btri) .AND. (ep1 .EQ. bedg)) THEN
                  btri = u
                  bedg = f
               END IF
               l = - (3*u + f - 1)
               tt = u
               ee = fm1
               DO WHILE (neigh(ee,tt) .GT. 0)
                  tt = neigh(ee,tt)
                  IF (ien(1,tt) .EQ. b) THEN
                     ee = 3
                  ELSE IF (ien(2,tt) .EQ. b) THEN
                     ee = 1
                  ELSE
                     ee = 2
                  END IF
               END DO
               neigh(ee,tt) = l
            END IF

         END IF ! swap

      END DO

      END SUBROUTINE SWAPEC

!**************************************************

      INTEGER FUNCTION DIAEDG(x0, y0, x1, y1, x2, y2, x3, y3)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: x0, y0, x1, y1, x2, y2, x3, y3

      REAL(KIND=8) :: ca,cb,dx10,dx12,dx30,dx32
      REAL(KIND=8) :: dy10,dy12,dy30,dy32
      REAL(KIND=8) :: s,tol,tola,tolb

      tol = 1D2 * EPSILON(tol)

      dx10 = x1 - x0
      dy10 = y1 - y0
      dx12 = x1 - x2
      dy12 = y1 - y2
      dx30 = x3 - x0
      dy30 = y3 - y0
      dx32 = x3 - x2
      dy32 = y3 - y2

      tola = tol * MAX(ABS(dx10), ABS(dy10), ABS(dx30), ABS(dy30))
      tolb = tol * MAX(ABS(dx12), ABS(dy12), ABS(dx32), ABS(dy32))

      ca = dx10 * dx30 + dy10 * dy30
      cb = dx12 * dx32 + dy12 * dy32

      IF ((tola .LT. ca) .AND. (tolb .LT. cb)) THEN
         DIAEDG = -1
      ELSE IF ((ca .LT. -tola) .AND. (cb .LT. -tolb)) THEN
         DIAEDG = 1
      ELSE
         tola = MAX(tola, tolb)
         s = (dx10*dy30 - dx30*dy10) * cb +
     2       (dx32*dy12 - dx12*dy32) * ca
         IF (tola .LT. s) THEN
            DIAEDG = -1
         ELSE IF (s .LT. -tola) THEN
            DIAEDG = 1
         ELSE
            DIAEDG = 0
         END IF
      END IF

      END FUNCTION DIAEDG

!**************************************************

      INTEGER FUNCTION I4_SIGN(x)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: x

      IF (x .LT. 0) THEN
         I4_SIGN = -1
      ELSE
         I4_SIGN = +1
      END IF

      END FUNCTION I4_SIGN

!**************************************************

      INTEGER FUNCTION I4_WRAP(ival, ilo, ihi)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ival,ilo,ihi
      INTEGER :: jhi,jlo,value,wide
      INTEGER, EXTERNAL :: I4_MODP

      jlo = MIN(ilo, ihi)
      jhi = MAX(ilo, ihi)
      wide = jhi - jlo + 1

      IF (wide .EQ. 1) THEN
         value = jlo
      ELSE
         value = jlo + I4_MODP(ival-jlo, wide)
      END IF

      I4_WRAP = value

      END FUNCTION I4_WRAP

!**************************************************

      INTEGER FUNCTION I4_MODP(i, j)
      USE COMMOD, ONLY : err
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i,j
      INTEGER :: value

      IF (j .EQ. 0) THEN
         err = "DELAU2D I4_MODP Fatal: divide by 0"
      END IF

      value = MOD(i, j)

      IF (value .LT. 0) THEN
         value = value + ABS (j)
      END IF
      I4_MODP = value

      END FUNCTION I4_MODP

!**************************************************

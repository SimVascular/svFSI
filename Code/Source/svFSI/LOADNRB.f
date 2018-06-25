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
!     This routine is desinged to read NURBS patches and if neccessary
!     refine it by knot insertion
!      
!--------------------------------------------------------------------
!     This is for reading a single NURBS patch
      SUBROUTINE READNRB(list, lM)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      
      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list
      TYPE(mshType), INTENT(INOUT) :: lM

      LOGICAL flag, isZ, dirSet(3)
      INTEGER n, i, j, rep, fid, nDir
      REAL(KIND=8) tmp, tmpO
      TYPE(listType), POINTER :: lPtr, lPD
      TYPE(fileType) ftmp

      lPtr => list%get(ftmp,"NURBS data file path")
      IF (.NOT.ASSOCIATED(lPtr)) RETURN
      dbg = "Reading a NURBS patch"
      fid = ftmp%open()

      lM%nNo   = 1
      lM%nFa   = 2*nsd
      lM%eType = eType_NRB
      ALLOCATE(lM%bs(nsd), lM%fa(lM%nFa))
!     These are the default faces      
      lM%fa(1)%name = "XN_"//TRIM(lM%name)
      lM%fa(2)%name = "XP_"//TRIM(lM%name)
      lM%fa(3)%name = "YN_"//TRIM(lM%name)
      lM%fa(4)%name = "YP_"//TRIM(lM%name)
      IF (nsd .EQ. 3) THEN
         lM%fa(5)%name = "ZN_"//TRIM(lM%name)
         lM%fa(6)%name = "ZP_"//TRIM(lM%name)
      END IF
      DO i=1, nsd
         lM%fa(2*i-1)%d = i
         lM%fa(2*i)%d   = i

         CALL GTBLK(fid,'knotV',n)
         flag = .TRUE.
         lM%bs(i)%n   = n
         lM%bs(i)%p   = 0
         lM%bs(i)%nEl = 0
         ALLOCATE(lM%bs(i)%xi(n))
         DO j=1, n
            READ(fid,*) tmp
            IF (j .GT. 1) THEN 
               isZ = ISZERO(tmpO,tmp)
               IF (tmpO .GT. tmp) err = "Knot vector must be"
     2            //" non-decreasing"
               IF (.NOT.isZ) lM%bs(i)%nEl = lM%bs(i)%nEl + 1
               IF (flag) THEN
                  IF (isZ) THEN
                     lM%bs(i)%p = lM%bs(i)%p + 1
                  ELSE
                     IF (j .EQ. 2) err = "Violating p > 0"
                     flag = .FALSE.
                  END IF
               END IF
               IF (j .GT. n-lM%bs(i)%p) THEN
                  IF (.NOT.isZ) err = "Expecting p+1 similar"//
     2                " knots at the end of knot vector"
               END IF
            END IF
            lM%bs(i)%xi(j) = tmp
            tmpO = tmp
         END DO
         lM%bs(i)%nNo = lM%bs(i)%n - lM%bs(i)%p - 1
         lM%nNo       = lM%nNo*lM%bs(i)%nNo
!     Default values for number of guass points and sample points         
         lM%bs(i)%nG  = MIN(lM%bs(i)%p+1,5)
         lM%bs(i)%nSl = lM%bs(i)%p + 1
      END DO
      CALL GTBLK(fid,'ctrlPts',n)
      IF (n .NE. lM%nNo) err = "Unexpected number of "//
     2   " control points"
      
      ALLOCATE(lM%nW(n), x(nsd,n))
      DO i=1, n
         READ(fid,*) x(:,i), lM%nW(i)
      END DO
      CLOSE(fid)
         
      dirSet = .FALSE.
      nDir   = list%srch("Set direction")
      DO j=1, nDir
!     This would be the direction of refinement      
         lPD => list%get(i,"Set direction",j,ul=3,ll=1)
         IF (dirSet(i)) err = TRIM(lPD%ping("Set direction"))//
     2      " is repeated more than once"
         dirSet(i) = .TRUE.

!     This would be the total number of inserted knots in each element
         n = 0
         lPtr => lPD%get(n,"Number of knot insertion",ll=0)
         IF (n .NE. 0) THEN
!     This would be the times each knot is repeated
            rep = 1
            lPtr => lPD%get(rep,"Inserted knots repetition",
     2         ll=1,ul=lM%bs(i)%p)

            IF (MOD(n,rep) .NE. 0) err = TRIM(lPD%ping("Inserted"//
     2         "knots repetition",lPtr))//" Number of repetition is"//
     3         " not compatible with the total number of inserted knots"
            CALL KNOTINS(lM, lM%bs(i), i, n, rep)
         END IF
         lPtr => lPD%get(lM%bs(i)%nG,"Number of Gauss points",ll=1,ul=5)
!     Number of sample points per element for display
         lPtr => lPD%get(lM%bs(i)%nSl,"Number of sample points",ll=2)
         lPtr => lPD%get(lM%fa(2*i-1)%name,"Start face name")
         lPtr => lPD%get(lM%fa(2*i)%name,"End face name")
      END DO
      CALL CONSTNRB(lM)

      dbg = " Number of DOF/element: "//lM%eNoN
      dbg = " Number of Gauss points/element: "//lM%nG

      RETURN
      END SUBROUTINE READNRB

!####################################################################
!     Inserts a knot into a NURBS. "i" is the direction, "n" is total
!     number of inserted knts, and "r" is the repetition of each knot
      SUBROUTINE KNOTINS(lM, bs, dir, n, rep)
      
      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(bsType), INTENT(INOUT) :: bs
      INTEGER, INTENT(IN) :: dir, n, rep
   
      INTEGER l, i, j, k, ns, tn, p, ir, m
      TYPE(bsType) tbs
      REAL(KIND=8) u, tmp, tmpO, c, alpha
      REAL(KIND=8), ALLOCATABLE :: xo(:,:), xn(:,:)

!     For convinience
      p = bs%p
!     Number of segments in each element      
      ns = n/rep + 1
!     The total number of inserted knots
      tn = n*bs%nEl
!     Using tbs, x, w as temporary containers
      ALLOCATE(tbs%xi(bs%n), xo(nsd+1,lM%nNo))
      tbs = bs

      DO i=1, lM%nNo
         xo(1:nsd,i) = x(:,i)
         xo(nsd+1,i) = lM%nW(i)
      END DO
      DEALLOCATE(bs%xi, x, lM%nW)
!     For one knot insertion
      bs%n   = tbs%n + tn
      bs%nNo = bs%nNo + tn
      bs%nEl = bs%nEl*ns
      lM%nNo = lM%nNo*bs%nNo/tbs%nNo
      ALLOCATE(bs%xi(bs%n), x(nsd,lM%nNo), lM%nW(lM%nNo),
     2   xn(nsd+1,lM%nNo))

      DO i=1, p+1
         bs%xi(i) = tbs%xi(i)
         CALL SETROW(i,i)
      END DO
      DO i=tbs%n-p,tbs%n
         j = i + tn
         bs%xi(j) = tbs%xi(i)
      END DO

      j = p
      DO i=p+1, tbs%n - p
         tmp = tbs%xi(i)
         IF (i .EQ. p+1) THEN
            tmpO     = tmp
            j        = j + 1
            bs%xi(j) = tmp
            CYCLE
         END IF
         IF (ISZERO(tmpO,tmp)) THEN
            j = j + 1
            bs%xi(j) = tmp
            CALL SETROW(j-1,i-1)
            CYCLE
         END IF
!     Found an element. Now I will break it down to ns pieces
         c = (tmp - tmpO)/REAL(ns,8)
         DO k=1, ns
            DO ir=1, rep
               IF (k.EQ.ns .AND. ir.GT.1) EXIT
               j = j + 1
               u = tmpO + c*REAL(k,8)
               CALL SETROW(j-1,i-1)
               bs%xi(j) = u
               IF (k .EQ. ns) EXIT
               DO l=p, 1, -1
                  m     = l + j - p - 1
                  alpha = (u - bs%xi(m))/(tbs%xi(l+i-1) - bs%xi(m))
                  CALL INTROW(m, alpha)
               END DO
            END DO
         END DO
         tmpO = tmp
      END DO
      DO i=1, lM%nNo
         x(:,i)   = xn(1:nsd,i)
         lM%nW(i) = xn(nsd+1,i)
      END DO

      RETURN
      CONTAINS 
!--------------------------------------------------------------------      
!     This routine sets %x and %nW
      SUBROUTINE INTROW(a1, alpha)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: a1
      REAL(KIND=8), INTENT(IN) :: alpha

      INTEGER i, j, k, c, p, jmp(nsd), a2, a3
      REAL(KIND=8) w1, w2

      IF (ISZERO(alpha,1D0)) RETURN

      i = dir
      IF (nsd .EQ. 2) THEN
         jmp(1) = lM%bs(2)%nNo
         jmp(2) = 1
         j      = 1
         IF (i .EQ. 1) j = 2

         DO a2=1, lM%bs(j)%nNo
            c = jmp(i)*(a1-1) + jmp(j)*(a2-1) + 1
            p = c - jmp(i)
 
            w1      = xn(3,c)*alpha
            w2      = xn(3,p)*(1D0-alpha)
            xn(3,c) = w1 + w2
            w1      = w1/xn(3,c)
            w2      = w2/xn(3,c)
            xn(1,c) = w1*xn(1,c) + w2*xn(1,p)
            xn(2,c) = w1*xn(2,c) + w2*xn(2,p)
         END DO
      ELSE
         jmp(1) = lM%bs(3)%nNo*lM%bs(2)%nNo
         jmp(2) = lM%bs(3)%nNo
         jmp(3) = 1
         IF (i .EQ. 1) THEN
            j = 2
            k = 3
         ELSEIF (i .EQ. 2) THEN
            j = 1
            k = 3
         ELSE
            j = 1
            k = 2
         END IF

         DO a2=1, lM%bs(j)%nNo
            DO a3=1, lM%bs(k)%nNo
               c = jmp(i)*(a1-1) + jmp(j)*(a2-1) + jmp(k)*(a3-1) + 1
               p = c - jmp(i)

               w1      = xn(4,c)*alpha
               w2      = xn(4,p)*(1D0-alpha)
               xn(4,c) = w1 + w2
               w1      = w1/xn(4,c)
               w2      = w2/xn(4,c)
               xn(1,c) = w1*xn(1,c) + w2*xn(1,p)
               xn(2,c) = w1*xn(2,c) + w2*xn(2,p)
               xn(3,c) = w1*xn(3,c) + w2*xn(3,p)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE INTROW
!--------------------------------------------------------------------      
!     This routine sets x and %nW
      SUBROUTINE SETROW(a1n,a1o)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: a1n, a1o

      INTEGER i, j, k, cn, co, jmpn(nsd), jmpo(nsd), a2, a3

      i = dir
      IF (nsd .EQ. 2) THEN
         a2      = lM%bs(2)%nNo
         jmpn(1) = a2
         jmpn(2) = 1
         IF (i .EQ. 1) THEN
            j = 2
         ELSE
            j  = 1
!     If this is direction 2, nNo in this direction has changed            
            a2 = tbs%nNo
         END IF
         jmpo(1) = a2
         jmpo(2) = 1

         DO a2=1, lM%bs(j)%nNo
            cn = jmpn(i)*(a1n-1) + jmpn(j)*(a2-1) + 1
            co = jmpo(i)*(a1o-1) + jmpo(j)*(a2-1) + 1
            xn(:,cn) = xo(:,co)
         END DO
      ELSE
         a2      = lM%bs(2)%nNo
         a3      = lM%bs(3)%nNo
         jmpn(1) = a3*a2
         jmpn(2) = a3
         jmpn(3) = 1
         IF (i .EQ. 1) THEN
            j = 2
            k = 3
         ELSEIF (i .EQ. 2) THEN
            j  = 1
            k  = 3
            a2 = tbs%nNo
         ELSE
            j  = 1
            k  = 2
            a3 = tbs%nNo
         END IF
         jmpo(1) = a3*a2
         jmpo(2) = a3
         jmpo(3) = 1

         DO a2=1, lM%bs(j)%nNo
            DO a3=1, lM%bs(k)%nNo
               cn = jmpn(i)*(a1n-1) + jmpn(j)*(a2-1) + jmpn(k)*(a3-1) +1
               co = jmpo(i)*(a1o-1) + jmpo(j)*(a2-1) + jmpo(k)*(a3-1) +1
               xn(:,cn) = xo(:,co)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE SETROW

      END SUBROUTINE KNOTINS

!####################################################################
!     This routine constructs a NURBS based on %bs sub-structures
      SUBROUTINE CONSTNRB(lM)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      
      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM

      TYPE dType  
         INTEGER, ALLOCATABLE :: INN(:)
         INTEGER, ALLOCATABLE :: IEN(:,:)
      END TYPE dType 
      
      INTEGER i, j, k, ex, ey, ez, e, Ec, ax, ay, az, a, Ac, iFa,
     2   shN, shE, jn(nsd), je(nsd), e1, e2, a1, a2
      REAL(KIND=8) tmp, tmpO
      TYPE(dType) d(nsd)

      IF (lM%eType .NE. eType_NRB) err = "Wrong call to CONSTNRB"

      lM%lShpF = .FALSE.
      lM%eNoN  = 1
      lM%gnNo  = 1
      lM%gnEl  = 1
      lM%nG    = 1
      lM%nSl   = 1
      DO i=1, nsd
         lM%bs(i)%nEl = 0
         DO j=1, lM%bs(i)%n
            tmp = lM%bs(i)%xi(j)
            IF (j .GT. 1) THEN 
               IF (.NOT.ISZERO(tmpO,tmp)) THEN
                  lM%bs(i)%nEl = lM%bs(i)%nEl + 1
               END IF
            END IF
            tmpO = tmp
         END DO
         lM%bs(i)%nNo = lM%bs(i)%n - lM%bs(i)%p - 1
         lM%eNoN = lM%eNoN*(lM%bs(i)%p + 1)
         lM%gnNo = lM%gnNo*lM%bs(i)%nNo
         lM%gnEl = lM%gnEl*lM%bs(i)%nEl
         lM%nG   = lM%nG  *lM%bs(i)%nG
         lM%nSl  = lM%nSl *lM%bs(i)%nSl
      END DO
      ALLOCATE(lM%gIEN(lM%eNoN,lM%gnEl), lM%INN(nsd,lM%gnEl))
 
!     Setting the faces
      CALL SELECTELE(lM)
      DO iFa=1, lM%nFa
         CALL SELECTELEB(lM, lM%fa(iFa))
      END DO

!     First constructing local IEN and INN arrays
      DO i=1, nsd
         ALLOCATE(d(i)%INN(lM%bs(i)%nEl),
     2            d(i)%IEN(lM%bs(i)%p+1,lM%bs(i)%nEl))
         e = 0
         DO j=lM%bs(i)%p+1, lM%bs(i)%n - lM%bs(i)%p - 1
            IF (.NOT.ISZERO(lM%bs(i)%xi(j),lM%bs(i)%xi(j+1))) THEN
               e = e + 1
               d(i)%INN(e) = j
               DO a=1, lM%bs(i)%p+1
                  d(i)%IEN(a,e) = j + a - lM%bs(i)%p - 1
               END DO
            END IF
         END DO
      END DO

!     Now based on the local IEN, I will construct global IEN         
      e = 0
      IF (nsd .EQ. 2) THEN
         DO ex=1, lM%bs(1)%nEl
            DO ey=1, lM%bs(2)%nEl
               e = e + 1
               a = 0
               DO ax=1, lM%bs(1)%p + 1
                  DO ay=1, lM%bs(2)%p + 1
                     Ac = lM%bs(2)%nNo*(d(1)%IEN(ax,ex)-1)
     2                  +  d(2)%IEN(ay,ey)
                     a = a + 1
                     lM%gIEN(a,e) = Ac
                  END DO
               END DO
               lM%INN(1,e) = d(1)%INN(ex)
               lM%INN(2,e) = d(2)%INN(ey)
            END DO
         END DO
      ELSE
         DO ex=1, lM%bs(1)%nEl
            DO ey=1, lM%bs(2)%nEl
               DO ez=1, lM%bs(3)%nEl
                  e = e + 1
                  a = 0
                  DO ax=1, lM%bs(1)%p + 1
                     DO ay=1, lM%bs(2)%p + 1
                        DO az=1, lM%bs(3)%p + 1
                           Ac = lM%bs(3)%nNo*(lM%bs(2)%nNo
     3                        *(d(1)%IEN(ax,ex) - 1) 
     4                        + d(2)%IEN(ay,ey) - 1) 
     5                        + d(3)%IEN(az,ez)
                           a = a + 1
                           lM%gIEN(a,e) = Ac
                        END DO
                     END DO
                  END DO
                  lM%INN(1,e) = d(1)%INN(ex)
                  lM%INN(2,e) = d(2)%INN(ey)
                  lM%INN(3,e) = d(3)%INN(ez)
               END DO
            END DO
         END DO
      END IF
!     Now working on faces
      DO iFa=1, lM%nFa
         i = lM%fa(iFa)%d
         lM%fa(iFa)%nEl = lM%gnEl/lM%bs(i)%nEl
         lM%fa(iFa)%nNo = lM%gnNo/lM%bs(i)%nNo
         ALLOCATE(lM%fa(iFa)%gE(lM%fa(iFa)%nEl), 
     4      lM%fa(iFa)%IEN(lM%fa(iFa)%eNoN,lM%fa(iFa)%nEl))

!     If control point is at the beginnig of the vector, starting ele%gl
!     element is always 1. Otherwise it is calculated based on
!     (ex-1)*(nEl_2*nEl_3) + (ey-1)*nEl_3 + ez. "pon" denotes moving in
!     positive or negative direction. This is to have normal allways
!     pointing outward. 
         shN = 0
         shE = 0
         IF (nsd .EQ. 2) THEN
            jn(1) = lM%bs(2)%nNo
            jn(2) = 1
            je(1) = lM%bs(2)%nEl
            je(2) = 1

            IF (i .EQ. 1) THEN
               j = 2
               IF (MOD(iFa,2) .EQ. 0) THEN
                  shN = lM%gnNo - jn(i)
                  shE = lM%gnEl - je(i)
               END IF
            ELSE
               j = 1
               IF (MOD(iFa,2) .EQ. 1) THEN
                  shN = jn(i) - jn(j)
                  shE = je(i) - je(j)
               END IF
            END IF

            DO e=1, lM%fa(iFa)%nEl
               Ec = shE + je(j)*e
               lM%fa(iFa)%gE(e) = Ec
               DO a=1, lM%bs(j)%p + 1
                  Ac = shN + jn(j)*d(j)%IEN(a,e)
                  lM%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         ELSE ! 3D NURBS
            jn(1) = lM%bs(3)%nNo*lM%bs(2)%nNo
            jn(2) = lM%bs(3)%nNo
            jn(3) = 1
            je(1) = lM%bs(3)%nEl*lM%bs(2)%nEl
            je(2) = lM%bs(3)%nEl
            je(3) = 1

            IF (i .EQ. 1) THEN
               j = 2
               k = 3
               IF (MOD(iFa,2) .EQ. 0) THEN
                  shN = lM%gnNo - jn(1)
                  shE = lM%gnEl - je(1)
               END IF
            ELSE IF (i .EQ. 2) THEN
               j = 1
               k = 3
               IF (MOD(iFa,2) .EQ. 0) THEN
                  shN = jn(1) - jn(2)
                  shE = je(1) - je(2)
               END IF
            ELSE
               j = 1
               k = 2
               IF (MOD(iFa,2) .EQ. 1) THEN
                  shN = jn(3) - jn(2)
                  shE = je(3) - je(2)
               END IF
            END IF

            e = 0
            DO e1=1, lM%bs(j)%nEl
               DO e2=1, lM%bs(k)%nEl
                  a  = 0
                  e  = e + 1
                  Ec = shE + je(j)*(e1-1) + je(k)*e2
                  lM%fa(iFa)%gE(e) = Ec
                  DO a1=1, lM%bs(j)%p + 1
                     DO a2=1, lM%bs(k)%p + 1
                        a  = a + 1
                        Ac = shN + jn(j)*(d(j)%IEN(a1,e1) - 1) 
     2                           + jn(k)* d(k)%IEN(a2,e2)
                        lM%fa(iFa)%IEN(a,e) = Ac
                     END DO
                  END DO
               END DO
            END DO
         END IF
         CALL CALCNBC(lM, lM%fa(iFa))
      END DO

      RETURN
      END SUBROUTINE CONSTNRB


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
!     All the routines that require an interface are included here.
!     This Mainly involves small routines with a very well-defined
!     functionality. To use these routines, just add "USE ALLFUN" to
!     your routine.
!
!--------------------------------------------------------------------

      MODULE ALLFUN
      
      USE MATFUN
      
      IMPLICIT NONE

      INTERFACE Integ
         MODULE PROCEDURE IntegS, IntegV, IntegG, vInteg
      END INTERFACE Integ

      INTERFACE COMMU
         MODULE PROCEDURE COMMUS, COMMUV
      END INTERFACE COMMU

      INTERFACE MKC
         MODULE PROCEDURE MKCS, MKCV
      END INTERFACE MKC
      INTERFACE MKCI
         MODULE PROCEDURE MKCIS, MKCIV
      END INTERFACE MKCI

      INTERFACE GLOBAL
         MODULE PROCEDURE GLOBALS, GLOBALV
      END INTERFACE GLOBAL

      INTERFACE LOCAL
         MODULE PROCEDURE LOCALIS, LOCALRV
      END INTERFACE LOCAL

      INTERFACE DESTROY
         MODULE PROCEDURE DESTROYFACE, DESTROYMSH, DESTROYBC, DESTROYEQ,
     2      DESTROYSTACK, DESTROYQUEUE
      END INTERFACE DESTROY

      CONTAINS

!####################################################################
!     This routine integrate s over the surface faId.
      FUNCTION IntegS(lFa, s)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: s(:)
      REAL(KIND=8) IntegS

      INTEGER a, e, g, Ac, nNo
      REAL(KIND=8) sHat, Jac, n(nsd)

      nNo  = SIZE(s)
      IF (nNo .NE. tnNo) err = "Incompatible vector size in Integ"

      IntegS = 0D0
      DO e=1, lFa%nEl
!     Updating the shape functions, if this is a NURB
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(lFa%iM), lFa, e)
         DO g=1, lFa%nG
            CALL GNNB(lFa, e, g, n)
            Jac = SQRT(NORM(n))
!     Calculating the function value
            sHat = 0D0
            DO a=1, lFa%eNoN
               Ac   = lFa%IEN(a,e)
               sHat = sHat + s(Ac)*lFa%N(a,g)
            END DO
!     Now integrating
            IntegS = IntegS + Jac*lFa%w(g)*sHat
         END DO
      END DO
      IntegS = cm%reduce(IntegS)

      RETURN
      END FUNCTION IntegS

!--------------------------------------------------------------------
!     This routine integrate s over the surface faId.
      FUNCTION IntegV(lFa, s)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: s(:,:)
      REAL(KIND=8) IntegV

      INTEGER a, i, e, Ac, g, ierr, nNo
      REAL(KIND=8) sHat, n(nsd)

      nNo = SIZE(s,2)
      IF (SIZE(s,1).NE.nsd .OR. nNo.NE.tnNo)
     2   err = "Incompatible vector size in IntegV"

      IntegV = 0D0
      DO e=1, lFa%nEl
!     Updating the shape functions, if this is a NURB
         IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(lFa%iM), lFa, e)
         DO g=1, lFa%nG
!     Since sHat must be multiplied by Jac, I don't normalize it here
            CALL GNNB(lFa, e, g, n)
!     Calculating the function value
            sHat = 0D0
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               DO i=1, nsd
                  sHat = sHat + lFa%N(a,g)*s(i,Ac)*n(i)
               END DO
            END DO
!     Now integrating
            IntegV = IntegV + lFa%w(g)*sHat
         END DO
      END DO

      IF (cm%seq()) RETURN
      CALL MPI_ALLREDUCE(IntegV, sHat, 1, mpreal, MPI_SUM, cm%com(),
     2   ierr)
      IntegV = sHat

      RETURN
      END FUNCTION IntegV

!--------------------------------------------------------------------
!     This routine integrate s(l,u,:) over the surface faId.
      FUNCTION IntegG(lFa, s, l, uo)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=8), INTENT(IN) :: s(:,:)
      INTEGER, INTENT(IN) :: l
      INTEGER, INTENT(IN), OPTIONAL :: uo

      INTEGER a, u, nNo
      REAL(KIND=8) IntegG
      REAL(KIND=8), ALLOCATABLE :: sclr(:), vec(:,:)

      u = l
      IF (PRESENT(uo)) u = uo

      nNo = SIZE(s,2)
      IF (nNo .NE. tnNo)
     2   err = "Incompatible vector size in Integ"

      IntegG = 0D0
      IF (u-l+1 .EQ. nsd) THEN
         ALLOCATE (vec(nsd,tnNo))
         DO a=1, tnNo
            vec(:,a) = s(l:u,a)
         END DO
         IntegG = IntegV(lFa,vec)
      ELSE IF (l .EQ. u) THEN
         ALLOCATE (sclr(tnNo))
         DO a=1, tnNo
            sclr(a) = s(l,a)
         END DO
         IntegG = IntegS(lFa,sclr)
      ELSE
         err = "Unexpected dof in IntegG"
      END IF

      RETURN
      END FUNCTION IntegG

!--------------------------------------------------------------------
!     This routine integrate an equation over a particular domain
      FUNCTION vInteg(dId, s, l, u)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: s(:,:)
      INTEGER, INTENT(IN) :: dId
      INTEGER, INTENT(IN) :: l
      INTEGER, INTENT(IN) :: u
      REAL(KIND=8) vInteg

      INTEGER a, e, g, Ac, i, iM, eNoN
      REAL(KIND=8) Jac, tmp(nsd,nsd), sHat, rtmp
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), Nx(:,:)

      IF (SIZE(s,2) .NE. tnNo)
     2   err = "Incompatible vector size in Integ"

      vInteg = 0D0
      DO iM=1, nMsh
         eNoN = msh(iM)%eNoN
         IF (iM .NE. 1) DEALLOCATE(xl, Nx)
         ALLOCATE(xl(nsd,eNoN), Nx(nsd,eNoN))
         DO e=1, msh(iM)%nEl
            IF (dId.GT.0 .AND. ALLOCATED(msh(iM)%eId)) THEN
               IF (.NOT.BTEST(msh(iM)%eId(e),dId)) CYCLE
            END IF
!     Updating the shape functions, if this is a NURB
            IF (msh(iM)%eType .EQ. eType_NRB) CALL NRBNNX(msh(iM), e)
            DO a=1, eNoN
               Ac      = msh(iM)%IEN(a,e)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)
            END DO

            DO g=1, msh(iM)%nG
               IF (g.EQ.1 .OR. .NOT.msh(iM)%lShpF)
     2            CALL GNN(eNoN, msh(iM)%Nx(:,:,g), xl, Nx, Jac, tmp)
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               sHat = 0D0
               DO a=1, eNoN
                  Ac = msh(iM)%IEN(a,e)
                  IF (l .EQ. u) THEN
                     rtmp = s(l,Ac)
                  ELSE
                     rtmp = SQRT(NORM(s(l:u,Ac)))
                  END IF
                  sHat = sHat + rtmp*msh(iM)%N(a,g)
               END DO
               vInteg = vInteg + msh(iM)%w(g)*Jac*sHat
            END DO
         END DO
      END DO

      IF (cm%seq()) RETURN
      CALL MPI_ALLREDUCE(vInteg, sHat, 1, mpreal, MPI_SUM, cm%com(), i)
      vInteg = sHat

      RETURN
      END FUNCTION vInteg

!####################################################################
!     This a both way communication with three main part:
      SUBROUTINE COMMUS(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(INOUT) :: U(:)

      IF (cm%seq()) RETURN

      IF (SIZE(U,1) .NE. lhs%nNo) THEN
         err = "COMMU is only specified for vector with size nNo"
      END IF
      U = MKC(U)
      CALL FSILS_COMMUS(lhs, U)
      CALL MKCI(U)

      RETURN
      END SUBROUTINE COMMUS
!--------------------------------------------------------------------
      SUBROUTINE COMMUV(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(INOUT) :: U(:,:)

      INTEGER m

      IF (cm%seq()) RETURN

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lhs%nNo) THEN
         err = "COMMU is only specified for vector with size nNo"
      END IF

      U = MKC(U)
      CALL FSILS_COMMUV(lhs, m, U)
      CALL MKCI(U)

      RETURN
      END SUBROUTINE COMMUV

!####################################################################

      FUNCTION MKCS(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: U(:)
      REAL(KIND=8), ALLOCATABLE :: MKCS(:)

      INTEGER a

      IF (SIZE(U,1) .NE. lhs%nNo) THEN
         err = "MKC is only specified for vector with size nNo"
      END IF
      ALLOCATE(MKCS(lhs%nNo))
      IF (cm%seq()) THEN
         MKCS = U
      ELSE
         DO a=1, lhs%nNo
            MKCS(lhs%map(a)) = U(a)
         END DO
      END IF

      RETURN
      END FUNCTION MKCS
!--------------------------------------------------------------------
      FUNCTION MKCV(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), ALLOCATABLE :: MKCV(:,:)

      INTEGER m, a

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lhs%nNo) THEN
         err = "MKC is only specified for vector with size nNo"
      END IF

      ALLOCATE(MKCV(m,lhs%nNo))
      IF (cm%seq()) THEN
         MKCV = U
      ELSE
         DO a=1, lhs%nNo
            MKCV(:,lhs%map(a)) = U(:,a)
         END DO
      END IF

      RETURN
      END FUNCTION MKCV
!--------------------------------------------------------------------
      SUBROUTINE MKCIS(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(INOUT) :: U(:)

      INTEGER a
      REAL(KIND=8), ALLOCATABLE :: tmp(:)

      IF (cm%seq()) RETURN

      IF (SIZE(U,1) .NE. lhs%nNo) THEN
         err = "MKC is only specified for vector with size nNo"
      END IF
      ALLOCATE(tmp(lhs%nNo))
      tmp = U
      DO a=1, lhs%nNo
         U(a) = tmp(lhs%map(a))
      END DO

      RETURN
      END SUBROUTINE MKCIS
!--------------------------------------------------------------------
      SUBROUTINE MKCIV(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(INOUT) :: U(:,:)

      INTEGER m, a
      REAL(KIND=8), ALLOCATABLE :: tmp(:,:)

      IF (cm%seq()) RETURN

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lhs%nNo) THEN
         err = "MKC is only specified for vector with size nNo"
      END IF

      ALLOCATE(tmp(m,lhs%nNo))
      tmp = U
      DO a=1, lhs%nNo
         U(:,a) = tmp(:,lhs%map(a))
      END DO

      RETURN
      END SUBROUTINE MKCIV

!####################################################################

      FUNCTION GLOBALS(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: U(:)
      REAL(KIND=8), ALLOCATABLE :: GLOBALS(:)

      INTEGER Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      IF (SIZE(U,1) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALS(lM%gnNo))
         GLOBALS = U
         RETURN
      END IF

      ALLOCATE(ienU(lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(lM%eNoN,lM%gnEl), GLOBALS(lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALS(0))
      END IF

      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            ienU(a,e) = U(Ac)
         END DO
      END DO

      DO i=1, cm%np()
         disp(i)   = lM%eDist(i-1)*lM%eNoN
         sCount(i) = lM%eDist(i)*lM%eNoN - disp(i)
      END DO

      CALL MPI_GATHERV(ienU, lM%nEl*lM%eNoN, mpreal, gienU, sCount,
     2   disp, mpreal, master, cm%com(), ierr)

      IF (cm%slv()) RETURN

      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            GLOBALS(Ac) = gienU(a,e)
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALS
!--------------------------------------------------------------------
      FUNCTION GLOBALV(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), ALLOCATABLE :: GLOBALV(:,:)

      INTEGER m, Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALV(m,lM%gnNo))
         GLOBALV = U
         RETURN
      END IF

      ALLOCATE(ienU(m*lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(m*lM%eNoN,lM%gnEl), GLOBALV(m,lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALV(0,0))
      END IF

      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            DO i=1, m
               ienU(m*(a-1)+i,e) = U(i,Ac)
            END DO
         END DO
      END DO

      a = lM%eNoN*m
      DO i=1, cm%np()
         disp(i)   = lM%eDist(i-1)*a
         sCount(i) = lM%eDist(i)*a - disp(i)
      END DO

      CALL MPI_GATHERV(ienU, lM%nEl*a, mpreal, gienU, sCount, disp,
     2   mpreal, master, cm%com(), ierr)

      IF (cm%slv()) RETURN

      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            DO i=1, m
               GLOBALV(i,Ac) = gienU(m*(a-1)+i,e)
            END DO
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALV

!####################################################################
!     These routine use ltg to go from gtnNo to tnNo
      FUNCTION LOCALIS(U)
      USE COMMOD
      USE UTILMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: U(:)
      INTEGER, ALLOCATABLE :: LOCALIS(:)

      INTEGER Ac, a
      INTEGER, ALLOCATABLE :: tmpU(:)

      IF (.NOT.ALLOCATED(ltg)) err = "ltg is not set yet"
      IF (cm%mas()) THEN
         IF (SIZE(U,1) .NE. gtnNo) err = "LOCAL is only"//
     2      " specified for vector with size gnNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(LOCALIS(gtnNo))
         LOCALIS = U
         RETURN
      END IF

      ALLOCATE(LOCALIS(tnNo), tmpU(gtnNo))
      IF (cm%mas()) tmpU = U
      CALL cm%bcast(tmpU)
      DO a=1, tnNo
         Ac = ltg(a)
         LOCALIS(a) = tmpU(Ac)
      END DO

      RETURN
      END FUNCTION LOCALIS
!--------------------------------------------------------------------
      FUNCTION LOCALRV(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), ALLOCATABLE :: LOCALRV(:,:)

      INTEGER m, s, e, a, Ac
      REAL(KIND=8), ALLOCATABLE :: tmpU(:)

      IF (.NOT.ALLOCATED(ltg)) err = "ltg is not set yet"
      IF (cm%mas()) THEN
         m = SIZE(U,1)
         IF (SIZE(U,2) .NE. gtnNo) err = "LOCAL is only"//
     2      " specified for vector with size gtnNo"
      END IF
      CALL cm%bcast(m)

      IF (cm%seq()) THEN
         ALLOCATE(LOCALRV(m,gtnNo))
         LOCALRV = U
         RETURN
      END IF

      ALLOCATE(LOCALRV(m,tnNo), tmpU(m*gtnNo))
      IF (cm%mas()) THEN
         DO a=1, gtnNo
            s = m*(a-1) + 1
            e = m*a
            tmpU(s:e) = U(:,a)
         END DO
      END IF
      CALL cm%bcast(tmpU)
      DO a=1, tnNo
         Ac = ltg(a)
         s  = m*(Ac-1) + 1
         e  = m*Ac
         LOCALRV(:,a) = tmpU(s:e)
      END DO

      RETURN
      END FUNCTION LOCALRV

!####################################################################
!     This function returns the domain that a set of nodes belongs to
      FUNCTION DOMAIN(lM, iEq, e)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, e
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER DOMAIN

      INTEGER iDmn

      DOMAIN = 0
!     Domain Id of -1 counts for the entire domain
      DO iDmn=1, eq(iEq)%nDmn
         DOMAIN = iDmn
         IF (eq(iEq)%dmn(iDmn)%Id .EQ. -1) RETURN
      END DO
      IF (.NOT. ALLOCATED(lM%eId)) err = "eId is not allocated,"//
     2   " while Id.NE.-1"
      DO iDmn=1, eq(iEq)%nDmn
         DOMAIN = iDmn
         IF (BTEST(lM%eId(e),eq(iEq)%dmn(iDmn)%Id)) RETURN
      END DO
      err = "Unable to find the domain ID of element "//e

      RETURN
      END FUNCTION DOMAIN

!--------------------------------------------------------------------
!     This function is true if "phys" is solved on "node"
      PURE FUNCTION ISDOMAIN(iEq, node, phys)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, node, phys
      LOGICAL ISDOMAIN

      INTEGER iDmn

      ISDOMAIN = .FALSE.
      IF (ALLOCATED(dmnId)) THEN
         DO iDmn=1, eq(iEq)%nDmn
            IF (eq(iEq)%dmn(iDmn)%phys .EQ. phys) THEN
               IF (BTEST(dmnId(node),eq(iEq)%dmn(iDmn)%Id)) THEN
                  ISDOMAIN = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      ELSE
!     No domain partitioning exists, so single domain is assumed and we
!     only need to check that
            IF (eq(iEq)%dmn(1)%phys .EQ. phys) ISDOMAIN = .TRUE.
      END IF

      RETURN
      END FUNCTION ISDOMAIN

!####################################################################
!     This routine (get to block) searches for "kwd" in file "fid"
!     and returns the position at the next line. If "m" and/or "n"
!     are present, the size of matrix is checked to be compatible
      SUBROUTINE GTBLK(fid, kwd, n)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: fid
      CHARACTER(LEN=*), INTENT(IN) :: kwd
      INTEGER, INTENT(OUT) :: n

      INTEGER l
      CHARACTER(LEN=stdL) rLine

      l = LEN(kwd)
      DO
         READ(fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(2:l+1) .EQ. kwd) EXIT
      END DO
      rLine = rLine(l+2:stdL)
      CALL GET(rLine,n)
      RETURN

 001  err = "Block with keyword <"//kwd//"> was not found"

      END SUBROUTINE GTBLK

!####################################################################
!     These set of routines destroy an object.
      PURE SUBROUTINE DESTROYFACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(OUT) :: lFa

      IF (ALLOCATED(lFa%IEN))  DEALLOCATE(lFa%IEN)
      IF (ALLOCATED(lFa%nV))   DEALLOCATE(lFa%nV)
      IF (ALLOCATED(lFa%gN))   DEALLOCATE(lFa%gN)
      IF (ALLOCATED(lFa%gE))   DEALLOCATE(lFa%gE)
      IF (ALLOCATED(lFa%Nx))   DEALLOCATE(lFa%Nx)
      IF (ALLOCATED(lFa%w))    DEALLOCATE(lFa%w)
      IF (ALLOCATED(lFa%N))    DEALLOCATE(lFa%N)
      IF (ALLOCATED(lFa%gebc)) DEALLOCATE(lFa%gebc)
      IF (ALLOCATED(lFa%x))    DEALLOCATE(lFa%x)

      lFa%eType = eType_NA
      lFa%nEl   = 0
      lFa%nNo   = 0
      lFa%gnEl  = 0

      RETURN
      END SUBROUTINE DESTROYFACE
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYMSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(OUT) :: lM

      INTEGER iFa

      IF (ALLOCATED(lM%eDist))  DEALLOCATE(lM%eDist)
      IF (ALLOCATED(lM%eId))    DEALLOCATE(lM%eId)
      IF (ALLOCATED(lM%gIEN))   DEALLOCATE(lM%gIEN)
      IF (ALLOCATED(lM%IEN))    DEALLOCATE(lM%IEN)
      IF (ALLOCATED(lM%INN))    DEALLOCATE(lM%INN)
      IF (ALLOCATED(lM%gN))     DEALLOCATE(lM%gN)
      IF (ALLOCATED(lM%lN))     DEALLOCATE(lM%lN)
      IF (ALLOCATED(lM%Nx))     DEALLOCATE(lM%Nx)
      IF (ALLOCATED(lM%nW))     DEALLOCATE(lM%nW)
      IF (ALLOCATED(lM%bs))     DEALLOCATE(lM%bs)
      IF (ALLOCATED(lM%w))      DEALLOCATE(lM%w)
      IF (ALLOCATED(lM%N))      DEALLOCATE(lM%N)
      IF (ALLOCATED(lM%x))      DEALLOCATE(lM%x)
      IF (ALLOCATED(lM%otnIEN)) DEALLOCATE(lM%otnIEN)
      IF (ALLOCATED(lM%gpN))    DEALLOCATE(lM%gpN)
      IF (ALLOCATED(lM%fa)) THEN
         DO iFa=1, lM%nFa
            CALL DESTROYFACE(lM%fa(iFa))
         END DO
         DEALLOCATE(lM%fa)
      END IF

      lM%eType = eType_NA
      lM%gnEl  = 0
      lM%gnNo  = 0
      lM%nEl   = 0
      lM%nFa   = 0
      lM%nNo   = 0

      RETURN
      END SUBROUTINE DESTROYMSH
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYBC(lBc)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(OUT) :: lBc

      IF (ALLOCATED(lBc%gt)) THEN
         IF (ALLOCATED(lBc%gt%r)) DEALLOCATE(lBc%gt%r)
         IF (ALLOCATED(lBc%gt%i)) DEALLOCATE(lBc%gt%i)
         DEALLOCATE(lBc%gt)
      END IF
      IF (ALLOCATED(lBc%gm)) THEN
         IF (ALLOCATED(lBc%gm%t))    DEALLOCATE(lBc%gm%t)
         IF (ALLOCATED(lBc%gm%d))    DEALLOCATE(lBc%gm%d)
         DEALLOCATE(lBc%gm)
      END IF
      IF (ALLOCATED(lBc%gx))    DEALLOCATE(lBc%gx)
      IF (ALLOCATED(lBc%eDrn)) DEALLOCATE(lBc%eDrn)
      
      lBc%bType    = 0
      lBc%cplBCptr = 0
      lBc%eDrn     = 0
      lBc%g        = 0D0
      lBc%r        = 0D0

      RETURN
      END SUBROUTINE DESTROYBC
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYEQ(lEq)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(OUT) :: lEq

      INTEGER iBc

      IF (ALLOCATED(lEq%bc)) THEN
         DO iBc=1, lEq%nBc
            CALL DESTROYBC(lEq%bc(iBc))
         END DO
         DEALLOCATE(lEq%bc)
      END IF
      IF (ALLOCATED(lEq%output)) DEALLOCATE(lEq%output)
      IF (ALLOCATED(lEq%dmn))    DEALLOCATE(lEq%dmn)

      lEq%coupled = .TRUE.
      lEq%dof     = 0
      lEq%maxItr  = 5
      lEq%minItr  = 1
      lEq%nBc     = 0
      lEq%nDmn    = 0
      lEq%nOutput = 0
      lEq%dBr     = -4D1
      lEq%tol     = 1D64

      RETURN
      END SUBROUTINE DESTROYEQ
!--------------------------------------------------------------------
      SUBROUTINE DESTROYSTACK(stk)
      USE COMMOD
      IMPLICIT NONE
      TYPE(stackType), INTENT(OUT) :: stk

      IF (ALLOCATED(stk%v)) DEALLOCATE(stk%v)

      stk%n    = 0
      stk%maxN = 0

      RETURN
      END SUBROUTINE DESTROYSTACK
!--------------------------------------------------------------------
      SUBROUTINE DESTROYQUEUE(que)
      USE COMMOD
      IMPLICIT NONE
      TYPE(queueType), INTENT(OUT) :: que

      IF (ALLOCATED(que%v)) DEALLOCATE(que%v)

      que%n    = 0
      que%maxN = 0

      RETURN
      END SUBROUTINE DESTROYQUEUE

!####################################################################
!     Spliting "m" jobs between "n" workers. "b" contains amount of jobs
!     and "A" will store the distribution of jobs
      RECURSIVE SUBROUTINE SPLITJOBS(m,n,A,b)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m, n
      REAL(KIND=8), INTENT(IN) :: b(m)
      REAL(KIND=8), INTENT(OUT) :: A(m,n)

      INTEGER i, j, k, ml, nl, mr, nr
      REAL(KIND=8) sb, sbl, sl, optsl
      REAL(KIND=8), ALLOCATABLE :: Al(:,:), Ar(:,:), bl(:), br(:)

      IF (m.LE.0 .OR. n.LE.0) RETURN

!     Single worker, but multiple jobs.
      IF (n .EQ. 1) THEN
         j = 1
         DO i=1, m
            A(i,j) = b(i)
         END DO
         RETURN
      END IF

!     Multiple workers, but a single job.
      IF (m .EQ. 1) THEN
         i = 1
         DO j=1, n
            A(i,j) = b(i)/REAL(n,8)
         END DO
         RETURN
      END IF

!     Multiple workers and multiple jobs
!     This is the initial guess for nl, nr
      nl  = n/2
      nr  = n - nl
!     This is the total amount of work
      sb  = SUM(b)
!     The work that suppose to be done by "l"
      sbl = sb*REAL(nl,8)/REAL(n,8)

      sl = 0D0
      DO i=1, m
         IF (sl+b(i) .GT. sbl) EXIT
         sl = sl + b(i)
      END DO

      optsl = ABS(sl - sbl)
      ml = i - 1
      j  = 0
      DO i=ml+1, m
         IF (ABS(sl + b(i) - sbl) .LT. optsl) THEN
            j = i
            optsl = ABS(sl + b(i) - sbl)
         END IF
      END DO
      IF (j .NE. 0) ml = ml + 1
      mr = m - ml
      ALLOCATE (bl(ml), br(mr))
      IF (j .NE. 0) THEN
         DO i=1, ml-1
            bl(i) = b(i)
         END DO
         bl(ml) = b(j)

         k = 0
         DO i=ml, m
            IF (i .EQ. j) CYCLE
            k = k + 1
            br(k) = b(i)
         END DO
      ELSE
         bl = b(1:ml)
         br = b(ml+1:m)
      END IF
      nl = NINT(REAL(n,8)*SUM(bl)/sb)
      IF (nl .EQ. 0) nl = 1
      IF (nl .EQ. n) nl = n - 1
      nr = n - nl
      ALLOCATE (Al(ml,nl), Ar(mr,nr))

      CALL SPLITJOBS(ml,nl,Al,bl)
      CALL SPLITJOBS(mr,nr,Ar,br)

      A = 0D0
      IF (j .NE. 0) THEN
         DO i=1, ml-1
            A(i,1:nl) = Al(i,:)
         END DO
         A(j,1:nl) = Al(ml,:)

         k = 0
         DO i=ml, m
            IF (i .EQ. j) CYCLE
            k = k + 1
            A(i,nl+1:n) = Ar(k,:)
         END DO
      ELSE
         A(1:ml,1:nl) = Al
         A(ml+1:m,nl+1:n) = Ar
      END IF

      RETURN
      END SUBROUTINE SPLITJOBS

!####################################################################
!     Set domain ID to a given number for the entire or a range of
!     elements in a mesh
      SUBROUTINE SETDMNID(lM, iDmn, ifirst, ilast)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iDmn
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(IN), OPTIONAL :: ifirst, ilast

      INTEGER e, first, last

      first = 1
      IF (PRESENT(ifirst)) first = ifirst
      IF (first.LE.0 .OR. first.GT.lM%gnEl) err = "Out of range"

      last = lM%gnEl
      IF (PRESENT(ilast)) last = ilast
      IF (last.LT.first .OR. last.GT.lM%gnEl) err = "Out of range"

      IF (.NOT.ALLOCATED(lM%eId)) THEN
         ALLOCATE(lM%eId(lM%gnEl))
         lM%eId = 0
      END IF
      DO e=first, last
         lM%eId(e) = IBSET(lM%eId(e),iDmn)
      END DO

      RETURN
      END SUBROUTINE SETDMNID

!####################################################################
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE FINDFACE(faceName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: faceName
      INTEGER, INTENT(OUT) :: iM, iFa

      iFa = 0
      MY_LOOP : DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            IF (msh(iM)%fa(iFa)%name .EQ. faceName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP

      IF (iM .GT. nMsh) err="Unable to find face <"//TRIM(faceName)//">"

      RETURN
      END SUBROUTINE FINDFACE

!####################################################################
!     Computes the JACOBIAN of an element
      FUNCTION JACOBIAN(nDim,eNoN,x,Nxi) result(Jac)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nDim,eNoN
      REAL(KIND=8), INTENT(IN) :: x(nDim,eNoN), Nxi(nDim,eNoN)

      INTEGER :: a
      REAL(KIND=8) :: Jac, xXi(nDim,nDim)

      xXi = 0D0
      DO a=1, eNoN
         xXi(:,1) = xXi(:,1) + (x(:,a) * Nxi(1,a))
         xXi(:,2) = xXi(:,2) + (x(:,a) * Nxi(2,a))
         IF (nDim .EQ. 3) THEN
            xXi(:,3) = xXi(:,3) + (x(:,a) * Nxi(3,a))
         END IF
      END DO
      Jac = MAT_DET(xXi,nDim)

      RETURN
      END FUNCTION JACOBIAN

!####################################################################
!     Computes the Skewness of an element
      FUNCTION SKEWNESS(nDim,eNoN,x) result(SkF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nDim,eNoN
      REAL(KIND=8), INTENT(IN) :: x(nDim,eNoN)

      INTEGER :: i,j,a,cnt
      REAL(KIND=8) :: circumRad, integ_eq, integ_el
      REAL(KIND=8) :: SkF, Dmat(eNoN,nDim+2), Dsub(eNoN,nDim+1)
      REAL(KIND=8), DIMENSION(nDim+2) :: coeff,detD

      IF (nDim .EQ. 2) THEN
         coeff = (/1D0, -1D0, 1D0, -1D0/)
      ELSE
         coeff = (/1D0, 1D0, -1D0, 1D0, 1D0/)
      END IF

      Dmat(:,:) = 1D0
      Dsub(:,:) = 0D0

      DO a=1, eNoN
         Dmat(a,1) = sum(x(:,a)**2)
         Dmat(a,2:nDim+1) = x(:,a)
      END DO

      DO j=1, nDim+2
         cnt = 0
         DO i=1, nDim+2
            IF (i .EQ. j) CYCLE
            cnt = cnt+1
            Dsub(:,cnt) = Dmat(:,i)
         END DO
         detD(j) = coeff(j) * MAT_DET(Dsub,nDim+1)
      END DO

      circumRad = SQRT( sum(detD(2:nDim+1)**2) -
     2                  4D0*detD(1)*detD(nDim+2) ) /
     3                ( 2D0*ABS(detD(1)) )

      IF (nDim .EQ. 2) THEN
         integ_eq = 2.5D-1 * SQRT(27D0) * circumRad**2
         integ_el =   5D-1 * ABS(detD(1))
      ELSE IF (nDim .EQ. 3) THEN
         integ_eq = ( 8D0*circumRad**3 ) / SQRT(243D0)
         integ_el = ABS(detD(1)) / 6D0
      END IF

      SkF = ABS(integ_eq - integ_el) / integ_eq

      RETURN
      END FUNCTION SKEWNESS

!####################################################################
!     Computes the Aspect Ratio of an element
      FUNCTION ASPECTRATIO(nDim,eNoN,x) result(AR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nDim,eNoN
      REAL(KIND=8), INTENT(IN) :: x(nDim,eNoN)

      INTEGER :: i,j,ap,a,b,icol,irow
      INTEGER :: rowM(eNoN,eNoN-1), colM(nDim,nDim-1)
      REAL(KIND=8) :: AR, s(eNoN)
      REAL(KIND=8) :: Dsub(nDim,nDim), detD(nDim)
      
      IF (nDim .EQ. 2) THEN
         s(:) = 0D0
         DO a=1, eNoN
            ap = a+1
            IF (a .EQ. eNoN) ap = 1
            s(a) = SQRT( sum( (x(:,a)-x(:,ap))**2 ) )
         END DO

      ELSE IF (nDim .EQ. 3) THEN
         rowM = reshape( (/1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4/),
     2                   shape(rowM) )
         colM = reshape( (/1, 2, 3, 2, 3, 1/), shape(colM) )

         DO a=1, eNoN
            DO b=1, nDim
               Dsub(:,:) = 1D0
               DO i=1, eNoN-1
                  irow = rowM(a,i)
                  DO j=1, nDim-1
                     icol = colM(b,j)
                     Dsub(i,j) = x(icol,irow)
                  END DO ! j
               END DO ! i

               detD(b) = MAT_DET(Dsub,nDim)
            END DO ! b
            s(a) = 5D-1 * SQRT( sum( detD(:)**2 ) )
         END DO ! a
      END IF
      
      AR = MAXVAL(s)/MINVAL(s)

      RETURN
      END FUNCTION ASPECTRATIO

!####################################################################
      
      END MODULE ALLFUN


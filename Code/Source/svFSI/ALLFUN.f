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
         MODULE PROCEDURE GLOBALIS, GLOBALRS, GLOBALIV, GLOBALRV
      END INTERFACE GLOBAL

      INTERFACE LOCAL
         MODULE PROCEDURE LOCALIS, LOCALRV
      END INTERFACE LOCAL

      INTERFACE DESTROY
         MODULE PROCEDURE DESTROYFACE, DESTROYMSH, DESTROYBC, DESTROYEQ,
     2      DESTROYBS, DESTROYMB, DESTROYDATA, DESTROYADJ, DESTROYSTACK,
     3      DESTROYQUEUE, DESTROYTRACE, DESTROYFCELL, DESTROYIBCM
      END INTERFACE DESTROY

      INTERFACE GETNADJCNCY
         MODULE PROCEDURE GETNADJ_MSH, GETNADJ_FACE
      END INTERFACE

      INTERFACE GETEADJCNCY
         MODULE PROCEDURE GETEADJ_MSH, GETEADJ_FACE
      END INTERFACE

      INTERFACE IB_SYNC
         MODULE PROCEDURE IB_SYNCS, IB_SYNCV
      END INTERFACE

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

      INTEGER a, e, g, Ac, iM, eNoN, insd, ibl
      REAL(KIND=8) Jac, nV(nsd), sHat, tmp(nsd,nsd)

      REAL(KIND=8), ALLOCATABLE :: xl(:,:), sl(:), Nxi(:,:), Nx(:,:),
     2   tmps(:,:)

      TYPE(fCellType) :: fCell

      IF (SIZE(s,2) .NE. tnNo)
     2   err = "Incompatible vector size in Integ"

      vInteg = 0D0
      DO iM=1, nMsh
         eNoN = msh(iM)%eNoN
         insd = nsd
         IF (msh(iM)%lShl) insd = nsd-1

         ALLOCATE(xl(nsd,eNoN), Nxi(insd,eNoN), Nx(insd,eNoN),
     2      sl(eNoN), tmps(nsd,insd))

         DO e=1, msh(iM)%nEl
            IF (dId.GT.0 .AND. ALLOCATED(msh(iM)%eId)) THEN
               IF (.NOT.BTEST(msh(iM)%eId(e),dId)) CYCLE
            END IF
!     Updating the shape functions, if this is a NURB
            IF (msh(iM)%eType .EQ. eType_NRB) CALL NRBNNX(msh(iM), e)
            ibl = 0
            DO a=1, eNoN
               Ac      = msh(iM)%IEN(a,e)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)
               IF (l .EQ. u) THEN
                  sl(a) = s(l,Ac)
               ELSE
                  sl(a) = SQRT(NORM(s(l:u,Ac)))
               END IF
               ibl = ibl + iblank(Ac)
            END DO
            IF (ibl .EQ. eNoN) CYCLE

            IF (msh(iM)%iGC(e) .EQ. 1) THEN
               IF (msh(iM)%lShl) CYCLE
               IF (ib%fcFlag) THEN
                  CALL FC_INIT(fCell, msh(iM), xl)
                  CALL FC_SET(fCell, msh(iM), xl, ib%Uo)
                  CALL FC_vInteg(fCell, xl, sl, vInteg)
                  CALL DESTROY(fCell)
                  CYCLE
               END IF
            END IF

            DO g=1, msh(iM)%nG
               Nxi(:,:) = msh(iM)%Nx(:,:,g)
               IF (g.EQ.1 .OR. .NOT.msh(iM)%lShpF) THEN
                  IF (msh(iM)%lShl) THEN
                     CALL GNNS(eNoN, Nxi, xl, nV, tmps, tmps)
                     Jac = SQRT(NORM(nV))
                  ELSE
                     CALL GNN(eNoN, Nxi, xl, Nx, Jac, tmp)
                  END IF
               END IF
               IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

               sHat = 0D0
               DO a=1, eNoN
                  Ac = msh(iM)%IEN(a,e)
                  sHat = sHat + sl(a)*msh(iM)%N(a,g)
               END DO
               vInteg = vInteg + msh(iM)%w(g)*Jac*sHat
            END DO
         END DO

         DEALLOCATE(xl, Nxi, Nx, sl, tmps)
      END DO

      IF (cm%seq()) RETURN
      CALL MPI_ALLREDUCE(vInteg, sHat, 1, mpreal, MPI_SUM, cm%com(), g)
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
      FUNCTION GLOBALIS(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: U(:)
      INTEGER, ALLOCATABLE :: GLOBALIS(:)

      INTEGER Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      INTEGER, ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      IF (SIZE(U,1) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALIS(lM%gnNo))
         GLOBALIS = U
         RETURN
      END IF

      ALLOCATE(ienU(lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(lM%eNoN,lM%gnEl), GLOBALIS(lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALIS(0))
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

      CALL MPI_GATHERV(ienU, lM%nEl*lM%eNoN, mpint, gienU, sCount,
     2   disp, mpint, master, cm%com(), ierr)

      IF (cm%slv()) RETURN

      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            GLOBALIS(Ac) = gienU(a,e)
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALIS
!--------------------------------------------------------------------
      FUNCTION GLOBALRS(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: U(:)
      REAL(KIND=8), ALLOCATABLE :: GLOBALRS(:)

      INTEGER Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      IF (SIZE(U,1) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALRS(lM%gnNo))
         GLOBALRS = U
         RETURN
      END IF

      ALLOCATE(ienU(lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(lM%eNoN,lM%gnEl), GLOBALRS(lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALRS(0))
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
            GLOBALRS(Ac) = gienU(a,e)
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALRS
!--------------------------------------------------------------------
      FUNCTION GLOBALIV(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(IN) :: U(:,:)
      INTEGER, ALLOCATABLE :: GLOBALIV(:,:)

      INTEGER m, Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      INTEGER, ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALIV(m,lM%gnNo))
         GLOBALIV = U
         RETURN
      END IF

      ALLOCATE(ienU(m*lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(m*lM%eNoN,lM%gnEl), GLOBALIV(m,lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALIV(0,0))
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

      CALL MPI_GATHERV(ienU, lM%nEl*a, mpint, gienU, sCount, disp,
     2   mpint, master, cm%com(), ierr)

      IF (cm%slv()) RETURN

      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            DO i=1, m
               GLOBALIV(i,Ac) = gienU(m*(a-1)+i,e)
            END DO
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALIV
!--------------------------------------------------------------------
      FUNCTION GLOBALRV(lM, U)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=8), INTENT(IN) :: U(:,:)
      REAL(KIND=8), ALLOCATABLE :: GLOBALRV(:,:)

      INTEGER m, Ac, a, e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: ienU(:,:), gienU(:,:)

      ALLOCATE(sCount(cm%np()), disp(cm%np()))
      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. lM%nNo) THEN
         err = "GLOBAL is only specified for vector with size nNo"
      END IF

      IF (cm%seq()) THEN
         ALLOCATE(GLOBALRV(m,lM%gnNo))
         GLOBALRV = U
         RETURN
      END IF

      ALLOCATE(ienU(m*lM%eNoN,lM%nEl))
      IF (cm%mas()) THEN
         ALLOCATE(gienU(m*lM%eNoN,lM%gnEl), GLOBALRV(m,lM%gnNo))
      ELSE
         ALLOCATE(gienU(0,0), GLOBALRV(0,0))
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
               GLOBALRV(i,Ac) = gienU(m*(a-1)+i,e)
            END DO
         END DO
      END DO

      RETURN
      END FUNCTION GLOBALRV
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
!     This function returns the domain that an element of a mesh
!     belongs to
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
      DOMAIN = 0
!      err = "Unable to find the domain ID of element "//e//", mesh <"//
!     2   TRIM(lM%name)//"> and equation "//iEq

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
!     This function returns the IB domain that an element of an IB mesh
!     belongs to
      FUNCTION IB_DOMAIN(lM, iEq, e)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, e
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER IB_DOMAIN

      INTEGER iDmn

      IB_DOMAIN = 0
!     Domain Id of -1 counts for the entire domain
      DO iDmn=1, eq(iEq)%nDmnIB
         IB_DOMAIN = iDmn
         IF (eq(iEq)%dmnIB(iDmn)%Id .EQ. -1) RETURN
      END DO
      IF (.NOT. ALLOCATED(lM%eId)) err = "eId is not allocated,"//
     2   " while Id.NE.-1"
      DO iDmn=1, eq(iEq)%nDmnIB
         IB_DOMAIN = iDmn
         IF (BTEST(lM%eId(e),eq(iEq)%dmnIB(iDmn)%Id)) RETURN
      END DO
      IB_DOMAIN = 0
!      err = "Unable to find the domain ID of element "//e//", mesh <"//
!     2   TRIM(lM%name)//"> and equation "//iEq

      RETURN
      END FUNCTION IB_DOMAIN
!--------------------------------------------------------------------
!     This function is true if "phys" is solved on "node"
      PURE FUNCTION IB_ISDOMAIN(iEq, node, phys)
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: iEq, node, phys
      LOGICAL IB_ISDOMAIN

      INTEGER iDmn

      IB_ISDOMAIN = .FALSE.
      IF (ALLOCATED(ib%dmnId)) THEN
         DO iDmn=1, eq(iEq)%nDmnIB
            IF (eq(iEq)%dmnIB(iDmn)%phys .EQ. phys) THEN
               IF (BTEST(ib%dmnId(node),eq(iEq)%dmnIB(iDmn)%Id)) THEN
                  IB_ISDOMAIN = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      ELSE
!     No domain partitioning exists, so single domain is assumed and we
!     only need to check that
            IF (eq(iEq)%dmnIB(1)%phys .EQ. phys) IB_ISDOMAIN = .TRUE.
      END IF

      RETURN
      END FUNCTION IB_ISDOMAIN
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

      IF (ALLOCATED(lFa%gE))   DEALLOCATE(lFa%gE)
      IF (ALLOCATED(lFa%gN))   DEALLOCATE(lFa%gN)
      IF (ALLOCATED(lFa%lN))   DEALLOCATE(lFa%lN)
      IF (ALLOCATED(lFa%IEN))  DEALLOCATE(lFa%IEN)
      IF (ALLOCATED(lFa%gebc)) DEALLOCATE(lFa%gebc)
      IF (ALLOCATED(lFa%w))    DEALLOCATE(lFa%w)
      IF (ALLOCATED(lFa%x))    DEALLOCATE(lFa%x)
      IF (ALLOCATED(lFa%xi))   DEALLOCATE(lFa%xi)
      IF (ALLOCATED(lFa%N))    DEALLOCATE(lFa%N)
      IF (ALLOCATED(lFa%nV))   DEALLOCATE(lFa%nV)
      IF (ALLOCATED(lFa%Nx))   DEALLOCATE(lFa%Nx)
      IF (ALLOCATED(lFa%Nxx))  DEALLOCATE(lFa%Nxx)

      CALL DESTROYTRACE(lFa%trc)

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

      INTEGER iFa, i, insd

      insd = nsd
      IF (lM%lShl) insd = nsd - 1

      IF (ALLOCATED(lM%eDist))  DEALLOCATE(lM%eDist)
      IF (ALLOCATED(lM%eId))    DEALLOCATE(lM%eId)
      IF (ALLOCATED(lM%gN))     DEALLOCATE(lM%gN)
      IF (ALLOCATED(lM%gpN))    DEALLOCATE(lM%gpN)
      IF (ALLOCATED(lM%gIEN))   DEALLOCATE(lM%gIEN)
      IF (ALLOCATED(lM%IEN))    DEALLOCATE(lM%IEN)
      IF (ALLOCATED(lM%otnIEN)) DEALLOCATE(lM%otnIEN)
      IF (ALLOCATED(lM%INN))    DEALLOCATE(lM%INN)
      IF (ALLOCATED(lM%lN))     DEALLOCATE(lM%lN)
      IF (ALLOCATED(lM%eIEN))   DEALLOCATE(lM%eIEN)
      IF (ALLOCATED(lM%sbc))    DEALLOCATE(lM%sbc)
      IF (ALLOCATED(lM%iGC))    DEALLOCATE(lM%iGC)
      IF (ALLOCATED(lM%nW))     DEALLOCATE(lM%nW)
      IF (ALLOCATED(lM%w))      DEALLOCATE(lM%w)
      IF (ALLOCATED(lM%xi))     DEALLOCATE(lM%xi)
      IF (ALLOCATED(lM%x))      DEALLOCATE(lM%x)
      IF (ALLOCATED(lM%N))      DEALLOCATE(lM%N)
      IF (ALLOCATED(lM%nV))     DEALLOCATE(lM%nV)
      IF (ALLOCATED(lM%Nx))     DEALLOCATE(lM%Nx)
      IF (ALLOCATED(lM%Nxx))    DEALLOCATE(lM%Nxx)

      IF (ALLOCATED(lM%bs)) THEN
         DO i=1, insd
            CALL DESTROY(lM%bs(i))
         END DO
         DEALLOCATE(lM%bs)
      END IF

      IF (ALLOCATED(lM%fa)) THEN
         DO iFa=1, lM%nFa
            CALL DESTROYFACE(lM%fa(iFa))
         END DO
         DEALLOCATE(lM%fa)
      END IF

      CALL DESTROYADJ(lM%nAdj)
      CALL DESTROYADJ(lM%eAdj)
      CALL DESTROYTRACE(lM%trc)
      IF (ALLOCATED(lM%bf)) THEN
         CALL DESTROYMB(lM%bf)
         DEALLOCATE(lM%bf)
      END IF

      lM%lShl  = .FALSE.
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
      IF (ALLOCATED(lBc%h))    DEALLOCATE(lBc%h)
      IF (ALLOCATED(lBc%gx))   DEALLOCATE(lBc%gx)
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
      IF (ALLOCATED(lEq%bcIB)) THEN
         DO iBc=1, lEq%nBcIB
            CALL DESTROYBC(lEq%bcIB(iBc))
         END DO
         DEALLOCATE(lEq%bcIB)
      END IF
      IF (ALLOCATED(lEq%output)) DEALLOCATE(lEq%output)
      IF (ALLOCATED(lEq%dmn))    DEALLOCATE(lEq%dmn)
      IF (ALLOCATED(lEq%dmnIB))  DEALLOCATE(lEq%dmnIB)

      lEq%coupled = .TRUE.
      lEq%dof     = 0
      lEq%maxItr  = 5
      lEq%minItr  = 1
      lEq%nBc     = 0
      lEq%nDmn    = 0
      lEq%nOutput = 0
      lEq%dBr     = -4D1
      lEq%tol     = 1D64
      lEq%nBcIB   = 0
      lEq%nDmnIB  = 0

      RETURN
      END SUBROUTINE DESTROYEQ
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYMB(lMB)
      USE COMMOD
      IMPLICIT NONE
      TYPE(MBType), INTENT(OUT) :: lMB

      IF (ALLOCATED(lMB%t)) DEALLOCATE(lMB%t)
      IF (ALLOCATED(lMB%d)) DEALLOCATE(lMB%d)

      lMB%nTP = 0

      RETURN
      END SUBROUTINE DESTROYMB
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYBS(bs)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bsType), INTENT(OUT) :: bs

      IF (ALLOCATED(bs%xi)) DEALLOCATE(bs%xi)

      RETURN
      END SUBROUTINE DESTROYBS
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYDATA(d)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dataType), INTENT(OUT) :: d

      IF (ALLOCATED(d%IEN))   DEALLOCATE(d%IEN)
      IF (ALLOCATED(d%xe))    DEALLOCATE(d%xe)
      IF (ALLOCATED(d%gx))    DEALLOCATE(d%gx)
      IF (ALLOCATED(d%x))     DEALLOCATE(d%x)

      RETURN
      END SUBROUTINE DESTROYDATA
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYADJ(adj)
      USE COMMOD
      IMPLICIT NONE
      TYPE(adjType), INTENT(OUT) :: adj

      IF (ALLOCATED(adj%prow)) DEALLOCATE(adj%prow)
      IF (ALLOCATED(adj%pcol)) DEALLOCATE(adj%pcol)
      adj%nnz = 0

      RETURN
      END SUBROUTINE DESTROYADJ
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYSTACK(stk)
      USE COMMOD
      IMPLICIT NONE
      TYPE(stackType), INTENT(OUT) :: stk

      IF (ALLOCATED(stk%v)) DEALLOCATE(stk%v)

      stk%n    = 0
      stk%maxN = 0

      RETURN
      END SUBROUTINE DESTROYSTACK
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYQUEUE(que)
      USE COMMOD
      IMPLICIT NONE
      TYPE(queueType), INTENT(OUT) :: que

      IF (ALLOCATED(que%v)) DEALLOCATE(que%v)

      que%n    = 0
      que%maxN = 0

      RETURN
      END SUBROUTINE DESTROYQUEUE
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYTRACE(trc)
      USE COMMOD
      IMPLICIT NONE
      TYPE(traceType), INTENT(OUT) :: trc

      IF (ALLOCATED(trc%gE))  DEALLOCATE(trc%gE)
      IF (ALLOCATED(trc%ptr)) DEALLOCATE(trc%ptr)
      trc%n = 0

      RETURN
      END SUBROUTINE DESTROYTRACE
!--------------------------------------------------------------------
      RECURSIVE SUBROUTINE DESTROYFCELL(lFc)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fCellType), INTENT(OUT) :: lFc
      INTEGER :: i

      IF (lFc%nsub .GT. 0) THEN
         DO i=1, lFc%nSub
            CALL DESTROY(lFc%sub(i))
         END DO
         NULLIFY(lFc%sub)
      END IF

      lFc%eType = eType_NA
      lFc%ilev  = 0
      lFc%nsub  = 0
      IF (ALLOCATED(lFc%w))    DEALLOCATE(lFc%w)
      IF (ALLOCATED(lFc%x))    DEALLOCATE(lFc%x)
      IF (ALLOCATED(lFc%xi))   DEALLOCATE(lFc%xi)
      IF (ALLOCATED(lFc%xiGP)) DEALLOCATE(lFc%xiGP)
      IF (ALLOCATED(lFc%incG)) DEALLOCATE(lFc%incG)

      RETURN
      END SUBROUTINE DESTROYFCELL
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYIBCM(ibCm)
      USE COMMOD
      IMPLICIT NONE
      TYPE(ibCommType), INTENT(OUT) :: ibCm

      IF (ALLOCATED(ibCm%n))  DEALLOCATE(ibCm%n)
      IF (ALLOCATED(ibCm%gE)) DEALLOCATE(ibCm%gE)

      RETURN
      END SUBROUTINE DESTROYIBCM
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
         lM%eId(e) = IBSET(lM%eId(e), iDmn)
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
      FUNCTION JACOBIAN(nDim, eNoN, x, Nxi) result(Jac)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nDim, eNoN
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
!--------------------------------------------------------------------
!     Computes the Skewness of an element
      FUNCTION SKEWNESS(nDim, eNoN, x) result(SkF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nDim, eNoN
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
!--------------------------------------------------------------------
!     Computes the Aspect Ratio of an element
      FUNCTION ASPECTRATIO(nDim, eNoN, x) result(AR)
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
!     Find nodal adjacency of a given mesh. Computes list of all nodes
!     around a given node of a mesh.
      SUBROUTINE GETNADJ_MSH(lM)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER :: a, b, e, Ac, Bc, i, j, maxAdj
      LOGICAL :: flag

      INTEGER, ALLOCATABLE :: incNd(:), adjL(:,:)

      ALLOCATE(incNd(lM%nNo))
      incNd = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            DO b=1, lM%eNoN
               IF (b .EQ. a) CYCLE
               incNd(Ac) = incNd(Ac) + 1
            END DO
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(adjL(maxAdj, lM%nNo))
      incNd = 0
      adjL  = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            DO b=1, lM%eNoN
               IF (b .EQ. a) CYCLE
               Bc = lM%IEN(b,e)
               Bc = lM%lN(Bc)
               flag = .TRUE.
               DO i=1, incNd(Ac)
                  IF (adjL(i,Ac) .EQ. Bc) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(Ac) = incNd(Ac) + 1
                  adjL(incNd(Ac),Ac) = Bc
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE(incNd)

      lM%nAdj%nnz = 0
      DO a=1, lM%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            lM%nAdj%nnz = lM%nAdj%nnz + 1
         END DO
      END DO

      ALLOCATE(lM%nAdj%prow(lM%nNo+1), lM%nAdj%pcol(lM%nAdj%nnz))
      j = 0
      lM%nAdj%prow(1) = j + 1
      DO a=1, lM%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            j = j + 1
            lM%nAdj%pcol(j) = adjL(i,a)
         END DO
         lM%nAdj%prow(a+1) = j + 1
      END DO
      DEALLOCATE(adjL)

      RETURN
      END SUBROUTINE GETNADJ_MSH
!--------------------------------------------------------------------
!     Find nodal adjacency of a given face. Computes list of all nodes
!     around a given node of a face.
      SUBROUTINE GETNADJ_FACE(lFa)
      USE COMMOD
      IMPLICIT NONE

      TYPE(faceType),  INTENT(INOUT) :: lFa

      INTEGER :: a, b, e, Ac, Bc, i, j, maxAdj
      LOGICAL :: flag

      INTEGER, ALLOCATABLE :: incNd(:), adjL(:,:)

      ALLOCATE(incNd(lFa%nNo))
      incNd = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO b=1, lFa%eNoN
               IF (b .EQ. a) CYCLE
               incNd(Ac) = incNd(Ac) + 1
            END DO
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(adjL(maxAdj, lFa%nNo))
      incNd = 0
      adjL  = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO b=1, lFa%eNoN
               IF (b .EQ. a) CYCLE
               Bc = lFa%IEN(b,e)
               Bc = lfa%lN(Bc)
               flag = .TRUE.
               DO i=1, incNd(Ac)
                  IF (adjL(i,Ac) .EQ. Bc) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(Ac) = incNd(Ac) + 1
                  adjL(incNd(Ac),Ac) = Bc
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE(incNd)

      lFa%nAdj%nnz = 0
      DO a=1, lFa%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            lFa%nAdj%nnz = lFa%nAdj%nnz + 1
         END DO
      END DO

      ALLOCATE(lFa%nAdj%prow(lFa%nNo+1), lFa%nAdj%pcol(lFa%nAdj%nnz))
      j = 0
      lFa%nAdj%prow(1) = j + 1
      DO a=1, lFa%nNo
         DO i=1, maxAdj
            b = adjL(i,a)
            IF (b .EQ. 0) EXIT
            j = j + 1
            lFa%nAdj%pcol(j) = adjL(i,a)
         END DO
         lFa%nAdj%prow(a+1) = j + 1
      END DO
      DEALLOCATE(adjL)

      RETURN
      END SUBROUTINE GETNADJ_FACE
!####################################################################
!     Find element adjacency of a given mesh. Computes list of all
!     elements around a given element of a mesh.
      SUBROUTINE GETEADJ_MSH(lM)
      USE COMMOD
      IMPLICIT NONE

      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER, ALLOCATABLE :: incNd(:), adjList(:,:), tmpI(:,:)

      ALLOCATE(incNd(lM%nNo))
      incNd = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(tmpI(maxAdj, lM%nNo))
      incNd = 0
      tmpI  = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
            tmpI(incNd(Ac), Ac) = e
         END DO
      END DO
      b = 2*maxAdj

 001  b = b + maxAdj
      DEALLOCATE(incNd)
      ALLOCATE(incNd(lM%nEl))
      IF (ALLOCATED(adjList)) DEALLOCATE(adjList)
      ALLOCATE(adjList(b, lM%nEl))
      adjList = 0
      incNd   = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            DO i=1, maxAdj
               IF (tmpI(i,Ac) .EQ. 0) EXIT
               flag = .TRUE.
               DO j=1, incNd(e)
                  IF (adjList(j,e) .EQ. tmpI(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(e) = incNd(e) + 1
                  IF (incNd(e) .GE. b) GOTO 001
                  adjList(incNd(e),e) = tmpI(i,Ac)
               END IF
            END DO
         END DO
      END DO
      maxAdj = MAXVAL(incNd)
      DEALLOCATE(tmpI, incNd)

      lM%eAdj%nnz = 0
      DO e=1, lM%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               lM%eAdj%nnz = lM%eAdj%nnz + 1
            ELSE
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE(lM%eAdj%prow(lM%nEl+1), lM%eAdj%pcol(lM%eAdj%nnz))
      j = 0
      lM%eAdj%prow(1) = j + 1
      DO e=1, lM%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               j = j + 1
               lM%eAdj%pcol(j) = adjList(i,e)
            ELSE
               EXIT
            END IF
         END DO
         lM%eAdj%prow(e+1) = j + 1
      END DO
      DEALLOCATE(adjList)

      RETURN
      END SUBROUTINE GETEADJ_MSH
!--------------------------------------------------------------------
!     Find element adjacency of a given face. Computes list of all
!     elements around a given element of a face.
      SUBROUTINE GETEADJ_FACE(lFa)
      USE COMMOD
      IMPLICIT NONE

      TYPE(faceType),  INTENT(INOUT) :: lFa

      INTEGER :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER, ALLOCATABLE :: incNd(:), adjList(:,:), tmpI(:,:)

      ALLOCATE(incNd(lFa%nNo))
      incNd = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
         END DO
      END DO

      maxAdj = MAXVAL(incNd)
      ALLOCATE(tmpI(maxAdj, lFa%nNo))
      incNd = 0
      tmpI  = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            incNd(Ac) = incNd(Ac) + 1
            tmpI(incNd(Ac), Ac) = e
         END DO
      END DO
      b = 2*maxAdj

 001  b = b + maxAdj
      DEALLOCATE(incNd)
      ALLOCATE(incNd(lFa%nEl))
      IF (ALLOCATED(adjList)) DEALLOCATE(adjList)
      ALLOCATE(adjList(b, lFa%nEl))
      adjList = 0
      incNd   = 0
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            DO i=1, maxAdj
               IF (tmpI(i,Ac) .EQ. 0) EXIT
               flag = .TRUE.
               DO j=1, incNd(e)
                  IF (adjList(j,e) .EQ. tmpI(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incNd(e) = incNd(e) + 1
                  IF (incNd(e) .GE. b) GOTO 001
                  adjList(incNd(e),e) = tmpI(i,Ac)
               END IF
            END DO
         END DO
      END DO
      maxAdj = MAXVAL(incNd)
      DEALLOCATE(tmpI, incNd)

      lFa%eAdj%nnz = 0
      DO e=1, lFa%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               lFa%eAdj%nnz = lFa%eAdj%nnz + 1
            ELSE
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE(lFa%eAdj%prow(lFa%nEl+1), lFa%eAdj%pcol(lFa%eAdj%nnz))
      j = 0
      lFa%eAdj%prow(1) = j + 1
      DO e=1, lFa%nEl
         DO i=1, maxAdj
            IF (adjList(i,e) .NE. 0) THEN
               j = j + 1
               lFa%eAdj%pcol(j) = adjList(i,e)
            ELSE
               EXIT
            END IF
         END DO
         lFa%eAdj%prow(e+1) = j + 1
      END DO
      DEALLOCATE(adjList)

      RETURN
      END SUBROUTINE GETEADJ_FACE
!####################################################################
!     Finds an element of a mesh for any probe. Returns 0 if not found.
!     Uses Newton method to compute the parametric coordinate (xi) of
!     the probe with respect to an element for the given physical
!     coordinate (xp). Elements are searched prescribed by eList.
      SUBROUTINE FINDE(xp, lM, xg, Dg, nNo, ne, eList, Ec, xi)
      USE COMMOD
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nNo, ne, eList(ne)
      REAL(KIND=8), INTENT(IN) :: xp(nsd), xg(nsd,nNo), Dg(nsd,nNo)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER, INTENT(OUT) :: Ec
      REAL(KIND=8), INTENT(OUT) :: xi(nsd)

      LOGICAL :: flag
      INTEGER :: a, e, i, Ac, eNoN
      REAL(KIND=8) :: rt

      LOGICAL, ALLOCATABLE :: eChck(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:)

      IF (lM%lShl) err = " Finding traces is not applicable for shells"

      Ec   = 0
      eNoN = lM%eNoN
      ALLOCATE(eChck(lM%nEl), xl(nsd,eNoN), N(eNoN), Nxi(nsd,eNoN))

c      WRITE(1000+cm%tF(),'(A)') "=================================="
c      WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') "Finding trace for xp: "
c      DO i=1,  nsd
c         WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xp(i))
c      END DO
c      WRITE(1000+cm%tF(),'(A)')
c      WRITE(1000+cm%tF(),'(A)')
c      WRITE(1000+cm%tF(),'(A)') "List of elems: <"//STR(ne)//">"
c      DO e=1, ne
c         WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(eList(e))
c      END DO
c      WRITE(1000+cm%tF(),'(A)')
c      WRITE(1000+cm%tF(),'(A)') "=================================="

      eChck = .FALSE.
      DO e=1, ne
         Ec = eList(e)
         IF (Ec .EQ. 0) CYCLE
         IF (eChck(Ec)) CYCLE
         eChck(Ec) = .TRUE.

         DO a=1, eNoN
            Ac = lM%IEN(a,Ec)
            xi(:) = xi(:) + lM%xi(:,a)
            xl(:,a) = xg(:,Ac) + Dg(:,Ac)
         END DO

         xi = 0D0
         DO i=1, lM%nG
            xi = xi + lM%xi(:,i)
         END DO
         xi(:) = xi(:) / REAL(lM%nG,KIND=8)

c         WRITE(1000+cm%tF(),'(A)') "----------------------------"
c         WRITE(1000+cm%tF(),'(A)') "Probe el: "//STR(Ec)
c         DO a=1, eNoN
c            Ac = lM%IEN(a,Ec)
c            WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') STR(Ac)
c            DO i=1, nsd
c               WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xl(i,a))
c            END DO
c            WRITE(1000+cm%tF(),'(A)')
c         END DO

         CALL GETXI(lM%eType, eNoN, xl, xp, xi, flag)
!        If Newton's method fails, continue
         IF (.NOT.flag) CYCLE

c         WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') "xi: "
c         DO i=1, nsd
c            WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xi(i))
c         END DO
c         WRITE(1000+cm%tF(),'(A)')

!        Check for shape function even if the Newton's method returns OK
         CALL GETGNN(nsd, lM%eType, eNoN, xi, N, Nxi)

c         WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') "N: "
c         DO a=1, eNoN
c            WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(N(a))
c         END DO
c         WRITE(1000+cm%tF(),'(A)')
c         CALL FLUSH(1000+cm%tF())

         i  = 0
         rt = 0.0D0
         DO a=1, eNoN
            rt = rt + N(a)
            IF (N(a).GT.-1E-4 .AND. N(a).LT.1.0001D0) i = i + 1
         END DO
         flag = rt.GE.0.9999D0 .AND. rt.LE.1.0001D0
         IF (i.EQ.eNoN .AND. flag) RETURN
      END DO

      Ec = 0
      xi = 0D0

      RETURN
      END SUBROUTINE FINDE
!####################################################################
      SUBROUTINE IB_SYNCS(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: U(:)

      INTEGER i, a, e, Ac, iM, iFa, nl, ng, ierr, tag, sReq

      INTEGER, ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=8), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      IF (SIZE(U,1) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%n)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            DO i=1, ib%msh(iM)%trc%n
               e = ib%msh(iM)%trc%gE(1,i)
               DO a=1, ib%msh(iM)%eNoN
                  Ac = ib%msh(iM)%IEN(a,e)
                  incNd(Ac) = 1
               END DO
            END DO
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               DO i=1, ib%msh(iM)%fa(iFa)%trc%n
                  e = ib%msh(iM)%fa(iFa)%trc%gE(1,i)
                  DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                     Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                     incNd(Ac) = 1
                  END DO
               END DO
            END DO
         END IF
      END DO

      nl = SUM(incNd)
      ALLOCATE(lU(nl))
      nl = 0
      DO a=1, ib%tnNo
         IF (incNd(a) .EQ. 1) THEN
            nl = nl + 1
            lU(nl) = U(a)
         END IF
      END DO
      DEALLOCATE(incNd)

      ng = SUM(ib%cm%n)
      ALLOCATE(gU(ng))
      IF (cm%mas()) THEN
         ALLOCATE(rReq(cm%np()))
         rReq = 0
         gU   = 0D0
         ng   = 0
         DO i=1, cm%np()
            nl = ib%cm%n(i)
            IF (nl .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               gU(ng+1:ng+nl) = lU(1:nl)
            ELSE
               tag = i * 100
               CALL MPI_IRECV(gU(ng+1:ng+nl), nl, mpreal, i-1, tag,
     2            cm%com(), rReq(i), ierr)
            END IF
            ng = ng + nl
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ib%cm%n(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO

         U = 0D0
         DO a=1, ng
            Ac    = ib%cm%gE(a)
            U(Ac) = U(Ac) + gU(a)
         END DO
      ELSE IF (nl .NE. 0) THEN
         tag = cm%tF() * 100
         CALL MPI_SEND(lU, nl, mpreal, master, tag, cm%com(), ierr)
      END IF

      DEALLOCATE(lU, gU)

      CALL cm%bcast(U)

      RETURN
      END SUBROUTINE IB_SYNCS
!--------------------------------------------------------------------
      SUBROUTINE IB_SYNCV(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=8), INTENT(INOUT) :: U(:,:)

      INTEGER m, i, a, e, s, Ac, iM, iFa, nl, ng, ierr, tag, sReq

      INTEGER, ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=8), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%n)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            DO i=1, ib%msh(iM)%trc%n
               e = ib%msh(iM)%trc%gE(1,i)
               DO a=1, ib%msh(iM)%eNoN
                  Ac = ib%msh(iM)%IEN(a,e)
                  incNd(Ac) = 1
               END DO
            END DO
         ELSE
            DO iFa=1, ib%msh(iM)%nFa
               DO i=1, ib%msh(iM)%fa(iFa)%trc%n
                  e = ib%msh(iM)%fa(iFa)%trc%gE(1,i)
                  DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                     Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                     incNd(Ac) = 1
                  END DO
               END DO
            END DO
         END IF
      END DO

      nl = SUM(incNd)
      ALLOCATE(lU(m*nl))
      nl = 0
      DO a=1, ib%tnNo
         IF (incNd(a) .EQ. 1) THEN
            nl = nl + 1
            s  = (nl-1)*m + 1
            e  = (nl-1)*m + m
            lU(s:e) = U(1:m,a)
         END IF
      END DO
      DEALLOCATE(incNd)

      ng = SUM(ib%cm%n)
      ALLOCATE(gU(m*ng))
      IF (cm%mas()) THEN
         ALLOCATE(rReq(cm%np()))
         rReq = 0
         gU   = 0D0
         ng   = 0
         DO i=1, cm%np()
            nl = ib%cm%n(i)
            IF (nl .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               s = 1
               e = nl*m
               gU(s:e) = lU(:)
            ELSE
               tag = i * 100
               s = ng*m + 1
               e = ng*m + nl*m
               CALL MPI_IRECV(gU(s:e), nl*m, mpreal, i-1, tag,
     2            cm%com(), rReq(i), ierr)
            END IF
            ng = ng + nl
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ib%cm%n(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO

         U = 0D0
         DO a=1, ng
            Ac = ib%cm%gE(a)
            s  = (a-1)*m + 1
            e  = (a-1)*m + m
            U(1:m,Ac) = U(1:m,Ac) + gU(s:e)
         END DO
      ELSE IF (nl .NE. 0) THEN
         tag = cm%tF() * 100
         CALL MPI_SEND(lU, nl*m, mpreal, master, tag, cm%com(), ierr)
      END IF

      DEALLOCATE(lU, gU)

      CALL cm%bcast(U)

      RETURN
      END SUBROUTINE IB_SYNCV
!####################################################################
      END MODULE ALLFUN
!####################################################################

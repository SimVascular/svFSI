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
         MODULE PROCEDURE DESTROYFACE, DESTROYMSH, DESTROYBC, DESTROYBF,
     2      DESTROYDMN, DESTROYEQ, DESTROYBS, DESTROYMB, DESTROYDATA,
     3      DESTROYADJ, DESTROYSTACK, DESTROYQUEUE, DESTROYTRACE,
     4      DESTROYIBCM, DESTROYFS
      END INTERFACE DESTROY

      INTERFACE GETNADJCNCY
         MODULE PROCEDURE GETNADJ_MSH, GETNADJ_FACE
      END INTERFACE

      INTERFACE GETEADJCNCY
         MODULE PROCEDURE GETEADJ_MSH, GETEADJ_FACE
      END INTERFACE

      INTERFACE IB_SYNCN
         MODULE PROCEDURE IB_SYNCNS, IB_SYNCNV
      END INTERFACE

      INTERFACE IB_SYNCG
         MODULE PROCEDURE IB_SYNCGS, IB_SYNCGV
      END INTERFACE

      CONTAINS
!####################################################################
!     This routine integrate s over the surface faId.
      FUNCTION IntegS(lFa, s, pflag)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: s(:)
      LOGICAL, INTENT(IN), OPTIONAL :: pflag
      REAL(KIND=RKIND) IntegS

      LOGICAL isIB, flag
      INTEGER(KIND=IKIND) a, e, g, Ac, nNo, insd
      REAL(KIND=RKIND) sHat, Jac, n(nsd)
      TYPE(fsType) :: fs

      flag = .FALSE.
      IF (PRESENT(pflag)) flag = pFlag

      insd = nsd - 1
      IF (msh(lFa%iM)%lShl) insd = insd - 1
      IF (msh(lFa%iM)%lFib) insd = 0

      nNo = SIZE(s)
      IF (nNo .NE. tnNo) THEN
         IF (ibFlag) THEN
            IF (nNo .NE. ib%tnNo) err =
     2         "Incompatible vector size in Integ"
         ELSE
            err = "Incompatible vector size in vInteg"
         END IF
      END IF

      isIB = .FALSE.
      IF (ibFlag) THEN
         IF (nNo .EQ. ib%tnNo) isIB = .TRUE.
      END IF

!     Update pressure function space for Taylor-Hood element
      IF (flag) THEN
         IF (lFa%nFs.NE.2) err = "Incompatible boundary "//
     2      "integral function call and face element type"
         fs%nG    = lFa%fs(2)%nG
         fs%eType = lFa%fs(2)%eType
         fs%lShpF = lFa%fs(2)%lShpF
         fs%eNoN  = lFa%fs(2)%eNoN
         ALLOCATE(fs%w(fs%nG), fs%N(fs%eNoN,fs%nG),
     2      fs%Nx(insd,fs%eNoN,fs%nG))
         IF (fs%eType .NE. eType_NRB) THEN
            fs%w  = lFa%fs(2)%w
            fs%N  = lFa%fs(2)%N
            fs%Nx = lFa%fs(2)%Nx
         END IF
      ELSE
         fs%nG    = lFa%fs(1)%nG
         fs%eType = lFa%fs(1)%eType
         fs%lShpF = lFa%fs(1)%lShpF
         fs%eNoN  = lFa%fs(1)%eNoN
         ALLOCATE(fs%w(fs%nG), fs%N(fs%eNoN,fs%nG),
     2      fs%Nx(insd,fs%eNoN,fs%nG))
         IF (fs%eType .NE. eType_NRB) THEN
            fs%w  = lFa%fs(1)%w
            fs%N  = lFa%fs(1)%N
            fs%Nx = lFa%fs(1)%Nx
         END IF
      END IF

      IntegS = 0._RKIND
      DO e=1, lFa%nEl
!     Updating the shape functions, if this is a NURB
         IF (lFa%eType .EQ. eType_NRB) THEN
            IF (.NOT.isIB) THEN
               CALL NRBNNXB(msh(lFa%iM), lFa, e)
            ELSE
               CALL NRBNNXB(ib%msh(lFa%iM), lFa, e)
            END IF
            fs%w  = lFa%w
            fs%N  = lFa%N
            fs%Nx = lFa%Nx
         END IF

         DO g=1, fs%nG
            IF (.NOT.isIB) THEN
               CALL GNNB(lFa, e, g, insd, fs%eNoN, fs%Nx(:,:,g), n)
            ELSE
               CALL GNNIB(lFa, e, g, n)
            END IF
            Jac = SQRT(NORM(n))

!     Calculating the function value
            sHat = 0._RKIND
            DO a=1, fs%eNoN
               Ac   = lFa%IEN(a,e)
               sHat = sHat + s(Ac)*fs%N(a,g)
            END DO
!     Now integrating
            IntegS = IntegS + Jac*fs%w(g)*sHat
         END DO
      END DO

      IF (cm%seq() .OR. isIB) RETURN
      IntegS = cm%reduce(IntegS)

      RETURN
      END FUNCTION IntegS
!--------------------------------------------------------------------
!     This routine integrate s over the surface faId.
      FUNCTION IntegV(lFa, s)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: s(:,:)
      REAL(KIND=RKIND) IntegV

      LOGICAL isIB
      INTEGER(KIND=IKIND) a, i, e, Ac, g, nNo
      REAL(KIND=RKIND) sHat, n(nsd)

      IF (SIZE(s,1) .NE. nsd) err = "Incompatible vector size in IntegV"

      nNo = SIZE(s,2)
      IF (nNo .NE. tnNo) THEN
         IF (ibFlag) THEN
            IF (nNo .NE. ib%tnNo) err =
     2         "Incompatible vector size in IntegV"
         ELSE
            err = "Incompatible vector size in IntegV"
         END IF
      END IF

      isIB = .FALSE.
      IF (ibFlag) THEN
         IF (nNo .EQ. ib%tnNo) isIB = .TRUE.
      END IF

      IntegV = 0._RKIND
      DO e=1, lFa%nEl
!     Updating the shape functions, if this is a NURB
         IF (lFa%eType .EQ. eType_NRB) THEN
            IF (.NOT.isIB) THEN
               CALL NRBNNXB(msh(lFa%iM), lFa, e)
            ELSE
               CALL NRBNNXB(ib%msh(lFa%iM), lFa, e)
            END IF
         END IF

         DO g=1, lFa%nG
            IF (.NOT.isIB) THEN
               CALL GNNB(lFa, e, g, nsd-1, lFa%eNoN, lFa%Nx(:,:,g), n)
            ELSE
               CALL GNNIB(lFa, e, g, n)
            END IF

!     Calculating the function value
            sHat = 0._RKIND
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

      IF (cm%seq() .OR . isIB) RETURN
      IntegV = cm%reduce(IntegV)

      RETURN
      END FUNCTION IntegV
!--------------------------------------------------------------------
!     This routine integrate s(l:u,:) over the surface faId.
      FUNCTION IntegG(lFa, s, l, uo, THflag)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: s(:,:)
      INTEGER(KIND=IKIND), INTENT(IN) :: l
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: uo
      LOGICAL, INTENT(IN), OPTIONAL :: THflag

      LOGICAL flag
      INTEGER(KIND=IKIND) a, u, nNo
      REAL(KIND=RKIND) IntegG
      REAL(KIND=RKIND), ALLOCATABLE :: sclr(:), vec(:,:)

      u = l
      IF (PRESENT(uo)) u = uo

      flag = .FALSE.
      IF (PRESENT(THflag)) flag = THflag

      nNo = SIZE(s,2)
      IF (nNo .NE. tnNo) THEN
         IF (ibFlag) THEN
            IF (nNo .NE. ib%tnNo) err =
     2         "Incompatible vector size in IntegG"
         ELSE
            err = "Incompatible vector size in IntegG"
         END  IF
      END IF

      IntegG = 0._RKIND
      IF (u-l+1 .EQ. nsd) THEN
         ALLOCATE (vec(nsd,nNo))
         DO a=1, nNo
            vec(:,a) = s(l:u,a)
         END DO
         IntegG = IntegV(lFa,vec)
      ELSE IF (l .EQ. u) THEN
         ALLOCATE (sclr(nNo))
         DO a=1, nNo
            sclr(a) = s(l,a)
         END DO
         IntegG = IntegS(lFa,sclr,flag)
      ELSE
         err = "Unexpected dof in IntegG"
      END IF

      RETURN
      END FUNCTION IntegG
!--------------------------------------------------------------------
!     This routine integrate an equation over a particular domain
      FUNCTION vInteg(dId, s, l, u, pFlag)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: s(:,:)
      INTEGER(KIND=IKIND), INTENT(IN) :: dId
      INTEGER(KIND=IKIND), INTENT(IN) :: l
      INTEGER(KIND=IKIND), INTENT(IN) :: u
      LOGICAL, INTENT(IN), OPTIONAL :: pFlag
      REAL(KIND=RKIND) vInteg

      LOGICAL isIB, flag
      INTEGER(KIND=IKIND) a, e, g, Ac, iM, eNoN, insd, ibl, nNo
      REAL(KIND=RKIND) Jac, nV(nsd), sHat, tmp(nsd,nsd)
      TYPE(fsType) :: fs

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), sl(:), Nxi(:,:),
     2   Nx(:,:), tmps(:,:)

      nNo = SIZE(s,2)
      IF (nNo .NE. tnNo) THEN
         IF (ibFlag) THEN
            IF (nNo .NE. ib%tnNo) err =
     2         "Incompatible vector size in vInteg"
         ELSE
            err = "Incompatible vector size in vInteg"
         END IF
      END IF

      flag = .FALSE.
      IF (PRESENT(pFlag)) THEN
         flag = pFlag
         IF (l .NE. u) err = "Incompatible spatial output setting "//
     2      "and element type"
      END IF

      isIB = .FALSE.
      IF (ibFlag) THEN
         IF (nNo .EQ. ib%tnNo) isIB = .TRUE.
      END IF

      vInteg = 0._RKIND
      IF (.NOT.isIB) THEN
         DO iM=1, nMsh
            insd = nsd
            IF (msh(iM)%lShl) insd = nsd-1
            IF (msh(iM)%lFib) insd = 1

!           Update pressure function space for Taylor-Hood type element
            IF (flag .AND. (msh(iM)%nFs.EQ.2)) THEN
               fs%nG    = msh(iM)%fs(2)%nG
               fs%eType = msh(iM)%fs(2)%eType
               fs%lShpF = msh(iM)%fs(2)%lShpF
               fs%eNoN  = msh(iM)%fs(2)%eNoN
               ALLOCATE(fs%w(fs%nG), fs%N(fs%eNoN,fs%nG),
     2            fs%Nx(nsd,fs%eNoN,fs%nG))
               IF (fs%eType .NE. eType_NRB) THEN
                  fs%w  = msh(iM)%fs(2)%w
                  fs%N  = msh(iM)%fs(2)%N
                  fs%Nx = msh(iM)%fs(2)%Nx
               END IF
            ELSE
               fs%nG    = msh(iM)%fs(1)%nG
               fs%eType = msh(iM)%fs(1)%eType
               fs%lShpF = msh(iM)%fs(1)%lShpF
               fs%eNoN  = msh(iM)%fs(1)%eNoN
               ALLOCATE(fs%w(fs%nG), fs%N(fs%eNoN,fs%nG),
     2            fs%Nx(nsd,fs%eNoN,fs%nG))
               IF (fs%eType .NE. eType_NRB) THEN
                  fs%w  = msh(iM)%fs(1)%w
                  fs%N  = msh(iM)%fs(1)%N
                  fs%Nx = msh(iM)%fs(1)%Nx
               END IF
            END IF
            eNoN = fs%eNoN

            ALLOCATE(xl(nsd,eNoN), Nxi(insd,eNoN), Nx(insd,eNoN),
     2         sl(eNoN), tmps(nsd,insd))

            DO e=1, msh(iM)%nEl
               IF (dId.GT.0 .AND. ALLOCATED(msh(iM)%eId)) THEN
                  IF (.NOT.BTEST(msh(iM)%eId(e),dId)) CYCLE
               END IF
!           Updating the shape functions, if this is a NURB
               IF (msh(iM)%eType .EQ. eType_NRB) THEN
                  CALL NRBNNX(msh(iM), e)
                  fs%w  = msh(iM)%w
                  fs%N  = msh(iM)%N
                  fs%Nx = msh(iM)%Nx
               END IF

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

               DO g=1, fs%nG
                  Nxi(:,:) = fs%Nx(:,:,g)
                  IF (g.EQ.1 .OR. .NOT.fs%lShpF) THEN
                     IF (msh(iM)%lShl) THEN
                        CALL GNNS(eNoN, Nxi, xl, nV, tmps, tmps)
                        Jac = SQRT(NORM(nV))
                     ELSE
                        CALL GNN(eNoN, insd, Nxi, xl, Nx, Jac, tmp)
                     END IF
                  END IF
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

                  sHat = 0._RKIND
                  DO a=1, eNoN
                     Ac = msh(iM)%IEN(a,e)
                     sHat = sHat + sl(a)*fs%N(a,g)
                  END DO
                  vInteg = vInteg + fs%w(g)*Jac*sHat
               END DO
            END DO

            DEALLOCATE(xl, Nxi, Nx, sl, tmps)
            CALL DESTROY(fs)
         END DO
      ELSE
         DO iM=1, ib%nMsh
            eNoN = ib%msh(iM)%eNoN
            insd = nsd

            ALLOCATE(xl(nsd,eNoN), Nxi(insd,eNoN), Nx(insd,eNoN),
     2         sl(eNoN), tmps(nsd,insd))

            DO e=1, ib%msh(iM)%nEl
               IF (dId.GT.0 .AND. ALLOCATED(ib%msh(iM)%eId)) THEN
                  IF (.NOT.BTEST(ib%msh(iM)%eId(e),dId)) CYCLE
               END IF

!           Updating the shape functions, if this is a NURB
               IF (ib%msh(iM)%eType .EQ. eType_NRB)
     2            CALL NRBNNX(ib%msh(iM), e)
               DO a=1, eNoN
                  Ac      = ib%msh(iM)%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Ubo(:,Ac)
                  IF (l .EQ. u) THEN
                     sl(a) = s(l,Ac)
                   ELSE
                     sl(a) = SQRT(NORM(s(l:u,Ac)))
                  END IF
               END DO

               DO g=1, ib%msh(iM)%nG
                  Nxi(:,:) = ib%msh(iM)%Nx(:,:,g)
                  IF (g.EQ.1 .OR. .NOT.ib%msh(iM)%lShpF) THEN
                     CALL GNN(eNoN, insd, Nxi, xl, Nx, Jac, tmp)
                  END IF
                  IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

                  sHat = 0._RKIND
                  DO a=1, eNoN
                     Ac = ib%msh(iM)%IEN(a,e)
                     sHat = sHat + sl(a)*ib%msh(iM)%N(a,g)
                  END DO
                  vInteg = vInteg + ib%msh(iM)%w(g)*Jac*sHat
               END DO
            END DO

            DEALLOCATE(xl, Nxi, Nx, sl, tmps)
         END DO
      END IF

      IF (cm%seq() .OR. isIB) RETURN
      vInteg = cm%reduce(vInteg)

      RETURN
      END FUNCTION vInteg
!####################################################################
!     This a both way communication with three main part:
      SUBROUTINE COMMUS(U)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: U(:)

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
      REAL(KIND=RKIND), INTENT(INOUT) :: U(:,:)

      INTEGER(KIND=IKIND) m

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
      REAL(KIND=RKIND), INTENT(IN) :: U(:)
      REAL(KIND=RKIND), ALLOCATABLE :: MKCS(:)

      INTEGER(KIND=IKIND) a

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
      REAL(KIND=RKIND), INTENT(IN) :: U(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: MKCV(:,:)

      INTEGER(KIND=IKIND) m, a

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
      REAL(KIND=RKIND), INTENT(INOUT) :: U(:)

      INTEGER(KIND=IKIND) a
      REAL(KIND=RKIND), ALLOCATABLE :: tmp(:)

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
      REAL(KIND=RKIND), INTENT(INOUT) :: U(:,:)

      INTEGER(KIND=IKIND) m, a
      REAL(KIND=RKIND), ALLOCATABLE :: tmp(:,:)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: U(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: GLOBALIS(:)

      INTEGER(KIND=IKIND) Ac, a, e, i, ierr
      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disp(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ienU(:,:), gienU(:,:)

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
      REAL(KIND=RKIND), INTENT(IN) :: U(:)
      REAL(KIND=RKIND), ALLOCATABLE :: GLOBALRS(:)

      INTEGER(KIND=IKIND) Ac, a, e, i, ierr
      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: ienU(:,:), gienU(:,:)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: U(:,:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: GLOBALIV(:,:)

      INTEGER(KIND=IKIND) m, Ac, a, e, i, ierr
      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disp(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: ienU(:,:), gienU(:,:)

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
      REAL(KIND=RKIND), INTENT(IN) :: U(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: GLOBALRV(:,:)

      INTEGER(KIND=IKIND) m, Ac, a, e, i, ierr
      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: ienU(:,:), gienU(:,:)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: U(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: LOCALIS(:)

      INTEGER(KIND=IKIND) Ac, a
      INTEGER(KIND=IKIND), ALLOCATABLE :: tmpU(:)

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
      REAL(KIND=RKIND), INTENT(IN) :: U(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: LOCALRV(:,:)

      INTEGER(KIND=IKIND) m, s, e, a, Ac
      REAL(KIND=RKIND), ALLOCATABLE :: tmpU(:)

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, e
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND) DOMAIN

      INTEGER(KIND=IKIND) iDmn

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
      FUNCTION ISDOMAIN(iEq, node, phys)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, node, phys
      LOGICAL ISDOMAIN

      INTEGER(KIND=IKIND) iDmn

      ISDOMAIN = .FALSE.
      IF (eq(iEq)%nDmn .EQ. 1) THEN
!        Single domain is assumed and we only need to check that
         IF (eq(iEq)%dmn(1)%phys .EQ. phys) ISDOMAIN = .TRUE.
      ELSE
!        Domain partition is expected
         IF (.NOT.ALLOCATED(dmnId)) err = "Domain partitioning info "//
     2      "is not provided"
         DO iDmn=1, eq(iEq)%nDmn
            IF (eq(iEq)%dmn(iDmn)%phys .EQ. phys) THEN
               IF (BTEST(dmnId(node),eq(iEq)%dmn(iDmn)%Id)) THEN
                  ISDOMAIN = .TRUE.
                  RETURN
               END IF
            END IF
         END DO
      END IF

      RETURN
      END FUNCTION ISDOMAIN
!####################################################################
!     This function returns the IB domain that an element of an IB mesh
!     belongs to
      FUNCTION IB_DOMAIN(lM, iEq, e)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, e
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND) IB_DOMAIN

      INTEGER(KIND=IKIND) iDmn

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
      INTEGER(KIND=IKIND), INTENT(IN) :: iEq, node, phys
      LOGICAL IB_ISDOMAIN

      INTEGER(KIND=IKIND) iDmn

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
      INTEGER(KIND=IKIND), INTENT(IN) :: fid
      CHARACTER(LEN=*), INTENT(IN) :: kwd
      INTEGER(KIND=IKIND), INTENT(OUT) :: n

      INTEGER(KIND=IKIND) l
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
      PURE SUBROUTINE DESTROYFS(fs)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fsType), INTENT(OUT) :: fs

      fs%eType = eType_NA

      IF (ALLOCATED(fs%w))      DEALLOCATE(fs%w)
      IF (ALLOCATED(fs%xi))     DEALLOCATE(fs%xi)
      IF (ALLOCATED(fs%xib))    DEALLOCATE(fs%xib)
      IF (ALLOCATED(fs%N))      DEALLOCATE(fs%N)
      IF (ALLOCATED(fs%Nb))     DEALLOCATE(fs%Nb)
      IF (ALLOCATED(fs%Nx))     DEALLOCATE(fs%Nx)
      IF (ALLOCATED(fs%Nxx))    DEALLOCATE(fs%Nxx)

      RETURN
      END SUBROUTINE DESTROYFS
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYFACE(lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(OUT) :: lFa

      INTEGER(KIND=IKIND) i

      IF (ALLOCATED(lFa%gE))     DEALLOCATE(lFa%gE)
      IF (ALLOCATED(lFa%gN))     DEALLOCATE(lFa%gN)
      IF (ALLOCATED(lFa%lN))     DEALLOCATE(lFa%lN)
      IF (ALLOCATED(lFa%IEN))    DEALLOCATE(lFa%IEN)
      IF (ALLOCATED(lFa%gebc))   DEALLOCATE(lFa%gebc)
      IF (ALLOCATED(lFa%w))      DEALLOCATE(lFa%w)
      IF (ALLOCATED(lFa%x))      DEALLOCATE(lFa%x)
      IF (ALLOCATED(lFa%xi))     DEALLOCATE(lFa%xi)
      IF (ALLOCATED(lFa%N))      DEALLOCATE(lFa%N)
      IF (ALLOCATED(lFa%nV))     DEALLOCATE(lFa%nV)
      IF (ALLOCATED(lFa%Nx))     DEALLOCATE(lFa%Nx)
      IF (ALLOCATED(lFa%Nxx))    DEALLOCATE(lFa%Nxx)

      CALL DESTROYADJ(lFa%nAdj)
      CALL DESTROYADJ(lFa%eAdj)

      IF (ALLOCATED(lFa%fs)) THEN
         DO i=1, lFa%nFs
            CALL DESTROY(lFa%fs(i))
         END DO
         DEALLOCATE(lFa%fs)
      END IF

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

      INTEGER(KIND=IKIND) iFa, i, insd

      insd = nsd
      IF (lM%lShl) insd = nsd - 1

      IF (ALLOCATED(lM%eDist))   DEALLOCATE(lM%eDist)
      IF (ALLOCATED(lM%eId))     DEALLOCATE(lM%eId)
      IF (ALLOCATED(lM%gN))      DEALLOCATE(lM%gN)
      IF (ALLOCATED(lM%gpN))     DEALLOCATE(lM%gpN)
      IF (ALLOCATED(lM%gIEN))    DEALLOCATE(lM%gIEN)
      IF (ALLOCATED(lM%IEN))     DEALLOCATE(lM%IEN)
      IF (ALLOCATED(lM%otnIEN))  DEALLOCATE(lM%otnIEN)
      IF (ALLOCATED(lM%INN))     DEALLOCATE(lM%INN)
      IF (ALLOCATED(lM%lN))      DEALLOCATE(lM%lN)
      IF (ALLOCATED(lM%eIEN))    DEALLOCATE(lM%eIEN)
      IF (ALLOCATED(lM%sbc))     DEALLOCATE(lM%sbc)
      IF (ALLOCATED(lM%iGC))     DEALLOCATE(lM%iGC)
      IF (ALLOCATED(lM%nW))      DEALLOCATE(lM%nW)
      IF (ALLOCATED(lM%w))       DEALLOCATE(lM%w)
      IF (ALLOCATED(lM%xib))     DEALLOCATE(lM%xib)
      IF (ALLOCATED(lM%xi))      DEALLOCATE(lM%xi)
      IF (ALLOCATED(lM%x))       DEALLOCATE(lM%x)
      IF (ALLOCATED(lM%N))       DEALLOCATE(lM%N)
      IF (ALLOCATED(lM%Nb))      DEALLOCATE(lM%Nb)
      IF (ALLOCATED(lM%nV))      DEALLOCATE(lM%nV)
      IF (ALLOCATED(lM%fN))      DEALLOCATE(lM%fN)
      IF (ALLOCATED(lM%Nx))      DEALLOCATE(lM%Nx)
      IF (ALLOCATED(lM%Nxx))     DEALLOCATE(lM%Nxx)

      IF (ALLOCATED(lM%fs)) THEN
         DO i=1, lM%nFs
            CALL DESTROY(lM%fs(i))
         END DO
         DEALLOCATE(lM%fs)
      END IF

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

      lM%lShl  = .FALSE.
      lM%eType = eType_NA
      lM%gnEl  = 0
      lM%gnNo  = 0
      lM%nEl   = 0
      lM%nFa   = 0
      lM%nNo   = 0
      lM%nFn   = 0

      RETURN
      END SUBROUTINE DESTROYMSH
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYBC(lBc)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(OUT) :: lBc

      IF (ALLOCATED(lBc%eDrn)) DEALLOCATE(lBc%eDrn)
      IF (ALLOCATED(lBc%h))    DEALLOCATE(lBc%h)
      IF (ALLOCATED(lBc%gx))   DEALLOCATE(lBc%gx)

      IF (ALLOCATED(lBc%gt)) THEN
         IF (ALLOCATED(lBc%gt%qi)) DEALLOCATE(lBc%gt%qi)
         IF (ALLOCATED(lBc%gt%qs)) DEALLOCATE(lBc%gt%qs)
         IF (ALLOCATED(lBc%gt%r))  DEALLOCATE(lBc%gt%r)
         IF (ALLOCATED(lBc%gt%i))  DEALLOCATE(lBc%gt%i)
         DEALLOCATE(lBc%gt)
      END IF

      IF (ALLOCATED(lBc%gm)) THEN
         CALL DESTROYMB(lBc%gm)
         DEALLOCATE(lBc%gm)
      END IF

      lBc%weakDir  = .FALSE.
      lBc%flwP     = .FALSE.
      lBc%bType    = 0
      lBc%cplBCptr = 0
      lBc%g        = 0._RKIND
      lBc%r        = 0._RKIND
      lBc%k        = 0._RKIND
      lBc%k        = 0._RKIND
      lBc%tauB     = 0._RKIND

      RETURN
      END SUBROUTINE DESTROYBC
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYBF(lBf)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bfType), INTENT(OUT) :: lBf

      IF (ALLOCATED(lBf%b))  DEALLOCATE(lBf%b)
      IF (ALLOCATED(lBf%bx)) DEALLOCATE(lBf%bx)
      IF (ALLOCATED(lBf%bt)) THEN
         IF (ALLOCATED(lBf%bt%qi)) DEALLOCATE(lBf%bt%qi)
         IF (ALLOCATED(lBf%bt%qs)) DEALLOCATE(lBf%bt%qs)
         IF (ALLOCATED(lBf%bt%r))  DEALLOCATE(lBf%bt%r)
         IF (ALLOCATED(lBf%bt%i))  DEALLOCATE(lBf%bt%i)
         DEALLOCATE(lBf%bt)
      END IF
      IF (ALLOCATED(lBf%bm)) THEN
         CALL DESTROYMB(lBf%bm)
         DEALLOCATE(lBf%bm)
      END IF
      lBf%bType = 0

      RETURN
      END SUBROUTINE DESTROYBF
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYDMN(lDmn)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(OUT) :: lDmn

      lDmn%Id   = -1
      lDmn%v    = 0._RKIND
      lDmn%prop = 0._RKIND

      ! lDmn%stM
      lDmn%stM%volType = stVol_NA
      lDmn%stM%Kpen    = 0._RKIND
      lDmn%stM%isoType = stIso_NA
      lDmn%stM%C10     = 0._RKIND
      lDmn%stM%C01     = 0._RKIND
      lDmn%stM%a       = 0._RKIND
      lDmn%stM%b       = 0._RKIND
      lDmn%stM%aff     = 0._RKIND
      lDmn%stM%bff     = 0._RKIND
      lDmn%stM%ass     = 0._RKIND
      lDmn%stM%bss     = 0._RKIND
      lDmn%stM%afs     = 0._RKIND
      lDmn%stM%bfs     = 0._RKIND
      lDmn%stM%khs     = 100._RKIND

      lDmn%stM%Tf%g     = 0._RKIND
      lDmn%stM%Tf%fType = 0
      IF (ALLOCATED(lDmn%stM%Tf%gt%qi)) DEALLOCATE(lDmn%stM%Tf%gt%qi)
      IF (ALLOCATED(lDmn%stM%Tf%gt%qs)) DEALLOCATE(lDmn%stM%Tf%gt%qs)
      IF (ALLOCATED(lDmn%stM%Tf%gt%r))  DEALLOCATE(lDmn%stM%Tf%gt%r)
      IF (ALLOCATED(lDmn%stM%Tf%gt%i))  DEALLOCATE(lDmn%stM%Tf%gt%i)

      ! lDmn%cep
      lDmn%cep%cepType  = cepModel_NA
      lDmn%cep%Diso     = 0._RKIND
      IF (ALLOCATED(lDmn%cep%Dani)) DEALLOCATE(lDmn%cep%Dani)

      lDmn%cep%Istim%Ts = 0._RKIND
      lDmn%cep%Istim%Td = 0._RKIND
      lDmn%cep%Istim%CL = 0._RKIND
      lDmn%cep%Istim%A  = 0._RKIND

      lDmn%cep%odeS%tIntType = tIntType_NA
      lDmn%cep%odeS%maxItr   = 5
      lDmn%cep%odeS%absTol   = 1.E-8_RKIND
      lDmn%cep%odeS%relTol   = 1.E-4_RKIND

      RETURN
      END SUBROUTINE DESTROYDMN
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYEQ(lEq)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(OUT) :: lEq

      INTEGER(KIND=IKIND) iBc, iBf, iDmn

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
      IF (ALLOCATED(lEq%dmn)) THEN
         DO iDmn=1, lEq%nDmn
            CALL DESTROYDMN(lEq%dmn(iDmn))
         END DO
      END IF
      IF (ALLOCATED(lEq%dmnIB)) THEN
         DO iDmn=1, lEq%nDmnIB
            CALL DESTROYDMN(lEq%dmnIB(iDmn))
         END DO
      END IF
      IF (ALLOCATED(lEq%bf)) THEN
         DO iBf=1, lEq%nBf
            CALL DESTROYBF(lEq%bf(iBf))
         END DO
         DEALLOCATE(lEq%bf)
      END IF
      IF (ALLOCATED(lEq%dmn))      DEALLOCATE(lEq%dmn)
      IF (ALLOCATED(lEq%dmnIB))    DEALLOCATE(lEq%dmnIB)
      IF (ALLOCATED(lEq%output))   DEALLOCATE(lEq%output)
      IF (ALLOCATED(lEq%outIB))    DEALLOCATE(lEq%outIB)

      lEq%coupled = .TRUE.
      lEq%dof     = 0
      lEq%maxItr  = 5
      lEq%minItr  = 1
      lEq%nOutput = 0
      lEq%nOutIB  = 0
      lEq%nDmn    = 0
      lEq%nDmnIB  = 0
      lEq%nBc     = 0
      lEq%nBcIB   = 0
      lEq%nBf     = 0
      lEq%useTLS  = .FALSE.
      lEq%assmTLS = .FALSE.

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

      IF (ALLOCATED(trc%gN))    DEALLOCATE(trc%gN)
      IF (ALLOCATED(trc%gE))    DEALLOCATE(trc%gE)
      IF (ALLOCATED(trc%nptr))  DEALLOCATE(trc%nptr)
      IF (ALLOCATED(trc%gptr))  DEALLOCATE(trc%gptr)
      IF (ALLOCATED(trc%xi))    DEALLOCATE(trc%xi)
      IF (ALLOCATED(trc%xiG))   DEALLOCATE(trc%xiG)

      trc%n  = 0
      trc%nG = 0

      RETURN
      END SUBROUTINE DESTROYTRACE
!--------------------------------------------------------------------
      PURE SUBROUTINE DESTROYIBCM(ibCm)
      USE COMMOD
      IMPLICIT NONE
      TYPE(ibCommType), INTENT(OUT) :: ibCm

      IF (ALLOCATED(ibCm%n))  DEALLOCATE(ibCm%n)
      IF (ALLOCATED(ibCm%nG)) DEALLOCATE(ibCm%nG)
      IF (ALLOCATED(ibCm%gN)) DEALLOCATE(ibCm%gN)
      IF (ALLOCATED(ibCm%gE)) DEALLOCATE(ibCm%gE)

      RETURN
      END SUBROUTINE DESTROYIBCM
!####################################################################
!     Spliting "m" jobs between "n" workers. "b" contains amount of jobs
!     and "A" will store the distribution of jobs
      RECURSIVE SUBROUTINE SPLITJOBS(m,n,A,b)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: m, n
      REAL(KIND=RKIND), INTENT(IN) :: b(m)
      REAL(KIND=RKIND), INTENT(OUT) :: A(m,n)

      INTEGER(KIND=IKIND) i, j, k, ml, nl, mr, nr
      REAL(KIND=RKIND) sb, sbl, sl, optsl
      REAL(KIND=RKIND), ALLOCATABLE :: Al(:,:), Ar(:,:), bl(:), br(:)

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
            A(i,j) = b(i)/REAL(n, KIND=RKIND)
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
      sbl = sb*REAL(nl, KIND=RKIND)/REAL(n, KIND=RKIND)

      sl = 0._RKIND
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
      nl = NINT(REAL(n, KIND=RKIND)*SUM(bl)/sb, KIND=IKIND)
      IF (nl .EQ. 0) THEN
         nl = 1
         nr = n - 1
      ELSE IF (nl .EQ. n) THEN
         nl = n - 1
         nr = 1
      ELSE
         nr = n - nl
      END IF
      ALLOCATE (Al(ml,nl), Ar(mr,nr))

      CALL SPLITJOBS(ml,nl,Al,bl)
      CALL SPLITJOBS(mr,nr,Ar,br)

      A = 0._RKIND
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
      INTEGER(KIND=IKIND), INTENT(IN) :: iDmn
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: ifirst, ilast

      INTEGER(KIND=IKIND) e, first, last

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
!     Finding the mesh ID based on the mesh name
      SUBROUTINE FINDMSH(mshName, iM)
      USE COMMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: mshName
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM

      DO iM=1, nMsh
         IF (msh(iM)%name .EQ. mshName) EXIT
      END DO
      IF (iM .GT. nMsh) err="Unable to find msh <"//TRIM(mshName)//">"

      RETURN
      END SUBROUTINE FINDMSH
!--------------------------------------------------------------------
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE FINDFACE(faceName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE
      CHARACTER(LEN=stdL) :: faceName
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM, iFa

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
      INTEGER(KIND=IKIND), INTENT(IN) :: nDim, eNoN
      REAL(KIND=RKIND), INTENT(IN) :: x(nDim,eNoN), Nxi(nDim,eNoN)

      INTEGER(KIND=IKIND) :: a
      REAL(KIND=RKIND) :: Jac, xXi(nDim,nDim)

      xXi = 0._RKIND
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
      INTEGER(KIND=IKIND), INTENT(IN) :: nDim, eNoN
      REAL(KIND=RKIND), INTENT(IN) :: x(nDim,eNoN)

      INTEGER(KIND=IKIND) :: i,j,a,cnt
      REAL(KIND=RKIND) :: circumRad, integ_eq, integ_el
      REAL(KIND=RKIND) :: SkF, Dmat(eNoN,nDim+2), Dsub(eNoN,nDim+1)
      REAL(KIND=RKIND), DIMENSION(nDim+2) :: coeff,detD

      IF (nDim .EQ. 2) THEN
         coeff = (/1._RKIND, -1._RKIND, 1._RKIND, -1._RKIND/)
      ELSE
         coeff = (/1._RKIND, 1._RKIND, -1._RKIND, 1._RKIND, 1._RKIND/)
      END IF

      Dmat(:,:) = 1._RKIND
      Dsub(:,:) = 0._RKIND

      DO a=1, eNoN
         Dmat(a,1) = SUM(x(:,a)**2._RKIND)
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

      circumRad = SQRT(SUM(detD(2:nDim+1)**2._RKIND) -
     2   4._RKIND*detD(1)*detD(nDim+2)) /(2._RKIND*ABS(detD(1)))

      IF (nDim .EQ. 2) THEN
         integ_eq = 0.25_RKIND*SQRT(27._RKIND)*circumRad**2._RKIND
         integ_el = 0.5_RKIND*ABS(detD(1))
      ELSE IF (nDim .EQ. 3) THEN
         integ_eq = 80._RKIND*circumRad**3._RKIND /SQRT(243._RKIND)
         integ_el = ABS(detD(1))/6._RKIND
      END IF

      SkF = ABS(integ_eq - integ_el) / integ_eq

      RETURN
      END FUNCTION SKEWNESS
!--------------------------------------------------------------------
!     Computes the Aspect Ratio of an element
      FUNCTION ASPECTRATIO(nDim, eNoN, x) result(AR)
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nDim,eNoN
      REAL(KIND=RKIND), INTENT(IN) :: x(nDim,eNoN)

      INTEGER(KIND=IKIND) :: i,j,ap,a,b,icol,irow
      INTEGER(KIND=IKIND) :: rowM(eNoN,eNoN-1), colM(nDim,nDim-1)
      REAL(KIND=RKIND) :: AR, s(eNoN)
      REAL(KIND=RKIND) :: Dsub(nDim,nDim), detD(nDim)

      IF (nDim .EQ. 2) THEN
         s(:) = 0._RKIND
         DO a=1, eNoN
            ap = a+1
            IF (a .EQ. eNoN) ap = 1
            s(a) = SQRT(SUM((x(:,a)-x(:,ap))**2._RKIND))
         END DO
      ELSE IF (nDim .EQ. 3) THEN
         rowM = reshape( (/1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4/),
     2                   shape(rowM) )
         colM = reshape( (/1, 2, 3, 2, 3, 1/), shape(colM) )

         DO a=1, eNoN
            DO b=1, nDim
               Dsub(:,:) = 1._RKIND
               DO i=1, eNoN-1
                  irow = rowM(a,i)
                  DO j=1, nDim-1
                     icol = colM(b,j)
                     Dsub(i,j) = x(icol,irow)
                  END DO ! j
               END DO ! i

               detD(b) = MAT_DET(Dsub,nDim)
            END DO ! b
            s(a) = 0.5_RKIND*SQRT(SUM(detD(:)**2._RKIND))
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

      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, i, j, k, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjL(:,:)

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
               j = 1
               DO i=1, incNd(Ac)
                  IF (Bc .EQ. adjL(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  ELSE IF (Bc .GT. adjL(i,Ac)) THEN
                     j = i + 1
                  END IF
               END DO
               IF (flag) THEN
                  IF (incNd(Ac) .EQ. 0) THEN
                     incNd(Ac)  = 1
                     adjL(1,Ac) = Bc
                  ELSE
                     DO k=incNd(Ac), j, -1
                        adjL(k+1,Ac) = adjL(k,Ac)
                     END DO
                     adjL(j,Ac) = Bc
                     incNd(Ac) = incNd(Ac) + 1
                  END IF
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

      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, i, j, k, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjL(:,:)

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
               j = 1
               DO i=1, incNd(Ac)
                  IF (Bc .EQ. adjL(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  ELSE IF (Bc .GT. adjL(i,Ac)) THEN
                     j = i + 1
                  END IF
               END DO
               IF (flag) THEN
                  IF (incNd(Ac) .EQ. 0) THEN
                     incNd(Ac)  = 1
                     adjL(1,Ac) = Bc
                  ELSE
                     DO k=incNd(Ac), j, -1
                        adjL(k+1,Ac) = adjL(k,Ac)
                     END DO
                     adjL(j,Ac) = Bc
                     incNd(Ac)  = incNd(Ac) + 1
                  END IF
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

      INTEGER(KIND=IKIND) :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjList(:,:),
     2   tmpI(:,:)

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

      INTEGER(KIND=IKIND) :: a, b, e, Ac, i, j, maxAdj
      LOGICAL :: flag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), adjList(:,:),
     2   tmpI(:,:)

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
      SUBROUTINE FINDE(xp, lM, xg, Dg, nNo, ne, eList, Ec, xi, lDebug)
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, ne, eList(ne)
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd), xg(nsd,nNo), Dg(nsd,nNo)
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
      REAL(KIND=RKIND), INTENT(OUT) :: xi(nsd)
      LOGICAL, INTENT(IN), OPTIONAL :: lDebug

      LOGICAL :: ldbg, l1, l2, l3, l4
      INTEGER(KIND=IKIND) :: a, e, i, Ac, eNoN
      REAL(KIND=RKIND) :: rt, xi0(nsd)

      LOGICAL, ALLOCATABLE :: eChck(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), N(:), Nxi(:,:)

      IF (lM%lShl) err = " Finding traces is not applicable for shells"

      ldbg = .FALSE.
      IF (PRESENT(lDebug)) ldbg = lDebug
      Ec   = 0
      eNoN = lM%eNoN
      ALLOCATE(eChck(lM%nEl), xl(nsd,eNoN), N(eNoN), Nxi(nsd,eNoN))

!     Initialize parameteric coordinate for Newton's iterations
      xi0 = 0._RKIND
      DO i=1, lM%nG
         xi0 = xi0 + lM%xi(:,i)
      END DO
      xi0 = xi0 / REAL(lM%nG, KIND=RKIND)

      IF (ldbg) THEN
         WRITE(1000+cm%tF(),'(A)') "=================================="
         WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') "Finding trace for xp: "
         DO i=1,  nsd
            WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xp(i))
         END DO
         WRITE(1000+cm%tF(),'(A)')
         WRITE(1000+cm%tF(),'(A)')
         WRITE(1000+cm%tF(),'(A)') "List of elems: <"//STR(ne)//">"
         DO e=1, ne
            WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(eList(e))
         END DO
         WRITE(1000+cm%tF(),'(A)')
         WRITE(1000+cm%tF(),'(A)') "=================================="
      END IF

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

         IF (ldbg) THEN
            WRITE(1000+cm%tF(),'(A)') "----------------------------"
            WRITE(1000+cm%tF(),'(A)') "Probe el: "//STR(Ec)
            DO a=1, eNoN
               Ac = lM%IEN(a,Ec)
               WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') STR(Ac)
               DO i=1, nsd
                  WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//
     2               STR(xl(i,a))
               END DO
               WRITE(1000+cm%tF(),'(A)')
            END DO
         END IF

         xi = xi0
         CALL GETXI(lM%eType, eNoN, xl, xp, xi, l1)

         IF (ldbg) THEN
            WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') "xi: "
            DO i=1, nsd
               WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(xi(i))
            END DO
            WRITE(1000+cm%tF(),'(A)')
         END IF

!        Check if parameteric coordinate is within bounds
         a = 0
         DO i=1, nsd
            IF (xi(i).GE.lM%xib(1,i) .AND. xi(i).LE.lM%xib(2,i))
     2         a = a + 1
         END DO
         l2 = a .EQ. nsd

!        Check for shape function even if the Newton's method returns OK
         CALL GETGNN(nsd, lM%eType, eNoN, xi, N, Nxi)

         IF (ldbg) THEN
            WRITE(1000+cm%tF(),'(4X,A)',ADVANCE='NO') "N: "
            DO a=1, eNoN
               WRITE(1000+cm%tF(),'(A)',ADVANCE='NO') " "//STR(N(a))
            END DO
            WRITE(1000+cm%tF(),'(A)')
            CALL FLUSH(1000+cm%tF())
         END IF

!        Check if shape functions are within bounds and sum to unity
         i  = 0
         rt = 0._RKIND
         DO a=1, eNoN
            rt = rt + N(a)
            IF (N(a).GT.lM%Nb(1,a) .AND. N(a).LT.lM%Nb(2,a)) i = i + 1
         END DO
         l3 = i .EQ. eNoN
         l4 = rt.GE.0.9999_RKIND .AND. rt.LE.1.0001_RKIND

         l1 = ALL((/l1, l2, l3, l4/))
         IF (l1) RETURN
      END DO

      Ec = 0
      xi = 0._RKIND

      RETURN
      END SUBROUTINE FINDE
!####################################################################
      SUBROUTINE IB_SYNCNS(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(INOUT) :: U(:)

      INTEGER(KIND=IKIND) i, a, Ac, iM, nl, ng, ierr, tag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      IF (SIZE(U,1) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%n)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         DO i=1, ib%msh(iM)%trc%n
            a  = ib%msh(iM)%trc%gN(i)
            Ac = ib%msh(iM)%gN(a)
            incNd(Ac) = 1
         END DO
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
         gU   = 0._RKIND
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

         U = 0._RKIND
         DO a=1, ng
            Ac    = ib%cm%gN(a)
            U(Ac) = U(Ac) + gU(a)
         END DO
      ELSE IF (nl .NE. 0) THEN
         tag = cm%tF() * 100
         CALL MPI_SEND(lU, nl, mpreal, master, tag, cm%com(), ierr)
      END IF

      DEALLOCATE(lU, gU)

      CALL cm%bcast(U)

      RETURN
      END SUBROUTINE IB_SYNCNS
!--------------------------------------------------------------------
      SUBROUTINE IB_SYNCNV(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(INOUT) :: U(:,:)

      INTEGER(KIND=IKIND) m, i, a, e, s, Ac, iM, nl, ng, ierr, tag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%n)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         DO i=1, ib%msh(iM)%trc%n
            a  = ib%msh(iM)%trc%gN(i)
            Ac = ib%msh(iM)%gN(a)
            incNd(Ac) = 1
         END DO
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
         gU   = 0._RKIND
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

         U = 0._RKIND
         DO a=1, ng
            Ac = ib%cm%gN(a)
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
      END SUBROUTINE IB_SYNCNV
!####################################################################
      SUBROUTINE IB_SYNCGS(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(INOUT) :: U(:)

      INTEGER(KIND=IKIND) i, a, e, Ac, iM, nl, ng, ierr, tag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      IF (SIZE(U,1) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%nG)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         DO i=1, ib%msh(iM)%trc%nG
            e = ib%msh(iM)%trc%gE(1,i)
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               incNd(Ac) = 1
            END DO
         END DO
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

      ng = SUM(ib%cm%nG)
      ALLOCATE(gU(ng))
      IF (cm%mas()) THEN
         ALLOCATE(rReq(cm%np()))
         rReq = 0
         gU   = 0._RKIND
         ng   = 0
         DO i=1, cm%np()
            nl = ib%cm%nG(i)
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
            IF (i.EQ.1 .OR. ib%cm%nG(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO

         U = 0._RKIND
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
      END SUBROUTINE IB_SYNCGS
!--------------------------------------------------------------------
      SUBROUTINE IB_SYNCGV(U)
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND), INTENT(INOUT) :: U(:,:)

      INTEGER(KIND=IKIND) m, i, a, e, s, Ac, iM, nl, ng, ierr, tag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), rReq(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lU(:), gU(:)

      IF (cm%seq()) RETURN

      m = SIZE(U,1)
      IF (SIZE(U,2) .NE. ib%tnNo) err = " Inconsistent vector size "//
     2   "to synchronize IB data"

      IF (.NOT.ALLOCATED(ib%cm%nG)) err = " IB comm structure not "//
     2   "initialized. Correction necessary"

      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         DO i=1, ib%msh(iM)%trc%nG
            e = ib%msh(iM)%trc%gE(1,i)
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               incNd(Ac) = 1
            END DO
         END DO
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

      ng = SUM(ib%cm%nG)
      ALLOCATE(gU(m*ng))
      IF (cm%mas()) THEN
         ALLOCATE(rReq(cm%np()))
         rReq = 0
         gU   = 0._RKIND
         ng   = 0
         DO i=1, cm%np()
            nl = ib%cm%nG(i)
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
            IF (i.EQ.1 .OR. ib%cm%nG(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO

         U = 0._RKIND
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
      END SUBROUTINE IB_SYNCGV
!####################################################################
      END MODULE ALLFUN
!####################################################################

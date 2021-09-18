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
!     This routine remeshes any deteriorated mesh during run-time
!     and restarts the simulation.
!
!--------------------------------------------------------------------

      SUBROUTINE REMESHRESTART(timeP)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: timeP(3)

      TYPE(mshType)  :: tMsh

      INTEGER(KIND=IKIND), PARAMETER :: fid=1
      INTEGER(KIND=IKIND) :: iM, i, e, a, Ac, Ec, ierr, iEq, lDof
      REAL(KIND=RKIND) :: t1, t2
      CHARACTER(LEN=stdL) :: sTmp, fTmp

      REAL(KIND=RKIND), ALLOCATABLE, DIMENSION(:,:) :: tempX, gX, gtX,
     2   tempD, gD, gnD, gtD

      INTERFACE
         SUBROUTINE INTERP(lDof, iM, tMsh, sD, tgD)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         TYPE(mshType), INTENT(INOUT) :: tMsh
         INTEGER(KIND=IKIND), INTENT(IN) :: lDof, iM
         REAL(KIND=RKIND), DIMENSION(:,:), INTENT(IN)  :: sD
         REAL(KIND=RKIND), DIMENSION(:,:), INTENT(OUT) :: tgD
         END SUBROUTINE INTERP
      END INTERFACE

      t1 = CPUT()

      sTmp = TRIM(stFileName)//"_last.bin"
      fTmp = TRIM(stFileName)//"_"//STR(rmsh%rTS)//".bin"
      IF (cm%mas()) THEN
         OPEN(fid, FILE=TRIM(sTmp))
         CLOSE(fid, STATUS='DELETE')
      END IF
!     This call is to block all processors
      CALL cm%bcast(rmsh%rTS)
      OPEN(fid, FILE=TRIM(fTmp), ACCESS='DIRECT', RECL=recLn)
      IF (dFlag) THEN
         WRITE(fid, REC=cm%tF()) stamp, rmsh%rTS, time,
     2      CPUT()-timeP(1), eq%iNorm, cplBC%xn,
     3      rmsh%Y0, rmsh%A0, rmsh%D0
      ELSE
         WRITE(fid, REC=cm%tF()) stamp, rmsh%rTS, rmsh%time,
     2      CPUT()-timeP(1), rmsh%iNorm, cplBC%xn, rmsh%Y0, rmsh%A0
      END IF
      CLOSE(fid)
      IF (cm%mas()) CALL SYSTEM("ln -f "//TRIM(fTmp)//" "//TRIM(sTmp))

      gtnNo = 0
      lDof = 3*tDof
      DO iM=1, nMsh
         IF (rmsh%flag(iM)) THEN
            std = " "
            std = "cccccccccccccccccccccccccccccccccccccccccccccccccc"//
     2         "cccccccc"
            std = " "

            IF (rmsh%method .EQ. RMSH_TETGEN) THEN
               std = " Remeshing <"//CLR(TRIM(msh(iM)%name))//
     2            "> using <"//CLR("Tetgen")//"> library at time "//
     3            STR(rmsh%rTS)
            ELSE
               err = "Unexpected behavior in Remesher"
            END IF

            ALLOCATE(tempX(nsd,msh(iM)%nNo))
            ALLOCATE(gD(lDof,msh(iM)%nNo))
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               tempX(:,a) = x(:,Ac) + rmsh%D0(nsd+2:2*nsd+1,Ac)
               gD(1:tDof,a) = rmsh%A0(:,Ac)
               gD(tDof+1:2*tDof,a) = rmsh%Y0(:,Ac)
               gD(2*tDof+1:3*tDof,a) = rmsh%D0(:,Ac)
            END DO

            tMsh%nFa  = 1
            ALLOCATE(tMsh%fa(tMsh%nFa))
            IF (cm%mas()) THEN
               tMsh%gnNo = msh(iM)%gnNo
               tMsh%gnEl = msh(iM)%gnEl
               tMsh%eNoN = msh(iM)%eNoN
            ELSE
               tMsh%gnNo = 0
               tMsh%gnEl = 0
               tMsh%eNoN = 0
            END IF
            ALLOCATE(gX(nsd,tMsh%gnNo))
            gX = GLOBAL(msh(iM), tempX)
            DEALLOCATE(tempX)

            IF (cm%mas()) THEN
               ALLOCATE(tMsh%gIEN(tMsh%eNoN,tMsh%gnEl))
               tMsh%gIEN = msh(iM)%gIEN
               DO e=1, tMsh%gnEl
                  Ec = msh(iM)%otnIEN(e)
                  DO a=1, tMsh%eNoN
                     msh(iM)%gIEN(a,e) = tMsh%gIEN(a,Ec)
                  END DO
               END DO
            END IF

            CALL INTMSHSRF(msh(iM), tMsh%fa(1))

            IF (cm%mas()) THEN
               ALLOCATE(tMsh%fa(1)%x(nsd,tMsh%fa(1)%nNo))
               DO a=1, tMsh%fa(1)%nNo
                  Ac = tMsh%fa(1)%gN(a)
                  tMsh%fa(1)%x(:,a) = gX(:,Ac)
               END DO

               IF (nsd .EQ. 2) THEN
                  err = "Remesher not yet developed for 2D objects"
               ELSE
                  CALL REMESHER_3D(iM, tMsh%fa(1), tMsh)
               END IF

               ALLOCATE(gnD(lDof,tMsh%gnNo))
            ELSE
               ALLOCATE(tMsh%fa(1)%x(0,0), gnD(0,0))
            END IF

            CALL cm%bcast(tMsh%gnNo)
            CALL cm%bcast(tMsh%eNoN)
            CALL cm%bcast(tMsh%gnEl)

            IF (cm%slv()) THEN
               ALLOCATE(tMsh%x(nsd,tMsh%gnNo))
               ALLOCATE(tMsh%gIEN(tMsh%eNoN,tMsh%gnEl))
            END IF

            CALL MPI_BCAST(tMsh%x, nsd*tMsh%gnNo, mpreal, master,
     2         cm%com(), ierr)
            CALL MPI_BCAST(tMsh%gIEN, tMsh%eNoN*tMsh%gnEl, mpint,
     2         master, cm%com(), ierr)
            CALL MPI_BCAST(tMsh%fa(1)%IEN, tMsh%fa(1)%eNoN *
     2         tMsh%fa(1)%nEl, mpint, master, cm%com(), ierr)
            CALL MPI_BCAST(tMsh%fa(1)%gN, tMsh%fa(1)%nNo, mpint,
     2         master, cm%com(), ierr)

            CALL INTERP(lDof, iM, tMsh, gD, gnD)

            IF (cm%mas()) THEN
               msh(iM)%gnNo = tMsh%gnNo
               a = gtnNo + msh(iM)%gnNo
               IF (iM .GT. 1) THEN
                  ALLOCATE(tempX(nsd,gtnNo), tempD(lDof,gtnNo))
                  tempX = gtX
                  tempD = gtD
                  DEALLOCATE(gtX, gtD)
                  ALLOCATE(gtX(nsd,a), gtD(lDof,a))
                  gtX(:,1:gtnNo) = tempX(:,:)
                  gtD(:,1:gtnNo) = tempD(:,:)
                  DEALLOCATE(tempX, tempD)
               ELSE
                  ALLOCATE(gtX(nsd,a))
                  ALLOCATE(gtD(lDof,a))
               END IF
               IF (ALLOCATED(msh(iM)%x)) DEALLOCATE(msh(iM)%x)
               ALLOCATE(msh(iM)%x(nsd,msh(iM)%gnNo))
               gtX(:,gtnNo+1:a) = tMsh%x(:,:) -
     2            gnD(2*tDof+nsd+2:2*tDof+2*nsd+1,:)
               gtD(:,gtnNo+1:a) = gnD(:,:)
               msh(iM)%x(:,:) = gtX(:,gtnNo+1:a)
               gtnNo = a

               CALL SETFACEEBC(tMsh%fa(1), tMsh)

               msh(iM)%eNoN = tMsh%eNoN
               msh(iM)%gnEl = tMsh%gnEl
               IF (ALLOCATED(msh(iM)%gIEN)) DEALLOCATE(msh(iM)%gIEN)
               ALLOCATE(msh(iM)%gIEN(msh(iM)%eNoN,msh(iM)%gnEl))
               msh(iM)%gIEN(:,:) = tMsh%gIEN(:,:)

               CALL DISTMSHSRF(tMsh%fa(1), msh(iM), 1)

               sTmp = TRIM(appPath)//".remesh_tmp.dir"
               fTmp = TRIM(sTmp)//"/"//TRIM(msh(iM)%name)//
     2            "_"//STR(rmsh%rTS)//".vtu"
               msh(iM)%gIEN(:,:) = msh(iM)%gIEN(:,:) - 1
               CALL WRITEVTU(msh(iM), fTmp)
               msh(iM)%gIEN(:,:) = msh(iM)%gIEN(:,:) + 1
            ELSE
               IF (.NOT.ALLOCATED(gtX)) THEN
                  ALLOCATE(gtX(0,0), gtD(0,0))
               END IF
            END IF

            CALL MPI_BARRIER(cm%com(), ierr)

            CALL DESTROY(tMsh)
            DEALLOCATE(gX, gD, gnD)
         ELSE
            ALLOCATE(tempX(nsd,msh(iM)%nNo))
            ALLOCATE(tempD(lDof,msh(iM)%nNo))
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               tempX(1:nsd,a) = x(:,Ac)
               tempD(1:tDof,a) = rmsh%A0(:,Ac)
               tempD(tDof+1:2*tDof,a) = rmsh%Y0(:,Ac)
               tempD(2*tDof+1:3*tDof,a) = rmsh%D0(:,Ac)
            END DO

            IF (cm%mas()) THEN
               a = msh(iM)%gnNo
            ELSE
               a = 0
            END IF
            ALLOCATE(gX(nsd,a), gD(lDof,a))
            gX = GLOBAL(msh(iM), tempX)
            gD = GLOBAL(msh(iM), tempD)
            IF (cm%mas()) THEN
               ALLOCATE(tMsh%gIEN(msh(iM)%eNoN,msh(iM)%gnEl))
               tMsh%gIEN = msh(iM)%gIEN
               DO e=1, msh(iM)%gnEl
                  Ec = msh(iM)%otnIEN(e)
                  msh(iM)%gIEN(:,e) = tMsh%gIEN(:,Ec)
               END DO
               DEALLOCATE(tMsh%gIEN)
            ELSE
            END IF
            DEALLOCATE(tempX, tempD)
            IF (ALLOCATED(msh(iM)%x)) DEALLOCATE(msh(iM)%x)
            ALLOCATE(msh(iM)%x(nsd,msh(iM)%gnNo))

            tMsh%nFa = 1
            ALLOCATE(tMsh%fa(tMsh%nFa))
            CALL INTMSHSRF(msh(iM), tMsh%fa(1))

            IF (cm%mas()) THEN
               ALLOCATE(tMsh%fa(1)%x(nsd,tMsh%fa(1)%nNo))
               DO i=1, tMsh%fa(1)%nNo
                  Ac = tMsh%fa(1)%gN(i)
                  tMsh%fa(1)%x(:,i) = gX(:,Ac)
               END DO
               msh(iM)%x = gX

               CALL DISTMSHSRF(tMsh%fa(1), msh(iM), 2)

               a = gtnNo + a
               IF (iM .GT. 1) THEN
                  ALLOCATE(tempX(nsd,gtnNo))
                  ALLOCATE(tempD(lDof,gtnNo))
                  tempX = gtX
                  tempD = gtD
                  IF (ALLOCATED(gtX)) DEALLOCATE(gtX)
                  IF (ALLOCATED(gtD)) DEALLOCATE(gtD)
                  ALLOCATE(gtX(nsd,a), gtD(lDof,a))
                  gtX(:,1:gtnNo) = tempX(:,:)
                  gtD(:,1:gtnNo) = tempD(:,:)
                  DEALLOCATE(tempX, tempD)
               ELSE
                  ALLOCATE(gtX(nsd,a), gtD(lDof,a))
               END IF
               gtX(:,gtnNo+1:a) = gX(:,:)
               gtD(:,gtnNo+1:a) = gD(:,:)
               gtnNo = a
            ELSE
               ALLOCATE(tMsh%fa(1)%x(0,0))
               IF (ALLOCATED(gtX)) DEALLOCATE(gtX)
               IF (ALLOCATED(gtD)) DEALLOCATE(gtD)
               ALLOCATE(gtX(0,0), gtD(0,0))
            END IF
            DEALLOCATE(gX, gD)
            CALL DESTROY(tMsh)
         END IF ! reMesh flag
      END DO
      CALL cm%bcast(gtnNo)
      DEALLOCATE(x, rmsh%A0, rmsh%Y0, rmsh%D0)

      IF (cm%mas()) THEN
         ALLOCATE(x(nsd,gtnNo))
         ALLOCATE(rmsh%A0(tDof,gtnNo))
         ALLOCATE(rmsh%Y0(tDof,gtnNo))
         ALLOCATE(rmsh%D0(tDof,gtnNo))
         x = gtX
         DO a=1, gtnNo
            rmsh%A0(:,a) = gtD(1:tDof,a)
            rmsh%Y0(:,a) = gtD(tDof+1:2*tDof,a)
            rmsh%D0(:,a) = gtD(2*tDof+1:3*tDof,a)
         END DO
      ELSE
         ALLOCATE(rmsh%A0(0,0),rmsh%Y0(0,0),rmsh%D0(0,0))
      END IF
      DEALLOCATE(gtX, gtD)

      DO iM=1, nMsh
         IF (cm%mas()) THEN
            IF (ALLOCATED(msh(iM)%eDist))  DEALLOCATE(msh(iM)%eDist)
            IF (ALLOCATED(msh(iM)%eId))    DEALLOCATE(msh(iM)%eId)
            IF (ALLOCATED(msh(iM)%gN))     DEALLOCATE(msh(iM)%gN)
            IF (ALLOCATED(msh(iM)%gpN))    DEALLOCATE(msh(iM)%gpN)
            IF (ALLOCATED(msh(iM)%IEN))    DEALLOCATE(msh(iM)%IEN)
            IF (ALLOCATED(msh(iM)%otnIEN)) DEALLOCATE(msh(iM)%otnIEN)
            IF (ALLOCATED(msh(iM)%INN))    DEALLOCATE(msh(iM)%INN)
            IF (ALLOCATED(msh(iM)%lN))     DEALLOCATE(msh(iM)%lN)
            IF (ALLOCATED(msh(iM)%eIEN))   DEALLOCATE(msh(iM)%eIEN)
            IF (ALLOCATED(msh(iM)%sbc))    DEALLOCATE(msh(iM)%sbc)
            IF (ALLOCATED(msh(iM)%iGC))    DEALLOCATE(msh(iM)%iGC)
            IF (ALLOCATED(msh(iM)%nW))     DEALLOCATE(msh(iM)%nW)
            IF (ALLOCATED(msh(iM)%w))      DEALLOCATE(msh(iM)%w)
            IF (ALLOCATED(msh(iM)%xi))     DEALLOCATE(msh(iM)%xi)
            IF (ALLOCATED(msh(iM)%xib))    DEALLOCATE(msh(iM)%xib)
            IF (ALLOCATED(msh(iM)%x))      DEALLOCATE(msh(iM)%x)
            IF (ALLOCATED(msh(iM)%N))      DEALLOCATE(msh(iM)%N)
            IF (ALLOCATED(msh(iM)%Nb))     DEALLOCATE(msh(iM)%Nb)
            IF (ALLOCATED(msh(iM)%nV))     DEALLOCATE(msh(iM)%nV)
            IF (ALLOCATED(msh(iM)%fN))     DEALLOCATE(msh(iM)%fN)
            IF (ALLOCATED(msh(iM)%Nx))     DEALLOCATE(msh(iM)%Nx)
            IF (ALLOCATED(msh(iM)%Nxx))    DEALLOCATE(msh(iM)%Nxx)
            CALL DESTROY(msh(iM)%nAdj)
            CALL DESTROY(msh(iM)%eAdj)
            DO i=1, msh(iM)%nFs
               CALL DESTROY(msh(iM)%fs(i))
            END DO
            IF (ALLOCATED(msh(iM)%fs))     DEALLOCATE(msh(iM)%fs)
         ELSE
            CALL DESTROY(msh(iM))
         END IF
      END DO
      IF (cm%slv()) DEALLOCATE(msh)

      IF (ALLOCATED(eq)) THEN
         DO iEq=1, nEq
            CALL DESTROY(eq(iEq))
         END DO
         DEALLOCATE(eq)
      END IF

      IF (ALLOCATED(colPtr))   DEALLOCATE(colPtr)
      IF (ALLOCATED(dmnID))    DEALLOCATE(dmnID)
      IF (ALLOCATED(ltg))      DEALLOCATE(ltg)
      IF (ALLOCATED(rowPtr))   DEALLOCATE(rowPtr)
      IF (ALLOCATED(idMap))    DEALLOCATE(idMap)
      IF (ALLOCATED(cmmBdry))  DEALLOCATE(cmmBdry)
      IF (ALLOCATED(iblank))   DEALLOCATE(iblank)
      IF (ALLOCATED(Ao))       DEALLOCATE(Ao)
      IF (ALLOCATED(An))       DEALLOCATE(An)
      IF (ALLOCATED(Do))       DEALLOCATE(Do)
      IF (ALLOCATED(Dn))       DEALLOCATE(Dn)
      IF (ALLOCATED(R))        DEALLOCATE(R)
      IF (ALLOCATED(Val))      DEALLOCATE(Val)
      IF (ALLOCATED(Yo))       DEALLOCATE(Yo)
      IF (ALLOCATED(Yn))       DEALLOCATE(Yn)
      IF (ALLOCATED(Bf))       DEALLOCATE(Bf)
      cplBC%nFa = 0

!     Additional physics based variables to be deallocated
      IF (ALLOCATED(Ad))       DEALLOCATE(Ad)
      IF (ALLOCATED(Rd))       DEALLOCATE(Rd)
      IF (ALLOCATED(Kd))       DEALLOCATE(Kd)
      IF (ALLOCATED(pS0))      DEALLOCATE(pS0)
      IF (ALLOCATED(pSn))      DEALLOCATE(pSn)
      IF (ALLOCATED(pSa))      DEALLOCATE(pSa)

!     Varwall properties -----------------------------------------------
      IF (ALLOCATED(vWP0))      DEALLOCATE(vWP0)
!     ------------------------------------------------------------------

      t2 = CPUT()

      std = " Time taken for remeshing: "//STR(t2-t1)//" (s)"
      std = " "
      std = "cccccccccccccccccccccccccccccccccccccccccccccccccc"//
     2   "cccccccc"
      std = " "

      RETURN
      END SUBROUTINE REMESHRESTART
!--------------------------------------------------------------------
      SUBROUTINE INTMSHSRF(lM, lFa)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      INTEGER(KIND=IKIND) :: a, e, Ac, iFa, eNoNb, eoff

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:)

      eNoNb = lM%fa(1)%eNoN
      lFa%eNoN = eNoNb
      lFa%gnEl = SUM(lM%fa(:)%gnEl)
      lFa%nNo  = 0
      lFa%name = "old_mesh_surface"
      lFa%nEl = lFa%gnEl

      ALLOCATE(lFa%IEN(lFa%eNoN,lFa%gnEl), lFa%gE(lFa%gnEl))

      IF(eNoNb .NE. lM%fa(lM%nFa)%eNoN) THEN
         wrn = "    Remeshing not formulated for different face types"
         RETURN
      END IF

      IF(cm%mas()) THEN
         DO iFa=1, lM%nFa
            eoff = 0
            IF(iFa .GT. 1) eoff = SUM(lM%fa(1:iFa-1)%gnEl)
            DO e=eoff+1, eoff+lM%fa(iFa)%gnEl
               lFa%gE(e) = lM%fa(iFa)%gebc(1,e-eoff)
               lFa%IEN(:,e) = lM%fa(iFa)%gebc(2:1+eNoNb,e-eoff)
            END DO
         END DO

         CALL CALCNBC(lM, lFa)

         ALLOCATE(incNd(lM%gnNo))
         incNd = 0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            incNd(Ac) = a
         END DO
         DO e=1, lFa%gnEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               lFa%IEN(a,e) = incNd(Ac)
            END DO
         END DO
         DEALLOCATE(incNd)
      END IF

      CALL cm%bcast(lFa%nNo)
      IF (cm%slv()) ALLOCATE(lFa%gN(lFa%nNo))

      END SUBROUTINE INTMSHSRF
!--------------------------------------------------------------------
      SUBROUTINE REMESHER_3D(iM, lFa, lM)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iM
      TYPE(faceType), INTENT(INOUT) :: lFa
      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: e, i, Ac, fid, iOK
      REAL(KIND=RKIND) :: rparams(3)
      TYPE(fileType) :: fTmp

      rparams(1) = rmsh%maxRadRatio
      rparams(2) = rmsh%minDihedAng
      rparams(3) = rmsh%maxEdgeSize(iM)
      iOK = 0
      IF (rmsh%method .EQ. RMSH_TETGEN) THEN
         CALL remesh3d_tetgen(lFa%nNo, lFa%nEl, lFa%x, lFa%IEN,
     2      rparams, iOK)
         IF (iOK .LT. 0)
     2      err = "Fatal! TetGen returned with error. Check log"
      ELSE
         err = "Unknown remesher choice."
      END IF

      dbg = " Reading new mesh data.."
      fTmp%fname  = "new-vol-mesh.ele"
      fid = fTmp%open()
      lM%gnEl = 0
      DO
         READ (fid,*,END=111)
         lM%gnEl = lM%gnEl + 1
      END DO
 111  REWIND(fid)
      IF (rmsh%method .EQ. RMSH_TETGEN) lM%gnEl = lM%gnEl - 1
      std = "    Number of elements after remesh "//STR(lM%gnEl)

      IF (ALLOCATED(lM%gIEN)) DEALLOCATE(lM%gIEN)
      ALLOCATE(lM%gIEN(lM%eNoN,lM%gnEl))
      IF (rmsh%method .EQ. RMSH_TETGEN) READ(fid,*)
      DO e=1, lM%gnEl
         READ(fid,*) i, lM%gIEN(:,e)
      END DO
      CLOSE(fid, STATUS='DELETE')

      fTmp%fname  = "new-vol-mesh.node"
      fid = fTmp%open()
      lM%gnNo = 0
      DO
         READ (fid,*,END=112)
         lM%gnNo = lM%gnNo + 1
      END DO
 112  REWIND(fid)
      IF (rmsh%method .EQ. RMSH_TETGEN) lM%gnNo = lM%gnNo - 1
      std = "    Number of vertices after remesh "//STR(lM%gnNo)

      IF (ALLOCATED(lM%x)) DEALLOCATE(lM%x)
      ALLOCATE(lM%x(nsd,lM%gnNo))
      IF (rmsh%method .EQ. RMSH_TETGEN) READ(fid,*)
      DO Ac=1, lM%gnNo
         READ(fid,*) i, lM%x(:,Ac)
      END DO
      CLOSE(fid, STATUS='DELETE')

      CALL SELECTELE(lM)

      RETURN
      END SUBROUTINE REMESHER_3D
!--------------------------------------------------------------------
      SUBROUTINE SETFACEEBC(lFa, lM)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      TYPE(mshType),  INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: i, j, e, a, b
      INTEGER(KIND=IKIND) :: Ac, maxAssocEl
      REAL(KIND=RKIND) :: sn, Jac

      INTEGER(KIND=IKIND), ALLOCATABLE :: nAssocEl(:), bin(:,:),
     2   assocEl(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), v(:,:)

      dbg = " Finding global surface elems"
      ALLOCATE(nAssocEl(lM%gnNo))
      nAssocEl(:) = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            nAssocEl(Ac) = nAssocEl(Ac) + 1
         END DO
      END DO

      maxAssocEl = MAXVAL(nAssocEl)
      ALLOCATE(assocEl(maxAssocEl,lM%gnNo))
      nAssocEl(:) = 0
      assocEl(:,:) = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            nAssocEl(Ac) = nAssocEl(Ac) + 1
            assocEl(nAssocEl(Ac),Ac) = e
         END DO
      END DO

      IF (ALLOCATED(lFa%gE)) DEALLOCATE(lFa%gE)
      ALLOCATE(lFa%gE(lFa%nEl))
      ALLOCATE(bin(maxAssocEl,2))
      DO e=1, lFa%nEl
         bin(:,:) = 0
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            IF (a .EQ. 1) THEN
               bin(1:nAssocEl(Ac),1) = assocEl(1:nAssocEl(Ac),Ac)
               bin(1:nAssocEl(Ac),2) = 1
            ELSE
               DO i=1, nAssocEl(Ac)
                  DO j=1, maxAssocEl
                     IF (bin(j,1) .EQ. 0) THEN
                        bin(j,1) = assocEl(i,Ac)
                        bin(j,2) = 1
                        EXIT
                     ELSE IF (bin(j,1) .EQ. assocEl(i,Ac)) THEN
                        bin(j,2) = bin(j,2) + 1
                     END IF
                  END DO ! j
               END DO ! i
            END IF
         END DO ! a

         DO j=1, maxAssocEl
            IF (bin(j,2) .EQ. lFa%eNoN) THEN
               lFa%gE(e) = bin(j,1)
               EXIT
            END IF
         END DO ! j
      END DO ! e
      DEALLOCATE(nAssocEl,assocEl,bin)

      ALLOCATE(xl(nsd,lM%eNoN), v(nsd,lM%eNoN))
c     Check mesh quality and reset IEN if necessary
      std = " Checking IEN array"
      e = lFa%gE(1)
      xl = lM%x(:,lM%gIEN(:,e))
      v(:,1) = xl(:,2) - xl(:,1)
      v(:,2) = xl(:,3) - xl(:,2)
      v(:,3) = xl(:,4) - xl(:,3)
      v(:,4) = CROSS(v(:,1:2))
      sn = SGN(SUM(v(:,3)*v(:,4)))
      IF (sn .EQ. 1) THEN
         a=1; b=2
         std = "    Reordering element connectivity"
         DO e=1, lM%gnEl
            Ac = lM%gIEN(a,e)
            lM%gIEN(a,e) = lM%gIEN(b,e)
            lM%gIEN(b,e) = Ac
         END DO
      ELSE IF (sn .EQ. 0) THEN
         err = "Surface element "//STR(e)//" is distorted"
      END IF

      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            xl(:,a) = lM%x(:,Ac)
         END DO
         Jac = JACOBIAN(nsd, lM%eNoN, xl, lM%nX(:,:,1))
         IF (Jac .LT. 0) err = "Remeshing didn't improve mesh quality"
      END DO
      DEALLOCATE(xl,v)

      RETURN
      END SUBROUTINE SETFACEEBC
!--------------------------------------------------------------------
      SUBROUTINE DISTMSHSRF(lFa, lM, iOpt)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      TYPE(mshType),  INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: iOpt

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: e, a, Ac, iFa, eoff
      CHARACTER(LEN=stdL) :: sTmp, fTmp

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:)

      dbg = " Distributing surface mesh"
      sTmp = TRIM(appPath)//".remesh_tmp.dir"
      INQUIRE(FILE=TRIM(sTmp)//"/.", EXIST=flag)
      IF (.NOT. flag) THEN
         CALL SYSTEM("mkdir  -p  "//TRIM(sTmp))
      END IF

      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%gN(Ac)
            lFa%IEN(a,e) = Ac
         END DO
      END DO

      DO iFa=1, lM%nFa
         lM%fa(iFa)%nEl = lM%fa(iFa)%gnEl
         IF (ALLOCATED(lM%fa(iFa)%IEN)) DEALLOCATE(lM%fa(iFa)%IEN)
         IF (ALLOCATED(lM%fa(iFa)%gE)) DEALLOCATE(lM%fa(iFa)%gE)
         ALLOCATE(lM%fa(iFa)%IEN(lM%fa(iFa)%eNoN,lM%fa(iFa)%nEl))
         ALLOCATE(lM%fa(iFa)%gE(lM%fa(iFa)%nEl))

         eoff = 0
         IF (iFa .GT. 1) eoff = SUM(lM%fa(1:iFa-1)%gnEl)
         DO e=1, lM%fa(iFa)%gnEl
            lM%fa(iFa)%gE(e) = lFa%gE(eoff+e)
            lM%fa(iFa)%IEN(:,e) = lFa%IEN(:,eoff+e)
         END DO

         IF (ALLOCATED(lM%fa(iFa)%gN)) DEALLOCATE(lM%fa(iFa)%gN)
         CALL CALCNBC(lM, lM%fa(iFa))

         IF (ALLOCATED(lM%fa(iFa)%x)) DEALLOCATE(lM%fa(iFa)%x)
         ALLOCATE(lM%fa(iFa)%x(nsd, lM%fa(iFa)%nNo))
         DO a=1, lM%fa(iFa)%nNo
            Ac = lM%fa(iFa)%gN(a)
            lM%fa(iFa)%x(:,a) = lM%x(:,Ac)
         END DO

         lM%fa(iFa)%gebc(1,:) = lM%fa(iFa)%gE(:)
         lM%fa(iFa)%gebc(2:1+lM%fa(iFa)%eNoN,:) = lM%fa(iFa)%IEN(:,:)

         IF (iOpt.EQ.1) THEN
            fTmp = TRIM(sTmp)//"/"//TRIM(lM%fa(iFa)%name)//"_"//
     2         STR(rmsh%rTS)//".vtp"
            ALLOCATE(incNd(lM%gnNo))
            incNd = 0
            DO a=1, lM%fa(iFa)%nNo
               Ac = lM%fa(iFa)%gN(a)
               incNd(Ac) = a
            END DO
            DO e=1, lM%fa(iFa)%nEl
               DO a=1, lM%fa(iFa)%eNoN
                  Ac = lM%fa(iFa)%IEN(a,e)
                  lM%fa(iFa)%IEN(a,e) = incNd(Ac) - 1
               END DO
            END DO
            CALL WRITEVTP(lM%fa(iFa), fTmp)
            DO e=1, lM%fa(iFa)%nEl
               DO a=1, lM%fa(iFa)%eNoN
                  Ac = lM%fa(iFa)%IEN(a,e) + 1
                  Ac = lM%fa(iFa)%gN(Ac)
                  lM%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
            DEALLOCATE(incNd)
         END IF
      END DO

      RETURN
      END SUBROUTINE DISTMSHSRF
!--------------------------------------------------------------------
!     Interpolation of data variables from source mesh to target mesh
      SUBROUTINE INTERP(lDof, iM, tMsh, sD, tgD)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: tMsh
      INTEGER(KIND=IKIND), INTENT(IN) :: lDof, iM
      REAL(KIND=RKIND), DIMENSION(:,:), INTENT(IN)  :: sD
      REAL(KIND=RKIND), DIMENSION(:,:), INTENT(OUT) :: tgD(:,:)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: a, e, b, i, iFa, ierr, Ac, Bc, Ec, nn, ne,
     2   try, nbnd, nNo, nEl, eNon, gnNo, gnEl, maxKNE, maxKNN, probe,
     3   bTag
      REAL(KIND=RKIND) :: dS
      TYPE(queueType) :: rootNdQ

      LOGICAL, ALLOCATABLE :: chckNp(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:), gE(:), tagNd(:),
     2   masEList(:), srcAdjEl(:,:), tgtAdjNd(:,:), eList(:), rootEl(:),
     3   srfNds(:), sCount(:), disp(:), tmpL(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Xp(:), Nsf(:), gNsf(:,:),
     2   Dg(:,:), tmpX(:,:), vec(:), gvec(:)

      INTERFACE
         SUBROUTINE DISTRN(iM, lM, Dg, nNo, gN)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         TYPE(mshType), INTENT(INOUT) :: lM
         INTEGER(KIND=IKIND), INTENT(IN) :: iM
         INTEGER(KIND=IKIND), INTENT(OUT) :: nNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
         REAL(KIND=RKIND), INTENT(IN) :: Dg(nsd,tnNo)
         END SUBROUTINE DISTRN
      END INTERFACE

      INTERFACE
         SUBROUTINE DISTRE(lM, nEl, gE)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         TYPE(mshType), INTENT(INOUT) :: lM
         INTEGER(KIND=IKIND), INTENT(OUT) :: nEl
         INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:)
         END SUBROUTINE DISTRE
      END INTERFACE

      INTERFACE
         SUBROUTINE GETADJESRC(lM, kneList)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         TYPE(mshType), INTENT(IN) :: lM
         INTEGER(KIND=IKIND), ALLOCATABLE :: kneList(:,:)
         END SUBROUTINE GETADJESRC
      END INTERFACE

      INTERFACE
         SUBROUTINE GETADJNTGT(lM, nNo, nEl, gN, gE, knnList)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         TYPE(mshType), INTENT(IN) :: lM
         INTEGER(KIND=IKIND), INTENT(IN) :: nNo, nEl, gN(nNo), gE(nEl)
         INTEGER(KIND=IKIND), ALLOCATABLE :: knnList(:,:)
         END SUBROUTINE GETADJNTGT
      END INTERFACE

      INTERFACE
         SUBROUTINE FINDN(xp, iM, Dg, eList, Ec, Nsf)
         USE COMMOD
         USE UTILMOD
         USE ALLFUN
         IMPLICIT NONE
         INTEGER(KIND=IKIND), INTENT(IN) :: iM, eList(:)
         INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
         REAL(KIND=RKIND), INTENT(IN) :: xp(nsd+1), Dg(nsd,tnNo)
         REAL(KIND=RKIND), INTENT(OUT) :: Nsf(msh(iM)%eNoN)
         END SUBROUTINE FINDN
      END INTERFACE

      std = " Interpolating.."
      IF (cm%seq()) wrn = "Interpolation not optimized for serial job"

      IF (nsd+1 .NE. msh(iM)%eNoN) err = "Inconsistent element type "//
     2   "for interpolation. Can support 2D Tri or 3D Tet elements only"

      gnNo = tMsh%gnNo
      eNoN = tMsh%eNoN
      gnEl = tMsh%gnEl

!     Distribute the new mesh nodes among all the processors
!     Need to transfer mesh displacement
      ALLOCATE(Dg(nsd,tnNo))
      i = nsd+1
      Dg(:,:) = rmsh%D0(i+1:i+nsd,:)
      CALL DISTRN(iM, tMsh, Dg, nNo, gN)

!     Distribute elements of the new mesh to all processors
      CALL DISTRE(tMsh, nEl, gE)

!     Setup data structures for octree search
!     Get adjacent cells for source (old) mesh
      CALL GETADJESRC(msh(iM), srcAdjEl)
      maxKNE = SIZE(srcAdjEl,1)

!     Get adjacent nodes for each node on the new mesh
      CALL GETADJNTGT(tMsh, nNo, nEl, gN, gE, tgtAdjNd)
      maxKNN = SIZE(tgtAdjNd,1)
      DEALLOCATE(gE)

      ALLOCATE(Xp(nsd+1), Nsf(eNoN), gNsf(eNoN,nNo), tagNd(gnNo),
     2 gE(nNo), rootEl(nNo), chckNp(nNo), masEList(msh(iM)%nEl))

      gNsf   = 0._RKIND
      gE     = 0
      tagNd  = 0
      rootEl = 0
      chckNp = .FALSE.
      DO e=1, msh(iM)%nEl
         masEList(e) = e
      END DO

!     Determine boundary nodes on the new mesh, where interpolation is
!     not needed, or boundary search is performed
      ALLOCATE(tmpL(gnNo))
      tmpL = 0
      DO e=1, tMsh%fa(1)%nEl
         DO a=1, tMsh%fa(1)%eNoN
            Ac = tMsh%fa(1)%IEN(a,e)
            IF (tmpL(Ac) .EQ. 0) tmpL(Ac) = 1
         END DO
      END DO

      ALLOCATE(srfNds(nNo))
      srfNds = 0
      DO a=1, nNo
         Ac = gN(a)
        IF (tmpL(Ac) .GT. 0) srfNds(a) = 1
      END DO
      DEALLOCATE(tmpL)

      bTag = 2*cm%np()
      DO a=1, nNo
         Ac = gN(a)
         IF (srfNds(a) .GT. 0) THEN
            chckNp(a) = .TRUE.
            tagNd(Ac) = bTag
            gE(a)     = 0
            gNsf(:,a) = 0._RKIND
         END IF
      END DO

!     Find a starting element through brute search and the starting node
!     is used to initialize the queue (local numbering)
      DO a=1, nNo
         IF (chckNp(a)) CYCLE
         Ac = gN(a)
         Xp = 1._RKIND
         Xp(1:nsd) = tMsh%x(:,Ac)
         CALL FINDN(Xp, iM, Dg, masEList, Ec, Nsf)
         chckNp(a) = .TRUE.
!     Once an element in the source mesh is found, set 'this' element as
!     'root' element for all the neighbor nodes of the target mesh. Then
!     add these neighbor nodes to the queue to form the 'front'
         IF (Ec .GT. 0) THEN
            gE(a) = Ec
            tagNd(Ac) = cm%tF()
            gNsf(:,a) = Nsf(:)
            DO nn=1, maxKNN
               b = tgtAdjNd(nn,a)
               IF (b .GT. 0) THEN
                  CALL ENQUEUE(rootNdQ, b)
                  rootEl(b) = Ec
               END IF
            END DO
            IF (rootNdQ%n .GT. 1) EXIT
         END IF
      END DO

!     Node-Cell search begins here. Uses fastest grid-to-grid algorithm
!     based on advancing front in the nearest vicinity
      ALLOCATE(eList(maxKNE)); eList = 0
      DO WHILE (DEQUEUE(rootNdQ, probe))
         IF (ALL(chckNp(:))) EXIT
         IF (chckNp(probe)) CYCLE
         Ac = gN(probe)
         Xp = 1._RKIND
         Xp(1:nsd) = tMsh%x(:,Ac)
         eList = 0
         DO e=1, maxKNE
            Ec = srcAdjEl(e,rootEl(probe))
            IF (Ec .GT. 0) eList(e) = Ec
         END DO

         DO try=1, 2
            CALL FINDN(Xp, iM, Dg, eList, Ec, Nsf)
!      IF an element is found, the neighbors of the node are initialized
!      with this element as root element and are added to the queue
!      (front propagation)
            IF (Ec .GT. 0) THEN
              gE(probe) = Ec
               tagNd(Ac) = cm%tF()
               gNsf(:,probe) = Nsf(:)
               DO nn=1, maxKNN
                  a = tgtAdjNd(nn,probe)
                  IF (a .GT. 0) THEN
                     CALL ENQUEUE(rootNdQ, a)
                     rootEl(a) = Ec
                  END IF
               END DO
               DEALLOCATE(eList)
               ALLOCATE(eList(maxKNE)); eList = 0
               EXIT
            END IF
!      If failed in first attempt, try again with increased search area
!      but add neighbors uniquely without repetition
            IF (try .EQ. 1) THEN
               ALLOCATE(tmpL(nEl))
               tmpL = 0
               ne = 0
               DO e=1, maxKNE
                  Ec = eList(e)
                  IF (Ec .GT. 0) THEN
                     DO a=1, maxKNE
                        Bc = srcAdjEl(a,Ec)
                        IF (Bc.NE.0 .AND. Bc.NE.Ec) THEN
                           flag = .TRUE.
                           DO i=1, ne
                              IF (Bc .EQ. tmpL(i)) THEN
                                 flag = .FALSE.
                                 EXIT
                              END IF
                           END DO
                           IF (flag) THEN
                              ne = ne + 1
                              tmpL(ne) = Bc
                           END IF
                        END IF
                     END DO
                  END IF
               END DO
               DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               eList = tmpL(1:ne)
               DEALLOCATE(tmpL)
            ELSE IF (try .EQ. 2) THEN
               DEALLOCATE(eList)
               ALLOCATE(eList(maxKNE)); eList = 0
            END IF
         END DO ! try
         chckNp(probe) = .TRUE.
      END DO
      CALL DESTROY(rootNdQ)

      ALLOCATE(tmpL(gnNo))
      tmpL = 0
      CALL MPI_ALLREDUCE(tagNd, tmpL, gnNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      nn   = 0
      nbnd = 0
      DO a=1, gnNo
         IF (tmpL(a) .GT. 0) nn = nn + 1
         IF (tmpL(a) .EQ. 2*cm%np()) nbnd = nbnd + 1
      END DO

!     Assign tag for nodes that were interpolated in other processors
      tagNd = 0
      DO a=1, nNo
         Ac = gN(a)
         tagNd(Ac) = tmpL(Ac)
      END DO

!     Now use brute force to find left over non-interpolated nodes
      DO a=1, nNo
         Ac = gN(a)
         IF (tagNd(Ac) .EQ. 0) THEN
            Xp = 1._RKIND
            Xp(1:nsd) = tMsh%x(:,Ac)
            CALL FINDN(Xp, iM, Dg, masElist, Ec, Nsf)
            IF (Ec .GT. 0) THEN
               gE(a) = Ec
               tagNd(Ac) = cm%tF()
               gNsf(:,a) = Nsf
            END IF
         END IF
      END DO

!     Nodes belonging to other procs are reassigned 0
      DO a=1, nNo
         Ac = gN(a)
         IF (tagNd(Ac).NE.cm%tF() .AND. tagNd(Ac).NE.bTag) THEN
            gE(a) = 0
            tagNd(Ac) = 0
            gNsf(:,a) = 0._RKIND
         END IF
      END DO
      DEALLOCATE(tmpL)

!     Now that all the elements have been found, data is interpolated
!     from the source to the target mesh
      ALLOCATE(tmpX(lDof,nNo))
      tmpX = 0._RKIND
      DO a=1, nNo
         Ac = gN(a)
         IF (tagNd(Ac) .EQ. cm%tF()) THEN
            Ec = gE(a)
            Nsf = gNsf(:,a)
            DO i=1, eNoN
               Bc = msh(iM)%IEN(i,Ec)
               Bc = msh(iM)%lN(Bc)
               tmpX(:,a) = tmpX(:,a) + Nsf(i)*sD(:,Bc)
            END DO
         END IF
      END DO

!     Since there is no direct mapping for face data, we use L2 norm
!     to find the nearest face node and copy its solution. This requires
!     face node/IEN structure to NOT be changed during remeshing.
      ALLOCATE(tmpL(nNo))
      tmpL = 0
      DO a=1, nNo
         Ac = gN(a)
         IF (srfNds(a) .GT. 0) THEN
            flag = .FALSE.
            DO iFa=1, msh(iM)%nFa
               DO b=1, msh(iM)%fa(iFa)%nNo
                  Bc = msh(iM)%fa(iFa)%gN(b)
                  dS = SQRT(SUM( (x(:,Bc)+Dg(:,Bc) -
     2               tMsh%x(:,Ac))**2._RKIND ))
                  IF (dS .LT. 1.E-12_RKIND) THEN
                     tmpL(a) = Bc
                     flag = .TRUE.
                     EXIT
                  END IF
               END DO
               IF (flag) EXIT
            END DO
            IF (flag) THEN
               Bc = msh(iM)%lN(tmpL(a))
               tmpX(:,a) = sD(:,Bc)
            ELSE
               tagNd(Ac) = 0
            END IF
         END IF
      END DO
      DEALLOCATE(tmpL)

!     Map the tagged nodes and solution to local vector within a proc,
!     including boundary nodes. Since the boundary nodes can be overlapping
!     across different procs, these are repeated. But this will not cause
!     problem as the solution is simply overwritten depending on the face pointer.
      nn = 0
      DO a=1, nNo
         Ac = gN(a)
         IF (tagNd(Ac) .GT. 0) nn = nn + 1
      END DO
      ALLOCATE(tmpL(nn))
      nn = 0
      DO a=1, nNo
         Ac = gN(a)
         IF (tagNd(Ac) .GT. 0) THEN
            nn = nn + 1
            tmpL(nn) = a
         END IF
      END DO

      DEALLOCATE(gNsf)
      ALLOCATE(gNsf(lDof,nn))
      DO i=1, nn
         a = tmpL(i)
         Ac = gN(a)
         gNsf(:,i) = tmpX(:,a)
         tmpL(i) = Ac
      END DO
      DEALLOCATE(tmpX)

      IF (cm%mas()) THEN
         ALLOCATE(disp(cm%np()))
         disp = 0
      ELSE
         ALLOCATE(disp(0))
      END IF

      CALL MPI_GATHER(nn, 1, mpint, disp, 1, mpint, master,
     2   cm%com(), ierr)

      IF (cm%mas()) THEN
         i = SUM(disp(:))
         i = i*(1 + lDof)
         ALLOCATE(gvec(i))
         ALLOCATE(sCount(cm%np()))
         DO i=1, cm%np()
            sCount(i) = disp(i)*(1+lDof)
         END DO
         disp(1) = 0
         DO i=2, cm%np()
            disp(i) = disp(i-1) + sCount(i-1)
         END DO
      ELSE
         ALLOCATE(gvec(0), sCount(0))
      END IF
      ALLOCATE(vec((lDof+1)*nn))
      e = 0
      DO a=1, nn
         e = e + 1
         vec(e) = REAL(tmpL(a), KIND=RKIND)
         DO b=1, lDof
            e = e + 1
            vec(e) = gNsf(b,a)
         END DO
      END DO
      DEALLOCATE(gNsf)

      CALL MPI_GATHERV(vec, (1+lDof)*nn, mpreal, gvec, sCount, disp,
     2   mpreal, master, cm%com(), ierr)

      IF (cm%mas()) THEN
         tgD = 0._RKIND
         nn = SUM(sCount(:))
         i = 0
         DO
            i = i + 1
            Ac = NINT(gvec(i), KIND=IKIND)
            DO b=1, lDof
               i = i + 1
               tgD(b,Ac) = gvec(i)
            END DO
            IF (i .EQ. nn) EXIT
         END DO
      END IF
      DEALLOCATE(vec, gvec, sCount, disp)

      RETURN
      END SUBROUTINE INTERP
!--------------------------------------------------------------------
!     Distribute the new mesh nodes among all the processors
      SUBROUTINE DISTRN(iM, lM, Dg, nNo, gN)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: iM
      INTEGER(KIND=IKIND), INTENT(OUT) :: nNo
      INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
      REAL(KIND=RKIND), INTENT(IN) :: Dg(nsd,tnNo)

      INTEGER(KIND=IKIND) :: i, a, b, Ac, gnNo, ierr
      REAL(KIND=RKIND) :: f, tol, dS, minS

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), tmpI(:)

      gnNo = lM%gnNo
      i = 0
      f = 2.5E-2_RKIND
      ALLOCATE(part(gnNo), tmpI(gnNo))
 003  part  = 0
      tmpI  = 0
      nNo   = 0
      f     = 2._RKIND*f
      tol   = (1._RKIND+f) * rmsh%maxEdgeSize(iM)
      i = i+1
      DO a=1, gnNo
         IF (part(a) .NE. 0) CYCLE
         minS  = HUGE(minS)
         DO b=1, msh(iM)%nNo
            Ac = msh(iM)%gN(b)
            dS = SQRT(SUM((x(:,Ac)+Dg(:,Ac)-lM%x(:,a))**2._RKIND))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            nNo = nNo + 1
            part(a) = cm%tF()
         END IF
      END DO

      CALL MPI_ALLREDUCE(part, tmpI, gnNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      b = 0
      DO a=1, gnNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. gnNo) THEN
         wrn = "Found only "//STR(b)//" nodes in pass "//STR(i)//
     2      " out of "//STR(gnNo)//" nodes"
         IF (i .GT. 5) err = "Could not distribute all nodes in "//
     2      STR(i)//" passes. Try changing tolerance."
         GOTO 003
      END IF

      ALLOCATE(gN(nNo), lM%gpN(gnNo))
      nNo = 0
      DO a=1, gnNo
         IF (part(a) .EQ. cm%tF()) THEN
            nNo = nNo + 1
            gN(nNo) = a
         END IF
      END DO
!      WRITE(*,'(A)') "CPU "//STR(cm%tF())//" nNo = "//STR(nNo)
      lM%gpN = part
      DEALLOCATE(part, tmpI)

      RETURN
      END SUBROUTINE DISTRN
!--------------------------------------------------------------------
!     Distribute the new mesh elements to all processors
      SUBROUTINE DISTRE(lM, nEl, gE)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(OUT) :: nEl
      INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:)

      INTEGER(KIND=IKIND) :: e, a, i, Ac, gnEl, eNoN

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:)

      gnEL = lM%gnEl
      eNoN = lM%eNoN

      ALLOCATE(part(gnEl))
      part = 0
      nEl  = 0
      DO e=1, gnEl
         i = 0
         DO a=1, eNoN
            Ac = lM%gIEN(a,e)
            IF (lM%gpN(Ac) .EQ. cm%tF()) i = i + 1
         END DO
         IF (i .EQ. eNoN) THEN
            nEl = nEl + 1
            part(e) = 1
         END IF
      END DO
      DEALLOCATE(lM%gpN)
!      WRITE(*,'(A)') "CPU "//STR(cm%tF())//" nEl = "//STR(nEl)

      ALLOCATE(gE(nEl))
      nEl = 0
      DO e=1, gnEl
         IF (part(e) .GT. 0) THEN
            nEl = nEl + 1
            gE(nEl) = e
         END IF
      END DO
      DEALLOCATE(part)

      RETURN
      END SUBROUTINE DISTRE
!--------------------------------------------------------------------
!     Create list of connected/adjacent elements for old/source mesh
      SUBROUTINE GETADJESRC(lM, kneList)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), ALLOCATABLE :: kneList(:,:)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: a, b, e, Ac, i, j, maxKNE

      INTEGER(KIND=IKIND), ALLOCATABLE :: nL(:), tmpList(:,:)

!     First get elements around all nodes
      ALLOCATE(nL(lM%nNo))
      nL = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            nL(Ac) = nL(Ac) + 1
         END DO
      END DO

      maxKNE = MAX(MAXVAL(nL), 0)
      ALLOCATE(tmpList(maxKNE, lM%nNo))
      nL(:)  = 0
      tmpList(:,:) = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            nL(Ac) = nL(Ac) + 1
            tmpList(nL(Ac), Ac) = e
         END DO
      END DO
      b = 2*maxKNE

!     Now get elements around each element and avoid redundancy
 001  b = b + maxKNE
      DEALLOCATE(nL)
      ALLOCATE(nL(lM%nEl))
      IF (ALLOCATED(kneList)) DEALLOCATE(kneList)
      ALLOCATE(kneList(b, lM%nEl))
      kneList = 0
      nL = 0
      DO e=1, lM%nEl
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            DO i=1, maxKNE
               IF (tmpList(i,Ac) .NE. 0) THEN
                  flag = .TRUE.
                  DO j=1, nL(e)
                     IF (kneList(j,e) .EQ. tmpList(i,Ac)) THEN
                        flag = .FALSE.
                        EXIT
                     END IF
                  END DO
                  IF (flag) THEN
                     nL(e) = nL(e) + 1
                     IF (nL(e) .GE. b) GOTO 001
                     kneList(nL(e), e) = tmpList(i,Ac)
                  END IF
               ELSE
                  EXIT
               END IF
            END DO
         END DO
      END DO
      DEALLOCATE(tmpList)
      maxKNE = MAX(MAXVAL(nL), 0)
!      WRITE(*,'(A)') "CPU "//STR(cm%tF())//", maxKNE = "//STR(maxKNE)

      ALLOCATE(tmpList(b,lM%nEl))
      tmpList = kneList
      DEALLOCATE(kneList)
      ALLOCATE(kneList(maxKNE, lM%nEl))
      DO e=1, lM%nEl
         kneList(:, e) = tmpList(1:maxKNE, e)
      END DO
      DEALLOCATE(nL, tmpList)

      RETURN
      END SUBROUTINE GETADJESRC
!--------------------------------------------------------------------
!     Create list of connected/adjacent nodes for new/target mesh
      SUBROUTINE GETADJNTGT(lM, nNo, nEl, gN, gE, knnList)
      USE COMMOD
      USE UTILMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: nNo, nEl, gN(nNo), gE(nEl)
      INTEGER(KIND=IKIND), ALLOCATABLE :: knnList(:,:)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, Ec, i, j, n, maxKNN

      INTEGER(KIND=IKIND), ALLOCATABLE :: lN(:), nL(:), tmpList(:,:)

      ALLOCATE(lN(lM%gnNo))
      lN = 0
      DO a=1, nNo
         Ac = gN(a)
         lN(Ac) = a
      END DO

      ALLOCATE(nL(nNo))
      nL = 0
      DO e=1, nEl
         Ec = gE(e)
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,Ec)
            i = lN(Ac)
            DO b=1, lM%eNoN
               IF (b .EQ. a) CYCLE
               Bc = lM%gIEN(b,Ec)
               nL(i) = nL(i) + 1
            END DO
         END DO
      END DO
      maxKNN = MAX(MAXVAL(nL), 0)

      ALLOCATE(tmpList(maxKNN, nNo))
      tmpList = 0
      nL = 0
      DO e=1, nEl
         Ec = gE(e)
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,Ec)
            i = lN(Ac)
            DO b=1, lM%eNoN
               IF (a .EQ. b) CYCLE
               Bc = lM%gIEN(b,Ec)
               j = lN(Bc)
               flag = .TRUE.
               DO n=1, nL(i)
                  IF (tmpList(n,i) .EQ. j) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  nL(i) = nL(i) + 1
                  tmpList(nL(i), i) = j
               END IF
            END DO
         END DO
      END DO
      maxKNN = MAX(MAXVAL(nL), 0)
!      WRITE(*,'(A)') "CPU "//STR(cm%tF())//", maxKNN = "//STR(maxKNN)

      ALLOCATE(knnList(maxKNN, nNo))
      DO a=1, nNo
         knnList(:,a) = tmpList(1:maxKNN,a)
      END DO
      DEALLOCATE(tmpList, nL, lN)

      RETURN
      END SUBROUTINE GETADJNTGT
!--------------------------------------------------------------------
!     Find element in the old mesh for each node in the new mesh
      SUBROUTINE FINDN(xp, iM, Dg, eList, Ec, Nsf)
      USE COMMOD
      USE MATFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iM, eList(:)
      INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd+1), Dg(nsd,tnNo)
      REAL(KIND=RKIND), INTENT(OUT) :: Nsf(msh(iM)%eNoN)

      INTEGER(KIND=IKIND) :: e, a, i, j, ne, Ac

      REAL(KIND=RKIND), ALLOCATABLE :: Amat(:,:)

      ne = size(eList)
      ALLOCATE(Amat(nsd+1,msh(iM)%eNoN))
      Nsf(:) = 0._RKIND
      DO e=1, ne
         Ec = eList(e)
         IF (Ec .EQ. 0) EXIT
         Amat = 1._RKIND
         DO a=1, msh(iM)%eNoN
            Ac = msh(iM)%IEN(a,Ec)
            Amat(1:nsd,a) = x(:,Ac) + Dg(:,Ac)
         END DO
         Amat = MAT_INV(Amat, msh(iM)%eNoN)
         a = 0
         DO i=1, nsd+1
            Nsf(i) = 0._RKIND
            DO j=1, msh(iM)%eNoN
               Nsf(i) = Nsf(i) + Amat(i,j)*Xp(j)
            END DO
            IF (Nsf(i).GT.-1.E-14_RKIND .AND.
     2          Nsf(i).LT.(1._RKIND+1.E-14_RKIND)) THEN
               a = a + 1
            END IF
         END DO
         IF (a .EQ. nsd+1) THEN
            RETURN
         END IF
      END DO
      Ec = 0
      Nsf(:) = 0._RKIND

      RETURN
      END SUBROUTINE FINDN
!--------------------------------------------------------------------


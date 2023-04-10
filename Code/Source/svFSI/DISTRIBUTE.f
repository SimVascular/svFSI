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
!     This routine partitions the mesh, distributes values from
!     single processor to multiple processors and prepares the problem
!     to be lunched with several processors.
!
!--------------------------------------------------------------------

      SUBROUTINE DISTRIBUTE
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: a, e, i, Ac, iEq, iM, iFa, iBf

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), gmtl(:)
      REAL(KIND=RKIND4), ALLOCATABLE :: iWgt(:)
      REAL(KIND=RKIND), ALLOCATABLE :: wgt(:,:), wrk(:), tmpX(:,:),
     2   tmpD(:,:,:)
      TYPE(mshType), ALLOCATABLE :: tMs(:)

!     Preparing IO incase of error or warning. I'm keeping dbg channel
!     closed is slave processors. Warning is closed only if it is
!     closed in master
      IF (.NOT.resetSim) THEN
         CALL cm%bcast(pClr)
         CALL cm%bcast(appPath)
         CALL cm%bcast(wrn%oTS)
         wrn%oTF = wrn%oTS

!     Constructing data structures one by one
         CALL cm%bcast(nMsh)
         CALL cm%bcast(nsd)
         CALL cm%bcast(rmsh%isReqd)
      END IF
      CALL cm%bcast(gtnNo)

      IF (cm%slv()) ALLOCATE(msh(nMsh))

!     tMs is a temporary variable to keep fa%gN of the old meshes.
!     wgt and wrk are the assigned portion of each mesh to the each
!     processor.
      ALLOCATE(tMs(nMsh), wgt(nMsh,cm%np()), wrk(nMsh), iWgt(cm%np()),
     2   gmtl(gtnNo))

!     Here is rough estimation of how each mesh should be splited
!     between processors
      wrk = REAL(msh%gnNo, KIND=RKIND)/REAL(gtnNo, KIND=RKIND)
      CALL cm%bcast(wrk)
      CALL SPLITJOBS(nMsh, cm%np(), wgt, wrk)

!     First partitioning the meshes
!     gmtl:  gtnNo --> tnNo. What is gmtl? I think gmtl is a mapping from
!     global nodal index to nodal index local to this proc. If this proc does
!     not own an element with global nodal index Ac, then gmtl(Ac) = 0.
      tnNo = 0
      gmtl = 0
      IF (cm%seq()) THEN
         tnNo = gtnNo
         ALLOCATE(ltg(tnNo))
         DO a=1, tnNo
            ltg(a) = a
         END DO
      END IF
      DO iM=1, nMsh
         dbg = "Partitioning mesh "//iM
         iWgt = REAL(wgt(iM,:)/SUM(wgt(iM,:)), KIND=RKIND4)
!        what is gmtl??
         CALL PARTMSH(msh(iM), gmtl, cm%np(), iWgt)
      END DO

!     Setting gtl pointer in case that it is needed and mapping IEN
      DO iM=1, nMsh
         IF (ALLOCATED(msh(iM)%lN)) DEALLOCATE(msh(iM)%lN)
         ALLOCATE(msh(iM)%lN(tnNo))
         msh(iM)%lN = 0
         DO a=1, msh(iM)%nNo
            Ac             = msh(iM)%gN(a)
            msh(iM)%lN(Ac) = a
         END DO
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac               = msh(iM)%IEN(a,e)
               msh(iM)%IEN(a,e) = msh(iM)%gN(Ac)
            END DO
         END DO
      END DO
      IF (cm%seq()) THEN
!        Rearrange body force structure, if necessary
         DO iEq=1, nEq
            DO iBf=1, eq(iEq)%nBf
               IF (ALLOCATED(eq(iEq)%bf(iBf)%bx)) THEN
                  i  = eq(iEq)%bf(iBf)%dof
                  iM = eq(iEq)%bf(iBf)%iM
                  ALLOCATE(tmpX(i,gtnNo))
                  tmpX = eq(iEq)%bf(iBf)%bx
                  DEALLOCATE(eq(iEq)%bf(iBf)%bx)
                  ALLOCATE(eq(iEq)%bf(iBf)%bx(i,msh(iM)%nNo))
                  DO a=1, msh(iM)%nNo
                     Ac = msh(iM)%gN(a)
                     eq(iEq)%bf(iBf)%bx(:,a) = tmpX(:,Ac)
                  END DO
                  DEALLOCATE(tmpX)
               ELSE IF (ALLOCATED(eq(iEq)%bf(iBf)%bm)) THEN
                  i  = eq(iEq)%bf(iBf)%bm%dof
                  a  = eq(iEq)%bf(iBf)%bm%nTP
                  iM = eq(iEq)%bf(iBf)%iM
                  ALLOCATE(tmpD(i,gtnNo,a))
                  tmpD = eq(iEq)%bf(iBf)%bm%d
                  DEALLOCATE(eq(iEq)%bf(iBf)%bm%d)
                  ALLOCATE(eq(iEq)%bf(iBf)%bm%d(i,msh(iM)%nNo,a))
                  eq(iEq)%bf(iBf)%bm%d = 0._RKIND
                  DO i=1, eq(iEq)%bf(iBf)%bm%nTP
                     DO a=1, msh(iM)%nNo
                        Ac = msh(iM)%gN(a)
                        eq(iEq)%bf(iBf)%bm%d(:,a,i) = tmpD(:,Ac,i)
                     END DO
                  END DO
                  DEALLOCATE(tmpD)
               END IF
            END DO
         END DO
         RETURN
      END IF

!     Partitioning the faces
      DO iM=1, nMsh
         ALLOCATE(tMs(iM)%fa(msh(iM)%nFa))
         DO iFa=1, msh(iM)%nFa
            CALL PARTFACE(msh(iM), msh(iM)%fa(iFa), tMs(iM)%fa(iFa),
     2         gmtl)
         END DO
      END DO

!     Sending data from read by master in READFILES to slaves
      IF (.NOT.resetSim) THEN
         CALL cm%bcast(nsymd)
         CALL cm%bcast(stopTrigName)
         CALL cm%bcast(iniFilePath)
         CALL cm%bcast(stFileName)
         CALL cm%bcast(stFileFlag)
         CALL cm%bcast(stFileIncr)
         CALL cm%bcast(stFileRepl)
         CALL cm%bcast(saveIncr)
         CALL cm%bcast(saveATS)
         CALL cm%bcast(saveAve)
         CALL cm%bcast(saveVTK)
         CALL cm%bcast(bin2VTK)
         CALL cm%bcast(mvMsh)
         CALL cm%bcast(nITS)
         CALL cm%bcast(nTS)
         CALL cm%bcast(startTS)
         CALL cm%bcast(nEq)
         CALL cm%bcast(dt)
         CALL cm%bcast(zeroAve)
         CALL cm%bcast(cmmInit)
         CALL cm%bcast(cmmVarWall)
         CALL cm%bcast(shlEq)
         CALL cm%bcast(pstEq)
         CALL cm%bcast(sstEq)
         CALL cm%bcast(cepEq)
         CALL cm%bcast(ecCpld)
         IF (rmsh%isReqd) THEN
            CALL cm%bcast(rmsh%method)
            CALL cm%bcast(rmsh%freq)
            CALL cm%bcast(rmsh%cpVar)
            IF (cm%slv()) THEN
               ALLOCATE(rmsh%maxEdgeSize(nMsh))
               rmsh%minDihedAng = 0._RKIND
               rmsh%maxRadRatio = 0._RKIND
            END IF
            call cm%bcast(rmsh%maxEdgeSize)
         END IF
         CALL cm%bcast(iCntct)
         IF (iCntct) THEN
            CALL cm%bcast(cntctM%cType)
            CALL cm%bcast(cntctM%k)
            CALL cm%bcast(cntctM%c)
            CALL cm%bcast(cntctM%h)
            CALL cm%bcast(cntctM%al)
            CALL cm%bcast(cntctM%p)
            CALL cm%bcast(cntctM%Rin)
            CALL cm%bcast(cntctM%Rout)
            CALL cm%bcast(cntctM%gap)
         END IF
         CALL cm%bcast(ibFlag)
         IF (ibFlag) CALL DISTIB()
         CALL cm%bcast(nXion)
      END IF

!     Distributing X to processors
      IF (cm%mas()) THEN
         ALLOCATE(tmpX(nsd,gtnNo))
         tmpX = x
         DEALLOCATE(x)
      ELSE
         ALLOCATE(tmpX(0,0))
      END IF
      ALLOCATE(x(nsd,tnNo))
      x = LOCAL(tmpX)
      DEALLOCATE(tmpX)

!     Distributing lM%dmnId if present to processors
      flag = ALLOCATED(dmnId)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(part(gtnNo))
            part = dmnId
            DEALLOCATE(dmnId)
         ELSE
            ALLOCATE(part(0))
         END IF
         ALLOCATE(dmnId(tnNo))
         dmnId = LOCAL(part)
         DEALLOCATE(part)
      END IF

!     Distribute prestress (pS0) to processors
      flag = ALLOCATED(pS0)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(nsymd,gtnNo))
            tmpX = pS0
            DEALLOCATE(pS0)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         ALLOCATE(pS0(nsymd,tnNo))
         pS0 = LOCAL(tmpX)
         DEALLOCATE(tmpX)
      END IF

!     Distribute initial flow quantities to processors
      flag = ALLOCATED(Pinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(1,gtnNo))
            tmpX(1,:) = Pinit
            DEALLOCATE(Pinit)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         IF (ALLOCATED(wgt)) DEALLOCATE(wgt)
         ALLOCATE(wgt(1,tnNo), Pinit(tnNo))
         wgt = LOCAL(tmpX)
         Pinit(:) = wgt(1,:)
         DEALLOCATE(tmpX, wgt)
      END IF

      flag = ALLOCATED(Vinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(nsd,gtnNo))
            tmpX = Vinit
            DEALLOCATE(Vinit)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         ALLOCATE(Vinit(nsd,tnNo))
         Vinit = LOCAL(tmpX)
         DEALLOCATE(tmpX)
      END IF

      flag = ALLOCATED(Dinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(nsd,gtnNo))
            tmpX = Dinit
            DEALLOCATE(Dinit)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         ALLOCATE(Dinit(nsd,tnNo))
         Dinit = LOCAL(tmpX)
         DEALLOCATE(tmpX)
      END IF

!     And distributing eq to processors
      IF (cm%slv()) ALLOCATE(eq(nEq))
      DO iEq=1, nEq
         CALL DISTEQ(eq(iEq), tMs, gmtl)
         dbg = "Distributed equation "//iEq
      END DO

!     For CMM initialization
      flag = ALLOCATED(cmmBdry)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(part(gtnNo))
            part = cmmBdry
            DEALLOCATE(cmmBdry)
         ELSE
            ALLOCATE(part(0))
         END IF
         ALLOCATE(cmmBdry(tnNo))
         cmmBdry = LOCAL(part)
         DEALLOCATE(part)
      END IF

!     For CMM variable wall properties
      flag = ALLOCATED(varWallProps)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(2,gtnNo))
            tmpX = varWallProps
            DEALLOCATE(varWallProps)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         ALLOCATE(varWallProps(2,tnNo))
         varWallProps = LOCAL(tmpX)
         DEALLOCATE(tmpX)
      END IF

!     Communicating cplBC data
      CALL cm%bcast(cplBC%nFa)
      CALL cm%bcast(cplBC%schm)
      CALL cm%bcast(cplBC%useGenBC)
      IF (cplBC%useGenBC) THEN
         IF (cm%slv()) THEN
            cplBC%nX = 0
            ALLOCATE(cplBC%xo(cplBC%nX))
         END IF
      ELSE
         CALL cm%bcast(cplBC%nX)
         IF (.NOT.ALLOCATED(cplBC%xo)) ALLOCATE(cplBC%xo(cplBC%nX))
         IF (cplBC%nX .NE. 0) CALL cm%bcast(cplBC%xo)
      END IF
      CALL cm%bcast(cplBC%initRCR)

      DO iM=1, nMsh
         CALL DESTROY(tMs(iM))
      END DO
      DEALLOCATE(tMs)

      RETURN
      END SUBROUTINE DISTRIBUTE
!####################################################################
!     This routine distributes immersed boundary data structures
      SUBROUTINE DISTIB()
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iM, iFa

      IF (cm%slv()) ALLOCATE(ib)

      CALL cm%bcast(ib%mthd)
      CALL cm%bcast(ib%cpld)
      CALL cm%bcast(ib%intrp)

      CALL cm%bcast(ib%nMsh)
      CALL cm%bcast(ib%tnNo)
      IF (cm%slv()) THEN
         ALLOCATE(ib%msh(ib%nMsh))
         ALLOCATE(ib%x(nsd,ib%tnNo))
      END IF
      CALL cm%bcast(ib%x)

      DO iM=1, ib%nMsh
         CALL DISTIBMSH(ib%msh(iM))
         DO iFa=1, ib%msh(iM)%nFa
            CALL DISTIBFa(ib%msh(iM), ib%msh(iM)%fa(iFa))
         END DO
      END DO

      RETURN
      END SUBROUTINE DISTIB
!--------------------------------------------------------------------
      SUBROUTINE DISTIBMSH(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM

      LOGICAL fnFlag
      INTEGER(KIND=IKIND) i, insd

      CALL cm%bcast(lM%lShpF)
      CALL cm%bcast(lM%lShl)
      CALL cm%bcast(lM%lFib)
      CALL cm%bcast(lM%eType)
      CALL cm%bcast(lM%eNoN)
      CALL cm%bcast(lM%gnEl)
      CALL cm%bcast(lM%gnNo)
      CALL cm%bcast(lM%nFa)
      CALL cm%bcast(lM%nFs)
      CALL cm%bcast(lM%nG)
      CALL cm%bcast(lM%vtkType)
      CALL cm%bcast(lM%nFn)
      CALL cm%bcast(lM%scF)
      CALL cm%bcast(lM%dx)
      CALL cm%bcast(lM%name)

      fnFlag = ALLOCATED(lM%fN)
      CALL cm%bcast(fnFlag)

      IF (cm%slv()) THEN
         lM%nNo = lM%gnNo
         lM%nEl = lM%gnEl
         ALLOCATE(lM%gN(lM%nNo))
         ALLOCATE(lM%lN(ib%tnNo))
         ALLOCATE(lM%IEN(lM%eNoN, lM%nEl))
         ALLOCATE(lM%eId(lM%nEl))
         ALLOCATE(lM%fa(lM%nFa))
         IF (fnFlag) ALLOCATE(lM%fN(lM%nFn*nsd,lM%nEl))
         CALL SELECTELE(lM)
      END IF
      CALL cm%bcast(lM%gN)
      CALL cm%bcast(lM%lN)
      CALL cm%bcast(lM%IEN)
      CALL cm%bcast(lM%eId)
      IF (fnFlag) CALL cm%bcast(lM%fN)

      IF (lM%eType .EQ. eType_NRB) THEN
         CALL cm%bcast(lM%nSl)
         insd = nsd
         IF (lM%lShl) insd = nsd - 1
         IF (cm%slv()) THEN
            ALLOCATE(lM%nW(lM%gnNo))
            ALLOCATE(lM%INN(insd,lM%gnEl))
            ALLOCATE(lM%bs(insd))
         END IF
         CALL cm%bcast(lM%nW)
         CALL cm%bcast(lM%INN)
         DO i=1, insd
            CALL cm%bcast(lM%bs(i)%n)
            CALL cm%bcast(lM%bs(i)%nG)
            CALL cm%bcast(lM%bs(i)%nEl)
            CALL cm%bcast(lM%bs(i)%nSl)
            CALL cm%bcast(lM%bs(i)%p)
            IF (cm%slv()) ALLOCATE(lM%bs(i)%xi(lM%bs(i)%n))
            CALL cm%bcast(lM%bs(i)%xi)
            lM%bs(i)%nNo = lM%bs(i)%n - lM%bs(i)%p - 1
         END DO
      END IF

      RETURN
      END SUBROUTINE DISTIBMSH
!--------------------------------------------------------------------
      SUBROUTINE DISTIBFa(lM, lFa)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      CALL cm%bcast(lFa%d)
      CALL cm%bcast(lFa%eNoN)
      CALL cm%bcast(lFa%nEl)
      CALL cm%bcast(lFa%nNo)
      CALL cm%bcast(lFa%name)
      IF (cm%slv()) THEN
         ALLOCATE(lFa%IEN(lFa%eNoN,lFa%nEl))
         ALLOCATE(lFa%gN(lFa%nNo))
         ALLOCATE(lFa%lN(ib%tnNo))
         ALLOCATE(lFa%gE(lFa%nEl))
         CALL SELECTELEB(lM, lFa)
      END IF
      CALL cm%bcast(lFa%IEN)
      CALL cm%bcast(lFa%gN)
      CALL cm%bcast(lFa%lN)
      CALL cm%bcast(lFa%gE)

      RETURN
      END SUBROUTINE DISTIBFa
!####################################################################
!     This routine distributes equations between processors
      SUBROUTINE DISTEQ(lEq, tMs, gmtl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: gmtl(gtnNo)
      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(mshType), INTENT(IN) :: tMs(nMsh)

      INTEGER(KIND=IKIND) iDmn, iOut, iBc, iBf

!     Distribute equation parameters
      CALL cm%bcast(lEq%nOutput)
      CALL cm%bcast(lEq%coupled)
      CALL cm%bcast(lEq%maxItr)
      CALL cm%bcast(lEq%minItr)
      CALL cm%bcast(lEq%roInf)
      CALL cm%bcast(lEq%phys)
      CALL cm%bcast(lEq%nDmn)
      CALL cm%bcast(lEq%nBc)
      CALL cm%bcast(lEq%nBf)
      CALL cm%bcast(lEq%tol)
      CALL cm%bcast(lEq%absTol)
      CALL cm%bcast(lEq%useTLS)
      CALL cm%bcast(lEq%assmTLS)
      IF (ibFlag) THEN
         CALL cm%bcast(lEq%nDmnIB)
         CALL cm%bcast(lEq%nBcIB)
      END IF

!     Distribute linear solver settings
      CALL cm%bcast(lEq%FSILS%foC)
      CALL cm%bcast(lEq%FSILS%LS_type)
      CALL cm%bcast(lEq%FSILS%RI%relTol)
      CALL cm%bcast(lEq%FSILS%GM%relTol)
      CALL cm%bcast(lEq%FSILS%CG%relTol)
      CALL cm%bcast(lEq%FSILS%RI%absTol)
      CALL cm%bcast(lEq%FSILS%GM%absTol)
      CALL cm%bcast(lEq%FSILS%CG%absTol)
      CALL cm%bcast(lEq%FSILS%RI%mItr)
      CALL cm%bcast(lEq%FSILS%GM%mItr)
      CALL cm%bcast(lEq%FSILS%CG%mItr)
      CALL cm%bcast(lEq%FSILS%RI%sD)
      CALL cm%bcast(lEq%FSILS%GM%sD)
      CALL cm%bcast(lEq%FSILS%CG%sD)

      CALL cm%bcast(lEq%ls%LS_Type)
      CALL cm%bcast(lEq%ls%PREC_Type)
      CALL cm%bcast(lEq%ls%relTol)
      CALL cm%bcast(lEq%ls%absTol)
      CALL cm%bcast(lEq%ls%mItr)
      CALL cm%bcast(lEq%ls%sD)

!     Distribute domain properties
      IF (cm%slv()) ALLOCATE(lEq%dmn(lEq%nDmn))
      DO iDmn=1, lEq%nDmn
         CALL cm%bcast(lEq%dmn(iDmn)%phys)
         CALL cm%bcast(lEq%dmn(iDmn)%Id)
         CALL cm%bcast(lEq%dmn(iDmn)%prop)
         IF (lEq%dmn(iDmn)%phys .EQ. phys_CEP) THEN
            CALL cm%bcast(lEq%dmn(iDmn)%cep%cepType)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%fpar_in)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%nX)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%nG)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%nFn)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%imyo)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%dt)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Ksac)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Diso)
            IF (cm%slv()) THEN
               ALLOCATE(lEq%dmn(iDmn)%cep%Dani(lEq%dmn(iDmn)%cep%nFn))
            END IF
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Dani)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Istim%Ts)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Istim%Td)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Istim%CL)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%Istim%A)
            CALL cm%bcast(lEq%dmn(iDmn)%cep%odes%tIntType)
            IF (lEq%dmn(iDmn)%cep%odes%tIntType .EQ. tIntType_CN2) THEN
               CALL cm%bcast(lEq%dmn(iDmn)%cep%odes%maxItr)
               CALL cm%bcast(lEq%dmn(iDmn)%cep%odes%absTol)
               CALL cm%bcast(lEq%dmn(iDmn)%cep%odes%relTol)
            END IF
         END IF

         IF ((lEq%dmn(iDmn)%phys .EQ. phys_struct)  .OR.
     2       (lEq%dmn(iDmn)%phys .EQ. phys_ustruct) .OR.
     3       (lEq%dmn(iDmn)%phys .EQ. phys_shell)) THEN
            CALL DIST_MATCONSTS(lEq%dmn(iDmn)%stM)
         END IF

         IF (ecCpld) THEN
            CALL DIST_ECMODEL(lEq%dmn(iDmn)%ec)
         END IF

         IF ((lEq%dmn(iDmn)%phys .EQ. phys_fluid)  .OR.
     2       (lEq%dmn(iDmn)%phys .EQ. phys_stokes) .OR.
     3       (lEq%dmn(iDmn)%phys .EQ. phys_CMM .AND. .NOT.cmmInit))THEN
            CALL DIST_VISCMODEL(lEq%dmn(iDmn)%visc)
         END IF
      END DO

      IF (ibFlag) THEN
         IF (cm%slv()) ALLOCATE(lEq%dmnIB(lEq%nDmnIB))
         DO iDmn=1, lEq%nDmnIB
            CALL cm%bcast(lEq%dmnIB(iDmn)%phys)
            CALL cm%bcast(lEq%dmnIB(iDmn)%Id)
            CALL cm%bcast(lEq%dmnIB(iDmn)%prop)
            CALL DIST_MATCONSTS(lEq%dmnIB(iDmn)%stM)
         END DO
      END IF

!     Distribute output parameters
      IF (cm%slv()) ALLOCATE(lEq%output(lEq%nOutput))
      DO iOut=1, lEq%nOutput
         CALL cm%bcast(lEq%output(iOut)%wtn)
         CALL cm%bcast(lEq%output(iOut)%grp)
         CALL cm%bcast(lEq%output(iOut)%o)
         CALL cm%bcast(lEq%output(iOut)%l)
         CALL cm%bcast(lEq%output(iOut)%name)
      END DO

!     Distribute BC information
      IF (cm%slv()) ALLOCATE(lEq%bc(lEq%nBc))
      DO iBc=1, lEq%nBc
         CALL DISTBC(lEq%bc(iBc), tMs, gmtl)
      END DO

      IF (ibFlag) THEN
         IF (cm%slv()) ALLOCATE(lEq%bcIB(lEq%nBcIB))
         DO iBc=1, lEq%nBcIB
            CALL DISTBCIB(lEq%bcIB(iBc))
         END DO
      END IF

!     Distribute BF information
      IF (cm%slv()) ALLOCATE(lEq%bf(lEq%nBf))
      DO iBf=1, lEq%nBf
         CALL DISTBF(lEq%bf(iBf))
      END DO

      RETURN
      END SUBROUTINE DISTEQ
!--------------------------------------------------------------------
!     This routine distributes the BCs between processors
      SUBROUTINE DISTBC(lBc, tMs, gmtl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: gmtl(gtnNo)
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(mshType), INTENT(IN) :: tMs(nMsh)

      LOGICAL flag
      INTEGER(KIND=IKIND) i, j, iDof, nTp, nNo, a, b, Ac, iM, iFa

      INTEGER(KIND=IKIND), ALLOCATABLE :: tmpI(:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmp(:)

      CALL cm%bcast(lBc%cplBCptr)
      CALL cm%bcast(lBc%bType)
      IF (cm%slv()) ALLOCATE(lBc%eDrn(nsd), lBc%h(nsd))
      CALL cm%bcast(lBc%eDrn)
      CALL cm%bcast(lBc%iFa)
      CALL cm%bcast(lBc%iM)
      CALL cm%bcast(lBc%r)
      CALL cm%bcast(lBc%g)
      CALL cm%bcast(lBc%k)
      CALL cm%bcast(lBc%c)
      CALL cm%bcast(lBc%h)
      CALL cm%bcast(lBc%weakDir)
      CALL cm%bcast(lBc%tauB)
      CALL cm%bcast(lBc%flwP)
      CALL cm%bcast(lBc%rbnN)
      IF (BTEST(lBc%bType,bType_RCR)) THEN
         CALL cm%bcast(lBc%RCR%Rp)
         CALL cm%bcast(lBc%RCR%C)
         CALL cm%bcast(lBC%RCR%Rd)
         CALL cm%bcast(lBC%RCR%Pd)
         CALL cm%bcast(lBC%RCR%Xo)
      END IF

!     Communicating time-dependent BC data
      flag = ALLOCATED(lBc%gt)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBc%gt)
         CALL cm%bcast(lBc%gt%lrmp)
         CALL cm%bcast(lBc%gt%d)
         CALL cm%bcast(lBc%gt%n)
         j = lBc%gt%d
         i = lBc%gt%n
         IF (cm%slv()) THEN
            ALLOCATE(lBc%gt%qi(j))
            ALLOCATE(lBc%gt%qs(j))
            ALLOCATE(lBc%gt%r(j,i))
            ALLOCATE(lBc%gt%i(j,i))
         END IF
         CALL cm%bcast(lBc%gt%ti)
         CALL cm%bcast(lBc%gt%T)
         CALL cm%bcast(lBc%gt%qi)
         CALL cm%bcast(lBc%gt%qs)
         CALL cm%bcast(lBc%gt%r)
         CALL cm%bcast(lBc%gt%i)
      END IF

!     Communicating moving BC data
      flag = ALLOCATED(lBc%gm)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBc%gm)
         CALL cm%bcast(lBc%gm%period)
!     Communication the %t data
         CALL cm%bcast(lBc%gm%nTP)
         CALL cm%bcast(lBc%gm%dof)
         nTp  = lBc%gm%nTP
         iDof = lBc%gm%dof
         IF (cm%slv()) ALLOCATE(lBc%gm%t(nTp))
         CALL cm%bcast(lBc%gm%t)

         nNo  = tMs(lBc%iM)%fa(lBc%iFa)%nNo
         a    = nTp*iDof*nNo
!     Allocating the container and copying the nodes which belong to
!     this processor
         ALLOCATE(tmp(a))
         IF (cm%mas()) THEN
            tmp = RESHAPE(lBc%gm%d,(/a/))
            DEALLOCATE(lBc%gm%d)
         END IF

         CALL cm%bcast(tmp)
!     This is the new number of nodes
         a = msh(lBc%iM)%fa(lBc%iFa)%nNo
         ALLOCATE(lBc%gm%d(iDof,a,nTp))
         b = 0
         DO a=1, nNo
            Ac = tMs(lBc%iM)%fa(lBc%iFa)%gN(a)
            Ac = gmtl(Ac)
            IF (Ac .NE. 0) THEN
               b = b + 1
               DO i=1, nTP
                  j = iDof*((i-1)*nNo + a - 1)
                  lBc%gm%d(:,b,i) = tmp(j+1:j+iDof)
               END DO
            END IF
         END DO
      END IF

!     Communicating profile data
      flag = ALLOCATED(lBc%gx)
      CALL cm%bcast(flag)
      IF (flag) THEN
         nNo = tMs(lBc%iM)%fa(lBc%iFa)%nNo
         IF (ALLOCATED(tmp)) DEALLOCATE(tmp)
         ALLOCATE(tmp(nNo))
         IF (cm%mas()) THEN
            tmp = lBc%gx
            DEALLOCATE(lBc%gx)
         END IF
         CALL cm%bcast(tmp)
!     This is the new number of nodes
         a = msh(lBc%iM)%fa(lBc%iFa)%nNo
         ALLOCATE(lBc%gx(a))
         b = 0
         DO a=1, nNo
            Ac = tMs(lBc%iM)%fa(lBc%iFa)%gN(a)
            Ac = gmtl(Ac)
            IF (Ac .NE. 0) THEN
               b = b + 1
               lBc%gx(b) = tmp(a)
            END IF
         END DO
         DEALLOCATE(tmp)
      END IF

!     Communicating and reordering master node data for clamped Neu BC
      IF (BTEST(lBc%bType,bType_clmpd)) THEN
         CALL cm%bcast(lBc%masN)
         iM   = lBc%iM
         iFa  = lBc%iFa
         nNo  = msh(iM)%fa(iFa)%nNo
         flag = nNo .NE. 0
!        Action performed only on the processes that share the face
         IF (flag) THEN
            Ac = lBc%masN
            a  = gmtl(Ac)
!        The process that owns the node is set to be master. For other
!        processes that share the face but do not own the node, we add
!        this node as a ghost master node. ltg pointer is reset.
            IF (a .NE. 0) THEN
               lBc%masN = a
            ELSE
               nNo  = tnNo
               tnNo = tnNo + 1
!              Remap ltg
               IF (ALLOCATED(tmpI)) DEALLOCATE(tmpI)
               ALLOCATE(tmpI(nNo))
               tmpI = ltg
               DEALLOCATE(ltg)
               ALLOCATE(ltg(tnNo))
               ltg(1:nNo) = tmpI(:)
               DEALLOCATE(tmpI)
               ltg(tnNo) = Ac
               lBc%masN  = tnNo

!              Add the ghost master node to the face data structure
               nNo = msh(iM)%fa(iFa)%nNo
               msh(iM)%fa(iFa)%nNo = nNo + 1
               ALLOCATE(tmpI(nNo))
               tmpI(:) = msh(iM)%fa(iFa)%gN(:)
               DEALLOCATE(msh(iM)%fa(iFa)%gN)
               ALLOCATE(msh(iM)%fa(iFa)%gN(msh(iM)%fa(iFa)%nNo))
               msh(iM)%fa(iFa)%gN(1:nNo) = tmpI(:)
               msh(iM)%fa(iFa)%gN(nNo+1) = tnNo
            END IF
         ELSE
!        Zero out master node if not part of the face
            lBc%masN = 0
         END IF
      ELSE
         lBc%masN = 0
      END IF

      RETURN
      END SUBROUTINE DISTBC
!--------------------------------------------------------------------
!     This routine distributes the BCs on immersed surfaces
      SUBROUTINE DISTBCIB(lBc)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc

      LOGICAL flag
      INTEGER(KIND=IKIND) i, j, iDof, nTp, nNo, a

      REAL(KIND=RKIND), ALLOCATABLE :: tmp(:)

      CALL cm%bcast(lBc%bType)
      IF (cm%slv()) ALLOCATE(lBc%eDrn(nsd))
      CALL cm%bcast(lBc%eDrn)
      CALL cm%bcast(lBc%iFa)
      CALL cm%bcast(lBc%iM)
      CALL cm%bcast(lBc%r)
      CALL cm%bcast(lBc%g)
      CALL cm%bcast(lBc%weakDir)
      CALL cm%bcast(lBc%flwP)

!     Communicating time-dependant BC data
      flag = ALLOCATED(lBc%gt)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBc%gt)
         CALL cm%bcast(lBc%gt%lrmp)
         CALL cm%bcast(lBc%gt%d)
         CALL cm%bcast(lBc%gt%n)
         j = lBc%gt%d
         i = lBc%gt%n
         IF (cm%slv()) THEN
            ALLOCATE(lBc%gt%qi(j))
            ALLOCATE(lBc%gt%qs(j))
            ALLOCATE(lBc%gt%r(j,i))
            ALLOCATE(lBc%gt%i(j,i))
         END IF
         CALL cm%bcast(lBc%gt%ti)
         CALL cm%bcast(lBc%gt%T)
         CALL cm%bcast(lBc%gt%qi)
         CALL cm%bcast(lBc%gt%qs)
         CALL cm%bcast(lBc%gt%r)
         CALL cm%bcast(lBc%gt%i)
      END IF

!     Communicating moving BC data
      flag = ALLOCATED(lBc%gm)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBc%gm)
         CALL cm%bcast(lBc%gm%period)
!     Communication the %t data
         CALL cm%bcast(lBc%gm%nTP)
         CALL cm%bcast(lBc%gm%dof)
         nTp  = lBc%gm%nTP
         iDof = lBc%gm%dof
         IF (cm%slv()) ALLOCATE(lBc%gm%t(nTp))
         CALL cm%bcast(lBc%gm%t)

         nNo  = ib%msh(lBc%iM)%fa(lBc%iFa)%nNo
         a    = nTp*iDof*nNo
!     Allocating the container and copying the nodes which belong to
!     this processor
         ALLOCATE(tmp(a))
         IF (cm%mas()) THEN
            tmp = RESHAPE(lBc%gm%d,(/a/))
            DEALLOCATE(lBc%gm%d)
         ELSE
            ALLOCATE(lBc%gm%d(iDof,nNo,nTp))
         END IF

         CALL cm%bcast(tmp)
         DO a=1, nNo
            DO i=1, nTP
               j = iDof*((i-1)*nNo + a - 1)
               lBc%gm%d(:,a,i) = tmp(j+1:j+iDof)
            END DO
         END DO
         DEALLOCATE(tmp)
      END IF

!     Communicating profile data
      flag = ALLOCATED(lBc%gx)
      CALL cm%bcast(flag)
      IF (flag) THEN
         nNo = ib%msh(lBc%iM)%fa(lBc%iFa)%nNo
         IF (cm%slv()) ALLOCATE(lBc%gx(nNo))
         CALL cm%bcast(lBc%gx)
      END IF

      RETURN
      END SUBROUTINE DISTBCIB
!--------------------------------------------------------------------
!     This routine distributes the BF between processors
      SUBROUTINE DISTBF(lBf)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bfType), INTENT(INOUT) :: lBf

      LOGICAL flag
      INTEGER(KIND=IKIND) a, i, j, Ac, iM, idof, nTP

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:), tmpD(:,:,:)

      CALL cm%bcast(lBf%bType)
      CALL cm%bcast(lBf%dof)
      CALL cm%bcast(lBf%iM)
      iM = lBf%iM

!     Steady value
      flag = ALLOCATED(lBf%b)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBf%b(lBf%dof))
         CALL cm%bcast(lBf%b)
      END IF

!     Communicating spatially dependent BF
      flag = ALLOCATED(lBf%bx)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpX(lBf%dof,gtnNo))
            tmpX = lBf%bx
            DEALLOCATE(lBf%bx)
         ELSE
            ALLOCATE(tmpX(0,0))
         END IF
         ALLOCATE(lBf%bx(lBf%dof,tnNo))
         lBf%bx = LOCAL(tmpX)
         DEALLOCATE(tmpX)
         ALLOCATE(tmpX(lBf%dof,tnNo))
         tmpX = lBf%bx
         DEALLOCATE(lBf%bx)
         ALLOCATE(lBf%bx(lBf%dof,msh(iM)%nNo))
         lBf%bx = 0._RKIND
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            lBf%bx(:,a) = tmpX(:,Ac)
         END DO
         DEALLOCATE(tmpX)
      END IF

!     Communicating time-dependent BF data
      flag = ALLOCATED(lBf%bt)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBf%bt)
         CALL cm%bcast(lBf%bt%lrmp)
         CALL cm%bcast(lBf%bt%d)
         CALL cm%bcast(lBf%bt%n)
         j = lBf%bt%d
         i = lBf%bt%n
         IF (cm%slv()) THEN
            ALLOCATE(lBf%bt%qi(j))
            ALLOCATE(lBf%bt%qs(j))
            ALLOCATE(lBf%bt%r(j,i))
            ALLOCATE(lBf%bt%i(j,i))
         END IF
         CALL cm%bcast(lBf%bt%ti)
         CALL cm%bcast(lBf%bt%T)
         CALL cm%bcast(lBf%bt%qi)
         CALL cm%bcast(lBf%bt%qs)
         CALL cm%bcast(lBf%bt%r)
         CALL cm%bcast(lBf%bt%i)
      END IF

!     Communicating moving BF data
      flag = ALLOCATED(lBf%bm)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (cm%slv()) ALLOCATE(lBf%bm)
         CALL cm%bcast(lBf%bm%period)
!     Communication the %t data
         CALL cm%bcast(lBf%bm%nTP)
         CALL cm%bcast(lBf%bm%dof)
         nTP  = lBf%bm%nTP
         idof = lBf%bm%dof
         IF (cm%slv()) ALLOCATE(lBf%bm%t(nTP))
         CALL cm%bcast(lBf%bm%t)

         IF (cm%mas()) THEN
            ALLOCATE(tmpD(lBf%bm%dof,gtnNo,lBf%bm%nTP))
            tmpD = lBf%bm%d
            DEALLOCATE(lBf%bm%d)
         ELSE
            ALLOCATE(tmpD(0,0,0))
         END IF
         ALLOCATE(tmpX(lBf%bm%dof,tnNo),
     2      lBf%bm%d(lBf%bm%dof,msh(iM)%nNo,lBf%bm%nTP))
         lBf%bm%d = 0._RKIND
         DO i=1, lBf%bm%nTP
            tmpX = LOCAL(tmpD(:,:,i))
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               lBf%bm%d(:,a,i) = tmpX(:,Ac)
            END DO
         END DO
         DEALLOCATE(tmpX, tmpD)
      END IF

      RETURN
      END SUBROUTINE DISTBF
!--------------------------------------------------------------------
!     This subroutine distributes constants and parameters of the
!     structural constitutive model to all processes
      SUBROUTINE DIST_MATCONSTS(lStM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(stModelType), INTENT(INOUT) :: lStM

      INTEGER(KIND=IKIND) i, j

      CALL cm%bcast(lStM%volType)
      CALL cm%bcast(lStM%Kpen)
      CALL cm%bcast(lStM%isoType)
      CALL cm%bcast(lStM%C01)
      CALL cm%bcast(lStM%C10)
      CALL cm%bcast(lStM%a)
      CALL cm%bcast(lStM%b)
      CALL cm%bcast(lStM%aff)
      CALL cm%bcast(lStM%bff)
      CALL cm%bcast(lStM%ass)
      CALL cm%bcast(lStM%bss)
      CALL cm%bcast(lStM%afs)
      CALL cm%bcast(lStM%bfs)
      CALL cm%bcast(lStM%kap)
      CALL cm%bcast(lStM%khs)
      CALL cm%bcast(lStM%a0)
      CALL cm%bcast(lStM%b1)
      CALL cm%bcast(lStM%b2)
      CALL cm%bcast(lStM%mu0)

!     Distribute fiber stress
      CALL cm%bcast(lStM%Tf%fType)
      IF (BTEST(lStM%Tf%fType, bType_std)) THEN
         CALL cm%bcast(lStM%Tf%g)
      ELSE IF (BTEST(lStM%Tf%fType, bType_ustd)) THEN
         CALL cm%bcast(lStM%Tf%gt%lrmp)
         CALL cm%bcast(lStM%Tf%gt%d)
         CALL cm%bcast(lStM%Tf%gt%n)
         j = lStM%Tf%gt%d
         i = lStM%Tf%gt%n
         IF (cm%slv()) THEN
            ALLOCATE(lStM%Tf%gt%qi(j))
            ALLOCATE(lStM%Tf%gt%qs(j))
            ALLOCATE(lStM%Tf%gt%r(j,i))
            ALLOCATE(lStM%Tf%gt%i(j,i))
         END IF
         CALL cm%bcast(lStM%Tf%gt%ti)
         CALL cm%bcast(lStM%Tf%gt%T)
         CALL cm%bcast(lStM%Tf%gt%qi)
         CALL cm%bcast(lStM%Tf%gt%qs)
         CALL cm%bcast(lStM%Tf%gt%r)
         CALL cm%bcast(lStM%Tf%gt%i)
      END IF
      CALL cm%bcast(lStM%Tf%eta_s)

      RETURN
      END SUBROUTINE DIST_MATCONSTS
!--------------------------------------------------------------------
!     This subroutine distributes constants and parameters of the
!     excitation-contraction coupling model
      SUBROUTINE DIST_ECMODEL(lEc)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eccModelType), INTENT(INOUT) :: lEc

      CALL cm%bcast(lEc%astress)
      CALL cm%bcast(lEc%astrain)
      CALL cm%bcast(lEc%asnType)
      CALL cm%bcast(lEc%k)
      CALL cm%bcast(lEc%caCpld)
      CALL cm%bcast(lEc%fpar_in)
      IF (.NOT.lEc%caCpld) THEN
         CALL cm%bcast(lEc%odes%tIntType)
         CALL cm%bcast(lEc%dt)
         CALL cm%bcast(lEc%dType)
         IF (BTEST(lEc%dType, bType_std)) THEN
            CALL cm%bcast(lEc%Ya)
         ELSE IF (BTEST(lEc%dType, bType_ustd)) THEN
            CALL cm%bcast(lEc%Yat%lrmp)
            CALL cm%bcast(lEc%Yat%d)
            CALL cm%bcast(lEc%Yat%n)
            IF (cm%slv()) THEN
               ALLOCATE(lEc%Yat%qi(lEc%Yat%d))
               ALLOCATE(lEc%Yat%qs(lEc%Yat%d))
               ALLOCATE(lEc%Yat%r(lEc%Yat%d,lEc%Yat%n))
               ALLOCATE(lEc%Yat%i(lEc%Yat%d,lEc%Yat%n))
            END IF
            CALL cm%bcast(lEc%Yat%ti)
            CALL cm%bcast(lEc%Yat%T)
            CALL cm%bcast(lEc%Yat%qi)
            CALL cm%bcast(lEc%Yat%qs)
            CALL cm%bcast(lEc%Yat%r)
            CALL cm%bcast(lEc%Yat%i)
         END IF
      END IF
      CALL cm%bcast(lEc%eta_s)

      RETURN
      END SUBROUTINE DIST_ECMODEL
!--------------------------------------------------------------------
!     This subroutine distributes constants and parameters of the
!     fluid viscosity constitutive model to all processes
      SUBROUTINE DIST_VISCMODEL(lVis)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(viscModelType), INTENT(INOUT) :: lVis

      CALL cm%bcast(lVis%viscType)
      CALL cm%bcast(lVis%mu_i)
      CALL cm%bcast(lVis%mu_o)
      CALL cm%bcast(lVis%lam)
      CALL cm%bcast(lVis%a)
      CALL cm%bcast(lVis%n)

      RETURN
      END SUBROUTINE DIST_VISCMODEL
!####################################################################
!     This is for partitioning a single mesh
      SUBROUTINE PARTMSH(lM, gmtl, nP, wgt)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: nP
      REAL(KIND=RKIND4), INTENT(IN) :: wgt(nP)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: gmtl(gtnNo) ! gtnNo = global total number of nodes
      TYPE(mshType), INTENT(INOUT) :: lM

      LOGICAL :: flag, fnFlag
      INTEGER(KIND=MPI_OFFSET_KIND) :: idisp
      INTEGER(KIND=IKIND) :: i, a, Ac, e, Ec, edgecut, nEl, nNo, eNoN,
     2   eNoNb, ierr, fid, SPLIT, insd, nFn
      CHARACTER(LEN=stdL) fTmp

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), gPart(:),
     2   tempIEN(:,:), gtlPtr(:), sCount(:), disp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:), tmpFn(:,:)

!     If sequential, copy global parameters as local parameters and return
      IF (cm%seq()) THEN
         lM%nEl = lM%gnEl ! number of elements in this mesh
         lM%nNo = lM%gnNo ! number of nodes in this mesh
         ALLOCATE(lM%IEN(lM%eNoN,lM%nEl), lM%eDist(0:cm%np()))
         lM%IEN      = lM%gIEN ! IEN array (maps element e (of mesh) and node a to global node A (of mesh))
         lM%eDist(0) = 0
         lM%eDist(1) = lM%gnEl ! element distribution array, lM%eDist(i) represents first element which belong to cm%id()=i
         ALLOCATE(lM%otnIEN(lM%nEl))! old to new IEN array (maps ??)
         DO e=1, lM%nEl
            lM%otnIEN(e) = e
         END DO
         ALLOCATE(lM%iGC(lM%nEl))
         lM%iGC = 0
         RETURN
      END IF

!     Sending data from read by master in READFILES to slaves
      CALL cm%bcast(lM%lShpF)
      CALL cm%bcast(lM%lShl)
      CALL cm%bcast(lM%lFib)
      CALL cm%bcast(lM%eType)
      CALL cm%bcast(lM%eNoN)
      CALL cm%bcast(lM%nFa)
      CALL cm%bcast(lM%nFs)
      CALL cm%bcast(lM%nG)
      CALL cm%bcast(lM%gnEl)
      CALL cm%bcast(lM%gnNo)
      CALL cm%bcast(lM%name)
      CALL cm%bcast(lM%nFn)
      CALL cm%bcast(lM%scF)
      nFn = lM%nFn

      insd = nsd
      IF (lM%lShl) insd = nsd - 1
      IF (lM%lFib) insd = 1

      eNoN = lM%eNoN ! number of nodes in the element
      IF (cm%slv()) THEN
         CALL SELECTELE(lM)
         ALLOCATE(lM%gIEN(0,0), lM%fa(lM%nFa))
      END IF
      ALLOCATE(sCount(cm%np()), disp(cm%np()), lM%eDist(0:cm%np()))

!     And distributing bs for NURBS
      IF (lM%eType .EQ. eType_NRB) THEN
         IF (cm%slv()) ALLOCATE(lM%bs(insd))
         CALL cm%bcast(lM%nSl)
         DO i=1, insd
            CALL cm%bcast(lM%bs(i)%n)
            CALL cm%bcast(lM%bs(i)%nG)
            CALL cm%bcast(lM%bs(i)%nEl)
            CALL cm%bcast(lM%bs(i)%nSl)
            CALL cm%bcast(lM%bs(i)%p)
            IF (cm%slv()) ALLOCATE(lM%bs(i)%xi(lM%bs(i)%n))
            CALL cm%bcast(lM%bs(i)%xi)
            lM%bs(i)%nNo = lM%bs(i)%n - lM%bs(i)%p - 1
         END DO

         a = lM%bs(2)%nEl
         IF (insd .EQ. 3) a = a*lM%bs(3)%nEl
         DO i=0, cm%np()
            lM%eDist(i) = a*NINT(SUM(wgt(1:i))*lM%bs(1)%nEl, KIND=IKIND)
            IF (lM%eDist(i) .GT. lM%gnEl) lM%eDist(i) = lM%gnEl
         END DO
      ELSE
!     A draft of splitting the mesh between processors
!     lM%eDist(i) represents first element which belong to cm%id()=i
         DO i=0, cm%np() ! Loop over procs
            lM%eDist(i) = NINT(SUM(wgt(1:i))*lM%gnEl, KIND=IKIND)
            IF (lM%eDist(i) .GT. lM%gnEl) lM%eDist(i) = lM%gnEl
         END DO
      END IF
      lM%eDist(cm%np()) = lM%gnEl

!     Count nodes assigned to each proc
      DO i=1, cm%np()
         disp(i)   = lM%eDist(i-1)*eNoN
         sCount(i) = lM%eDist(i)*eNoN - disp(i)
      END DO
!     Count elements assigned to each proc and allocate space for part array
!     part array maps an element to the proc it belongs to
      nEl = lM%eDist(cm%id() + 1) - lM%eDist(cm%id())
      idisp = lM%eDist(cm%id())*SIZEOF(nEl)
      ALLOCATE(part(nEl))

!     Read domain partition from file if desired (if remesh is used)
      fTmp = TRIM(appPath)//".partitioning_"//TRIM(lM%name)//".bin"
      flag = .FALSE.
      IF (rmsh%isReqd) INQUIRE(FILE=TRIM(fTmp), EXIST=flag)
      IF (lM%eType .EQ. eType_NRB) THEN
         part = cm%id()
      ELSE IF (flag .AND. .NOT.resetSim) THEN
         std = " Reading partition data from file"
         CALL MPI_FILE_OPEN(cm%com(), TRIM(fTmp), MPI_MODE_RDONLY,
     2      MPI_INFO_NULL, fid, ierr)
         CALL MPI_FILE_SET_VIEW(fid, idisp, mpint, mpint, 'native',
     2      MPI_INFO_NULL, ierr)
         CALL MPI_FILE_READ(fid, part, nEl, mpint,
     2      MPI_STATUS_IGNORE, ierr)
         CALL MPI_FILE_CLOSE(fid, ierr)
      ELSE
!        lM%IEN maps node a in element e (on this proc) to global node Ac on this proc?
         ALLOCATE(lM%IEN(eNoN,nEl))
!     Scattering the lM%gIEN array to processors. Splits IEN array according to
!     rough partition computed above. Places data in lM%IEN. I think this step
!     is needed to do the final partition in SPLIT() below
         CALL MPI_SCATTERV(lM%gIEN, sCount, disp, mpint, lM%IEN,
     2      nEl*eNoN, mpint, master, cm%com(), ierr)

!     This is to get eNoNb
         SELECT CASE (lM%eType)
         CASE(eType_TET4)
            eNoNb = 3
         CASE(eType_TET10)
            eNoNb = 6
         CASE(eType_HEX8)
            eNoNb = 4
         CASE(eType_HEX20)
            eNoNb = 8
         CASE(eType_HEX27)
            eNoNb = 9
         CASE(eType_WDG)
            eNoNb = 3
         CASE(eType_TRI3)
            eNoNb = 2
         CASE(eType_TRI6)
            eNoNb = 3
         CASE(eType_QUD4)
            eNoNb = 2
         CASE(eType_QUD8)
            eNoNb = 3
         CASE(eType_QUD9)
            eNoNb = 3
         CASE(eType_LIN1)
            eNoNb = 1
         CASE(eType_LIN2)
            eNoNb = 1
         CASE DEFAULT
            err = "Undefined element type"
         END SELECT
!     The output of this process is "part" array. part(i) says
!     which processor element "i" belongs to
!     Doing partitioning, using ParMetis
         edgecut = SPLIT(nEl, eNoN, eNoNb, lM%IEN, cm%np(), lM%eDist,
     2      wgt, part)
         IF (edgecut .EQ. 0) THEN
c            wrn = " ParMETIS failed to partition the mesh"
            part = cm%id()
         ELSE IF (edgecut .GT. 0) THEN
            std = " ParMETIS partitioned the mesh by cutting "//
     2         STR(edgecut)//" elements"
!     LT 0 is for the case that all elements reside in one processor
         END IF
         DEALLOCATE(lM%IEN) ! Destroy lM%IEN
         IF (rmsh%isReqd) THEN
            std = " Writing partition data to file"
            CALL MPI_FILE_OPEN(cm%com(), TRIM(fTmp), MPI_MODE_WRONLY +
     2         MPI_MODE_CREATE, MPI_INFO_NULL, fid, ierr)
            CALL MPI_FILE_SET_VIEW(fid, idisp, mpint, mpint, 'native',
     2         MPI_INFO_NULL, ierr)
            CALL MPI_FILE_WRITE(fid, part, nEl, mpint,
     2         MPI_STATUS_IGNORE, ierr)
            CALL MPI_FILE_CLOSE(fid, ierr)
         END IF
      END IF

!     Compute first element belonging to each proc (disp(i) and number 
!     of elements belonging to each proc (sCount(i))
      DO i=1, cm%np() ! Loop over procs
         disp(i)   = lM%eDist(i-1)
         sCount(i) = lM%eDist(i) - disp(i)
      END DO

!     Gathering the parts inside master, part(e) is equal to the
!     cm%id() that the element e belong to
      IF (cm%mas()) THEN
         ALLOCATE(gPart(lM%gnEl))
      ELSE
         ALLOCATE(gPart(0))
      END IF
!     gpart is a global version of part in which processor p = gpart(e)
!     is the owner of element "e"
      CALL MPI_GATHERV(part, nEl, mpint, gPart, sCount, disp, mpint,
     2   master, cm%com(), ierr)

!     Destroy part array in each proc. master holds gPart, which contains all
!     the info from part array in each proc. 
      DEALLOCATE(part)

!     From gPart, fill in eDist array, where eDist(i) = global index of first 
!     element that belongs to proc i
      IF (cm%mas()) THEN
         sCount = 0
         DO e=1, lM%gnEl ! Loop over global elements in mesh
            sCount(gPart(e) + 1) = sCount(gPart(e) + 1) + 1
         END DO
!        sCount(i) = num elements on proc i
         DO i=1, cm%np()
            lM%eDist(i) = lM%eDist(i-1) + sCount(i)
         END DO
!        eDist(i) = global index of first element that belong to proc i

         ALLOCATE(tempIEN(eNoN,lM%gnEl), lM%otnIEN(lM%gnEl))
!     Making the lM%IEN array in order, based on the cm%id() number in
!     master. lM%otnIEN maps old IEN order to new IEN order.
!     Why does the IEN order change? I think to facilitate scattering later.
         disp = 0
         DO e=1, lM%gnEl ! Loop over global elements in mesh
!           gPart(e) = proc that owns e
!           eDist(i) = first element that belong to proc i (minus 1?)
            Ec = lM%eDist(gPart(e)) + 1 
            lM%eDist(gPart(e)) = Ec
            tempIEN(:,Ec) = lM%gIEN(:,e)
            lM%otnIEN(e) = Ec
         END DO
         lM%gIEN = tempIEN
         lM%eDist(0) = 0
         DO i=1, cm%np()
            lM%eDist(i) = lM%eDist(i-1) + sCount(i)
         END DO

!     This it to distribute eId, if allocated
         flag = .FALSE.
         IF (ALLOCATED(lM%eId)) THEN
            flag = .TRUE.
            ALLOCATE(part(lM%gnEl))
            DO e=1, lM%gnEl
               Ec = lM%otnIEN(e)
               part(Ec) = lM%eId(e)
            END DO
            DEALLOCATE(lM%eId)
         END IF

!     This it to distribute fN, if allocated
         fnFlag = .FALSE.
         IF (ALLOCATED(lM%fN)) THEN
            fnFlag = .TRUE.
            ALLOCATE(tmpFn(nFn*nsd,lM%gnEl))
            tmpFn = 0._RKIND
            DO e=1, lM%gnEl
               Ec = lM%otnIEN(e)
               tmpFn(:,Ec) = lM%fN(:,e)
            END DO
            DEALLOCATE(lM%fN)
         END IF
      ELSE
         ALLOCATE(lM%otnIEN(0))
      END IF
      DEALLOCATE(gPart)

      CALL cm%bcast(flag)
      CALL cm%bcast(fnFlag)
      CALL cm%bcast(lM%eDist)

      nEl = lM%eDist(cm%id() + 1) - lM%eDist(cm%id())
      lM%nEl = nEl
      ALLOCATE(lM%IEN(eNoN,nEl), lM%iGC(nEl))
      lM%iGC = 0

!     Communicating eId, if neccessary
      IF (flag) THEN
         ALLOCATE(lM%eId(nEl))
         IF (.NOT.ALLOCATED(part)) ALLOCATE(part(0))
         DO i=1, cm%np()
            disp(i)   = lM%eDist(i-1)
            sCount(i) = lM%eDist(i) - disp(i)
         END DO
         CALL MPI_SCATTERV(part, sCount, disp, mpint, lM%eId, nEl,
     2      mpint, master, cm%com(), ierr)
         DEALLOCATE(part)
      END IF

!     Communicating fN, if neccessary
      IF (fnFlag) THEN
         ALLOCATE(lM%fN(nFn*nsd,nEl))
         IF (.NOT.ALLOCATED(tmpFn)) ALLOCATE(tmpFn(0,0))
         DO i=1, cm%np()
            disp(i)   = lM%eDist(i-1)*nFn*nsd
            sCount(i) = lM%eDist(i)*nFn*nsd - disp(i)
         END DO
         CALL MPI_SCATTERV(tmpFn, sCount, disp, mpreal, lM%fN,
     2      nEl*nFn*nsd, mpreal, master, cm%com(), ierr)
         DEALLOCATE(tmpFn)
      END IF

!     Now scattering the sorted lM%IEN to all processors. Data placed in
!     lM%IEN. Now, lM%IEN maps node a of element e (local to this proc and mesh)
!     to global node A on mesh
      IF (.NOT.ALLOCATED(tempIEN)) ALLOCATE(tempIEN(0,0))
      DO i=1, cm%np()
         disp(i)   = lM%eDist(i-1)*eNoN
         sCount(i) = lM%eDist(i)*eNoN - disp(i)
      END DO
      CALL MPI_SCATTERV(tempIEN, sCount, disp, mpint, lM%IEN, nEl*eNoN,
     2   mpint, master, cm%com(), ierr)
      DEALLOCATE(tempIEN)

!     Constructing the initial global to local pointer
!     lM%IEN: eNoN,nEl --> gnNo. Maps node a of element e (local to proc and mesh) to global node A on mesh
!     gtlPtr: gnNo     --> nNo. Maps global node A on mesh to local node a on proc and mesh
!     lM%IEN: eNoN,nEl --> nNo. Reconstructing IEN to map node a of element e (local to proc and mesh) to local node a on proc and mesh
      ALLOCATE(gtlPtr(lM%gnNo)) ! global number of nodes on this mesh
      nNo    = 0
      gtlPtr = 0
      DO e=1, nEl ! Loop over num elements on this proc on this mesh
         DO a=1, eNoN ! Loop over nodes in element
            Ac = lM%IEN(a,e) ! Get global node Ac on this mesh. This is the global node ID in Paraview
!           If Ac has not be seen yet, add it to this proc. Increment the number
!           of nodes on this proc, and update gtlPtr
            IF (gtlPtr(Ac) .EQ. 0) THEN
               nNo = nNo + 1
               gtlPtr(Ac) = nNo ! Maps global node Ac on mesh to local node a on this proc (on same mesh)
            END IF
!           Now, lM%IEN maps node a of element e (local to this proc and mesh) 
!           to local node a on this proc (on same mesh)
            lM%IEN(a,e) = gtlPtr(Ac)
         END DO
      END DO
      lM%nNo = nNo
      IF (cm%slv()) ALLOCATE(lM%gN(lM%gnNo))
      CALL cm%bcast(lM%gN)

!     Use gtlptr to distribute lM%tmX, if allocated
      flag = ALLOCATED(lM%tmX)
      CALL cm%bcast(flag)
      IF (flag) THEN
         ALLOCATE(tmpR(lM%gnNo))
         IF (cm%mas()) THEN
            tmpR = lM%tmX
            DEALLOCATE(lM%tmX)
         END IF

         CALL cm%bcast(tmpR)

         ALLOCATE(lM%tmX(lM%nNo))
         lM%tmX = 0._RKIND
         DO Ac=1, lM%gnNo
            a = gtlptr(Ac)
            IF (a .NE. 0) THEN
               lM%tmX(a) = tmpR(Ac)
            END IF
         END DO
         DEALLOCATE(tmpR)
      END IF

!     lM%gN: gnNo --> gtnNo. Maps global node index on this mesh to global node index across all meshes
!     part:  nNo  --> gtnNo. Maps local node index on this proc on this mesh to global
      ALLOCATE(part(nNo))
      DO Ac=1, lM%gnNo ! Loop over nodes on this mesh
         a = gtlPtr(Ac) ! get local node index on this proc and mesh
         ! If a is non-zero, then this proc owns global node Ac, so assign the
         ! corresponding global total node gN(Ac) to part(a).
         ! At this point, gN(Ac) maps global node index on this mesh to
         ! global node index across all meshes
         IF (a .NE. 0) part(a) = lM%gN(Ac)
      END DO
!     mapping and converting other parameters.
!     I will use an upper bound for gPart as a container for ltg,
!     since there can be repeated nodes. gPart is just a temp variable.
!     gmtl:  gtnNo --> tnNo. Maps global node index across all meshes to local node index (across all meshes or parts of meshes belong to this proc). gmtl is zero global node Ac does not belong to this proc
!     gPart: tnNo  --> gtnNo. Maps local node index on this proc (across multiple parts of meshes) to global node index across all meshes. Inverse of gmtl
!     ltg:   tnNo  --> gtnNo. How is this different from gPart, if at all? I think they contain the same info, gPart is just temporary as we assemble the info
!     lM%gN: nNo   --> tnNo. Changing gN to map local node index on this proc and this mesh to local node index on this proc across all meshes or parts of meshes belonging to this proc
      DEALLOCATE(lM%gN)
      ALLOCATE(gPart(tnNo+nNo), lM%gN(nNo))
!     Loop over nodes on this proc (potentially multiple meshes or parts of meshes).
!     If this is the first mesh, then tnNo = 0 and ltg is unfilled. So this loop
!     is only meaningful if we are partitioning the second mesh or later
      DO a=1, tnNo
         Ac       = ltg(a)
         gPart(a) = Ac
         gmtl(Ac) = a   ! Ac is a node index across all meshes and procs. a is a node index in this proc, which may be responsible for parts of multiple meshes
      END DO
      DO a=1, nNo ! Loop over nodes on this proc on this mesh
         Ac = part(a) ! part(a) gives the global node index across all meshes of node a on this proc on this mesh
         IF (gmtl(Ac) .EQ. 0) THEN ! If this node has not be reached yet, add it to this procs structures
            tnNo        = tnNo + 1 ! increment total nodes on this proc (across all meshes)
            gmtl(Ac)    = tnNo     ! Map the global node index across all meshes (Ac) to the recently added local node index on this proc across all meshes
            lM%gN(a)    = tnNo     ! Map the local node on this proc on this mesh to the recently added local node index on this proc across all meshes
            gPart(tnNo) = Ac       ! Map the local node index on this proc across all meshes to the global node index across all meshes (inverse of gmtl)
         ELSE
            lM%gN(a) = gmtl(Ac)
         END IF
      END DO
      IF (ALLOCATED(ltg)) DEALLOCATE(ltg)
      ALLOCATE(ltg(tnNo))
      ltg = gPart(1:tnNo) ! Set ltg array, which maps the local node index on this proc across all meshes to the global index across all meshes
      DEALLOCATE(gPart)

!     If neccessary communicate NURBS
      IF (lM%eType .EQ. eType_NRB) THEN
         ALLOCATE(tmpR(lM%gnNo))
         IF (cm%mas()) THEN
            tmpR = lM%nW
            DEALLOCATE(lM%nW)
         END IF
         CALL cm%bcast(tmpR)
         ALLOCATE(lM%nW(lM%nNo))
         DO Ac=1, lM%gnNo
            a = gtlPtr(Ac)
            IF (a .NE. 0) THEN
               lM%nW(a) = tmpR(Ac)
            END IF
         END DO
         DEALLOCATE(tmpR)

!     Distributing INN, using tempIEN as tmp array
         IF (cm%mas()) THEN
            ALLOCATE(tempIEN(insd,lM%gnEl))
            DO e=1, lM%gnEl
               Ec = lM%otnIEN(e)
               tempIEN(:,Ec) = lM%INN(:,e)
            END DO
            DEALLOCATE(lM%INN)
         ELSE
            ALLOCATE(tempIEN(0,0))
         END IF
         DO i=1, cm%np()
            disp(i)   = lM%eDist(i-1)*insd
            sCount(i) = lM%eDist(i)*insd - disp(i)
         END DO
         ALLOCATE(lM%INN(insd,nEl))
!     Now scattering the sorted lM%INN to all processors
         CALL MPI_SCATTERV(tempIEN, sCount, disp, mpint, lM%INN,
     2      nEl*insd, mpint, master, cm%com(), ierr)
      END IF

      RETURN
      END SUBROUTINE PARTMSH
!--------------------------------------------------------------------
!     This routine partitions the face based on the already partitioned
!     mesh
      SUBROUTINE PARTFACE(lM, lFa, gFa, gmtl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa, gFa
      INTEGER(KIND=IKIND), INTENT(INOUT) :: gmtl(gtnNo)

      INTEGER(KIND=IKIND) eNoNb, e, a, Ac, Ec, i, j, iM

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), ePtr(:)

!     Broadcasting the number of nodes and elements off to slaves and
!     populating gFa to all procs. gFa is a global face object that contains
!     info about the entire face, accessible by all procs.
      IF (cm%mas()) THEN
         gFa%d    = lFa%d
         gFa%eNoN = lFa%eNoN ! number of nodes in an element on this face
         gFa%iM   = lFa%iM   ! mesh index
         gFa%nEl  = lFa%nEl  ! number of elements on this face
         gFa%gnEl = lFa%gnEl ! global number of elements on this face
!        gFa%gnEL = gFa%nEl = lFa%gnEl, right?
         gFa%nNo  = lFa%nNo  ! number of nodes on this face
         gFa%virtual = lFa%virtual ! Is this face virtual
         gFa%capFaceName = lFa%capFaceName ! Name of capping face
         gFa%capFaceID = lFa%capFaceID ! ID number of capping face
         IF (rmsh%isReqd) ALLOCATE(gFa%gebc(1+gFa%eNoN,gFa%gnEl))
      ELSE
         IF (rmsh%isReqd) ALLOCATE(gFa%gebc(0,0))
      END IF
      CALL cm%bcast(gFa%d)
      CALL cm%bcast(gFa%eNoN)
      CALL cm%bcast(gFa%iM)
      CALL cm%bcast(gFa%nEl)
      CALL cm%bcast(gFa%gnEl)
      CALL cm%bcast(gFa%nNo)
      CALL cm%bcast(gFa%virtual)
      CALL cm%bcast(gFa%capFaceName)
      CALL cm%bcast(gFa%capFaceID)
!     Selects face mesh element type (linear triangle, quadratic quad, etc.)
!     Also sets up Gauss points and shape functions
!     All procs call this and individually populate their global face object gFa
      CALL SELECTELEB(lM, gFa) 

      eNoNb = gFa%eNoN ! number of nodes in an element on this face
      iM = gFa%iM      ! index of mesh that this face belongs to

!     Allocate connectivity array (IEN), global element IDs array (gE), 
!     global node IDs array (gN), and element pointer array (ePtr)
!     The IEN maps element e and element node a to global node A: IEN(a,e) = A
!     gE maps element e on this face to the global element index on the volume mesh
!     gN maps node a on this face to the global node index on the volume mesh
!     ePtr maps 
      ALLOCATE(gFa%IEN(eNoNb,gFa%nEl), gFa%gE(gFa%nEl), gFa%gN(gFa%nNo),
     2   ePtr(gFa%nEl))

!     If master proc, transfer local face structure to global face structure
!     and destroy local face structure
!     At this point, only master's gFa contains detailed info about the face
      IF (cm%mas()) THEN
         gFa = lFa
         CALL DESTROY(lFa)
      END IF

!     Set local face structure info according to broadcasted global face structure
      CALL cm%bcast(gFa%name)
      lFa%name = gFa%name ! Name of face
      lFa%d    = gFa%d
      lFa%eNoN = eNoNb ! Number of nodes in a face element
!     Selects face mesh element type and set up Gauss points and shape functions
!     for local face object (lFa). Before, we did it for gFa
      CALL SELECTELEB(lM, lFa)
      lFa%iM   = iM
      lFa%virtual = gFa%virtual
      lFa%capFaceName = gFa%capFaceName
      lFa%capFaceID = gFa%capFaceID

!     AB 7/7/22: If face is virtual, cannot partition it according to the already 
!     partitioned mesh. Instead, treat this separately.
      IF (lFa%virtual) THEN
         CALL PARTFACEV(lM, lFa, gFa, gmtl)
         RETURN
      END IF

!     Create and fill a structure (part) to communicate global face data among
!     all procs
      i = gFa%nEl*(2+eNoNb) + gFa%nNo
      ALLOCATE(part(i)) 
      IF (cm%mas()) THEN ! if master, fill part array with global face data
         DO e=1, gFa%nEl
            j  = (e-1)*(2+eNoNb) + 1
            Ec = gFa%gE(e)             ! Global element index. For virtual, this is 0 for all e. See READVTP.f
!           gIEN mapper from old to new. Was set in PARTMESH(). 
!           lM%otnIEN maps old IEN order to new IEN order. Indexing should start
!           at 1, so for virtual face (for which Ec = 0), otnIEN(Ec) = trash
            ePtr(e)   = lM%otnIEN(Ec)  
            part(j)   = Ec       
            part(j+1) = ePtr(e)        ! New element index for element e
            part(j+2:j+1+eNoNb) = gFa%IEN(:,e)  ! IEN array for this element
         END DO
         DO a=1, gFa%nNo
            j = gFa%nEl*(2+eNoNb) + a
            part(j) = gFa%gN(a) ! gN(a) maps node a on this face to global node Ac across all meshes
         END DO
      END IF

      CALL cm%bcast(part)
!     If slave proc, extract face data from part(:) and place into proper structures
      IF (cm%slv()) THEN
         DO e=1, gFa%nEl
            j = (e-1)*(2+eNoNb) + 1
            gFa%gE(e)    = part(j)     ! (old?) global element index
            ePtr(e)      = part(j+1)   ! New global element index 
            gFa%IEN(:,e) = part(j+2:j+1+eNoNb) ! IEN array for element e
         END DO
         DO a=1, gFa%nNo
            j = gFa%nEl*(2+eNoNb) + a
            gFa%gN(a) = part(j)        !
         END DO
      END IF
      DEALLOCATE(part)


!     Finding the number of lM%fas to allocate required space, also
!     maping global element number to processor element number
      lFa%nEl = 0
      DO e=1, gFa%nEl ! Loop over global number of elements on face
         Ec = ePtr(e) ! New global element index
         gFa%gE(e) = Ec ! Set gE as new global element index (changing gE(e) in general)
!        If this proc is responsible for this element, increment the number of 
!        elements in the local face object
         IF (Ec.LE.lM%eDist(cm%id()+1) .AND.
     2       Ec.GT.lM%eDist(cm%id()) ) THEN
            lFa%nEl = lFa%nEl + 1
         END IF
      END DO


      ALLOCATE(lFa%gE(lFa%nEl), lFa%IEN(eNoNb,lFa%nEl))

      lFa%nNo = 0
      DO a=1, gFa%nNo ! Loop over global number of nodes on face
!     gmtl maps global node index across all meshes to local index on 
!     this proc (potentially across multiple meshes)
!     For a global node A, gmtl is non-zero if node A lies in the mesh parition of this proc
         Ac = gmtl(gFa%gN(a)) 
         IF (Ac .NE. 0) THEN
            lFa%nNo = lFa%nNo + 1
         END IF
      END DO

      ALLOCATE(lFa%gN(lFa%nNo))

!     Time to form "face" structure in each processor
!     Only copying the element which belong to this processors
      j = 0
      DO e=1, gFa%nEl ! Loop over global number of elements
         Ec = gFa%gE(e) ! Get global element index (trash for virtual surface)
!        If the element index falls into the range of elements assigned to this
!        proc, add this face elements info to the local face structure on this proc
         IF (Ec.LE.lM%eDist(cm%id()+1) .AND.
     2       Ec.GT.lM%eDist(cm%id())) THEN
            j = j + 1
            lFa%gE(j) = Ec - lM%eDist(cm%id())
            DO a=1, eNoNb
               lFa%IEN(a,j) = gmtl(gFa%IEN(a,e))
            END DO
         END IF
      END DO

!     Analogously copying the nodes which belong to this processor
      j = 0
      DO a=1, gFa%nNo
         Ac = gmtl(gFa%gN(a))
         IF (Ac .NE. 0) THEN ! If Ac is 0, then the node does not belong to this proc
            j = j + 1
            lFa%gN(j) = Ac
         END IF
      END DO

      lFa%gnEl = gFa%gnEl
      IF (rmsh%isReqd) THEN
         IF(cm%mas()) THEN
            ALLOCATE(lFa%gebc(1+eNoNb,lFa%gnEl))
            DO e=1, gFa%gnEl
               lFa%gebc(1,e) = gFa%gebc(1,e)
               lFa%gebc(2:1+eNoNb,e) = gFa%gebc(2:1+eNoNb,e)
            END DO
         ELSE
            ALLOCATE(lFa%gebc(0,0))
         END IF
      END IF

      RETURN
      END SUBROUTINE PARTFACE
!--------------------------------------------------------------------
!     This routine partitions a virtual face. Since a virtual face does not
!     lie on the computational domain, it must be treated separately.
      SUBROUTINE PARTFACEV(lM, lFa, gFa, gmtl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa, gFa
      INTEGER(KIND=IKIND), INTENT(IN) :: gmtl(gtnNo)

      INTEGER(KIND=IKIND) eNoNb, e, a, Ac, Ec, i, j, iM

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), ePtr(:)


!     Define some useful variables
      eNoNb = gFa%eNoN ! number of nodes in an element on this face
      iM = gFa%iM      ! index of mesh that this face belongs to

      ALLOCATE(ePtr(gFa%nEl))

!     Create and fill a structure (part) to communicate global face data among
!     all procs. This synchronizes gFa across all procs
      i = gFa%nEl*(2+eNoNb) + gFa%nNo
      ALLOCATE(part(i)) 
!     If master proc, fill part array with global face data
      IF (cm%mas()) THEN 
         !  Add element information to part
         DO e=1, gFa%nEl
            j  = (e-1)*(2+eNoNb) + 1
            Ec = gFa%gE(e)             ! Global element index. For virtual, this is 0 for all e. See READVTP.f
!           gIEN mapper from old to new. Was set in PARTMESH(). 
!           lM%otnIEN maps old IEN order to new IEN order. Indexing should start
!           at 1, so for virtual face (for which Ec = 0), otnIEN(Ec) = trash
!           Instead, just set ePtr(e) = 0
!            ePtr(e)   = lM%otnIEN(Ec)  
            ePtr(e)   = 0
            part(j)   = Ec             ! 0 for virtual face 
            part(j+1) = ePtr(e)        ! New element index for element e
            part(j+2:j+1+eNoNb) = gFa%IEN(:,e)  ! IEN array for this element
         END DO

         ! Add nodal information to part
         DO a=1, gFa%nNo
            j = gFa%nEl*(2+eNoNb) + a
            part(j) = gFa%gN(a) ! gN(a) maps node a on this face to global node Ac across all meshes
         END DO
      END IF

!     Broadcast part array to all processors
      CALL cm%bcast(part)

!     If slave proc, extract face data from part(:) and place into proper structures
      IF (cm%slv()) THEN
         ! Extract element information from part
         DO e=1, gFa%nEl
            j = (e-1)*(2+eNoNb) + 1
            gFa%gE(e)    = part(j)     ! (old?) global element index. 0 for virtual face
            ePtr(e)      = part(j+1)   ! New global element index. 0 for virtual face
            gFa%IEN(:,e) = part(j+2:j+1+eNoNb) ! IEN array for element e
         END DO
         ! Extract nodal information from part
         DO a=1, gFa%nNo
            j = gFa%nEl*(2+eNoNb) + a
            gFa%gN(a) = part(j) 
         END DO
      END IF
      DEALLOCATE(part)

!     At this point, gFa is the same across all procs, and contains all
!     the face information

!     If face is virtual, then we can't partition the face according
!     to the element partition of the mesh, since the element of the
!     virtual face do not lie on the mesh. In this case, set lFa%nEl
!     equal to gFa%nEl for all procs. We'll deal with the parallelization
!     and communication when we need to integrate (in IntegS and IntegV)
!     Also, set gFa%gE(e) = -1 so we don't get confused
!     and behavior will be defined.
      lFa%nEl = gFa%nEl
!     This loop is not necessary, but including to match PARTFACE()
      DO e=1, gFa%nEl ! Loop over global number of elements on face
         Ec = ePtr(e) ! New global element index. 0 for virtual face
         gFa%gE(e) = Ec ! Set gE as new global element index (changing gE(e) in general).
         ! This was 0 before, and we're setting it to 0 again
      END DO

!     Allocate local face global element list and IEN array
      ALLOCATE(lFa%gE(lFa%nEl), lFa%IEN(eNoNb,lFa%nEl))

!     Compute the number of nodes on face that belong to this proc
      lFa%nNo = 0
      DO a=1, gFa%nNo ! Loop over global number of nodes on face
!     gmtl maps global node index across all meshes to local index on 
!     this proc (potentially across multiple meshes)
!     For a global node A, gmtl is non-zero if node A lies in the mesh parition of this proc
         Ac = gmtl(gFa%gN(a)) 
         IF (Ac .NE. 0) THEN
            lFa%nNo = lFa%nNo + 1
         END IF
      END DO

!     Allocate local face global node list
      ALLOCATE(lFa%gN(lFa%nNo))

!     Time to form "face" structure in each processor

!     To integrate later on, we need to fill in the lFa%IEN(a,e) array, which
!     gives the local (to processor) node index of the ath node of 
!     element e on the virtual face. To do this, we loop through
!     all elements and fill in what we can of IEN on each processor. On a single
!     processor, IEN will not be completely populated. For example, if a global
!     node A, which corresponds to node a on global face element e, does not
!     belong to processor 1 (this partition was done in PARTMESH and embodied in
!     gmtl), then lFa%IEN(a,e) = 0 on processor 1.
      DO e=1, lFa%nEl ! Loop over all elements of face (lFa%nEl = gFa%nEl)
!        Set lFa%gE(e) = 0, since for a virtual face, element e does not
!        lie on any element of the mesh
         lFa%gE(e) = 0
         DO a=1, eNoNb ! Loop over nodes on element
!           gmtl maps global node index (across all meshes and procs) to local
!           node index on this proc (across all parts of meshes belong to this proc)          
!           gmtl is zero if global node Ac does not belong to this proc
            lFa%IEN(a,e) = gmtl(gFa%IEN(a,e)) 
         END DO
      END DO

!     Analogously copying the nodes which belong to this processor
      j = 0
      DO a=1, gFa%nNo
         Ac = gmtl(gFa%gN(a))
         IF (Ac .NE. 0) THEN ! If Ac is 0, then the node does not belong to this proc
            j = j + 1
            lFa%gN(j) = Ac
         END IF
      END DO

      lFa%gnEl = gFa%gnEl
      IF (rmsh%isReqd) THEN
         IF(cm%mas()) THEN
            ALLOCATE(lFa%gebc(1+eNoNb,lFa%gnEl))
            DO e=1, gFa%gnEl
               lFa%gebc(1,e) = gFa%gebc(1,e)
               lFa%gebc(2:1+eNoNb,e) = gFa%gebc(2:1+eNoNb,e)
            END DO
         ELSE
            ALLOCATE(lFa%gebc(0,0))
         END IF
      END IF


      RETURN
      END SUBROUTINE PARTFACEV
!####################################################################

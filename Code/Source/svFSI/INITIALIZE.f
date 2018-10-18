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
!     This is to initialize or finalize svFSI variables/structures.
!
!--------------------------------------------------------------------

      SUBROUTINE INITIALIZE(timeP)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      REAL(KIND=8), INTENT(OUT) :: timeP(3)

      LOGICAL :: flag
      INTEGER :: i, iEq, iDmn, iM, ierr, nnz, gnnz
      CHARACTER(LEN=stdL) :: fTmp, sTmp
      REAL(KIND=8) :: am
      TYPE(FSILS_commuType) :: communicator

      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: s

      tDof     = 0
      flag     = .FALSE.
      dFlag    = .FALSE.
      pstEq    = .FALSE.
      cmmEq    = .FALSE.
      nFacesLS = SUM(eq%nBc)
      ibLSptr  = 0
      IF (ibFlag) THEN
         nFacesLS = nFacesLS + 1
         ibLSptr  = nFacesLS
      END IF
      IF (mvMsh) nFacesLS = nFacesLS + 1
      DO iEq=1, nEq
         DO i=1, eq(iEq)%nBc
            IF (eq(iEq)%bc(i)%weakDir) nFacesLS = nFacesLS - 1
         END DO
      END DO

      IF (ANY(eq(1)%bc%cplBCptr .NE. 0)) cplBC%coupled = .TRUE.

      DO iEq=1, nEq
!     This would be the default value of am for first order equations
         eq(iEq)%am = 5D-1*(3D0 - eq(iEq)%roInf) /
     2                (1D0 + eq(iEq)%roInf)
!     If the equation is second order, am is calculated as follows
         am = (2D0 - eq(iEq)%roInf)/(1D0 + eq(iEq)%roInf)
         SELECT CASE (eq(iEq)%phys)
         CASE (phys_fluid)
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'NS'
         CASE (phys_heatF)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'HF'
         CASE (phys_heatS)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'HS'
         CASE (phys_lElas)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'LE'
         CASE (phys_struct)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'ST'
         CASE (phys_ustruct)
            dFlag = .TRUE.
            flag  = .TRUE.
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'ST'
         CASE (phys_preSt)
            dFlag = .TRUE.
            pstEq = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'PS'
         CASE (phys_CMM)
            dFlag = .TRUE.
            cmmEq = .TRUE.
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'CM'
         CASE (phys_shell)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'SH'
         CASE (phys_FSI)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'FS'
         CASE (phys_mesh)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd
            eq(iEq)%am  = am
            eq(iEq)%sym = 'MS'
         CASE (phys_BBO)
            eq(iEq)%dof = nsd
            eq(iEq)%sym = 'BB'
         CASE (phys_CEP)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'EP'
         CASE DEFAULT
            err = "Equation type "//eq(iEq)%sym//" is not defined"
         END SELECT

         eq(iEq)%pNorm = HUGE(eq(iEq)%pNorm)
         eq(iEq)%af    = 1D0/(1D0 + eq(iEq)%roInf)
         eq(iEq)%beta  = 25D-2*(1D0 + eq(iEq)%am - eq(iEq)%af)**2D0
         eq(iEq)%gam   = 5D-1 + eq(iEq)%am - eq(iEq)%af
         eq(iEq)%s     = tDof + 1
         eq(iEq)%e     = tDof + eq(iEq)%dof
         tDof          = eq(iEq)%e
      END DO

      ierr = 0
      IF (dFlag .AND. .NOT.cmmEq) ierr = 1

      i = 0
      IF (cplBC%coupled) i = cplBC%nX

      stamp = (/cm%np(), nEq, nMsh, tnNo, i, tDof, ierr, version/)

!     Calculating the record length
      i = 2
      IF (dFlag .AND. .NOT.cmmEq) i = 3
      IF (pstEq) i = i + nstd
      i = 4*(1+SIZE(stamp)) + 8*(2+nEq+cplBC%nX+i*tDof*tnNo)
      IF (cm%seq()) THEN
         recLn = i
      ELSE
         CALL MPI_ALLREDUCE(i, recLn, 1, mpint, MPI_MAX, cm%com(), ierr)
      END IF

!     Initialize shell eIEN data structure. Used later in LHSA.
      DO iM=1, nMsh
         IF (msh(iM)%lShl) CALL SHLINI(msh(iM))
      END DO

!     Initialize tensor operations
      CALL TEN_INIT(nsd)

      std = " Constructing stiffness matrix sparse structure"
      CALL LHSA(nnz)

      gnnz = nnz
      CALL MPI_ALLREDUCE(nnz, gnnz, 1, mpint, MPI_SUM, cm%com(), ierr)
      std = " Total number of non-zeros in the LHS matrix: "//gnnz

!     Initialize FSILS structures
      IF (resetSim) THEN
         IF (communicator%foC) CALL FSILS_COMMU_FREE(communicator)
         IF (lhs%foC) CALL FSILS_LHS_FREE(lhs)
      END IF ! resetSim

      dbg = "Calling FSILS_COMMU_CREATE"
      CALL FSILS_COMMU_CREATE(communicator, cm%com())

      dbg = "Calling FSILS_LHS_CREATE"
      CALL FSILS_LHS_CREATE(lhs, communicator, gtnNo, tnNo, nnz, ltg,
     2   rowPtr, colPtr, nFacesLS)

!     Variable allocation and initialization
      ALLOCATE(Ao(tDof,tnNo), An(tDof,tnNo), Yo(tDof,tnNo),
     2   Yn(tDof,tnNo), Do(tDof,tnNo), Dn(tDof,tnNo))

!     Additional physics dependent variables
!     USTRUCT phys
      IF (flag) THEN
         ALLOCATE(Rd(nsd,tnNo), Kd((nsd+1)*nsd,nnz))
         IF (ustRd) THEN
            ALLOCATE(ADo(nsd,tnNo), ADn(nsd,tnNo))
         END IF
      END IF

!     PRESTRESS phys
      IF (pstEq) THEN
         IF (ALLOCATED(pS0)) err = "Prestress already allocated. "//
     2      "Correction needed"
         ALLOCATE(pS0(nstd,tnNo), pSn(nstd,tnNo), pSa(tnNo))
         pS0 = 0D0
         pSn = 0D0
         pSa = 0D0
      END IF

      IF (.NOT.resetSim) THEN
         IF (.NOT.ALLOCATED(rmsh%flag)) ALLOCATE(rmsh%flag(nMsh))
         rmsh%flag(:) = .FALSE.
         rmsh%fTS = rmsh%freq
         IF (rmsh%isReqd) THEN
            ALLOCATE(rmsh%A0(tDof,tnNo))
            ALLOCATE(rmsh%Y0(tDof,tnNo))
            ALLOCATE(rmsh%D0(tDof,tnNo))
            ALLOCATE(rmsh%iNorm(nEq))
         END IF

         INQUIRE (FILE=iniFilePath, EXIST=flag)
         IF (flag) THEN
            i = LEN(TRIM(iniFilePath))
            IF (iniFilePath(i-2:i) .EQ. "bin") THEN
               CALL INITFROMBIN(iniFilePath)
            ELSE
               CALL INITFROMVTU(iniFilePath)
            END IF
         ELSE
            IF (stFileFlag) THEN
               fTmp = TRIM(stFileName)//"_last.bin"
               INQUIRE (FILE=fTmp, EXIST=flag)
               IF (flag) THEN
                  CALL INITFROMBIN(fTmp)
               ELSE
                  IF (cm%mas()) wrn = TRIM(fTmp)//" can not be opened"
                  CALL ZEROINIT
               END IF
               IF (rmsh%isReqd) THEN
                  rmsh%fTS = (cTS/rmsh%fTS + 1)*rmsh%freq
                  rmsh%rTS = cTS
                  rmsh%time = time
                  rmsh%iNorm(:) = eq(:)%iNorm
                  rmsh%A0(:,:) = Ao(:,:)
                  rmsh%Y0(:,:) = Yo(:,:)
                  rmsh%D0(:,:) = Do(:,:)
               END IF
            ELSE
               CALL ZEROINIT
            END IF ! stFileFlag
         END IF

!        For CMM Equation
         IF (cmmEq) THEN
!           Reset time step
            cTS  = 0
            time = 0D0

!           Reset dof written to restart files
            i = 3
            IF (pstEq) i = i + nstd
            i = 4*(1+SIZE(stamp)) + 8*(2+nEq+cplBC%nX+i*tDof*tnNo)
            IF (cm%seq()) THEN
               recLn = i
            ELSE
               CALL MPI_ALLREDUCE(i, recLn, 1, mpint, MPI_MAX, cm%com(),
     2            ierr)
            END IF
         END IF

!        Load any explicitly provided solution variables
         CALL INITSOLVAR()

         rsTS = cTS
      ELSE
         cTS  = rmsh%rTS
         time = rmsh%time
         eq(:)%iNorm = rmsh%iNorm(:)
         Ao = LOCAL(rmsh%A0)
         Yo = LOCAL(rmsh%Y0)
         Do = LOCAL(rmsh%D0)
         DEALLOCATE(rmsh%A0,rmsh%Y0,rmsh%D0)
         ALLOCATE(rmsh%A0(tDof,tnNo)); rmsh%A0(:,:) = Ao(:,:)
         ALLOCATE(rmsh%Y0(tDof,tnNo)); rmsh%Y0(:,:) = Yo(:,:)
         ALLOCATE(rmsh%D0(tDof,tnNo)); rmsh%D0(:,:) = Do(:,:)
      END IF ! resetSim

      DO iM=1, nMsh
         IF (cm%mas()) THEN
            fTmp = TRIM(appPath)//".partitioning_"//
     2         TRIM(msh(iM)%name)//".bin"
            sTmp = TRIM(appPath)//".partitioning_"//
     2         TRIM(msh(iM)%name)//"_"//STR(cTS)//".bin"
            INQUIRE(FILE=TRIM(fTmp), EXIST=flag)
            IF (flag) THEN
               sTmp = "cp  "//TRIM(fTmp)//" "//TRIM(sTmp)
               CALL SYSTEM(TRIM(sTmp))
            END IF
         END IF
      END DO

!     Initialize Immersed Boundary data structures
      ALLOCATE(iblank(tnNo), ighost(tnNo))
      iblank = 0
      ighost = 0
      IF (ibFlag) THEN
         CALL IB_INIT(Do)
!        For all immersed shells, reset ibLSptr and nFacesLS
         i = SUM(iblank)
         i = cm%reduce(i)
         IF (i.EQ.0 .OR. ib%mthd.EQ.ibMthd_IFEM) THEN
            ibLSptr = 0
            nFacesLS = nFacesLS - 1

!        Reset FSILS structures as nFacesLS has been changed
            CALL FSILS_COMMU_FREE(communicator)
            CALL FSILS_LHS_FREE(lhs)
            CALL FSILS_COMMU_CREATE(communicator, cm%com())
            CALL FSILS_LHS_CREATE(lhs, communicator, gtnNo, tnNo, nnz,
     2         ltg, rowPtr, colPtr, nFacesLS)
         END IF
      END IF

!     Calculating the volume of each domain
      ALLOCATE(s(1,tnNo))
      s = 1D0
      DO iEq=1, nEq
         DO iDmn=1, eq(iEq)%nDmn
            eq(iEq)%dmn(iDmn)%v = Integ(eq(iEq)%dmn(iDmn)%Id, s, 1, 1)
            IF (cm%mas()) THEN
               WRITE(*,'(A,F12.6)') " Volume of domain ",
     2            eq(iEq)%dmn(iDmn)%v
            END IF
            IF (ISZERO(eq(iEq)%dmn(iDmn)%v)) wrn = "Volume of "//
     2         "domain "//iDmn//" of equation "//iEq//" is zero"
         END DO
      END DO

!     Predicting new variables
      CALL PICP

!     Preparing faces and BCs
      CALL BAFINI()

!     Making sure the old solution satisfies BCs
      CALL SETBCDIR(Ao, Yo, Do)

!     Preparing TXT files
      CALL TXT(.TRUE.)

!     Printing the first line and initializing timeP
      CALL OUTRESULT(timeP, 1, 1)
      rmsh%flag(:) = .FALSE.
      resetSim = .FALSE.

      RETURN
      CONTAINS
!--------------------------------------------------------------------
!     Initializing accelaration, velocity and displacement to zero
      SUBROUTINE ZEROINIT
      IMPLICIT NONE

      std = " Initializing state variables to zero"

!     This cTS corresponds to old variables. As soon as incrementing it
!     by one, it will be associated to new variables.
      cTS      = 0
      time     = 0D0
      timeP(1) = 0D0
      eq%iNorm = 0D0

      Ao = 0D0
      Yo = 0D0
      Do = 0D0
      IF (ALLOCATED(ADo)) ADo = 0D0

      IF (rmsh%isReqd) THEN
         rmsh%A0 = 0D0
         rmsh%Y0 = 0D0
         rmsh%D0 = 0D0
      END IF

      RETURN
      END SUBROUTINE ZEROINIT
!--------------------------------------------------------------------
!     Using the saved VTU files for initialization
      SUBROUTINE INITFROMVTU(fName)
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER l

      REAL(KIND=8), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:)

      l = LEN(TRIM(fName))
      IF (fName(l-2:l) .NE. "vtu") err = "Format of <"//
     2   TRIM(fName)//"> is not recognized"

      IF (ANY(msh%eType .EQ. eType_NRB)) err = "Only isoparametric"
     2   //" meshes can be read from VTU files"

      IF (nMsh .GT. 1) wrn = "Reading from VTU will fail in"//
     2   " presence of projected boundaries"
      std = " Initializing from "//fName

      cTS      = 0
      time     = 0D0
      timeP(1) = 0D0
      eq%iNorm = 0D0

      IF (cm%mas()) THEN
         ALLOCATE(tmpA(tDof,gtnNo), tmpY(tDof,gtnNo), tmpD(tDof,gtnNo))
         CALL READVTUS(tmpA, tmpY, tmpD, fName)
      ELSE
         ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
      END IF

      Ao = LOCAL(tmpA)
      Yo = LOCAL(tmpY)
      Do = LOCAL(tmpD)

      RETURN
      END SUBROUTINE INITFROMVTU
!--------------------------------------------------------------------
!     Using the svFSI specific format binary file for initialization
      SUBROUTINE INITFROMBIN(fName)
      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER, PARAMETER :: fid = 1

      INTEGER tStamp(SIZE(stamp))

      std = " Initializing from "//fName
      OPEN(fid, FILE=fName, ACCESS='DIRECT', RECL=recLn)
      IF (dFlag .AND. .NOT.cmmEq) THEN
         IF (pstEq) THEN
            READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm,
     2         cplBC%xo, Yo, Ao, Do, pS0
         ELSE
            READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm,
     2         cplBC%xo, Yo, Ao, Do
         END IF
      ELSE
         READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm,
     2      cplBC%xo, Yo, Ao
      END IF
      CLOSE(fid)

!     First checking all variables on master processor, since on the
!     other processor data will be shifted due to any change on the
!     sizes
      IF (cm%mas()) THEN
         IF (tStamp(1).NE.stamp(1)) err = "Number of processors <"//
     2      tStamp(1)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(1)//">"
         IF (tStamp(2).NE.stamp(2)) err = "Number of equations <"//
     2      tStamp(2)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(2)//">"
         IF (tStamp(3).NE.stamp(3)) err = "Number of meshes <"//
     2      tStamp(3)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(3)//">"
         IF (tStamp(4).NE.stamp(4)) err = "Number of nodes <"//
     2      tStamp(4)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(4)//">"
         IF (tStamp(5).NE.stamp(5)) err = "Number of cplBC%x <"//
     2      tStamp(5)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(5)//">"
         IF (tStamp(6).NE.stamp(6)) err = "Number of dof <"//
     2      tStamp(6)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(6)//">"
         IF (tStamp(7).NE.stamp(7)) err = "dFlag specification"//
     2      " <"//tStamp(7)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(7)//">"
         IF (tStamp(8).NE.stamp(8)) err = "Version of solver <"//
     2      tStamp(8)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(8)//">"
      END IF

      CALL MPI_BARRIER(cm%com(), ierr)
      IF (ANY(tStamp.NE.stamp)) err = "Simulation stamp"//
     2   " does not match with "//fName

      RETURN
      END SUBROUTINE INITFROMBIN
!--------------------------------------------------------------------
      SUBROUTINE INITSOLVAR()
      IMPLICIT NONE

      INTEGER a
      REAL(KIND=8), ALLOCATABLE :: tmpR(:,:)

!     Transfer initial velocity field
      flag = ALLOCATED(Vinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (tDof .LT. nsd) err = "Not enough dof to initialize "//
     2      "velocity field"
         IF (cm%mas()) THEN
            ALLOCATE(tmpR(nsd,gtnNo))
            tmpR = Vinit
            DEALLOCATE(Vinit)
         ELSE
            ALLOCATE(tmpR(0,0))
         END IF
         ALLOCATE(Vinit(nsd,tnNo))
         Vinit = LOCAL(tmpR)
         DEALLOCATE(tmpR)
         DO a=1, tnNo
            Yo(1:nsd,a) = Vinit(:,a)
         END DO
         DEALLOCATE(Vinit)
      END IF

!     Transfer initial pressure field
      flag = ALLOCATED(Pinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (tDof .LT. nsd+1) err = "Not enough dof to initialize "//
     2      "pressure field"
         IF (cm%mas()) THEN
            ALLOCATE(tmpR(1,gtnNo))
            tmpR(1,:) = Pinit
            DEALLOCATE(Pinit)
         ELSE
            ALLOCATE(tmpR(0,0))
         END IF
         ALLOCATE(Pinit(tnNo))
         tmpR = LOCAL(tmpR)
         Pinit(:) = tmpR(1,:)
         DEALLOCATE(tmpR)
         DO a=1, tnNo
            Yo(nsd+1,a) = Pinit(a)
         END DO
         DEALLOCATE(Pinit)
      END IF

!     Transfer initial displacement field
      flag = ALLOCATED(Dinit)
      CALL cm%bcast(flag)
      IF (flag) THEN
         IF (.NOT.dFlag) err = "Incorrect phys to initialize "//
     2      "displacement field"
         IF (cm%mas()) THEN
            ALLOCATE(tmpR(nsd,gtnNo))
            tmpR = Dinit
            DEALLOCATE(Dinit)
         ELSE
            ALLOCATE(tmpR(0,0))
         END IF
         ALLOCATE(Dinit(nsd,tnNo))
         Dinit = LOCAL(tmpR)
         DEALLOCATE(tmpR)
         DO a=1, tnNo
            Do(1:nsd,a) = Dinit(:,a)
         END DO
         DEALLOCATE(Dinit)
      END IF

      RETURN
      END SUBROUTINE INITSOLVAR
!--------------------------------------------------------------------
      END SUBROUTINE INITIALIZE
!####################################################################
      SUBROUTINE FINALIZE
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER iM, iEq

!     Deallocating meshes
      IF (ALLOCATED(msh)) THEN
         DO iM=1, nMsh
            CALL DESTROY(msh(iM))
         END DO
         DEALLOCATE(msh)
      END IF

!     Deallocating equations
      IF (ALLOCATED(eq)) THEN
         DO iEq=1, nEq
            CALL DESTROY(eq(iEq))
         END DO
         DEALLOCATE(eq)
      END IF

!     Deallocating sparse matrix structures
      IF(lhs%foc) CALL FSILS_LHS_FREE(lhs)
#ifdef WITH_TRILINOS
      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
         CALL TRILINOS_LHS_FREE()
      END IF
#endif
      IF (.NOT.useTrilinosAssemAndLS) THEN
         IF (ALLOCATED(Val))   DEALLOCATE(Val)
      END IF

      IF (ALLOCATED(colPtr))   DEALLOCATE(colPtr)
      IF (ALLOCATED(dmnId))    DEALLOCATE(dmnId)
      IF (ALLOCATED(ltg))      DEALLOCATE(ltg)
      IF (ALLOCATED(rowPtr))   DEALLOCATE(rowPtr)
      IF (ALLOCATED(iblank))   DEALLOCATE(iblank)
      IF (ALLOCATED(ighost))   DEALLOCATE(ighost)

      IF (ALLOCATED(Ao))       DEALLOCATE(Ao)
      IF (ALLOCATED(An))       DEALLOCATE(An)
      IF (ALLOCATED(Yo))       DEALLOCATE(Yo)
      IF (ALLOCATED(Yn))       DEALLOCATE(Yn)
      IF (ALLOCATED(Do))       DEALLOCATE(Do)
      IF (ALLOCATED(Dn))       DEALLOCATE(Dn)
      IF (ALLOCATED(R))        DEALLOCATE(R)
      IF (ALLOCATED(fN))       DEALLOCATE(fN)

      IF (ALLOCATED(ADo))      DEALLOCATE(ADo)
      IF (ALLOCATED(ADn))      DEALLOCATE(ADn)
      IF (ALLOCATED(Rd))       DEALLOCATE(Rd)
      IF (ALLOCATED(Kd))       DEALLOCATE(Kd)

      IF (ALLOCATED(pS0))      DEALLOCATE(pS0)
      IF (ALLOCATED(pSn))      DEALLOCATE(pSn)
      IF (ALLOCATED(pSa))      DEALLOCATE(pSa)

      IF (ALLOCATED(Vinit))    DEALLOCATE(Vinit)
      IF (ALLOCATED(Pinit))    DEALLOCATE(Pinit)
      IF (ALLOCATED(Dinit))    DEALLOCATE(Dinit)
      IF (ALLOCATED(Rib))      DEALLOCATE(Rib)

      IF (ALLOCATED(cplBC%fa)) DEALLOCATE(cplBC%fa)
      IF (ALLOCATED(cplBC%xn)) DEALLOCATE(cplBC%xn)
      IF (ALLOCATED(cplBC%xo)) DEALLOCATE(cplBC%xo)

!     IB structures
      IF (ibFlag) THEN
         IF (ALLOCATED(ib%dmnId))  DEALLOCATE(ib%dmnId)
         IF (ALLOCATED(ib%rowPtr)) DEALLOCATE(ib%rowPtr)
         IF (ALLOCATED(ib%colPtr)) DEALLOCATE(ib%colPtr)
         IF (ALLOCATED(ib%x))      DEALLOCATE(ib%x)
         IF (ALLOCATED(ib%fN))     DEALLOCATE(ib%fN)
         IF (ALLOCATED(ib%An))     DEALLOCATE(ib%An)
         IF (ALLOCATED(ib%Ao))     DEALLOCATE(ib%Ao)
         IF (ALLOCATED(ib%Yn))     DEALLOCATE(ib%Yn)
         IF (ALLOCATED(ib%Yo))     DEALLOCATE(ib%Yo)
         IF (ALLOCATED(ib%Un))     DEALLOCATE(ib%Un)
         IF (ALLOCATED(ib%Uo))     DEALLOCATE(ib%Uo)
         IF (ALLOCATED(ib%R))      DEALLOCATE(ib%R)

         DO iM=1, ib%nMsh
            CALL DESTROY(ib%msh(iM))
         END DO
         DEALLOCATE(ib%msh)

         IF (ib%lhs%foc) CALL FSILS_LHS_FREE(ib%lhs)

         DEALLOCATE(ib)
      END IF

!     Closing the output channels
      CALL std%close()
      CALL wrn%close()
      CALL err%close()
      CALL dbg%close()

      RETURN
      END SUBROUTINE FINALIZE
!####################################################################

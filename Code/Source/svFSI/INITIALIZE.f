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
      REAL(KIND=RKIND), INTENT(OUT) :: timeP(3)

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: i, a, iEq, iDmn, iM, iFa, ierr, nnz, gnnz
      CHARACTER(LEN=stdL) :: fTmp, sTmp
      REAL(KIND=RKIND) :: am
      TYPE(FSILS_commuType) :: communicator

      REAL(KIND=RKIND), ALLOCATABLE, DIMENSION(:,:) :: s

      tDof  = 0
      dFlag = .FALSE.

!     Set faces for linear solver
      nFacesLS = SUM(eq%nBc)

!     Remove LS pointer for faces with weakly applied Dir. BC
      DO iEq=1, nEq
         DO i=1, eq(iEq)%nBc
            IF (eq(iEq)%bc(i)%weakDir) nFacesLS = nFacesLS - 1
         END DO
      END DO

!     For FSI simulations, LS pointer to structure nodes
      IF (mvMsh) nFacesLS = nFacesLS + 1

!     For initializing CMM, LS pointer to fixed edge nodes
      IF (cmmInit) nFacesLS = nFacesLS + 1

      IF (ANY(eq(1)%bc%cplBCptr .NE. 0)) cplBC%coupled = .TRUE.

      flag = .FALSE.
      DO iEq=1, nEq
!     This would be the default value of am for first order equations
         eq(iEq)%am = 0.5_RKIND*(3._RKIND - eq(iEq)%roInf) /
     2                (1._RKIND + eq(iEq)%roInf)
!     If the equation is second order, am is calculated as follows
         am = (2._RKIND - eq(iEq)%roInf)/(1._RKIND + eq(iEq)%roInf)
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
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'ST'
         CASE (phys_CMM)
            dFlag = .TRUE.
            eq(iEq)%dof = nsd + 1
            IF (cmmInit) eq(iEq)%dof = nsd
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
         CASE (phys_CEP)
            eq(iEq)%dof = 1
            eq(iEq)%sym = 'EP'
         CASE (phys_stokes)
            eq(iEq)%dof = nsd + 1
            eq(iEq)%sym = 'SS'
         CASE DEFAULT
            err = "Equation type "//eq(iEq)%sym//" is not defined"
         END SELECT

         eq(iEq)%pNorm = HUGE(eq(iEq)%pNorm)
         eq(iEq)%af    = 1._RKIND/(1._RKIND + eq(iEq)%roInf)
         eq(iEq)%beta  = 0.25_RKIND*(1._RKIND + eq(iEq)%am -
     2      eq(iEq)%af)**2._RKIND
         eq(iEq)%gam   = 0.5_RKIND + eq(iEq)%am - eq(iEq)%af
         eq(iEq)%s     = tDof + 1
         eq(iEq)%e     = tDof + eq(iEq)%dof
         tDof          = eq(iEq)%e
         IF (eq(iEq)%useTLS) flag = .TRUE.
      END DO

      ierr = 0
      IF (dFlag) ierr = 1

      i = 0
      IF (cplBC%coupled) i = cplBC%nX

      stamp = (/cm%np(), nEq, nMsh, tnNo, i, tDof, ierr/)

!     Calculating the record length
      i = 2*tDof
      IF (dFlag) i = 3*tDof
      IF (pstEq) i = i + nsymd
      IF (sstEq) i = i + nsd
!      IF (useVarWall) i = i + nvwp
      IF (cepEq) THEN
         i = i + nXion
         IF (cem%cpld) i = i + 1
      END IF
      i = IKIND*(1+SIZE(stamp)) + RKIND*(2+nEq+cplBC%nX+i*tnNo)

      IF (ibFlag) i = i + RKIND*(3*nsd+1)*ib%tnNo
      IF (cm%seq()) THEN
         recLn = i
      ELSE
         CALL MPI_ALLREDUCE(i, recLn, 1, mpint, MPI_MAX, cm%com(), ierr)
      END IF

!     Initialize shell eIEN data structure. Used later in LHSA.
      IF (shlEq) THEN
         DO iM=1, nMsh
            IF (msh(iM)%lShl) THEN
               IF (msh(iM)%eType .EQ. eType_NRB) THEN
                  ALLOCATE(msh(iM)%eIEN(0,0))
                  ALLOCATE(msh(iM)%sbc(msh(iM)%eNoN,msh(iM)%nEl))
                  msh(iM)%sbc = 0
               ELSE IF (msh(iM)%eType .EQ. eType_TRI3) THEN
                  CALL SETSHLXIEN(msh(iM))
               END IF
            END IF
         END DO
      END IF

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

!     Initialize Trilinos data structure
      IF (flag) THEN
         ALLOCATE(tls)
         ALLOCATE(tls%ltg(tnNo))
         tls%ltg = 0
         DO a=1, tnNo
           tls%ltg(lhs%map(a)) = ltg(a)
         END DO
      END IF

!     Variable allocation and initialization
      ALLOCATE(Ao(tDof,tnNo), An(tDof,tnNo), Yo(tDof,tnNo),
     2   Yn(tDof,tnNo), Do(tDof,tnNo), Dn(tDof,tnNo), Bf(nsd,tnNo))
      IF (ibFlag) CALL IB_MEMALLOC()

!     Additional physics dependent variables
!     USTRUCT phys
      IF (sstEq) THEN
         ALLOCATE(Ad(nsd,tnNo), Rd(nsd,tnNo), Kd((nsd+1)*nsd,nnz))
         Ad = 0._RKIND
         Rd = 0._RKIND
         Kd = 0._RKIND
      END IF

!     PRESTRESS
      IF (pstEq) THEN
         IF (ALLOCATED(pS0)) err = "Prestress already allocated. "//
     2      "Correction needed"
         ALLOCATE(pS0(nsymd,tnNo), pSn(nsymd,tnNo), pSa(tnNo))
         pS0 = 0._RKIND
         pSn = 0._RKIND
         pSa = 0._RKIND
      END IF

!     Electrophysiology
      IF (cepEq) THEN
         ALLOCATE(Xion(nXion,tnNo))
         Xion(:,:) = 0._RKIND

         CALL CEPINIT()

!        Electro-Mechanics
         IF (cem%cpld) THEN
            ALLOCATE(cem%Ya(tnNo))
            cem%Ya = 0._RKIND
         END IF
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
               CALL INITFROMBIN(iniFilePath, timeP)
            ELSE
               CALL INITFROMVTU(iniFilePath, timeP)
            END IF
         ELSE
            IF (stFileFlag) THEN
               fTmp = TRIM(stFileName)//"_last.bin"
               INQUIRE (FILE=fTmp, EXIST=flag)
               IF (flag) THEN
                  CALL INITFROMBIN(fTmp, timeP)
               ELSE
                  IF (cm%mas()) THEN
                     wrn = TRIM(fTmp)//" can not be opened. "//
     2                  "Resetting restart flag to false"
                  END IF
                  stFileFlag = .FALSE.
                  CALL ZEROINIT(timeP)
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
               CALL ZEROINIT(timeP)
            END IF ! stFileFlag
         END IF

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

!     Initialize new variables
      An = Ao
      Yn = Yo
      Dn = Do

      DO iM=1, nMsh
         IF (cm%mas()) THEN
            fTmp = TRIM(appPath)//".partitioning_"//TRIM(msh(iM)%name)//
     2         ".bin"
            sTmp = TRIM(appPath)//".partitioning_"//TRIM(msh(iM)%name)//
     2         "_"//STR(cTS)//".bin"
            INQUIRE(FILE=TRIM(fTmp), EXIST=flag)
            IF (flag) THEN
               sTmp = "cp  "//TRIM(fTmp)//" "//TRIM(sTmp)
               CALL SYSTEM(TRIM(sTmp))
            END IF
         END IF
      END DO

!     Initialize function spaces
      DO iM=1, nMsh
         CALL INITFSMSH(msh(iM))
         DO iFa=1, msh(iM)%nFa
            CALL INITFSFACE(msh(iM), msh(iM)%fa(iFa))
         END DO
      END DO

!     Initialize Immersed Boundary data structures
      ALLOCATE(iblank(tnNo))
      iblank = 0
      IF (ibFlag) CALL IB_INIT(Do)

!     Calculating the volume of each domain
      ALLOCATE(s(1,tnNo))
      s = 1._RKIND
      DO iEq=1, nEq
         IF (.NOT.shlEq .AND. .NOT.cmmInit) std = " Eq. <"//
     2      CLR(eq(iEq)%sym, iEq)//">"
         DO iDmn=1, eq(iEq)%nDmn
            i = eq(iEq)%dmn(iDmn)%Id
            eq(iEq)%dmn(iDmn)%v = Integ(i, s, 1, 1)
            IF (.NOT.shlEq .AND. .NOT.cmmInit) THEN
               std = "    Volume of domain <"//STR(i)//"> is "//
     2            STR(eq(iEq)%dmn(iDmn)%v)
               IF (ISZERO(eq(iEq)%dmn(iDmn)%v)) wrn = "<< Volume of "//
     2            "domain "//iDmn//" of equation "//iEq//" is zero >>"
            END IF
         END DO
      END DO

!     Preparing faces and BCs
      CALL BAFINI()

!     As all the arrays are allocated, call BIN to VTK for conversion
      IF (bin2VTK) CALL PPBIN2VTK()

!     Making sure the old solution satisfies BCs
      CALL SETBCDIR(Ao, Yo, Do)

!     Preparing TXT files
      CALL TXT(.TRUE.)

!     Printing the first line and initializing timeP
      CALL OUTRESULT(timeP, 1, 1)
      rmsh%flag(:) = .FALSE.
      resetSim = .FALSE.

      RETURN
      END SUBROUTINE INITIALIZE
!--------------------------------------------------------------------
!     Initializing accelaration, velocity and displacement to zero
      SUBROUTINE ZEROINIT(timeP)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(OUT) :: timeP(3)

      INTEGER(KIND=IKIND) a

      std = " Initializing state variables to zero"

!     This cTS corresponds to old variables. As soon as incrementing it
!     by one, it will be associated to new variables.
      cTS      = startTS
      time     = 0._RKIND
      timeP(1) = 0._RKIND
      eq%iNorm = 0._RKIND

      Ao = 0._RKIND
      Yo = 0._RKIND
      Do = 0._RKIND

      IF (rmsh%isReqd) THEN
         rmsh%A0 = 0._RKIND
         rmsh%Y0 = 0._RKIND
         rmsh%D0 = 0._RKIND
      END IF

      IF (ibFlag) THEN
         ib%Yb  = 0._RKIND
         ib%Auo = 0._RKIND
         ib%Ubo = 0._RKIND

         IF (ib%cpld .EQ. ibCpld_I) THEN
            ib%Aun = 0._RKIND
            ib%Auk = 0._RKIND
            ib%Ubn = 0._RKIND
            ib%Ubk = 0._RKIND
            ib%Uo  = 0._RKIND
            ib%Un  = 0._RKIND
         END IF
      END IF

!     Load any explicitly provided solution variables
      IF (ALLOCATED(Vinit)) THEN
         DO a=1, tnNo
            Yo(1:nsd,a) = Vinit(:,a)
         END DO
      END IF

      IF (ALLOCATED(Pinit)) THEN
         DO a=1, tnNo
            Yo(nsd+1,a) = Pinit(a)
         END DO
      END IF

      IF (ALLOCATED(Dinit)) THEN
         DO a=1, tnNo
            Do(1:nsd,a) = Dinit(:,a)
         END DO
      END IF

      RETURN
      END SUBROUTINE ZEROINIT
!--------------------------------------------------------------------
!     Using the saved VTU files for initialization
      SUBROUTINE INITFROMVTU(fName, timeP)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      REAL(KIND=RKIND), INTENT(OUT) :: timeP(3)

      INTEGER(KIND=IKIND) l
      REAL(KIND=RKIND), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:)

      l = LEN(TRIM(fName))
      IF (fName(l-2:l) .NE. "vtu") err = "Format of <"//
     2   TRIM(fName)//"> is not recognized"

      IF (ANY(msh%eType .EQ. eType_NRB)) err = "Only isoparametric"
     2   //" meshes can be read from VTU files"

      IF (nMsh .GT. 1) wrn = "Reading from VTU will fail in"//
     2   " presence of projected boundaries"
      std = " Initializing from "//fName

      cTS      = 0
      time     = 0._RKIND
      timeP(1) = 0._RKIND
      eq%iNorm = 0._RKIND

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
      SUBROUTINE INITFROMBIN(fName, timeP)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      REAL(KIND=RKIND), INTENT(OUT) :: timeP(3)

      INTEGER(KIND=IKIND), PARAMETER :: fid = 1
      INTEGER(KIND=IKIND) tStamp(SIZE(stamp)), i

      i = 0
      IF (.NOT.bin2VTK) THEN
         std = " Initializing from "//fName
      END IF

      OPEN(fid, FILE=fName, ACCESS='DIRECT', RECL=recLn)
      IF (.NOT.ibFlag) THEN
         IF (dFlag) THEN
            IF (sstEq) THEN
               IF (pstEq) THEN
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do, pS0, Ad
               ELSE IF (cepEq) THEN
                  IF (.NOT.cem%cpld) err = "Incorrect equation "//
     2               "combination. Cannot load restart files"
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do, Ad, Xion, cem%Ya
               ELSE
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do, Ad
               END IF
            ELSE
               IF (pstEq) THEN
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do, pS0
               ELSE IF (cepEq) THEN
                  IF (.NOT.cem%cpld) err = "Incorrect equation "//
     2               "combination. Cannot load restart files"
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do, Xion, cem%Ya
               ELSE
                  READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2               eq%iNorm, cplBC%xo, Yo, Ao, Do
               END IF
            END IF
         ELSE
            IF (cepEq) THEN
               READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2          eq%iNorm, cplBC%xo, Yo, Ao, Xion
            ELSE
               READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2            eq%iNorm, cplBC%xo, Yo, Ao
            END IF
         END IF
      ELSE
         IF (dFlag) THEN
            IF (pstEq) THEN
               READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2            eq%iNorm, cplBC%xo, Yo, Ao, Do, pS0, ib%Yb, ib%Auo,
     3            ib%Ubo
            ELSE
               READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1),
     2            eq%iNorm, cplBC%xo, Yo, Ao, Do, ib%Yb, ib%Auo, ib%Ubo
            END IF
         ELSE
            READ(fid,REC=cm%tF()) tStamp, cTS, time, timeP(1), eq%iNorm,
     2         cplBC%xo, Yo, Ao, ib%Yb, ib%Auo, ib%Ubo
         END IF
         IF (ib%cpld .EQ. ibCpld_I) THEN
            ib%Aun = ib%Auo
            ib%Ubn = ib%Ubo
         END IF
      END IF
      CLOSE(fid)

!     First checking all variables on master processor, since on the
!     other processor data will be shifted due to any change on the
!     sizes
      IF (cm%mas()) THEN
         IF (tStamp(1) .NE. stamp(1)) err = "Number of processors <"//
     2      tStamp(1)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(1)//">"
         IF (tStamp(2) .NE. stamp(2)) err = "Number of equations <"//
     2      tStamp(2)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(2)//">"
         IF (tStamp(3) .NE. stamp(3)) err = "Number of meshes <"//
     2      tStamp(3)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(3)//">"
         IF (tStamp(4) .NE. stamp(4)) err = "Number of nodes <"//
     2      tStamp(4)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(4)//">"
         IF (tStamp(5) .NE. stamp(5)) err = "Number of cplBC%x <"//
     2      tStamp(5)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(5)//">"
         IF (tStamp(6) .NE. stamp(6)) err = "Number of dof <"//
     2      tStamp(6)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(6)//">"
         IF (tStamp(7) .NE. stamp(7)) err = "dFlag specification"//
     2      " <"//tStamp(7)//"> does not match with "//
     3      TRIM(fName)//" <"//stamp(7)//">"
      END IF

      CALL cm%bcast(i)
      IF (ANY(tStamp .NE. stamp)) err = "Simulation stamp"//
     2   " does not match with "//fName

      RETURN
      END SUBROUTINE INITFROMBIN
!####################################################################
      SUBROUTINE FINALIZE
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iM, iEq

!     Deallocating meshes
      IF (ALLOCATED(msh)) THEN
         DO iM=1, nMsh
            CALL DESTROY(msh(iM))
         END DO
         DEALLOCATE(msh)
      END IF
      IF (ALLOCATED(rmsh%flag)) DEALLOCATE(rmsh%flag)

!     Deallocating equations
      IF (ALLOCATED(eq)) THEN
         DO iEq=1, nEq
            CALL DESTROY(eq(iEq))
         END DO
         DEALLOCATE(eq)
      END IF

!     Deallocating sparse matrix structures
      IF(lhs%foc) CALL FSILS_LHS_FREE(lhs)
      IF (ALLOCATED(tls)) THEN
         IF (ALLOCATED(tls%W))   DEALLOCATE(tls%W)
         IF (ALLOCATED(tls%R))   DEALLOCATE(tls%R)
         IF (ALLOCATED(tls%ltg)) DEALLOCATE(tls%ltg)
         DEALLOCATE(tls)
#ifdef WITH_TRILINOS
         CALL TRILINOS_LHS_FREE()
#endif
      ELSE
         IF (ALLOCATED(Val))   DEALLOCATE(Val)
      END IF

      IF (ALLOCATED(colPtr))   DEALLOCATE(colPtr)
      IF (ALLOCATED(dmnId))    DEALLOCATE(dmnId)
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
      IF (ALLOCATED(x))        DEALLOCATE(x)
      IF (ALLOCATED(Yo))       DEALLOCATE(Yo)
      IF (ALLOCATED(Yn))       DEALLOCATE(Yn)
      IF (ALLOCATED(Bf))       DEALLOCATE(Bf)

      IF (ALLOCATED(Ad))       DEALLOCATE(Ad)
      IF (ALLOCATED(Rd))       DEALLOCATE(Rd)
      IF (ALLOCATED(Kd))       DEALLOCATE(Kd)

      IF (ALLOCATED(pS0))      DEALLOCATE(pS0)
      IF (ALLOCATED(pSn))      DEALLOCATE(pSn)
      IF (ALLOCATED(pSa))      DEALLOCATE(pSa)

!     Varwall properties -----------------------------------------------
      IF (ALLOCATED(vWP0))      DEALLOCATE(vWP0)
!     ------------------------------------------------------------------

      IF (ALLOCATED(Pinit))    DEALLOCATE(Pinit)
      IF (ALLOCATED(Vinit))    DEALLOCATE(Vinit)
      IF (ALLOCATED(Dinit))    DEALLOCATE(Dinit)

      IF (ALLOCATED(cplBC%fa)) DEALLOCATE(cplBC%fa)
      IF (ALLOCATED(cplBC%xn)) DEALLOCATE(cplBC%xn)
      IF (ALLOCATED(cplBC%xo)) DEALLOCATE(cplBC%xo)
      IF (ALLOCATED(cplBC%xp)) DEALLOCATE(cplBC%xp)

      IF (ALLOCATED(varWallProps)) DEALLOCATE(varWallProps)

!     Electrophysiology and Electromechanics
      IF (cepEq) THEN
         IF (ALLOCATED(Xion))  DEALLOCATE(Xion)
         IF (cem%cpld) THEN
            IF (ALLOCATED(cem%Ya))  DEALLOCATE(cem%Ya)
         END IF
      END IF

!     IB structures
      IF (ibFlag) THEN
         IF (ALLOCATED(ib%dmnId))  DEALLOCATE(ib%dmnId)
         IF (ALLOCATED(ib%rowPtr)) DEALLOCATE(ib%rowPtr)
         IF (ALLOCATED(ib%colPtr)) DEALLOCATE(ib%colPtr)
         IF (ALLOCATED(ib%x))      DEALLOCATE(ib%x)
         IF (ALLOCATED(ib%Yb))     DEALLOCATE(ib%Yb)
         IF (ALLOCATED(ib%Auo))    DEALLOCATE(ib%Auo)
         IF (ALLOCATED(ib%Aun))    DEALLOCATE(ib%Aun)
         IF (ALLOCATED(ib%Auk))    DEALLOCATE(ib%Auk)
         IF (ALLOCATED(ib%Ubo))    DEALLOCATE(ib%Ubo)
         IF (ALLOCATED(ib%Ubn))    DEALLOCATE(ib%Ubn)
         IF (ALLOCATED(ib%Ubk))    DEALLOCATE(ib%Ubk)
         IF (ALLOCATED(ib%Uo))     DEALLOCATE(ib%Uo)
         IF (ALLOCATED(ib%Un))     DEALLOCATE(ib%Un)
         IF (ALLOCATED(ib%R))      DEALLOCATE(ib%R)
         IF (ALLOCATED(ib%Ru))     DEALLOCATE(ib%Ru)
         IF (ALLOCATED(ib%Rub))    DEALLOCATE(ib%Rub)
         IF (ALLOCATED(ib%Ku))     DEALLOCATE(ib%Ku)
         IF (ALLOCATED(ib%cm%n))   DEALLOCATE(ib%cm%n)
         IF (ALLOCATED(ib%cm%gN))  DEALLOCATE(ib%cm%gN)

         DO iM=1, ib%nMsh
            CALL DESTROY(ib%msh(iM))
         END DO
         DEALLOCATE(ib%msh)

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

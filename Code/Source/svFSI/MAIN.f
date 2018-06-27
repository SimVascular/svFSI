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
!     Main routine that contains the calls to all major routines and
!     general structure of the code.
!
!--------------------------------------------------------------------

      PROGRAM MAIN

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      LOGICAL l1, l2, l3, l4, l5
      INTEGER i, a, Ac, e, ierr, iEqOld, iBc, eNoN, iM, j, fid
c      INTEGER OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
      REAL(KIND=8) timeP(3)
      CHARACTER(LEN=stdL) fName, tmpS

      LOGICAL, ALLOCATABLE :: isS(:)
      INTEGER, ALLOCATABLE :: ptr(:), incL(:), ltgReordered(:)
      REAL(KIND=8), ALLOCATABLE :: xl(:,:), Ag(:,:), al(:,:), Yg(:,:),
     2   yl(:,:), Dg(:,:), dl(:,:), dol(:,:), fNl(:,:), res(:),
     3   RTrilinos(:,:), dirW(:,:)

!--------------------------------------------------------------------
      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      l4 = .FALSE.
      l5 = .FALSE.

      fid       = 27
      savedOnce = .FALSE.
      CALL MPI_INIT(i)
      CALL cm%new(MPI_COMM_WORLD)

!     Initiating the exception tracing
      CALL EXCEPTIONS

      resetSim  = .FALSE.
      rmsh%cntr = 0

!     Reading the user-defined parameters from foo.inp
 101  CALL READFILES

!     Doing the partitioning and distributing the data to the all
!     Processors
      CALL DISTRIBUTE

!     Initializing the solution vectors and constructing LHS matrix
!     format
      CALL INITIALIZE(timeP)

!     only compute once
      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
         ALLOCATE(ltgReordered(tnNo))
         DO i = 1, tnNo
           ltgReordered(lhs%map(i)) = ltg(i)
         END DO
      END IF

      dbg = 'Allocating intermediate variables'
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo),
     3   res(nFacesLS), incL(nFacesLS), isS(tnNo))
      isS = .FALSE.
      DO Ac=1, tnNo
         IF (ISDOMAIN(1, Ac, phys_struct)) isS(Ac) = .TRUE.
      END DO

#ifdef WITH_TRILINOS
#else
      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
           write(*,*) "svFSI is not compiled with Trilinos:
     2         currently FSILS preconditioning will be used instead.
     3         Use FSILS in the input file or set it default."
      ENDIF
#endif
!--------------------------------------------------------------------
!     Outer loop for marching in time. When entring this loop, all old
!     variables are completely set and satisfy BCs.
      IF (cTS .LE. nITS) dt = dt/1D1
      DO
!     Adjusting the time step size once initialization stage is over
         IF (cTS .EQ. nITS) THEN
            dt = dt*1D1
            std = "New time step size: "//dt
         END IF
!     Incrementing time step, hence cTS will be associated with new
!     variables, i.e. An, Yn, and Dn
         cTS    = cTS + 1
         time   = time + dt
         cEq    = 1
         eq%itr = 0
         eq%ok  = .FALSE.

!     Compute mesh properties to check if remeshing is required
         IF (mvMsh .AND. rmsh%isReqd) THEN
            CALL CALCMESHPROPS(nMsh, msh)
            IF (resetSim) EXIT
         END IF

!     Predictor step
         CALL PICP
         CALL SETBCDIR(An, Yn, Dn)

!     Inner loop for iteration
         DO
            iEqOld = cEq

            IF (cplBC%coupled .AND. cEq.EQ.1) THEN
               CALL SETBCCPL
               CALL SETBCDIR(An, Yn, Dn)
            END IF
!     Initiator step
            CALL PICI(Ag, Yg, Dg)

            dbg = 'Allocating the RHS and LHS'
            IF (ALLOCATED(R)) THEN
               IF (SIZE(R,1) .NE. dof) THEN
                  DEALLOCATE(R)
                  ALLOCATE (R(dof,tnNo))
                  IF (.NOT. useTrilinosAssemAndLS) THEN
                     DEALLOCATE(Val)
                     ALLOCATE (Val(dof*dof,lhs%nnz))
                  END IF

#ifdef WITH_TRILINOS
                  IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
                     DEALLOCATE(dirW, RTrilinos)
                     ALLOCATE(dirW(dof,tnNo), RTrilinos(dof,tnNo))
                     CALL TRILINOS_LHS_FREE() !free K and R in C++
                     CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo,
     2                  lhs%nnz, ltgReordered, ltg, rowPtr, colPtr, dof)
                  END IF
#endif
               END IF
            ELSE
               ALLOCATE (R(dof,tnNo))
               IF (.NOT. useTrilinosAssemAndLS) THEN
                  ALLOCATE(Val(dof*dof,lhs%nnz))
               END IF
#ifdef WITH_TRILINOS
               IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
                  CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo,
     2               lhs%nnz, ltgReordered, ltg, rowPtr, colPtr, dof)
                  ALLOCATE(dirW(dof, tnNo), RTrilinos(dof,tnNo))
               END IF
#endif
            END IF

            IF (.NOT. useTrilinosAssemAndLS) THEN
!$OMP DO SCHEDULE(GUIDED,mpBs)
               DO i=1, lhs%nnz
                  Val(:,i) = 0D0
               END DO
!$OMP END DO
            END IF
!$OMP DO SCHEDULE(GUIDED,mpBs)
            DO a=1, tnNo
               R(:,a) = 0D0
            END DO
!$OMP END DO
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, e, a, Ac, al, yl, dl,
!$OMP&   xl, fNl, ptr, cDmn)
            dbg = "Assembling equation <"//eq(cEq)%sym//">"
            DO iM=1, nMsh
               eNoN = msh(iM)%eNoN
               i    = nsd + 2
               j    = 2*nsd + 1
               ALLOCATE(al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
     2            dol(nsd,eNoN), xl(nsd,eNoN), fNl(nsd,eNoN),
     3            ptr(eNoN))
!$OMP DO SCHEDULE(GUIDED,mpBs)
               DO e=1, msh(iM)%nEl
                  fNl = 0D0
                  DO a=1, eNoN
                     Ac      = msh(iM)%IEN(a,e)
                     ptr(a)  = Ac
                     xl(:,a) = x (:,Ac)
                     al(:,a) = Ag(:,Ac)
                     yl(:,a) = Yg(:,Ac)
                     dl(:,a) = Dg(:,Ac)
                     IF (mvMsh) THEN
                        dol(:,a) = Do(i:j,Ac)
                     END IF
                     IF (ALLOCATED(fN)) fNl(:,a) = fN(:,Ac)
                  END DO
!     Add contribution of current equation to the LHS/RHS
                  CALL CONSTRUCT(msh(iM), al, yl, dl, dol, xl, fNl,
     2               ptr, e)
               END DO
!$OMP END DO
               DEALLOCATE(al, yl, dl, xl, dol, fNl, ptr)
               dbg = "Mesh "//iM//" is assembled"
            END DO

!     Constructing the element stiffness matrix for boundaries
            CALL SETBCNEU(Yg, Dg)
            incL = 0
            IF (eq(cEq)%phys .EQ. phys_mesh) incL(nFacesLS) = 1
            DO iBc=1, eq(cEq)%nBc
               i = eq(cEq)%bc(iBc)%lsPtr
               IF (i .NE. 0) THEN
                  res(i) = eq(cEq)%gam*dt*eq(cEq)%bc(iBc)%r
                  incL(i) = 1
               END IF
            END DO
!$OMP END PARALLEL

            dbg = "Solving equation <"//eq(cEq)%sym//">"
!     Initialize Dirichlet and coupled Neumann resistnace BC
!     for Trilinos
#ifdef WITH_TRILINOS
            IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
               CALL INIT_DIR_AND_COUPNEU_BC(incL, res, dirW)
            END IF
!     Solving the linear system of equations
            IF (useTrilinosAssemAndLS) THEN
               eq(cEq)%FSILS%RI%suc = .FALSE.
!     Calls C++ solver and stores solution vector in RTrilinos
               CALL TRILINOS_SOLVE(RTrilinos, dirW,
     2            eq(cEq)%FSILS%RI%fNorm, eq(cEq)%FSILS%RI%iNorm,
     3            eq(cEq)%FSILS%RI%itr, eq(cEq)%FSILS%RI%callD,
     4            eq(cEq)%FSILS%RI%dB, eq(cEq)%FSILS%RI%suc,
     5            eq(cEq)%ls%LS_type, eq(cEq)%FSILS%RI%reltol,
     6            eq(cEq)%FSILS%RI%mItr, eq(cEq)%FSILS%RI%sD,
     7            eq(cEq)%ls%PREC_Type)
            ELSE IF(useTrilinosLS) THEN
               CALL TRILINOS_GLOBAL_SOLVE(Val, R, RTrilinos,
     2            dirW, eq(cEq)%FSILS%RI%fNorm, eq(cEq)%FSILS%RI%iNorm,
     3            eq(cEq)%FSILS%RI%itr, eq(cEq)%FSILS%RI%callD,
     4            eq(cEq)%FSILS%RI%dB, eq(cEq)%FSILS%RI%suc,
     5            eq(cEq)%ls%LS_type, eq(cEq)%FSILS%RI%reltol,
     6            eq(cEq)%FSILS%RI%mItr, eq(cEq)%FSILS%RI%sD,
     7            eq(cEq)%ls%PREC_Type)
            ELSE
#endif
               CALL FSILS_SOLVE(lhs, eq(cEq)%FSILS, dof, R, Val, isS,
     2                          incL, res)
#ifdef WITH_TRILINOS
            END IF
!     For trilinos case need to do the reordering
            IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
               DO i=1, tnNo
!                 Rtrilinos is reordered so use map to convert back
                  R(:,i) = RTrilinos(:,lhs%map(i))
               END DO
            END IF
#endif
!     Solution is obtained, now updating (Corrector)
            CALL PICC

!     Checking for exceptions
            CALL EXCEPTIONS

!     Writing out the time passed, residual, and etc.
            IF (ALL(eq%ok)) EXIT
            CALL OUTRESULT(timeP, 2, iEqOld)
         END DO
!     End of inner loop

!     Saving the TXT files containing average and fluxes
         CALL TXT(.FALSE.)

         IF (rmsh%isReqd) THEN
            l1 = MOD(cTS,rmsh%cpVar) .EQ. 0
            IF (l1) THEN
               rmsh%rTS = cTS-1
               rmsh%time = time-dt
               rmsh%iNorm(:) = eq(:)%iNorm
               rmsh%A0(:,:) = Ao(:,:)
               rmsh%Y0(:,:) = Yo(:,:)
               rmsh%D0(:,:) = Do(:,:)
            END IF
         END IF

         IF (cm%mas()) INQUIRE(FILE=stopTrigName, EXIST=l1)
         CALL cm%bcast(l1)
         l2 = cTS .GE. nTS
!     Saving the result to restart simulation in the future
         IF (.NOT. stFileIncr .EQ. 0) THEN
            l3 = MOD(cTS,stFileIncr) .EQ. 0
            IF (l1 .OR. l2 .OR. l3) THEN
               fName = TRIM(stFileName)//"_last.bin"
               tmpS  = fName
               IF (.NOT.stFileRepl) THEN
                  WRITE(fName,'(I3.3)') cTS
                  IF (cTS .GE. 1000) fName = STR(cTS)
                  fName = TRIM(stFileName)//"_"//TRIM(fName)//".bin"
               END IF
               IF (cm%mas()) THEN
                  OPEN(fid, FILE=TRIM(fName))
                  CLOSE(fid, STATUS='DELETE')
               END IF
!     This call is to block all processors
               CALL cm%bcast(l1)
               OPEN(fid, FILE=TRIM(fName), ACCESS='DIRECT', RECL=recLn)
               IF (dFlag) THEN
                  IF (rmsh%isReqd .AND. saveAve) THEN
                     WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2                  CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An, Dn,
     3                  rmsh%Aav, rmsh%Yav, rmsh%Dav
                  ELSE
                     WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2                  CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An, Dn
                  END IF
               ELSE
                  WRITE(fid, REC=cm%tF()) stamp, cTS, time,
     2               CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An
               END IF
               CLOSE(fid)
               IF (.NOT.stFileRepl .AND. cm%mas()) THEN
                  CALL SYSTEM("ln -f "//TRIM(fName)//" "//TRIM(tmpS))
               END IF
            END IF
         END IF

         l3 = MOD(cTS,saveIncr) .EQ. 0
         l4 = saveFormat .NE. saveF_none
         l5 = cTS .GE. saveATS
!     Writing results into the disk with VTU format
         IF (l3 .AND. l4 .AND. l5) THEN
            CALL OUTRESULT(timeP, 3, iEqOld)
            CALL WRITEVTUS(An, Yn, Dn)
            IF (rmsh%isReqd .AND. saveAve) THEN
               rmsh%Aav = rmsh%Aav + An
               rmsh%Yav = rmsh%Yav + Yn
               rmsh%Dav = rmsh%Dav + Dn
            END IF
         ELSE
            CALL OUTRESULT(timeP, 2, iEqOld)
         END IF

!     Exiting outer loop if l1 or l2 happens
         IF (l1 .OR. l2) EXIT

!     Solution is stored here before replacing it at next time step
         Ao = An
         Yo = Yn
         IF (dFlag) Do = Dn
         cplBC%xo = cplBC%xn
      END DO
!     End of outer loop

      IF (resetSim) THEN
         CALL REMESHRESTART(timeP)
         DEALLOCATE(Ag, Yg, Dg, incL, res, isS)
         IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
            DEALLOCATE(ltgReordered, dirW, RTrilinos)
         END IF
         GOTO 101
      END IF

      IF (l2 .AND. saveAve) CALL CALCAVE

      CALL FINALIZE

#ifdef WITH_TRILINOS
      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
         DEALLOCATE(ltgReordered, dirW, RTrilinos)
      END IF
#endif
      DEALLOCATE(Ag, Yg, Dg, incL, res, isS)

      CALL MPI_FINALIZE(ierr)

      END PROGRAM MAIN

!--------------------------------------------------------------------

      SUBROUTINE STOPSIM()

      CALL FINALIZE
c      CALL cm%fStop()
      STOP "MPI is forced to stop by a fatal error"

      END SUBROUTINE STOPSIM

!--------------------------------------------------------------------
      SUBROUTINE INIT_DIR_AND_COUPNEU_BC(incL, res, dirW)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE
      !define input arguments
      INTEGER, INTENT(IN) :: incL(lhs%nFaces)
      REAL(KIND=8), INTENT(IN) :: res(lhs%nFaces)
      REAL(KIND=8), INTENT(INOUT) :: dirW(dof,tnNo)
      REAL(KIND=8), ALLOCATABLE :: v(:,:)

      INTEGER faIn, i, faDof, Ac
      LOGICAL flag, isCoupledBC

!     First initialize boundary flags as done in FSILS/SOLVE.f
      IF (lhs%nFaces .NE. 0) THEN
        lhs%face%incFlag = .TRUE.
        DO faIn=1, lhs%nFaces
!         used in dir bc condition
          IF (incL(faIn) .EQ. 0) lhs%face(faIn)%incFlag = .FALSE.
        END DO
!       used in coupled Neumann BC
        flag = ANY(lhs%face%bGrp.EQ.BC_TYPE_Neu)
        DO faIn=1, lhs%nFaces
          lhs%face(faIn)%coupledFlag = .FALSE.
          IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
          flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
          IF (flag .AND. res(faIn).NE.0D0) THEN
            lhs%face(faIn)%res = res(faIn)
            lhs%face(faIn)%coupledFlag = .TRUE. !first set the coupled flag
          END IF
        END DO
      END IF

!     First construct Dirichlet BC
!     dirW is matrix used to multiply in the Dirichlet boundary conditions
      dirW = 1D0 !will be used to be multiplied by as in PRECOND implementation
      DO faIn = 1, lhs%nFaces
        IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
        !copy over to global space
        faDof = MIN(lhs%face(faIn)%dof,dof)
        !check if the boundary face has a dirichlet bc on it
        IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
          DO i = 1, lhs%face(faIn)%nNo
!           uses sorted ordering Ac = lhs%map(gNodes(a))
            Ac = lhs%face(faIn)%glob(i)
            dirW(1:faDof, Ac) = dirW(1:faDof, Ac) *
     2                                 lhs%face(faIn)%val(1:faDof,i)
          END DO
        END IF
      END DO

      !Now implement coupled Neumann resistance BC for C
      ALLOCATE(v(dof,tnNo))
      v = 0D0
      isCoupledBC = .FALSE. !only copy over values for outer product if nonzeros exist
      !add flag if true and v is nonzero otherwise skip rotuine
      DO faIn = 1, lhs%nFaces
        IF (lhs%face(faIn)%coupledFlag) THEN !check for coupled condition
          isCoupledBC = .TRUE. !there are some nonzeros in v
          !copy over to global space
          faDof = MIN(lhs%face(faIn)%dof,dof)
          !check if the boundary face has a dirichlet bc on it
          DO i = 1, lhs%face(faIn)%nNo
            Ac = lhs%face(faIn)%glob(i) !uses sorted ordering Ac = lhs%map(gNodes(a))
            v(1:faDof,Ac) = v(1:faDof,Ac) + SQRT(res(faIn))* !will add on res*v*v'
     2                        lhs%face(faIn)%val(1:faDof,i)
          END DO
        END IF
      END DO
#ifdef WITH_TRILINOS
      CALL TRILINOS_BC_CREATE(v, isCoupledBC) !add flag argumentul
#endif
      DEALLOCATE(v)

      END SUBROUTINE INIT_DIR_AND_COUPNEU_BC

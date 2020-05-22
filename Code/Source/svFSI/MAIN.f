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

      LOGICAL l1, l2, l3, l4
      INTEGER(KIND=IKIND) i, iM, iBc, ierr, iEqOld
      REAL(KIND=RKIND) timeP(3)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incL(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Ag(:,:), Yg(:,:), Dg(:,:), res(:)

      IF (IKIND.NE.LSIP .OR. RKIND.NE.LSRP) THEN
         STOP "Incompatible datatype precision between solver and FSILS"
      END IF

      l1 = .FALSE.
      l2 = .FALSE.
      l3 = .FALSE.
      l4 = .FALSE.

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

      dbg = 'Allocating intermediate variables'
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo),
     2   res(nFacesLS), incL(nFacesLS))

!--------------------------------------------------------------------
!     Outer loop for marching in time. When entring this loop, all old
!     variables are completely set and satisfy BCs.
      IF (cTS .LE. nITS) dt = dt/10._RKIND
      DO
!     Adjusting the time step size once initialization stage is over
         IF (cTS .EQ. nITS) THEN
            dt = dt*10._RKIND
            std = " New time step size: "//dt
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

!     Apply Dirichlet BCs strongly
         CALL SETBCDIR(An, Yn, Dn)

!     Inner loop for iteration
         DO
            iEqOld = cEq

            IF (cplBC%coupled .AND. cEq.EQ.1) THEN
               CALL SETBCCPL
               CALL SETBCDIR(An, Yn, Dn)
            END IF

!     Initiator step (quantities at n+am, n+af)
            CALL PICI(Ag, Yg, Dg)
            IF (ALLOCATED(Rd)) THEN
               Rd = 0._RKIND
               Kd = 0._RKIND
            END IF

            dbg = 'Allocating the RHS and LHS'
            CALL LSALLOC(eq(cEq))

!           Compute body forces. If phys is shells or CMM (init),
!           apply contribution from body forces (pressure) to residue
            CALL SETBF(Dg)

            dbg = "Assembling equation <"//eq(cEq)%sym//">"
            DO iM=1, nMsh
               CALL GLOBALEQASSEM(msh(iM), Ag, Yg, Dg)
               dbg = "Mesh "//iM//" is assembled"
            END DO

!     Treatment of boundary conditions on faces
!     Apply Neumman or Traction boundary conditions
            CALL SETBCNEU(Yg, Dg)

!     Apply CMM BC conditions
            IF (.NOT.cmmInit) CALL SETBCCMM(Ag, Dg)

!     Apply weakly applied Dirichlet BCs
            CALL SETBCDIRW(Yg, Dg)

!     Apply contact model and add its contribution to residue
            IF (iCntct) CALL CONTACTFORCES(Dg)

!     Synchronize R across processes. Note: that it is important to
!     synchronize residue, R before treating immersed boundaries
            IF (.NOT.eq(cEq)%assmTLS) CALL COMMU(R)

!     Update residue in displacement equation for USTRUCT physics. Note
!     that this step is done only first iteration. The residue will be 0
!     for subsequent iterations
            IF (sstEq) CALL USTRUCTR(Yg)

            CALL SETBCUNDEFNEU()

!     Add contribution from IB to residue
            IF (ibFlag) CALL IB_CONSTRUCT()

            incL = 0
            IF (eq(cEq)%phys .EQ. phys_mesh) incL(nFacesLS) = 1
            IF (cmmInit) incL(nFacesLS) = 1
            DO iBc=1, eq(cEq)%nBc
               i = eq(cEq)%bc(iBc)%lsPtr
               IF (i .NE. 0) THEN
                  res(i) = eq(cEq)%gam*dt*eq(cEq)%bc(iBc)%r
                  incL(i) = 1
               END IF
            END DO
            IF (ibLSptr .NE. 0) incL(ibLSptr) = 1

            dbg = "Solving equation <"//eq(cEq)%sym//">"
            CALL LSSOLVE(eq(cEq), incL, res)

!     Solution is obtained, now updating (Corrector)
            CALL PICC

!     Checking for exceptions
            CALL EXCEPTIONS

!     Writing out the time passed, residual, and etc.
            IF (ALL(eq%ok)) EXIT
            CALL OUTRESULT(timeP, 2, iEqOld)
         END DO
!     End of inner loop

!     Immersed body treatment: project flow variables from fluid mesh
!     to IB solid mesh and enforce Dirichlet BCs
         IF (ibFlag) THEN
            CALL IB_PROJFVAR(An, Yn, Do, ib%An, ib%Yn, ib%Un)
            CALL IB_SETBCDIR(ib%An, ib%Yn, ib%Un)
         END IF

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
         l3 = MOD(cTS,stFileIncr) .EQ. 0
!     Saving the result to restart bin file
         IF (l1 .OR. l2 .OR. l3) CALL WRITERESTART(timeP)

!     Writing results into the disk with VTU format
         IF (saveVTK) THEN
            l3 = MOD(cTS,saveIncr) .EQ. 0
            l4 = cTS .GE. saveATS
            IF (l3 .AND. l4) THEN
               CALL OUTRESULT(timeP, 3, iEqOld)
               CALL WRITEVTUS(An, Yn, Dn)
               IF (ibFlag) CALL IB_WRITEVTUS(ib%An, ib%Yn, ib%Un)
            ELSE
               CALL OUTRESULT(timeP, 2, iEqOld)
            END IF
         ELSE
            CALL OUTRESULT(timeP, 2, iEqOld)
         END IF
         IF (pstEq) CALL OUTDNORM()

         IF (ibFlag) CALL IB_OUTR()

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
         DEALLOCATE(Ag, Yg, Dg, incL, res)
         IF (ALLOCATED(tls)) THEN
            DEALLOCATE(tls%ltg, tls%W, tls%R)
            DEALLOCATE(tls)
         END IF
         GOTO 101
      END IF

      IF (l2 .AND. saveAve) CALL CALCAVE

      DEALLOCATE(Ag, Yg, Dg, incL, res)
      CALL FINALIZE()
      CALL MPI_FINALIZE(ierr)

      END PROGRAM MAIN
!####################################################################
      SUBROUTINE STOPSIM()

      CALL FINALIZE
      STOP "MPI is forced to stop by a fatal error"

      END SUBROUTINE STOPSIM
!####################################################################

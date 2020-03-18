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
!     Communicator related module
!
!--------------------------------------------------------------------

      MODULE CMMOD
      USE TYPEMOD
      IMPLICIT NONE
      INCLUDE "mpif.h"

!     Size of blocks for openMP communications
      INTEGER(KIND=IKIND), PARAMETER :: mpBs = 1000
!     master is assumed to have zero ID
      INTEGER(KIND=IKIND), PARAMETER :: master = 0

!     Abstracted MPI names
      INTEGER(KIND=IKIND), PARAMETER :: mplog  = MPI_LOGICAL
      INTEGER(KIND=IKIND), PARAMETER :: mpint  = MPI_INTEGER
      INTEGER(KIND=IKIND), PARAMETER :: mpreal = MPI_DOUBLE_PRECISION
      INTEGER(KIND=IKIND), PARAMETER :: mpchar = MPI_CHARACTER

      TYPE cmType
         PRIVATE
!        Communicator handle
         INTEGER(KIND=IKIND) cHndl
!        Processors ID
         INTEGER(KIND=IKIND) taskId
!        Number of openMP threads in this cm
         INTEGER(KIND=IKIND) nThreads
!        Number of processors
         INTEGER(KIND=IKIND) nProcs
      CONTAINS
!        Create a new Communicator
         PROCEDURE, PUBLIC :: new => NEWCM
!        Returns commu handle
         PROCEDURE, PUBLIC :: com
!        Returns processor ID
         PROCEDURE, PUBLIC :: id => IDCM
!        Returns number of processors
         PROCEDURE, PUBLIC :: np => NUMPROC
!        Returns number of threads
         PROCEDURE, PUBLIC :: nT
!        Returns processor ID in fortran indexing
         PROCEDURE, PUBLIC :: tF
!        Returns true if this is master
         PROCEDURE, PUBLIC :: mas
!        Returns true if this is slave
         PROCEDURE, PUBLIC :: slv
!        Returns true if there is only one processor
         PROCEDURE, PUBLIC :: seq
!        Forces the MPI communicator to abort
         PROCEDURE, PUBLIC :: fStop
!        Broadcasting scaler/vector of logic/integer/real/character
         PROCEDURE :: BCASTLS
         PROCEDURE :: BCASTLV
         PROCEDURE :: BCASTIS
         PROCEDURE :: BCASTIV
         PROCEDURE :: BCASTIA2
         PROCEDURE :: BCASTRS
         PROCEDURE :: BCASTRV
         PROCEDURE :: BCASTRA2
         PROCEDURE :: BCASTSS
         PROCEDURE :: BCASTSV
         GENERIC :: bcast => BCASTLS, BCASTLV, BCASTIS, BCASTIV,
     2      BCASTIA2, BCASTRS, BCASTRV, BCASTRA2, BCASTSS, BCASTSV
!        Blocking MPI send
         PROCEDURE :: send => SENDRV
!        Blocking MPI recv
         PROCEDURE :: recv => RECVRV
!        Non-blocking MPI send
         PROCEDURE :: isend => ISENDRV
!        Non-blocking MPI recv
         PROCEDURE :: irecv => IRECVRV
!        Doing MPI wait on a set of requests
         PROCEDURE :: WAITS
         PROCEDURE :: WAITV
         GENERIC :: wait => WAITS, WAITV
!        Doing MPI reduce on a set of data types
         PROCEDURE :: REDUCEIS
         PROCEDURE :: REDUCEIV
         PROCEDURE :: REDUCERS
         PROCEDURE :: REDUCERV
         GENERIC :: reduce => REDUCEIS, REDUCEIV, REDUCERS, REDUCERV
!        Creating a communication handle from a set of processors
         PROCEDURE :: createCH
      END TYPE cmType

      INTERFACE ASSIGNMENT(=)
         MODULE PROCEDURE cmAssignCm
      END INTERFACE ASSIGNMENT(=)

      CONTAINS

!####################################################################
      SUBROUTINE NEWCM(cm, comHandle)
      IMPLICIT NONE
      CLASS(cmType) cm
      INTEGER(KIND=IKIND), INTENT(IN) :: comHandle

      LOGICAL ierr
      INTEGER(KIND=IKIND) i

      cm%cHndl  = comHandle
      cm%taskId = 0
      cm%nProcs = 1
      CALL MPI_COMM_RANK(comHandle, cm%taskId, ierr)
      CALL MPI_COMM_SIZE(comHandle, cm%nProcs, ierr)
      cm%nThreads = 1

      IF (cm%nProcs .EQ. 1) RETURN

      CALL MPI_TYPE_SIZE(mpint, i, ierr)
      IF (i .NE. IKIND) STOP "SIZE(MPI_INTEGER) .NE. default value"
      CALL MPI_TYPE_SIZE(mpreal, i, ierr)
      IF (i .NE. RKIND) STOP "SIZE(MPI_REAL) .NE. default value"
      CALL MPI_TYPE_SIZE(mplog, i, ierr)
      IF (i .NE. 4) STOP "SIZE(MPI_LOG) .NE. 4"
      CALL MPI_TYPE_SIZE(mpchar, i, ierr)
      IF (i .NE. 1) STOP "SIZE(MPI_CHAR) .NE. 1"

      RETURN
      END SUBROUTINE NEWCM
!--------------------------------------------------------------------
      FUNCTION COM(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND) COM

      COM = cm%cHndl

      RETURN
      END FUNCTION COM
!--------------------------------------------------------------------
      FUNCTION IDCM(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND) IDcm

      IDCM = cm%taskId

      RETURN
      END FUNCTION IDCM
!--------------------------------------------------------------------
      FUNCTION NUMPROC(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND) NUMPROC

      NUMPROC = cm%nProcs

      RETURN
      END FUNCTION NUMPROC
!--------------------------------------------------------------------
      FUNCTION NT(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND) NT

      NT = cm%nThreads

      RETURN
      END FUNCTION NT
!--------------------------------------------------------------------
      FUNCTION TF(cm)
      IMPLICIT NONE

      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND) TF

      TF = cm%taskId + 1

      RETURN
      END FUNCTION TF
!--------------------------------------------------------------------
      FUNCTION MAS(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL MAS

      IF (cm%taskId .NE. master) THEN
         MAS = .FALSE.
      ELSE
         MAS = .TRUE.
      END IF

      RETURN
      END FUNCTION MAS
!--------------------------------------------------------------------
      FUNCTION SLV(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL SLV

      IF (cm%taskId .NE. master) THEN
         SLV = .TRUE.
      ELSE
         SLV = .FALSE.
      END IF

      RETURN
      END FUNCTION SLV
!--------------------------------------------------------------------
      FUNCTION SEQ(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL SEQ

      SEQ = .FALSE.
      IF (cm%nProcs .EQ. 1) SEQ = .TRUE.

      RETURN
      END FUNCTION SEQ
!--------------------------------------------------------------------
      SUBROUTINE FSTOP(cm)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm

      LOGICAL ierr

      CALL MPI_ABORT(cm, MPI_ERR_OTHER, ierr)

      RETURN
      END SUBROUTINE FSTOP
!####################################################################
      SUBROUTINE BCASTLS(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL, INTENT(INOUT) :: u

      INTEGER(KIND=IKIND) ierr

      CALL MPI_BCAST(u, 1, mplog, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTLS
!--------------------------------------------------------------------
      SUBROUTINE BCASTLV(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      LOGICAL, INTENT(INOUT) :: u(:)

      INTEGER(KIND=IKIND) m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mplog, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTLV
!--------------------------------------------------------------------
      SUBROUTINE BCASTIS(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(INOUT) :: u

      INTEGER(KIND=IKIND) ierr

      CALL MPI_BCAST(u, 1, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIS
!--------------------------------------------------------------------
      SUBROUTINE BCASTIV(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(INOUT) :: u(:)

      INTEGER(KIND=IKIND) m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIV
!--------------------------------------------------------------------
      SUBROUTINE BCASTIA2(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(INOUT) :: u(:,:)

      INTEGER(KIND=IKIND) m, n, ierr

      m = SIZE(u,1)
      n = SIZE(u,2)
      CALL MPI_BCAST(u, m*n, mpint, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTIA2
!--------------------------------------------------------------------
      SUBROUTINE BCASTRS(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(INOUT) :: u

      INTEGER(KIND=IKIND) ierr

      CALL MPI_BCAST(u, 1, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRS
!--------------------------------------------------------------------
      SUBROUTINE BCASTRV(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(INOUT) :: u(:)

      INTEGER(KIND=IKIND) m, ierr

      m = SIZE(u)
      CALL MPI_BCAST(u, m, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRV
!--------------------------------------------------------------------
      SUBROUTINE BCASTRA2(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(INOUT) :: u(:,:)

      INTEGER(KIND=IKIND) m, n, ierr

      m = SIZE(u,1)
      n = SIZE(u,2)
      CALL MPI_BCAST(u, m*n, mpreal, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTRA2
!--------------------------------------------------------------------
      SUBROUTINE BCASTSS(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      CHARACTER(LEN=*), INTENT(INOUT) :: u

      INTEGER(KIND=IKIND) l, ierr

      l = LEN(u)
      CALL MPI_BCAST(u, l, mpchar, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTSS
!--------------------------------------------------------------------
      SUBROUTINE BCASTSV(cm, u)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      CHARACTER(LEN=*), INTENT(INOUT) :: u(:)

      INTEGER(KIND=IKIND) m, l, ierr

      m = SIZE(u)
      l = LEN(u)
      CALL MPI_BCAST(u, l*m, mpchar, master, cm%com(), ierr)

      RETURN
      END SUBROUTINE BCASTSV
!####################################################################
      SUBROUTINE SENDRV(cm, u, to, tag)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(IN) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN) :: to
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: tag

      INTEGER(KIND=IKIND) m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag

      ftag = ftag + cm%id()

      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_SEND(u, m, mpreal, to, ftag, cm%com(), ierr)

      RETURN
      END SUBROUTINE SENDRV
!--------------------------------------------------------------------
      SUBROUTINE RECVRV(cm, u, from, tag)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(OUT) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN) :: from
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: tag

      INTEGER(KIND=IKIND) m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag

      ftag = ftag + cm%id()

      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_RECV(u, m, mpreal, from, ftag, cm%com(),
     2   MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE RECVRV
!--------------------------------------------------------------------
      FUNCTION ISENDRV(cm, u, to, tag) RESULT(req)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(IN) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN) :: to
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: tag
      INTEGER(KIND=IKIND) req

      INTEGER(KIND=IKIND) m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag

      ftag = ftag + cm%id()

      req = MPI_REQUEST_NULL
      IF (to .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_ISEND(u, m, mpreal, to, ftag, cm%com(), req, ierr)

      RETURN
      END FUNCTION ISENDRV
!--------------------------------------------------------------------
      FUNCTION IRECVRV(cm, u, from, tag) RESULT(req)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(OUT) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN) :: from
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: tag
      INTEGER(KIND=IKIND) req

      INTEGER(KIND=IKIND) m, ierr, ftag

      ftag = 0
      IF (PRESENT(tag)) ftag = cm%np()*tag

      ftag = ftag + cm%id()

      req = MPI_REQUEST_NULL
      IF (from .EQ. cm%id()) RETURN
      m = SIZE(u)
      CALL MPI_IRECV(u, m, mpreal, from, ftag, cm%com(), req, ierr)

      RETURN
      END FUNCTION IRECVRV
!--------------------------------------------------------------------
      SUBROUTINE WAITS(cm, iReq)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(IN) :: iReq

      INTEGER(KIND=IKIND) ierr

      CALL MPI_WAIT(iReq, MPI_STATUS_IGNORE, ierr)

      RETURN
      END SUBROUTINE WAITS
!--------------------------------------------------------------------
      SUBROUTINE WAITV(cm, iReq)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(IN) :: iReq(:)

      INTEGER(KIND=IKIND) i, n

      n = SIZE(iReq)
      DO i=1, n
         CALL cm%wait(iReq(i))
      END DO

      RETURN
      END SUBROUTINE WAITV
!####################################################################
      FUNCTION REDUCERS(cm, u, iOp) RESULT(gU)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(IN) :: u
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: iOp
      REAL(KIND=RKIND) gU

      INTEGER(KIND=IKIND) ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, 1, mpreal, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCERS
!--------------------------------------------------------------------
      FUNCTION REDUCEIS(cm, u, iOp) RESULT(gU)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(IN) :: u
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: iOp
      INTEGER(KIND=IKIND) gU

      INTEGER(KIND=IKIND) ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, 1, mpint, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCEIS
!--------------------------------------------------------------------
      FUNCTION REDUCEIV(cm, u, iOp) RESULT(gU)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(IN) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: iOp
      INTEGER(KIND=IKIND), ALLOCATABLE :: gU(:)

      INTEGER(KIND=IKIND) n, ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      n = SIZE(u)
      ALLOCATE(gU(n))
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, SIZE(u), mpint, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCEIV
!--------------------------------------------------------------------
      FUNCTION REDUCERV(cm, u, iOp) RESULT(gU)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      REAL(KIND=RKIND), INTENT(IN) :: u(:)
      INTEGER(KIND=IKIND), INTENT(IN), OPTIONAL :: iOp
      REAL(KIND=RKIND), ALLOCATABLE :: gU(:)

      INTEGER(KIND=IKIND) n, ierr, op

      op = MPI_SUM
      IF (PRESENT(iOp)) op = iOp
      n = SIZE(u)
      ALLOCATE(gU(n))
      IF (cm%seq()) THEN
         gU = u
         RETURN
      END IF
      CALL MPI_ALLREDUCE(u, gU, n, mpreal, op, cm%com(), ierr)

      RETURN
      END FUNCTION REDUCERV
!####################################################################
      FUNCTION createCH(cm, pid) RESULT(ch)
      IMPLICIT NONE
      CLASS(cmType), INTENT(IN) :: cm
      INTEGER(KIND=IKIND), INTENT(IN) :: pid(:)
      INTEGER(KIND=IKIND) ch

      INTEGER(KIND=IKIND) n, cmGrp, newGrp, ierr

      n = SIZE(pid)
      CALL MPI_COMM_GROUP(cm%com(), cmGrp, ierr)
      CALL MPI_GROUP_INCL(cmGrp, n, pid, newGrp, ierr)
      CALL MPI_GROUP_FREE(cmGrp,ierr)
      CALL MPI_COMM_CREATE(cm%com(), newGrp, ch, ierr)
      CALL MPI_GROUP_FREE(newGrp,ierr)

      RETURN
      END FUNCTION createCH
!####################################################################
      SUBROUTINE cmAssignCm(s,r)
      IMPLICIT NONE
      CLASS(cmType), INTENT(OUT) :: s
      TYPE(cmType), INTENT(IN) :: r

      s%cHndl    = r%cHndl
      s%taskId   = r%taskId
      s%nThreads = r%nThreads
      s%nProcs   = r%nProcs

      RETURN
      END SUBROUTINE cmAssignCm
!####################################################################
      END MODULE CMMOD
!####################################################################

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
!     This routine initializes required structure for boundaries,
!     faces, those that interface with FSILS and cplBC.
!
!--------------------------------------------------------------------

      SUBROUTINE BAFINI()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) iM, iFa, iEq, i, a, iBc, lsPtr, faIn, faInCap

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:)

!     Compute face normals and area
      DO iM=1, nMsh
         DO iFa=1, msh(iM)%nFa
            IF (msh(iM)%lFib) CYCLE
            CALL FACEINI(msh(iM), msh(iM)%fa(iFa))
         END DO
         IF (msh(iM)%lShl) CALL SHLINI(msh(iM))
      END DO

!     Initialize face BC profile
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            CALL BCINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa))
            IF (msh(iM)%lShl .AND. (msh(iM)%eType.EQ.eType_TRI3)) THEN
               CALL SHLBCINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), msh(iM))
            END IF
         END DO
      END DO

!     cplBC faces are initialized here. Assuming that cplBC is always
!     in the first equation
      iEq = 1
      IF (ALLOCATED(cplBC%fa)) DEALLOCATE(cplBC%fa)
      IF (ALLOCATED(cplBC%xn)) DEALLOCATE(cplBC%xn)
      ALLOCATE(cplBC%fa(cplBC%nFa), cplBC%xn(cplBC%nX))
      IF (cplBC%coupled) THEN
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_cpl) .OR.
     2          BTEST(eq(iEq)%bc(iBc)%bType,bType_RCR)) THEN
               i = eq(iEq)%bc(iBc)%cplBCPtr
               cplBC%fa(i)%name = TRIM(msh(iM)%fa(iFa)%name)
               cplBC%fa(i)%y    = 0._RKIND
               IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Dir)) THEN
!                 Set cplBC internal flag to Dirichlet
                  cplBC%fa(i)%bGrp = cplBC_Dir 
               ELSE IF (BTEST(eq(iEq)%bc(iBc)%bType,bType_Neu)) THEN
!                 Set cplBC internal flag to Neumann
                  cplBC%fa(i)%bGrp = cplBC_Neu 
!                 If coupled scheme is implicit or semi-implicit, set bType to resistance
                  IF (cplBC%schm .NE. cplBC_E) eq(iEq)%bc(iBc)%bType=
     2               IBSET(eq(iEq)%bc(iBc)%bType,bType_res)

!                 Copy RCR structure from bc() to cplBC()
                  cplBC%fa(i)%RCR%Rp = eq(iEq)%bc(iBc)%RCR%Rp
                  cplBC%fa(i)%RCR%C  = eq(iEq)%bc(iBc)%RCR%C
                  cplBC%fa(i)%RCR%Rd = eq(iEq)%bc(iBc)%RCR%Rd
                  cplBC%fa(i)%RCR%Pd = eq(iEq)%bc(iBc)%RCR%Pd
                  cplBC%fa(i)%RCR%Xo = eq(iEq)%bc(iBc)%RCR%Xo
               ELSE
                  err = " Not a compatible cplBC_type"
               END IF
            END IF
         END DO
         IF (.NOT.stFileFlag) CALL RCRINIT()
         IF (cplBC%useGenBC) CALL genBC_Integ_X('I')
         IF (cplBC%schm .NE. cplBC_E) CALL CALCDERCPLBC()
      END IF

!     Setting up FSILS
      lsPtr = 0
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            eq(iEq)%bc(iBc)%lsPtr = 0
            CALL FSILSINI(eq(iEq)%bc(iBc), msh(iM)%fa(iFa), lsPtr)
!           Store mesh and face index in corresponding lhs%face(i). This
!           explicitly maps mesh and face indices between lhs object
!           (used in the linear solver for coupled surfaces) and the
!           msh object. Used in MATCHFACE below.
            IF (eq(iEq)%bc(iBc)%lsPtr .NE. 0) THEN 
               lhs%face(eq(iEq)%bc(iBc)%lsPtr)%iM = iM
               lhs%face(eq(iEq)%bc(iBc)%lsPtr)%iFa = iFa
            END IF
         END DO
      END DO
!     If any faces in msh(:)%fa(:) are capped, share information with
!     lhs%face(:)
      DO iEq=1,nEq
         DO iBc=1, eq(iEq)%nBc
            iFa = eq(iEq)%bc(iBc)%iFa
            iM  = eq(iEq)%bc(iBc)%iM
            IF (msh(iM)%fa(iFa)%capFaceID .NE. 0) THEN ! If face is capped
!              Find lhs%face(:) index of face being capped
               CALL MATCHFACE(iM, iFa, faIn)
!              Find lhs%face(:) index of capping face
               CALL MATCHFACE(iM, msh(iM)%fa(iFa)%capFaceID, faInCap)
!              Store capping relation in lhs%face(faIn)
               lhs%face(faIn)%faInCap = faInCap
            END IF
         END DO
      END DO

      IF (mvMsh) THEN
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct) .OR.
     2          ISDOMAIN(1, a, phys_ustruct) .OR.
     3          ISDOMAIN(1, a, phys_lElas)) i = i + 1
         END DO
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (ISDOMAIN(1, a, phys_struct) .OR.
     2          ISDOMAIN(1, a, phys_ustruct) .OR.
     3          ISDOMAIN(1, a, phys_lElas)) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         lsPtr = nFacesLS
         CALL FSILS_BC_CREATE(lhs, lsPtr, i, nsd, BC_TYPE_Dir, gNodes)
         DEALLOCATE(gNodes)
      END IF

      IF (cmmInit) THEN
         i = SUM(cmmBdry)
         ALLOCATE(gNodes(i))
         i = 0
         DO a=1, tnNo
            IF (cmmBdry(a) .EQ. 1) THEN
               i = i + 1
               gNodes(i) = a
            END IF
         END DO
         lsPtr = nFacesLS
         CALL FSILS_BC_CREATE(lhs, lsPtr, i, nsd, BC_TYPE_Dir, gNodes)
         DEALLOCATE(gNodes)
      END IF

      RETURN
      END SUBROUTINE BAFINI
!####################################################################
!     Initializing faces
      SUBROUTINE FACEINI(lM, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, g, i, Ac, Bc, Ec
      REAL(KIND=RKIND) area, sln, nV(nsd), v(nsd), xXi(nsd,nsd-1)
      TYPE(fsType) :: fs

      LOGICAL, ALLOCATABLE :: setIt(:)
      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), sA(:), sV(:,:)

!     Calculating face area
      ALLOCATE(sA(tnNo)) ! tnNo: Total number of nodes on this proc across all meshes
      sA   = 1._RKIND
      area = Integ(lFa, sA)
      std  = "    Area of face <"//TRIM(lFa%name)//"> is "//STR(area)
      IF (ISZERO(area)) THEN
         IF (cm%mas()) wrn = " <"//TRIM(lFa%name)//"> area is zero"
      END IF
      lFa%area = area
      DEALLOCATE(sA)

!     Compute face normals at nodes, stored in nV
      IF (ALLOCATED(lFa%nV)) DEALLOCATE(lFa%nV)
      ALLOCATE(lFa%nV(nsd,lFa%nNo), sV(nsd,tnNo))
      sV = 0._RKIND

      flag = .FALSE.
      IF (lM%eType.EQ.eType_TRI6  .OR. lM%eType.EQ.eType_QUD8  .OR.
     2    lM%eType.EQ.eType_QUD9  .OR. lM%eType.EQ.eType_TET10 .OR.
     3    lM%eType.EQ.eType_HEX20 .OR. lM%eType .EQ. eType_HEX27) THEN
         flag =.TRUE.
      END IF

      IF (.NOT.flag) THEN
!        For linear elements or NURBS, we simply project element normals
!        to nodes
         DO e=1, lFa%nEl
            IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(lM, lFa, e)
            DO g=1, lFa%nG
               CALL GNNB(lFa, e, g, nsd-1, lFa%eNoN, lFa%Nx(:,:,g), nV)
               DO a=1, lFa%eNoN
                  Ac       = lFa%IEN(a,e)
                  IF (Ac .NE. 0) THEN ! Ac could equal zero if virtual face
!                    Compute integral of normal vector over face
                     sV(:,Ac) = sV(:,Ac) + nV*lFa%N(a,g)*lFa%w(g)
                  END IF
               END DO
            END DO
         END DO

      ELSE
!        For higher order elements, use reduced order basis on mesh
!        to project element normals. Lumping method is used to project
!        to face corners. Normals at edge nodes are computed by simple
!        interpolation from reduced basis. Standard lumping using higher
!        order basis could lead to spurious errors
         CALL SETTHOODFS(fs, lFa%eType)
         CALL INITFS(fs, nsd-1)

         ALLOCATE(xl(nsd,lM%eNoN), ptr(lM%eNoN), setIt(lM%eNoN))
         DO e=1, lFa%nEl
            Ec = lFa%gE(e)
            setIt = .TRUE.
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               DO b=1, lM%eNoN
                  IF (setIt(b)) THEN
                     Bc = lM%IEN(b,Ec)
                     IF (Bc .EQ. Ac) EXIT
                  END IF
               END DO
               IF (b .GT. lM%eNoN) THEN
                  WRITE(*,'(A)') " ERROR: could not find matching "//
     2               "face node on higher order mesh"
                  CALL STOPSIM()
               END IF
               ptr(a) = b
               setIt(b) = .FALSE.
            END DO

            a = lFa%eNoN
            DO b=1, lM%eNoN
               IF (setIt(b)) THEN
                  a = a + 1
                  ptr(a) = b
               END IF
            END DO

            DO a=1, lM%eNoN
               Ac = lM%IEN(a,Ec)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Do(nsd+2:2*nsd+1,Ac)
            END DO

            DO g=1, fs%nG
               xXi = 0._RKIND
               DO a=1, fs%eNoN
                  b = ptr(a)
                  DO i=1, nsd-1
                     xXi(:,i) = xXi(:,i) + fs%Nx(i,a,g)*xl(:,b)
                  END DO
               END DO
               nV = CROSS(xXi)

               a  = ptr(1)
               b  = ptr(lFa%eNoN+1)
               v  = xl(:,a) - xl(:,b)
               IF (NORM(nV,v) .LT. 0._RKIND) nV = -nV

               DO a=1, fs%eNoN
                  Ac = lFa%IEN(a,e)
                  IF (Ac .NE. 0) THEN
!                    Compute integral of normal vector over face
                     sV(:,Ac) = sV(:,Ac) + fs%w(g)*fs%N(a,g)*nV(:)
                  END IF
               END DO
            END DO

            IF (MOD(lFa%eNoN,2) .EQ. 0) THEN
               g = lFa%eNoN
            ELSE
               g  = lFa%eNoN - 1
               Ac = lFa%IEN(lFa%eNoN,e)
               IF (Ac .NE. 0) THEN
                  DO b=1, fs%eNoN
                     Bc = lFa%IEN(b,e)
                     IF (Bc .NE. 0) THEN 
!                       Averaging sV? For higher order elements
                        sV(:,Ac) = sV(:,Ac) + sV(:,Bc)
                     END IF
                  END DO
                  sV(:,Ac) = sV(:,Ac)/REAL(fs%eNoN,KIND=RKIND)
               END IF
            END IF

            DO a=fs%eNoN+1, g
               b  = a - fs%eNoN
               Ac = lFa%IEN(a,e)
               Bc = lFa%IEN(b,e)
               IF ((Ac .NE. 0) .AND. (Bc .NE. 0)) THEN
                  nV = sV(:,Bc)
                  IF (b .EQ. fs%eNoN) THEN
                     Bc = lFa%IEN(1,e)
                  ELSE
                     Bc = lFa%IEN(b+1,e)
                  END IF
                  sV(:,Ac) = (nV + sV(:,Bc))*0.5_RKIND
               END IF
            END DO
         END DO
         DEALLOCATE(xl, ptr, setIt)
      END IF

      CALL COMMU(sV)
      flag = .TRUE.
      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         sln = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(sln)) THEN
            IF (flag) THEN
               wrn = "Skipping normal calculation of node "//a//
     2            " in face <"//TRIM(lFa%name)//">"
               flag = .FALSE.
            END IF
            lFa%nV(:,a) = 0._RKIND
            lFa%nV(1,a) = 1._RKIND
            CYCLE
         END IF
         lFa%nV(:,a) = sV(:,Ac)/sln
      END DO
      DEALLOCATE(sV)

      RETURN
      END SUBROUTINE FACEINI
!--------------------------------------------------------------------
      SUBROUTINE BCINI(lBc, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) iM, iFa, jFa, i, a, b, Ac, j, ierr
      REAL(KIND=RKIND) tmp, nV(nsd), center(nsd), maxN

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:), sCount(:), disp(:)
      REAL(KIND=RKIND), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = "Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      IF (BTEST(lBc%bType,bType_Robin) .AND. .NOT.dFlag) err =
     2    "Robin BC can be set for a displacement-based eqn only"

      iM  = lFa%iM
      iFa = lBc%iFa
      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lFa%nNo))

      ALLOCATE(s(tnNo), sCount(cm%np()), disp(cm%np()))
      s = 0._RKIND
      IF (BTEST(lBc%bType,bType_flat)) THEN
!     Just a constant value for Flat profile
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = 1._RKIND
         END DO
      ELSE IF (BTEST(lBc%bType,bType_para)) THEN
!     Here is the method that is used for imposing parabolic profile:
!     1- Find the coordinate of the points on the boundary 2- find unit
!     vector from center to each of points on the boundary: ew
!     3- maximize ew(i).e where e is the unit vector from current
!     point to the center 4- Use the point i as the diam here
         DO i=1, nsd
            center(i) = Integ(lFa, x, i)/lFa%area
         END DO
         ALLOCATE(gNodes(tnNo), sVl(nsd,lFa%nNo), sV(nsd,tnNo))
!     gNodes is one if a node located on the boundary (beside iFa)
         gNodes = 0
         DO jFa=1, msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, msh(iM)%fa(jFa)%nNo
               Ac         = msh(iM)%fa(jFa)%gN(a)
               gNodes(Ac) = 1
            END DO
         END DO
!     "j" is a counter for the number of nodes that are located on the
!     boundary of lFa and sVl contains the list of their coordinates
         j   = 0
         sVl = 0._RKIND
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (gNodes(Ac) .EQ. 1) THEN
               j = j + 1
               sVl(:,j) = x(:,Ac)
            END IF
         END DO
         IF (.NOT.cm%seq()) THEN
!     Getting the length data that is going to be received at each proc
            i = j*nsd
            CALL MPI_ALLGATHER(i, 1, mpint, sCount, 1, mpint, cm%com(),
     2         ierr)
            disp(1) = 0
            DO i=2, cm%np()
               disp(i) = disp(i-1) + sCount(i-1)
            END DO
            CALL MPI_ALLGATHERV(sVl, j*nsd, mpreal, sV, sCount, disp,
     2         mpreal, cm%com(), ierr)
            j = SUM(sCount)/nsd
         ELSE
            sV(:,1:j) = sVl(:,1:j)
         END IF
         IF (cm%mas() .AND. j.EQ.0) err = "No perimeter"//
     2      " found for face "//lFa%name
!     sVl will keep the normal unit vector from center to perimeter
         DEALLOCATE(sVl)
         ALLOCATE(sVl(nsd,j))
         DO a=1, j
            sV(:,a)  = sV(:,a) - center
            sVl(:,a) = sV(:,a)/SQRT(NORM(sV(:,a)))
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            nV = x(:,Ac) - center
            maxN = NORM(nV, sVl(:,1))
            i = 1
            DO b=2, j
               tmp = NORM(nV, sVl(:,b))
               IF (tmp .GT. maxN) THEN
                  maxN = tmp
                  i = b
               END IF
            END DO
            s(Ac) = 1._RKIND - NORM(nV)/NORM(sV(:,i))
         END DO
      ELSE IF (BTEST(lBc%bType,bType_ud)) THEN
         DO a=1, lFa%nNo
            Ac    = lFa%gN(a)
            s(Ac) = lBc%gx(a)
         END DO
      END IF

!     Now correcting the inlet BC for the inlet ring
      IF (BTEST(lBc%bType,bType_zp)) THEN
         DO jFa=1, msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, msh(iM)%fa(jFa)%nNo
               Ac    = msh(iM)%fa(jFa)%gN(a)
               s(Ac) = 0._RKIND
            END DO
         END DO
      END IF

!     Normalizing the profile for flux
      tmp = 1._RKIND
      IF (BTEST(lBc%bType,bType_flx)) THEN
         tmp = Integ(lFa, s)
         IF (ISZERO(tmp)) THEN
            tmp = 1._RKIND
            wrn = "Using face <"//TRIM(lFa%name)//
     2         "> to impose BC led to no non-zero node."
         END IF
      END IF

      DO a=1, lFa%nNo
         Ac        = lFa%gN(a)
         lBc%gx(a) = s(Ac)/tmp
      END DO

      RETURN
      END SUBROUTINE BCINI
!####################################################################
      SUBROUTINE FSILSINI(lBc, lFa, lsPtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(INOUT) :: lsPtr
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) a, e, Ac, g, iM, i, nNo
      REAL(KIND=RKIND) n(nsd)
      LOGICAL :: eDrn

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:)
      REAL(KIND=RKIND), ALLOCATABLE :: sV(:,:), sVl(:,:)

      iM  = lFa%iM
      nNo = lFa%nNo
      ALLOCATE(sVl(nsd,nNo), sV(nsd,tnNo), gNodes(nNo))

!     Copy mesh node id corresponding to face node id to gNodes
!     a in 1:nNo (nNo is the number of nodes on this face on this proc)
!     gN(a) in 1:tnNo (tnNo is the total number of nodes on this proc)
      DO a=1, nNo
         gNodes(a) = lFa%gN(a)
      END DO

      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (lBc%weakDir) THEN
            lBc%lsPtr = 0
         ELSE
            lsPtr     = lsPtr + 1
            lBc%lsPtr = lsPtr
            sVl = 0._RKIND
            eDrn = .FALSE.
            DO i=1, nsd
               IF (lBc%eDrn(i) .NE. 0) THEN
                  eDrn = .TRUE.
                  EXIT
               END IF
            END DO
            IF (eDrn) THEN
               sVl = 1._RKIND
               DO i=1, nsd
                  IF (lBc%eDrn(i) .NE. 0) sVl(i,:) = 0._RKIND
               END DO
            END IF
            CALL FSILS_BC_CREATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Dir,
     2         gNodes, sVl)
         END IF
      ELSE IF (BTEST(lBc%bType,bType_Neu)) THEN
!        Compute integral of normal vector over face in Moghadam et al. eq. 27. 
!        Note that this function is only computed once at initialization
         IF (BTEST(lBc%bType,bType_res)) THEN ! If resistance BC (or cpl BC)
            sV = 0._RKIND     ! The value of the integral
            DO e=1, lFa%nEl   ! Loop over elements on face
               IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(msh(iM),lFa,e) ! If NURBS
               DO g=1, lFa%nG ! Loop over Gauss point
                  CALL GNNB(lFa, e, g, nsd-1, lFa%eNoN, lFa%Nx(:,:,g),n) ! get weighted normal vector in ref config
                  DO a=1, lFa%eNoN     ! Loop over nodes in element
                     Ac = lFa%IEN(a,e) ! Extract global nodal index
!                    For a virtual face, Ac can be 0 if a proc does not own node Ac
                     IF (Ac .NE. 0) THEN 
!                       Integral of shape function times weighted normal
                        sV(:,Ac) = sV(:,Ac) + lFa%N(a,g)*lFa%w(g)*n 
                     END IF
                  END DO
               END DO
            END DO
            DO a=1, lFa%nNo
!              For a virtual face, Ac is obtained correctly as the index of the 
!              node local to this proc, Ac in [1, tnNo].
!              If lFa%nNo = 0, we do not enter loop and sVl is allocated with no space
               Ac       = lFa%gN(a)
               IF (Ac .EQ. 0) THEN
                  err = "Ac = 0"
               END IF
               sVl(:,a) = sV(:,Ac)
            END DO
            lsPtr     = lsPtr + 1
            lBc%lsPtr = lsPtr
!           Fills lhs%face(i) variables, including val if sVl exists
            CALL FSILS_BC_CREATE(lhs, lsPtr, lFa%nNo, nsd, BC_TYPE_Neu,
     2         gNodes, sVl, lFa%virtual)
         ELSE
            lBc%lsPtr = 0
         END IF
      ELSE IF (BTEST(lBc%bType,bType_trac)) THEN
         lBc%lsPtr = 0
      ELSE IF (BTEST(lBc%bType,bType_CMM)) THEN
         lsPtr     = lsPtr + 1
         lBc%lsPtr = lsPtr

         nNo = 0
         DO a=1, lFa%nNo
            IF (ISZERO(lBc%gx(a))) nNo = nNo + 1
         END DO
         DEALLOCATE(gNodes, sVl)
         ALLOCATE(sVl(nsd,nNo), gNodes(nNo))
         sVl  = 0._RKIND

         eDrn = .FALSE.
         DO i=1, nsd
            IF (lBc%eDrn(i) .NE. 0) THEN
               eDrn = .TRUE.
               EXIT
            END IF
         END DO
         IF (eDrn) THEN
            sVl = 1._RKIND
            DO i=1, nsd
               IF (lBc%eDrn(i) .NE. 0) sVl(i,:) = 0._RKIND
            END DO
         END IF

         nNo = 0
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (ISZERO(lBc%gx(a))) THEN
               nNo = nNo + 1
               gNodes(nNo) = Ac
            END IF
         END DO

         CALL FSILS_BC_CREATE(lhs, lsPtr, nNo, nsd, BC_TYPE_Dir, gNodes,
     2      sVl)
      ELSE
         err = "Unxpected bType in FSILSINI"
      END IF

      DEALLOCATE(sVl, sV, gNodes)

      RETURN
      END SUBROUTINE FSILSINI
!####################################################################
!     Compute shell extended IEN for triangular elements. Reqd. to
!     resolve bending moments.
      SUBROUTINE SETSHLXIEN(lM)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: a, b, e, f, Ac, Bc, nEl, eNoN, ep(2,3)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incN(:)

      eNoN = lM%eNoN
      nEl  = lM%nEl

      ep   = RESHAPE((/2,3,3,1,1,2/), SHAPE(ep))
      ALLOCATE(incN(eNoN), lM%eIEN(eNoN,nEl), lM%sbc(eNoN,nEl))

      lM%eIEN = 0
      lM%sbc  = 0
      DO e=1, nEl
         DO a=1, eNoN
            Ac = lM%IEN(ep(1,a),e)
            Bc = lM%IEN(ep(2,a),e)
            DO f=1, nEl
               IF (e .EQ. f) CYCLE
               incN = 0
               DO b=1, eNoN
                  IF (lM%IEN(b,f).EQ.Ac .OR. lM%IEN(b,f).EQ.Bc)
     2               incN(b) = incN(b) + 1
               END DO
               IF (SUM(incN) .EQ. 2) THEN
                  DO b=1, eNoN
                     IF (incN(b) .EQ. 0) THEN
                        lM%eIEN(a,e) = lM%IEN(b,f)
                        EXIT
                     END IF
                  END DO
                  EXIT
               END IF
            END DO
            IF (lM%eIEN(a,e) .EQ. 0) THEN
               lM%sbc(a,e) = IBSET(lM%sbc(a,e), bType_free)
            END IF
         END DO
      END DO
      DEALLOCATE(incN)

      RETURN
      END SUBROUTINE SETSHLXIEN
!--------------------------------------------------------------------
!     Initializing shell normals and area
      SUBROUTINE SHLINI(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: a, e, g, Ac, nNo, nEl, eNoN
      LOGICAL :: flag
      REAL(KIND=RKIND) :: Jac, area, nV(nsd), tmpR(nsd,nsd-1)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), sV(:,:)

      nNo  = lM%nNo
      nEl  = lM%nEl
      eNoN = lM%eNoN

!     Compute shell director (normal)
      ALLOCATE(xl(nsd,eNoN), sV(nsd,tnNo), lM%nV(nsd,nNo))
      sV   = 0._RKIND
      area = 0._RKIND
      DO e=1, nEl
         IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)
         DO a=1, eNoN
            Ac = lM%IEN(a,e)
            xl(:,a) = x(:,Ac)
         END DO

         DO g=1, lM%nG
            CALL GNNS(eNoN, lM%Nx(:,:,g), xl, nV, tmpR, tmpR)
            Jac = SQRT(NORM(nV))
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               sV(:,Ac) = sV(:,Ac) + lM%w(g)*lM%N(a,g)*nV
               area     = area     + lM%w(g)*lM%N(a,g)*Jac
            END DO
         END DO
      END DO
      area = cm%reduce(area)
      CALL COMMU(sV)

      std = " Area of the shell surface <"//TRIM(lM%name)//"> is "//
     2   STR(area)

      flag = .TRUE.
      DO a=1, nNo
         Ac  = lM%gN(a)
         Jac = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(Jac)) THEN
            IF (flag) THEN
               wrn = "Skipping normal calculation of node "//a//
     2            " in mesh <"//TRIM(lM%name)//">"
               flag = .FALSE.
            END IF
            lM%nV(:,a) = 0._RKIND
            lM%nV(1,a) = 1._RKIND
            CYCLE
         END IF
         lM%nV(:,a) = sV(:,Ac)/Jac
      END DO

      DEALLOCATE(xl, sV)

      RETURN
      END SUBROUTINE SHLINI
!--------------------------------------------------------------------
!     Initializing boundary condition variables for CST shells - this
!     BC is set on the interior node adjacent to the boundary
      SUBROUTINE SHLBCINI(lBc, lFa, lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(IN) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER(KIND=IKIND) :: a, b, e, Ac, Bc, Ec
      LOGICAL :: bFlag

      DO e=1, lFa%nEl
         Ec = lFa%gE(e)
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,Ec)

!           Find the interior node of the mesh element that is part of
!           the boundary. [bFlag=='T' => node is on the boundary]
            bflag = .FALSE.
            DO b=1, lFa%eNoN
               Bc = lFa%IEN(b,e)
               IF (Ac .EQ. Bc) THEN
                  bFlag = .TRUE.
                  EXIT
               END IF
            END DO

!           Set the BC on the interior node [bFlag == 'F']
            IF (.NOT.bFlag) THEN
               IF (.NOT.BTEST(lM%sbc(a,Ec),bType_free)) THEN
                  err = "BC detected on a non-boundary shell element."//
     2               " Correction needed"
               END IF
               lM%sbc(a,Ec) = IBCLR(lM%sbc(a,Ec), bType_free)
               IF (BTEST(lBc%bType,bType_free)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_free)
               ELSE IF (BTEST(lBc%bType,bType_fix)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_fix)
               ELSE IF (BTEST(lBc%bType,bType_hing)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_hing)
               ELSE IF (BTEST(lBc%bType,bType_symm)) THEN
                  lM%sbc(a,Ec) = IBSET(lM%sbc(a,Ec),bType_symm)
               END IF
               EXIT
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE SHLBCINI
!--------------------------------------------------------------------
!     Find lhs%face(:) index faIn corresponding iM and iFa, which
!     index msh(:)%fa(:)
      SUBROUTINE MATCHFACE(iM, iFa, faIn)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: iM, iFa
      INTEGER(KIND=IKIND), INTENT(OUT) :: faIn

      INTEGER(KIND=IKIND) a
!     Loop over lhs%faces to find corresonding face
      DO a=1, lhs%nFaces
!        If lhs%face matches mesh and face index, 
         IF ((lhs%face(a)%iFa .EQ. iFa)
     2      .AND. (lhs%face(a)%iM .EQ. iM)) THEN
            faIn = a
         END IF
      END DO
      END SUBROUTINE
!####################################################################


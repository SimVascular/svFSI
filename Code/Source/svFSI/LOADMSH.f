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
!     This routine is desinged to read isoparametric meshes and
!     retype it upon request
!
!--------------------------------------------------------------------
!     This reads coordinate, connectivity, and ebc/vtk files
      SUBROUTINE READCCNE(list, lM)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list
      TYPE(mshType), INTENT(INOUT) :: lM

      INTEGER i, fid, Ac, e, iFa
      CHARACTER(LEN=stdL) stmp
      TYPE(listType), POINTER :: lPtr, lPBC
      TYPE(fileType) ftmp

!     Reading connectivity file
      lPtr => list%get(ftmp,"Connectivity file path")
!     First testing to find out whether CCNE is the mesh format
      IF (.NOT.ASSOCIATED(lPtr)) RETURN
      fid = ftmp%open()

      dbg = " Reading coordinate/connectivity based mesh"
      READ (fid,"(A)") stmp
      lM%eNoN = CheckNoNumbers(stmp) - 1
!     Based on eNoN and nsd, finding the eType
      CALL SELECTELE(lM)

      lM%gnEl = 1
      DO
         READ (fid,*,END=113)
         lM%gnEl = lM%gnEl + 1
      END DO
 113  REWIND(fid)

      ALLOCATE (lM%gIEN(lM%eNoN,lM%gnEl))
      DO e=1, lM%gnEl
         READ (fid,*) i, lM%gIEN(:,e)
      END DO
      CLOSE (fid)

!     Reading coordinates file
      lPtr => list%get(ftmp,"Coordinates file path",1)
      fid = ftmp%open()
      READ (fid,"(A)") stmp
      IF (CheckNoNumbers(stmp) .NE. nsd+1) THEN
         err = "nsd is not consistent with coordinate file"
      END IF
      lM%gnNo = 1
      DO
         READ (fid,*,END=111)
         lM%gnNo = lM%gnNo + 1
      END DO
 111  REWIND(fid)

!     allocating the required space for solver
      ALLOCATE (x(nsd,lM%gnNo))
      DO Ac=1, lM%gnNo
         READ (fid,*) i, x(:,Ac)
      END DO
      CLOSE (fid)

!     Satisfying right-handed property for the internal elements and
!     making sure everything is compatible
      IF (ichckIEN) CALL CHECKIEN(lM, .FALSE.)

!     Reading faces names and setting the required parameters
      lM%nFa = list%srch("Add face")
      std = " Number of available faces: "//lM%nFa
      ALLOCATE (lM%fa(lM%nFa))

      DO iFa=1, lM%nFa
         lPBC => list%get(lM%fa(iFa)%name,"Add face",iFa)

!     First searching for a vtk file
         lPtr => lPBC%get(ftmp,"vtk file path")
         IF (ASSOCIATED(lPtr)) THEN
            CALL RDBCVTK(lM, lM%fa(iFa), ftmp%open())
         ELSE
!     Reading face connectivity file
            lPtr => lPBC%get(ftmp,"Connectivity file (ebc) path",1)
            fid = ftmp%open()
            READ (fid,"(A)",END=110) stmp
            lM%fa(iFa)%eNoN = CheckNoNumbers(stmp) - 2
            CALL SELECTELEB(lM, lM%fa(iFa))
            lM%fa(iFa)%nEl = 1
            DO
               READ (fid,*,END=110)
               lM%fa(iFa)%nEl = lM%fa(iFa)%nEl + 1
            END DO
 110        ALLOCATE (lM%fa(iFa)%gE(lM%fa(iFa)%nEl),
     2         lM%fa(iFa)%IEN(lM%fa(iFa)%eNoN,lM%fa(iFa)%nEl))
            REWIND(fid)
            DO e=1, lM%fa(iFa)%nEl
               READ (fid,*) lM%fa(iFa)%gE(e),i,lM%fa(iFa)%IEN(:,e)
            END DO
            CLOSE (fid)
            CALL CALCNBC(lM, lM%fa(iFa))
         END IF
         
         lM%fa(iFa)%gnEl = lM%fa(iFa)%nEl
         ALLOCATE(lM%fa(iFa)%gebc(1+lM%fa(iFa)%eNoN,lM%fa(iFa)%gnEl))
         lM%fa(iFa)%gebc(1,:) = lM%fa(iFa)%gE(:)
         lM%fa(iFa)%gebc(2:1+lM%fa(iFa)%eNoN,:) = lM%fa(iFa)%IEN(:,:)
      END DO
      
      RETURN
      END SUBROUTINE READCCNE

!####################################################################
!     This reads msh file produced by Gambit
      SUBROUTINE READGAMBIT(list, lM)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(listType), INTENT(INOUT) :: list
      TYPE(mshType), INTENT(INOUT) :: lM

      TYPE blkType
         INTEGER nEl
         INTEGER id
         INTEGER facePtr
      END TYPE blkType

      TYPE ldmnType
         INTEGER iDmn
         INTEGER :: blk = 0
         CHARACTER(LEN=stdL) name
      END TYPE ldmnType

      INTEGER, PARAMETER :: maxNblk = 1024

      LOGICAL flag, nonL
      INTEGER A, Ac, Ec, b, i, j, k, l, nBlk, cntr, nCon, iFa, e, fid,
     2   eNoNb, nlDmn
      INTEGER(KIND=8) tmp
      CHARACTER(LEN=stdL) fName, stmp, ctmp
      TYPE(blkType) blk(maxNblk)
      TYPE(listType),POINTER :: lPtr, lPD
      TYPE(fileType) ftmp

      INTEGER, ALLOCATABLE :: adj(:,:), ef(:)
      REAL(KIND=8), ALLOCAtABLE :: Xtmp(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: names(:)
      TYPE(ldmnType), ALLOCATABLE :: ldmn(:)

      lPtr => list%get(ftmp,"Gambit mesh file path")
      IF (.NOT.ASSOCIATED(lPtr)) RETURN
!     If that is not the case, we continue with this format of mesh
      dbg = " Reading a Gambit formate mesh"
      fid = ftmp%open()

      nonL     = .FALSE.
      lM%eNoN  = 4
      eNoNb    = 2
      nonL = .FALSE.
      lPtr => list%get(nonL,"Convert elements to biquadratic elements")
      IF (nonL) THEN
         lM%eNoN = 9
         eNoNb   = 3
      END IF
      CALL SELECTELE(lM)

      nlDmn = list%srch("Set domain")
      ALLOCATE (ldmn(nlDmn))
      DO i=1, nlDmn
         lPD => list%get(ldmn(i)%name,"Set domain",i)
         lPtr => lPD%get(ldmn(i)%iDmn,"Domain",1,
     2      ll=0,ul=BIT_SIZE(dmnId)-1)
      END DO
!     This is used for quick reading through the stuff that I don't need
!     later on
      cntr = 0
!     First reading coordinate data
      DO
         cntr = cntr + 1
         READ(fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(5:13) .EQ. "Dimension") EXIT
      END DO
      cntr = cntr + 1
      READ(fid,"(A)") stmp
      CALL GET(stmp, i)
      CALL GET(stmp, i)
      IF (i .NE. nsd) err = "nsd does not match with msh file"
      DO
         cntr = cntr + 1
         READ (fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(1:3) .EQ. "(10") EXIT
      END DO
      cntr = cntr + 1
      READ(fid,"(A)") stmp
      DO i=1, 4
         CALL GET(stmp, tmp)
      END DO
      lM%gnNo = INT(tmp)
!     allocating the required space for solver
      ALLOCATE (x(nsd,lM%gnNo))
      DO A=1, lM%gnNo
         cntr = cntr + 1
         READ (fid,*) x(:,A)
      END DO

!     Now reading data block by block
      DO
         cntr = cntr + 1
         READ(fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(5:9) .EQ. "Faces") EXIT
      END DO
      l    = 0
      nBlk = 0
      DO
         READ(fid,"(A)") stmp
         READ(fid,"(A)") stmp
         IF (stmp(1:1) .NE. "(") EXIT
         nBlk = nBlk + 1
         IF (nBlk .GT. maxnBlk) err = "Increase maxnBlk"
         CALL GET(stmp, tmp)
         CALL GET(stmp, tmp)
         blk(nBlk)%id = INT(tmp)
         CALL GET(stmp, tmp)
         i = INT(tmp)
         CALL GET(stmp, tmp)
         i = INT(tmp) - i + 1
         blk(nBlk)%nEl = i
         READ(fid,"(A)") stmp
         l = CheckNoNumbers(stmp) - 1
         IF (l .NE. 4) err = "GMBT - Unrecognized element type"
         DO j=1, blk(nBlk)%nEl-1
            READ(fid,*)
         END DO
      END DO
      IF (nBlk .EQ. 0) err = "No block found in .msh file"
      DO
         READ(fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(5:9) .EQ. "Zones") EXIT
      END DO
      ALLOCATE (names(nBlk))
      nCon   = 0
      lM%nFa = 0
      DO
         READ(fid,'(A)',END=015) stmp
         CALL GET(stmp, tmp)
         CALL GET(stmp, j)
         READ(stmp,*) ctmp, fName
         fName = ADJUSTR(fName)
         fName = ADJUSTL(fName(:stdL-4))
!     Corresponding the name of the face to its identifier
         DO i=1, nBlk
            IF (j .EQ. blk(i)%id) THEN
!     Identifying the type of this block, it can be either connectivity
!     or surfaces connectivity data
               IF (ctmp .EQ. "interior") THEN
!     Using number of connectivity as the upper limit for number of
!     elements
                  nCon = blk(nBlk)%nEl
                  ALLOCATE (lM%IEN(lM%eNoN,nCon))
                  lM%IEN = 0
                  blk(i)%facePtr = 0
                  IF (nonL) THEN
                     ALLOCATE (adj(lM%nEf,nCon))
                     adj = 0
                  END IF
               ELSE
                  names(i)       = fName
                  lM%nFa         = lM%nFa + 1
                  blk(i)%facePtr = lM%nFa
               END IF
               EXIT
            END IF
         END DO
!     This is a domain
         IF (i .GT. nBlk) THEN
            DO i=1, nlDmn
               IF (ldmn(i)%name .EQ. fName) ldmn(i)%blk = j
            END DO
         END IF
      END DO
 015  REWIND(fid)

      DO i=1, nlDmn
         IF (ldmn(i)%blk .EQ. 0) err = "Failed to find domain <"//
     2      TRIM(ldmn(i)%name)//">"
      END DO
      IF (nCon .EQ. 0) err = "No iterior found in .msh file"
      std = " Number of external surfaces: "//lM%nFa
      ALLOCATE (lM%fa(lM%nFa), ef(l))
      lM%fa%eNoN = eNoNb
!     Passing through coordinate data
      DO i=1, cntr
         READ(fid,*)
      END DO

      lM%gnEl = 0
!     Reading connectivity and ebc data
      DO i=1, nBlk
         iFa = blk(i)%facePtr
         IF (iFa .NE. 0) THEN
            lM%fa(iFa)%name = names(i)
            lM%fa(iFa)%nEl = blk(i)%nEl
            ALLOCATE (lM%fa(iFa)%gE(lM%fa(iFa)%nEl),
     2         lM%fa(iFa)%IEN(eNoNb,lM%fa(iFa)%nEl))
            lM%fa(iFa)%gE = 0
         END IF
         READ(fid,"(A)") stmp
         READ(fid,"(A)") stmp
         DO e=1, blk(i)%nEl
            READ(fid,"(A)") stmp
            CALL GET(stmp,j) ! This is garbage
            DO a=1, l
               CALL GET(stmp, tmp)
               ef(a) = INT(tmp)
            END DO
!        Producing adjacency data in case of nonlinear elements
            IF (nonL .AND. ef(l-1)*ef(l).NE.0) THEN
               DO j=1, lM%nEf
                  IF (adj(j,ef(l)) .EQ. 0) THEN
                     adj(j,ef(l)) = ef(l-1)
                     DO k=1, lM%nEf
                        IF (adj(k,ef(l-1)) .EQ. 0) THEN
                           adj(k,ef(l-1)) = ef(l)
                           EXIT
                        END IF
                     END DO
                     EXIT
                  ELSE IF (adj(j,ef(l)) .EQ. ef(l-1)) THEN
                     EXIT
                  END IF
               END DO
            END IF
!     Each face is connected to two elements
            DO j=l-1, l
               Ec = ef(j)
               IF (Ec .EQ. 0) CYCLE
               IF (Ec .GT. lM%gnEl) THEN
                  lM%gnEl = Ec
                  IF (lM%gnEl .GT. nCon) THEN
                     err = "nCon is surprisingly less than gnEl"
                  END IF
               END IF
               DO k=1, lM%eNoN
                  IF (lM%IEN(k,Ec) .EQ. 0) THEN
                     DO a=1, l-2
                        Ac = ef(a)
                        flag = .TRUE.
                        DO b=1, k-1
                           IF (lM%IEN(b,Ec) .EQ. Ac) flag = .FALSE.
                        END DO
                        IF (flag) THEN
                           lM%IEN(k,Ec) = Ac
                           EXIT
                        END IF
                     END DO
                  END IF
               END DO
!     If this is a face, we set
               IF (iFa .NE. 0) THEN
                  DO a=1, l-2
                     lM%fa(iFa)%IEN(a,e) = ef(a)
                  END DO
                  IF (Ec .GT. lM%fa(iFa)%gE(e)) lM%fa(iFa)%gE(e) = Ec
               END IF
            END DO
         END DO
      END DO
      DO
         READ(fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(5:9) .EQ. "Cells") EXIT
      END DO
      DO
         READ(fid,'(A)') stmp
         stmp = ADJUSTL(stmp)
         IF (stmp(1:1) .NE. "(") EXIT
         CALL GET(stmp, tmp)
         CALL GET(stmp, j)
         DO i=1, nlDmn
            IF (ldmn(i)%blk .EQ. j) EXIT
         END DO
         IF (i .LE. nlDmn) THEN
!     Saving the range of elements in each domain so we can set the dmn
!     after constructing gIEN
            CALL GET(stmp, tmp)
            j = INT(tmp)
            CALL GET(stmp, tmp)
            k = INT(tmp)
            CALL SETDMNID(lM, ldmn(i)%iDmn, j, k)
         END IF
      END DO
      CLOSE(fid)

!     Constructing IEN array
      ALLOCATE (lM%gIEN(lM%eNoN,lM%gnEl))
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            lM%gIEN(a,e) = lM%IEN(a,e)
         END DO
      END DO
      DEALLOCATE (lM%IEN)

!     Satisfying right-handed property for the internal elements and
!     making sure everything is compatible
      IF (ichckIEN) CALL CHECKIEN(lM, nonL)

!     Adding intermediate nodes to the elements
      IF (nonL) THEN
!     Upper estimation for gnNo
         ALLOCATE (Xtmp(nsd,lM%eNoN*lM%gnEl))
         DO a=1, lM%gnNo
            Xtmp(:,a) = x(:,a)
         END DO
         DEALLOCATE (x)
         IF (lM%eType .EQ. eType_BIQ) THEN
!     We will start adding intermediate nodes from faces. Here I assume
!     there is no overlapping between faces
            DO iFa=1, lM%nFa
               DO e=1, lM%fa(iFa)%nEl
                  i  = lM%fa(iFa)%IEN(1,e)
                  j  = lM%fa(iFa)%IEN(2,e)
                  Ec = lM%fa(iFa)%gE(e)
                  CALL ADDBETWEEN(lM, i, j, Ec, .TRUE., Xtmp)
                  lM%fa(iFa)%IEN(3,e) = lM%gnNo
                  DO k=1, lM%nEf
                     Ac = adj(k,Ec)
                     IF (Ac .NE. 0) CALL ADDBETWEEN(lM, i, j, Ac,
     2                  .FALSE., Xtmp)
                  END DO
               END DO
            END DO
!     Now adding nodes for interior elements
            DO e=1, lM%gnEl
               DO a=1, 4
                  i = lM%gIEN(a,e)
                  j = lM%gIEN(a+1,e)
                  IF (a .EQ. 4) j = lM%gIEN(1,e)
                  CALL ADDBETWEEN(lM, i, j, e, .TRUE., Xtmp)
                  DO k=1, lM%nEf
                     Ec = adj(k,e)
                     IF (Ec .NE. 0) CALL ADDBETWEEN(lM, i, j, Ec,
     2                  .FALSE., Xtmp)
                  END DO
               END DO
            END DO
!     And now the node for the middle of elements
            DO e=1, lM%gnEl
               lM%gnNo = lM%gnNo + 1
               lM%gIEN(9,e) = lM%gnNo
               i = lM%gIEN(1,e)
               j = lM%gIEN(3,e)
               Xtmp(:,lM%gnNo) = (Xtmp(:,i) + Xtmp(:,j))/2D0
            END DO
            ALLOCATE (x(nsd,lM%gnNo))
            DO a=1, lM%gnNo
               x(:,a) = Xtmp(:,a)
            END DO
            DEALLOCATE (Xtmp)
         ELSE
            err = "READFILE correction is required"
         END IF
      END IF

      DO iFa=1, lM%nFa
         CALL CALCNBC(lM, lM%fa(iFa))
         CALL SELECTELEB(lM, lM%fa(iFa))
         
         lM%fa(iFa)%gnEl = lM%fa(iFa)%nEl
         ALLOCATE(lM%fa(iFa)%gebc(1+lM%fa(iFa)%eNoN,lM%fa(iFa)%gnEl))
         lM%fa(iFa)%gebc(1,:) = lM%fa(iFa)%gE(:)
         lM%fa(iFa)%gebc(2:1+lM%fa(iFa)%eNoN,:) = lM%fa(iFa)%IEN(:,:)
      END DO

      RETURN
      END SUBROUTINE READGAMBIT

!####################################################################
!     Adds node "a" between two given nodes "a1" and "a2" in element "e"
      SUBROUTINE ADDBETWEEN(lM, a1, a2, e, flag, Xtmp)

      USE COMMOD

      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      LOGICAL, INTENT(IN) :: flag
      INTEGER, INTENT(IN) :: a1, a2, e
      REAL(KIND=8), INTENT(INOUT) :: Xtmp(nsd,lM%eNoN*lM%gnEl)

      INTEGER b1, b2, i

      b1 = 0
      b2 = 0
      DO i=1, lM%eNoN
         IF (lM%gIEN(i,e) .EQ. a1) THEN
            b1 = i
         ELSEIF(lM%gIEN(i,e) .EQ. a2) THEN
            b2 = i
         END IF
      END DO
      IF (b1 .GT. b2) THEN
         i = b1
         b1 = b2
         b2 = i
      END IF
      IF (b1 .EQ. 0) RETURN

      IF (lM%eType .EQ. eType_BIQ) THEN
         IF (b1.EQ.1 .AND. b2.EQ.2) THEN
            i = 5
         ELSE IF (b1.EQ.2 .AND. b2.EQ.3) THEN
            i = 6
         ELSE IF (b1.EQ.3 .AND. b2.EQ.4) THEN
            i = 7
         ELSE IF (b1.EQ.1 .AND. b2.EQ.4) THEN
            i = 8
         ELSE
            err = "You should not be here"
            RETURN
         END IF
         IF (flag .AND. lM%gIEN(i,e).EQ.0) THEN
            lM%gnNo = lM%gnNo + 1
            Xtmp(:,lM%gnNo) = (Xtmp(:,a1) + Xtmp(:,a2))/2D0
         END IF
         IF (lM%gIEN(i,e) .EQ. 0) lM%gIEN(i,e) = lM%gnNo
      ELSE
         err = "eType insertion needs to be specified"
      END IF

      RETURN
      END SUBROUTINE ADDBETWEEN

!####################################################################
!     Reading BC (ebc/nbc) from a VTK file
      SUBROUTINE RDBCVTK(lM, lFa, fid)

      USE COMMOD
      USE ALLFUN

      IMPLICIT NONE

      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa
      INTEGER, INTENT(IN) :: fid

      INTEGER e, a, b, Ac, ie, maxnEtN
      CHARACTER(LEN=stdL) rLine

      INTEGER, ALLOCATABLE :: nbc(:), ndEs(:,:), bcNdEs(:,:), gnID(:)

      DO
         READ (fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:10) .EQ. 'POINT_DATA') EXIT
      END DO
      rLine = rLine(11:)
      READ (rLine,*) a
      ALLOCATE (nbc(a))
      DO
         READ (fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
!     Continuing until to get to a block of nodal data
         IF (rLine(1:7) .EQ. 'SCALARS') THEN
            READ (fid,'(a)',END=001) rLine
            EXIT
         ELSE
            CYCLE
         END IF
      END DO
!     Reading global ID of the nodes
      READ(fid,*) nbc

      REWIND(fid)
      DO
         READ (fid,'(a)',END=001) rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:8) .EQ. 'POLYGONS') EXIT
      END DO
      rLine = rLine(9:)
      READ (rLine,*) lFa%nEl
      READ (fid,'(a)',END=001) rLine
      lFa%eNoN = CheckNoNumbers(rLine) - 1
      CALL SELECTELEB(lM, lFa)
      ALLOCATE (lFa%gE(lFa%nEl), lFa%IEN(lFa%eNoN,lFa%nEl),
     2   gnID(lFa%eNoN))
      READ (rLine,*) a, gnID
      IF (a .NE. lFa%eNoN) err = "Inconsistent number of columns"
      DO a=1, lFa%eNoN
         lFa%IEN(a,1) = nbc(gnID(a)+1)
      END DO
!     And reading the ebc data
      DO e=2, lFa%nEl
         READ (fid,*,END=001) a, gnID
         DO a=1, lFa%eNoN
            lFa%IEN(a,e) = nbc(gnID(a)+1)
         END DO
      END DO
      CLOSE (fid)
!     Since vtk might not contain ele%gl values, I am going to
!     construct it myself based on gIEN and face%IEN
!     To do this, I am going to construct an array, ndEs, with size
!     gnNo*(max no element touching a node) that row "a" contains
!     elements that have "a" as a node. Here "nbc" array keeps track of
!     the length ndEs
      DEALLOCATE (nbc)
      ALLOCATE (nbc(lM%gnNo))
      maxnEtN = 3*lM%eNoN**3
 003  maxnEtN = maxnEtN + 5
      IF (ALLOCATED(ndEs)) DEALLOCATE (ndEs)
      ALLOCATE (ndEs(lM%gnNo,maxnEtN))
      nbc = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            nbc(Ac) = nbc(Ac) + 1
            IF (nbc(Ac) .GT. maxnEtN) THEN
               wrn = "Increasing maxnEtN to "//maxnEtN
               GOTO 003
            END IF
            ndEs(Ac,nbc(Ac)) = e
         END DO
      END DO
!     Now I am going to construct another array that contains the
!     elements that touch boundary nodes in the first index and the
!     number of repetition of those elements in the second index
      ALLOCATE (bcNdEs(maxnEtN*lFa%eNoN,2))
      DO e=1, lFa%nEl
         bcNdEs = 0
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            DO ie=1, nbc(Ac)
               DO b=1, maxnEtN*lFa%eNoN
                  IF (bcNdEs(b,1) .EQ. 0) THEN
                     bcNdEs(b,1) = ndEs(Ac,ie)
                     bcNdEs(b,2) = 1
                     EXIT
                  ELSE IF (bcNdEs(b,1) .EQ. ndEs(Ac,ie)) THEN
                     bcNdEs(b,2) = bcNdEs(b,2) + 1
                     EXIT
                  END IF
               END DO
            END DO
         END DO
!     Now if the number of repetition of an element is equal to face
!     eNoN, then that element should have a face shared with the
!     boundary element
         DO b=1, maxnEtN*lFa%eNoN
            IF (bcNdEs(b,2) .EQ. lFa%eNoN) THEN
               lFa%gE(e) = bcNdEs(b,1)
               EXIT
            ELSE IF (bcNdEs(b,2) .EQ. 0) THEN
               err = "global element number for element <"//e
     2            //"> of face <"//TRIM(lFa%name)//"> was not found"
            END IF
         END DO
      END DO
!     And calculating nbc
      CALL CALCNBC(lM, lFa)

      RETURN
 001  err = "A block of data is missing from BC vtk file"

      RETURN
      END SUBROUTINE RDBCVTK



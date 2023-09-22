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
!     These set of routines load fiber directions for a mesh and
!     computes fiber direction at either a node, an element, or
!     at an integration point depending on how it is stored.
!
!--------------------------------------------------------------------

!####################################################################
!     Read fiber orientation. Here we will try to load fibers in
!     multiple ways. First we will attempt to read fibers from a single
!     file. The code will look for variables named as "FIB_DIR1",
!     "FIB_DIR2", etc. in a single vtu file.
!     If not found, the code will then search for multiple fiber files
!     (vtu format) with fiber direction variable named as "FIB_DIR".
!     Otherwise, we will search for a uniform fiber vector with the
!     keyword "Fiber direction".
      SUBROUTINE READ_FIBERS(list, lM, fib)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(fibDirType), INTENT(INOUT) :: fib

      LOGICAL :: fib_alloc
      INTEGER(KIND=IKIND) :: a, b, i, j
      REAL(KIND=RKIND) :: fibN(nsd)
      CHARACTER(LEN=stdL) :: ctmp
      TYPE(listType), POINTER :: lPtr

      std = " ------------------------------------"
      std = " Searching for fiber data.. "

      fib_alloc = .FALSE.
      fib%nFn = 0

!     Check if num fiber directions are provided by the user
!     If not, the code will try to compute nFn from the variable names
!     in the vtu file, or from a text file
      lPtr => list%get(fib%nFn, "Number of fiber directions")

!     First attempt to read fibers from a single vtu file.
      lPtr => list%get(ctmp, "Fiber directions file path")
      IF (ASSOCIATED(lPtr)) THEN
         IF (rmsh%isReqd) err = "Fiber directions read from a file"//
     2      " is not allowed with remeshing"
         i = LEN(TRIM(ctmp))
         IF (ctmp(i-3:i) .EQ. '.vtu') THEN
            CALL READFIBNSVTU(lM, fib, ctmp)
         ELSE
            CALL READFIBNSTXT(lM, fib, ctmp)
         END IF
         IF (ALLOCATED(fib%fN)) fib_alloc = .TRUE.
      END IF

!     If the fibers are not found, we will try searching for multiple
!     vtu-based fiber files with keyword FIB_DIR. We will also determine
!     if the fibers are located at nodes or elements.
      IF (.NOT.fib_alloc) THEN
         fib%nFn = list%srch("Fiber direction file path")

         IF (fib%nFn .NE. 0) THEN
            fib_alloc = .TRUE.
            IF (rmsh%isReqd) err = "Fiber directions read from file "//
     2         "is not allowed with remeshing"

            DO i=1, fib%nFn
               lPtr => list%get(ctmp, "Fiber direction file path", i)
               CALL READFIBNFF(lM, fib, ctmp, "FIB_DIR", i)
            END DO
         END IF ! msh%nFn
      END IF ! fib_alloc

!     If fibers are still not found, look for a prescribed constant
!     vector and store them at element nodes
      IF (.NOT.fib_alloc) THEN
         fib%nFn = list%srch("Fiber direction")

         IF (fib%nFn .NE. 0) THEN
            fib%locNd = .TRUE.
            fib_alloc = .TRUE.
            ALLOCATE(fib%fN(nsd*fib%nFn,1,lM%gnNo))
            fib%fN = 0._RKIND

            DO i=1, fib%nFn
               lPtr => list%get(fibN, "Fiber direction", i)
               b = (i-1)*nsd
               DO a=1, lM%gnNo
                  DO j=1, nsd
                     fib%fN(b+j,1,a) = fibN(j)
                  END DO
               END DO
            END DO
         END IF ! fib%nFn
      END IF ! fib_alloc

      IF (fib_alloc) THEN
         std = "   Found "//STR(fib%nFn)//" fiber directions"//
     2      " for <"//TRIM(lM%name)//">"
         IF (fib%locNd) THEN
            std = "   Fibers are located at mesh nodes"
         ELSE IF (fib%locEl) THEN
            std = "   Fibers are located at mesh elements"
         ELSE IF (fib%locGP) THEN
            std = "   Fibers are located at the integration points"
         END IF

!        Normalizing fiber directions
         CALL FIBER_NORMALIZE(lM, fib)
      ELSE
         std = " ..fibers not found"
      END IF
      std = " ------------------------------------"

!     For 3-noded (tri3) shells, fibers should be at elements only
      IF (lM%lShl .AND. lM%eType.EQ.eType_TRI3
     2    .AND. .NOT.fib%locEl) THEN
         err = " Fibers are required to be element centroids for "//
     2      "triangular shell surfaces"
      END IF

      RETURN
      END SUBROUTINE READ_FIBERS
!--------------------------------------------------------------------
!     Read multiple fiber directions from a single vtu file with
!     variable names FIB_DIR1, FIB_DIR2, etc.
!     This can be used if fibers are stored at nodes or at element
!     centers but not integration points.
      SUBROUTINE READFIBNSVTU(lM, fib, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(fibDirType), INTENT(INOUT) :: fib
      CHARACTER(LEN=*) :: fName

      TYPE(vtkXMLType) :: vtu
      INTEGER(KIND=IKIND) :: a, e, i, is, ie, nvar, nFn, istat
      CHARACTER(LEN=stdL) :: stmp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: varNames(:)

      istat = 0
      std = "   <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, istat)
      IF (istat .LT. 0) err = " VTU file read error (init)"

      CALL getVTK_numPoints(vtu, a, istat)
      IF (a .NE. lM%gnNo) THEN
         err = " Mismatch in num nodes while loading fibers"
      END IF

      CALL getVTK_numElems(vtu, e, istat)
      IF (e .NE. lM%gnEl) THEN
         err = " Mismatch in num elems while loading fibers"
      END IF

!     We will first check if the fibers are stored at elements
      nvar = 0
      nFn  = 0
      CALL getVTK_numElemData(vtu, nvar, istat)
      IF (nvar .NE. 0) THEN
         ALLOCATE(varNames(nvar))
         CALL getVTK_elemDataNames(vtu, varNames, istat)
         IF (istat .LT. 0) err = " VTU file read error (elemDataNames)"

         DO i=1, nvar
            stmp = varNames(i)
            IF (stmp(1:7) .EQ. 'FIB_DIR') THEN
               READ(stmp(8:8),*,IOSTAT=istat) e
               IF (istat .NE. 0) err = " Cannot find fiber directions"
               IF (e .GT. nFn) nFn = e
            END IF
         END DO
         DEALLOCATE(varNames)

!        If found, allocate appropriate arrays
         IF (nFn .NE. 0) THEN
            fib%locEl = .TRUE.
            IF (fib%nFn .EQ. 0) THEN
               fib%nFn = nFn
            ELSE IF (fib%nFn .NE. nFn) THEN
               err = " Inconsistent number of fibers detected"
            END IF
            ALLOCATE(fib%fN(nsd*nFn,1,lM%gnEl))
            fib%fN = 0._RKIND
         END IF
      END IF

!     If the fibers are not found at elements, we will look for fibers
!     at nodal locations
      IF (nFn .EQ. 0) THEN
         nvar  = 0
         istat = 0
         CALL getVTK_numPointData(vtu, nvar, istat)
         IF (istat .LT. 0) THEN
            err = " VTU file read error (numPointData). Could not "//
     2         "find fibers in the vtu file <"//TRIM(fName)//">"
         END IF

         ALLOCATE(varNames(nvar))
         CALL getVTK_pointDataNames(vtu, varNames, istat)
         IF (istat .LT. 0) err = " VTU file read error (pointDataNames)"

         DO i=1, nvar
            stmp = varNames(i)
            IF (stmp(1:7) .EQ. 'FIB_DIR') THEN
               READ(stmp(8:8),*,IOSTAT=istat) a
               IF (a .GT. nFn) nFn = a
            END IF
         END DO
         DEALLOCATE(varNames)

!        Throw error if no fiber directions are found in the given file
         IF (nFn .EQ. 0) THEN
            err = " Could not find fiber directions in <"//TRIM(fName)//
     2         ">. Check variable names (FIB_DIR)."
         ELSE
!        If fibers are found at nodes, allocate appropriate arrays
            fib%locNd = .TRUE.
            IF (fib%nFn .EQ. 0) THEN
               fib%nFn = nFn
            ELSE IF (fib%nFn .NE. nFn) THEN
               err = " Inconsistent number of fibers detected"
            END IF
            ALLOCATE(fib%fN(nsd*nFn,1,lM%gnNo))
            fib%fN = 0._RKIND
         END IF
      END IF

      IF (fib%locEl) THEN
         ALLOCATE(tmpR(maxNSD,lM%gnEl))
         DO i=1, fib%nFn
            WRITE(stmp,'(A)') "FIB_DIR"//STR(i)
            CALL getVTK_elemData(vtu, TRIM(stmp), tmpR, istat)
            IF (istat .LT. 0) err = " VTU file read error "//TRIM(stmp)

            is = (i-1)*nsd + 1
            ie = i*nsd
            DO e=1, lM%gnEl
               fib%fN(is:ie,1,e) = tmpR(1:nsd,e)
            END DO
         END DO
         DEALLOCATE(tmpR)

      ELSE IF (fib%locNd) THEN
         ALLOCATE(tmpR(maxNSD,lM%gnNo))
         DO i=1, fib%nFn
            WRITE(stmp,'(A)') "FIB_DIR"//STR(i)
            CALL getVTK_pointData(vtu, TRIM(stmp), tmpR, istat)
            IF (istat .LT. 0) err = " VTU file read error "//TRIM(stmp)

            is = (i-1)*nsd + 1
            ie = i*nsd
            DO a=1, lM%gnNo
               fib%fN(is:ie,1,a) = tmpR(1:nsd,a)
            END DO
         END DO
         DEALLOCATE(tmpR)

      ELSE
         err = " Unexpected behavior while loading fibers from vtu file"
      END IF

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READFIBNSVTU
!--------------------------------------------------------------------
!     Read multiple fiber directions from a txt file.
      SUBROUTINE READFIBNSTXT(lM, fib, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(fibDirType), INTENT(INOUT) :: fib
      CHARACTER(LEN=*) :: fName

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: a, e, g, i, j, fid, nFn, nnL, slen, nToks,
     2   istat
      CHARACTER(LEN=1024) :: sline
      CHARACTER(LEN=stdL), DIMENSION(250) :: tokList


      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:,:)

      flag = .TRUE.
      INQUIRE (FILE=TRIM(fName), EXIST=flag)
      IF (.NOT.flag) THEN
         err = "Fiber directions file <"//TRIM(fName)//"> does not"//
     2      " exist or can not be opened"
      END IF

      fid = 1367
      OPEN(fid,FILE=TRIM(fName))

!     Determine the size of nFn
      IF (fib%nFn .EQ. 0) THEN
         nFn = 0
         DO
            READ(fid,'(A)',IOSTAT=istat) sline
            IF (istat .NE. 0) THEN
               err = " Error reading fiber file <"//TRIM(fName)//">"
            END IF
            sline = ADJUSTC(sline)
            slen  = LEN(TRIM(sline))
            IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

            CALL parseString(sline, tokList, nToks)
            nFn = (nToks - (nsd+2))/nsd
            IF (nFn.LT.1 .OR. nFn.GT.2) THEN
               err = " Inconsistent format in fiber directions file <"//
     2            TRIM(fName)//">"
            END IF
            fib%nFn = nFn
            EXIT
         END DO
      ELSE
         nFn = fib%nFn
      END IF

      std = " Num fibers: "//STR(nFn)

!     Find the number of lines to determing if the fibers are located
!     at nodes, elements, or at the integration points
      nnL = 0
      REWIND(fid)
      DO
         READ(fid,'(A)',IOSTAT=istat) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         IF (istat .GT. 0) THEN
            err = " Error reading fiber file <"//TRIM(fName)//">"
         ELSE IF (istat .LT. 0) THEN ! End of File !
            EXIT
         END IF
         nnL = nnL + 1
      END DO

      std = " Num Lines: "//STR(nnL)

!     Read all the fibers into a temporaary array
      ALLOCATE(tmpR(nsd*nFn,nnL))
      tmpR = 0._RKIND
      nnL  = 0
      REWIND(fid)
      DO
         READ(fid,'(A)',IOSTAT=istat) sline
         sline = ADJUSTC(sline)
         slen  = LEN(TRIM(sline))
         IF (sline(1:1).EQ.'#' .OR. slen.EQ.0) CYCLE

         IF (istat .GT. 0) THEN
            err = " Error reading fiber file <"//TRIM(fName)//">"
         ELSE IF (istat .LT. 0) THEN ! End of File !
            EXIT
         END IF
         nnL = nnL + 1

         CALL parseString(sline, tokList, nToks)
         IF (nToks .NE. (2+(nFn+1)*nsd)) THEN
            err = " Inconsistent format in fiber directions file <"//
     2         TRIM(fName)//">"
         END IF

         DO i=1, nFn
            a = (i-1)*nsd
            e = (2+nsd) + a
            DO j=1, nsd !is, ie
               READ(tokList(e+j),*,IOSTAT=istat) tmpR(a+j,nnL)
               IF (istat .NE. 0) THEN
                  err = " Error copying fiber direction to arrays"
               END IF
            END DO
         END DO

      END DO
      CLOSE(fid)

!     Determine where the fibers are located based on no. of lines
      IF (nnL .EQ. lM%gnNo) THEN
         fib%locNd = .TRUE.
         ALLOCATE(fib%fN(nsd*nFn,1,lM%gnNo))
         fib%fN = 0._RKIND
         DO a=1, lM%gnNo
            fib%fN(:,1,a) = tmpR(:,a)
         END DO

      ELSE IF (nnL .EQ. lM%gnEl) THEN
         fib%locEl = .TRUE.
         ALLOCATE(fib%fN(nsd*nFn,1,lM%gnEl))
         fib%fN = 0._RKIND
         DO e=1, lM%gnEl
            fib%fN(:,1,e) = tmpR(:,e)
         END DO

      ELSE IF (nnL .EQ. (lM%nG*lM%gnEl)) THEN
         fib%locGP = .TRUE.
         ALLOCATE(fib%fN(nsd*nFn,lM%nG,lM%gnEl))
         fib%fN = 0._RKIND
         i = 0
         DO e=1, lM%gnEl
            DO g=1, lM%nG
               i = i + 1
               fib%fN(:,g,e) = tmpR(:,i)
            END DO
         END DO

      ELSE
         err = " Error reading fiber file <"//TRIM(fName)//">. "//
     2      "Could not determine where fibers are stored."

      END IF

      DEALLOCATE(tmpR)

      RETURN
      END SUBROUTINE READFIBNSTXT
!--------------------------------------------------------------------
!     Read fiber direction from a vtu file with variable name FIB_DIR
!     This can be used if fibers are stored at nodes or at element
!     centers but not integration points.
      SUBROUTINE READFIBNFF(lM, fib, fName, kwrd, idx)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(fibDirType), INTENT(INOUT) :: fib
      CHARACTER(LEN=*) :: fName, kwrd
      INTEGER(KIND=IKIND), INTENT(IN) :: idx

      TYPE(vtkXMLType) :: vtu
      INTEGER(KIND=IKIND) :: a, e, i, is, ie, nvar, istat
      CHARACTER(LEN=stdL) :: stmp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: varNames(:)

      istat = 0
      std = "   <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = " VTU file read error (init)"

      CALL getVTK_numPoints(vtu, a, istat)
      IF (a .NE. lM%gnNo) THEN
         err = " Mismatch in num nodes while loading fibers"
      END IF

      CALL getVTK_numElems(vtu, e, iStat)
      IF (e .NE. lM%gnEl) err = " Mismatch in num elems for "//
     2   TRIM(kwrd)

!     If fib%fN is not allocated, we will first check if the fibers
!     are stored at elements. If not, check for nodes.
      IF (.NOT.ALLOCATED(fib%fN)) THEN
         nvar = 0
         CALL getVTK_numElemData(vtu, nvar, istat)
         IF (nvar .NE. 0) THEN
            ALLOCATE(varNames(nvar))
            CALL getVTK_elemDataNames(vtu, varNames, istat)
            IF (istat .LT. 0) THEN
               err = " VTU file read error (elemDataNames)"
            END IF

            DO i=1, nvar
               stmp = varNames(i)
               IF (stmp(1:7) .EQ. 'FIB_DIR') THEN
                  fib%locEl = .TRUE.
                  ALLOCATE(fib%fN(nsd*fib%nFn,1,lM%gnEl))
                  fib%fN = 0._RKIND
                  EXIT
               END IF
            END DO
            DEALLOCATE(varNames)
         END IF

!        If fibers are not lcoated at elements, we assume that the
!        fibers are located at points. If not found at points, an error
!        is thrown later.
         IF (.NOT.fib%locEl) THEN
            istat = 0
            fib%locNd = .TRUE.
            ALLOCATE(fib%fN(nsd*fib%nFn,1,lM%gnNo))
            fib%fN = 0._RKIND
         END IF
      END IF ! fib%fN

      IF (fib%locEl) THEN
         ALLOCATE(tmpR(maxNSD,lM%gnEl))
         tmpR = 0._RKIND
         istat = 0
         CALL getVTK_elemData(vtu, TRIM(kwrd), tmpR, iStat)
         IF (istat .LT. 0) err = " VTU file read error "//TRIM(kwrd)

         is = (idx-1)*nsd + 1
         ie = idx*nsd
         DO e=1, lM%gnEl
            fib%fN(is:ie,1,e) = tmpR(1:nsd,e)
         END DO
         DEALLOCATE(tmpR)

      ELSE IF (fib%locNd) THEN
         ALLOCATE(tmpR(maxNSD,lM%gnNo))
         tmpR = 0._RKIND
         istat = 0
         CALL getVTK_pointData(vtu, TRIM(kwrd), tmpR, iStat)
         IF (istat .LT. 0) THEN
            err = " VTU file read error "//TRIM(kwrd)//
     2         " Fibers are not found in <"//TRIM(fName)//"> at "//
     3         " nodes and cell centers"
         END IF

         is = (idx-1)*nsd + 1
         ie = idx*nsd
         DO a=1, lM%gnNo
            fib%fN(is:ie,1,a) = tmpR(1:nsd,a)
         END DO
         DEALLOCATE(tmpR)

      END IF

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READFIBNFF
!####################################################################
!     Normalizes fiber directions to unit magnitude
      SUBROUTINE FIBER_NORMALIZE(lM, fib)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(fibDirType), INTENT(INOUT) :: fib

      INTEGER(KIND=IKIND) :: a, e, g, i, is, ie
      REAL(KIND=RKIND) :: rtmp, fibN(nsd)

      IF (fib%locNd) THEN
         DO a=1, lM%gnNo
            DO i=1, fib%nFn
               is = (i-1)*nsd + 1
               ie = i*nsd
               fibN = fib%fN(is:ie,1,a)
               rtmp = SQRT(NORM(fibN))
               IF (.NOT.ISZERO(rtmp)) THEN
                  fib%fN(is:ie,1,a) = fibN(:) / rtmp
               END IF
            END DO
         END DO

      ELSE IF (fib%locEl) THEN
         DO e=1, lM%gnEl
            DO i=1, fib%nFn
               is = (i-1)*nsd + 1
               ie = i*nsd
               fibN = fib%fN(is:ie,1,e)
               rtmp = SQRT(NORM(fibN))
               IF (.NOT.ISZERO(rtmp)) THEN
                  fib%fN(is:ie,1,e) = fibN(:) / rtmp
               END IF
            END DO
         END DO

      ELSE IF (fib%locGP) THEN
         DO e=1, lM%gnEl
            DO g=1, lM%nG
               DO i=1, fib%nFn
                  is = (i-1)*nsd + 1
                  ie = i*nsd
                  fibN = fib%fN(is:ie,g,e)
                  rtmp = SQRT(NORM(fibN))
                  IF (.NOT.ISZERO(rtmp)) THEN
                     fib%fN(is:ie,g,e) = fibN(:) / rtmp
                  END IF
               END DO
            END DO
         END DO

      END IF

      RETURN
      END SUBROUTINE FIBER_NORMALIZE
!####################################################################
!     Computes fiber direction at the integration point depending on
!     whether the fiber are stored directly at the integration point,
!     elements-centroids, or at nodal locations.
      SUBROUTINE GET_FIBN(lM, fib, e, g, eNoN, N, fN)
      USE COMMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(fibDirType), INTENT(IN) :: fib
      INTEGER(KIND=IKIND), INTENT(IN) :: e, g, eNoN
      REAL(KIND=RKIND), INTENT(IN) :: N(eNoN)
      REAL(KIND=RKIND), INTENT(OUT) :: fN(nsd,fib%nFn)

      INTEGER(KIND=IKIND) :: a, b, i, j

      REAL(KIND=RKIND), ALLOCATABLE :: fNl(:,:)

!     Initialize
      fN = 0._RKIND

!     Fibers located at integration points
      IF (fib%locGP) THEN
         DO i=1, fib%nFn
            a = (i-1)*nsd + 1
            b = i*nsd
            fN(1:nsd,i) = fib%fN(a:b,g,e)
         END DO

!     Fibers located at element centroids
      ELSE IF (fib%locEl) THEN
         DO i=1, fib%nFn
            a = (i-1)*nsd + 1
            b = i*nsd
            fN(1:nsd,i) = fib%fN(a:b,1,e)
         END DO

!     Fibers located at mesh nodes
      ELSE IF (fib%locNd) THEN
         ALLOCATE(fNl(nsd*fib%nFn,eNoN))
         fNl = 0._RKIND
         DO a=1, eNoN
            i = lM%IEN(a,e)
            b = lM%lN(i)
            fNl(:,a) = fib%fN(:,1,b)
         END DO

         DO i=1, fib%nFn
            b = (i-1)*nsd
            DO j=1, nsd
               DO a=1, eNoN
                  fN(j,i) = fN(j,i) + N(a)*fNl(b+j,a)
               END DO
            END DO
         END DO
         DEALLOCATE(fNl)
      END IF

      RETURN
      END SUBROUTINE GET_FIBN
!####################################################################

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
!     Here data (mesh, solution) input/output is handled interfacing
!     with fortran-based VTK module.
!
!--------------------------------------------------------------------

      SUBROUTINE READVTU(lM, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND) :: iStat
      TYPE(vtkXMLType) :: vtu

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (init)"

      CALL getVTK_numPoints(vtu, lM%gnNo, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num points)"

      CALL getVTK_numElems(vtu, lM%gnEl, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num cells)"

      CALL getVTK_nodesPerElem(vtu, lM%eNoN, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (nodes per cell)"

      ALLOCATE(lM%x(nsd,lM%gnNo),tmpX(maxNSD,lM%gnNo))
      CALL getVTK_pointCoords(vtu, tmpX, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (coords)"
      lM%x(:,:) = tmpX(1:nsd,:)
      DEALLOCATE(tmpX)

      ALLOCATE(lM%gIEN(lM%eNoN,lM%gnEl))
      CALL getVTK_elemIEN(vtu, lM%gIEN, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (ien)"
      lM%gIEN = lM%gIEN + 1

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READVTU
!--------------------------------------------------------------------
      SUBROUTINE READVTP(lFa, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND) :: iStat, a, e, Ac
      TYPE(vtkXMLType) :: vtp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtp, fName, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (init)"

      CALL getVTK_numPoints(vtp, lFa%nNo, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (num points)"

      CALL getVTK_numElems(vtp, lFa%nEl, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (num cells)"

      CALL getVTK_nodesPerElem(vtp, lFa%eNoN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (nodes per cell)"

      ALLOCATE(lFa%x(nsd,lFa%nNo), tmpX(maxNSD,lFa%nNo))
      CALL getVTK_pointCoords(vtp, tmpX, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (coords)"
      lFa%x(:,:) = tmpX(1:nsd,:)
      DEALLOCATE(tmpX)

      ALLOCATE(lFa%IEN(lFa%eNoN,lFa%nEl))
      CALL getVTK_elemIEN(vtp, lFa%IEN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (ien)"

      ALLOCATE(lFa%gN(lFa%nNo))
      CALL getVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      IF (iStat .LT. 0) THEN
         DEALLOCATE(lFa%gN)
         wrn = " Could not find GlobalNodeID in the vtp file"
      ELSE
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)+1
               Ac = lFa%gN(Ac)
               lFa%IEN(a,e) = Ac
            END DO
         END DO
      END IF

      ALLOCATE(lFa%gE(lFa%nEl))
      CALL getVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      IF (iStat .LT. 0) THEN
         DEALLOCATE(lFa%gE)
         wrn = " Could not find GlobalElementID in the vtp file"
      ELSE
         lFa%gnEl = lFa%nEl
         ALLOCATE(lFa%gebc(1+lFa%eNoN,lFa%gnEl))
         lFa%gebc(1,:) = lFa%gE(:)
         lFa%gebc(2:1+lFa%eNoN,:) = lFa%IEN(:,:)
      END IF

      CALL flushVTK(vtp)

      RETURN
      END SUBROUTINE READVTP
!--------------------------------------------------------------------
      SUBROUTINE READVTUS(lA, lY, lD, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(tDof,gtnNo), lY(tDof,gtnNo),
     2   lD(tDof,gtnNo)
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, l, s, e, a, b, Ac,
     2   nNo, oGrp
      CHARACTER(LEN=stdL) :: varName
      TYPE(vtkXMLType) :: vtu

      REAL(KIND=RKIND), ALLOCATABLE :: tmpS(:,:), tmpGS(:,:)

      iStat = 0
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (init)"

      CALL getVTK_numPoints(vtu, nNo, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (num points)"
      IF (nNo .NE. SUM(msh(:)%gnNo))
     2   err = "Incompatible mesh and "//TRIM(fName)

      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutput
            IF (.NOT.eq(iEq)%output(iOut)%wtn(1)) CYCLE
            l = eq(iEq)%output(iOut)%l
            s = eq(iEq)%s + eq(iEq)%output(iOut)%o
            e = s + l - 1
            varName = TRIM(eq(iEq)%output(iOut)%name)

            oGrp = eq(iEq)%output(iOut)%grp
            SELECT CASE(oGrp)
            CASE (outGrp_A, outGrp_Y, outGrp_D)
               IF (l .GT. 1) THEN
                  ALLOCATE(tmpGS(maxNSD,nNo))
               ELSE
                  ALLOCATE(tmpGS(1,nNo))
               END IF
               CALL getVTK_pointData(vtu, TRIM(varName), tmpGS, iStat)
               IF (iStat .LT. 0) err ="VTU file read error (point data)"

               IF (nNo .NE. gtnNo) THEN
                  IF (l .GT. 1) THEN
                     ALLOCATE(tmpS(maxNSD,gtnNo))
                  ELSE
                     ALLOCATE(tmpS(1,gtnNo))
                  END IF
                  IF (.NOT.ALLOCATED(msh(1)%gpN)) err = "Unexpected "//
     2               "behavior. VTU file read error"
                  b = 0
                  DO iM=1, nMsh
                     DO a=1, msh(iM)%gnNo
                        b = b + 1
                        Ac = msh(iM)%gpN(a)
                        tmpS(:,Ac) = tmpGS(:,b)
                     END DO
                  END DO
                  DEALLOCATE(tmpGS)
                  ALLOCATE(tmpGS(SIZE(tmpS,1),gtnNo))
                  tmpGS = tmpS
                  DEALLOCATE(tmpS)
               END IF
            CASE DEFAULT
               CYCLE
            END SELECT

            SELECT CASE(oGrp)
            CASE (outGrp_A)
               lA(s:e,:) = tmpGS(1:l,:)
            CASE (outGrp_Y)
               lY(s:e,:) = tmpGS(1:l,:)
            CASE (outGrp_D)
               lD(s:e,:) = tmpGS(1:l,:)
            END SELECT
            DEALLOCATE(tmpGS)
         END DO
      END DO

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READVTUS
!--------------------------------------------------------------------
!     Read a particular point dataset from a vtu file. Point data will
!     be read and stored in lM%x array
      SUBROUTINE READVTUPDATA(lM, fName, kwrd, m, idx)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: m, idx
      CHARACTER(LEN=*) :: fName, kwrd

      INTEGER(KIND=IKIND) :: iStat, a
      TYPE(vtkXMLType) :: vtu

      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtu, fName, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (init)"

      CALL getVTK_numPoints(vtu, a, iStat)
      IF (a .NE. lM%gnNo) err = "Mismatch in num points for "//
     2   TRIM(kwrd)

      IF (m .EQ. nsd) THEN
         ALLOCATE(tmpR(maxNSD,lM%gnNo))
         tmpR = 0._RKIND
         CALL getVTK_pointData(vtu, TRIM(kwrd), tmpR, iStat)
         IF (iStat .LT. 0) err = "VTU file read error "//TRIM(kwrd)
         DO a=1, lM%gnNo
            lM%x((idx-1)*nsd+1:idx*nsd,a) = tmpR(1:nsd,a)
         END DO
         DEALLOCATE(tmpR)
      ELSE
         CALL getVTK_pointData(vtu, TRIM(kwrd), lM%x, iStat)
         IF (iStat .LT. 0) err = "VTU file read error "//TRIM(kwrd)
      END IF

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE READVTUPDATA
!--------------------------------------------------------------------
!     Read a particular point dataset from a vtp file. Point data will
!     be read and stored in lFa%x array
      SUBROUTINE READVTPPDATA(lFa, fName, kwrd, m, idx)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      INTEGER(KIND=IKIND), INTENT(IN) :: m, idx
      CHARACTER(LEN=*) :: fName, kwrd

      INTEGER(KIND=IKIND) :: iStat, a
      TYPE(vtkXMLType) :: vtp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpR(:,:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtp, fName, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (init)"

      CALL getVTK_numPoints(vtp, a, iStat)
      IF (a .NE. lFa%nNo) err = "Mismatch in num points for "//
     2   TRIM(kwrd)

      IF (m .EQ. nsd) THEN
         ALLOCATE(tmpR(maxNSD,lFa%nNo))
         tmpR = 0._RKIND
         CALL getVTK_pointData(vtp, TRIM(kwrd), tmpR, iStat)
         IF (iStat .LT. 0) err = "VTP file read error "//TRIM(kwrd)
         DO a=1, lFa%nNo
            lFa%x((idx-1)*nsd+1:idx*nsd,a) = tmpR(1:nsd,a)
         END DO
         DEALLOCATE(tmpR)
      ELSE
         CALL getVTK_pointData(vtp, TRIM(kwrd), lFa%x, iStat)
         IF (iStat .LT. 0) err = "VTU file read error "//TRIM(kwrd)
      END IF

      CALL flushVTK(vtp)

      RETURN
      END SUBROUTINE READVTPPDATA
!####################################################################
      SUBROUTINE WRITEVTU(lM, fName)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER(KIND=IKIND) :: iStat
      TYPE(vtkXMLType) :: vtu

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

      CALL putVTK_pointCoords(vtu, lM%x, iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

      CALL putVTK_elemIEN(vtu, lM%gIEN, lM%vtkType, iStat)
      IF (iStat .LT. 0) err = "VTU file write error (ien)"

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"
      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE WRITEVTU
!--------------------------------------------------------------------
      SUBROUTINE WRITEVTP(lFa, fName)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      INTEGER(KIND=IKIND) :: iStat, vtkType
      TYPE(vtkXMLType) :: vtp

      IF (nsd .EQ. 3) THEN
         SELECT CASE (lFa%eNoN)
         CASE (3)
            vtkTYpe = 5
         CASE (4)
            vtkType = 9
         END SELECT
      ELSE
         SELECT CASE (lFa%eNoN)
         CASE (2)
            vtkType = 3
         CASE (3)
            vtkType = 21
         END SELECT
      END IF

      CALL vtkInitWriter(vtp, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTP file write error (init)"

      CALL putVTK_pointCoords(vtp, lFa%x, iStat)
      IF (iStat .LT. 0) err = "VTP file write error (coords)"

      CALL putVTK_elemIEN(vtp, lFa%IEN, vtkType, iStat)
      IF (iStat .LT. 0) err = "VTP file write error (ien)"

      IF (ALLOCATED(lFa%gN)) THEN
         CALL putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
         IF (iStat .LT. 0) err = "VTP file write error (point data)"
      END IF

      IF (ALLOCATED(lFa%gE)) THEN
         CALL putVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
         IF (iStat .LT. 0) err = "VTP file write error (cell data)"
      END IF

      CALL vtkWriteToFile(vtp, iStat)
      IF (iStat .LT. 0) err = "VTP file write error"
      CALL flushVTK(vtp)

      RETURN
      END SUBROUTINE WRITEVTP
!--------------------------------------------------------------------
      SUBROUTINE WRITEVTUS(lA, lY, lD, lAve)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lA(tDof,tnNo), lY(tDof,tnNo),
     2   lD(tDof,tnNo)
      LOGICAL, INTENT(IN) :: lAve

      LOGICAL :: lIbl, lD0
      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo,
     2   nEl, s, l, ie, is, nSh, oGrp, outDof, nOut, cOut, ne, iFn, nFn,
     3   nOute
      CHARACTER(LEN=stdL) :: fName
      TYPE(dataType) :: d(nMsh)
      TYPE(vtkXMLType) :: vtu

      INTEGER(KIND=IKIND), ALLOCATABLE :: outS(:), tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:), tmpVe(:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:), outNamesE(:)

      lIbl = .FALSE.
      lD0  = .FALSE.

      a  = SUM(iblank(:))
      a  = cm%reduce(a)
      IF (a .GT. 0) lIbl = .TRUE.

      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_CMM .AND. ALLOCATED(Dinit)) THEN
            lD0 = .TRUE.
            EXIT
         END IF
      END DO

      nOut   = 1
      nOute  = 0
      outDof = nsd
      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutput
            IF (.NOT.eq(iEq)%output(iOut)%wtn(1)) CYCLE
            oGrp   = eq(iEq)%output(iOut)%grp
            IF (oGrp .EQ. outGrp_fN) THEN
               nFn = 1
               DO iM=1, nMsh
                  nFn = MAX(nFn, msh(iM)%nFn)
               END DO
               nOut   = nOut + nFn
               outDof = outDof + eq(iEq)%output(iOut)%l * nFn
            ELSE
               nOut   = nOut + 1
               outDof = outDof + eq(iEq)%output(iOut)%l
            END IF

            IF (oGrp.EQ.outGrp_J .OR. oGrp.EQ.outGrp_Mises)
     2         nOute = nOute + 1
         END DO
      END DO

!     iblank array for immersed bodies
      IF (lIbl) THEN
         nOut   = nOut + 1
         outDof = outDof + 1
      END IF

!     Initial displacements for CMM equation
      IF (lD0) THEN
         nOut = nOut + 1
         outDof = outDof + nsd
      END IF

      ALLOCATE(outNames(nOut), outS(nOut+1), outNamesE(nOute))

!     Prepare all solultions in to dataType d
      DO iM=1, nMsh
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE(d(iM)%x(outDof,msh(iM)%nNo))
         ALLOCATE(d(iM)%xe(nOute,msh(iM)%nEl))
         ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = x(:,Ac)/msh(iM)%scF
         END DO

         IF (lD0) THEN
            cOut           = cOut + 1
            is             = outS(cOut)
            ie             = is + nsd - 1
            outS(cOut+1)   = ie + 1
            outNames(cOut) = "Initial_displacement"
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               d(iM)%x(is:ie,a) = Dinit(:,Ac)
            END DO
         END IF

         nOute = 0
         outNamesE = ""

         DO iEq=1, nEq
            DO iOut=1, eq(iEq)%nOutput
               IF (.NOT.eq(iEq)%output(iOut)%wtn(1)) CYCLE
               l = eq(iEq)%output(iOut)%l
               s = eq(iEq)%s + eq(iEq)%output(iOut)%o
               e = s + l - 1

               cOut = cOut + 1
               is   = outS(cOut)
               ie   = is + l - 1
               outS(cOut+1) = ie + 1

               oGrp = eq(iEq)%output(iOut)%grp
               outNames(cOut) = TRIM(eq(iEq)%output(iOut)%name)

               SELECT CASE (oGrp)
                  CASE (outGrp_NA)
                  err = "Undefined output grp in VTK"

               CASE (outGrp_A)
                  DO a=1, msh(iM)%nNo
                     Ac = msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lA(s:e,Ac)
                  END DO

               CASE (outGrp_Y)
                  DO a=1, msh(iM)%nNo
                     Ac = msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lY(s:e,Ac)
                  END DO

               CASE (outGrp_D)
                  DO a=1, msh(iM)%nNo
                     Ac = msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lD(s:e,Ac)/msh(iM)%scF
                  END DO

               CASE (outGrp_WSS, outGrp_trac)
                  CALL BPOST(msh(iM), tmpV, lY, lD, oGrp)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is:ie,a) = tmpV(1:l,a)
                  END DO

               CASE (outGrp_vort, outGrp_eFlx, outGrp_hFlx,
     2            outGrp_stInv, outGrp_vortex, outGrp_Visc)
                  CALL POST(msh(iM), tmpV, lY, lD, oGrp, iEq)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is:ie,a) = tmpV(1:l,a)
                  END DO

               CASE (outGrp_absV)
                  DO a=1, msh(iM)%nNo
                     Ac = msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lY(1:nsd,Ac)
     2                                - lY(nsd+2:2*nsd+1,Ac)
                  END DO

               CASE (outGrp_fN)
                  cOut = cOut - 1
                  IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(nFn*nsd,msh(iM)%nNo))
                  tmpV = 0._RKIND
                  IF (msh(iM)%nFn .NE. 0)
     2               CALL FIBDIRPOST(msh(iM), nFn, tmpV, lD, iEq)
                  DO iFn=1, nFn
                     cOut = cOut + 1
                     is   = outS(cOut)
                     ie   = is + l - 1
                     outS(cOut+1)   = ie + 1
                     outNames(cOut) = TRIM(eq(iEq)%output(iOut)%name)//
     2                  STR(iFn)
                     DO a=1, msh(iM)%nNo
                        d(iM)%x(is:ie,a) =
     2                     tmpV((iFn-1)*nsd+1:iFn*nsd,a)
                     END DO
                  END DO
                  DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

               CASE (outGrp_fA)
                  IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(1,msh(iM)%nNo))
                  tmpV = 0._RKIND
                  IF (msh(iM)%nFn .EQ. 2)
     2               CALL FIBALGNPOST(msh(iM), tmpV, lD, iEq)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is,a) = tmpV(1,a)
                  END DO
                  DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

               CASE (outGrp_stress, outGrp_cauchy, outGrp_mises)
                  IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(l,msh(iM)%nNo), tmpVe(msh(iM)%nEl))
                  tmpV  = 0._RKIND
                  tmpVe = 0._RKIND
                  IF (pstEq) THEN
                     DO a=1, msh(iM)%nNo
                        Ac = msh(iM)%gN(a)
                        tmpV(:,a) = pS0(:,Ac)
                     END DO
                  END IF

                  IF (.NOT.cmmInit) CALL TPOST(msh(iM), l, tmpV, tmpVe,
     2               lD, lY, iEq, oGrp)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is:ie,a) = tmpV(1:l,a)
                  END DO

                  IF (oGrp .EQ. outGrp_mises) THEN
                     nOute = nOute + 1
                     outNamesE(nOute) = "E_VonMises"
                     DO a=1, msh(iM)%nEl
                        d(iM)%xe(nOute,a) = tmpVe(a)
                     END DO
                  END IF

                  DEALLOCATE(tmpV, tmpVe)
                  ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

               CASE (outGrp_J, outGrp_F, outGrp_strain)
                  IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(l,msh(iM)%nNo), tmpVe(msh(iM)%nEl))
                  tmpV  = 0._RKIND
                  tmpVe = 0._RKIND

                  CALL TPOST(msh(iM), l, tmpV, tmpVe, lD, lY, iEq, oGrp)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is:ie,a) = tmpV(1:l,a)
                  END DO

                  IF (oGrp .EQ. outGrp_J) THEN
                     nOute = nOute + 1
                     outNamesE(nOute) = "E_Jacobian"
                     DO a=1, msh(iM)%nEl
                        d(iM)%xe(nOute,a) = tmpVe(a)
                     END DO
                  END IF

                  DEALLOCATE(tmpV, tmpVe)
                  ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))

               CASE (outGrp_divV)
                  IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(1,msh(iM)%nNo))
                  tmpV = 0._RKIND
                  CALL DIVPOST(msh(iM), tmpV, lY, lD, iEq)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is,a) = tmpV(1,a)
                  END DO
                  DEALLOCATE(tmpV)
                  ALLOCATE(tmpV(maxnsd,msh(iM)%nNo))
               CASE DEFAULT
                  err = "Undefined output"
               END SELECT
            END DO
         END DO

         IF (lIbl) THEN
            cOut = cOut + 1
            is   = outS(cOut)
            ie   = is
            outS(cOut+1)   = ie + 1
            outNames(cOut) = "IBLANK"
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               d(iM)%x(is:ie,a) = REAL(iblank(Ac), KIND=RKIND)
            END DO
         END IF
      END DO

!     Integrate data from all processors
      nNo = 0
      nEl = 0
      DO iM=1, nMsh
         IF (msh(iM)%eType .EQ. eType_NRB) THEN
            CALL INTNRBDATA(msh(iM), d(iM), outDof, nOute)
         ELSE
            CALL INTMSHDATA(msh(iM), d(iM), outDof, nOute)
         END IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
      END DO
      IF (cm%slv()) THEN
         savedOnce = .TRUE.
         RETURN
      END IF

      DEALLOCATE(tmpV)
      ALLOCATE(tmpV(maxnsd, nNo))

!     Writing to vtu file (master only)
      IF (cTS.GE.1000 .OR. lAve) THEN
         fName = STR(cTS)
      ELSE
         WRITE(fName,'(I3.3)') cTS
      END IF

      fName = TRIM(saveName)//"_"//TRIM(ADJUSTL(fName))//".vtu"
      dbg = "Writing VTU"

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

!     Writing the position data
      iOut = 1
      s   = outS(iOut)
      e   = outS(iOut+1)-1
      nSh  = 0
      tmpV = 0._RKIND
      DO iM=1, nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%gx(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO
      CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

!     Writing the connectivity data
      nSh = -1
      DO iM=1, nMsh
         ALLOCATE(tmpI(d(iM)%eNoN, d(iM)%nEl))
         DO e=1, d(iM)%nEl
            tmpI(:,e) = d(iM)%IEN(:,e) + nSh
         END DO
         CALL putVTK_elemIEN(vtu, tmpI, d(iM)%vtkType, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (ien)"
         DEALLOCATE(tmpI)
         nSh = nSh + d(iM)%nNo
      END DO

!     Writing all solutions
      DO iOut=2, nOut
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         s = outS(iOut)
         e = outS(iOut+1) - 1
         l = e - s + 1
         ALLOCATE(tmpV(l, nNo))
         nSh = 0
         DO iM=1, nMsh
            DO a=1, d(iM)%nNo
               tmpV(:,a+nSh) = d(iM)%gx(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (point data)"
      END DO

!     Write element-based variables
      ne = 0
      IF (.NOT.savedOnce .OR. nMsh.GT.1) THEN
         ALLOCATE(tmpI(1,nEl))
!     Write the domain ID
         ne = 1
         IF (ALLOCATED(dmnID)) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = INT(d(iM)%xe(e,ne), KIND=IKIND)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (dom id)"
         END IF

         IF (.NOT.savedOnce) THEN
            savedOnce = .TRUE.
            ne = ne + 1

!        Write partition data
            IF (.NOT.cm%seq()) THEN
               Ec = 0
               DO iM=1, nMsh
                  DO e=1, d(iM)%nEl
                     Ec = Ec + 1
                     tmpI(1,Ec) = INT(d(iM)%xe(e,ne), KIND=IKIND)
                  END DO
               END DO
               CALL putVTK_elemData(vtu, 'Proc_ID', tmpI, iStat)
               IF (iStat .LT. 0) err = "VTU file write error (proc id)"
            END IF
         END IF

!     Write the mesh ID
         IF (nMsh .GT. 1) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = iM
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Mesh_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (mesh id)"
         END IF
         DEALLOCATE(tmpI)
      END IF

!     Write element Jacobian and von Mises stress if necessary
      DO l=1, nOute
         ne = ne + 1
         ALLOCATE(tmpVe(nEl))
         tmpVe = 0._RKIND
         Ec = 0
         DO iM=1, nMsh
            IF (ALLOCATED(d(iM)%xe)) THEN
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpVe(Ec) = d(iM)%xe(e,ne)
               END DO
            END IF
         END DO
         CALL putVTK_elemData(vtu, outNamesE(l), tmpVe, iStat)
         IF (iStat .LT. 0) err = "VTU file write error ("//
     2      TRIM(outNamesE(l))//")"
         DEALLOCATE(tmpVe)
      END DO

!     Write element ghost cells if necessary
      IF (lIbl) THEN
         ne = ne + 1
         ALLOCATE(tmpI(1,nEl))
         tmpI = 0
         Ec = 0
         DO iM=1, nMsh
            IF (ALLOCATED(d(iM)%xe)) THEN
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = INT(d(iM)%xe(e,ne), KIND=IKIND)
               END DO
            END IF
         END DO
         CALL putVTK_elemData(vtu, 'EGHOST', tmpI, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (EGHOST)"
         DEALLOCATE(tmpI)
      END IF

      DO iM=1, nMsh
         CALL DESTROY(d(iM))
      END DO

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE WRITEVTUS
!####################################################################
!     This routine prepares data array of a regular mesh
      SUBROUTINE INTMSHDATA(lM, d, outDof, nOute)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(dataType), INTENT(INOUT) :: d
      INTEGER(KIND=IKIND), INTENT(IN) :: outDof, nOute

      INTEGER(KIND=IKIND) :: e, i, ierr, Ec, m

      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disps(:), tmpI(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lDe(:,:), gDe(:,:)

      d%eNoN    = lM%eNoN
      d%vtkType = lM%vtkType
      IF (cm%mas()) THEN
         d%nNo  = lM%gnNo
         d%nEl  = lM%gnEl
         d%eNoN = lM%eNoN
      ELSE
         d%nNo  = 0
         d%nEl  = 0
      END IF
      ALLOCATE(d%IEN(d%eNoN,d%nEl), d%gx(outDof,d%nNo))

      d%gx = GLOBAL(lM, d%x)
      DEALLOCATE(d%x)

      IF (cm%mas()) THEN
         DO e=1, lM%gnEl
            Ec = lM%otnIEN(e)
            d%IEN(:,e) = lM%gIEN(:,Ec)
         END DO
      END IF

!     Element variables
      m = nOute
      IF (.NOT.savedOnce .OR. nMsh.GT.1) THEN
         IF (savedOnce) THEN
            m = m + 1
         ELSE
            m = m + 2
         END IF
      END IF
      IF (ALLOCATED(d%xe)) THEN
         ALLOCATE(lDe(nOute, lM%nEl))
         DO e=1, lM%nEl
            DO i=1, nOute
               lDe(i,e) = d%xe(i,e)
            END DO
         END DO
         DEALLOCATE(d%xe)
      END IF
      IF (ALLOCATED(lM%iGC)) m = m + 1
      IF (m .NE. 0) ALLOCATE(d%xe(d%nEl,m))

      m = 0
      IF (.NOT.savedOnce .OR. nMsh.GT.1) THEN
         m = 1
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sCount(cm%np()), tmpI(d%nEl))
            DO i=1, cm%np()
               sCount(i) = lM%eDist(i) - lM%eDist(i-1)
            END DO
            CALL MPI_GATHERV(lM%eId, lM%nEl, mpint, tmpI, sCount,
     2         lM%eDist, mpint, master, cm%com(), ierr)

            IF (cm%mas()) THEN
               DO e=1, lM%gnEl
                  Ec = lM%otnIEN(e)
                  d%xe(e,m) = REAL(tmpI(Ec), KIND=RKIND)
               END DO
            END IF
            DEALLOCATE(sCount, tmpI)
         ELSE
            d%xe(:,m) = 1._RKIND
         END IF

         IF (.NOT.savedOnce) THEN
            m = m + 1
            ALLOCATE(tmpI(d%nEl))
            IF (cm%mas()) THEN
               IF (.NOT.cm%seq()) THEN
                  i = 0
                  DO e=1, lM%gnEl
                     IF (e .GT. lM%eDist(i)) THEN
                        DO
                           i = i + 1
                           IF (e .LE. lM%eDist(i)) EXIT
                        END DO
                     END IF
                     tmpI(e) = i
                  END DO
                  ALLOCATE(sCount(lM%gnEl))
                  sCount = tmpI(:)
                  DO e=1, lM%gnEl
                     Ec = lM%otnIEN(e)
                     d%xe(e,m) = REAL(sCount(Ec), KIND=RKIND)
                  END DO
                  DEALLOCATE(sCount)
               END IF
            END IF
            DEALLOCATE(tmpI)
         END IF
      END IF

      IF (ALLOCATED(lDe)) THEN
         IF (cm%mas()) THEN
            ALLOCATE(disps(cm%np()), sCount(cm%np()), gDe(nOute,d%nEl))
            DO i=1, cm%np()
               disps(i)  = lM%eDist(i-1)*nOute
               sCount(i) = lM%eDist(i)*nOute - disps(i)
            END DO
         ELSE
            ALLOCATE(disps(0), sCount(0), gDe(0,0))
         END IF
         CALL MPI_GATHERV(lDe, nOute*lM%nEl, mpreal, gDe, sCount, disps,
     2      mpreal, master, cm%com(), ierr)

         IF (cm%mas()) THEN
            DO e=1, lM%gnEl
               Ec = lM%otnIEN(e)
               DO i=1, nOute
                  d%xe(e,m+i) = gDe(i,Ec)
               END DO
            END DO
         END IF
         DEALLOCATE(lDe, gDe, sCount, disps)
         m = m + nOute
      END IF

      IF (ALLOCATED(lM%iGC)) THEN
         m = m + 1
         ALLOCATE(sCount(cm%np()), tmpI(d%nEl))
         DO i=1, cm%np()
            sCount(i) = lM%eDist(i) - lM%eDist(i-1)
         END DO
         CALL MPI_GATHERV(lM%iGC, lM%nEl, mpint, tmpI, sCount,
     2      lM%eDist, mpint, master, cm%com(), ierr)

         IF (cm%mas()) THEN
            DO e=1, lM%gnEl
               Ec = lM%otnIEN(e)
               d%xe(e,m) = REAL(tmpI(Ec), KIND=RKIND)
            END DO
         END IF
         DEALLOCATE(sCount, tmpI)
      END IF

      RETURN
      END SUBROUTINE INTMSHDATA
!--------------------------------------------------------------------
!     This routine will interpolates NURBS into a QUD4/HEX8 mesh, so you
!     can write it into a VTK file.
!     TODO: Element data are not written to VTK file for NURBS
      SUBROUTINE INTNRBDATA(lM, d, outDof, nOute)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(dataType), INTENT(INOUT) :: d
      INTEGER(KIND=IKIND), INTENT(IN) :: outDof, nOute

      LOGICAL flag, clcDmn
      INTEGER(KIND=IKIND) ie, e, ex, ey, ez, s, sx, sy, sz, ia, Ac, a,
     2   ax, ay, az, jx, jy, jz, ierr, i, nShft, insd, m

      INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:), tDmnId(:), sCe(:),
     2   dise(:), sCn(:), disn(:)
      REAL(KIND=RKIND), ALLOCATABLE :: N(:,:), tmpX(:,:), tmpI(:),
     2   lDe(:,:), gDe(:,:)

      insd = nsd
      IF (lM%lShl) insd = nsd - 1

      jy = lM%bs(2)%nEl*(lM%bs(2)%nSl - 1) + 1
      IF (insd .EQ. 2) THEN
         jz     = 0
         d%eNoN = 4
      ELSE IF (insd .EQ. 3) THEN
         jz     = lM%bs(3)%nEl*(lM%bs(3)%nSl - 1) + 1
         d%eNoN = 8
      END IF
      ALLOCATE(sCe(cm%np()), dise(0:cm%np()), sCn(cm%np()),
     2   disn(0:cm%np()))

!     Preparing sCn, sCe for the future communications
      dise = 0
      disn = 0
      DO i=1, cm%np()
         e = lM%bs(1)%nEl*(lM%eDist(i) - lM%eDist(i-1))/lM%gnEl
         jx = e*(lM%bs(1)%nSl - 1) + 1
         IF (insd .EQ. 2) jz = 2
         sCe(i) = (jx - 1)*(jy - 1)*(jz - 1)

         IF (insd .EQ. 2) jz = 1
         IF (lM%eDist(i).NE.lM%gnEl .OR. e.EQ.0) jx = jx - 1
         sCn(i) = jx*jy*jz

         dise(i) = dise(i-1) + sCe(i)
         disn(i) = disn(i-1) + sCn(i)
      END DO
      IF (cm%mas()) THEN
         d%nEl = dise(cm%np())
         d%nNo = disn(cm%np())
      ELSE
         d%nEl = 0
         d%nNo = 0
      END IF

      ALLOCATE(d%IEN(d%eNoN,d%nEl), d%gx(outDof,d%nNo),
     2   IEN(d%eNoN,sCe(cm%tF())), tmpX(outDof,sCn(cm%tF())),
     3   N(lM%eNoN,lM%nSl))

!     I am calculating the last column only if this is the last
!     partition that has a section of this mesh. This is to prevent
!     repetition associated with the node on the boundary of partitions
      flag = .FALSE.
      IF (lM%eDist(cm%tF()).NE.lM%gnEl .OR. sCe(cm%tF()).EQ.0) THEN
         flag = .TRUE.
      END IF
!     This would be the number of nodes in the lower procesors
      nShft = disn(cm%id())

!     Element variables
      clcDmn = .FALSE.
      m = nOute
      IF (.NOT.savedOnce .OR. nMsh.GT.1) THEN
         IF (savedOnce) THEN
            m = m + 1
         ELSE
            m = m + 2
         END IF
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(tDmnId(sCe(cm%tF())))
            IF (ALLOCATED(lM%eId)) THEN
               clcDmn = .TRUE.
            ELSE
               tDmnId = 1
            END IF
         END IF
      END IF

      IF (ALLOCATED(d%xe)) THEN
         ALLOCATE(lDe(nOute,lM%nEl))
         DO e=1, lM%nEl
            DO i=1, nOute
               lDe(i,e) = d%xe(i,e)
            END DO
         END DO
         DEALLOCATE(d%xe)
      END IF
      IF (m .NE. 0) ALLOCATE(d%xe(d%nEl,m))

      ie = 0
      IF (insd .EQ. 2) THEN
!     I am using d%eType = eType_BIL
         d%vtkType = 9
         DO e=1, lM%nEl
            ex =    (e-1)/lM%bs(2)%nEl  + 1
            ey = MOD(e-1, lM%bs(2)%nEl) + 1
            ax = (ex - 1)*(lM%bs(1)%nSl - 1)
            DO sx=1, lM%bs(1)%nSl-1
               ax = ax + 1
               ay = (ey - 1)*(lM%bs(2)%nSl - 1)
               DO sy=1, lM%bs(2)%nSl-1
                  ie = ie + 1
                  ay = ay + 1
                  ia = (ax - 1)*jy + ay + nShft

                  IEN(1,ie) = ia + jy + 1
                  IEN(2,ie) = ia      + 1
                  IEN(3,ie) = ia
                  IEN(4,ie) = ia + jy
                  IF (clcDmn) tDmnId(ie) = lM%eId(e)
               END DO
            END DO

            s  = 0
            CALL NRBNNS(lM, N, e)
            ax = (ex - 1)*(lM%bs(1)%nSl - 1)
            DO sx=1, lM%bs(1)%nSl
               IF (flag .AND. sx.EQ.lM%bs(1)%nSl) CYCLE
               ax = ax + 1
               ay = (ey - 1)*(lM%bs(2)%nSl - 1)
               DO sy=1, lM%bs(2)%nSl
                  s  = s + 1
                  ay = ay + 1
                  ia = (ax - 1)*jy + ay
                  tmpX(:,ia) = 0._RKIND
                  DO a=1, lM%eNoN
                     Ac = lM%IEN(a,e)
                     Ac = lM%lN(Ac)
                     tmpX(:,ia) = tmpX(:,ia) + d%x(:,Ac)*N(a,s)
                  END DO
               END DO
            END DO
         END DO
      ELSE
!     I am using d%eType = eType_BRK
         d%vtkType = 12
         DO e=1, lM%nEl
            ex =   (e-1)/(lM%bs(2)%nEl*lM%bs(3)%nEl) + 1
            ey = (MOD(e-1,lM%bs(2)%nEl*lM%bs(3)%nEl))/lM%bs(3)%nEl + 1
            ez = MOD(e-1, lM%bs(3)%nEl) + 1
            ax = (ex - 1)*(lM%bs(1)%nSl - 1)
            DO sx=1, lM%bs(1)%nSl-1
               ax = ax + 1
               ay = (ey - 1)*(lM%bs(2)%nSl - 1)
               DO sy=1, lM%bs(2)%nSl-1
                  ay = ay + 1
                  az = (ez - 1)*(lM%bs(3)%nSl - 1)
                  DO sz=1, lM%bs(3)%nSl-1
                     az = az + 1
                     ie = ie + 1
                     ia = ((ax - 1)*jy + ay - 1)*jz + az + nShft

                     IEN(1,ie) = ia + jy*jz + jz + 1
                     IEN(2,ie) = ia         + jz + 1
                     IEN(3,ie) = ia         + jz
                     IEN(4,ie) = ia + jy*jz + jz
                     IEN(5,ie) = ia + jy*jz      + 1
                     IEN(6,ie) = ia              + 1
                     IEN(7,ie) = ia
                     IEN(8,ie) = ia + jy*jz
                     IF (clcDmn) tDmnId(ie) = lM%eId(e)
                  END DO
               END DO
            END DO

            s  = 0
            CALL NRBNNS(lM, N, e)
            ax = (ex - 1)*(lM%bs(1)%nSl - 1)
            DO sx=1, lM%bs(1)%nSl
               IF (flag .AND. sx.EQ.lM%bs(1)%nSl) CYCLE
               ax = ax + 1
               ay = (ey - 1)*(lM%bs(2)%nSl - 1)
               DO sy=1, lM%bs(2)%nSl
                  ay = ay + 1
                  az = (ez - 1)*(lM%bs(3)%nSl - 1)
                  DO sz=1, lM%bs(3)%nSl
                     az = az + 1
                     s  = s + 1
                     ia = ((ax - 1)*jy + ay - 1)*jz + az
                     tmpX(:,ia) = 0._RKIND
                     DO a=1, lM%eNoN
                        Ac = lM%IEN(a,e)
                        Ac = lM%lN(Ac)
                        tmpX(:,ia) = tmpX(:,ia) + d%x(:,Ac)*N(a,s)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END IF
      DEALLOCATE(d%x)

      IF (cm%seq()) THEN
         d%gx  = tmpX
         d%IEN = IEN
      ELSE
!     First collecting all the elements connectivity
         sCe  = sCe*d%eNoN
         sCn  = sCn*outDof
         dise = dise*d%eNoN
         disn = disn*outDof
         CALL MPI_GATHERV(IEN, sCe(cm%tF()), mpint, d%IEN, sCe, dise,
     2      mpint, master, cm%com(), ierr)
!     Now collecting the solutions
         CALL MPI_GATHERV(tmpX, sCn(cm%tF()), mpreal, d%gx, sCn, disn,
     2      mpreal, master, cm%com(), ierr)
         sCe  = sCe/d%eNoN
         dise = dise/d%eNoN
      END IF

!     Default element variables (Domain_ID, Proc_ID)
      m = 0
      IF (.NOT.savedOnce .OR. nMsh.GT.1) THEN
!        Domain ID
         m = 1
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(tmpI(d%nEl))
            CALL MPI_GATHERV(tDmnId, sCe(cm%tF()), mpint, tmpI, sCe,
     2         dise, mpint, master, cm%com(), ierr)
            IF (cm%mas()) THEN
               DO e=1, d%nEl
                  d%xe(e,1) = REAL(tmpI(e), KIND=RKIND)
               END DO
            END IF
            DEALLOCATE(tmpI)
         ELSE
            d%xe(:,m) = 1._RKIND
         END IF

         IF (.NOT.savedOnce) THEN
!           Proc_ID for parallel run
            m = m + 1
            IF (cm%mas()) THEN
               IF (.NOT.cm%seq()) THEN
                  i = 0
                  DO e=1, d%nEl
                     IF (e .GT. dise(i)) THEN
                        DO
                           i = i + 1
                           IF (e .LE. dise(i)) EXIT
                        END DO
                     END IF
                     d%xe(e,m) = REAL(i, KIND=RKIND)
                  END DO
               END IF
            END IF
         END IF
      END IF

      IF (ALLOCATED(lDe)) THEN
         sCe = sCe*nOute
         dise = dise*nOute
         IF (cm%mas()) THEN
            ALLOCATE(gDe(nOute,d%nEl))
         ELSE
            ALLOCATE(gDe(0,0))
         END IF

         CALL MPI_GATHERV(lDe, sCe(cm%tF()), mpreal, gDe, sCe, dise,
     2      mpreal, master, cm%com(), ierr)

         IF (cm%mas()) THEN
            DO e=1, d%nEl
               DO i=1, nOute
                  d%xe(i,m+1) = gDe(i,e)
               END DO
            END DO
         END IF
         DEALLOCATE(lDe, gDe)
      END IF

      RETURN
      END SUBROUTINE INTNRBDATA
!####################################################################
!     This routine calculates average of VTU files.
      SUBROUTINE CALCAVE
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      LOGICAL :: lExist
      INTEGER(KIND=IKIND) :: fid, n, lTS, fTS, iTS
      REAL(KIND=RKIND) :: cntr
      CHARACTER(LEN=stdL) fName

      REAL(KIND=RKIND), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:),
     2   Ag(:,:), Yg(:,:), Dg(:,:)

      fid  = 1
      fTS  = nTS
      lTS  = 0
      cntr = 0._RKIND
      iTS  = 0
      IF (zeroAve) iTS = rsTS
      IF (cm%mas()) THEN
         ALLOCATE(tmpA(tDof,gtnNo), tmpY(tDof,gtnNo),
     2      tmpD(tDof,gtnNo))
      ELSE
         ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
      END IF
      ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo))

      Ag   = 0._RKIND
      Yg   = 0._RKIND
      Dg   = 0._RKIND
      std = " Computing average quantities"
      DO n=iTS, nTS, saveIncr
         IF (n .LT. saveATS) CYCLE
         IF (n .LT. fTS) fTS = n
         IF (n .GT. lTS) lTS = n
         cntr = cntr + 1._RKIND
         IF (cm%mas()) THEN
            IF (n .GE. 1000) THEN
               fName = STR(n)
            ELSE
               WRITE(fName,'(I3.3)') n
            END IF
            fName = TRIM(saveName)//"_"//TRIM(ADJUSTL(fName))//".vtu"
            INQUIRE (FILE=TRIM(fName), EXIST=lExist)
            IF (lExist) THEN
               CALL READVTUS(tmpA, tmpY, tmpD, fName)
            ELSE
               err = "File "//TRIM(fName)//" is missing"
            END IF
         END IF
         Ag = Ag + LOCAL(tmpA)
         Yg = Yg + LOCAL(tmpY)
         Dg = Dg + LOCAL(tmpD)
      END DO
      n = NINT(cntr, KIND=IKIND)
      IF (n .EQ. 0) RETURN
      Ag = Ag/cntr
      Yg = Yg/cntr
      Dg = Dg/cntr
      cTS = lTS
      fName = saveName
      saveName = TRIM(saveName)//"_ave_"//STR(fTS)//"_"//
     2   STR(saveIncr)
      CALL WRITEVTUS(Ag, Yg, Dg, .TRUE.)
      saveName = fName

      RETURN
      END SUBROUTINE CALCAVE
!####################################################################

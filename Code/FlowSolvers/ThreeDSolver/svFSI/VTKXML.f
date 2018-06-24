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

!**********************************************************************

      SUBROUTINE READVTU(lM, fName)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod

      IMPLICIT NONE

      TYPE(MSHTYPE), INTENT(INOUT) :: lM
      CHARACTER(LEN=STDL) :: fName

      TYPE(vtkXMLType) :: vtu
      INTEGER :: iStat
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: tmpX

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

      CALL SELECTELE(lM)

      ALLOCATE(x(nsd,lM%gnNo),tmpX(maxNSD,lM%gnNo))
      CALL getVTK_pointCoords(vtu, tmpX, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (coords)"
      x(:,:) = tmpX(1:nsd,:)
      deALLOCATE(tmpX)
      
      ALLOCATE(lM%gIEN(lM%eNoN,lM%gnEl))
      CALL getVTK_elemIEN(vtu, lM%gIEN, iStat)
      IF (iStat .LT. 0) err = "VTU file read error (ien)"
      lM%gIEN = lM%gIEN + 1
      
      IF (ichckIEN) CALL CHECKIEN(lM, .FALSE.)

      CALL flushVTK(vtu)
      
      RETURN
      END SUBROUTINE READVTU

!**********************************************************************

      SUBROUTINE READVTP(lM, lFa, fName)

      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod

      IMPLICIT NONE

      TYPE(MSHTYPE), INTENT(INOUT) :: lM
      TYPE(FACETYPE), INTENT(INOUT) :: lFa
      CHARACTER(LEN=STDL) :: fName

      TYPE(vtkXMLType) :: vtp
      INTEGER :: iStat, a, e, Ac

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

      CALL SELECTELEB(lM, lFa)

      ALLOCATE(lFa%IEN(lFa%eNoN,lFa%nEl))
      CALL getVTK_elemIEN(vtp, lFa%IEN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (ien)"

      ALLOCATE(lFa%gN(lFa%nNo))
      CALL getVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (point data)"

      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)+1
            Ac = lFa%gN(Ac)
            lFa%IEN(a,e) = Ac
         END DO
      END DO

      ALLOCATE(lFa%gE(lFa%nEl))
      CALL getVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (cell data)"

      CALL flushVTK(vtp)
      
      lFa%gnEl = lFa%nEl
      ALLOCATE(lFa%gebc(1+lFa%eNoN,lFa%gnEl))
      lFa%gebc(1,:) = lFa%gE(:)
      lFa%gebc(2:1+lFa%eNoN,:) = lFa%IEN(:,:)
      
      RETURN
      END SUBROUTINE READVTP

!**********************************************************************

      SUBROUTINE WRITEVTU(lM, fName)
      
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      
      IMPLICIT NONE
      
      TYPE(mshType), INTENT(IN) :: lM
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      
      TYPE(vtkXMLType) :: vtu
      INTEGER :: iStat
      
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
      
!**********************************************************************

      SUBROUTINE WRITEVTP(lFa, fName)
      
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      
      IMPLICIT NONE
      
      TYPE(faceType), INTENT(IN) :: lFa
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      
      TYPE(vtkXMLType) :: vtp
      INTEGER :: iStat, vtkType
      
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
      
      CALL putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
      IF (iStat .LT. 0) err = "VTP file write error (point data)"
      
      CALL putVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
      IF (iStat .LT. 0) err = "VTP file write error (cell data)"
      
      CALL vtkWriteToFile(vtp, iStat)
      IF (iStat .LT. 0) err = "VTP file write error"
      CALL flushVTK(vtp)
      
      RETURN
      END SUBROUTINE WRITEVTP
      
!**********************************************************************

      SUBROUTINE WRITEVTUS(lA, lY, lD)
      
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(IN) :: lA(tDof, tnNo), lY(tDof, tnNo),
     2   lD(tDof, tnNo)
      
      TYPE(dataType) :: d(nMsh)
      TYPE(vtkXMLType) :: vtu
      
      CHARACTER(LEN=stdL) :: fName
      INTEGER :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo, nEl
      INTEGER :: s, l, ie, is, nSh, oGrp, outDof, nOut, cOut
      
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)
      INTEGER, ALLOCATABLE :: outS(:), tmpI(:,:)
      REAL(KIND=8), ALLOCATABLE :: tmpV(:,:)
      
      nOut   = 1
      outDof = nsd
      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutput
            IF (.NOT.eq(iEq)%output(iOut)%wtn(1)) CYCLE
            nOut   = nOut + 1
            outDof = outDof + eq(iEq)%output(iOut)%l
         END DO
      END DO
      ALLOCATE(outNames(nOut), outS(nOut+1))
      
!     Prepare all solultions in to dataType d
      DO iM=1, nMsh
         IF (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
         ALLOCATE(d(iM)%x(outDof,msh(iM)%nNo), tmpV(maxnsd,msh(iM)%nNo))
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = x(:,Ac)
         END DO

         DO iEq=1, nEq
            DO iOut=1, eq(iEq)%nOutput
               IF (.NOT.eq(iEq)%output(iOut)%wtn(1)) CYCLE
               l  = eq(iEq)%output(iOut)%l
               s  = eq(iEq)%s + eq(iEq)%output(iOut)%o
               e  = s + l - 1

               cOut = cOut + 1
               is   = outS(cOut)
               ie   = is + l - 1
               outS(cOut+1)   = ie + 1
               outNames(cOut) = eq(iEq)%sym//"_"//
     2            TRIM(eq(iEq)%output(iOut)%name)

               oGrp = eq(iEq)%output(iOut)%grp
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
                     d(iM)%x(is:ie,a) = lD(s:e,Ac)
                  END DO
               CASE (outGrp_WSS)
                  CALL BPOST(msh(iM), tmpV, lY, lD, oGrp)
                  DO a=1, msh(iM)%nNo
                     d(iM)%x(is:ie,a) = tmpV(1:l,a) 
                  END DO
               CASE (outGrp_vort, outGrp_eFlx, outGrp_hFlx, 
     2            outGrp_stInv, outGrp_vortex)
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
               CASE DEFAULT
                  err = "Undefined output"
               END SELECT
            END DO
         END DO
      END DO
      
!     Integrate data from all processors
      nNo = 0
      nEl = 0
      DO iM=1, nMsh
         IF (msh(iM)%eType .EQ. eType_NRB) THEN
            CALL INTNRBDATA(msh(iM), d(iM), outDof)
         ELSE
            CALL INTMSHDATA(msh(iM), d(iM), outDof)
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
      IF (cTS .GE. 1000) THEN
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
      tmpV = 0D0
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
      IF (.NOT.savedOnce .OR. mvMsh) THEN
         savedOnce = .TRUE.
         ALLOCATE(tmpI(1,nEl))
         
!     Write partition data
         IF (.NOT.cm%seq()) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = d(iM)%xe(e,1)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Proc_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (cell data)"
         END IF

!     Write the domain ID
         IF (ALLOCATED(dmnID)) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = d(iM)%xe(e,2)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (cell data)"
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
            IF (iStat .LT. 0) err = "VTU file write error (cell data)"
         END IF
         DEALLOCATE(tmpI)
      END IF
      
      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"
      
      CALL flushVTK(vtu)
      
      RETURN
      END SUBROUTINE WRITEVTUS
      
!**********************************************************************
!     This routine prepares data array of a regular mesh
      SUBROUTINE INTMSHDATA(lM, d, outDof)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(dataType), INTENT(INOUT) :: d
      INTEGER, INTENT(IN) :: outDof

      INTEGER e, i, ierr, Ec
      INTEGER, ALLOCATABLE :: sCount(:)

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

      d%gx = GLOBAL(lM,d%x)
      DEALLOCATE(d%x)
      IF (cm%mas()) THEN
         DO e=1, lM%gnEl
            Ec = lM%otnIEN(e)
            d%IEN(:,e) = lM%gIEN(:,Ec)
         END DO
      END IF

      IF (.NOT.savedOnce .OR. mvMsh) THEN
         ALLOCATE(d%xe(d%nEl,2))
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
                  d%xe(e,1) = i
               END DO
               ALLOCATE(sCount(lM%gnEl))
               sCount = d%xe(:,1)
               DO e=1, lM%gnEl
                  Ec = lM%otnIEN(e)
                  d%xe(e,1) = sCount(Ec)
               END DO
               DEALLOCATE(sCount)
            END IF
         END IF
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sCount(cm%np()))
            DO i=1, cm%np()
               sCount(i) = lM%eDist(i) - lM%eDist(i-1)
            END DO
            CALL MPI_GATHERV(lM%eId, lM%nEl, mpint, d%xe(:,2), sCount, 
     2         lM%eDist, mpint, master, cm%com(), ierr)
            
            IF (cm%mas()) THEN
               DEALLOCATE(sCount)
               ALLOCATE(sCount(lM%gnEl))
               sCount = d%xe(:,2)
               DO e=1, lM%gnEl
                  Ec = lM%otnIEN(e)
                  d%xe(e,2) = sCount(Ec)
               END DO
               DEALLOCATE(sCount)
            END IF
         ELSE
            d%xe(:,2) = 1
         END IF
      END IF
      
      RETURN
      END SUBROUTINE INTMSHDATA

!**********************************************************************
!     This routine will interpolates NURBS into a BIL/BRK mesh, so you
!     can write it into a VTK file.
      SUBROUTINE INTNRBDATA(lM, d, outDof)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(dataType), INTENT(INOUT) :: d
      INTEGER, INTENT(IN) :: outDof

      LOGICAL flag, clcDmn
      INTEGER ie, e, ex, ey, ez, s, sx, sy, sz, ia, Ac, a, ax, ay, az,
     2   jx, jy, jz, ierr, i, nShft
      INTEGER, ALLOCATABLE :: IEN(:,:), tDmnId(:), sCe(:), dise(:),
     2   sCn(:), disn(:)
      REAL(KIND=8), ALLOCATABLE :: N(:,:), tmpX(:,:)

      jy = lM%bs(2)%nEl*(lM%bs(2)%nSl - 1) + 1
      IF (nsd .EQ. 2) THEN
         jz     = 0
         d%eNoN = 4
      ELSE
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
         IF (nsd .EQ. 2) jz = 2
         sCe(i) = (jx - 1)*(jy - 1)*(jz - 1)

         IF (nsd .EQ. 2) jz = 1
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

      clcDmn = .FALSE.
      IF (.NOT.savedOnce) THEN
         ALLOCATE(d%xe(d%nEl,2))
         IF (cm%mas()) THEN
!     Based on the two posible fields
            IF (.NOT.cm%seq()) THEN
               i = 0
               DO e=1, d%nEl
                  IF (e .GT. dise(i)) THEN
                     DO
                        i = i + 1
                        IF (e .LE. dise(i)) EXIT
                     END DO
                  END IF
                  d%xe(e,1) = i
               END DO
            END IF
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

      ie = 0
      IF (nsd .EQ. 2) THEN
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
                  tmpX(:,ia) = 0D0
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
                     tmpX(:,ia) = 0D0
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
         IF (ALLOCATED(dmnId)) CALL MPI_GATHERV(tDmnId, sCe(cm%tF()), 
     2      mpint, d%xe(:,2), sCe, dise, mpint, master, cm%com(), ierr)

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
      END IF

      RETURN
      END SUBROUTINE INTNRBDATA

!**********************************************************************
!     This routine calculates average of VTU files.
      SUBROUTINE CALCAVE
      
      USE COMMOD
      USE ALLFUN
      
      IMPLICIT NONE
      
      LOGICAL :: lExist
      INTEGER :: fid, n, lTS, fTS, iTS
      CHARACTER(LEN=stdL) fName
      REAL(KIND=8) :: cntr
      
      REAL(KIND=8), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:),
     2   Ag(:,:), Yg(:,:), Dg(:,:)
      
      fid  = 1
      fTS  = nTS
      lTS  = 0
      cntr = 0D0
      iTS  = 0
      IF (zeroAve) iTS = rsTS
      IF (.NOT.rmsh%isReqd) THEN
         IF (cm%mas()) THEN
            ALLOCATE(tmpA(tDof,gtnNo), tmpY(tDof,gtnNo), 
     2         tmpD(tDof,gtnNo))
         ELSE
            ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
         END IF
         ALLOCATE(Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo))
         
         Ag   = 0D0
         Yg   = 0D0
         Dg   = 0D0
         std = " Computing average quantities"
         DO n=iTS, nTS, saveIncr
            IF (n .LT. saveATS) CYCLE
            IF (n .LT. fTS) fTS = n
            IF (n .GT. lTS) lTS = n
            cntr = cntr + 1D0
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
         n = NINT(cntr)
         IF (n .EQ. 0) RETURN
         Ag = Ag/cntr
         Yg = Yg/cntr
         Dg = Dg/cntr
         cTS = lTS
         fName = saveName
         saveName = TRIM(saveName)//"_ave_"//STR(fTS)//"_"//
     2      STR(saveIncr)
         CALL WRITEVTUS(Ag, Yg, Dg)
         saveName = fName
      ELSE
         DO n=iTS, nTS, saveIncr
            IF (n .LT. saveATS) CYCLE
            IF (n .LT. fTS) fTS = n
            IF (n .GT. lTS) lTS = n
            cntr = cntr + 1D0
         END DO
         n = NINT(cntr)
         IF (n.EQ.0 .OR. .NOT.ALLOCATED(rmsh%Aav)) RETURN
         rmsh%Aav = rmsh%Aav/cntr
         rmsh%Yav = rmsh%Yav/cntr
         rmsh%Dav = rmsh%Dav/cntr
         cTS = lTS
         fName = saveName
         saveName = TRIM(saveName)//"_ave_"//STR(fTS)//"_"//
     2      STR(saveIncr)
         CALL WRITEVTUS(rmsh%Aav, rmsh%Yav, rmsh%Dav)
         saveName = fName
      END IF
      
      RETURN
      END SUBROUTINE CALCAVE
      
!**********************************************************************

      SUBROUTINE READVTUS(lA, lY, lD, fName)
      
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      
      IMPLICIT NONE
      
      REAL(KIND=8), INTENT(INOUT) :: lA(tDof,gtnNo), lY(tDof,gtnNo),
     2   lD(tDof,gtnNo)
      CHARACTER(LEN=STDL) :: fName
      
      TYPE(vtkXMLType) :: vtu
      
      INTEGER :: iStat, iEq, iOut, iM, l, s, e, a, b, Ac, nNo, oGrp
      CHARACTER(LEN=stdL) :: varName
      REAL(KIND=8), ALLOCATABLE :: tmpS(:,:), tmpGS(:,:)
      
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
            varName = eq(iEq)%sym//"_"//TRIM(eq(iEq)%output(iOut)%name)
            
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

!**********************************************************************

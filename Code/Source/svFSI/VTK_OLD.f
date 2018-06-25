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
!     This routine generates VTK/VTKB files.
!      
!--------------------------------------------------------------------

      SUBROUTINE VTK(lA, lY, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: lA(tDof,tnNo), lY(tDof,tnNo), 
     2   lD(tDof,tnNo)
      
      LOGICAL isBinary
      INTEGER s, e, l, ie, is, iEq, fid, iOut, a, nNo, nEl, iM, nSh, 
     2   cnt, outDof, nOut, cOut, oGrp, Ac, Ec
      CHARACTER(LEN=stdL) fileName
      TYPE(dataType) d(nMsh)
      INTEGER, ALLOCATABLE :: outS(:), tmpI(:)
      REAL(KIND=8), ALLOCATABLE :: tmpV(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)
 
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
!--------------------------------------------------------------------
!     Preparing all solutions, e.g. velocity and pressure, and writing
!     them into d
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
!     Don't write it when it doesn't suppose to be written
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
               CASE (outGrp_vort,outGrp_eFlx,outGrp_hFlx,outGrp_stInv)
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
!     Now it is time to transform the local data to global form before
!     outputting it
      nNo = 0
      nEl = 0
      cnt = 0
      DO iM=1, nMsh
         IF (msh(iM)%eType .EQ. eType_NRB) THEN
            CALL INTNRBDATA(msh(iM), d(iM), outDof)
         ELSE
            CALL INTMSHDATA(msh(iM), d(iM), outDof)
         END IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
         cnt = cnt + (d(iM)%eNoN + 1)*d(iM)%nEl
      END DO
      IF (cm%slv()) THEN
         savedOnce = .TRUE.
         RETURN
      END IF

      DEALLOCATE(tmpV)
      ALLOCATE(tmpV(maxnsd,nNo))
      
      IF (saveFormat .EQ. saveF_VTK) THEN
         isBinary = .FALSE.
      ELSE
         isBinary = .TRUE.
      END IF

      fid  = 1
      IF (cTS .GE. 1000) THEN
         fileName = STR(cTS)
      ELSE
         WRITE(fileName,'(I3.3)') cTS
      END IF
      
      fileName = TRIM(saveName)//"_"//TRIM(ADJUSTL(fileName))//'.vtk'
      dbg = "Writing VTK"
      OPEN (fid, FILE=filename)
      IF (isBinary) THEN
         WRITE (fid,*) 1
         CLOSE (fid)
         OPEN (fid, FILE=fileName, STATUS='UNKNOWN', ACCESS='STREAM',
     2      FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
         WRITE (fid) "# vtk DataFile Version 3.0"//eol
         WRITE (fid) "Simulation results"//eol
         WRITE (fid) "BINARY"//eol
         WRITE (fid) "DATASET UNSTRUCTURED_GRID"//eol
         WRITE (fid) "POINTS "//STR(nNo)//" float"//eol
      ELSE
         WRITE (fid,'(A)') "# vtk DataFile Version 3.0"
         WRITE (fid,'(A)') "Simulation results"
         WRITE (fid,'(A)') "ASCII"
         WRITE (fid,'(A)') "DATASET UNSTRUCTURED_GRID"
         WRITE (fid,'(A)') "POINTS "//STR(nNo)//" float"
      END IF
!--------------------------------------------------------------------
!     Writing the Position data
      iOut = 1
      s    = outS(iOut)
      e    = outS(iOut+1) - 1
      nSh  = 0
      DO iM=1, nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%gx(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO   
      CALL WRITEVTK(nsd,nNo,tmpV,'',fid)
!--------------------------------------------------------------------
!     Writing the connectivity data
      IF (isBinary) THEN
         WRITE(fid) "CELLS "//STR(nEl)//" "//STR(cnt)//eol
         nSh = -1
         DO iM=1, nMsh
            DO e=1, d(iM)%nEl
               WRITE (fid) d(iM)%eNoN, d(iM)%IEN(:,e) + nSh
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         WRITE(fid) eol//"CELL_TYPES "//STR(nEl)//eol
         DO iM=1, nMsh
            DO e=1, d(iM)%nEl
               WRITE (fid) d(iM)%vtkType
            END DO
         END DO
         WRITE(fid) eol//"POINT_DATA "//STR(nNo)//eol
      ELSE
         WRITE(fid,'(A)')"CELLS "//STR(nEl)//" "//STR(cnt)
         nSh = -1
         DO iM=1, nMsh
            DO e=1, d(iM)%nEl
               WRITE (fid,'(A)',ADVANCE='NO') STR(d(iM)%eNoN)
               DO a=1, d(iM)%eNoN
                  WRITE (fid,'(A)',ADVANCE='NO') " "//
     2               STR(d(iM)%IEN(a,e) + nSh)
               END DO
               WRITE (fid,'(A)')
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         WRITE(fid,'(A)') "CELL_TYPES "//STR(nEl)
         DO iM=1, nMsh
            DO e=1, d(iM)%nEl
               WRITE (fid,'(A)') STR(d(iM)%vtkType)
            END DO
         END DO
         WRITE(fid,'(A)') "POINT_DATA "//STR(nNo)
      END IF
!     Writing all solutions, e.g. velocity and pressure
      DO iOut=2, nOut
         s    = outS(iOut)
         e    = outS(iOut+1) - 1
         l    = e - s + 1
         nSh  = 0
         DO iM=1, nMsh
            DO a=1, d(iM)%nNo
               tmpV(1:l,a+nSh) = d(iM)%gx(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL WRITEVTK(l,nNo,tmpV,outNames(iOut),fid)
      END DO

!     Here goes the element based variables
      IF (.NOT.savedOnce .OR. mvMsh) THEN
         IF (isBinary) THEN
            WRITE(fid) "CELL_DATA "//STR(nEl)//eol
         ELSE
            WRITE(fid,'(A)') "CELL_DATA "//STR(nEl)
         END IF
      END IF

      IF (.NOT.savedOnce) THEN
         savedOnce = .TRUE.
         ALLOCATE(tmpI(nEl))

!     Writing the Partitioning data
         IF (.NOT.cm%seq()) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(Ec) = d(iM)%xe(e,1)
               END DO
            END DO
            CALL WRITEVTKELE(nEl,tmpI,'Proc_ID',fid)
         END IF
!     Writing the domain data
         IF (ALLOCATED(dmnId)) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(Ec) = d(iM)%xe(e,2)
               END DO
            END DO
            CALL WRITEVTKELE(nEl,tmpI,'Domain_ID',fid)
         END IF
!     Writing the mesh data
         IF (nMsh .GT. 1) THEN
            Ec = 0
            DO iM=1, nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(Ec) = iM
               END DO
            END DO
            CALL WRITEVTKELE(nEl,tmpI,'Mesh_ID',fid)
         END IF
      END IF

      CALL FLUSH(fid)
      CLOSE (fid)

      RETURN
      END SUBROUTINE VTK

!--------------------------------------------------------------------
      SUBROUTINE WRITEVTKELE(nEl, tbw, paramName, fid)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nEl, fid, tbw(nEl)
      CHARACTER(*), INTENT(IN) :: paramName
 
      INTEGER e

      IF (saveFormat .EQ. saveF_VTKB) THEN
         WRITE(fid) "SCALARS "//TRIM(paramName)//" int"//eol
         WRITE(fid) "LOOKUP_TABLE default"//eol
         DO e=1, nEl
            WRITE(fid) tbw(e)
         END DO
         WRITE(fid) eol
      ELSE
         WRITE(fid,'(A)') "SCALARS "//TRIM(paramName)//" int"
         WRITE(fid,'(A)') "LOOKUP_TABLE default"
         DO e=1, nEl
            WRITE(fid,'(A)') STR(tbw(e))
         END DO
      END IF

      RETURN
      END SUBROUTINE WRITEVTKELE

!--------------------------------------------------------------------
      SUBROUTINE WRITEVTK(m, nNo, tbw, paramName, fid)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: m, nNo, fid
      REAL(KIND=8), INTENT(INOUT) :: tbw(maxnsd,nNo)
      CHARACTER(*), INTENT(IN) :: paramName

      LOGICAL flag, isBinary
      INTEGER i, Ac
 
      dbg = "Writing "//TRIM(paramName)//" with dof: "//m
      IF (m.LT.1 .OR. m.GT.3) err = "Unexpected dof in WRITEVTK "//
     2 TRIM(paramName)//" dof: "//m
      
      IF (saveFormat .EQ. saveF_VTK) THEN
         isBinary = .FALSE.
      ELSE
         isBinary = .TRUE.
      END IF
      
      flag = .TRUE.
      DO Ac=1, nNo
         DO i=1, m
            IF ((tbw(i,Ac).NE.tbw(i,Ac)) .OR. 
     2          (tbw(i,Ac).GT.HUGE(tbw))) THEN
               IF (flag) THEN
                  wrn = "Exceptional number detected, replaced by zero"
                  flag = .FALSE.
               END IF
               tbw(i,Ac) = 0D0
            END IF
         END DO
      END DO

      IF (paramName .NE. '') THEN
         IF (isBinary) THEN
            IF (m .EQ. 1) THEN
               WRITE(fid)"SCALARS "//TRIM(paramName)//" float"//eol
               WRITE(fid) "LOOKUP_TABLE default"//eol
            ELSE
               WRITE(fid)"VECTORS "//TRIM(paramName)//" float"//eol
            END IF
         ELSE
            IF (m .EQ. 1) THEN
               WRITE(fid,'(A)') "SCALARS "//TRIM(paramName)//" float"
               WRITE(fid,'(A)') "LOOKUP_TABLE default"
            ELSE
               WRITE(fid,'(A)') "VECTORS "//TRIM(paramName)//" float"
            END IF
         END IF
      END IF
      
      IF (isBinary) THEN
         IF (m .EQ. 2) THEN
            DO Ac=1, nNo
               WRITE (fid) REAL(tbw(1:m,Ac)), 0.0
            END DO
         ELSE
            DO Ac=1, nNo
               WRITE (fid) REAL(tbw(1:m,Ac))
            END DO
         END IF
      ELSE
         IF (m .EQ. 1) THEN
            DO Ac=1, nNo
               WRITE (fid,'(A)') STR(tbw(1,Ac))
            END DO
         ELSE IF (m .EQ. 2) THEN
            DO Ac=1, nNo
               WRITE (fid,'(A)') STR(tbw(1,Ac))//" "//
     2            STR(tbw(2,Ac))//" 0.0"
            END DO
         ELSE
            DO Ac=1, nNo
               WRITE (fid,'(A)') STR(tbw(1,Ac))//" "//
     2            STR(tbw(2,Ac))//" "//STR(tbw(3,Ac))
            END DO
         END IF
      END IF
      IF (isBinary) WRITE(fid) eol

      RETURN
      END SUBROUTINE WRITEVTK

!####################################################################
!     This routine prepares data array of a regular mesh
      SUBROUTINE INTMSHDATA(lM, d, outDof)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      TYPE(dataType), INTENT(INOUT) :: d
      INTEGER, INTENT(IN) :: outDof

      INTEGER e, i, ierr
      INTEGER, ALLOCATABLE :: sCount(:), disp(:)
      REAL(KIND=8), ALLOCATABLE :: ienU(:,:), gienU(:,:)

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
      IF (cm%mas()) d%IEN = lM%gIEN

      IF (.NOT.savedOnce) THEN
         ALLOCATE(d%xe(d%nEl,2))
         IF (cm%mas()) THEN
!     Based on the two posible fields
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
            END IF
         END IF
         IF (ALLOCATED(dmnId)) THEN
            ALLOCATE(sCount(cm%np()))
            DO i=1, cm%np()
               sCount(i) = lM%eDist(i) - lM%eDist(i-1)
            END DO
            CALL MPI_GATHERV(lM%eId, lM%nEl, mpint, d%xe(:,2), sCount, 
     2         lM%eDist, mpint, master, cm%com(), ierr)
         ELSE
            d%xe(:,2) = 1
         END IF
      END IF

      RETURN
      END SUBROUTINE INTMSHDATA

!--------------------------------------------------------------------
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

!####################################################################
!     This routine calculates average of VTK files.
      SUBROUTINE CALCAVE
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      LOGICAL ierr
      INTEGER fid, n, fN, lN
      REAL(KIND=8) cntr
      CHARACTER(LEN=stdL) fName
      
      REAL(KIND=8), ALLOCATABLE :: tmpA(:,:), tmpY(:,:), tmpD(:,:),
     2   A(:,:), Y(:,:), D(:,:)

      IF (cm%mas()) THEN
         ALLOCATE(tmpA(tDof,gtnNo), tmpY(tDof,gtnNo), tmpD(tDof,gtnNo))
      ELSE
         ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
      END IF
      ALLOCATE(A(tDof,tnNo), Y(tDof,tnNo), D(tDof,tnNo))

      fid  = 1
      fN   = nTS
      lN   = 0
      cntr = 0D0
      A    = 0D0
      Y    = 0D0
      D    = 0D0
      DO n=0, nTS, saveIncr
         IF (n .LT. saveATS) CYCLE
         IF (n .LT. fN) fN = n
         IF (n .GT. lN) lN = n
         cntr = cntr + 1D0
         IF (cm%mas()) THEN
            IF (n .GE. 1000) THEN
               fName = STR(n)
            ELSE
               WRITE(fName,'(I3.3)') n
            END IF
            fName = TRIM(saveName)//"_"//TRIM(ADJUSTL(fName))//'.vtk'
            INQUIRE (FILE=fName, EXIST=ierr)
            IF (ierr) THEN
               std = " Reading "//TRIM(fName)//
     2            " to calculate average file"
               CALL READVTK(fName, tmpA, tmpY, tmpD)
            ELSE
               err = "File "//TRIM(fName)//" is missing"
            END IF
         END IF
         A = A + LOCAL(tmpA)
         Y = Y + LOCAL(tmpY)
         D = D + LOCAL(tmpD)
      END DO
!     For any reason if there is no file to average, we return
      n = NINT(cntr)
      IF (n .EQ. 0) RETURN
      
      A = A/cntr
      Y = Y/cntr
      D = D/cntr
      
      cTS = lN
      fName = saveName
      saveName = TRIM(saveName)//"_ave_"//STR(fN)//"_"//STR(saveIncr)
!     Now saving a single vtk file with average_ prefix
      CALL VTK(A, Y, D)
      saveName = fName

      RETURN
      END SUBROUTINE CALCAVE

!--------------------------------------------------------------------
!     To read data from a vtk file. Only master calls this routine and
!     data are stored in global arrays through the arguments
      SUBROUTINE READVTK(fName, tmpA, tmpY, tmpD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      CHARACTER(LEN=stdL), INTENT(IN) :: fName
      REAL(KIND=8), INTENT(OUT) :: tmpA(tDof,gtnNo), tmpY(tDof,gtnNo), 
     2   tmpD(tDof,gtnNo)

      LOGICAL isBinary
      INTEGER i, a, e, fid, iEq, iOut, l, s, m, iPos, nNo, iM
      CHARACTER c
      CHARACTER(LEN=stdL) rLine, tmp
      REAL, ALLOCATABLE :: tmpGS(:,:), tmpS(:)

      tmpA = 0D0
      tmpY = 0D0
      tmpD = 0D0

      fid = 1
      OPEN (fid, FILE=fName, STATUS='OLD')
      
      isBinary = .FALSE.
      DO i=1, 20
         READ (fid,'(a)') rLine
         rLine = ADJUSTL(rLine)
         IF (rLine(1:6) .EQ. 'BINARY') THEN
            isBinary = .TRUE.
            CLOSE(fid)
            OPEN (fid, FILE=fName, STATUS='UNKNOWN', ACCESS='STREAM',
     2         FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
            iPos = 1
            EXIT
         END IF
      END DO
      REWIND(fid)
!     Skipping stuff until I get to a point which contains data
      DO
         IF (isBinary) THEN
!     Since reading from STREAM, we need to find the end of a line
!     manually
            rLine = ''
            DO i=1, stdL
               READ (fid,END=001,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=001) rLine
         END IF
         rLine = ADJUSTL(rLine)
         IF (rLine(1:10) .EQ. 'POINT_DATA') EXIT
      END DO
      rLine = rLine(11:)
      READ (rLine,*) a
      nNo = SUM(msh%gnNo)
      IF (a .NE. nNo) err = "Incompatible mesh and "//fName
      DO
         IF (isBinary) THEN
            rLine = ''
            DO i=1, stdL
               READ (fid,END=003,POS=iPos) c
               iPos = iPos + 1
               IF (c .EQ. eol) EXIT
               rLine(i:i) = c
            END DO
         ELSE
            READ (fid,'(a)',END=003) rLine
         END IF
         rLine = ADJUSTL(rLine)
!     Continuing until to get to a block of nodal data
         IF (rLine(1:7) .EQ. 'SCALARS') THEN
            m = 1
!     To skip lookup table statment            
            IF (isBinary) THEN
               DO i=1, stdL
                  READ (fid,END=003,POS=iPos) c
                  iPos = iPos + 1
                  IF (c .EQ. eol) EXIT
               END DO
            ELSE
               READ (fid,'(a)',END=001) tmp
            END IF
         ELSE IF (rLine(1:7) .EQ. 'VECTORS') THEN
            m = maxnsd
         ELSE 
            CYCLE
         END IF
!     Skipping VECTOR/SCALAR
         rLine = rLine(9:)
         i = LEN(TRIM(rLine))
!     Skipping 'float'
         rLine = rLine(1:i-6)
!     Finding the equation based on prefix
         DO iEq=1, nEq
            IF (eq(iEq)%sym .EQ. rLine(1:2)) EXIT
         END DO
         IF (iEq .GT. nEq) THEN
            wrn ="Skipping equation data <"//rLine(1:2)//"> in "//fName
            CYCLE
         END IF
         rLine = rLine(4:)
!     Finding output based on the name of variable
         DO iOut=1, eq(iEq)%nOutput
            IF (eq(iEq)%output(iOut)%name .EQ. rLine) EXIT
         END DO
         IF (iOut .GT. eq(iEq)%nOutput) THEN
            wrn = "Skipping output <"//TRIM(rLine)//"> in "//fName
            CYCLE
         ELSE
            std = " Reading <"//TRIM(rLine)//"> from "//fName
         END IF
         IF (ALLOCATED(tmpGS)) DEALLOCATE(tmpGS, tmpS)
         ALLOCATE(tmpGS(m,gtnNo), tmpS(m))
         s = eq(iEq)%s + eq(iEq)%output(iOut)%o
         IF (nMsh .EQ. 1) THEN
!     Reading data as single precision
            IF (isBinary) THEN
               DO a=1, gtnNo
                  READ(fid,END=001) tmpGS(:,a)
                  iPos = iPos + m*KIND(tmpGS)
               END DO
            ELSE
               DO a=1, gtnNo
                  READ(fid,*,END=001) tmpGS(:,a)
               END DO
            END IF
         ELSE
            IF (.NOT.cm%seq()) err = "Unable to read vtk with multiple"
     2         //" procs and msh"
            IF (isBinary) THEN
               DO iM=1, nMsh
                  DO a=1, msh(iM)%gnNo
                     READ(fid,END=001) tmpS
                     tmpGS(:,msh(iM)%gN(a)) = tmpS
                     iPos = iPos + m*KIND(tmpGS)
                  END DO
               END DO
            ELSE
               DO iM=1, nMsh
                  DO a=1, msh(iM)%gnNo
                     READ(fid,*,END=001) tmpS
                     tmpGS(:,msh(iM)%gN(a)) = tmpS
                     iPos = iPos + m*KIND(tmpGS)
                  END DO
               END DO
            END IF
         END IF

         l = eq(iEq)%output(iOut)%l
         s = eq(iEq)%s + eq(iEq)%output(iOut)%o
         e = s + l - 1
         SELECT CASE (eq(iEq)%output(iOut)%grp)
         CASE (outGrp_A)
            tmpA(s:e,:) = tmpGS(1:l,:)
         CASE (outGrp_Y)
            tmpY(s:e,:) = tmpGS(1:l,:)
         CASE (outGrp_D)
            tmpD(s:e,:) = tmpGS(1:l,:)
         CASE DEFAULT
            CYCLE
         END SELECT
      END DO        

 003  RETURN     
 001  err = "A block of data is missing from "//fName

      END SUBROUTINE READVTK


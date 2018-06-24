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

!**************************************************

      module vtkXMLlib
      use dataTypeParams
      use stdParams, only : strL

      logical, parameter :: debug = .false.
       character(len=64), parameter :: b64List = &
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

      integer(IK), parameter :: nVTKElms=5
      integer(IK), parameter :: nPieceTyps=5
      integer(IK), parameter :: nPieceElms=6
      integer(IK), parameter :: nPieceData=5
      integer(IK), parameter :: nPieceAtts=9
      integer(IK), parameter :: nDataElms=7
      integer(IK), parameter :: nDataTyps=10
      integer(IK), parameter :: nDataFrmt=3
      integer(IK), parameter :: nDataEnc=2

      character(len=strL), dimension(nVTKElms)   :: libVTKElms
      character(len=strL), dimension(nPieceTyps) :: libVTKPcTyps
      character(len=strL), dimension(nPieceElms) :: libPieceElms
      character(len=strL), dimension(nPieceData) :: libPcPtClData
      character(len=strL), dimension(nPieceAtts) :: libPieceAtts
      character(len=strL), dimension(nDataElms)  :: libDataElms
      character(len=strL), dimension(nDataTyps)  :: libDataTyps
      character(len=strL), dimension(nDataFrmt)  :: libDataFrmt
      character(len=strL), dimension(nDataFrmt)  :: libDataEnc

      contains

         subroutine initVTKXMLlib
         implicit none

         libVTKElms=""
         libVTKPcTyps=""
         libPieceElms=""
         libPcPtClData=""
         libPieceAtts=""
         libDataElms=""
         libDataTyps=""
         libDataFrmt=""

         ! libVTKElms(5) !
         libVTKElms(1) = "type"
         libVTKElms(2) = "version"
         libVTKElms(3) = "byte_order"
         libVTKElms(4) = "header_type"
         libVTKElms(5) = "compressor"

         ! libVTKPcTyps(5) !
         libVTKPcTyps(1) = "ImageData"
         libVTKPcTyps(2) = "RectilinearGrid"
         libVTKPcTyps(3) = "StructuredGrid"
         libVTKPcTyps(4) = "PolyData"
         libVTKPcTyps(5) = "UnstructuredGrid"

         ! libPieceElms(6) !
         libPieceElms(1) = "NumberOfPoints"
         libPieceElms(2) = "NumberOfCells"
         libPieceElms(3) = "NumberOfVerts"
         libPieceElms(4) = "NumberOfLines"
         libPieceElms(5) = "NumberOfStrips"
         libPieceElms(6) = "NumberOfPolys"

         ! libPieceData(5) !
         libPcPtClData(1)  = "Scalars"
         libPcPtClData(2)  = "Vectors"
         libPcPtClData(3)  = "Normals"
         libPcPtClData(4)  = "Tensors"
         libPcPtClData(5)  = "Tcoords"

         ! libPieceAtts(9) !
         libPieceAtts(1) = "PointData"
         libPieceAtts(2) = "CellData"
         libPieceAtts(3) = "Points"
         libPieceAtts(4) = "Coords"
         libPieceAtts(5) = "Verts"
         libPieceAtts(6) = "Lines"
         libPieceAtts(7) = "Strips"
         libPieceAtts(8) = "Polys"
         libPieceAtts(9) = "Cells"

         ! libDataElms(7) !
         libDataElms(1) = "type"
         libDataElms(2) = "Name"
         libDataElms(3) = "NumberOfComponents"
         libDataElms(4) = "format"
         libDataElms(5) = "offset"
         libDataElms(6) = "RangeMin"
         libDataElms(7) = "RangeMax"

         ! libDataTyps(10) !
         libDataTyps(1)  = "Int8"
         libDataTyps(2)  = "UInt8"
         libDataTyps(3)  = "Int16"
         libDataTyps(4)  = "UInt16"
         libDataTyps(5)  = "Int32"
         libDataTyps(6)  = "UInt32"
         libDataTyps(7)  = "Int64"
         libDataTyps(8)  = "UInt64"
         libDataTyps(9)  = "Float32"
         libDataTyps(10) = "Float64"

         ! libDataFrmt(3) !
         libDataFrmt(1) = "ascii"
         libDataFrmt(2) = "binary"
         libDataFrmt(3) = "appended"

         ! libDataEnc(2) !
         libDataEnc(1) = "raw"
         libDataEnc(2) = "base64"

         end subroutine initVTKXMLlib

      end module vtkXMLlib

!**************************************************

      module vtkXMLMod
      use stdParams
      use dataTypeParams
      use genUtils
      use vtkXMLlib
      implicit none

      private
      public :: vtkXMLType
      public :: loadVTK, flushVTK
      public :: vtkInitWriter
      public :: vtkWriteToFile
      public :: getVTK_nodesPerElem
      public :: getVTK_numElems
      public :: getVTK_numPoints
      public :: getVTK_elemIEN
      public :: getVTK_pointCoords
      public :: getVTK_pointData
      public :: getVTK_elemData
      public :: putVTK_pointCoords
      public :: putVTK_elemIEN
      public :: putVTK_pointData
      public :: putVTK_elemData

      interface getVTK_pointData
         module procedure getVTK_pointDataIntS, getVTK_pointDataRealS, &
                          getVTK_pointDataIntV, getVTK_pointDataRealV
      end interface getVTK_pointData

      interface getVTK_elemData
         module procedure getVTK_elemDataIntS, getVTK_elemDataRealS, &
                          getVTK_elemDataIntV, getVTK_elemDataRealV
      end interface getVTK_elemData

      interface putVTK_pointData
         module procedure putVTK_pointDataIntS, putVTK_pointDataRealS, &
                          putVTK_pointDataIntV, putVTK_pointDataRealV
      end interface putVTK_pointData

      interface putVTK_elemData
         module procedure putVTK_elemDataIntS, putVTK_elemDataRealS, &
                          putVTK_elemDataIntV, putVTK_elemDataRealV
      end interface putVTK_elemData

      logical :: flag
      character(len=strL) :: rLine,stmp
      character(len=strL) :: startKwrd,stopKwrd
      character(len=strL), dimension(maxToks) :: tokenList
      character :: c
      integer(IK) :: itok,ntoks,slen,iPos
      integer(IK) :: rank,maxRank
      integer(IK) :: pcStPos,pcEndPos

      type dataArrType
         private
         character(len=strL), dimension(nDataElms)  :: dElms
         character(len=strL) :: dType,dName,dFrmt,hdrType
         logical :: isInt
         integer(IK) :: hdrKind
         integer(IK) :: iOffst,appRank,ikind,rank
         integer(IK) :: stPos,endPos,nBytes
         integer(IK) :: nElms,nComps,nVals
         integer(IK4), dimension(:), allocatable :: iarr
         real(RK), dimension(:), allocatable :: darr
      end type dataArrType

      type pieceAttType
         private
         integer(IK) :: n,stPos,endPos
         character(len=strL) :: pName
         character(len=strL), dimension(nPieceData) :: ptClField
         character(len=strL), dimension(nPieceData) :: ptClFieldName
         type(dataArrType), dimension(:), allocatable :: dataArr
      end type pieceAttType

      type vtkXMLType
         private
         logical :: isBinApp
         character(len=strL) :: fileName
         character(len=strL) :: dataFormat
         character(len=strL) :: dataEncdng
         character(len=strL) :: vtkPcType
         integer(IK) :: fid
         integer(IK) :: stAppendPos
         integer(IK) :: endAppendPos
         integer(IK) :: offsets(100)

         character(len=strL), dimension(nVTKElms)   :: vtkElms
         integer(IK), dimension(nPieceElms) :: pieceElms
         type(pieceAttType), dimension(nPieceAtts) :: pcAtt
      end type vtkXMLType

      contains

         !==========================================

         subroutine loadVTK(vtk,fName,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: fName
         integer(IK) :: istat

         istat = 0
         inquire(file=trim(fName), exist=flag)
         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: File "//trim(fName)//" does not exist"
            istat=-1; return
         end if

         slen = len(trim(fName))
         select case (fName(slen-2:slen))
         case ("vtu","vtp")
         case default
            write(stdout,ftab4) &
               "ERROR: unknown file extension &
               &(can only be vtu or vtp)"
            istat=-1; return
         end select

!         write(stdout,ftab1) "<VTK XML Parser> Loading file <-----"// &
!            "  "//trim(fName)
         call initVTKXMLlib

         call initVTKXMLstruct(vtk,fName)

         call readHeader(vtk,istat)
         if ( istat.lt.0 ) return

         call parseVTKKernel(vtk,istat)
         if ( istat.lt.0 ) return

         call vtkDataLoader(vtk,istat)
         if ( istat.lt.0 ) return

!         write(stdout,ftab2) "Success!"
!         write(stdout,ftab1)

         return

         end subroutine loadVTK

         !==========================================

         subroutine flushVTK(vtk)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK) :: iatt,i

         do iatt=1, nPieceAtts
            if ( vtk%pcAtt(iatt)%n.lt.1 .or. &
                .not.allocated(vtk%pcAtt(iatt)%dataArr) ) cycle
            do i=1, vtk%pcAtt(iatt)%n
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
         end do
         call flush(vtk%fid)
         close(vtk%fid)

         end subroutine flushVTK

         !==========================================

         subroutine vtkFlushData(dArr)
         implicit none
         type(dataArrType), intent(inout) :: dArr

         if ( allocated(dArr%iarr) ) deallocate(dArr%iarr)
         if ( allocated(dArr%darr) ) deallocate(dArr%darr)
         dArr%dElms(:) = ""
         dArr%dType = ""
         dArr%dName = ""
         dArr%dFrmt = ""
         dArr%hdrType = ""

         end subroutine vtkFlushData

         !==========================================

         subroutine vtkDeepCopyData(src, dest)
         implicit none
         type(dataArrType), intent(inout) :: src, dest

         dest%dElms(:) =  src%dElms(:)
         dest%dType    =  src%dType
         dest%dName    =  src%dName
         dest%dFrmt    =  src%dFrmt
         dest%hdrType  =  src%hdrType
         dest%isInt    =  src%isInt
         dest%hdrKind  =  src%hdrKind
         dest%ikind    =  src%ikind
         dest%nComps   =  src%nComps
         dest%nVals    =  src%nVals
         dest%nElms    =  src%nElms
         if ( allocated(src%iarr) ) then
            allocate(dest%iarr(dest%nElms))
            dest%iarr(:)  =  src%iarr(:)
         end if
         if ( allocated(src%darr) ) then
            allocate(dest%darr(dest%nElms))
            dest%darr(:)  =  src%darr(:)
         end if

         end subroutine vtkDeepCopyData

         !==========================================

         subroutine initVTKXMLstruct(vtk,fName)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: fName
         integer :: fid,iatt

         vtk%fileName        = trim(fName)
         vtk%isBinApp        = .false.
         vtk%dataFormat      = ""
         vtk%dataEncdng      = ""
         vtk%vtkPcType       = ""
         vtk%stAppendPos     = 0
         vtk%endAppendPos    = 0
         vtk%offsets(:)      = 0
         vtk%vtkElms(:)      = ""
         vtk%pieceElms(:)    = 0
         do iatt=1, nPieceAtts
            vtk%pcAtt(iatt)%pName  = ""
            vtk%pcAtt(iatt)%ptClField(:) = ""
            vtk%pcAtt(iatt)%ptClFieldName(:) = ""
         end do

         do fid=11, 1024
            inquire(unit=fid, opened=flag)
            if ( .not.flag ) exit
         end do

         vtk%fid = fid

         end subroutine initVTKXMLstruct

         !==========================================

         subroutine readHeader(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat

         rLine = ""
         stmp = ""
         tokenList(:) = ""
         iPos = 0

         open(vtk%fid,file=trim(vtk%fileName),status='unknown')
         read(vtk%fid,'(A)') rLine
         rLine = adjustl(rLine)
         if ( rLine(2:5).eq."?xml" ) then
            read(vtk%fid,'(A)') rLine
         end if

         ! get data format and encoding !
         call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat)
         if ( istat.lt.0 ) return
         iPos = 0
         call parseString(rLine,tokenList,ntoks)
         stmp = getTokenValue(tokenList,ntoks,"format")
         if ( trim(stmp).eq."binary" .or. &
             trim(stmp).eq."appended" ) then
            close(vtk%fid)
            open(vtk%fid,file=trim(vtk%fileName),form="unformatted", &
            access="stream",convert="big_endian")
            vtk%isBinApp = .true.
            iPos = 1
         end if
         vtk%dataFormat = trim(stmp)
         if ( debug ) write(stdout,ftab2) &
            "Data format: <"//trim(vtk%dataFormat)//">"

         vtk%dataEncdng = "base64"
         if ( vtk%dataFormat.eq."appended" ) then
            call findKwrdXML(vtk,"<AppendedData",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            do
               read(vtk%fid,pos=iPos,end=001) c
               iPos = iPos + 1
               if ( c.eq."_" ) exit
            end do
            vtk%stAppendPos = iPos
            call parseString(rLine,tokenList,ntoks)
            stmp = getTokenValue(tokenList,ntoks,"encoding")
            vtk%dataEncdng = trim(stmp)
            call findKwrdXML(vtk,"</AppendedData",iPos,rLine,istat)
            if ( istat.lt.0 ) return
            vtk%endAppendPos = iPos - len(trim(rLine))
         end if

         if ( debug ) then
            write(stdout,ftab2) &
            "Data encoding: <"//trim(vtk%dataEncdng)//">"
            if ( vtk%dataFormat.eq."appended" ) then
               write(stdout,ftab2) &
                  "Data begins at pos: "//trim(STR(vtk%stAppendPos))
               write(stdout,ftab2) &
                  "Data ends at pos: "// trim(STR(vtk%endAppendPos))
            end if
         end if

         iPos = 0
         if ( vtk%isBinApp ) iPos = 1
         rewind(vtk%fid)
         return

 001      write(stdout,ftab4) "ERROR: end of file reached.."
         istat = -1; return

         end subroutine readHeader

         !==========================================

         subroutine parseVTKKernel(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat

         rLine = ""
         stmp = ""
         tokenList(:) = ""

         call findKwrdXML(vtk,"<VTKFile",iPos,rLine,istat)
         if ( istat.lt.0 ) return
         if ( debug ) write(stdout,ftab2) trim(rLine)
         call parseString(rLine,tokenList,ntoks)

         do itok=1, nVTKElms
            stmp = getTokenValue(tokenList,ntoks,libVTKElms(itok))
            vtk%vtkElms(itok) = trim(stmp)
         end do

         do itok=1, nPieceTyps
            if ( trim(vtk%vtkElms(1)).eq.trim(libVTKPcTyps(itok)) ) exit
         end do
         if ( itok.le.3 ) then
            write(stdout,ftab4) &
               "ERROR: unknown piece type"
            write(stdout,ftab4) &
               "Piece <"//trim(tokenList(1))//'>'
            write(stdout,ftab4) &
               "Piece can be <UnstructuredGrid> or <PolyData> only"
            istat=-1; return
         end if
         vtk%vtkPcType = libVTKPcTyps(itok)
         if ( debug .and. len(trim(vtk%vtkElms(5))).gt.0 ) &
            write(stdout,ftab2) &
               "Data compression: "//trim(vtk%vtkElms(5))

         ! Piece elements parser !
         call findKwrdXML(vtk,"<Piece",iPos,rLine,istat)
         if ( istat.lt.0 ) return
         if ( debug ) write(stdout,ftab2) trim(rLine)
         call parseString(rLine,tokenList,ntoks)

         do itok=1, nPieceElms
            stmp = getTokenValue(tokenList,ntoks,libPieceElms(itok))
            slen = len(trim(stmp))
            if ( slen.gt.0 ) then
               read(stmp(1:slen),*) vtk%pieceElms(itok)
            end if
         end do

         pcStPos = iPos! - len(trim(rLine))
         call findKwrdXML(vtk,"</Piece",iPos,rLine,istat)
         if ( istat.lt.0 ) return
         if ( vtk%isBinApp ) then
            pcEndPos = iPos - len(trim(rLine))
         else
            pcEndPos = iPos - 1
         end if

         call readPieceAttributes(vtk,istat)

         return
         end subroutine parseVTKKernel

         !==========================================

         subroutine readPieceAttributes(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer :: cntr,iatt,i

         rLine = ""
         stmp = ""
         tokenList(:) = ""
         cntr = 0

         do iatt=1, nPieceAtts
            call resetFilePos(vtk%fid,vtk%fileName,pcStPos,vtk%isBinApp)
            iPos = pcStPos
            call findKwrdXML(vtk,"<"//libPieceAtts(iatt),iPos,rLine,istat,pcEndPos)
            if ( istat.lt.0 ) return
            if ( len(trim(rLine)).lt.1 ) cycle
            if ( debug ) write(stdout,ftab2) trim(rLine)

            call parseString(rLine,tokenList,ntoks)
            stmp = getToken(tokenList,ntoks,libPieceAtts(iatt))
            vtk%pcAtt(iatt)%pName = trim(stmp)

            if ( iatt.le.2 ) then
               do itok=1, nPieceData
                  stmp = getToken(tokenList,ntoks,libPcPtClData(itok))
                  vtk%pcAtt(iatt)%ptClField(itok) = trim(stmp)
                  stmp = getTokenValue(tokenList,ntoks,libPcPtClData(itok))
                  vtk%pcAtt(iatt)%ptClFieldName(itok) = trim(stmp)
               end do
            else
               vtk%pcAtt(iatt)%ptClField(:) = ""
               vtk%pcAtt(iatt)%ptClFieldName(:) = ""
            end if

            vtk%pcAtt(iatt)%stPos = iPos! - len(trim(rLine))
            call findKwrdXML(vtk,"</"//libPieceAtts(iatt),iPos,rLine,istat,pcEndPos)
            if ( istat.lt.0 ) return
            if ( vtk%isBinApp ) then
               vtk%pcAtt(iatt)%endPos= iPos - len(trim(rLine))
            else
               vtk%pcAtt(iatt)%endPos= iPos - 1
            end if

            vtk%pcAtt(iatt)%n = 0
            iPos = vtk%pcAtt(iatt)%stPos
            if ( .not.vtk%isBinApp ) &
               call resetFilePos(vtk%fid,vtk%fileName,iPos,vtk%isBinApp)

            ! Count DataArray elements within a Piece attribute !
            if ( debug ) then
               write(stdout,ftab3) &
                  "Counting <DataArray> elements in <"//&
                  trim(vtk%pcAtt(iatt)%pName)//">"
            end if
            do
               call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat,vtk%pcAtt(iatt)%endPos)
               if ( istat.lt.0 ) return
               if ( iPos.ge.vtk%pcAtt(iatt)%endPos ) exit
               if ( debug ) write(stdout,ftab3) trim(rLine)
               call parseString(rLine,tokenList,ntoks)
               if ( trim(tokenList(1)).eq."DataArray" ) &
                  vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            end do ! inner loop over DataArray elms !
            if ( vtk%pcAtt(iatt)%n.eq.0 ) then
               if ( debug ) write(stdout,ftab3) "None found.."
               cycle
            end if

            if ( debug ) then
               write(stdout,ftab3) "No. of <DataArray> elems: "//&
                  trim(STR(vtk%pcAtt(iatt)%n))
            end if

            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            iPos = vtk%pcAtt(iatt)%stPos
            call resetFilePos(vtk%fid,vtk%fileName,iPos,vtk%isBinApp)
            if ( debug ) then
               write(stdout,ftab2) &
                  "Reading <DataArray> elements in <"//&
                  trim(vtk%pcAtt(iatt)%pName)//">"
            end if
            do i=1, vtk%pcAtt(iatt)%n
               call findKwrdXML(vtk,"<DataArray",iPos,rLine,istat)
               if ( istat.lt.0 ) return
               vtk%pcAtt(iatt)%dataArr(i)%stPos = iPos+1

               if ( debug ) write(stdout,ftab3) trim(rLine)

               call parseString(rLine,tokenList,ntoks)
               do itok=1, nDataElms
                  stmp = getTokenValue(tokenList,ntoks,libDataElms(itok))
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(itok) = trim(stmp)
               end do

               if ( vtk%dataFormat.ne."appended" ) then
                  call findKwrdXML(vtk,"</DataArray>",iPos,rLine,istat)
                  if ( istat.lt.0 ) return
                  if ( vtk%isBinApp ) then
                     vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                        iPos - len(trim(rLine))
                  else
                     vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                        iPos
                  end if
               end if

               cntr = cntr+1
               vtk%pcAtt(iatt)%dataArr(i)%rank = cntr

               vtk%pcAtt(iatt)%dataArr(i)%dType = &
                  trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(1))
               vtk%pcAtt(iatt)%dataArr(i)%dName = &
                  trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(2))

               slen = len(trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(3)))
               if ( slen.gt.0 ) then
                  read(vtk%pcAtt(iatt)%dataArr(i)%dElms(3)(1:slen),*) &
                     vtk%pcAtt(iatt)%dataArr(i)%nComps
               else
                  vtk%pcAtt(iatt)%dataArr(i)%nComps = 0 !1
               end if

               vtk%pcAtt(iatt)%dataArr(i)%dFrmt = &
                  trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(4))

               slen = len(trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(5)))
               if ( slen.gt.0 ) then
                  read(vtk%pcAtt(iatt)%dataArr(i)%dElms(5)(1:slen),*) &
                     vtk%pcAtt(iatt)%dataArr(i)%iOffst
               else
                  vtk%pcAtt(iatt)%dataArr(i)%iOffst = 0
               end if

               if ( vtk%isBinApp ) then
                  if ( vtk%dataFormat.eq."binary" ) then
                     call adjustdataArray(vtk%fid, &
                     vtk%pcAtt(iatt)%dataArr(i)%stPos, &
                     vtk%pcAtt(iatt)%dataArr(i)%endPos )

                     iPos = vtk%pcAtt(iatt)%dataArr(i)%stPos
                     vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                        vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                        vtk%pcAtt(iatt)%dataArr(i)%stPos
                  else if ( vtk%dataFormat.eq."appended" ) then
                     vtk%pcAtt(iatt)%dataArr(i)%stPos = &
                        vtk%pcAtt(iatt)%dataArr(i)%iOffst + &
                        vtk%stAppendPos
                     vtk%offsets(vtk%pcAtt(iatt)%dataArr(i)%rank) = &
                        vtk%pcAtt(iatt)%dataArr(i)%stPos
                  end if
               else
                  vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                     vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                     vtk%pcAtt(iatt)%dataArr(i)%stPos
               end if

            end do ! idata
         end do ! outer loop over Piece atts !
         maxRank = cntr

         return
         end subroutine readPieceAttributes

         !==========================================

         subroutine vtkDataLoader(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,rank

         if ( debug ) write(stdout,ftab2) &
            "Piece Type: "//trim(vtk%vtkPcType)
         do iatt=1, nPieceAtts
            if ( vtk%pcAtt(iatt)%n.lt.1 .or. &
            .not.allocated(vtk%pcAtt(iatt)%dataArr) ) cycle
            if ( debug ) write(stdout,ftab3) &
               "Piece Attribute: "//trim(vtk%pcAtt(iatt)%pName)
            do i=1, vtk%pcAtt(iatt)%n
               if ( vtk%dataFormat.eq."appended" ) then
                  rank = vtk%pcAtt(iatt)%dataArr(i)%rank
                  do
                     rank = rank+1
                     if ( rank.gt.maxRank ) then
                        ipos = vtk%endAppendPos-1
                        do
                           read(vtk%fid,pos=iPos,end=001) c
                           if ( c.eq.' ' .or. c.eq.eol .or. &
                              c.eq.'   ') then
                              iPos = iPos-1
                              cycle
                           else
                              exit
                           end if
                        end do
                        vtk%pcAtt(iatt)%dataArr(i)%endPos = iPos
                        exit
                     end if

                     if ( vtk%offsets(rank).gt. &
                         vtk%pcAtt(iatt)%dataArr(i)%stPos ) then
                         vtk%pcAtt(iatt)%dataArr(i)%endPos = &
                           vtk%offsets(rank)
                        exit
                     end if
                  end do
                  vtk%pcAtt(iatt)%dataArr(i)%nbytes = &
                     vtk%pcAtt(iatt)%dataArr(i)%endPos - &
                        vtk%pcAtt(iatt)%dataArr(i)%stPos
               end if ! appended

               if ( vtk%dataFormat.eq."ascii" ) then
                  call resetFilePos(vtk%fid,vtk%fileName,&
                     vtk%pcAtt(iatt)%dataArr(i)%stPos-1,vtk%isBinApp)
               end if

               if ( debug ) write(stdout,ftab4) &
                  "DataArray name: "// &
                  trim(vtk%pcAtt(iatt)%dataArr(i)%dName)

               call vtkXMLDataParser(vtk,vtk%pcAtt(iatt)%dataArr(i),iatt,i,istat)
               if ( istat.lt.0 ) return
            end do
            if (debug) write(stdout,ftab1) repeat("*",76)
         end do

         return

 001     write(stdout,ftab4) "ERROR: end of file reached.."
         istat = -1; return

         end subroutine vtkDataLoader

         !==========================================

         subroutine vtkXMLDataParser(vtk,dataArr,iatt,idata,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt,idata
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr

         dPtr => dataArr

         if ( len(trim(vtk%vtkElms(4))).gt.0 ) then
            dPtr%hdrType = trim(vtk%vtkElms(4))
         else
            dPtr%hdrType = "UInt32"
         end if
         call selectDataType(dPtr%hdrType,dPtr%hdrKind,dPtr%isInt,istat)
         if ( istat.lt.0 ) return

         if ( debug ) then
            write(stdout,ftab1) repeat("*",76)
            write(stdout,ftab3) &
               "Piece Type: "//trim(vtk%vtkPcType)
            write(stdout,ftab3) &
               "Piece Attribute: "//trim(vtk%pcAtt(iatt)%pName)
            write(stdout,ftab3) &
               "Header type: "//trim(dPtr%hdrType)
            write(stdout,ftab3) &
               "DataArray type: "//trim(dPtr%dType)
            write(stdout,ftab3) &
               "DataArray name: "//trim(dPtr%dName)
            write(stdout,ftab3) &
               "DataArray frmt: "//trim(dPtr%dFrmt)
            write(stdout,ftab3) &
               "DataArray nComps: "//trim(STR(dPtr%nComps))
            write(stdout,ftab3) &
               "DataArray offset: "//trim(STR(dPtr%ioffst))
            write(stdout,ftab3) &
               "DataArray start pos: "//trim(STR(dPtr%stPos))
            write(stdout,ftab3) &
               "DataArray end pos:   "//trim(STR(dPtr%endPos))
            write(stdout,ftab3) &
               "DataArray nbytes:   "//trim(STR(dPtr%nbytes))
         end if

         select case (trim(vtk%vtkPcType))
         case ("UnstructuredGrid")
            call vtkParseUnstrucGrid(vtk,dataArr,iatt,idata,istat)
            if ( istat.lt.0 ) return

         case ("PolyData")
            call vtkParsePolyData(vtk,dataArr,iatt,idata,istat)
            if ( istat.lt.0 ) return

         end select

         if ( vtk%isBinApp ) then
            if ( vtk%vtkElms(5).eq."vtkZLibDataCompressor" ) then
               call readZlibBinaryData(vtk,dPtr,istat)
            else
               if ( dPtr%nElms.eq.0 ) &
                  return
               call readBinaryData(vtk,dPtr,istat)
            end if
         else
            if ( dPtr%nElms.eq.0 ) &
               return
            call readAsciiData(vtk,dPtr,istat)
         end if

         if ( istat.lt.0 ) return

         end subroutine vtkXMLDataParser

         !==========================================

         subroutine vtkParseUnstrucGrid(vtk,dataArr,iatt,idata,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt,idata
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nCells

         dPtr => dataArr

         nPoints = vtk%pieceElms(1)
         nCells  = vtk%pieceElms(2)
         if ( nPoints.lt.1 ) then
            write(stdout,ftab4) &
            "ERROR: VTK Piece element NumberOfPoints not defined.."
            istat=-1; return
         end if

         if ( nCells.lt.1 ) then
            write(stdout,ftab4) &
            "ERROR: VTK Piece element NumberOfCells not defined.."
            istat=-1; return
         end if

         dPtr%isInt = .false.
         call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)
         if ( istat.lt.0 ) return

         select case (trim(vtk%pcAtt(iatt)%pName))
         case ("PointData")
            if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nPoints
            dPtr%nElms = nPoints * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPoints "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("CellData")
            if ( nCells.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nCells)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nCells
            dPtr%nElms = nCells * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nCells "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Points")
            if ( dPtr%nComps.ne.3 ) then
               write(stdout,ftab4) &
               "WARNING: Element <NumberOfComponents> in <Points> &
                &attribute < 3"
            end if
            if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) &
               dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
            dPtr%nVals = nPoints
            dPtr%nElms = nPoints * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPoints "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Cells")
            if ( nCells.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nCells)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nCells
            dPtr%nElms = nCells * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nCells "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         end select

         end subroutine vtkParseUnstrucGrid

         !==========================================

         subroutine vtkParsePolyData(vtk,dataArr,iatt,idata,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt,idata
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nVerts,nLines,nStrips,nPolys

         dPtr => dataArr

         nPoints = vtk%pieceElms(1)
         nVerts  = vtk%pieceElms(3)
         nLines  = vtk%pieceElms(4)
         nStrips = vtk%pieceElms(5)
         nPolys  = vtk%pieceElms(6)

         if ( nPoints.lt.1 ) then
            write(stdout,ftab4) &
            "ERROR: VTK Piece element NumberOfPoints not defined.."
            istat=-1; return
         end if

         if ( nPolys.lt.1 ) then
            write(stdout,ftab4) &
            "ERROR: VTK Piece element NumberOfPolys not defined.."
            istat=-1; return
         end if

         dPtr%isInt = .false.
         call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)

         select case (trim(vtk%pcAtt(iatt)%pName))
         case ("PointData")
            if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nPoints
            dPtr%nElms = nPoints * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPoints "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("CellData")
            if ( nPolys.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nPolys)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nPolys
            dPtr%nElms = nPolys * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPolys "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Points")
            if ( nPoints.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nPoints)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            if ( dPtr%nComps.ne.3 ) then
               write(stdout,ftab4) &
               "WARNING: Element <NumberOfComponents> in <Points> &
                &attribute < 3"
            end if
            dPtr%nVals = nPoints
            dPtr%nElms = nPoints * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPoints "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Verts")
            if ( nVerts.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nVerts)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nVerts
            dPtr%nElms = nVerts * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nVerts "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Lines")
            if ( nLines.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nLines)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nLines
            dPtr%nElms = nLines * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nLines "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Strips")
            if ( nStrips.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nStrips)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nStrips
            dPtr%nElms = nStrips * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nStrips "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if

         case ("Polys")
            if ( nPolys.gt.0 .and. dPtr%nComps.eq.0 ) then
               dPtr%nComps = getNumComps(vtk,dataArr,nPolys)
               if ( .not.vtk%isBinApp ) dPtr%nComps = 1
            end if
            dPtr%nVals = nPolys
            dPtr%nElms = nPolys * dPtr%nComps
            if ( debug ) then
               write(stdout,ftab4) &
                  "nPolys "// trim(STR(dPtr%nVals)) //&
                  "; nComps "// trim(STR(dPtr%nComps)) //&
                  "; nElems "// trim(STR(dPtr%nElms))
            end if
         end select

         end subroutine vtkParsePolyData

         !==========================================

         function getNumComps(vtk,dArr,m) result(n)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         type(dataArrType), intent(inout) :: dArr
         integer(IK), intent(in) :: m
         integer(IK) :: n,npadd

         tokenList(:) = ""
         rLine = ""
         if ( vtk%isBinApp ) then
            n = dArr%nbytes
            npadd = 0
            if ( vtk%dataEncdng.eq."base64" ) then
               n = n*3_IK/4_IK
               if ( mod(n,3_IK).gt.0 ) &
                  npadd = 3_IK-mod(n,3_IK)
            end if
            n = n - npadd - dArr%hdrKind
            n = n / dArr%ikind / m
         else
            read(vtk%fid,'(A)') rLine
            call parseString(trim(rLine),tokenList,ntoks)
            if ( ntoks.ne.0 ) n = ntoks
            call resetFilePos(vtk%fid,vtk%fileName,dArr%stPos-1,vtk%isBinApp)
         end if

         end function getNumComps

         !==========================================

         subroutine readBinaryData(vtk,dPtr,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr

         character(len=:), allocatable :: code
         integer(IK) :: i,j,hdr,ikind
         integer(IK) :: np,np1,np2,npadd
         integer(IK1), allocatable, dimension(:) :: pIK1,p1,p2

         ikind = dPtr%ikind

         np  = dPtr%nBytes
         np1 = 1_IK*dPtr%hdrKind
         np2 = dPtr%nelms * ikind
         npadd = 0_IK

         allocate(pIK1(np),p1(np1),p2(np2))
         pIK1 = 0_IK1; p1 = 0_IK1; p2 = 0_IK1

         if ( vtk%dataEncdng.eq."base64" ) then
            code = repeat(" ",np)
            j = 0_IK
            do i=dPtr%stPos, dPtr%endPos-1
               j = j+1_IK
               read(vtk%fid,pos=i,end=001) code(j:j)
            end do

            np = np*3_IK/4_IK
            if ( mod(np1+np2,3_IK).ne.0 ) &
               npadd = 3_IK - mod(np1+np2,3_IK)
         end if

         if ( np.ne.(np1+np2+npadd) ) then
            write(stdout,ftab4) &
               "ERROR: inconsistent data array dimension.."
            write(stdout,ftab4) &
               "Expected size: "//trim(STR(np1+np2+npadd))
            write(stdout,ftab4) &
               "Actual size in file: "//trim(STR(np))
            istat=-1; return
         end if

         if ( vtk%dataEncdng.eq."base64" ) then
            call decode_bits(code,pIK1)
         else
            read(vtk%fid,pos=dPtr%stPos,end=001) pIK1(1:np)
         end if

         p1(1:np1) = pIK1(1:np1)
         p2(1:np2) = pIK1(np1+1:np1+np2)

         select case(dPtr%hdrkind)
         case(IK1)
            call transferBits(p1,np1,int(ikind,kind=IK1),hdr)
         case(IK2)
            call transferBits(p1,np1,int(ikind,kind=IK2),hdr)
         case(IK4)
            call transferBits(p1,np1,int(ikind,kind=IK4),hdr)
         case(IK8)
            call transferBits(p1,np1,int(ikind,kind=IK8),hdr)
         end select

         if ( dPtr%isInt ) then
            allocate(dPtr%iarr(dPtr%nElms))
            select case(ikind)
            case(IK1)
               call transferBits(p2,np2,int(ikind,kind=IK1),dPtr%iarr,dPtr%nelms)
            case(IK2)
               call transferBits(p2,np2,int(ikind,kind=IK2),dPtr%iarr,dPtr%nelms)
            case(IK4)
               call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%iarr,dPtr%nelms)
            case(IK8)
               call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%iarr,dPtr%nelms)
               if ( debug ) write(stdout,ftab4) &
                  "WARNING: Typecasting from INT64 to INT32 &
                  &could lead to errors.."
            end select
         else
            allocate(dPtr%darr(dPtr%nelms))
            select case(ikind)
            case(RK4)
               call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%darr,dPtr%nelms)
            case(RK8)
               call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%darr,dPtr%nelms)
            end select
         end if

         return

 001     write(stdout,ftab4) "ERROR: end of file reached.."
         istat=-1; return

         end subroutine readBinaryData

         !==========================================

         subroutine readZlibBinaryData(vtk,dPtr,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr

         character(len=:), allocatable :: code
         integer(IK) :: hdr(3),nlen,ikind
         integer(IK) :: i,j,ist,iend,iBlk
         integer(IK) :: np,np1,np2,npadd
         integer(IK), dimension(:), allocatable :: szBlk
         integer(IK1), allocatable, dimension(:) :: pIK1,p1,p2

         ikind = dPtr%ikind

         iPos = dPtr%stPos
         np1 = 3_IK * dPtr%hdrKind
         allocate(p1(np1)); p1=0_IK1
         npadd = 0_IK

         if ( vtk%dataEncdng.eq."base64" ) then
            nlen = 4_IK * dPTr%hdrKind
            code = repeat(" ",nlen)
            j = 0_IK
            do i=iPos, iPos+nlen-1_IK
               j = j+1_IK
               read(vtk%fid,pos=i,end=001) code(j:j)
            end do
            call decode_bits(code,p1)
            iPos = iPos + nlen
         else
            read(vtk%fid,pos=iPos,end=001) p1(1:np1)
            iPos = iPos + np1
         end if

         select case ( dPtr%hdrKind )
         case (IK1)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK1),hdr,3_IK)
         case (IK2)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK2),hdr,3_IK)
         case (IK4)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK4),hdr,3_IK)
         case (IK8)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK8),hdr,3_IK)
         end select

         if ( debug ) then
            write(stdout,ftab4) &
               "num blocks "//trim(STR(hdr(1)))
            write(stdout,ftab4) &
               "size of block (before compr.) "// &
               trim(STR(hdr(2)))
            write(stdout,ftab4) &
               "size of last block (before compr.) "// &
               trim(STR(hdr(3)))
         end if

         if ( hdr(1).lt.1_IK ) return
         if ( hdr(3) .eq. 0_IK ) hdr(3) = hdr(2)

         np1 = (hdr(2)*(hdr(1)-1) + hdr(3))/ikind

         if ( np1.ne.dPtr%nElms ) then
            dPtr%nElms = np1
            dPtr%nComps = dPtr%nElms / dPtr%nVals
            if ( debug ) then
               write(stdout,ftab4) "Reset : nComps "// &
                  trim(STR(dPtr%nComps))//" "//"nElms "//&
                  trim(STR(dPtr%nElms))
            end if
         end if
         deallocate(p1)

         np1 = hdr(1) * dPtr%hdrKind
         allocate(p1(np1)); p1(:) = 0_IK1
         allocate(szBlk(hdr(1))); szBlk(:) = 0_IK

         npadd = 0_IK
         if ( vtk%dataEncdng.eq."base64" ) then
            if ( mod(np1,3_IK).gt.0_IK ) &
               npadd = 3_IK-mod(np1,3_IK)
            nlen = (np1+npadd)*4_IK/3_IK
            j = 0_IK
            code = repeat(" ",nlen)
            do i=iPos, iPos+nlen-1_IK
               j = j+1_IK
               read(vtk%fid,pos=i,end=001) code(j:j)
            end do
            call decode_bits(code,p1)
            iPos = iPos + nlen
         else
            read(vtk%fid,pos=iPos,end=001) p1(1:np1)
            iPos = iPos + np1
         end if

         select case ( dPtr%hdrKind )
         case(IK1)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK1),szBlk,hdr(1))
         case(IK2)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK2),szBlk,hdr(1))
         case(IK4)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK4),szBlk,hdr(1))
         case(IK8)
            call transferBits(p1,np1,int(dPtr%hdrKind,kind=IK8),szBlk,hdr(1))
         end select

         if ( debug ) then
            do iBlk=1_IK, hdr(1)
               write(stdout,ftab4) &
                  "size of block#"//trim(STR(iBlk))// &
                  " (after compr.) "//trim(STR(szBlk(iBlk)))
            end do
         end if

         if ( dPtr%isInt ) then
            allocate(dPtr%iarr(dPtr%nElms))
         else
            allocate(dPtr%darr(dPtr%nElms))
         end if

         if ( vtk%dataEncdng.eq."base64" ) then
            npadd = 0_IK
            np = sum(szBlk)
            if ( allocated(pIK1) ) deallocate(pIK1)
            allocate(pIK1(np)); pIK1 = 0_IK1
            if ( mod(np,3_IK).gt.0_IK ) &
               npadd = 3_IK-mod(np,3_IK)
            nlen = (np + npadd)*4_IK/3_IK
            code = repeat(" ",nlen)
            j = 0_IK
            do i=iPos, iPos+nlen-1_IK
               j = j+1_IK
               read(vtk%fid,pos=i,end=001) code(j:j)
            end do
            call decode_bits(code,pIK1)
            iPos = iPos + nlen
         end if

         ist = 1_IK; np = 0_IK
         do iBlk=1_IK, hdr(1)
            np1 = szBlk(iBlk)
            np2 = hdr(2)
            if ( iBlk.eq.hdr(1) ) np2 = hdr(3)

            if ( allocated(p1) ) deallocate(p1)
            if ( allocated(p2) ) deallocate(p2)
            allocate(p1(np1)); p1(:) = 0_IK1
            allocate(p2(np2)); p2(:) = 0_IK1

            if ( vtk%dataEncdng.eq."base64" ) then
               p1(1:np1) = pIK1(np+1:np+np1)
               np = np+np1
            else
               read(vtk%fid,pos=ipos) p1(1:np1)
               iPos = iPos + np1
            end if

            call infZlibData(p1,np1,p2,np2,istat)
            if ( istat.lt.0 ) return

            nlen = np2/ikind
            iend = ist-1 + nlen

            if ( dPtr%isInt ) then
               select case ( ikind )
               case (IK1)
                  call transferBits(p2,np2,int(ikind,kind=IK1),dPtr%iarr(ist:iend),nlen)
               case (IK2)
                  call transferBits(p2,np2,int(ikind,kind=IK2),dPtr%iarr(ist:iend),nlen)
               case (IK4)
                  call transferBits(p2,np2,int(ikind,kind=IK4),dPtr%iarr(ist:iend),nlen)
               case (IK8)
                  call transferBits(p2,np2,int(ikind,kind=IK8),dPtr%iarr(ist:iend),nlen)
               end select
            else
               select case ( ikind )
               case (RK4)
                  call transferBits(p2,np2,int(ikind,kind=RK4),dPtr%darr(ist:iend),nlen)
               case (RK8)
                  call transferBits(p2,np2,int(ikind,kind=RK8),dPtr%darr(ist:iend),nlen)
               end select
               ist = iend+1
            end if
            ist = iend+1
         end do

         return

 001      write(stdout,ftab4) "ERROR: end of file reached.."
         istat=-1; return

         end subroutine readZlibBinaryData

         !==========================================

         subroutine decode_bits(code,bits)
         implicit none
         character(len=*), intent(in) :: code
         integer(IK1), intent(out) :: bits(:)

         integer(IK1) :: sixb(1:4)
         integer(IK) :: c,e,Nb,i

         Nb = size(bits,dim=1,kind=IK)
         e = 1_IK
         do c=1_IK, len(code), 4_IK
            sixb(:) = 0_IK1
            sixb(1) = index(b64List,code(c  :c  ),kind=IK1) - 1_IK1
            sixb(2) = index(b64List,code(c+1:c+1),kind=IK1) - 1_IK1
            sixb(3) = index(b64List,code(c+2:c+2),kind=IK1) - 1_IK1
            sixb(4) = index(b64List,code(c+3:c+3),kind=IK1) - 1_IK1

            call mvbits(sixb(1),0,6,bits(e),2)
            call mvbits(sixb(2),4,2,bits(e),0)

            if ( e+1.le.Nb ) then
               call mvbits(sixb(2),0,4,bits(e+1),4)
               call mvbits(sixb(3),2,4,bits(e+1),0)
            end if

            if ( e+2.le.Nb ) then
               call mvbits(sixb(3),0,2,bits(e+2),6)
               call mvbits(sixb(4),0,6,bits(e+2),0)
            end if
            e = e+3_IK
         end do

         end subroutine decode_bits

         !==========================================

         pure subroutine encode_bits(bits,npadd,code)
         implicit none
         integer(IK1), intent(in) :: bits(:)
         integer(IK4), intent(in) :: npadd
         character(len=*), intent(out) :: code

         integer(IK1) :: sixb(1:4)
         integer(IK ) :: c,e,Nb

         Nb = size(bits,dim=1,kind=IK )
         c = 1_IK
         do e=1_IK, Nb, 3_IK
            sixb(:) = 0_IK1
            call mvbits(bits(e),2,6,sixb(1),0)
            call mvbits(bits(e),0,2,sixb(2),4)
            if ( e+1_IK .le. Nb ) then
               call mvbits(bits(e+1),4,4,sixb(2),0)
               call mvbits(bits(e+1),0,4,sixb(3),2)
            end if
            if ( e+2_IK .le. Nb ) then
               call mvbits(bits(e+2),6,2,sixb(3),0)
               call mvbits(bits(e+2),0,6,sixb(4),0)
            end if
            sixb(:) = sixb(:) + 1_IK1
            code(c  :c  ) = b64List(sixb(1):sixb(1))
            code(c+1:c+1) = b64List(sixb(2):sixb(2))
            code(c+2:c+2) = b64List(sixb(3):sixb(3))
            code(c+3:c+3) = b64List(sixb(4):sixb(4))
            c = c+4_IK
         end do

         if ( npadd.gt.0 ) &
            code(len(code)-npadd+1:) = repeat('=',npadd)

         end subroutine encode_bits

         !==========================================

         subroutine readAsciiData(vtk,dPtr,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: i

         if ( dPtr%isInt ) then
            allocate(dPtr%iarr(dPtr%nElms))
            read(vtk%fid,*,end=001) (dPtr%iarr(i),i=1, dPtr%nElms)
            if ( debug ) then
               write(stdout,'(13X,A)',advance="no")
               do i=1, dPtr%nElms
                  write(stdout,'(A)',advance="no") " "//&
                     trim(STR(dPtr%iarr(i)))
               end do
               write(stdout,'(A)')
            end if
         else
            allocate(dPtr%darr(dPtr%nElms))
            read(vtk%fid,*,end=001) (dPtr%darr(i),i=1, dPtr%nElms)
            if ( debug ) then
               write(stdout,'(13X,A)',advance="no")
               do i=1, dPtr%nElms
                  write(stdout,'(A)',advance="no") " "//&
                     STR(dPtr%darr(i))
               end do
               write(stdout,'(A)')
            end if
         end if

         return

 001      write(stdout,ftab4) "ERROR: end of file reached.."
         istat=-1; return

         end subroutine readAsciiData

         !==========================================

         subroutine findKwrdXML(vtk,sKwrd,iPos,strng,istat,ePos)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: iPos
         character(len=*), intent(in) :: sKwrd
         character(len=strL), intent(out) :: strng
         integer(IK), intent(in), optional :: ePos
         integer(IK), intent(inout) :: istat

         integer(IK) :: i,kwrdL,cnt
         character :: c

         kwrdL = len(trim(sKwrd))
         do
            strng = ''
            if ( vtk%isBinApp ) then
               do i=1, strL
                  read(vtk%fid,pos=iPos,end=001) c
                  iPos = iPos + 1
                  if ( c.eq. '<' ) then
                     cnt = 1
                     strng(cnt:cnt) = c
                     do
                        read(vtk%fid,pos=iPos,end=001) c
                        iPos = iPos + 1
                        if ( present(ePos) .and. iPos.ge.ePos) &
                        then
                           strng = ''
                           return
                        end if
                        cnt = cnt+1
                        strng(cnt:cnt) = c
                        if ( cnt.lt.kwrdL ) then
                           if (strng(1:cnt).ne.sKwrd(1:cnt)) &
                              exit
                        else if ( cnt.eq.kwrdL ) then
                           if (strng(1:cnt).eq.trim(skwrd)) then
                              cycle
                           else
                              strng = ''
                              exit
                           end if
                        end if
                        if ( c.eq.'>' ) exit
                     end do
                     strng = adjustl(strng)
                     if ( strng(1:kwrdL).eq.trim(sKwrd) ) return
                  end if
                  if(c.eq.eol) exit
               end do
            else
               strng = ''
               read(vtk%fid,'(A)',end=001) strng
               iPos = iPos+1
               if ( present(ePos) .and. iPos.ge.ePos ) then
                  strng = ''
                  return
               end if
            end if
            strng = adjustl(strng)
            if ( strng(1:kwrdL).eq.trim(sKwrd) ) return
         end do

         return

 001     write(stdout,ftab4) "ERROR: end of file reached.."
         write(stdout,ftab4) 'Failed processing for keyword "'// &
            trim(skwrd)//'"'
         istat = -1; return

         end subroutine findKwrdXML

         !==========================================

         subroutine adjustDataArray(fid,spos,epos)
         implicit none
         integer(IK), intent(in) :: fid
         integer(IK), intent(inout) :: spos,epos
         integer(IK) :: npos
         character :: c

         write(stdout,'(A)') trim(STR(spos))//" "//trim(STR(epos))
         npos = spos
         do
            read(fid,pos=npos) c
            if ( c.eq.' ' .or. c.eq.eol .or. c.eq.'   ') then
               npos = npos+1
               cycle
            else
               exit
            end if
         end do
         spos = npos

         do
            read(fid,pos=npos) c
            write(stdout,'(A)') trim(STR(npos))//" '"// &
            c//"' "//trim(STR(ichar(c)))
            if ( c.eq.' ' .or. c.eq.eol .or. c.eq.'   ' ) exit
            npos = npos+1
         end do
         epos = npos-1

         end subroutine adjustDataArray

         !==========================================

         subroutine resetFilePos(fileId,fName,fPos,isBin)
         implicit none
         integer, intent(in) :: fileId,fPos
         logical, intent(in) :: isBin
         character(len=*), intent(in) :: fName
         integer(IK), parameter :: SEEK_SET=0
         integer(IK) :: i,ierr

         if ( isBin ) then
            call fseek(fileId,fPos,SEEK_SET,ierr)
         else
            close(fileId)
            open(fileId,file=trim(fName),status='old')
            do i=1, fPos
               read(fileId,*)
            end do
         end if

         end subroutine resetFilePos

         !==========================================

         subroutine selectDataType(varType,ikind,isint,istat)
         implicit none
         character(len=*), intent(in) :: varType
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: ikind
         logical, intent(out) :: isint

         isint = .true.
         ikind = 0
         select case (trim(varType))
         case ("Int8", "UInt8")
            ikind = IK1
         case ("Int16", "UInt16")
            ikind = IK2
         case ("Int32", "UInt32")
            ikind = IK4
         case ("Int64", "UInt64")
            ikind = IK8
         case ("Float32")
            ikind = RK4
            isint = .false.
         case ("Float64","Double")
            ikind = RK8
            isint = .false.
         case default
            write(stdout,ftab4) "ERROR: unknown data type.."
            istat=-1; return
         end select

         end subroutine selectDataType

         !==========================================

         subroutine getVTK_numPoints(vtk,nn,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: nn
         integer(IK) :: i

         nn = vtk%pieceElms(1)
         if ( nn.eq.0 ) then
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute"
            istat=-1; return
         end if

         end subroutine getVTK_numPoints

         !==========================================

         subroutine getVTK_numElems(vtk,ne,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: ne

         if ( vtk%pieceElms(2).gt.0 ) then
            ne = vtk%pieceElms(2)
         else if ( vtk%pieceElms(6).gt.0 ) then
            ne = vtk%pieceElms(6)
         end if

         if ( ne.eq.0 ) then
            write(stdout,ftab4) &
               "ERROR: could not find CELLS or POLYS attributes"
            istat=-1; return
         end if

         end subroutine getVTK_numElems

         !==========================================

         subroutine getVTK_nodesPerElem(vtk,eNoN,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: eNoN
         integer(IK) :: iatt,i,itmp

         if ( vtk%pieceElms(2).gt.0 ) then   ! nCells !
            iatt = 9
         else if (vtk%pieceElms(6).gt.0 ) then ! nPolys !
            iatt = 8
         else
            write(stdout,ftab4) &
               "ERROR: could not find CELLS or POLYS attributes"
            istat=-1; return
         end if

         eNoN = 0; itmp = 0
         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               "connectivity" ) then
               eNoN = vtk%pcAtt(iatt)%dataArr(i)%nComps
            end if
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               "offsets" ) then
               itmp = vtk%pcAtt(iatt)%dataArr(i)%iarr(2) - &
                     vtk%pcAtt(iatt)%dataArr(i)%iarr(1)
            end if
         end do

         eNoN = max(eNoN, itmp)
         if ( eNoN.eq.0 ) then
            write(stdout,ftab4) &
               "ERROR: unexpected VTK behavior (eNoN)"
            istat=-1; return
         end if

         end subroutine getVTK_nodesPerElem

         !==========================================

         subroutine getVTK_elemIEN(vtk,ien,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(out) :: ien(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,ne,eNoN,i,j,k,l

         eNoN = size(ien,1)
         ne   = size(ien,2)
         if ( vtk%pieceElms(2).gt.0 ) then   ! nCells !
            iatt = 9
         else if (vtk%pieceElms(6).gt.0 ) then ! nPolys !
            iatt = 8
         else
            write(stdout,ftab4) &
               "ERROR: could not find CELLS or POLYS attributes"
            istat=-1; return
         end if

         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               "connectivity" ) then
               if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. &
                   (eNoN*ne) ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in IEN params.."
                  write(stdout,ftab4) &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                  write(stdout,ftab4) &
                     trim(STR(eNoN))//" "//trim(STR(ne))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                  write(stdout,ftab4) &
                     "ERROR: connectivity array unallocated.."
                  istat=-1; return
               end if

               l = 0
               do j=1, ne
                  do k=1, eNoN
                     l = l+1
                     ien(k,j) = &
                        vtk%pcAtt(iatt)%dataArr(i)%iarr(l)
                  end do
               end do
               exit
            end if
         end do

         if ( i.gt.vtk%pcAtt(iatt)%n ) then
            write(stdout,ftab4) &
               "ERROR: could not find connectivity in Cells or &
               &Polys attributes"
            istat=-1; return
         end if

         end subroutine getVTK_elemIEN

         !==========================================

         subroutine getVTK_pointCoords(vtk,x,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         real(RK), intent(out) :: x(:,:)
         integer(IK) :: iatt,i,j,k,l,nd,nn

         nd = size(x,1)
         nn = size(x,2)
         if ( vtk%pieceElms(1).eq.nn ) then   ! nPoints !
            iatt = 3
         else
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute"
            istat=-1; return
         end if

         i=1
         if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .ne. &
            "Points" ) then
            write(stdout,ftab4) &
               "WARNING: <DataArray> name does not match Points.."
         end if

         if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. (nd*nn) ) then
            write(stdout,ftab4) &
               "ERROR: mismatch in POINTS params.."
            write(stdout,ftab4) &
               trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
            write(stdout,ftab4) &
               trim(STR(nd))//" "//trim(STR(nn))
            istat=-1; return
         end if

         if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
            write(stdout,ftab4) &
               "ERROR: Points array unallocated.."
            istat=-1; return
         end if

         l = 0
         do j=1, nn
            do k=1, nd
               l = l+1
               x(k,j) = vtk%pcAtt(iatt)%dataArr(i)%darr(l)
            end do
         end do

         end subroutine getVTK_pointCoords

         !==========================================

         subroutine getVTK_pointDataIntS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         n = size(u)
         if ( vtk%pieceElms(1).gt.0 ) then   ! nPoints !
            iatt = 1
         else
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute to read &
               &PointData"
            istat=-1; return
         end if

         flag = .false.
         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               trim(kwrd) ) then
               flag = .true.
               if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array size and numNodes"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                  write(stdout,ftab4) "Input size: "// &
                     trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               u(:) = vtk%pcAtt(iatt)%dataArr(i)%iarr(:)
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in PointData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_pointDataIntS

         !==========================================

         subroutine getVTK_pointDataIntV(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,j,k,m,n
         logical :: flag

         m = size(u,1)
         if ( m.eq.1 ) then
            call getVTK_pointDataIntS(vtk,kwrd,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(1).gt.0 ) then   ! nPoints !
            iatt = 1
         else
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute to read &
               &PointData"
            istat=-1; return
         end if

         flag = .false.
         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               trim(kwrd) ) then
               flag = .true.
               if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                   (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array sizes"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                  write(stdout,ftab4) "Input size: "// &
                     trim(STR(m))//", "//trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               do j=1, n
                  do k=1, m
                     u(k,j) = &
                       vtk%pcAtt(iatt)%dataArr(i)%iarr((j-1)*m+k)
                  end do
               end do
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in PointData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_pointDataIntV

         !==========================================

         subroutine getVTK_pointDataRealS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         REAL(RK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         n = size(u)
         if ( vtk%pieceElms(1).gt.0 ) then   ! nPoints !
            iatt = 1
         else
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute to read &
               &PointData"
            istat=-1; return
         end if

         flag = .false.
         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               trim(kwrd) ) then
               flag = .true.
               if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array size and numNodes"
                  write(stdout,ftab4) &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                  write(stdout,ftab4) trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               u(:) = vtk%pcAtt(iatt)%dataArr(i)%darr(:)
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in PointData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_pointDataRealS

         !==========================================

         subroutine getVTK_pointDataRealV(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         real(RK), intent(inout) :: u(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,j,k,m,n
         logical :: flag

         m = size(u,1)
         if ( m.eq.1 ) then
            call getVTK_pointDataRealS(vtk,kwrd,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(1).gt.0 ) then   ! nPoints !
            iatt = 1
         else
            write(stdout,ftab4) &
               "ERROR: could not find POINTS attribute to read &
               &PointData"
            istat=-1; return
         end if

         flag = .false.
         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
               trim(kwrd) ) then
               flag = .true.
               if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                   (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array sizes"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                  write(stdout,ftab4) "Input size: "// &
                     trim(STR(m))//", "//trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               do j=1, n
                  do k=1, m
                     u(k,j) = &
                       vtk%pcAtt(iatt)%dataArr(i)%darr((j-1)*m+k)
                  end do
               end do
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in PointData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_pointDataRealV

         !==========================================

         subroutine getVTK_elemDataIntS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         n = size(u)
         if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
            iatt = 2
         else
            write(stdout,ftab4) &
               "ERROR: could not find either CELLS or POLYS &
               &attributes to read CellData"
            istat=-1; return
         end if

         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                trim(kwrd) ) then
               flag = .true.
               if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array size and numElems"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                  write(stdout,ftab4) "Input size: "//&
                     trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               u(:) = vtk%pcAtt(iatt)%dataArr(i)%iarr(:)
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in CellData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_elemDataIntS

          !==========================================

         subroutine getVTK_elemDataIntV(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,j,k,m,n
         logical :: flag

         m = size(u,1)
         if ( m.eq.1 ) then
            call getVTK_elemDataIntS(vtk,kwrd,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
            iatt = 2
         else
            write(stdout,ftab4) &
               "ERROR: could not find either CELLS or POLYS &
               &attributes to read CellData"
            istat=-1; return
         end if

         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                trim(kwrd) ) then
               flag = .true.
               if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                   (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array sizes"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                  write(stdout,ftab4) "Input size: "// &
                     trim(STR(m))//", "//trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%iarr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               do j=1, n
                  do k=1, m
                     u(k,j) = &
                       vtk%pcAtt(iatt)%dataArr(i)%iarr((j-1)*m+k)
                  end do
               end do
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in CellData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_elemDataIntV

         !==========================================

         subroutine getVTK_elemDataRealS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         real(RK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         n = size(u)
         if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
            iatt = 2
         else
            write(stdout,ftab4) &
               "ERROR: could not find either CELLS or POLYS &
               &attributes to read CellData"
            istat=-1; return
         end if

         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                trim(kwrd) ) then
               flag = .true.
               if ( vtk%pcAtt(iatt)%dataArr(i)%nElms .ne. n ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array size and numElems"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nElms))
                  write(stdout,ftab4) "Input size: "//&
                     trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               u(:) = vtk%pcAtt(iatt)%dataArr(i)%darr(:)
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in CellData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_elemDataRealS

          !==========================================

         subroutine getVTK_elemDataRealV(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         real(RK), intent(inout) :: u(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,j,k,m,n
         logical :: flag

         m = size(u,1)
         if ( m.eq.1 ) then
            call getVTK_elemDataRealS(vtk,kwrd,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(2).gt.0 .or. vtk%pieceElms(6).gt.0 ) then
            iatt = 2
         else
            write(stdout,ftab4) &
               "ERROR: could not find either CELLS or POLYS &
               &attributes to read CellData"
            istat=-1; return
         end if

         do i=1, vtk%pcAtt(iatt)%n
            if ( trim(vtk%pcAtt(iatt)%dataArr(i)%dName) .eq. &
                trim(kwrd) ) then
               flag = .true.
               if ( (vtk%pcAtt(iatt)%dataArr(i)%nVals  .ne. n) .or. &
                   (vtk%pcAtt(iatt)%dataArr(i)%nComps .ne. m) ) then
                  write(stdout,ftab4) &
                     "ERROR: mismatch in array sizes"
                  write(stdout,ftab4) "Actual size: "//&
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nComps))//", "// &
                     trim(STR(vtk%pcAtt(iatt)%dataArr(i)%nVals))
                  write(stdout,ftab4) "Input size: "// &
                     trim(STR(m))//", "//trim(STR(n))
                  istat=-1; return
               end if
               if ( .not.allocated(vtk%pcAtt(iatt)%dataArr(i)%darr) ) then
                  write(stdout,ftab4) &
                     "ERROR: data not found.."
                  istat=-1; return
               end if

               do j=1, n
                  do k=1, m
                     u(k,j) = &
                       vtk%pcAtt(iatt)%dataArr(i)%darr((j-1)*m+k)
                  end do
               end do
               exit
            end if
         end do

         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: could not find <"//trim(kwrd)//"> &
               &in CellData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_elemDataRealV

         !==========================================

         subroutine vtkInitWriter(vtk,fName,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: fName
         integer, intent(inout) :: istat
         integer :: fid

         if ( debug ) &
            write(stdout,ftab1) "<VTK XML Writer> Initializing..."

         istat = 0
         call initVTKXMLlib

         vtk%fileName = trim(fName)
         vtk%isBinApp = .true.
         vtk%dataFormat = "appended"
         vtk%dataEncdng = "raw"
         slen = len(trim(fName))
         select case (fName(slen-2:slen))
         case ("vtu")
            vtk%vtkPcType = "UnstructuredGrid"
            call vtkInitUGridWriter(vtk, istat)
         case ("vtp")
            vtk%vtkPcType = "PolyData"
            call vtkInitPDataWriter(vtk, istat)
         case default
            write(stdout,ftab4) &
               "ERROR: unknown file extension &
               &(can only be vtu or vtp)"
            istat=-1; return
         end select

         vtk%vtkElms(1) = vtk%vtkPcType
         vtk%vtkElms(2) = "0.1"
         vtk%vtkElms(3) = "LittleEndian"
         vtk%vtkElms(4) = ""
         vtk%vtkElms(5) = "vtkZLibDataCompressor"
         vtk%pieceElms(:) = 0

         do fid=11, 1024
            inquire(unit=fid, opened=flag)
            if ( .not.flag ) exit
         end do
         vtk%fid = fid

         open(10, status='scratch', form='unformatted', &
            access='stream', convert='big_endian')

         end subroutine vtkInitWriter

         !==========================================

         subroutine vtkInitUGridWriter(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer, intent(inout) :: istat
         integer :: iatt,i

         istat = 0
         do iatt=1, nPieceAtts
            vtk%pcAtt(iatt)%pName = ""
            vtk%pcAtt(iatt)%n = 0
            vtk%pcAtt(iatt)%ptClField(:) = ""
            vtk%pcAtt(iatt)%ptClFieldName(:) = ""
            select case (trim(libPieceAtts(iatt)))
            case ('Cells')
               vtk%pcAtt(iatt)%n = 3
               vtk%pcAtt(iatt)%pName = "Cells"

               allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
               vtk%pcAtt(iatt)%dataArr(1:2)%dType = "Int64"
               vtk%pcAtt(iatt)%dataArr(3)%dType   = "UInt8"

               vtk%pcAtt(iatt)%dataArr(1)%dName   = "connectivity"
               vtk%pcAtt(iatt)%dataArr(2)%dName   = "offsets"
               vtk%pcAtt(iatt)%dataArr(3)%dName   = "types"

               vtk%pcAtt(iatt)%dataArr(:)%dFrmt   = "appended"

               vtk%pcAtt(iatt)%dataArr(:)%nComps  = 1
               vtk%pcAtt(iatt)%dataArr(:)%nVals   = 0
               vtk%pcAtt(iatt)%dataArr(:)%nElms   = 0

               do i=1, vtk%pcAtt(iatt)%n
                  vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"
                  vtk%pcAtt(iatt)%dataArr(i)%isInt   = .true.
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dType
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
                  vtk%pcAtt(iatt)%dataArr(i)%dName
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dFrmt
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = ""
                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
                     vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if

                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
                     vtk%pcAtt(iatt)%dataArr(i)%iKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if
               end do

            end select
         end do

         end subroutine vtkInitUGridWriter

         !==========================================

         subroutine vtkInitPDataWriter(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer, intent(inout) :: istat
         integer :: iatt,i

         istat = 0
         do iatt=1, nPieceAtts
            vtk%pcAtt(iatt)%pName = ""
            vtk%pcAtt(iatt)%n = 0
            vtk%pcAtt(iatt)%ptClField(:) = ""
            vtk%pcAtt(iatt)%ptClFieldName(:) = ""
            select case (trim(libPieceAtts(iatt)))
            case ('Verts')
               vtk%pcAtt(iatt)%n = 2
               vtk%pcAtt(iatt)%pName = "Verts"

               allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
               vtk%pcAtt(iatt)%dataArr(1:2)%dType = "Int64"

               vtk%pcAtt(iatt)%dataArr(1)%dName   = "connectivity"
               vtk%pcAtt(iatt)%dataArr(2)%dName   = "offsets"

               vtk%pcAtt(iatt)%dataArr(:)%dFrmt   = "appended"
               vtk%pcAtt(iatt)%dataArr(:)%nComps  = 1
               vtk%pcAtt(iatt)%dataArr(:)%nVals   = 0
               vtk%pcAtt(iatt)%dataArr(:)%nElms   = 0

               do i=1, vtk%pcAtt(iatt)%n
                  vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"
                  vtk%pcAtt(iatt)%dataArr(i)%isInt   = .true.
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dType
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
                  vtk%pcAtt(iatt)%dataArr(i)%dName
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dFrmt
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = ""
                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
                     vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if

                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
                     vtk%pcAtt(iatt)%dataArr(i)%iKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if
               end do
            case ('Lines')
               vtk%pcAtt(iatt)%n = 2
               vtk%pcAtt(iatt)%pName = "Lines"

               allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
               vtk%pcAtt(iatt)%dataArr(1:2)%dType = "Int64"

               vtk%pcAtt(iatt)%dataArr(1)%dName   = "connectivity"
               vtk%pcAtt(iatt)%dataArr(2)%dName   = "offsets"

               vtk%pcAtt(iatt)%dataArr(:)%dFrmt   = "appended"
               vtk%pcAtt(iatt)%dataArr(:)%nComps  = 1
               vtk%pcAtt(iatt)%dataArr(:)%nVals   = 0
               vtk%pcAtt(iatt)%dataArr(:)%nElms   = 0

               do i=1, vtk%pcAtt(iatt)%n
                  vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"
                  vtk%pcAtt(iatt)%dataArr(i)%isInt   = .true.
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dType
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
                  vtk%pcAtt(iatt)%dataArr(i)%dName
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dFrmt
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = ""
                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
                     vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if

                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
                     vtk%pcAtt(iatt)%dataArr(i)%iKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if
               end do
            case ('Strips')
               vtk%pcAtt(iatt)%n = 2
               vtk%pcAtt(iatt)%pName = "Strips"

               allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
               vtk%pcAtt(iatt)%dataArr(1:2)%dType = "Int64"

               vtk%pcAtt(iatt)%dataArr(1)%dName   = "connectivity"
               vtk%pcAtt(iatt)%dataArr(2)%dName   = "offsets"

               vtk%pcAtt(iatt)%dataArr(:)%dFrmt   = "appended"
               vtk%pcAtt(iatt)%dataArr(:)%nComps  = 1
               vtk%pcAtt(iatt)%dataArr(:)%nVals   = 0
               vtk%pcAtt(iatt)%dataArr(:)%nElms   = 0

               do i=1, vtk%pcAtt(iatt)%n
                  vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"
                  vtk%pcAtt(iatt)%dataArr(i)%isInt   = .true.
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dType
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
                  vtk%pcAtt(iatt)%dataArr(i)%dName
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dFrmt
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = ""
                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
                     vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if

                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
                     vtk%pcAtt(iatt)%dataArr(i)%iKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if
               end do
            case ('Polys')
               vtk%pcAtt(iatt)%n = 2
               vtk%pcAtt(iatt)%pName = "Polys"

               allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
               vtk%pcAtt(iatt)%dataArr(1:2)%dType = "Int64"

               vtk%pcAtt(iatt)%dataArr(1)%dName   = "connectivity"
               vtk%pcAtt(iatt)%dataArr(2)%dName   = "offsets"

               vtk%pcAtt(iatt)%dataArr(:)%dFrmt   = "appended"
               vtk%pcAtt(iatt)%dataArr(:)%nComps  = 1
               vtk%pcAtt(iatt)%dataArr(:)%nVals   = 0
               vtk%pcAtt(iatt)%dataArr(:)%nElms   = 0

               do i=1, vtk%pcAtt(iatt)%n
                  vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"
                  vtk%pcAtt(iatt)%dataArr(i)%isInt   = .true.
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dType
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
                  vtk%pcAtt(iatt)%dataArr(i)%dName
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
                     vtk%pcAtt(iatt)%dataArr(i)%dFrmt
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = ""
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = ""
                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
                     vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if

                  call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
                     vtk%pcAtt(iatt)%dataArr(i)%iKind, &
                     vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
                  if (istat .lt. 0) then
                     inquire(unit=10, opened=flag)
                     if (flag) close(10)
                     return
                  end if
               end do
            end select
         end do

         end subroutine vtkInitPDataWriter

         !==========================================

         subroutine vtkWriteToFile(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer, intent(inout) :: istat
         integer :: fid,i,iatt,ivar,tmpI(100),ios
         logical, dimension(nPieceAtts) :: dAttToW
         character c

         if ( debug )write(stdout,ftab1) &
            "<VTK XML Writer> Writing to file ----->  "// &
            trim(vtk%fileName)

         istat = 0
         fid = vtk%fid
         inquire(unit=10, opened=flag)
         if ( fid.lt.11 .OR. fid.gt.1024 .OR. .not.flag) then
            write(stdout,ftab4) "ERROR: VTK structure not properly "// &
            " initialized. Call vtkInitWriter()"
            istat = -1; return
         end if

         ! write data first !
         write(10) '  <AppendedData encoding="raw">'//newl//'  _'
         ivar = 0
         vtk%offsets = 0
         tmpI = 0
         do iatt=1, nPieceAtts
            if (vtk%pcAtt(iatt)%n .gt. 0) then
               do i=1, vtk%pcAtt(iatt)%n
                  if (.NOT.allocated(vtk%pcAtt(iatt)%dataArr)) then
                     write(stdout,ftab4) "ERROR: Unexpected "// &
                     "behavior. VTK data structure unallocated."
                     istat=-1; close(10); return
                  end if
                  call vtkWriteZlibBinaryData(vtk,vtk%pcAtt(iatt)%dataArr(i),istat)
                  ivar = ivar + 1
                  tmpI(ivar) = vtk%pcAtt(iatt)%dataArr(i)%nBytes
                  if (ivar .gt. 1) vtk%offsets(ivar) = &
                     vtk%offsets(ivar-1) + tmpI(ivar-1)
                  vtk%pcAtt(iatt)%dataArr(i)%iOffst = vtk%offsets(ivar)
                  vtk%pcAtt(iatt)%dataArr(i)%dElms(5) = &
                     STR(vtk%pcAtt(iatt)%dataArr(i)%iOffst)
               end do
            end if
         end do
         write(10) newl//"  </AppendedData>"//newl
         write(10) "</VTKFile>"

         open(fid,file=trim(vtk%fileName),form="unformatted", &
            access="stream",convert="big_endian")
         write(fid) '<VTKFile '
         do i=1, nVTKElms
            if (len(trim(vtk%vtkElms(i))) .gt. 0) then
               write(fid) trim(libVTKElms(i))//'="'// &
                  trim(vtk%vtkElms(i))//'" '
            end if
         end do
         write(fid) '>'//newl

         write(fid) '  <'//trim(vtk%vtkPcType)//'>'//newl
         write(fid) '    <Piece '
         do i=1, nPieceElms
            if (vtk%pieceElms(i) .gt. 0) then
               write(fid) trim(libPieceElms(i))//'="'// &
                  trim(STR(vtk%pieceElms(i)))//'"    '
            end if
         end do
         write(fid) '>'//newl

         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            dAttToW = .false.
            dAttToW(1) = .true.   ! PointData
            dAttToW(2) = .true.   ! CellData
            dAttToW(3) = .true.   ! Points
            dAttToW(9) = .true.   ! Cells
            call vtkWriteDataArr(vtk, dAttToW, istat)
            if (istat .lt. 0) return
         case("PolyData")
            dAttToW = .true.
            dAttToW(4) = .false.  ! Coords
            dAttToW(9) = .false.  ! Cells
            call vtkWriteDataArr(vtk, dAttToW, istat)
            if (istat .lt. 0) return
         end select
         write(fid) '    </Piece>'//newl
         write(fid) '  </'//trim(vtk%vtkPcType)//'>'//newl

         rewind(10)
         do
            read(10,end=001) c
            write(fid) c
         end do

 001     close(10)
         close(fid)

!         write(stdout,ftab2) "Success!"
!         write(stdout,ftab1)

         return

         end subroutine vtkWriteToFile

         !==========================================

         subroutine vtkWriteDataArr(vtk, dToW, istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         logical, dimension(nPieceAtts) :: dToW
         integer(IK), intent(inout) :: istat
         integer :: fid,iatt,i,j

         istat = 0
         fid = vtk%fid

         do iatt=1, nPieceAtts
            if ( dToW(iatt) .and. &
               len(trim(vtk%pcAtt(iatt)%pName)).gt.0 ) then
               write(fid) '      <'//trim(libPieceAtts(iatt))
               do i=1, nPieceData
                  if (len(trim(vtk%pcAtt(iatt)%ptClField(i))) .eq. 0) &
                     cycle
                  write(fid) ' '//trim(vtk%pcAtt(iatt)%ptClField(i))// &
                  '="'//trim(vtk%pcAtt(iatt)%ptClFieldName(i))//'" '
               end do
               write(fid) '>'//newl
               do i=1, vtk%pcAtt(iatt)%n
                  write(fid) '        <DataArray '
                  do j=1, nDataElms
                     if ( (len(trim(vtk%pcAtt(iatt)%dataArr(i)% &
                        dElms(j))).gt.0) .or. (j.gt.5) ) then
                        write(fid) trim(libDataElms(j))//'="'// &
                        trim(vtk%pcAtt(iatt)%dataArr(i)%dElms(j))//'" '
                     end if
                  end do
                  write(fid) '/>'//newl
               end do
               write(fid) '      </'//trim(libPieceAtts(iatt))// &
                  '>'//newl
            end if
         end do

         end subroutine vtkWriteDataArr

         !==========================================

         subroutine vtkWriteZlibBinaryData(vtk,dArr,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         type(dataArrType), intent(inout) :: dArr
         integer, parameter :: CHUNK=32768
         integer, parameter :: COMPR=6

         integer :: fid,ist,iend,iBlk
         integer :: np,np1,np2,nlen,ikind
         integer, dimension(:), allocatable  :: hdr
         integer(IK1), dimension(:), allocatable :: p1,p2,pIK1
         integer(IK1), dimension(:), allocatable :: iarrIK1
         integer(IK2), dimension(:), allocatable :: iarrIK2
         integer(IK4), dimension(:), allocatable :: iarrIK4
         integer(IK8), dimension(:), allocatable :: iarrIK8
         real(RK4), dimension(:), allocatable :: rarrRK4
         real(RK8), dimension(:), allocatable :: rarrRK8

         istat = 0
         fid = vtk%fid
         ikind = dArr%ikind
         np = dArr%nElms * ikind
         nlen = np / CHUNK
         if (mod(np, CHUNK).ne.0) nlen = nlen+1
         allocate(pIK1(2*np)); pIK1 = 0_IK1

         allocate(hdr(nlen+3))
         hdr(1) = nlen
         hdr(2) = CHUNK
         hdr(3) = mod(np, CHUNK)
         np1 = 3*dArr%hdrKind
         allocate(p1(np1)); p1 = 0_IK1
         p1 = transfer(hdr(1:3), p1)
         write(10) p1(1:np1)
         deallocate(p1)

         if ( debug ) then
            write(stdout,ftab3) "DataArray: "//trim(dArr%dName)//" ("//&
               trim(dArr%dType)//", kind="//trim(STR(dArr%iKind))//")"
            write(stdout,ftab4) &
               "num blocks "//trim(STR(hdr(1)))
            write(stdout,ftab4) &
               "size of block (before compr.) "// &
               trim(STR(hdr(2)))
            write(stdout,ftab4) &
               "size of last block (before compr.) "// &
               trim(STR(hdr(3)))
         end if

         ist=1; np=0; np1=0
         do iBlk=1, hdr(1)
            np1 = hdr(2)
            if (iBlk.eq.hdr(1) .and. hdr(3).ne.0) np1 = hdr(3)
            allocate(p1(np1)); p1 = 0_IK1
            nlen = np1 / ikind
            iend = ist-1 + nlen
            if (dArr%isInt) then
               select case (dArr%iKind)
               case(IK1)
                  allocate(iarrIK1(nlen))
                  iarrIK1 = int(dArr%iarr(ist:iend), kind=IK1)
                  p1 = transfer(iarrIK1, p1)
                  deallocate(iarrIK1)
               case(IK2)
                  allocate(iarrIK2(nlen))
                  iarrIK2 = int(dArr%iarr(ist:iend), kind=IK2)
                  p1 = transfer(iarrIK2, p1)
                  deallocate(iarrIK2)
               case(IK4)
                  allocate(iarrIK4(nlen))
                  iarrIK4 = int(dArr%iarr(ist:iend), kind=IK4)
                  p1 = transfer(iarrIK4, p1)
                  deallocate(iarrIK4)
               case(IK8)
                  allocate(iarrIK8(nlen))
                  iarrIK8 = int(dArr%iarr(ist:iend), kind=IK8)
                  p1 = transfer(iarrIK8, p1)
                  deallocate(iarrIK8)
               case default
                  write(stdout,ftab4) "ERROR: unknown data type. <"// &
                     trim(dArr%dType)//">"
                  istat=-1; return
               end select
            else
               select case (dArr%iKind)
               case(RK4)
                  allocate(rarrRK4(nlen))
                  rarrRK4 = real(dArr%darr(ist:iend), kind=RK4)
                  p1 = transfer(rarrRK4, p1)
                  deallocate(rarrRK4)
               case(RK8)
                  allocate(rarrRK8(nlen))
                  rarrRK8 = real(dArr%darr(ist:iend), kind=RK8)
                  p1 = transfer(rarrRK8, p1)
                  deallocate(rarrRK8)
               case default
                  write(stdout,ftab4) "ERROR: unknown data type. <"// &
                     trim(dArr%dType)//">"
                  istat=-1; return
               end select
            end if
            np2 = CHUNK
            allocate(p2(np2)); p2 = 0_IK1
            call defZlibData(p1,np1,p2,np2,COMPR,istat)
            if (istat .lt. 0) return
            if ( debug ) then
               write(stdout,ftab4) "size of block#"//trim(STR(iBlk))// &
                  " (after compr.) "//trim(STR(np2))
            end if
            hdr(3+iBlk) = np2
            pIK1(np+1:np+np2) = p2(1:np2)
            np = np + np2
            deallocate(p1,p2)

            ist = iend + 1
         end do

         nlen = hdr(1)
         if (nlen.gt.0) then
            np1 = nlen * dArr%hdrKind
            allocate(p1(np1)); p1 = 0_IK1
            p1 = transfer(hdr(4:nlen+3), p1)
            write(10) p1(1:np1)
            write(10) pIK1(1:np)
            deallocate(p1, pIK1)
         end if

         dArr%nBytes = np1 + np + (3*dArr%hdrKind)

         return

         end subroutine vtkWriteZlibBinaryData

         !==========================================

         subroutine putVTK_pointCoords(vtk,x,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         real(RK), intent(in) :: x(:,:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,j,k,nd,nn,ind

         istat = 0
         nd = size(x,1)
         nn = size(x,2)

         vtk%pieceElms(1) = nn     ! NumberOfPoints
         iatt = 3                  ! Points

         vtk%pcAtt(iatt)%n = 1
         vtk%pcAtt(iatt)%pName = "Points"

         allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))

         i=1
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Float32"
         vtk%pcAtt(iatt)%dataArr(i)%dName = "Points"
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .false.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = MAX(nd,maxNSD)
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = nn
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = MAX(nd,maxNSD)*nn

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = &
            STR(vtk%pcAtt(iatt)%dataArr(i)%nComps)
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(x))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(x))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%darr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0D0
         do j=1, nn
            do k=1, nd
               ind = (j-1)*maxNSD + k
               vtk%pcAtt(iatt)%dataArr(i)%darr(ind) = x(k,j)
            end do
         end do

         return

         end subroutine putVTK_pointCoords

         !==========================================

         subroutine putVTK_elemIEN(vtk,ien,vtkType,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(in) :: ien(:,:), vtkType
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,ne,eNoN,i,j,k,n
         integer(IK), allocatable :: tmpI(:),ioff(:)

         istat = 0
         eNoN = size(ien,1)
         ne   = size(ien,2)

         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            iatt = 9
            allocate(ioff(vtk%pcAtt(iatt)%n))
            ioff = 0
            if ( vtk%pieceElms(2).eq.0 ) then
               vtk%pieceElms(2) = ne
               vtk%pcAtt(iatt)%dataArr(1)%nVals  = ne*eNoN
               vtk%pcAtt(iatt)%dataArr(2)%nVals  = ne
               vtk%pcAtt(iatt)%dataArr(3)%nVals  = ne
               do i=1, vtk%pcAtt(iatt)%n
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  vtk%pcAtt(iatt)%dataArr(i)%nElms = n
                  allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(n))
               end do
            else
               vtk%pieceElms(2) = vtk%pieceElms(2) + ne
               vtk%pcAtt(iatt)%dataArr(1)%nVals  = ne*eNoN + &
                  vtk%pcAtt(iatt)%dataArr(1)%nVals
               vtk%pcAtt(iatt)%dataArr(2)%nVals  = ne + &
                  vtk%pcAtt(iatt)%dataArr(2)%nVals
               vtk%pcAtt(iatt)%dataArr(3)%nVals  = ne + &
                  vtk%pcAtt(iatt)%dataArr(3)%nVals
               do i=1, vtk%pcAtt(iatt)%n
                  n = vtk%pcAtt(iatt)%dataArr(i)%nElms
                  ioff(i) = n
                  allocate(tmpI(n))
                  tmpI = vtk%pcAtt(iatt)%dataArr(i)%iarr
                  deallocate(vtk%pcAtt(iatt)%dataArr(i)%iarr)
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(n))
                  n = vtk%pcAtt(iatt)%dataArr(i)%nElms
                  vtk%pcAtt(iatt)%dataArr(i)%iarr(1:n) = tmpI
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  vtk%pcAtt(iatt)%dataArr(i)%nElms = n
                  deallocate(tmpI)
               end do
            end if
         case("PolyData")
            iatt = 8
            allocate(ioff(vtk%pcAtt(iatt)%n))
            ioff = 0
            if ( vtk%pieceElms(6).eq.0 ) then
               vtk%pieceElms(6) = ne
               vtk%pcAtt(iatt)%n = 2
               vtk%pcAtt(iatt)%pName = "Polys"
               vtk%pcAtt(iatt)%dataArr(1)%nVals  = ne*eNoN
               vtk%pcAtt(iatt)%dataArr(2)%nVals  = ne
               do i=1, vtk%pcAtt(iatt)%n
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  vtk%pcAtt(iatt)%dataArr(i)%nElms = n
                  allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(n))
               end do
            else
               vtk%pieceElms(6) = vtk%pieceElms(6) + ne
               vtk%pcAtt(iatt)%dataArr(1)%nVals = ne*eNoN + &
                  vtk%pcAtt(iatt)%dataArr(1)%nVals
               vtk%pcAtt(iatt)%dataArr(2)%nVals = ne + &
                  vtk%pcAtt(iatt)%dataArr(2)%nVals
               do i=1, vtk%pcAtt(iatt)%n
                  n = vtk%pcAtt(iatt)%dataArr(i)%nElms
                  ioff(i) = n
                  allocate(tmpI(n))
                  tmpI = vtk%pcAtt(iatt)%dataArr(i)%iarr
                  deallocate(vtk%pcAtt(iatt)%dataArr(i)%iarr)
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(n))
                  n = vtk%pcAtt(iatt)%dataArr(i)%nElms
                  vtk%pcAtt(iatt)%dataArr(i)%iarr(1:n) = tmpI
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  vtk%pcAtt(iatt)%dataArr(i)%nElms = n
                  deallocate(tmpI)
               end do
            end if
         case default
            write(stdout,ftab4) &
               "ERROR: cannot find valid vtk piece type. "// &
               "Initialize vtk writer.."
            istat=-1
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end select
         vtk%pcAtt(iatt)%dataArr(:)%nComps = 1

         do i=1, vtk%pcAtt(iatt)%n
            select case(i)
            case(1)
               do j=1, ne
                  do k=1, eNoN
                    n = ioff(i) + (j-1)*eNoN + k
                    vtk%pcAtt(iatt)%dataArr(i)%iarr(n) = ien(k,j)
                  end do
               end do
            case(2)
               do j=1, ne
                  n = ioff(i) + j
                  vtk%pcAtt(iatt)%dataArr(i)%iarr(n) = ioff(1) + j*eNoN
               end do
            case(3)
               do j=1, ne
                  n = ioff(i) + j
                  vtk%pcAtt(iatt)%dataArr(i)%iarr(n) = vtkType
               end do
            end select
         end do

         deallocate(ioff)

         return
         end subroutine putVTK_elemIEN

         !==========================================

         subroutine putVTK_pointDataIntS(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         integer(IK), intent(in) :: u(:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,n
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         n = size(u)
         if ( vtk%pieceElms(1).eq.0 ) then
            vtk%pieceElms(1) = n
         else if ( vtk%pieceElms(1).gt.0 .and. &
                   vtk%pieceElms(1).ne.n ) then
            write(stdout,ftab4) "ERROR: Inconsistent numPoints "// &
               "between POINTS attribute and "//trim(dName)
            istat=-1
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if
         iatt = 1

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "PointData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (trim(vtk%pcAtt(iatt)%ptClField(1)).eq."") then
            vtk%pcAtt(iatt)%ptClField(1) = "Scalars"
            vtk%pcAtt(iatt)%ptClFieldName(1) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Int32"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = 1
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%iarr(:) = u(:)

         return

         end subroutine putVTK_pointDataIntS

         !==========================================

         subroutine putVTK_pointDataIntV(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         integer(IK), intent(in) :: u(:,:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,j,k,m,n,ind
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         m = size(u,1)
         if (m .eq. 1) then
            call putVTK_pointDataIntS(vtk,dName,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(1).eq.0 ) then
            vtk%pieceElms(1) = n
         else if ( vtk%pieceElms(1).gt.0 .and. &
                   vtk%pieceElms(1).ne.n ) then
            write(stdout,ftab4) "ERROR: Inconsistent numPoints "// &
               "between POINTS attribute and "//trim(dName)
            istat=-1
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if
         iatt = 1

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "PointData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (m.le.maxNSD .AND. trim(vtk%pcAtt(iatt)%ptClField(2)).eq."") then
            vtk%pcAtt(iatt)%ptClField(2) = "Vectors"
            vtk%pcAtt(iatt)%ptClFieldName(2) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Int32"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = MAX(m, maxNSD)
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = MAX(m,maxNSD)*n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         if (vtk%pcAtt(iatt)%dataArr(i)%nComps .gt. 1) &
            vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = &
               STR(vtk%pcAtt(iatt)%dataArr(i)%nComps)
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%iarr = 0
         do j=1, n
            do k=1, m
               ind = (j-1)*maxNSD + k
               vtk%pcAtt(iatt)%dataArr(i)%iarr(ind) = u(k,j)
            end do
         end do

         return

         end subroutine putVTK_pointDataIntV

         !==========================================

         subroutine putVTK_pointDataRealS(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         real(RK), intent(in) :: u(:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,n
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         n = size(u)
         if ( vtk%pieceElms(1).eq.0 ) then
            vtk%pieceElms(1) = n
         else if ( vtk%pieceElms(1).gt.0 .and. &
                   vtk%pieceElms(1).ne.n ) then
            write(stdout,ftab4) "ERROR: Inconsistent numPoints "// &
               "between POINTS attribute and "//trim(dName)
            istat=-1
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if
         iatt = 1

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "PointData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (trim(vtk%pcAtt(iatt)%ptClField(1)).eq."") then
            vtk%pcAtt(iatt)%ptClField(1) = "Scalars"
            vtk%pcAtt(iatt)%ptClFieldName(1) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Float64"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = 1
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%darr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%darr(:) = u(:)

         return

         end subroutine putVTK_pointDataRealS

         !==========================================

         subroutine putVTK_pointDataRealV(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         real(RK), intent(in) :: u(:,:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,j,k,m,n,ind
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         m = size(u,1)
         if (m .eq. 1) then
            call putVTK_pointDataRealS(vtk,dName,u(1,:),istat)
            return
         end if

         n = size(u,2)
         if ( vtk%pieceElms(1).eq.0 ) then
            vtk%pieceElms(1) = n
         else if ( vtk%pieceElms(1).gt.0 .and. &
                   vtk%pieceElms(1).ne.n ) then
            write(stdout,ftab4) "ERROR: Inconsistent numPoints "// &
               "between POINTS attribute and "//trim(dName)
            istat=-1
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if
         iatt = 1

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "PointData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (m.le.maxNSD .AND. trim(vtk%pcAtt(iatt)%ptClField(2)).eq."") then
            vtk%pcAtt(iatt)%ptClField(2) = "Vectors"
            vtk%pcAtt(iatt)%ptClFieldName(2) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Float64"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = MAX(m,maxNSD)
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = MAX(m,maxNSD)*n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         if (vtk%pcAtt(iatt)%dataArr(i)%nComps .gt. 1) &
            vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = &
               STR(vtk%pcAtt(iatt)%dataArr(i)%nComps)
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
               vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%darr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0D0
         do j=1, n
            do k=1, m
               ind = (j-1)*maxNSD + k
               vtk%pcAtt(iatt)%dataArr(i)%darr(ind) = u(k,j)
            end do
         end do

         return

         end subroutine putVTK_pointDataRealV

         !==========================================

         subroutine putVTK_elemDataIntS(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         integer(IK), intent(in) :: u(:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,n
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         n = size(u)

         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            if ( vtk%pieceElms(2).eq.0 ) then
               vtk%pieceElms(2) = n
            else if ( vtk%pieceElms(2).gt.0 .and. &
                      vtk%pieceElms(2).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numCells "// &
                  "between CELLS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         case("PolyData")
            if ( vtk%pieceElms(6).eq.0 ) then
               vtk%pieceElms(6) = n
            else if ( vtk%pieceElms(6).gt.0 .and. &
                      vtk%pieceElms(6).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numPolys "// &
                  "between POLYS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         end select
         iatt = 2

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "CellData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (trim(vtk%pcAtt(iatt)%ptClField(1)).eq."") then
            vtk%pcAtt(iatt)%ptClField(1) = "Scalars"
            vtk%pcAtt(iatt)%ptClFieldName(1) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Int32"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = 1
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%iarr(:) = u(:)

         return

         end subroutine putVTK_elemDataIntS

         !==========================================

         subroutine putVTK_elemDataIntV(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         integer(IK), intent(in) :: u(:,:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,j,k,m,n,ind
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         m = size(u,1)
         if (m .eq. 1) then
            call putVTK_elemDataIntS(vtk,dName,u(1,:),istat)
            return
         end if

         n = size(u,2)
         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            if ( vtk%pieceElms(2).eq.0 ) then
               vtk%pieceElms(2) = n
            else if ( vtk%pieceElms(2).gt.0 .and. &
                      vtk%pieceElms(2).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numCells "// &
                  "between CELLS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         case("PolyData")
            if ( vtk%pieceElms(6).eq.0 ) then
               vtk%pieceElms(6) = n
            else if ( vtk%pieceElms(6).gt.0 .and. &
                      vtk%pieceElms(6).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numPolys "// &
                  "between POLYS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         end select
         iatt = 2

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "CellData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (m.le.maxNSD .AND. trim(vtk%pcAtt(iatt)%ptClField(2)).eq."") then
            vtk%pcAtt(iatt)%ptClField(2) = "Vectors"
            vtk%pcAtt(iatt)%ptClFieldName(2) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Int32"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = MAX(m,maxNSD)
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = MAX(m,maxNSD)*n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         if (vtk%pcAtt(iatt)%dataArr(i)%nComps .gt. 1) &
            vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = &
               STR(vtk%pcAtt(iatt)%dataArr(i)%nComps)
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
               vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%iarr = 0
         do j=1, n
            do k=1, m
               ind = (j-1)*maxNSD + k
               vtk%pcAtt(iatt)%dataArr(i)%iarr(ind) = u(k,j)
            end do
         end do

         return

         end subroutine putVTK_elemDataIntV

         !==========================================

         subroutine putVTK_elemDataRealS(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         real(RK), intent(in) :: u(:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,n
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         n = size(u)

         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            if ( vtk%pieceElms(2).eq.0 ) then
               vtk%pieceElms(2) = n
            else if ( vtk%pieceElms(2).gt.0 .and. &
                      vtk%pieceElms(2).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numCells "// &
                  "between CELLS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         case("PolyData")
            if ( vtk%pieceElms(6).eq.0 ) then
               vtk%pieceElms(6) = n
            else if ( vtk%pieceElms(6).gt.0 .and. &
                      vtk%pieceElms(6).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numPolys "// &
                  "between POLYS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         end select
         iatt = 2

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "CellData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (trim(vtk%pcAtt(iatt)%ptClField(1)).eq."") then
            vtk%pcAtt(iatt)%ptClField(1) = "Scalars"
            vtk%pcAtt(iatt)%ptClFieldName(1) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Float64"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = 1
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = ""
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%darr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%darr(:) = u(:)

         return

         end subroutine putVTK_elemDataRealS

         !==========================================

         subroutine putVTK_elemDataRealV(vtk,dName,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: dName
         real(RK), intent(in) :: u(:,:)
         integer(IK), intent(inout) :: istat

         integer(IK) :: iatt,i,j,k,m,n,ind
         type(dataArrType), dimension(:), allocatable :: dArr

         istat = 0
         m = size(u,1)
         if (m .eq. 1) then
            call putVTK_elemDataRealS(vtk,dName,u(1,:),istat)
            return
         end if

         n = size(u,2)
         select case(trim(vtk%vtkPcType))
         case("UnstructuredGrid")
            if ( vtk%pieceElms(2).eq.0 ) then
               vtk%pieceElms(2) = n
            else if ( vtk%pieceElms(2).gt.0 .and. &
                      vtk%pieceElms(2).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numCells "// &
                  "between CELLS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         case("PolyData")
            if ( vtk%pieceElms(6).eq.0 ) then
               vtk%pieceElms(6) = n
            else if ( vtk%pieceElms(6).gt.0 .and. &
                      vtk%pieceElms(6).ne.n ) then
               write(stdout,ftab4) "ERROR: Inconsistent numPolys "// &
                  "between POLYS attribute and "//trim(dName)
               istat=-1
               inquire(unit=10, opened=flag)
               if (flag) close(10)
               return
            end if
         end select
         iatt = 2

         if ( allocated(vtk%pcAtt(iatt)%dataArr) ) then
            allocate(dArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n
               call vtkDeepCopyData(vtk%pcAtt(iatt)%dataArr(i), dArr(i))
               call vtkFlushData(vtk%pcAtt(iatt)%dataArr(i))
            end do
            deallocate(vtk%pcAtt(iatt)%dataArr)
            vtk%pcAtt(iatt)%n = vtk%pcAtt(iatt)%n + 1
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
            do i=1, vtk%pcAtt(iatt)%n-1
               call vtkDeepCopyData(dArr(i), vtk%pcAtt(iatt)%dataArr(i))
               call vtkFlushData(dArr(i))
            end do
            deallocate(dArr)
         else
            vtk%pcAtt(iatt)%n = 1
            vtk%pcAtt(iatt)%pName = "CellData"
            allocate(vtk%pcAtt(iatt)%dataArr(vtk%pcAtt(iatt)%n))
         end if

         if (m.le.maxNSD .AND. trim(vtk%pcAtt(iatt)%ptClField(2)).eq."") then
            vtk%pcAtt(iatt)%ptClField(2) = "Vectors"
            vtk%pcAtt(iatt)%ptClFieldName(2) = trim(dName)
         end if

         i = vtk%pcAtt(iatt)%n
         vtk%pcAtt(iatt)%dataArr(i)%dType = "Float64"
         vtk%pcAtt(iatt)%dataArr(i)%dName = trim(dName)
         vtk%pcAtt(iatt)%dataArr(i)%dFrmt = "appended"
         vtk%pcAtt(iatt)%dataArr(i)%hdrType = "UInt32"

         vtk%pcAtt(iatt)%dataArr(i)%isInt = .true.

         vtk%pcAtt(iatt)%dataArr(i)%nComps = MAX(m,maxNSD)
         vtk%pcAtt(iatt)%dataArr(i)%nVals  = n
         vtk%pcAtt(iatt)%dataArr(i)%nElms  = MAX(m,maxNSD)*n

         vtk%pcAtt(iatt)%dataArr(i)%dElms(1) = &
            vtk%pcAtt(iatt)%dataArr(i)%dType
         vtk%pcAtt(iatt)%dataArr(i)%dElms(2) = &
            vtk%pcAtt(iatt)%dataArr(i)%dName
         if (vtk%pcAtt(iatt)%dataArr(i)%nComps .gt. 1) &
            vtk%pcAtt(iatt)%dataArr(i)%dElms(3) = &
               STR(vtk%pcAtt(iatt)%dataArr(i)%nComps)
         vtk%pcAtt(iatt)%dataArr(i)%dElms(4) = &
            vtk%pcAtt(iatt)%dataArr(i)%dFrmt
         vtk%pcAtt(iatt)%dataArr(i)%dElms(6) = STR(MINVAL(u))
         vtk%pcAtt(iatt)%dataArr(i)%dElms(7) = STR(MAXVAL(u))

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%hdrType, &
            vtk%pcAtt(iatt)%dataArr(i)%hdrKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         call selectDataType(vtk%pcAtt(iatt)%dataArr(i)%dType, &
            vtk%pcAtt(iatt)%dataArr(i)%iKind, &
            vtk%pcAtt(iatt)%dataArr(i)%isInt, istat)
         if (istat .lt. 0) then
            inquire(unit=10, opened=flag)
            if (flag) close(10)
            return
         end if

         allocate(vtk%pcAtt(iatt)%dataArr(i)%darr(vtk%pcAtt(iatt)%dataArr(i)%nElms))
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0D0
         do j=1, n
            do k=1, m
               ind = (j-1)*maxNSD + k
               vtk%pcAtt(iatt)%dataArr(i)%darr(ind) = u(k,j)
            end do
         end do

         return

         end subroutine putVTK_elemDataRealV

         !==========================================

      end module vtkXMLMod

!**************************************************


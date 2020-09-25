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
!     Module to read/write VTU/VTP/VTI files.
!
!--------------------------------------------------------------------

      module vtkXMLlib
      use dataTypeParams
      use stdParams, only : strL

      logical, parameter :: debug = .false.
      character(len=64), parameter :: b64List = &
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"

      integer(IK), parameter :: nVTKElms=5
      integer(IK), parameter :: nVTKDataTyps=5
      integer(IK), parameter :: nVTKDataAtts=3
      integer(IK), parameter :: nPieceElms=7
      integer(IK), parameter :: nPieceData=5
      integer(IK), parameter :: nPieceAtts=9
      integer(IK), parameter :: nDataElms=7
      integer(IK), parameter :: nDataTyps=10
      integer(IK), parameter :: nDataFrmt=3
      integer(IK), parameter :: nDataEnc=2

      character(len=strL), dimension(nVTKElms)     :: libVTKElms
      character(len=strL), dimension(nVTKDataTyps) :: libVTKDataTyps
      character(len=strL), dimension(nVTKDataAtts) :: libVTKDataAtts
      character(len=strL), dimension(nPieceElms)   :: libPieceElms
      character(len=strL), dimension(nPieceData)   :: libPcPtClData
      character(len=strL), dimension(nPieceAtts)   :: libPieceAtts
      character(len=strL), dimension(nDataElms)    :: libDataElms
      character(len=strL), dimension(nDataTyps)    :: libDataTyps
      character(len=strL), dimension(nDataFrmt)    :: libDataFrmt
      character(len=strL), dimension(nDataFrmt)    :: libDataEnc

      contains

         subroutine initVTKXMLlib
         implicit none

         libVTKElms=""
         libVTKDataTyps=""
         libVTKDataAtts=""
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

         ! libVTKDataTyps(5) !
         libVTKDataTyps(1) = "ImageData"
         libVTKDataTyps(2) = "RectilinearGrid"
         libVTKDataTyps(3) = "StructuredGrid"
         libVTKDataTyps(4) = "PolyData"
         libVTKDataTyps(5) = "UnstructuredGrid"

         ! libVTKDataAtts(3) !
         libVTKDataAtts(1) = "WholeExtent"
         libVTKDataAtts(2) = "Origin"
         libVTKDataAtts(3) = "Spacing"

         ! libPieceElms(7) !
         libPieceElms(1) = "NumberOfPoints"
         libPieceElms(2) = "NumberOfCells"
         libPieceElms(3) = "NumberOfVerts"
         libPieceElms(4) = "NumberOfLines"
         libPieceElms(5) = "NumberOfStrips"
         libPieceElms(6) = "NumberOfPolys"
         libPieceElms(7) = "Extent"

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

      public :: getVTK_numElems
      public :: getVTK_numPoints
      public :: getVTK_nodesPerElem
      public :: putVTK_numElems
      public :: putVTK_numPoints

      public :: getVTK_pointCoords
      public :: getVTK_elemIEN
      public :: putVTK_pointCoords
      public :: putVTK_elemIEN

      public :: getVTK_numPointData
      public :: getVTK_pointDataNames
      public :: getVTK_pointData

      public :: getVTK_numElemData
      public :: getVTK_elemDataNames
      public :: getVTK_elemData
      public :: putVTK_pointData
      public :: putVTK_elemData

      public :: getVTK_imageExtent
      public :: getVTK_imageOrigin
      public :: getVTK_imageSpacing
      public :: getVTK_pieceExtent
      public :: putVTK_imageExtent
      public :: putVTK_imageOrigin
      public :: putVTK_imageSpacing
      public :: putVTK_pieceExtent

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
      character(len=strL), dimension(maxToks) :: tokenList
      character :: c
      integer(IK) :: itok,ntoks,slen,iPos
      integer(IK) :: maxRank
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

      type vtkPcElmType
         private
         integer(IK) :: nPoints
         integer(IK) :: nCells
         integer(IK) :: nVerts
         integer(IK) :: nLines
         integer(IK) :: nStrips
         integer(IK) :: nPolys
         integer(IK) :: extent(6)
      end type vtkPcElmType

      type vtkDataType
         private
         character(len=strL) :: str
         integer(IK) :: extent(6)
         real(RK) :: origin(3)
         real(RK) :: spacng(3)
      end type vtkDataType

      type vtkXMLType
         private
         logical :: isBinApp
         character(len=strL) :: fileName
         character(len=strL) :: dataFormat
         character(len=strL) :: dataEncdng
         integer(IK) :: fid
         integer(IK) :: stAppendPos
         integer(IK) :: endAppendPos
         integer(IK) :: offsets(100)

         character(len=strL), dimension(nVTKElms)   :: vtkElms

         type(vtkDataType)  :: dataType
         type(vtkPcElmType) :: pcElms
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
            write(stdout,ftab4) "ERROR: File "//trim(fName)// &
     &         " does not exist"
            istat=-1; return
         end if

         slen = len(trim(fName))
         select case (fName(slen-2:slen))
         case ("vtu","vtp","vti")
         case default
            write(stdout,ftab4) "ERROR: unknown file extension "//&
     &         "(can only be vtu or vtp or vti)"
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
         vtk%dataType%str    = ""
         vtk%stAppendPos     = 0
         vtk%endAppendPos    = 0
         vtk%offsets(:)      = 0

         vtk%vtkElms(:)      = ""

         vtk%dataType%str    = ""
         vtk%dataType%extent = 0
         vtk%dataType%origin = 0.0_RK
         vtk%dataType%spacng = 0.0_RK

         vtk%pcElms%nPoints  = 0
         vtk%pcElms%nCells   = 0
         vtk%pcElms%nVerts   = 0
         vtk%pcElms%nLines   = 0
         vtk%pcElms%nStrips  = 0
         vtk%pcElms%nPolys   = 0
         vtk%pcElms%extent   = 0

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

         select case (trim(vtk%vtkElms(1)))
         case ("ImageData")
            vtk%dataType%str = libVTKDataTyps(1)
            call loadVTKDataAtts(vtk,istat)

         case ("RectilinearGrid")
            vtk%dataType%str = libVTKDataTyps(2)
            call loadVTKDataAtts(vtk,istat)

         case ("StructuredGrid")
            vtk%dataType%str = libVTKDataTyps(3)
            call loadVTKDataAtts(vtk,istat)

         case ("PolyData")
            vtk%dataType%str = libVTKDataTyps(4)

         case ("UnstructuredGrid")
            vtk%dataType%str = libVTKDataTyps(5)
         end select

         if ( debug .and. len(trim(vtk%vtkElms(5))).gt.0 ) &
            write(stdout,ftab2) &
               "Data compression: "//trim(vtk%vtkElms(5))

         ! Piece elements parser !
         call loadVTKPieceElms(vtk,istat)

         return
         end subroutine parseVTKKernel

         !==========================================

         subroutine loadVTKDataAtts(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer :: i

         rLine = ""
         stmp = "<"//trim(vtk%dataType%str)
         tokenList(:) = ""

         call findKwrdXML(vtk,stmp,iPos,rLine,istat)
         if ( istat.lt.0 ) return
         if ( debug ) write(stdout,ftab2) trim(rLine)
         call parseString(rLine,tokenList,ntoks)

         stmp = ""
         do i=1, 6
            stmp = getTokenValue(tokenList,ntoks,libVTKDataAtts(1), i)
            slen = len(trim(stmp))
            if (slen.gt.0) then
               read(stmp(1:slen),*) vtk%dataType%extent(i)
            end if
         end do

         if (trim(vtk%dataType%str) .eq. libVTKDataTyps(2) .or. &
             trim(vtk%dataType%str) .eq. libVTKDataTyps(3)) return

         stmp = ""
         do i=1, 3
            stmp = getTokenValue(tokenList,ntoks,libVTKDataAtts(2), i)
            slen = len(trim(stmp))
            if (slen.gt.0) then
               read(stmp(1:slen),*) vtk%dataType%origin(i)
            end if
         end do

         stmp = ""
         do i=1, 3
            stmp = getTokenValue(tokenList,ntoks,libVTKDataAtts(3), i)
            slen = len(trim(stmp))
            if (slen.gt.0) then
               read(stmp(1:slen),*) vtk%dataType%spacng(i)
            end if
         end do

         return
         end subroutine loadVTKDataAtts

         !==========================================

         subroutine loadVTKPieceElms(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer :: i

         rLine = ""
         stmp = ""
         tokenList(:) = ""

         call findKwrdXML(vtk,"<Piece",iPos,rLine,istat)
         if ( istat.lt.0 ) return
         if ( debug ) write(stdout,ftab2) trim(rLine)
         call parseString(rLine,tokenList,ntoks)

         do itok=1, nPieceElms
            if (itok .lt. 7) then
               stmp = getTokenValue(tokenList,ntoks,libPieceElms(itok))
               slen = len(trim(stmp))
               if ( slen.gt.0 ) then
                  select case (trim(libPieceElms(itok)))
                  case ("NumberOfPoints")
                     read(stmp(1:slen),*) vtk%pcElms%nPoints
                  case ("NumberOfCells")
                     read(stmp(1:slen),*) vtk%pcElms%nCells
                  case ("NumberOfVerts")
                     read(stmp(1:slen),*) vtk%pcElms%nVerts
                  case ("NumberOfLines")
                     read(stmp(1:slen),*) vtk%pcElms%nLines
                  case ("NumberOfStrips")
                     read(stmp(1:slen),*) vtk%pcElms%nStrips
                  case ("NumberOfPolys")
                     read(stmp(1:slen),*) vtk%pcElms%nPolys
                  end select
               end if
            else
               do i=1, 6
                  stmp = getTokenValue(tokenList,ntoks,libPieceElms(itok),i)
                  slen = len(trim(stmp))
                  if ( slen.gt.0 ) then
                     read(stmp(1:slen),*) vtk%pcElms%extent(i)
                  end if
               end do
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

         call loadPieceAttributes(vtk,istat)

         return
         end subroutine loadVTKPieceElms

         !==========================================

         subroutine loadPieceAttributes(vtk,istat)
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

            end do ! i
         end do ! outer loop over Piece atts !
         maxRank = cntr

         return
         end subroutine loadPieceAttributes

         !==========================================

         subroutine vtkDataLoader(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,rank

         if ( debug ) write(stdout,ftab2) &
            "VTK Data Type: "//trim(vtk%dataType%str)
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

               call vtkXMLDataParser(vtk,vtk%pcAtt(iatt)%dataArr(i),iatt,istat)
               if ( istat.lt.0 ) return
            end do
            if (debug) write(stdout,ftab1) repeat("*",76)
         end do

         return

 001     write(stdout,ftab4) "ERROR: end of file reached.."
         istat = -1; return

         end subroutine vtkDataLoader

         !==========================================

         subroutine vtkXMLDataParser(vtk,dataArr,iatt,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt
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
               "Piece Type: "//trim(vtk%dataType%str)
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

         select case (trim(vtk%dataType%str))
         case ("ImageData")
            call vtkParseImageData(vtk,dataArr,iatt,istat)
            if ( istat.lt.0 ) return

         case ("RectilinearGrid")
            write(stdout,ftab4) "ERROR: Rectilinear grid parser is "//&
     &         "not available yet.."
!            call vtkParseRectGrid(vtk,dataArr,iatt,istat)
            if ( istat.lt.0 ) return

         case ("StructuredGrid")
            call vtkParseStrucGrid(vtk,dataArr,iatt,istat)
            if ( istat.lt.0 ) return

         case ("UnstructuredGrid")
            call vtkParseUnstrucGrid(vtk,dataArr,iatt,istat)
            if ( istat.lt.0 ) return

         case ("PolyData")
            call vtkParsePolyData(vtk,dataArr,iatt,istat)
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

         subroutine vtkParseImageData(vtk,dataArr,iatt,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nCells

         dPtr => dataArr

         dPtr%isInt = .false.
         call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)
         if ( istat.lt.0 ) return

         nPoints = (vtk%pcElms%extent(2)-vtk%pcElms%extent(1)+1) * &
                   (vtk%pcElms%extent(4)-vtk%pcElms%extent(3)+1) * &
                   (vtk%pcElms%extent(6)-vtk%pcElms%extent(5)+1)

         nCells  = (vtk%pcElms%extent(2)-vtk%pcElms%extent(1)) * &
                   (vtk%pcElms%extent(4)-vtk%pcElms%extent(3)) * &
                   (vtk%pcElms%extent(6)-vtk%pcElms%extent(5))

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

         end select

         end subroutine vtkParseImageData

         !==========================================

!         subroutine vtkParseRectGrid(vtk,dataArr,iatt,istat)
!         implicit none
!         type(vtkXMLType), intent(inout) :: vtk
!         type(dataArrType), intent(inout), target :: dataArr
!         integer(IK), intent(in) :: iatt
!         integer(IK), intent(inout) :: istat

!         end subroutine vtkParseRectGrid

         !==========================================

         subroutine vtkParseStrucGrid(vtk,dataArr,iatt,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nCells

         dPtr => dataArr

         dPtr%isInt = .false.
         call selectDataType(dPtr%dType,dPtr%ikind,dPtr%isInt,istat)
         if ( istat.lt.0 ) return

         nPoints = (vtk%pcElms%extent(2)-vtk%pcElms%extent(1)+1) * &
                   (vtk%pcElms%extent(4)-vtk%pcElms%extent(3)+1) * &
                   (vtk%pcElms%extent(6)-vtk%pcElms%extent(5)+1)

         nCells  = (vtk%pcElms%extent(2)-vtk%pcElms%extent(1)) * &
                   (vtk%pcElms%extent(4)-vtk%pcElms%extent(3)) * &
                   (vtk%pcElms%extent(6)-vtk%pcElms%extent(5))

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
               write(stdout,ftab4) "WARNING: Element "// &
     &            "<NumberOfComponents> in <Points> attribute < 3"
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

         end select

         end subroutine vtkParseStrucGrid

         !==========================================

         subroutine vtkParseUnstrucGrid(vtk,dataArr,iatt,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nCells

         dPtr => dataArr

         nPoints = vtk%pcElms%nPoints
         nCells  = vtk%pcElms%nCells
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
               write(stdout,ftab4) "WARNING: Element "// &
     &            "<NumberOfComponents> in <Points> attribute < 3"
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

         subroutine vtkParsePolyData(vtk,dataArr,iatt,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         type(dataArrType), intent(inout), target :: dataArr
         integer(IK), intent(in) :: iatt
         integer(IK), intent(inout) :: istat
         type(dataArrType), pointer :: dPtr
         integer(IK) :: nPoints,nVerts,nLines,nStrips,nPolys

         dPtr => dataArr

         nPoints = vtk%pcElms%nPoints
         nVerts  = vtk%pcElms%nVerts
         nLines  = vtk%pcElms%nLines
         nStrips = vtk%pcElms%nStrips
         nPolys  = vtk%pcElms%nPolys

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
               write(stdout,ftab4) "WARNING: Element "// &
     &            "<NumberOfComponents> in <Points> attribute < 3"
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
                  "WARNING: Typecasting from INT64 to INT32 "// &
                  "could lead to errors.."
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
         integer(IK) :: c,e,Nb

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

         nn = vtk%pcElms%nPoints
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

         if ( vtk%pcElms%nCells.gt.0 ) then
            ne = vtk%pcElms%nCells
         else if ( vtk%pcElms%nPolys.gt.0 ) then
            ne = vtk%pcElms%nPolys
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

         if ( vtk%pcElms%nCells.gt.0 ) then   ! nCells !
            iatt = 9
         else if (vtk%pcElms%nPolys.gt.0 ) then ! nPolys !
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
               itmp = 0
               if (vtk%pcAtt(iatt)%dataArr(i)%nelms .gt. 1) then
                  itmp = vtk%pcAtt(iatt)%dataArr(i)%iarr(2) - &
                        vtk%pcAtt(iatt)%dataArr(i)%iarr(1)
               end if
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

         subroutine getVTK_pointCoords(vtk,x,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         real(RK), intent(out) :: x(:,:)
         integer(IK) :: iatt,i,j,k,l,nd,nn

         nd = size(x,1)
         nn = size(x,2)
         if ( vtk%pcElms%nPoints.eq.nn ) then   ! nPoints !
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

         subroutine getVTK_elemIEN(vtk,ien,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(out) :: ien(:,:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,ne,eNoN,i,j,k,l

         eNoN = size(ien,1)
         ne   = size(ien,2)
         if ( vtk%pcElms%nCells.gt.0 ) then   ! nCells !
            iatt = 9
         else if (vtk%pcElms%nPolys.gt.0 ) then ! nPolys !
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
               "ERROR: could not find connectivity in "// &
               "Cells or Polys attributes"
            istat=-1; return
         end if

         end subroutine getVTK_elemIEN

         !==========================================

         subroutine getVTK_imageExtent(vtk,iLims,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: iLims(2,3)
         integer(IK) :: i
         logical :: l1

         l1 = .true.
         do i=1, 6
            if (vtk%dataType%extent(i) .ne. 0) then
               l1 = .false.
               exit
            end if
         end do
         if (l1) then
            write(stdout,ftab4) &
               "ERROR: could not find WholeExtent attribute"
            istat=-1; return
         end if
         iLims = reshape(vtk%dataType%extent,(/2,3/))

         return
         end subroutine getVTK_imageExtent

         !==========================================

         subroutine getVTK_imageOrigin(vtk,orig,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         real(RK), intent(out) :: orig(maxNSD)

         istat = 0
         orig = vtk%dataType%origin

         return
         end subroutine getVTK_imageOrigin

         !==========================================

         subroutine getVTK_imageSpacing(vtk,dx,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         real(RK), intent(out) :: dx(maxNSD)

         istat = 0
         dx = vtk%dataType%spacng

         return
         end subroutine getVTK_imageSpacing

         !==========================================

         subroutine getVTK_pieceExtent(vtk,pLims,istat)
         implicit none
         type(vtkXMLType), intent(in) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK), intent(out) :: pLims(2,3)
         integer(IK) :: i
         logical :: l1

         l1 = .true.
         do i=1, 6
            if (vtk%pcElms%extent(i) .ne. 0) then
               l1 = .false.
               exit
            end if
         end do
         if (l1) then
            write(stdout,ftab4) &
               "ERROR: could not find WholeExtent attribute"
            istat=-1; return
         end if
         pLims = reshape(vtk%pcElms%extent,(/2,3/))

         return
         end subroutine getVTK_pieceExtent

         !==========================================

         subroutine getVTK_numPointData(vtk,nPtData,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(out) :: nPtData
         integer(IK), intent(inout) :: istat

         nPtData = vtk%pcAtt(1)%n

         return
         end subroutine getVTK_numPointData

         !==========================================

         subroutine getVTK_pointDataNames(vtk,names,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(inout) :: names(:)
         integer(IK), intent(inout) :: istat
         integer(IK) i, n

         n = size(names)
         if (n .ne. vtk%pcAtt(1)%n) then
            write(stdout,ftab4) "ERROR: not enough size of input "//&
              "array for point data names"
            istat = -1; return
         end if

         do i=1, vtk%pcAtt(1)%n
            names(i) = trim(vtk%pcAtt(1)%dataArr(i)%dName)
         end do

         return
         end subroutine getVTK_pointDataNames

         !==========================================

         subroutine getVTK_pointDataIntS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         iatt = 1
         n    = size(u)

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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in PointData attribute"
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

         iatt = 1
         m    = size(u,1)
         n    = size(u,2)
         if ( m.eq.1 ) then
            call getVTK_pointDataIntS(vtk,kwrd,u(1,:),istat)
            return
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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in PointData attribute"
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

         iatt = 1
         n    = size(u)

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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in PointData attribute"
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

         iatt = 1
         m    = size(u,1)
         n    = size(u,2)
         if ( m.eq.1 ) then
            call getVTK_pointDataRealS(vtk,kwrd,u(1,:),istat)
            return
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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in PointData attribute"
            istat=-1; return
         end if

         end subroutine getVTK_pointDataRealV

         !==========================================

         subroutine getVTK_numElemData(vtk,nElData,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(out) :: nElData
         integer(IK), intent(inout) :: istat

         nElData = vtk%pcAtt(2)%n

         return
         end subroutine getVTK_numElemData

         !==========================================

         subroutine getVTK_elemDataNames(vtk,names,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(inout) :: names(:)
         integer(IK), intent(inout) :: istat
         integer(IK) i, n

         n = size(names)
         if (n .ne. vtk%pcAtt(2)%n) then
            write(stdout,ftab4) "ERROR: not enough size of input "//&
              "array for point data names"
            istat = -1; return
         end if

         do i=1, vtk%pcAtt(2)%n
            names(i) = trim(vtk%pcAtt(2)%dataArr(i)%dName)
         end do

         return
         end subroutine getVTK_elemDataNames

         !==========================================

         subroutine getVTK_elemDataIntS(vtk,kwrd,u,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         character(len=*), intent(in) :: kwrd
         integer(IK), intent(inout) :: u(:)
         integer(IK), intent(inout) :: istat
         integer(IK) :: iatt,i,n
         logical :: flag

         iatt = 2
         n    = size(u)

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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in CellData attribute"
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

         iatt = 2
         m    = size(u,1)
         n    = size(u,2)
         if ( m.eq.1 ) then
            call getVTK_elemDataIntS(vtk,kwrd,u(1,:),istat)
            return
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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in CellData attribute"
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

         iatt = 2
         n    = size(u)

         flag = .false.
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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in CellData attribute"
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

         iatt = 2
         m    = size(u,1)
         n    = size(u,2)
         if ( m.eq.1 ) then
            call getVTK_elemDataRealS(vtk,kwrd,u(1,:),istat)
            return
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
               "ERROR: could not find <"//trim(kwrd)// &
               "> in CellData attribute"
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
            vtk%dataType%str = "UnstructuredGrid"
            call vtkInitUGridWriter(vtk, istat)
         case ("vtp")
            vtk%dataType%str = "PolyData"
            call vtkInitPDataWriter(vtk, istat)
         case ("vti")
            vtk%dataType%str    = "ImageData"
            vtk%dataType%extent = 0
            vtk%dataType%origin = 0.0_RK
            vtk%dataType%spacng = 0.0_RK
            call vtkInitImageWriter(vtk, istat)
         case default
            write(stdout,ftab4) "ERROR: unknown file extension "//&
     &         "(can only be vtu or vtp or vti)"
            istat=-1; return
         end select

         vtk%vtkElms(1) = vtk%dataType%str
         vtk%vtkElms(2) = "0.1"
         vtk%vtkElms(3) = "LittleEndian"
         vtk%vtkElms(4) = ""
         vtk%vtkElms(5) = "vtkZLibDataCompressor"
         vtk%pcElms%nPoints  = 0
         vtk%pcElms%nCells   = 0
         vtk%pcElms%nVerts   = 0
         vtk%pcElms%nLines   = 0
         vtk%pcElms%nStrips  = 0
         vtk%pcElms%nPolys   = 0
         vtk%pcElms%extent   = 0

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

         subroutine vtkInitImageWriter(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer, intent(inout) :: istat
         integer :: iatt

         istat = 0
         do iatt=1, nPieceAtts
            vtk%pcAtt(iatt)%pName = ""
            vtk%pcAtt(iatt)%n = 0
            vtk%pcAtt(iatt)%ptClField(:) = ""
            vtk%pcAtt(iatt)%ptClFieldName(:) = ""
         end do

         end subroutine vtkInitImageWriter

         !==========================================

         subroutine vtkWriteToFile(vtk,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(inout) :: istat
         integer(IK) :: fid,i,iatt,ivar,tmpI(100)
         logical, dimension(nPieceAtts) :: dAttToW
         character c
         logical l1

         if ( debug )write(stdout,ftab1) &
            "<VTK XML Writer> Writing to file ----->  "// &
            trim(vtk%fileName)

         istat = 0
         fid   = vtk%fid
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

         write(fid) '  <'//trim(vtk%dataType%str)
         if (trim(vtk%dataType%str) .eq. "ImageData") then
            write(fid) ' '//'WholeExtent="'
            do i=1, 6
               if (i.ne.1) write(fid) ' '
               write(fid) trim(STR(vtk%dataType%extent(i)))
            end do
            write(fid) '" '//'Origin="'
            do i=1, 3
               if (i.ne.1) write(fid) ' '
               write(fid) trim(STR(vtk%dataType%origin(i)))
            end do
            write(fid) '" '//'Spacing="'
            do i=1, 3
               if (i.ne.1) write(fid) ' '
               write(fid) trim(STR(vtk%dataType%spacng(i)))
            end do
            write(fid) '"'
         end if
         write(fid) '>'//newl

         write(fid) '    <Piece '
         if (vtk%pcElms%nPoints .gt. 0) then
            write(fid) trim(libPieceElms(1))//'="'// &
               trim(STR(vtk%pcElms%nPoints))//'"    '
         end if
         if (vtk%pcElms%nCells .gt. 0) then
            write(fid) trim(libPieceElms(2))//'="'// &
               trim(STR(vtk%pcElms%nCells))//'"    '
         end if
         if (vtk%pcElms%nVerts .gt. 0) then
            write(fid) trim(libPieceElms(3))//'="'// &
               trim(STR(vtk%pcElms%nVerts))//'"    '
         end if
         if (vtk%pcElms%nLines .gt. 0) then
            write(fid) trim(libPieceElms(4))//'="'// &
               trim(STR(vtk%pcElms%nLines))//'"    '
         end if
         if (vtk%pcElms%nStrips .gt. 0) then
            write(fid) trim(libPieceElms(5))//'="'// &
               trim(STR(vtk%pcElms%nStrips))//'"    '
         end if
         if (vtk%pcElms%nPolys .gt. 0) then
            write(fid) trim(libPieceElms(6))//'="'// &
               trim(STR(vtk%pcElms%nPolys))//'"    '
         end if
         l1 = .false.
         do i=1, 6
            if (vtk%pcElms%extent(i) .ne. 0) then
               l1 = .true.
               exit
            end if
         end do
         if (l1) then
            write(fid) trim(libPieceElms(7))//'="'
            do i=1, 6
               if (i .eq. 1) then
                  write(fid) trim(STR(vtk%pcElms%extent(i)))
               else
                  write(fid) ' '//trim(STR(vtk%pcElms%extent(i)))
               end if
            end do
            write(fid) '"'
         end if
         write(fid) '>'//newl

         select case(trim(vtk%dataType%str))
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
         case("ImageData")
            dAttToW = .false.
            dAttToW(1) = .true.   ! PointData
            dAttToW(2) = .true.   ! CellData
            call vtkWriteDataArr(vtk, dAttToW, istat)
            if (istat .lt. 0) return
         end select
         write(fid) '    </Piece>'//newl
         write(fid) '  </'//trim(vtk%dataType%str)//'>'//newl

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

         vtk%pcElms%nPoints = nn     ! NumberOfPoints
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
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0.0_RK
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

         select case(trim(vtk%dataType%str))
         case("UnstructuredGrid")
            iatt = 9
            allocate(ioff(vtk%pcAtt(iatt)%n))
            ioff = 0
            if ( vtk%pcElms%nCells.eq.0 ) then
               vtk%pcElms%nCells = ne
               vtk%pcAtt(iatt)%dataArr(1)%nVals  = ne*eNoN
               vtk%pcAtt(iatt)%dataArr(2)%nVals  = ne
               vtk%pcAtt(iatt)%dataArr(3)%nVals  = ne
               do i=1, vtk%pcAtt(iatt)%n
                  n = vtk%pcAtt(iatt)%dataArr(i)%nVals
                  vtk%pcAtt(iatt)%dataArr(i)%nElms = n
                  allocate(vtk%pcAtt(iatt)%dataArr(i)%iarr(n))
               end do
            else
               vtk%pcElms%nCells = vtk%pcElms%nCells + ne
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
            if ( vtk%pcElms%nPolys.eq.0 ) then
               vtk%pcElms%nPolys = ne
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
               vtk%pcElms%nPolys = vtk%pcElms%nPolys + ne
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

         subroutine putVTK_numPoints(vtk,nn,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(in) :: nn
         integer(IK), intent(inout) :: istat

         istat = 0
         vtk%pcElms%nPoints = nn

         end subroutine putVTK_numPoints

         !==========================================

         subroutine putVTK_numElems(vtk,ne,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(in) :: ne
         integer(IK), intent(inout) :: istat

         istat = 0
         select case(trim(vtk%dataType%str))
         case("UnstructuredGrid")
            vtk%pcElms%nCells = ne

         case("PolyData")
            vtk%pcElms%nPolys = ne
         end select

         end subroutine putVTK_numElems

         !==========================================

         subroutine putVTK_imageExtent(vtk,iLims,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(in) :: iLims(2,3)
         integer(IK), intent(inout) :: istat

         istat = 0
         vtk%dataType%extent(:) = reshape(iLims,(/6/))

         return
         end subroutine putVTK_imageExtent

         !==========================================

         subroutine putVTK_imageOrigin(vtk,orig,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         real(RK), intent(in) :: orig(maxNSD)
         integer(IK), intent(inout) :: istat

         istat = 0
         vtk%dataType%origin = orig

         return
         end subroutine putVTK_imageOrigin

         !==========================================

         subroutine putVTK_imageSpacing(vtk,dx,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         real(RK), intent(in) :: dx(maxNSD)
         integer(IK), intent(inout) :: istat

         istat = 0
         vtk%dataType%spacng = dx

         return
         end subroutine putVTK_imageSpacing

         !==========================================

         subroutine putVTK_pieceExtent(vtk,pLims,istat)
         implicit none
         type(vtkXMLType), intent(inout) :: vtk
         integer(IK), intent(in) :: pLims(2,3)
         integer(IK), intent(inout) :: istat

         istat = 0
         vtk%pcElms%extent = reshape(pLims,(/6/))

         return
         end subroutine putVTK_pieceExtent

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
         iatt  = 1
         n     = size(u)

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
         iatt  = 1
         m     = size(u,1)
         n     = size(u,2)

         if (m .eq. 1) then
            call putVTK_pointDataIntS(vtk,dName,u(1,:),istat)
            return
         end if

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
               ind = (j-1)*vtk%pcAtt(iatt)%dataArr(i)%nComps + k
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
         iatt  = 1
         n     = size(u)

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
         iatt  = 1
         m     = size(u,1)
         n     = size(u,2)

         if (m .eq. 1) then
            call putVTK_pointDataRealS(vtk,dName,u(1,:),istat)
            return
         end if

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
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0.0_RK
         do j=1, n
            do k=1, m
               ind = (j-1)*vtk%pcAtt(iatt)%dataArr(i)%nComps + k
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
         iatt  = 2
         n     = size(u)

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
         iatt  = 2
         m     = size(u,1)
         n     = size(u,2)

         if (m .eq. 1) then
            call putVTK_elemDataIntS(vtk,dName,u(1,:),istat)
            return
         end if

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
               ind = (j-1)*vtk%pcAtt(iatt)%dataArr(i)%nComps + k
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
         iatt  = 2
         n     = size(u)

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
         iatt  = 2
         m     = size(u,1)
         n     = size(u,2)
         if (m .eq. 1) then
            call putVTK_elemDataRealS(vtk,dName,u(1,:),istat)
            return
         end if

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
         vtk%pcAtt(iatt)%dataArr(i)%darr = 0.0_RK
         do j=1, n
            do k=1, m
               ind = (j-1)*vtk%pcAtt(iatt)%dataArr(i)%nComps + k
               vtk%pcAtt(iatt)%dataArr(i)%darr(ind) = u(k,j)
            end do
         end do

         return
         end subroutine putVTK_elemDataRealV

         !==========================================

      end module vtkXMLMod

!**************************************************


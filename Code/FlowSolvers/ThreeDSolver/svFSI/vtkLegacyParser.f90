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

      module vtkLegacyMod
      use dataTypeParams
      use stdParams
      use genUtils

      integer(IK), parameter :: maxPtDOF=30
      integer(IK), parameter :: maxElDOF=12
      integer(IK), parameter :: scalarDOF=1
      integer(IK), parameter :: vectorDOF=maxNSD

      type vtkPtData
         integer(IK) :: nVar,tDOF
         integer(IK), dimension(maxPtDOF) :: ndof,ioff,ikind
         character(len=strL), dimension(maxPtDOF) :: varName
         character(len=strL), dimension(maxPtDOF) :: varType
         real(RK4), dimension(:), allocatable :: arr
      end type vtkPtData

      type vtkElData
         integer(IK) :: nVar,tDOF
         integer(IK), dimension(maxElDOF) :: ndof,ioff,ikind
         character(len=strL), dimension(maxElDOF) :: varName
         character(len=strL), dimension(maxElDOF) :: varType
         real(RK4), dimension(:), allocatable :: arr
      end type vtkElData

      type vtkUnstrucGridType
         logical :: isBinary
         character(len=strL) :: dAtt
         integer(IK) :: nsd,cellType
         integer(IK) :: nNo, nEl, eNoN
         integer(IK) :: scalarDOF,vectorDOF,nFields
         integer(IK), dimension(:,:), allocatable :: ien
         real(RK4), dimension(:,:), allocatable :: x
         type(vtkPtData) :: ptData
         type(vtkElData) :: elData
      end type vtkUnstrucGridType

      logical :: flag
      integer :: nn, ne, ivar
      save nn, ne, ivar

      contains

         !==========================================

         subroutine loadLegacyVTK(vtk,fName,istat)
         implicit none
         type(vtkUnstrucGridType), intent(inout) :: vtk
         character(len=*), intent(in) :: fName
         integer(IK) :: istat

         logical :: flag
         character(len=strL) :: rLine,stmp
         character(len=strL), dimension(maxToks) :: tokenList
         character :: c

         integer(IK) :: i,j,fid,iPos,n1,n2
         integer(IK) :: itmp,ntoks,slen
         real(RK4)   :: r4tmp

         istat = 0
         inquire(file=trim(fName), exist=flag)
         if ( .not.flag ) then
            write(stdout,ftab4) &
               "ERROR: File "//trim(fName)//" does not exist"
            istat=-1; return
         end if

         write(stdout,ftab1) &
            "Reading file: <"//trim(fName)//">"

         vtk%isBinary = .false.
         rLine = ''
         vtk%dAtt = ''
         stmp  = ''

         fid  = 10
         open(fid,file=trim(fName),status='old')
         read(fid,*) ! # vtk DataFile Version 3.0
         read(fid,*) ! Simulation Results < dummy hdr >
         read(fid,'(A)') rLine ! ASCII or BINARY

         rLine = adjustl(rLine)
         if ( rLine(1:6) .eq. "BINARY" ) then
            vtk%isBinary = .true.
            close(fid)
            open(fid,file=trim(fName),status="unknown", &
            access="stream",form="unformatted",convert="big_endian")
            iPos = 1
         end if
         rewind(fid)

         if (vtk%isBinary) then
            write(stdout,ftab2) "Data format: <BINARY>"
         else
            write(stdout,ftab2) "Data format: <ASCII>"
         end if

         call findKwrd(vtk,fid,"POINTS",iPos,rLine,istat)
         if (istat .eq. -1) return
         call parseString(rLine,tokenList,ntoks)
         stmp = tokenList(2)
         slen = len(trim(stmp))
         read(stmp(1:slen),*) vtk%nNo

         allocate(vtk%x(maxNSD,vtk%nNo)); vtk%x = 0.0
         if (vtk%isBinary) then
            do i=1, vtk%nNo
               read(fid,pos=iPos) vtk%x(:,i)
               iPos = iPos + maxNSD*kind(r4tmp)
            end do
         else
            do i=1, vtk%nNo
               read(fid,*) vtk%x(:,i)
            end do
         end if

         call findKwrd(vtk,fid,"CELL_TYPES",iPos,rLine,istat)
         if (istat .eq. -1) return
         if (vtk%isBinary) then
            read(fid,pos=iPos) vtk%cellType
         else
            read(fid,*) vtk%cellType
         end if
         select case(vtk%cellType)
         case(5)
            vtk%eNoN = 3
         case(10)
            vtk%eNoN = 4
         case default
            write(stdout,ftab3) &
               "ERROR: inconsistent cell type. "// &
               "Only tri or tet elements allowed"
            istat=-1; return
         end select
         iPos = 1
         rewind(fid)

         call findKwrd(vtk,fid,"CELLS",iPos,rLine,istat)
         if (istat .eq. -1) return
         rLine = adjustl(rLine(6:))
         read(rLine,*) vtk%nEl, i
         if (i .ne. vtk%nEl*(vtk%eNoN+1)) then
            write(stdout,ftab4) &
               "ERROR: inconsistent element type and connectivity.."
            istat=-1; return
         end if

         allocate(vtk%ien(vtk%eNoN,vtk%nEl)); vtk%ien = 0

         if (vtk%isBinary) then
            do i=1, vtk%nEl
               read(fid,pos=iPos) itmp, vtk%ien(:,i)
               iPos = iPos + (vtk%eNoN+1)*kind(itmp)
            end do
         else
            do i=1, vtk%nEl
               read(fid,*) itmp, vtk%ien(:,i)
            end do
         end if

         nn = vtk%nNo
         vtk%ptData%nvar = 0
         vtk%ptData%tDOF = 0
         vtk%ptData%varName(:) = ''
         vtk%ptData%varType(:) = ''

         ne = vtk%nEl
         vtk%elData%nvar = 0
         vtk%elData%tDOF = 0
         vtk%elData%varName(:) = ''
         vtk%elData%varType(:) = ''

         do
            rLine = ''
            if (vtk%isBinary) then
               do i=1, strL
                  read(fid,pos=iPos,end=001) c
                  iPos = iPos + 1
                  if (c .eq. newl) exit
                  rLine(i:i) = c
               end do
            else
               read(fid,'(A)',end=001) rLine
            end if

            rLine = adjustl(rLine)

            if (rLine(1:10) .eq. "POINT_DATA") then
               if (.not.allocated(vtk%ptData%arr)) then
                  allocate(vtk%ptData%arr(nn*maxPtDOF))
                  vtk%ptData%arr = 0.0
               end if
               vtk%dAtt = "POINT_DATA"
               cycle
            else if (rLine(1:9) .eq. "CELL_DATA") then
               if (.not.allocated(vtk%elData%arr)) then
                  allocate(vtk%elData%arr(ne*maxElDOF))
                  vtk%elData%arr = 0.0
               end if
               vtk%dAtt = "CELL_DATA"
               cycle
            else if (rLine(1:5).eq."FIELD") then
               vtk%dAtt = "FIELD"
            end if

            if (trim(vtk%dAtt) .eq. "POINT_DATA") then
               if (rLine(1:7) .eq. "SCALARS") then
                  vtk%ptData%nvar = vtk%ptData%nvar + 1
                  ivar = vtk%ptData%nvar
                  call parseString(rLine,tokenList,ntoks)
                  vtk%ptData%varName(ivar) = trim(tokenList(2))
                  vtk%ptData%varType(ivar) = trim(tokenList(3))

                  select case (trim(vtk%ptData%varType(ivar)))
                  case ("int")
                     vtk%ptData%ikind(ivar) = kind(itmp)
                  case ("float")
                     vtk%ptData%ikind(ivar) = kind(r4tmp)
                  end select

                  write(stdout,ftab3) "Data Attribute: "//trim(vtk%dAtt)
                  write(stdout,ftab4) "DataArray name: "// &
                  trim(vtk%ptData%varName(ivar))// &
                  " scalars ("// &
                  trim(vtk%ptData%varType(ivar))//")"

                  call vtkExtractData(vtk,scalarDOF,rLine,fid,iPos,istat)
                  if (istat .eq. -1) return
               else if ( rLine(1:7).eq."VECTORS" ) then
                  vtk%ptData%nvar = vtk%ptData%nvar + 1
                  ivar = vtk%ptData%nvar
                  call parseString(rLine,tokenList,ntoks)
                  vtk%ptData%varName(ivar) = trim(tokenList(2))
                  vtk%ptData%varType(ivar) = trim(tokenList(3))

                  select case (trim(vtk%ptData%varType(ivar)))
                  case ("int")
                     vtk%ptData%ikind(ivar) = kind(itmp)
                  case ("float")
                     vtk%ptData%ikind(ivar) = kind(r4tmp)
                  end select

                  write(stdout,ftab3) "Data Attribute: "//trim(vtk%dAtt)
                  write(stdout,ftab4) "DataArray name: "// &
                  trim(vtk%ptData%varName(ivar))// &
                  " vectors ("// &
                  trim(vtk%ptData%varType(ivar))//")"

                  call vtkExtractData(vtk,vectorDOF,rLine,fid,iPos,istat)
                  if (istat .eq. -1) return
               end if

            else if (trim(vtk%dAtt).eq."CELL_DATA") then
               if (rLine(1:7) .eq. "SCALARS") then
                  vtk%elData%nvar = vtk%elData%nvar + 1
                  ivar = vtk%elData%nvar
                  call parseString(rLine,tokenList,ntoks)
                  vtk%elData%varName(ivar) = trim(tokenList(2))
                  vtk%elData%varType(ivar) = trim(tokenList(3))

                  select case (trim(vtk%elData%varType(ivar)))
                  case ("int")
                     vtk%elData%ikind(ivar) = kind(itmp)
                  case ("float")
                     vtk%elData%ikind(ivar) = kind(r4tmp)
                  end select

                  write(stdout,ftab3) "Data Attribute: "//trim(vtk%dAtt)
                  write(stdout,ftab4) "DataArray name: "// &
                  trim(vtk%elData%varName(ivar))// &
                  " scalars ("// &
                  trim(vtk%elData%varType(ivar))//")"

                  call vtkExtractData(vtk,scalarDOF,rLine,fid,iPos,istat)
                  if (istat .eq. -1) return
               else if (rLine(1:7) .eq. "VECTORS") then
                  vtk%elData%nvar = vtk%elData%nvar + 1
                  ivar = vtk%elData%nvar
                  call parseString(rLine,tokenList,ntoks)
                  vtk%elData%varName(ivar) = trim(tokenList(2))
                  vtk%elData%varType(ivar) = trim(tokenList(3))

                  select case (trim(vtk%elData%varType(ivar)))
                  case ("int")
                     vtk%elData%ikind(ivar) = kind(itmp)
                  case ("float")
                     vtk%elData%ikind(ivar) = kind(r4tmp)
                  end select

                  write(stdout,ftab3) "Data Attribute: "//trim(vtk%dAtt)
                  write(stdout,ftab4) "DataArray name: "// &
                  trim(vtk%elData%varName(ivar))// &
                  " vectors ("// &
                  trim(vtk%elData%varType(ivar))//")"

                  call vtkExtractData(vtk,vectorDOF,rLine,fid,iPos,istat)
                  if (istat .eq. -1) return
               end if

            else if (trim(vtk%dAtt) .eq. "FIELD") then
               call parseString(rLine,tokenList,ntoks)
               stmp = tokenList(3)
               slen = len(trim(stmp))
               read(stmp(1:slen),*) vtk%nFields

               do j=1, vtk%nFields
                  rLine = ''
                  !if ( j.gt.1 .and. vtk%isBinary ) &
                  !      iPos = iPos+1 ! skip a new line char !

                  if (vtk%isBinary) then
                     do i=1, strL
                        read(fid,pos=iPos,end=001) c
                        iPos = iPos + 1
                        if ( c.eq.newl ) exit
                        rLine(i:i) = c
                     end do
                  else
                     read(fid,'(A)',end=001) rLine
                  end if
                  rLine = adjustl(rLine)

                  call parseString(rLine,tokenList,ntoks)
                  stmp = tokenList(2)
                  slen = len(trim(stmp))
                  read(stmp(1:slen),*) n1
                  stmp = tokenList(3)
                  slen = len(trim(stmp))
                  read(stmp(1:slen),*) n2

                  if ((n1.ne.scalarDOF .and. n1.ne.vectorDOF) .or. &
                     (n2.ne.nn .and. n2.ne.ne) ) then
                     write(stdout,ftab3) &
                     "ERROR (VTK): Unknown FIELD dimensions.."
                     istat=-1; return
                  end if

                  if (n2 .eq. nn) then
                     if (.not.allocated(vtk%ptData%arr)) then
                        allocate(vtk%ptData%arr(nn*maxPtDOF))
                        vtk%ptData%arr = 0.0
                     end if
                     vtk%dAtt = "POINT_DATA"
                     vtk%ptData%nvar = vtk%ptData%nvar + 1
                     ivar = vtk%ptData%nvar
                     vtk%ptData%varName(ivar) = trim(tokenList(1))
                     vtk%ptData%varType(ivar) = trim(tokenList(4))

                     select case (trim(vtk%ptData%varType(ivar)))
                     case ("int")
                        vtk%ptData%ikind(ivar) = kind(itmp)
                     case ("float")
                        vtk%ptData%ikind(ivar) = kind(r4tmp)
                     end select

                     if ( n1.eq.scalarDOF ) then
                        write(stdout,ftab3) "Data Attribute: (FIELD) "//&
                        trim(vtk%dAtt)
                        write(stdout,ftab4) "DataArray name: "// &
                        trim(vtk%ptData%varName(ivar))// &
                        " scalars ("// &
                        trim(vtk%ptData%varType(ivar))//")"
                     else if ( n1.eq.vectorDOF ) then
                        write(stdout,ftab3) "Data Attribute: (FIELD) "//&
                        trim(vtk%dAtt)
                        write(stdout,ftab4) "DataArray name: "// &
                        trim(vtk%ptData%varName(ivar))// &
                        " vectors ("// &
                        trim(vtk%ptData%varType(ivar))//")"
                     else
                        write(stdout,ftab3) &
                           "ERROR (VTK): Unknown FIELD dimensions.."
                        istat=-1; return
                     end if

                     call vtkExtractData(vtk,n1,rLine,fid,iPos,istat)
                     if (istat .eq. -1) return
                  else if (n2 .eq. ne) then
                     if (.not.allocated(vtk%elData%arr)) then
                        allocate(vtk%elData%arr(ne*maxElDOF))
                        vtk%elData%arr = 0.0
                     end if
                     vtk%dAtt = "CELL_DATA"
                     vtk%elData%nvar = vtk%elData%nvar + 1
                     ivar = vtk%elData%nvar
                     vtk%elData%varName(ivar) = trim(tokenList(1))
                     vtk%elData%varType(ivar) = trim(tokenList(4))

                     select case (trim(vtk%elData%varType(ivar)))
                     case ("int")
                        vtk%elData%ikind(ivar) = kind(itmp)
                     case ("float")
                        vtk%elData%ikind(ivar) = kind(r4tmp)
                     end select

                     if ( n1.eq.scalarDOF ) then
                        write(stdout,ftab3) "Data Attribute: (FIELD) "//&
                        trim(vtk%dAtt)
                        write(stdout,ftab4) "DataArray name: "// &
                        trim(vtk%elData%varName(ivar))// &
                        " scalars ("// &
                        trim(vtk%elData%varType(ivar))//")"
                     else if ( n1.eq.vectorDOF ) then
                        write(stdout,ftab3) "Data Attribute: (FIELD) "//&
                        trim(vtk%dAtt)
                        write(stdout,ftab4) "DataArray name: "// &
                        trim(vtk%elData%varName(ivar))// &
                        " vectors ("// &
                        trim(vtk%elData%varType(ivar))//")"
                     else
                        write(stdout,ftab3) &
                           "ERROR (VTK): Unknown FIELD dimensions.."
                        istat=-1; return
                     end if

                     call vtkExtractData(vtk,n1,rLine,fid,iPos,istat)
                     if (istat .eq. -1) return
                  end if ! n2

               end do
            end if
         end do

 001     close(fid)

         end subroutine loadLegacyVTK

         !==========================================

         subroutine findKwrd(vtk,fileId,sKwrd,iPos,sLine,istat)
         use genUtils
         implicit none
         type(vtkUnstrucGridType), intent(in) :: vtk
         integer, intent(in) :: fileId
         integer, intent(inout) :: iPos,istat
         character(len=*), intent(in) :: sKwrd
         character(len=strL), intent(out) :: sLine

         integer :: i,kwrdL
         character :: c

         kwrdL = len(trim(sKwrd))
         do
            if (vtk%isBinary) then
               sLine = ''
               do i=1, strL
                  read(fileId,pos=iPos,end=001) c
                  iPos = iPos + 1
                  if(c .eq. newl) exit
                  sLine(i:i) = c
               end do
            else
               read(fileId,'(A)',end=001) sLine
            end if
            sLine = adjustl(sLine)
            if (sLine(1:kwrdL) .eq. trim(sKwrd)) return
         end do

 001     write(stdout,ftab4) "ERROR: EOF reached while finding "// &
         "keyword <"//trim(sKwrd)//">"
         istat = -1
         return

         end subroutine findKwrd

         !==========================================

         subroutine vtkExtractData(vtk,nComps,rLine,fid,iPos,istat)
         implicit none
         type(vtkUnstrucGridType), intent(inout) :: vtk
         integer, intent(in) :: nComps
         integer, intent(inout) :: fid,iPos,istat
         character(len=*), intent(inout) :: rLine

         integer :: i,j
         integer :: ist,iend,vecl
         character(len=strL) :: vType
         integer, dimension(:), allocatable :: tmpI

         vType = ''
         select case (trim(vtk%dAtt))
         case ("POINT_DATA")

            vType = trim(vtk%ptData%varType(ivar))
            vtk%ptData%ndof(ivar) = nComps
            vtk%ptData%tdof = vtk%ptData%tdof + &
               vtk%ptData%ndof(ivar)

            ist  = ( vtk%ptData%tdof - &
                   vtk%ptData%ndof(ivar) )*nn + 1
            vecl =  vtk%ptData%ndof(ivar)*nn
            iend =  ist + vecl - 1
            vtk%ptData%ioff(ivar) = ist-1
            if (vtk%isBinary) then
               if (rLine(1:7) .eq. "SCALARS") &
                  call findKwrd(vtk,fid,"LOOKUP_TABLE",iPos,rLine,istat)
               if (trim(vType) .eq. "int") then
                  if (allocated(tmpI)) deallocate(tmpI)
                  allocate(tmpI(vecl))
                  read(fid,pos=iPos,end=001) tmpI(1:vecl)
                  vtk%ptData%arr(ist:iend) = &
                     real(tmpI(1:vecl))
               else if (trim(vType) .eq. "float") then
                  read(fid,pos=iPos,end=001) &
                     vtk%ptData%arr(ist:iend)
               else
                  write(stdout,ftab3) &
                  "ERROR (VTK): Unknown data type. <"//trim(vType)// &
                  "> for var "//trim(vtk%ptData%varName(ivar))
                  istat=-1; return
               end if
               iPos = iPos + vecl*vtk%ptData%ikind(ivar)
               !iPos = iPos + 1 ! newl !
            else
               if (rLine(1:7) .eq. "SCALARS") &
                  read(fid,'(A)',end=001) rLine ! LOOKUP_TABLE
               do i=1, nn
                  read(fid,*,end=001) &
                     ( vtk%ptData%arr( ist-1+ &
                     (nComps*(i-1)+j) ),j=1, nComps )
               end do
            end if

         case ("CELL_DATA")

            vType = trim(vtk%elData%varType(ivar))
            vtk%elData%ndof(ivar) = nComps
            vtk%elData%tdof = vtk%elData%tdof + &
               vtk%elData%ndof(ivar)

            ist  = ( vtk%elData%tdof - &
                   vtk%elData%ndof(ivar) )*ne + 1
            vecl =  vtk%elData%ndof(ivar)*ne
            iend =  ist + vecl - 1
            vtk%elData%ioff(ivar) = ist-1
            if (vtk%isBinary) then
               if (rLine(1:7) .eq. "SCALARS") &
                  call findKwrd(vtk,fid,"LOOKUP_TABLE",iPos,rLine,istat)
               if (trim(vType) .eq. "int") then
                  if (allocated(tmpI)) deallocate(tmpI)
                  allocate(tmpI(vecl))
                  read(fid,pos=iPos,end=001) tmpI(1:vecl)
                  vtk%elData%arr(ist:iend) = &
                     real(tmpI(1:vecl))
               else if (trim(vType) .eq. "float") then
                  read(fid,pos=iPos,end=001) &
                     vtk%elData%arr(ist:iend)
               else
                  write(stdout,ftab3) &
                  "ERROR (VTK): Unknown data type. <"//trim(vType)// &
                  "> for var "//trim(vtk%elData%varName(ivar))
                  istat=-1; return
               end if
               iPos = iPos + vecl*vtk%elData%ikind(ivar)
               !iPos = iPos + 1 ! newl !
            else
               if (rLine(1:7) .eq. "SCALARS") &
                  read(fid,'(A)',end=001) rLine ! LOOKUP_TABLE
               do i=1, ne
                  read(fid,*,end=001) &
                     ( vtk%elData%arr( ist-1+ &
                     (nComps*(i-1)+j) ),j=1, nComps )
               end do
            end if

         end select

         return

 001     write(stdout,ftab4) "ERROR: EOF reached while reading data"
         istat = -1
         return

         end subroutine vtkExtractData

         !==========================================

      end module vtkLegacyMod


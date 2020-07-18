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
!     Here, all functions for treating immersed boundaries are defined.
!
!--------------------------------------------------------------------

!     This routine reads IB mesh data
      SUBROUTINE IB_READMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL :: flag
      INTEGER(KIND=IKIND) :: i, j, iM, iFa, a, b, Ac, e
      REAL(KIND=RKIND) :: fibN(nsd), rtmp
      CHARACTER(LEN=stdL) :: ctmp, fExt
      TYPE(listType), POINTER :: lPtr, lPM
      TYPE(fileType) :: fTmp

      REAL(KIND=RKIND), ALLOCATABLE :: tmpX(:,:), gX(:,:)

      ib%nMsh  = list%srch("Add IB",ll=1)
      std = " Number of immersed boundaries: "//ib%nMsh
      ALLOCATE (ib%msh(ib%nMsh), gX(0,0))

      ib%tnNo = 0
      DO iM=1, ib%nMsh
         lPM => list%get(ib%msh(iM)%name, "Add IB", iM)
         lPtr => lPM%get(ib%msh(iM)%lShl, "Set mesh as shell")
         IF (ib%msh(iM)%lShl) err = "Immersed shells are not allowed"

         std  = " Reading IB mesh <"//CLR(TRIM(ib%msh(iM)%name))//">"
         CALL READSV(lPM, ib%msh(iM))
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READCCNE(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READNRB(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            CALL READGAMBIT(lPM, ib%msh(iM))
         END IF
         IF (ib%msh(iM)%eType .EQ. eType_NA) THEN
            err = " Failed to identify format of the IB mesh"
         END IF

         std = " Number of IB nodes: "//ib%msh(iM)%gnNo
         std = " Number of IB elements: "//ib%msh(iM)%gnEl

!     Making sure face names are unique
         DO iFa=1, ib%msh(iM)%nFa
            ib%msh(iM)%fa(iFa)%iM = iM
            ctmp = ib%msh(iM)%fa(iFa)%name
            DO i=1, iM
               DO j=1, ib%msh(i)%nFa
                  IF (ctmp .EQ. ib%msh(i)%fa(j)%name .AND.
     2               (i.NE.iM .OR. j.NE.iFa)) THEN
                     err = " Repeating face names is not allowed"
                  END IF
               END DO
            END DO
         END DO

!     To scale the mesh, while attaching x to gX
         ib%msh(iM)%scF = 1._RKIND
         lPtr => lPM%get(ib%msh(iM)%scF,"Mesh scale factor",lb=0._RKIND)
         a = ib%tnNo + ib%msh(iM)%gnNo
         IF (iM .GT. 1) THEN
            ALLOCATE(tmpX(nsd,ib%tnNo))
            tmpX = gX
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
            gX(:,1:ib%tnNo) = tmpX
            DEALLOCATE(tmpX)
         ELSE
            DEALLOCATE(gX)
            ALLOCATE(gX(nsd,a))
         END IF
         gX(:,ib%tnNo+1:a) = ib%msh(iM)%x * ib%msh(iM)%scF
         ib%tnNo           = a
         DEALLOCATE(ib%msh(iM)%x)

         lPtr => lPM%get(ib%msh(iM)%dx,"Mesh global edge size",1)
      END DO
      ALLOCATE(ib%x(nsd,ib%tnNo))
      ib%x = gX
      DEALLOCATE(gX)

!     Setting msh%gN, msh%lN parameter
      b = 0
      DO iM=1, ib%nMsh
         ib%msh(iM)%nNo = ib%msh(iM)%gnNo
         ALLOCATE(ib%msh(iM)%gN(ib%msh(iM)%nNo), ib%msh(iM)%lN(ib%tnNo))
         ib%msh(iM)%gN = 0
         ib%msh(iM)%lN = 0
         DO a=1, ib%msh(iM)%nNo
            b = b + 1
            ib%msh(iM)%gN(a) = b
            ib%msh(iM)%lN(b) = a
         END DO
      END DO
      IF (b .NE. ib%tnNo) err =
     2   " Mismatch in ib%tnNo. Correction needed"

!     Remap msh%gIEN array
      DO iM=1, ib%nMsh
         ib%msh(iM)%nEl = ib%msh(iM)%gnEl
         ALLOCATE(ib%msh(iM)%IEN(ib%msh(iM)%eNoN,ib%msh(iM)%nEl))
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%gIEN(a,e)
               Ac = ib%msh(iM)%gN(Ac)
               ib%msh(iM)%IEN(a,e) = Ac
            END DO
         END DO
         DEALLOCATE(ib%msh(iM)%gIEN)
      END DO

!     Re-arranging fa structure - %gN, %lN, %IEN
      b = 0
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            ALLOCATE(ib%msh(iM)%fa(iFa)%lN(ib%tnNo))
            ib%msh(iM)%fa(iFa)%lN = 0
            DO a=1, ib%msh(iM)%fa(iFa)%nNo
               Ac = ib%msh(iM)%fa(iFa)%gN(a)
               Ac = ib%msh(iM)%gN(Ac)
               ib%msh(iM)%fa(iFa)%gN(a)  = Ac
               ib%msh(iM)%fa(iFa)%lN(Ac) = a
            END DO
            DO e=1, ib%msh(iM)%fa(iFa)%nEl
               DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                  Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                  Ac = ib%msh(iM)%gN(Ac)
                  ib%msh(iM)%fa(iFa)%IEN(a,e) = Ac
               END DO
            END DO
         END DO
      END DO

!     Setting dmnId parameter here, if there is at least one mesh that
!     has defined eId.
      DO iM=1, nMsh
         lPM => list%get(ib%msh(iM)%name,"Add IB",iM)

         lPtr => lPM%get(i,"Domain (IB)",ll=0,ul=BIT_SIZE(ib%dmnId)-1)
         IF (ASSOCIATED(lPtr)) CALL SETDMNID(ib%msh(iM), i)

         lPtr => lPM%get(fTmp,"Domain (IB) file path")
         IF (ASSOCIATED(lPtr)) THEN
            i = LEN(TRIM(fTmp%fname))
            fExt = fTmp%fname(i-2:i)
            IF (TRIM(fExt).EQ."vtp" .OR. TRIM(fExt).EQ."vtu") THEN
               CALL SETDMNIDVTK(ib%msh(iM), fTmp%fname, "DOMAIN_ID")
            ELSE
               CALL SETDMNIDFF(ib%msh(iM), fTmp%open())
            END IF
         END IF

         IF (.NOT.ALLOCATED(ib%msh(iM)%eId)) err =
     2      " Immersed bodies require domain ID parameter to be set"
      END DO
      ALLOCATE(ib%dmnId(ib%tnNo))
      ib%dmnId = 0
      DO iM=1, ib%nMsh
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               ib%dmnId(Ac) = IOR(ib%dmnId(Ac),ib%msh(iM)%eId(e))
            END DO
         END DO
      END DO

!     Read fiber orientation
      flag = .FALSE.
      DO iM=1, ib%nMsh
         lPM => list%get(ib%msh(iM)%name,"Add IB",iM)
         j = lPM%srch("Fiber direction file path")
         IF (j .EQ. 0) j = lPM%srch("Fiber direction")
         IF (j .NE. 0) THEN
            flag = .TRUE.
            EXIT
         END IF
      END DO

      IF (flag) THEN
         DO iM=1, ib%nMsh
            lPM => list%get(ib%msh(iM)%name,"Add IB",iM)

            ib%msh(iM)%nFn = lPM%srch("Fiber direction file path")
            j = ib%msh(iM)%nFn
            IF (ib%msh(iM)%nFn .NE. 0) THEN
               ALLOCATE(ib%msh(iM)%fN(j*nsd,ib%msh(iM)%nEl))
               ib%msh(iM)%fN = 0._RKIND
               DO i=1, ib%msh(iM)%nFn
                  lPtr => lPM%get(cTmp, "Fiber direction file path", i)
                  CALL READFIBNFF(ib%msh(iM), cTmp, "FIB_DIR", i)
               END DO
            ELSE
               ib%msh(iM)%nFn = lPM%srch("Fiber direction")
               j = ib%msh(iM)%nFn
               IF (ib%msh(iM)%nFn .NE. 0) THEN
                  ALLOCATE(ib%msh(iM)%fN(j*nsd,ib%msh(iM)%nEl))
                  ib%msh(iM)%fN = 0._RKIND
                  DO i=1, ib%msh(iM)%nFn
                     lPtr => lPM%get(fibN, "Fiber direction", i)
                     rtmp = SQRT(NORM(fibN))
                     IF (.NOT.ISZERO(rtmp)) fibN(:) = fibN(:)/rtmp
                     DO e=1, ib%msh(iM)%nEl
                        ib%msh(iM)%fN((i-1)*nsd+1:i*nsd,e) = fibN(1:nsd)
                     END DO
                  END DO
               END IF
            END IF
         END DO
      ELSE
         ib%msh(:)%nFn = 0
      END IF

      IF (ib%nMsh .GT. 1) THEN
         std = " Total number of IB nodes: "//ib%tnNo
         std = " Total number of IB elements: "//SUM(ib%msh%nEl)
      END IF

      std = CLR(" IB mesh data imported successfully",3)

      RETURN
      END SUBROUTINE IB_READMSH
!####################################################################
!     This routine reads IB domain properties and BCs in a given Eq
      SUBROUTINE IB_READEQ(lEq, list, eqName)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=stdL), INTENT(IN) :: eqName

      INTEGER(KIND=IKIND), PARAMETER :: maxOutput = 5

      INTEGER(KIND=IKIND) i, iBc, propL(maxNProp), iDmn, iProp, prop,
     2   nDOP(4), outputs(maxOutput)
      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPD, lPBC

      IF (.NOT.ALLOCATED(ib%msh)) err = " No IB mesh is read yet"
      IF (.NOT.ALLOCATED(ib%dmnId)) err = " No IB domain is read yet"

      lEq%nDmnIB = list%srch("Domain (IB)")
      IF (lEq%nDmnIB .EQ. 0) THEN
         err = "IB domain properties not specified"
      ELSE
         ALLOCATE(lEq%dmnIB(lEq%nDmnIB))
      END IF

!--------------------------------------------------------------------
!     Searching through each IB domain properties
      DO iDmn=1, lEq%nDmnIB
         lPD => list%get(lEq%dmnIB(iDmn)%Id,"Domain (IB)",iDmn,
     2      ll=0,ul=(BIT_SIZE(ib%dmnId)-1))
         DO i=1, iDmn-1
            IF (lEq%dmnIB(iDmn)%Id .EQ. lEq%dmnIB(i)%Id) THEN
               err = TRIM(list%ping("Domain (IB)",lPD))//
     2            " Repeated IB domain ID"
            END IF
         END DO

!        IB Equation being solved: struct
!        Other potential equations include vms_struct/lElas/shell
         lPtr => lPD%get(ctmp, "Equation", 1)
         propL = prop_NA
         SELECT CASE(TRIM(ctmp))
         CASE("struct")
            lEq%dmnIB(iDmn)%phys = phys_struct
            propL(1) = solid_density
            propL(2) = solid_viscosity
            propL(3) = elasticity_modulus
            propL(4) = poisson_ratio
            propL(5) = f_x
            propL(6) = f_y
            IF (nsd .EQ. 3) propL(7) = f_z

            nDOP = (/4,1,0,0/)
            outPuts(1) = out_displacement
            outPuts(2) = out_integ
            outPuts(3) = out_velocity
            outPuts(4) = out_pressure

         CASE DEFAULT
            err = TRIM(lPD%ping("Equation (IB)",lPtr))//
     2         "IB must be solved using struct equation only"
         END SELECT

!        Domain properties
         DO iProp=1, maxNProp
            rtmp = 0._RKIND
            prop = propL(iProp)
            SELECT CASE (prop)
            CASE (prop_NA)
               EXIT
            CASE (solid_density)
               lPtr => lPD%get(rtmp,"Density",1,ll=0._RKIND)
            CASE (solid_viscosity)
               lPtr => lPD%get(rtmp,"Viscosity",ll=0._RKIND)
            CASE (elasticity_modulus)
               lPtr => lPD%get(rtmp,"Elasticity modulus",1,lb=0._RKIND)
            CASE (poisson_ratio)
               lPtr => lPD%get(rtmp,"Poisson ratio",1,ll=0._RKIND,
     2            ul=0.5_RKIND)
            CASE (f_x)
               lPtr => lPD%get(rtmp,"Force_X")
            CASE (f_y)
               lPtr => lPD%get(rtmp,"Force_Y")
            CASE (f_z)
               lPtr => lPD%get(rtmp,"Force_Z")
            CASE DEFAULT
               err = "Undefined properties (IB)"
            END SELECT
            lEq%dmnIB(iDmn)%prop(prop) = rtmp
         END DO

         IF (lEq%dmnIB(iDmn)%phys .EQ. phys_struct) THEN
            CALL READMATMODEL(lEq%dmnIB(iDmn), lPD)
         END IF
      END DO

!     Read IB outputs
      CALL IB_READOUTPUTS(lEq, nDOP, outPuts, list)

!     Set number of function spaces
      DO i=1, nMsh
         ib%msh(i)%nFs = 1
      END DO

!--------------------------------------------------------------------
!     Searching for BCs on immersed bodies
      lEq%nBcIB = list%srch("Add BC (IB)")
      ALLOCATE(lEq%bcIB(lEq%nBcIB))
      std = " Number of imposed BC for equation <"//TRIM(eqName)//
     2   ">: "//lEq%nBcIB
      DO iBc=1, lEq%nBcIB
         lPBC => list%get(ctmp,"Add BC (IB)",iBc)
         CALL IB_FINDFACE(ctmp, lEq%bcIB(iBc)%iM, lEq%bcIB(iBc)%iFa)
         CALL IB_READBC(lEq%bcIB(iBc), lPBC)
      END DO

      RETURN
      END SUBROUTINE IB_READEQ
!####################################################################
!     This subroutine is to read from input file on how to process the
!     output quantities: for VTK files or for surface/volume integrated
!     quantities. nDOP(1) is the total number of outputs, nDOP(2) is the
!     default number of outputs for VTK files, nDOP(3) is for boundaries
!     nDOP(4) is for volume
      SUBROUTINE IB_READOUTPUTS(lEq, nDOP, outputs, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: nDOP(4), outputs(nDOP(1))
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER(KIND=IKIND) nOut, iOut, i, j
      CHARACTER(LEN=stdL) ctmp, stmp
      TYPE(listType), POINTER :: lPtr, lPO

      lEq%nOutIB = nDOP(1)
      ALLOCATE(lEq%outIB(nDOP(1)))

      DO iOut=1, nDOP(1)
         SELECT CASE (outputs(iOut))
         CASE (out_displacement)
            lEq%outIB(iOut)%grp  = outGrp_D
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Displacement"
         CASE (out_integ)
            lEq%outIB(iOut)%grp  = outGrp_I
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = 1
            IF (nsd .EQ. 2) THEN
               lEq%outIB(iOut)%name = "Area"
            ELSE
               lEq%outIB(iOut)%name = "Volume"
            END IF
         CASE (out_velocity)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = 0
            lEq%outIB(iOut)%l    = nsd
            lEq%outIB(iOut)%name = "Velocity"
         CASE (out_pressure)
            lEq%outIB(iOut)%grp  = outGrp_Y
            lEq%outIB(iOut)%o    = nsd
            lEq%outIB(iOut)%l    = 1
            lEq%outIB(iOut)%name = "Pressure"
         CASE DEFAULT
            err = "IB output undefined"
         END SELECT
      END DO

!     These are the default values, we use the first nDef/nBDef outputs
      DO j=1, 3
         lEq%outIB(1:nDOP(j+1))%wtn(j) = .TRUE.
      END DO

!     First reading the outputs for VTK files and then for boundaries
!     and last for the volume
      nOut = list%srch("Output (IB)")
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output (IB)",iOut)
         SELECT CASE(TRIM(ctmp))
         CASE("Spatial")
            j = 1
         CASE("B_INT", "Boundary_integral")
            j = 2
         CASE("V_INT", "Volume_integral")
            j = 3
         CASE("Alias")
            CYCLE
         CASE DEFAULT
            j = -1
            err = TRIM(list%ping("Output (IB)",lPO))//
     2         " Undefined keyword"
         END SELECT
         DO i=1, lEq%nOutIB
            lPtr => lPO%get(lEq%outIB(i)%wtn(j),lEq%outIB(i)%name)
         END DO
      END DO

!     Read any alias names for outputs
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output (IB)",iOut)
         SELECT CASE (TRIM(ctmp))
         CASE ("Alias")
            DO i=1, lEq%nOutIB
               lPtr => lPO%get(stmp, TRIM(lEq%output(i)%name))
               IF (ASSOCIATED(lPtr)) lEq%outIB(i)%name = TRIM(stmp)
            END DO
            EXIT
         CASE DEFAULT
            CYCLE
         END SELECT
      END DO

      RETURN
      END SUBROUTINE IB_READOUTPUTS
!####################################################################
!     Finding the face ID and mesh ID based on the face name
      SUBROUTINE IB_FINDFACE(faName, iM, iFa)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: iM, iFa
      CHARACTER(LEN=stdL) :: faName

      iFa = 0
      MY_LOOP : DO iM=1, ib%nMsh
         IF (ib%msh(iM)%lShl) THEN
!        If the BC is applied on the shell itself
            IF (ib%msh(iM)%name .EQ. faName) THEN
               iFa = 0
               EXIT
            END IF
         END IF
         DO iFa=1, ib%msh(iM)%nFa
            IF (ib%msh(iM)%fa(iFa)%name .EQ. faName) EXIT MY_LOOP
         END DO
      END DO MY_LOOP

      IF (iM .GT. ib%nMsh) err =
     2   " Unable to find face (IB) <"//TRIM(faName)//">"

      RETURN
      END SUBROUTINE IB_FINDFACE
!####################################################################
!     This routine reads BC for immersed bodies
      SUBROUTINE IB_READBC(lBc, list)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL ltmp
      INTEGER(KIND=IKIND) a, b, i, j, Ac, iM, iFa, fid, nNo
      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) fTmp

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)

      iM  = lBc%iM
      iFa = lBc%iFa
      nNo = ib%msh(iM)%fa(iFa)%nNo

!     Reading the type: Dir/Neu
      lPtr => list%get(ctmp,"Type")
      SELECT CASE (ctmp)
      CASE ("Dirichlet","Dir")
         lBc%bType = IBSET(lBc%bType,bType_Dir)
      CASE ("Neumann","Neu")
         lBc%bType = IBSET(lBc%bType,bType_Neu)
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected IB BC type"
      END SELECT

      ALLOCATE(lBc%eDrn(nsd))
      lBc%eDrn = 0
      lPtr => list%get(lBc%eDrn,"Effective direction")

      ctmp = "Steady"
      lPtr => list%get(ctmp,"Time dependence")
      SELECT CASE (ctmp)
      CASE ('Steady')
         lBc%bType = IBSET(lBc%bType,bType_std)
         lPtr => list%get(lBc%g,"Value",1)
      CASE ('Unsteady')
         lBc%bType = IBSET(lBc%bType,bType_ustd)
         ALLOCATE(lBc%gt)
         lPtr => list%get(fTmp,"Temporal values file path")
         IF (ASSOCIATED(lPtr)) THEN
            ltmp = .FALSE.
            lPtr => list%get(ltmp,"Ramp function")
            lBc%gt%lrmp = ltmp
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lBc%gt%n = j
            IF (lBc%gt%lrmp) lBc%gt%n = 1
            ALLOCATE(lBc%gt%r(j))
            ALLOCATE(lBc%gt%i(j))
            CALL FFT(fid, i, lBc%gt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
            IF (.NOT.ASSOCIATED(lPtr)) err = "Undefined inputs for "//
     2         "unsteady type BC"
            lBc%gt%lrmp = .FALSE.
            fid = fTmp%open()
            READ (fid,*) lBc%gt%ti
            READ (fid,*) lBc%gt%T
            READ (fid,*) lBc%gt%qi
            READ (fid,*) lBc%gt%qs
            READ (fid,*) j
            lBc%gt%n = j
            ALLOCATE(lBc%gt%r(j))
            ALLOCATE(lBc%gt%i(j))
            DO i=1, j
               READ (fid,*) lBc%gt%r(i), lBc%gt%i(i)
            END DO
            CLOSE(fid)
         END IF
      CASE ('Coupled')
         err = " Cannot apply Coupled BCs for immersed bodies"
      CASE ('Resistance')
         err = " Cannot apply Resistance BCs for immersed bodies"
      CASE ('RCR', 'Windkessel')
         err = " Cannot apply RCR BCs for immersed bodies"
      CASE ('Spatial')
         err = " Cannot apply Spatial Neu BCs for immersed bodies"
      CASE ('General')
         lBc%bType = IBSET(lBc%bType,bType_gen)
         lPtr => list%get(ftmp,"BCT file path")
         IF (ASSOCIATED(lPtr)) THEN
            ALLOCATE(lBc%gm)
            CALL READBCT(lBc%gm, ib%msh(iM), ib%msh(iM)%fa(iFa),
     2         ftmp%fname)
         ELSE
            lPtr =>list%get(fTmp,
     2         "Temporal and spatial values file path",1)
            IF (.NOT.ASSOCIATED(lPtr)) err = "Input file for General"//
     2         " BC (bct_ib.vtp) is not provided"
            fid = fTmp%open()
            READ (fid,*) i, j, a
            IF (a .NE. nNo) THEN
               err = "Number of nodes does not match between "//
     2            TRIM(ib%msh(iM)%fa(iFa)%name)//" and "//fTmp%fname
            END IF
            IF (i.LT.1 .OR. i.GT.nsd) err = "0 < dof <= "//nsd//
     2         " is violated in "//fTmp%fname

            ALLOCATE(lBc%gm)
            ALLOCATE(lBc%gm%t(j), lBc%gm%d(i,a,j), ptr(ib%msh(iM)%nNo))
!     I am seting all the nodes to zero just in case a node is not set
            lBc%gm%d   = 0._RKIND
            lBc%gm%dof = i
            lBc%gm%nTP = j
            ptr        = 0
!     Preparing the pointer array
            DO a=1, nNo
               Ac = ib%msh(iM)%fa(iFa)%gN(a)
               Ac = ib%msh(iM)%lN(Ac)
               IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2            "detected for BC. Mesh: "//TRIM(ib%msh(iM)%name)//
     3            ", Face: "//TRIM(ib%msh(iM)%fa(iFa)%name)//
     4            ", Node: "//STR(a)//" gN: "//
     5            STR(ib%msh(iM)%fa(iFa)%gN(a))
               ptr(Ac) = a
            END DO
            DO i=1, j
               READ (fid,*) rtmp
               lBc%gm%t(i) = rtmp
               IF (i .EQ. 1) THEN
                  IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2               " should be zero in <"//TRIM(ftmp%fname)//">"
               ELSE
                  rtmp = rtmp - lBc%gm%t(i-1)
                  IF (ISZERO(rtmp) .OR. rtmp.LT.0._RKIND) err =
     2               "Non-increasing time trend is found in <"//
     3               TRIM(ftmp%fname)//">"
               END IF
            END DO

            lBc%gm%period = lBc%gm%t(j)
            DO b=1, nNo
               READ (fid,*) Ac
               IF (Ac.GT.ib%msh(iM)%nNo .OR. Ac.LE.0) THEN
                  err = "Entry "//b//" is out of bound in "//ftmp%fname
               END IF
               a = ptr(Ac)
               IF (a .EQ. 0) THEN
                  err = "Entry "//b//" not found in "//ftmp%fname
               END IF
               DO i=1, j
                  READ (fid,*) lBc%gm%d(:,a,i)
               END DO
            END DO
            CLOSE(fid)
            DEALLOCATE(ptr)
         END IF
      CASE DEFAULT
         err = TRIM(list%ping("Time dependence",lPtr))//
     2      " Unexpected type"
      END SELECT

!     To zero-out perimeter or not. Default is .true. for Dir
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) ltmp = .TRUE.
      lPtr => list%get(ltmp,"Zero out perimeter")
      lBc%bType = IBCLR(lBc%bType,bType_zp)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_zp)

!     Impose BC on the state variable or its integral
      ltmp = .TRUE.
      lPtr => list%get(ltmp,"Impose on state variable integral")
      lBc%bType = IBCLR(lBc%bType,bType_impD)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_impD)

!     Reading the spatial profile: flat/para/ud
      ctmp = "Flat"
      lPtr => list%get(ctmp,"Profile")
      SELECT CASE (ctmp)
      CASE ('Flat')
         lBc%bType = IBSET(lBc%bType,bType_flat)
      CASE ('Parabolic')
         lBc%bType = IBSET(lBc%bType,bType_para)
      CASE ('User_defined')
         lBc%bType = IBSET(lBc%bType,bType_ud)
         lPtr => list%get(fTmp,"Spatial profile file path",1)
         fid = fTmp%open()
         ALLOCATE(lBc%gx(nNo), ptr(ib%msh(iM)%nNo))
!     I am seting all the nodes to zero just in case a node is not set
         ptr = 0
!     Preparing the pointer array
         DO a=1, nNo
            lBc%gx(a) = 0._RKIND
            Ac = ib%msh(iM)%fa(iFa)%gN(a)
            Ac = ib%msh(iM)%lN(Ac)
            IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2         "detected for BC. Mesh: "//TRIM(ib%msh(iM)%name)//
     3         ", Face: "//TRIM(ib%msh(iM)%fa(iFa)%name)//
     4         ", Node: "//STR(a)//" gN: "//
     5         STR(ib%msh(iM)%fa(iFa)%gN(a))
            ptr(Ac) = a
         END DO

         DO b=1, nNo
            READ (fid,*) Ac, rtmp
            IF (Ac.GT.ib%msh(iM)%nNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//fTmp%fname
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry <"//b//" not found in face "//fTmp%fname
            END IF
            lBc%gx(a) = rtmp
         END DO
         CLOSE(fid)
         DEALLOCATE(ptr)
      CASE DEFAULT
         err = TRIM(list%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

      lBc%weakDir = .FALSE.
      lBc%flwP = .FALSE.

      RETURN
      END SUBROUTINE IB_READBC
!####################################################################
!     Allocates memory for IB data structures
      SUBROUTINE IB_MEMALLOC()
      USE COMMOD
      IMPLICIT NONE

      ALLOCATE(ib%R(nsd+1,tnNo))
      ALLOCATE(ib%Yn(nsd+1,ib%tnNo))
      ALLOCATE(ib%Un(nsd,ib%tnNo))

      RETURN
      END SUBROUTINE IB_MEMALLOC
!####################################################################
!     This routine initializes IB solution and FSILS data structures
      SUBROUTINE IB_INIT(Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, iM, iFa, iEq, iDmn, iBc, nnz
      REAL(KIND=RKIND) lD(nsd,tnNo), s(1,ib%tnNo)
      CHARACTER(LEN=stdL) sOut

      ib%callD = 0._RKIND
      ib%callD(1) = CPUT()

      std = " ================================="
      std = " Initializing IB data structures.."

      ib%cEq = 0
      DO iEq=1, nEq
         IF (eq(iEq)%phys .EQ. phys_fluid) THEN
            ib%cEq = iEq
            EXIT
         END IF
      END DO
      IF (ib%cEq .EQ. 0) err = "Fluid equation not detected "//
     2   "to treat immersed bodies"

      lD = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     Initialize function spaces
      DO iM=1, ib%nMsh
         CALL INITFSMSH(ib%msh(iM))
         DO iFa=1, ib%msh(iM)%nFa
            CALL INITFSFACE(ib%msh(iM), ib%msh(iM)%fa(iFa))
         END DO
      END DO

!     Calculating the volume of each domain
      s = 1._RKIND
      DO iEq=1, nEq
         DO iDmn=1, eq(iEq)%nDmnIB
            eq(iEq)%dmnIB(iDmn)%v =
     2         Integ(eq(iEq)%dmnIB(iDmn)%Id, s, 1, 1)
            std = " Volume of domain "//STR(eq(iEq)%dmnIB(iDmn)%v)
            IF (ISZERO(eq(iEq)%dmnIB(iDmn)%v)) wrn = "Volume of "//
     2         "domain "//iDmn//" of equation "//iEq//" is zero"
         END DO
      END DO

!     Initialize IB face normals
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            ib%msh(iM)%fa(iFa)%iM = iM
            CALL IB_FACEINI(ib%msh(iM), ib%msh(iM)%fa(iFa))
         END DO
      END DO

!     Initialize IB face BC profile
      DO iEq=1, nEq
         DO iBc=1, eq(iEq)%nBcIB
            iFa = eq(iEq)%bcIB(iBc)%iFa
            iM  = eq(iEq)%bcIB(iBc)%iM
            CALL IB_BCINI(eq(iEq)%bcIB(iBc), ib%msh(iM)%fa(iFa))
         END DO
      END DO

!     Create sparse matrix data structures
      CALL IB_LHSA(nnz)
      std = "    Non-zeros in LHS matrix (IB): "//nnz

!     Set IB Dirichlet BCs
      CALL IB_SETBCDIR(ib%Yn, ib%Un)

!     To compute IB traces, calculate node/element adjacency for
!     background lumen mesh
      DO iM=1, nMsh
         msh(iM)%iGC = 0
         CALL GETEADJCNCY(msh(iM))
         CALL GETNADJCNCY(msh(iM))
      END DO

!     Calculate node/element adjacency for immersed bodies
      DO iM=1, ib%nMsh
         CALL GETEADJCNCY(ib%msh(iM))
      END DO

!     Set iblank field for immersed boundaries
      CALL IB_SETIBLANK(lD)

!     Compute the traces of IB nodes on the background mesh
      DO iM=1, ib%nMsh
         CALL IB_INITTRACES(ib%msh(iM), lD)
      END DO

!     Initialize IB communication structure
      CALL IB_SETCOMMU()

!     Identify ghost cells for immersed boundaries
      CALL IB_SETIGHOST()

      ib%callD(1) = CPUT() - ib%callD(1)
      WRITE(sOut,'(F6.2)') ib%callD(1)
      WRITE(sOut,'(A)') "    Immersed boundary setting time: "//
     2   TRIM(sOut)//" sec"
      std = TRIM(sOut)

      std = " IB initialization complete.."
      std = " ================================="

      RETURN
      END SUBROUTINE IB_INIT
!####################################################################
!     Initializing immersed boundary faces
      SUBROUTINE IB_FACEINI(lM, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      TYPE(faceType), INTENT(INOUT) :: lFa

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, e, g, Ac, Bc, Ec
      REAL(KIND=RKIND) area, tmp, xi0(nsd), xi(nsd), xp(nsd), nV(nsd)
      TYPE(fsType) :: fs, fsb

      INTEGER, ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: sV(:,:), sVl(:,:), xbl(:,:),
     2   xl(:,:), N(:), Nxi(:,:)

!     Calculating face area
      ALLOCATE(N(ib%tnNo))
      N    = 1._RKIND
      area = Integ(lFa, N)
      std  = "    Area of face <"//TRIM(lFa%name)//"> is "//STR(area)
      IF (ISZERO(area)) THEN
         IF (cm%mas()) wrn = " <"//TRIM(lFa%name)//"> area is zero"
      END IF
      lFa%area = area
      DEALLOCATE(N)

!     Compute face normals at nodes
      IF (ALLOCATED(lFa%nV)) DEALLOCATE(lFa%nV)
      ALLOCATE(lFa%nV(nsd,lFa%nNo), sV(nsd,ib%tnNo))

      flag = .FALSE.
      IF (lM%eType.EQ.eType_QUD .OR. lM%eType.EQ.eType_QTR .OR.
     2    lM%eType.EQ.eType_BIQ .OR. lM%eType.EQ.eType_QTE) THEN
         flag =.TRUE.
      END IF

      IF (.NOT.flag) THEN
!        For linear elements or NURBS, we simply project element normals
!        to nodes
         DO e=1, lFa%nEl
            IF (lFa%eType .EQ. eType_NRB) CALL NRBNNXB(lM, lFa, e)
            DO g=1, lFa%nG
               CALL GNNIB(lFa, e, g, nV)
               DO a=1, lFa%eNoN
                  Ac       = lFa%IEN(a,e)
                  sV(:,Ac) = sV(:,Ac) + nV*lFa%N(a,g)*lFa%w(g)
               END DO
            END DO
         END DO

      ELSE
!        For higher order elements, use reduced order basis on mesh
!        to project element normals. Lumping method is used to project
!        to face corners. Normals at edge nodes are computed by simple
!        interpolation from reduced basis. Standard lumping using higher
!        order basis could lead to spurious errors
         CALL SETTHOODFS(fs, lM%eType)
         CALL ALLOCFS(fs, nsd)
         CALL SETTHOODFS(fsb, lFa%eType)
         CALL GETGIP(nsd, fs%eType, fs%nG, fs%w, fs%xi)
         DO g=1, fs%nG
            CALL GETGNN(nsd, fs%eType, fs%eNoN, fs%xi(:,g), fs%N(:,g),
     2         fs%Nx(:,:,g))
         END DO
         CALL GETNNBNDS(fs%eType, fs%eNoN, fs%xib, fs%Nb)

         xi0 = 0._RKIND
         DO g=1, fs%nG
            xi0 = xi0 + fs%xi(:,g)
         END DO
         xi0 = xi0 / REAL(fs%nG, KIND=RKIND)

         ALLOCATE(sVl(nsd,lFa%eNoN), xbl(nsd,lFa%eNoN), xl(nsd,fs%eNoN),
     2      N(fs%eNoN), Nxi(nsd,fs%eNoN), ptr(fs%eNoN))
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xbl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO

            Ec  = lFa%gE(e)
            ptr = 0
            DO a=1, fs%eNoN
               Ac = lM%IEN(a,Ec)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
               DO b=1, fsb%eNoN
                  Bc = lFa%IEN(b,e)
                  IF (Ac .EQ. Bc) THEN
                     ptr(a) = b
                     EXIT
                  END IF
               END DO
            END DO

            sVl(:,:) = 0._RKIND
            DO g=1, lFa%nG
               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + lFa%N(a,g)*xbl(:,a)
               END DO

               xi = xi0
               CALL GETNNX(fs%eType, fs%eNoN, xl, fs%xib, fs%Nb, xp,
     2            xi, N, Nxi)

               CALL GNNIB(lFa, e, g, nV)

               DO a=1, fs%eNoN
                  b = ptr(a)
                  IF (b .EQ. 0) CYCLE
                  Ac = lM%IEN(a,Ec)
                  sVl(:,b) = sVl(:,b) + lFa%w(g)*N(a)*nV(:)
                  sV(:,Ac) = sV(:,Ac) + sVl(:,b)
               END DO
            END DO

            DO b=fsb%eNoN+1, lFa%eNoN
               xp = xbl(:,b)
               xi = xi0
               CALL GETNNX(fs%eType, fs%eNoN, xl, fs%xib, fs%Nb, xp,
     2            xi, N, Nxi)

               DO a=1, fs%eNoN
                  IF (ptr(a) .EQ. 0) CYCLE
                  sVl(:,b) = sVl(:,b) + N(a)*sVl(:,ptr(a))
               END DO

               Ac = lFa%IEN(b,e)
               sV(:,Ac) = sV(:,Ac) + sVl(:,b)
            END DO
         END DO
         DEALLOCATE(sVl, xbl, xl, N, Nxi, ptr)
         CALL DESTROY(fs)
         CALL DESTROY(fsb)
      END IF

      flag = .TRUE.
      DO a=1, lFa%nNo
         Ac  = lFa%gN(a)
         tmp = SQRT(NORM(sV(:,Ac)))
         IF (ISZERO(tmp)) THEN
            IF (flag) THEN
               wrn = " Skipping normal calculation of node "//a//
     2            " in face <"//TRIM(lFa%name)//">"
               flag = .FALSE.
            END IF
            lFa%nV(:,a) = 0._RKIND
            lFa%nV(1,a) = 1._RKIND
            CYCLE
         END IF
         lFa%nV(:,a) = sV(:,Ac)/tmp
      END DO
      DEALLOCATE(sV)

      RETURN
      END SUBROUTINE IB_FACEINI
!--------------------------------------------------------------------
!     Set BC spatial profile on the face
      SUBROUTINE IB_BCINI(lBc, lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) iM, iFa, jFa, i, a, b, Ac, j
      REAL(KIND=RKIND) tmp, nV(nsd), center(nsd), maxN

      INTEGER(KIND=IKIND), ALLOCATABLE :: gNodes(:)
      REAL(KIND=RKIND), ALLOCATABLE :: s(:), sV(:,:), sVl(:,:)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (BTEST(lBc%bType,bType_Neu) .AND. lBc%gm%dof.NE.1) THEN
            err = " Only DOF=1 is accepted for Neu general BCs"
         END IF
         RETURN
      END IF

      iM  = lFa%iM
      iFa = lBc%iFa
      IF (.NOT.ALLOCATED(lBc%gx)) ALLOCATE(lBc%gx(lFa%nNo))

      ALLOCATE(s(ib%tnNo))
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
         center = 0._RKIND
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            center(:) = center(:) + ib%x(:,Ac)
         END DO
         center(:) = center(:)/REAL(lFa%nNo, KIND=RKIND)
         ALLOCATE(gNodes(ib%tnNo), sV(nsd,ib%tnNo))
!     gNodes is one if a node located on the boundary (beside iFa)
         gNodes = 0
         DO jFa=1, ib%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ib%msh(iM)%fa(jFa)%nNo
               Ac         = ib%msh(iM)%fa(jFa)%gN(a)
               gNodes(Ac) = 1
            END DO
         END DO
!     "j" is a counter for the number of nodes that are located on the
!     boundary of lFa and sVl contains the list of their coordinates
         j  = 0
         sV = 0._RKIND
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            IF (gNodes(Ac) .EQ. 1) THEN
               j = j + 1
               sV(:,j) = ib%x(:,Ac)
            END IF
         END DO
         IF (cm%mas() .AND. j.EQ.0) err = " No perimeter"//
     2      " found for face "//lFa%name
!     sVl will keep the normal unit vector from center to perimeter
         ALLOCATE(sVl(nsd,j))
         DO a=1, j
            sV(:,a)  = sV(:,a) - center
            sVl(:,a) = sV(:,a)/SQRT(NORM(sV(:,a)))
         END DO
!     "s" is going to keep the ew.e value
         DO a=1, lFa%nNo
            Ac = lFa%gN(a)
            nV = ib%x(:,Ac) - center
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
         DO jFa=1, ib%msh(iM)%nFa
            IF (jFa .EQ. iFa) CYCLE
            DO a=1, ib%msh(iM)%fa(jFa)%nNo
               Ac    = ib%msh(iM)%fa(jFa)%gN(a)
               s(Ac) = 0._RKIND
            END DO
         END DO
      END IF

      DO a=1, lFa%nNo
         Ac        = lFa%gN(a)
         lBc%gx(a) = s(Ac)
      END DO

      RETURN
      END SUBROUTINE IB_BCINI
!####################################################################
!     Form the LHS sparse data structures for immersed bodies
      SUBROUTINE IB_LHSA(nnz)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(OUT) :: nnz

      INTEGER(KIND=IKIND) i, j, a, b, e, mnnzeic, rowN, colN, iM

      INTEGER(KIND=IKIND), ALLOCATABLE :: uInd(:,:)

      mnnzeic = 10*MAXVAL(ib%msh(:)%eNoN)
 003  mnnzeic = mnnzeic + MAX(5,mnnzeic/5)

      IF (ALLOCATED(uInd)) DEALLOCATE(uInd)
      ALLOCATE (uInd(mnnzeic,ib%tnNo))
      uInd = 0
      DO iM=1, ib%nMsh
         DO e=1, ib%msh(iM)%nEl
            DO a=1, ib%msh(iM)%eNoN
               rowN = ib%msh(iM)%IEN(a,e)
               DO b=1, ib%msh(iM)%eNoN
                  colN = ib%msh(iM)%IEN(b,e)
                  DO i=1, mnnzeic
                     IF (uInd(i,rowN) .EQ. 0) THEN
                        uInd(i,rowN) = colN
                        EXIT
                     END IF
                     IF (colN .GT. uInd(i,rowN)) CYCLE
                     IF (colN .EQ. uInd(i,rowN)) EXIT
                     IF (uInd(mnnzeic,rowN) .NE. 0) GOTO 003
                     DO j=mnnzeic, i+1, -1
                        uInd(j,rowN) = uInd(j-1,rowN)
                     END DO
                     uInd(i,rowN) = colN
                     EXIT
                  END DO
                  IF (i .GT. mnnzeic) GOTO 003
               END DO
            END DO
         END DO
      END DO

!--------------------------------------------------------------------
!     Finding number of non-zeros in colPtr vector
      nnz = 0
      DO rowN=1, ib%tnNo
         IF (uInd(1,rowN) .EQ. 0) THEN
            err = "Node "//rowN//" is isolated"
         END IF
         DO i = 1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               nnz = nnz + 1
            END IF
         END DO
      END DO

!--------------------------------------------------------------------
!     Now constructing compact form of rowPtr and colPtr
      ALLOCATE (ib%colPtr(nnz), ib%rowPtr(ib%tnNo+1))
      j  = 1
      ib%rowPtr(1) = 1
      DO rowN=1, ib%tnNo
         DO i=1, mnnzeic
            IF (uInd(i,rowN) .NE. 0) THEN
               ib%colPtr(j) = uInd(i,rowN)
               j = j + 1
            END IF
         END DO
         ib%rowPtr(rowN+1) = j
      END DO
      DEALLOCATE (uInd)

      RETURN
      END SUBROUTINE IB_LHSA
!####################################################################
!     Set iblank field
!        iblank is set only immersed solids and not set for thin shells
!        iblank(A) = 1   =>   node A is inside the immersed solid
!        iblank(A) = 0   =>   node A is outside the immersed solid
      SUBROUTINE IB_SETIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: incL
      INTEGER(KIND=IKIND) :: a, Ac, i, iM, itmp
      REAL(KIND=RKIND) :: minb(nsd), maxb(nsd), xp(nsd), dx

      LOGICAL, ALLOCATABLE :: chck(:)

!     Initialize iblank field
      iblank(:) = 0

!     Fit a bounding box around IB and probe only those nodes lying
!     inside the box
      dx = TINY(dx)
      DO iM=1, ib%nMsh
         IF (dx .LT. ib%msh(iM)%dx) THEN
            dx = ib%msh(iM)%dx
         END IF
      END DO

      DO i=1, nsd
         minb(i) = MINVAL(ib%x(i,:) + ib%Un(i,:)) - dx
         maxb(i) = MAXVAL(ib%x(i,:) + ib%Un(i,:)) + dx
      END DO

      ALLOCATE(chck(tnNo))
      chck = .FALSE.
      DO Ac=1, tnNo
         xp(:) = x(:,Ac) + lD(:,Ac)
         itmp = 0
         DO i=1, nsd
            IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i)) THEN
               itmp = itmp + 1
            END IF
         END DO
         IF (itmp .EQ. nsd) chck(Ac) = .TRUE.
      END DO

      DO iM=1, nMsh
!        Probe each node if it is inside or outside the immersed body
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            IF (.NOT.chck(Ac)) CYCLE

            xp(:) = x(:,Ac) + lD(:,Ac)
            incL  = .FALSE.
            CALL IB_CHECKINOUT(xp, incL)
            IF (incL) iblank(Ac) = 1
            chck(Ac) = .FALSE.
         END DO
      END DO

      RETURN
      END SUBROUTINE IB_SETIBLANK
!--------------------------------------------------------------------
!     Update iblank field around the neighborhood of previous iblank. In
!     case iblank is 0 entirely on the current process, redo iblank
!     initialization to look for any fresh iblank node
      SUBROUTINE IB_UPDATEIBLANK(lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL flag
      INTEGER(KIND=IKIND) a, b, i, Ac, Bc, iswp, iM, nNo
      REAL(KIND=RKIND) xp(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), ptr(:)

!     Just in case a fresh iblank appears in the current process
      IF (SUM(iblank) .EQ. 0) CALL IB_SETIBLANK(lD)

      ALLOCATE(incNd(tnNo))
      incNd(:) = iblank(:)

!     Do multiple sweeps to get the neighborhood of iblank field. Probe
!     is performed only on these neighborhood nodes.
      iswp = 0
      DO WHILE (iswp .LE. 2)
         DO iM=1, nMsh
            DO a=1, msh(iM)%nNo
               Ac = msh(iM)%gN(a)
               IF (incNd(Ac) .EQ. 1) THEN
                  DO i=msh(iM)%nAdj%prow(a), msh(iM)%nAdj%prow(a+1)-1
                     b = msh(iM)%nAdj%pcol(i)
                     Bc = msh(iM)%gN(b)
                     incNd(Bc) = 1
                  END DO
               END IF
            END DO
         END DO
         iswp = iswp + 1
      END DO

      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) a = a + 1
      END DO
      nNo = a

!     Transfer all the nodes to be probed to a local pointer array
      ALLOCATE(ptr(nNo))
      a = 0
      DO Ac=1, tnNo
         IF (incNd(Ac) .EQ. 1) THEN
            a = a + 1
            ptr(a) = Ac
         END IF
      END DO

!     Probe each node if it is inside or outside the immersed body
      iblank = 0
      DO a=1, nNo
         Ac = ptr(a)
         xp = x(:,Ac) + lD(:,Ac)
         flag  = .FALSE.
         CALL IB_CHECKINOUT(xp, flag)
         IF (flag) iblank(Ac) = 1
      END DO
      DEALLOCATE(incNd, ptr)

      RETURN
      END SUBROUTINE IB_UPDATEIBLANK
!--------------------------------------------------------------------
!     Checks if a probe lies inside or outside an immersed boundary
      SUBROUTINE IB_CHECKINOUT(xp, flag)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd)
      LOGICAL, INTENT(INOUT) :: flag

      INTEGER(KIND=IKIND) :: a, e, Ac, Ec, iM, iFa, jM, jFa
      REAL(KIND=RKIND) :: dS, minS, Jac, nV(nsd), xb(nsd), dotP

!     Find the closest immersed face centroid from the probe
      minS = HUGE(minS)
      DO iM=1, ib%nMsh
         DO iFa=1, ib%msh(iM)%nFa
            DO e=1, ib%msh(iM)%fa(iFa)%nEl
               xb = 0._RKIND
               DO a=1, ib%msh(iM)%fa(iFa)%eNoN
                  Ac = ib%msh(iM)%fa(iFa)%IEN(a,e)
                  xb = xb + ib%x(:,Ac) + ib%Un(:,Ac)
               END DO
               xb = xb / REAL(ib%msh(iM)%fa(iFa)%eNoN, KIND=RKIND)
               dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ) )
               IF (dS .LT. minS) THEN
                  minS = dS
                  Ec   = e
                  jM   = iM
                  jFa  = iFa
               END IF
            END DO
         END DO
      END DO

!     Compute the normal of the face element
      CALL GNNIB(ib%msh(jM)%fa(jFa), Ec, 1, nV)
      Jac = SQRT(NORM(nV))
      nV  = nV(:)/Jac

!     Check the sign of op.n
      xb = 0._RKIND
      DO a=1, ib%msh(jM)%fa(jFa)%eNoN
         Ac = ib%msh(jM)%fa(jFa)%IEN(a,Ec)
         xb = xb + ib%x(:,Ac) + ib%Un(:,Ac)
      END DO
      xb   = xb / REAL(ib%msh(jM)%fa(jFa)%eNoN, KIND=RKIND)
      dotP = NORM(xp-xb, nV)

      IF (dotP .LT. -1.E-9_RKIND) THEN
!        probe lies inside IB
         flag = .TRUE.
      ELSE IF (dotP .GT. 1.E-9_RKIND) THEN
!        probe lies outside IB
         flag = .FALSE.
      ELSE
!        if the probe is along the tangent, perform sign check with one
!        of the vertices of the closest node instead of the centroid
         Ac   = ib%msh(jM)%fa(jFa)%IEN(1,Ec)
         xb   = ib%x(:,Ac) + ib%Un(:,Ac)
         dotP = NORM(xp-xb,nV)
         IF (ABS(dotP) .LT. 1.E-9_RKIND) THEN
            flag = .TRUE.
         ELSE
            flag = .FALSE.
         END IF
      END IF

      RETURN
      END SUBROUTINE IB_CHECKINOUT
!####################################################################
!     Initialize traces of IB mesh integration points on the
!     background mesh elements and tag those traces as ghost elements.
      SUBROUTINE IB_INITTRACES(lM, lD)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      INTEGER(KIND=IKIND), ALLOCATABLE :: incEl(:), ePtr(:,:,:)

!     For shells and IFEM method, all the mesh elements are included for
!     trace search. For solids, only boundary/face nodes are included
!     for trace search.
      ALLOCATE(incEl(lM%nEl), ePtr(2,lM%nG,lM%nEl))
      incEl = 1
      ePtr  = 0
      CALL IB_FINDMSHTRACES(lM, lD, incEl, ePtr)
      DEALLOCATE(incEl, ePtr)

      RETURN
      END SUBROUTINE IB_INITTRACES
!--------------------------------------------------------------------
!     Update IB tracers
      SUBROUTINE IB_UPDATE(Dg)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, iM
      REAL(KIND=RKIND), ALLOCATABLE :: lD(:,:)

      ALLOCATE(lD(nsd,tnNo))
      lD = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, tnNo
            lD(:,a) = Dg(nsd+2:2*nsd+1,a)
         END DO
      END IF

!     As the IB may have moved update iblank field
      CALL IB_UPDATEIBLANK(lD)

!     Update IB traces on background fluid mesh and update ghost cells
      DO iM=1, nMsh
         msh(iM)%iGC(:) = 0
      END DO
      DO iM=1, ib%nMsh
         CALL IB_UPDATETRACESMSH(ib%msh(iM), lD)
      END DO

!     Reset IB communicator
      CALL IB_SETCOMMU()

!     Now that we have updated iblank field and IB traces, time to
!     update ghost cells
      CALL IB_SETIGHOST()

      RETURN
      END SUBROUTINE IB_UPDATE
!--------------------------------------------------------------------
!     Update traces of IB mesh integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_UPDATETRACESMSH(lM, lD)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: ipass
      INTEGER(KIND=IKIND) :: a, e, g, i, j, Ac, Ec, iM, iswp, ne
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd)

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: incEl(:), eList(:),
     2   ePtr(:,:,:), gPtr(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:)

!     It is not necessary to find traces for all nodes. Instead, use
!     previously determined traces and update them in their neighborhood
!     elements. Likewise, neighborhood nodes of IB are also added for
!     trace search.
      ALLOCATE(ePtr(2,lM%nG,lM%nEl), incEl(lM%nEl))
      ePtr  = 0
      incEl = 0
      DO i=1, lM%trc%n
         e = lM%trc%gE(1,i)
         g = lM%trc%gE(2,i)
         incEl(e)    = 1
         ePtr(:,g,e) = lM%trc%ptr(:,i)
      END DO

!     Ignore if there are no tracers earlier. However, some elements
!     could move to the neighboring process and create fresh tracer
!     pointers. These will be tracked later.
      ipass = .TRUE.
      IF (SUM(incEl) .EQ. 0) ipass = .FALSE.

!     For elements that already have tracer pointers, search in the
!     neighborhood of previous element tracers
      ALLOCATE(gPtr(2,lM%nG,lM%nEl), ichk(lM%nEl), xl(nsd,lM%eNoN))
      gPtr = 0
      iswp = 1
      ichk = .FALSE.
      DO WHILE (iswp.LE.3 .AND. ipass)
         DO e=1, lM%nEl
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO

            DO g=1, lM%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. ichk(e)) CYCLE

               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
               END IF

               gPtr(1,g,e) = Ec
               gPtr(2,g,e) = iM
               msh(iM)%iGC(Ec) = 1

!              Initialize the neighbors using the current trace
               DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
                  j = lM%eAdj%pcol(i)
                  IF (j.EQ.e .OR. incEl(j).EQ.0 .OR. ichk(j)) CYCLE
                  ePtr(1,g,j) = Ec
                  ePtr(2,g,j) = iM
               END DO
            END DO
            ichk(e) = .TRUE.
         END DO
         iswp = iswp + 1
      END DO
      DEALLOCATE(ePtr)

!     Reset incEl for the elements whose traces are newly determined
      incEl = 1
      DO e=1, lM%nEl
         i = 0
         DO g=1, lM%nG
            IF (gPtr(1,g,e) .NE. 0) i = i + 1
         END DO
         IF (i .EQ. lM%nG) incEl(e) = 0
      END DO

!     Perform a general search for all the remaining nodes including
!     the ones that may have migrated into other process domains
      CALL IB_FINDMSHTRACES(lM, lD, incEl, gPtr)
      DEALLOCATE(incEl, gPtr)

      RETURN
      END SUBROUTINE IB_UPDATETRACESMSH
!--------------------------------------------------------------------
!     Update traces of IB face integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_UPDATETRACESFACE(lFa, lD)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)

      LOGICAL :: ipass
      INTEGER(KIND=IKIND) :: a, e, g, i, j, Ac, Ec, iM, iswp, ne
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd)

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: incEl(:), eList(:),
     2   ePtr(:,:,:), gPtr(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:)

!     It is not necessary to find traces for all nodes. Instead, use
!     previously determined traces and update them in their neighborhood
!     elements. Likewise, neighborhood nodes of IB are also added for
!     trace search.
      ALLOCATE(ePtr(2,lFa%nG,lFa%nEl), incEl(lFa%nEl))
      ePtr  = 0
      incEl = 0
      DO i=1, lFa%trc%n
         e = lFa%trc%gE(1,i)
         g = lFa%trc%gE(2,i)
         incEl(e)    = 1
         ePtr(:,g,e) = lFa%trc%ptr(:,i)
      END DO

!     Ignore if there are no tracers earlier. However, some elements
!     could move to the neighboring process and create fresh tracer
!     pointers. These will be tracked later.
      ipass = .TRUE.
      IF (SUM(incEl) .EQ. 0) ipass = .FALSE.

!     For elements that already have tracer pointers, search in the
!     neighborhood of previous element tracers
      ALLOCATE(gPtr(2,lFa%nG,lFa%nEl), ichk(lFa%nEl), xl(nsd,lFa%eNoN))
      gPtr = 0
      iswp = 0
      ichk = .FALSE.
      DO WHILE (iswp.LE.1 .AND. ipass)
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO

            DO g=1, lFa%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. ichk(e)) CYCLE

               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
               END IF

               gPtr(1,g,e) = Ec
               gPtr(2,g,e) = iM
               msh(iM)%iGC(Ec) = 1

!              Initialize the neighbors using the current trace
               DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
                  j = lFa%eAdj%pcol(i)
                  IF (j.EQ.e .OR. incEl(j).EQ.0 .OR. ichk(j)) CYCLE
                  ePtr(1,g,j) = Ec
                  ePtr(2,g,j) = iM
               END DO
            END DO
            ichk(e) = .TRUE.
         END DO
         iswp = iswp + 1
      END DO
      DEALLOCATE(ePtr)

!     Reset incEl for the elements whose traces are newly determined
      incEl = 1
      DO e=1, lFa%nEl
         i = 0
         DO g=1, lFa%nG
            IF (gPtr(1,g,e) .NE. 0) i = i + 1
         END DO
         IF (i .EQ. lFa%nG) incEl(e) = 0
      END DO

!     Perform a general search for all the remaining nodes including
!     the ones that may have migrated into other process domains
      CALL IB_FINDFACETRACES(lFa, lD, incEl, gPtr)
      DEALLOCATE(incEl, gPtr)

      RETURN
      END SUBROUTINE IB_UPDATETRACESFACE
!--------------------------------------------------------------------
!     Find traces of IB integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_FINDMSHTRACES(lM, lD, srchEl, ePtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: srchEl(lM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: ePtr(2,lM%nG,lM%nEl)

      INTEGER(KIND=IKIND) :: a, b, e, g, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd)
      TYPE(queueType) :: probeElQ

      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: incEl(:), rootEl(:), eList(:),
     2   sCount(:), disps(:), masEList(:), ptr(:), gptr(:), tmpI(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), xpL(:,:)

!     Create a list of all the mesh nodes
      ALLOCATE(xpL(nsd,lM%nNo))
      xpL = 0._RKIND
      DO a=1, lM%nNo
         Ac = lM%gN(a)
         xpL(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - lM%dx
      maxb(:) = maxb(:) + lM%dx

!     Loop over all possible background mesh and find traces of each IB
!     integration point onto the elements of the background mesh.
      DO iM=1, nMsh
!        Identify the IB elements which are owned by the current
!        process. Search for traces is performed only these elements
         ALLOCATE(incEl(lM%nEl))
         incEl(:) = srchEl(:)
         DO e=1, lM%nEl
            b = 0
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               Ac = lM%lN(Ac)
               xp = xpL(:,Ac)
               j = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               j = j + 1
               END DO
               IF (j .EQ. nsd) b = b + 1
            END DO
            IF (b .LT. lM%eNoN) incEl(e) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           Reset eList if the there is no overlap with iblank field
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b  = b + iblank(Ac)
            END DO
            IF (b .EQ. 0) eList(e) = 0
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            CALL IB_PARTMSH(lM, msh(iM), lD, eList, incEl)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incEl) .EQ. 0) THEN
            DEALLOCATE(incEl, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IB nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lM%nEl), rootEl(lM%nEl))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed element using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed element and form a queue for
!        subsequent search. The root element for the seed element is
!        assigned as a seed search element for the queued neighboring
!        elements
         ALLOCATE(xl(nsd,lM%eNoN))
         Ec = 0
         E_LOOP: DO e=1, lM%nEl
            ichk(e) = .TRUE.
            IF (incEl(e) .EQ. 0) CYCLE
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lM%nG
               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               ichk(e) = .FALSE.
               EXIT E_LOOP
            END DO
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
            CYCLE
         END IF

         ePtr(1,g,e) = Ec
         ePtr(2,g,e) = iM
         rootEl(e)   = Ec
         DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
            j = lM%eAdj%pcol(i)
            IF (incEl(j) .EQ. 1) THEN
               CALL ENQUEUE(probeElQ, j)
               rootEl(j) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeElQ, e))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(e)) CYCLE
            ichk(e) = .TRUE.

            DO g=1, lM%nG
               DO a=1, lM%eNoN
                  Ac = lM%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
               END DO

               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO

               Ec = rootEl(e)
               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!              If a trace is not found, then include the neighborhood
!              elements for search. Repeat this twice. If not found yet,
!              search using master element list. If master list search
!              also fails, continue
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) THEN
                     CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2                  Ec, xi)
                     IF (Ec .EQ. 0) CYCLE
                  END IF
               END IF

!              At this point, the trace is definitely found. Save it.
               ePtr(1,g,e) = Ec
               ePtr(2,g,e) = iM

!              Use the trace as a seed search element for the neigboring
!              IB elements in the following trace search
               rootEl(e)  = Ec
               DO i=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
                  j = lM%eAdj%pcol(i)
                  IF (incEl(j).EQ.1 .AND. .NOT.ichk(j)) THEN
                     CALL ENQUEUE(probeElQ, j)
                     rootEl(j) = Ec
                  END IF
               END DO
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeElQ)

!        Perform a brute search on any missed elements
         DO e=1, lM%nEl
            IF (.NOT.ichk(e) .AND. incEl(e).EQ.1) THEN
               DO g=1, lM%nG
                  DO a=1, lM%eNoN
                     Ac = lM%IEN(a,e)
                     xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
                  END DO

                  xp = 0._RKIND
                  DO a=1, lM%eNoN
                     xp = xp + xl(:,a)*lM%N(a,g)
                  END DO

                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END DO
            END IF
         END DO

         DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
      END DO

      DEALLOCATE(xpL)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lM%trc)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%ptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)  = e
               lM%trc%gE(2,i)  = g
               lM%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Check if all traces are found and return if yes
      i = cm%reduce(lM%trc%n)
      j = lM%nEl*lM%nG
      IF (i .GE. j) RETURN

!     If all traces are not found, add a fool-proof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lM%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      j = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(ptr(4*lM%trc%n), gptr(4*j))
      DO i=1, lM%trc%n
         ptr(4*i-3) = lM%trc%gE(1,i)
         ptr(4*i-2) = lM%trc%gE(2,i)
         ptr(4*i-1) = lM%trc%ptr(1,i)
         ptr(4*i)   = lM%trc%ptr(2,i)
      END DO

      i = 4*lM%trc%n
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lM%nG,lM%nEl), ptr(lM%nEl))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            e  = gptr(4*i-3)
            g  = gptr(4*i-2)
            Ec = gptr(4*i-1)
            iM = gptr(4*i)
            tmpI(1,g,e) = Ec
            tmpI(2,g,e) = iM
         END DO

         DO e=1, lM%nEl
            DO g=1, lM%nG
               Ec = tmpI(1,g,e)
               iM = tmpI(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(e) = 1
            END DO
         END DO

         ne = SUM(ptr)
         ALLOCATE(incEl(ne))
         i = 0
         DO e=1, lM%nEl
            IF (ptr(e) .EQ. 1) THEN
               i = i + 1
               incEl(i) = e
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF

!     Share the element list to all processes
      CALL cm%bcast(ne)
      IF (cm%slv()) ALLOCATE(incEl(ne))
      CALL cm%bcast(incEl)

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(xl(nsd,lM%eNoN))
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        IB bounding box position
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, ne
            e = incEl(i)
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lM%nG
               xp = 0._RKIND
               DO a=1, lM%eNoN
                  xp = xp + xl(:,a)*lM%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .NE. 0) THEN
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END IF
            END DO
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(xl)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lM%trc)
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lM%trc%n = i
      ALLOCATE(lM%trc%gE(2,i), lM%trc%ptr(2,i))
      i = 0
      DO e=1, lM%nEl
         DO g=1, lM%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lM%trc%gE(1,i)  = e
               lM%trc%gE(2,i)  = g
               lM%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Abort simulation if all traces are still not found
      i = cm%reduce(lM%trc%n)
      j = lM%nEl * lM%nG
      IF (i .LT. j) CALL DEBUGIBMSHTRC(lM)

      RETURN
      END SUBROUTINE IB_FINDMSHTRACES
!--------------------------------------------------------------------
!     Find traces of IB integration points on the background mesh
!     elements and those elements are tagged as ghost elements.
      SUBROUTINE IB_FINDFACETRACES(lFa, lD, srchEl, ePtr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(INOUT) :: lFa
      INTEGER(KIND=IKIND), INTENT(IN) :: srchEl(lFa%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: ePtr(2,lFa%nG,lFa%nEl)

      LOGICAL lShl
      INTEGER(KIND=IKIND) :: a, b, e, g, i, j, iM, Ac, Ec, nNe, ne, ierr
      REAL(KIND=RKIND) :: xp(nsd), xi(nsd), minb(nsd), maxb(nsd)

      TYPE(queueType) :: probeElQ
      LOGICAL, ALLOCATABLE :: ichk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: incEl(:), rootEl(:), eList(:),
     2   sCount(:), disps(:), masEList(:), ptr(:), gptr(:), tmpI(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), xpL(:,:)

!     Create a list of all the face nodes
      ALLOCATE(xpL(nsd,lFa%nNo))
      xpL = 0._RKIND
      DO a=1, lFa%nNo
         Ac = lFa%gN(a)
         xpL(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
      END DO

!     Create a bounding box around the intersection of immersed body
!     and background mesh
      minb = HUGE(minb)
      maxb = TINY(maxb)
      DO iM=1, nMsh
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            xp(:) = x(:,Ac) + lD(:,Ac)
            DO i=1, nsd
               IF (minb(i) .GT. xp(i)) minb(i) = xp(i)
               IF (maxb(i) .LT. xp(i)) maxb(i) = xp(i)
            END DO
         END DO
      END DO

      DO i=1, nsd
         minb(i) = MAX(MINVAL(xpL(i,:)), minb(i))
         maxb(i) = MIN(MAXVAL(xpL(i,:)), maxb(i))
      END DO
      minb(:) = minb(:) - ib%msh(lFa%iM)%dx
      maxb(:) = maxb(:) + ib%msh(lFa%iM)%dx

!     Determine if there are any shell surfaces
      lShl = .FALSE.
      DO i=1, ib%nMsh
         IF (ib%msh(i)%lShl) THEN
            lShl = .TRUE.
            EXIT
         END IF
      END DO

!     Loop over all possible background mesh and find traces of each IB
!     integration point onto the elements of the background mesh.
      DO iM=1, nMsh
!        Identify the IB elements which are owned by the current
!        process. Search for traces is performed only these elements
         ALLOCATE(incEl(lFa%nEl))
         incEl(:) = srchEl(:)
         DO e=1, lFa%nEl
            b = 0
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               Ac = lFa%lN(Ac)
               xp = xpL(:,Ac)
               j = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               j = j + 1
               END DO
               IF (j .EQ. nsd) b = b + 1
            END DO
            IF (b .LT. lFa%eNoN) incEl(e) = 0
         END DO

!        Create a master list of elements of the background mesh
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
!           For solids, reset eList if the element has no overlap with
!           IB based on iblank
            IF (.NOT.lShl) THEN
               b = 0
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  b  = b + iblank(Ac)
               END DO
               IF (b .EQ. 0) eList(e) = 0
            END IF
         END DO

!        Save to master list
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         ne = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               ne = ne + 1
               masElist(ne) = e
            END IF
         END DO

         IF (.NOT.cm%seq()) THEN
            CALL IB_PARTFACE(lFa, msh(iM), lD, eList, incEl)
         END IF
         DEALLOCATE(eList)
         IF (SUM(incEl) .EQ. 0) THEN
            DEALLOCATE(incEl, masEList)
            CYCLE
         END IF

!        At this stage, only the processes involved in search continue.
!        We have also identified the IB nodes local to the process and
!        prepared a master list of background mesh elements involved in
!        the trace search process.
         ALLOCATE(ichk(lFa%nEl), rootEl(lFa%nEl))
         ichk   = .FALSE.
         rootEl = 0

!        Identify a seed element using master element trace search and
!        the corresponding trace element is chosen as root element.
!        Identify the neighbors of the seed element and form a queue for
!        subsequent search. The root element for the seed element is
!        assigned as a seed search element for the queued neighboring
!        elements
         ALLOCATE(xl(nsd,lFa%eNoN))
         Ec = 0
         E_LOOP: DO e=1, lFa%nEl
            ichk(e) = .TRUE.
            IF (incEl(e) .EQ. 0) CYCLE
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lFa%nG
               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .EQ. 0) CYCLE
               ichk(e) = .FALSE.
               EXIT E_LOOP
            END DO
         END DO E_LOOP

         IF (ALL(ichk(:))) THEN
            DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
            CYCLE
         END IF

         ePtr(1,g,e) = Ec
         ePtr(2,g,e) = iM
         rootEl(e)   = Ec
         DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
            j = lFa%eAdj%pcol(i)
            IF (incEl(j) .EQ. 1) THEN
               CALL ENQUEUE(probeElQ, j)
               rootEl(j) = Ec
            END IF
         END DO

!        The nonlinear grid-grid search begins here.
         DO WHILE (DEQUEUE(probeElQ, e))
            IF (ALL(ichk(:))) EXIT
            IF (ichk(e)) CYCLE
            ichk(e) = .TRUE.

            DO g=1, lFa%nG
               DO a=1, lFa%eNoN
                  Ac = lFa%IEN(a,e)
                  xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
               END DO

               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO

               Ec = rootEl(e)
               ne = msh(iM)%eAdj%prow(Ec+1) - msh(iM)%eAdj%prow(Ec)
               IF (ALLOCATED(eList)) DEALLOCATE(eList)
               ALLOCATE(eList(ne))
               j = 0
               DO i=msh(iM)%eAdj%prow(Ec), msh(iM)%eAdj%prow(Ec+1)-1
                  j = j + 1
                  eList(j) = msh(iM)%eAdj%pcol(i)
               END DO

               CALL FINDE(xp, msh(iM), x, lD, tnNo, ne, eList, Ec, xi)

!              If a trace is not found, then include the neighborhood
!              elements for search. Repeat this twice. If not found yet,
!              search using master element list. If master list search
!              also fails, continue
               IF (Ec .EQ. 0) THEN
                  CALL IB_FPSRCH(xp, msh(iM), lD, ne, eList, 2, Ec, xi)
                  IF (Ec .EQ. 0) THEN
                     CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2                  Ec, xi)
                     IF (Ec .EQ. 0) CYCLE
                  END IF
               END IF

!              At this point, the trace is definitely found. Save it.
               ePtr(1,g,e) = Ec
               ePtr(2,g,e) = iM

!              Use the trace as a seed search element for the neigboring
!              IB elements in the following trace search
               rootEl(e)  = Ec
               DO i=lFa%eAdj%prow(e), lFa%eAdj%prow(e+1)-1
                  j = lFa%eAdj%pcol(i)
                  IF (incEl(j) .EQ. 1) THEN
                     CALL ENQUEUE(probeElQ, j)
                     rootEl(j) = Ec
                  END IF
               END DO
            END DO
         END DO
         IF (ALLOCATED(eList)) DEALLOCATE(eList)
         CALL DESTROY(probeElQ)

!        Perform a brute search on any missed elements
         DO e=1, lFa%nEl
            IF (.NOT.ichk(e) .AND. incEl(e).EQ.1) THEN
               DO g=1, lFa%nG
                  DO a=1, lFa%eNoN
                     Ac = lFa%IEN(a,e)
                     xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
                  END DO

                  xp = 0._RKIND
                  DO a=1, lFa%eNoN
                     xp = xp + xl(:,a)*lFa%N(a,g)
                  END DO

                  CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist,
     2               Ec, xi)
                  IF (Ec .EQ. 0) CYCLE
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END DO
            END IF
         END DO

         DEALLOCATE(incEl, masElist, ichk, rootEl, xl)
      END DO

      DEALLOCATE(xpL)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lFa%trc)
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lFa%trc%n = i
      ALLOCATE(lFa%trc%gE(2,i), lFa%trc%ptr(2,i))
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lFa%trc%gE(1,i)  = e
               lFa%trc%gE(2,i)  = g
               lFa%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Check if all traces are found and return if yes
      i = cm%reduce(lFa%trc%n)
      j = lFa%nEl*lFa%nG
      IF (i .GE. j) RETURN

!     If all traces are not found, add a fool-proof search on master
!     list. If the search fails here, the element is perhaps distorted
!     and hence, the simulation aborts.
!     First, create a list of all successfully found traces on master
      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lFa%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      j = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(ptr(4*lFa%trc%n), gptr(4*j))
      DO i=1, lFa%trc%n
         ptr(4*i-3) = lFa%trc%gE(1,i)
         ptr(4*i-2) = lFa%trc%gE(2,i)
         ptr(4*i-1) = lFa%trc%ptr(1,i)
         ptr(4*i)   = lFa%trc%ptr(2,i)
      END DO

      i = 4*lFa%trc%n
      CALL MPI_GATHERV(ptr, i, mpint, gptr, sCount, disps, mpint,
     2   master, cm%com(), ierr)

      DEALLOCATE(ptr, disps, sCount)

      IF (cm%mas()) THEN
         ALLOCATE(tmpI(2,lFa%nG,lFa%nEl), ptr(lFa%nEl))
         tmpI = 0
         ptr  = 0
         DO i=1, j
            e  = gptr(4*i-3)
            g  = gptr(4*i-2)
            Ec = gptr(4*i-1)
            iM = gptr(4*i)
            tmpI(1,g,e) = Ec
            tmpI(2,g,e) = iM
         END DO

         DO e=1, lFa%nEl
            DO g=1, lFa%nG
               Ec = tmpI(1,g,e)
               iM = tmpI(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) ptr(e) = 1
            END DO
         END DO

         ne = SUM(ptr)
         ALLOCATE(incEl(ne))
         i = 0
         DO e=1, lFa%nEl
            IF (ptr(e) .EQ. 1) THEN
               i = i + 1
               incEl(i) = e
            END IF
         END DO
         DEALLOCATE(tmpI, ptr)
      END IF

!     Share the element list to all processes
      CALL cm%bcast(ne)
      IF (cm%slv()) ALLOCATE(incEl(ne))
      CALL cm%bcast(incEl)

!     Loop over all the background mesh and find traces of each left out
!     integration point onto the elements of the background mesh
      ALLOCATE(xl(nsd,lFa%eNoN))
      DO iM=1, nMsh
!        Create a master list of elements of the background mesh based
!        IB bounding box position
         ALLOCATE(eList(msh(iM)%nEl))
         eList = 0
         DO e=1, msh(iM)%nEl
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               xp(:) = x(:,Ac) + lD(:,Ac)
               b = 0
               DO i=1, nsd
                  IF (xp(i).GE.minb(i) .AND. xp(i).LE.maxb(i))
     2               b = b + 1
               END DO
               IF (b .EQ. nsd) THEN
                  eList(e) = 1
                  EXIT
               END IF
            END DO
         END DO
         nNe = SUM(eList)
         ALLOCATE(masElist(nNe))
         i = 0
         DO e=1, msh(iM)%nEl
            IF (eList(e) .EQ. 1) THEN
               i = i + 1
               masElist(i) = e
            END IF
         END DO

!        Now perform search for each integration point of an element
!        whose trace was not determined earlier
         DO i=1, ne
            e = incEl(i)
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lFa%nG
               xp = 0._RKIND
               DO a=1, lFa%eNoN
                  xp = xp + xl(:,a)*lFa%N(a,g)
               END DO
               CALL FINDE(xp, msh(iM), x, lD, tnNo, nNe, masElist, Ec,
     2            xi)
               IF (Ec .NE. 0) THEN
                  ePtr(1,g,e) = Ec
                  ePtr(2,g,e) = iM
               END IF
            END DO
         END DO
         DEALLOCATE(eList, masElist)
      END DO
      DEALLOCATE(xl)

!     Transfer ePtr to trace data structure
      CALL DESTROY(lFa%trc)
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) i = i + 1
         END DO
      END DO
      lFa%trc%n = i
      ALLOCATE(lFa%trc%gE(2,i), lFa%trc%ptr(2,i))
      i = 0
      DO e=1, lFa%nEl
         DO g=1, lFa%nG
            IF (ePtr(1,g,e) .NE. 0) THEN
               i = i + 1
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               lFa%trc%gE(1,i)  = e
               lFa%trc%gE(2,i)  = g
               lFa%trc%ptr(:,i) = ePtr(:,g,e)
               msh(iM)%iGC(Ec) = 1
            END IF
         END DO
      END DO

!     Abort simulation if all traces are still not found
      i = cm%reduce(lFa%trc%n)
      j = lFa%nEl * lFa%nG
      IF (i .LT. j) CALL DEBUGIBFATRC(lFa)

      RETURN
      END SUBROUTINE IB_FINDFACETRACES
!--------------------------------------------------------------------
      SUBROUTINE IB_PARTMSH(lM, gM, lD, eList, incEl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM, gM
      INTEGER(KIND=IKIND), INTENT(IN) :: eList(gM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: incEl(lM%nEl)

      INTEGER(KIND=IKIND) :: itr, a, b, e, Ac, Bc, ierr
      REAL(KIND=RKIND) :: f, tol, dS, minS, xp(nsd), xb(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), tmpI(:), gN(:)

      ALLOCATE(gN(tnNo))
      gN = 0
      DO e=1, gM%nEl
         IF (eList(e) .EQ. 0) CYCLE
         DO a=1, gM%eNoN
            Ac = gM%IEN(a,e)
            gN(Ac) = 1
         END DO
      END DO

      itr = 0
      f   = 0.05_RKIND
      ALLOCATE(part(lM%nNo), tmpI(lM%nNo))
 001  part = 0
      tmpI = 0
      itr  = itr + 1
      f    = 2._RKIND*f
      tol  = (1._RKIND + f)*lM%dx
      DO a=1, lM%nNo
         IF (part(a) .NE. 0) CYCLE
         Ac   = lM%gN(a)
         xp   = ib%x(:,Ac) + ib%Un(:,Ac)
         minS = HUGE(minS)
         DO b=1, gM%nNo
            Bc = gM%gN(b)
            IF (gN(Bc) .EQ. 0) CYCLE
            xb = x(:,Bc) + lD(:,Bc)
            dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            part(a) = cm%tF()
         END IF
      END DO

      CALL MPI_ALLREDUCE(part, tmpI, lM%nNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      b = 0
      DO a=1, lM%nNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. lM%nNo) THEN
         wrn = "Found only "//STR(b)//" nodes in pass "//STR(itr)//
     2      " out of "//STR(lM%nNo)//" nodes"
         IF (itr .GT. 5) err = "Could not distribute all nodes in "//
     2      STR(itr)//" passes. Try changing mesh edge size."
         GOTO 001
      END IF

!     While all the IB nodes are partitioned, the ones which do not lie
!     within the current process are ignored
      DO e=1, lM%nEl
         b = 0
         DO a=1, lM%eNoN
            Ac = lM%IEN(a,e)
            Ac = lM%lN(Ac)
            IF (part(Ac) .NE. cm%tF()) b = b + 1
         END DO
         IF (b .EQ. lM%eNoN) incEl(e) = 0
      END DO

      DEALLOCATE(gN, part, tmpI)

      RETURN
      END SUBROUTINE IB_PARTMSH
!--------------------------------------------------------------------
      SUBROUTINE IB_PARTFACE(lFa, gM, lD, eList, incEl)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa
      TYPE(mshType), INTENT(IN) :: gM
      INTEGER(KIND=IKIND), INTENT(IN) :: eList(gM%nEl)
      REAL(KIND=RKIND), INTENT(IN) :: lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(INOUT) :: incEl(lFa%nEl)

      INTEGER(KIND=IKIND) :: itr, a, b, e, Ac, Bc, ierr
      REAL(KIND=RKIND) :: f, tol, dS, minS, xp(nsd), xb(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: part(:), tmpI(:), gN(:)

      ALLOCATE(gN(tnNo))
      gN = 0
      DO e=1, gM%nEl
         IF (eList(e) .EQ. 0) CYCLE
         DO a=1, gM%eNoN
            Ac = gM%IEN(a,e)
            gN(Ac) = 1
         END DO
      END DO

      itr = 0
      f   = 0.05_RKIND
      ALLOCATE(part(lFa%nNo), tmpI(lFa%nNo))
 001  part = 0
      tmpI = 0
      itr  = itr + 1
      f    = 2._RKIND*f
      tol  = (1._RKIND + f)*ib%msh(lFa%iM)%dx
      DO a=1, lFa%nNo
         IF (part(a) .NE. 0) CYCLE
         Ac   = lFa%gN(a)
         xp   = ib%x(:,Ac) + ib%Un(:,Ac)
         minS = HUGE(minS)
         DO b=1, gM%nNo
            Bc = gM%gN(b)
            IF (gN(Bc) .EQ. 0) CYCLE
            xb = x(:,Bc) + lD(:,Bc)
            dS = SQRT( SUM( (xp(:)-xb(:))**2._RKIND ))
            IF (minS .GT. dS) minS = dS
         END DO
         IF (minS .LT. tol) THEN
            part(a) = cm%tF()
         END IF
      END DO

      CALL MPI_ALLREDUCE(part, tmpI, lFa%nNo, mpint, MPI_MAX, cm%com(),
     2   ierr)

      b = 0
      DO a=1, lFa%nNo
         IF (tmpI(a) .GT. 0) b = b + 1
      END DO

      IF (b .NE. lFa%nNo) THEN
         wrn = "Found only "//STR(b)//" nodes in pass "//STR(itr)//
     2      " out of "//STR(lFa%nNo)//" nodes"
         IF (itr .GT. 5) err = "Could not distribute all nodes in "//
     2      STR(itr)//" passes. Try changing mesh edge size."
         GOTO 001
      END IF

!     While all the IB nodes are partitioned, the ones which do not lie
!     within the current process are ignored
      DO e=1, lFa%nEl
         b = 0
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            Ac = lFa%lN(Ac)
            IF (part(Ac) .NE. cm%tF()) b = b + 1
         END DO
         IF (b .EQ. lFa%eNoN) incEl(e) = 0
      END DO

      DEALLOCATE(gN, part, tmpI)

      RETURN
      END SUBROUTINE IB_PARTFACE
!--------------------------------------------------------------------
      SUBROUTINE IB_FPSRCH(xp, lM, lD, ne, eSrch, itMax, Ec, xi)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: xp(nsd), lD(nsd,tnNo)
      INTEGER(KIND=IKIND), INTENT(IN)  :: ne, eSrch(ne), itMax
      INTEGER(KIND=IKIND), INTENT(OUT) :: Ec
      REAL(KIND=RKIND), INTENT(OUT) :: xi(nsd)
      TYPE(mshType), INTENT(IN) :: lM

      LOGICAL flag
      INTEGER(KIND=IKIND) :: a, e, i, El, En, nEl, nEn, iter, is, ie

      LOGICAL, ALLOCATABLE :: eChk(:)
      INTEGER(KIND=IKIND), ALLOCATABLE :: eList(:), tmpI(:)

      nEn = ne
      ALLOCATE(eList(nEn), eChk(lM%nEl))
      eList = eSrch
      eChk  = .FALSE.
      iter  = 0

      DO WHILE (iter .LE. itMax)
         ALLOCATE(tmpI(lM%nEl))
         tmpI = 0
         nEl  = nEn
         nEn  = 0
         DO e=1, nEl
            El = eList(e)
            eChk(El) = .TRUE.
         END DO
         DO e=1, nEl
            El = eList(e)
            is = lM%eAdj%prow(El)
            ie = lM%eAdj%prow(El+1) - 1
            DO i=is, ie
               En = lM%eAdj%pcol(i)
               IF (Ec.EQ.En .OR. eChk(En)) CYCLE
               flag = .TRUE.
               DO a=1, nEn
                  IF (En .EQ. tmpI(a)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  nEn = nEn + 1
                  tmpI(nEn) = En
               END IF
            END DO
         END DO
         DEALLOCATE(eList)
         ALLOCATE(eList(nEn))
         eList(:) = tmpI(1:nEn)
         CALL FINDE(xp, lM, x, lD, tnNo, nEn, eList, Ec, xi)
         DEALLOCATE(tmpI)
         IF (Ec .NE. 0) THEN
            DEALLOCATE(eChk, eList)
            RETURN
         ELSE
            iter = iter + 1
         END IF
      END DO

      DEALLOCATE(eChk, eList)

      RETURN
      END SUBROUTINE IB_FPSRCH
!####################################################################
!     Communication structure for IB is initialized here. Here we
!     create a list of nodal traces on master that are local to other
!     processes. The master then gathers all the data, projects flow
!     variables (acceleration, velocity and pressure) and broadcasts to
!     all the processes.
      SUBROUTINE IB_SETCOMMU()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) i, a, e, n, Ac, iM, ierr, tag

      INTEGER(KIND=IKIND), ALLOCATABLE :: incNd(:), ptr(:), rA(:,:),
     2   rReq(:)

!     Free memory if already allocated
      CALL DESTROY(ib%cm)

!     Map trace pointers local to a process into a nodal vector. Note,
!     however, that the traces point to integration point of an element.
!     Therefore, we set all the nodes of an element that get contribution
!     from a valid trace.
      ALLOCATE(incNd(ib%tnNo))
      incNd = 0
      DO iM=1, ib%nMsh
         DO i=1, ib%msh(iM)%trc%n
            e = ib%msh(iM)%trc%gE(1,i)
            DO a=1, ib%msh(iM)%eNoN
               Ac = ib%msh(iM)%IEN(a,e)
               incNd(Ac) = 1
            END DO
         END DO
      END DO

!     All the included IB nodes are now mapped to a local vector
      n = SUM(incNd)
      ALLOCATE(ptr(n))
      i = 0
      DO a=1, ib%tnNo
         IF (incNd(a) .NE. 0) THEN
            i = i + 1
            ptr(i) = a
         END IF
      END DO

!     Set IB comm data structures for sequential run
      IF (cm%seq()) THEN
         ALLOCATE(ib%cm%n(1), ib%cm%gE(n))
         ib%cm%n(1)  = n
         ib%cm%gE(:) = ptr
         DEALLOCATE(incNd, ptr)
         RETURN
      END IF

!     Gather no of included nodes from each process on master. Data at
!     these nodes will later be gathered by master from other processes
      ALLOCATE(ib%cm%n(cm%np()))
      ib%cm%n = 0
      CALL MPI_GATHER(n, 1, mpint, ib%cm%n, 1, mpint, master, cm%com(),
     2   ierr)

!     The processes that do not have any nodal traces return
      IF (.NOT.cm%mas() .AND. n.EQ.0) THEN
         DEALLOCATE(incNd, ptr)
         ALLOCATE(ib%cm%gE(0))
         RETURN
      END IF

!     Master receives list of all nodal traces from other processes
      IF (cm%mas()) THEN
         n = SUM(ib%cm%n)
         a = MAXVAL(ib%cm%n)
         ALLOCATE(ib%cm%gE(n), rA(a,cm%np()), rReq(cm%np()))
         ib%cm%gE  = 0
         rA(:,:)   = 0
         DO i=1, cm%np()
            n = ib%cm%n(i)
            IF (n .EQ. 0) CYCLE
            IF (i .EQ. 1) THEN
               rA(1:n,i) = ptr(:)
            ELSE
               tag = i*100
               CALL MPI_IRECV(rA(1:n,i), n, mpint, i-1, tag, cm%com(),
     2            rReq(i), ierr)
            END IF
         END DO

         DO i=1, cm%np()
            IF (i.EQ.1 .OR. ib%cm%n(i).EQ.0) CYCLE
            CALL MPI_WAIT(rReq(i), MPI_STATUS_IGNORE, ierr)
         END DO
         DEALLOCATE(rReq)

         a = 0
         DO i=1, cm%np()
            n = ib%cm%n(i)
            ib%cm%gE(a+1:a+n) = rA(1:n,i)
            a = a + n
         END DO
      ELSE
         ALLOCATE(ib%cm%gE(0), rA(0,0))
         tag = cm%tF() * 100
         CALL MPI_SEND(ptr, n, mpint, master, tag, cm%com(), ierr)
      END IF

      DEALLOCATE(incNd, ptr, rA)

      RETURN
      END SUBROUTINE IB_SETCOMMU
!####################################################################
!     Set ighost field
!        ighost is set for both solids and shells
!        For solids, ighost(A) = 1   =>  iblank(A)=1 and is connected to
!           atleast one iblank=0 node
!        For shells, ighost(A) = 1   =>  node belongs to a ghost cell
      SUBROUTINE IB_SETIGHOST()
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER(KIND=IKIND) :: a, b, e, Ac, iM

      REAL(KIND=RKIND), ALLOCATABLE :: rG(:)

      ALLOCATE(rG(tnNo))
      rG = 0._RKIND

!     Update mesh and nodal ghost cell pointers based on iblank field
      DO iM=1, nMsh
         DO e=1, msh(iM)%nEl
            b = 0
            DO a=1, msh(iM)%eNoN
               Ac = msh(iM)%IEN(a,e)
               b = b + iblank(Ac)
            END DO
            IF (b.GT.0 .AND. b.LT.msh(iM)%eNoN) THEN
               msh(iM)%iGC(e) = 1
               DO a=1, msh(iM)%eNoN
                  Ac = msh(iM)%IEN(a,e)
                  rG(Ac) = REAL(iblank(Ac), KIND=RKIND)
               END DO
            END IF
         END DO
      END DO

c!     For Shells, use IB traces to identify ghost cells
c      DO iM=1, ib%nMsh
c         IF (ib%msh(iM)%lShl) THEN
c            DO b=1, ib%msh(iM)%trc%n
c               e  = ib%msh(iM)%trc%ptr(1,b)
c               jM = ib%msh(iM)%trc%ptr(2,b)
c               DO a=1, msh(jM)%eNoN
c                  Ac = msh(jM)%IEN(a,e)
c                  rG(Ac) = 1._RKIND
c               END DO
c            END DO
c         END IF
c      END DO
      CALL COMMU(rG)

      ighost = 0
      DO a=1, tnNo
         IF (rG(a) .GT. 1.E-6_RKIND) ighost(a) = 1
      END DO

      RETURN
      END SUBROUTINE IB_SETIGHOST
!####################################################################
!     Write IB solution to a vtu file
      SUBROUTINE IB_WRITEVTUS(lY, lU)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: lY(nsd+1,ib%tnNo), lU(nsd,ib%tnNo)

      TYPE(dataType) :: d(ib%nMsh)
      TYPE(vtkXMLType) :: vtu

      INTEGER(KIND=IKIND) :: iStat, iEq, iOut, iM, a, e, Ac, Ec, nNo,
     2   nEl, s, l, ie, is, nSh, oGrp, outDof, nOut, cOut
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: outS(:), tmpI(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: tmpV(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: outNames(:)

      IF (cm%slv()) THEN
         ib%savedOnce = .TRUE.
         RETURN
      END IF

      nOut   = 1
      outDof = nsd
      DO iEq=1, nEq
         DO iOut=1, eq(iEq)%nOutIB
            IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
            nOut   = nOut + 1
            outDof = outDof + eq(iEq)%outIB(iOut)%l
         END DO
      END DO

      ALLOCATE(outNames(nOut), outS(nOut+1))

!     Prepare all solultions in to dataType d
      nNo = 0
      nEl = 0
      DO iM=1, nMsh
         cOut           = 1
         outS(cOut)     = 1
         outS(cOut+1)   = nsd + 1
         outNames(cOut) = ""

         IF (ib%msh(iM)%eType .EQ. eType_NRB) err =
     2      " Outputs for NURBS data is under development"

         d(iM)%nNo     = ib%msh(iM)%nNo
         d(iM)%nEl     = ib%msh(iM)%nEl
         d(iM)%eNoN    = ib%msh(iM)%eNoN
         d(iM)%vtkType = ib%msh(iM)%vtkType

         ALLOCATE(d(iM)%x(outDof,d(iM)%nNo),
     2      d(iM)%IEN(d(iM)%eNoN,d(iM)%nEl))
         DO a=1, ib%msh(iM)%nNo
            Ac = ib%msh(iM)%gN(a)
            d(iM)%x(1:nsd,a) = ib%x(:,Ac)
         END DO

         DO e=1, ib%msh(iM)%nEl
            d(iM)%IEN(:,e) = ib%msh(iM)%IEN(:,e)
         END DO

         DO iEq=1, nEq
            DO iOut=1, eq(iEq)%nOutIB
               IF (.NOT.eq(iEq)%outIB(iOut)%wtn(1)) CYCLE
               l  = eq(iEq)%outIB(iOut)%l
               s  = eq(iEq)%s + eq(iEq)%outIB(iOut)%o
               e  = s + l - 1

               cOut = cOut + 1
               is   = outS(cOut)
               ie   = is + l - 1
               outS(cOut+1)   = ie + 1
               outNames(cOut) = "IB_"//TRIM(eq(iEq)%outIB(iOut)%name)

               oGrp = eq(iEq)%outIB(iOut)%grp
               SELECT CASE (oGrp)
                  CASE (outGrp_NA)
                  err = "Undefined output grp in VTK"
               CASE (outGrp_Y)
                  DO a=1, ib%msh(iM)%nNo
                     Ac = ib%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lY(s:e,Ac)
                  END DO
               CASE (outGrp_D)
                  DO a=1, ib%msh(iM)%nNo
                     Ac = ib%msh(iM)%gN(a)
                     d(iM)%x(is:ie,a) = lU(s:e,Ac)
                  END DO
               CASE DEFAULT
                  err = "Undefined output "//
     2               TRIM(eq(iEq)%outIB(iOut)%name)
               END SELECT
            END DO
         END DO

         ALLOCATE(d(iM)%xe(d(iM)%nEl,1))
         IF (.NOT.savedOnce) THEN
            IF (ALLOCATED(ib%dmnID)) THEN
               d(iM)%xe(:,1) = REAL(ib%msh(iM)%eId(:), KIND=RKIND)
            ELSE
               d(iM)%xe(:,1) = 1._RKIND
            END IF
         END  IF
         nNo = nNo +  d(iM)%nNo
         nEl = nEl +  d(iM)%nEl
      END DO

      ALLOCATE(tmpV(maxnsd,nNo))

!     Writing to vtu file (master only)
      IF (cTS .GE. 1000) THEN
         fName = STR(cTS)
      ELSE
         WRITE(fName,'(I3.3)') cTS
      END IF

      fName = TRIM(saveName)//"_ib_"//TRIM(ADJUSTL(fName))//".vtu"
      dbg = "Writing VTU"

      CALL vtkInitWriter(vtu, TRIM(fName), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (init)"

!     Writing the position data
      iOut = 1
      s    = outS(iOut)
      e    = outS(iOut+1)-1
      nSh  = 0
      tmpV = 0._RKIND
      DO iM=1, ib%nMsh
         DO a=1, d(iM)%nNo
            tmpV(1:nsd,a+nSh) = d(iM)%x(s:e,a)
         END DO
         nSh = nSh + d(iM)%nNo
      END DO
      CALL putVTK_pointCoords(vtu, tmpV(1:nsd,:), iStat)
      IF (iStat .LT. 0) err = "VTU file write error (coords)"

!     Writing the connectivity data
      nSh = -1
      DO iM=1, ib%nMsh
         ALLOCATE(tmpI(d(iM)%eNoN,d(iM)%nEl))
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
         DO iM=1, ib%nMsh
            DO a=1, d(iM)%nNo
               tmpV(:,a+nSh) = d(iM)%x(s:e,a)
            END DO
            nSh = nSh + d(iM)%nNo
         END DO
         CALL putVTK_pointData(vtu, outNames(iOut), tmpV, iStat)
         IF (iStat .LT. 0) err = "VTU file write error (point data)"
      END DO

!     Write element-based variables
      IF (.NOT.savedOnce .OR. mvMsh) THEN
         ib%savedOnce = .TRUE.
         ALLOCATE(tmpI(1,nEl))
!     Write the domain ID
         IF (ALLOCATED(ib%dmnID)) THEN
            Ec = 0
            DO iM=1, ib%nMsh
               DO e=1, d(iM)%nEl
                  Ec = Ec + 1
                  tmpI(1,Ec) = INT(d(iM)%xe(e,1), KIND=IKIND)
               END DO
            END DO
            CALL putVTK_elemData(vtu, 'Domain_ID', tmpI, iStat)
            IF (iStat .LT. 0) err = "VTU file write error (dom id)"
         END IF

!     Write the mesh ID
         IF (ib%nMsh .GT. 1) THEN
            Ec = 0
            DO iM=1, ib%nMsh
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

      DO iM=1, nMsh
         CALL DESTROY(d(iM))
      END DO

      CALL vtkWriteToFile(vtu, iStat)
      IF (iStat .LT. 0) err = "VTU file write error"

      CALL flushVTK(vtu)

      RETURN
      END SUBROUTINE IB_WRITEVTUS
!####################################################################
!     This is to check/create the txt file
      SUBROUTINE IB_CCTXT(lEq, fName, wtn)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      LOGICAL, INTENT(IN) :: wtn(2)

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      LOGICAL flag
      INTEGER(KIND=IKIND) iM, iFa, fid, iDmn, i
      CHARACTER(LEN=stdL) stmp

      fid = 1
      IF (cm%slv()) RETURN

!     i=1 are the boundary values and i=2 are volume values
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         INQUIRE(FILE=TRIM(fName(i)), EXIST=flag)
         IF (cTS.NE.0 .AND. flag) THEN
            CALL TRIMFILE(cTS+3,fName(i))
            CYCLE
         END IF

         OPEN(fid, FILE=TRIM(fName(i)))
         IF (i .EQ. 1) THEN
            DO iM=1, ib%nMsh
               DO iFa=1, ib%msh(iM)%nFa
                  stmp = ib%msh(iM)%fa(iFa)%name
                  IF (LEN(TRIM(stmp)) .LE. prL) THEN
                     WRITE(fid,'(A)', ADVANCE='NO')
     2                  ADJUSTR(stmp(1:prL))//" "
                  ELSE
                     WRITE(fid,'(A)',ADVANCE='NO') TRIM(stmp)//" "
                  END IF
               END DO
            END DO
            WRITE(fid,*)
            DO iM=1, ib%nMsh
               DO iFa=1, ib%msh(iM)%nFa
                  stmp = STR(ib%msh(iM)%fa(iFa)%area,prL)
                  WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmnIB
               stmp = "DOMAIN-"//STR(lEq%dmnIB(iDmn)%Id)
               IF (lEq%dmnIB(iDmn)%Id .EQ. -1) stmp = "ENTIRE"
               WRITE(fid,'(A)', ADVANCE='NO') ADJUSTR(stmp(1:prL))//" "
            END DO
            WRITE(fid,*)
            DO iDmn=1, lEq%nDmnIB
               stmp = STR(lEq%dmnIB(iDmn)%v,prL)
               WRITE(fid,'(A)', ADVANCE='NO') stmp(1:prL+1)
            END DO
         END IF
         WRITE(fid,*)
         WRITE(fid,*)
         CLOSE(fid)
      END DO

      RETURN
      END SUBROUTINE IB_CCTXT
!####################################################################
!     This is to write to txt file
      SUBROUTINE IB_WTXT(lEq, m, fName, tmpV, wtn, div)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: m
      CHARACTER(LEN=stdL), INTENT(IN) :: fName(2)
      REAL(KIND=RKIND), INTENT(IN) :: tmpV(maxnsd,ib%tnNo)
      LOGICAL, INTENT(IN) :: wtn(2), div

      INTEGER(KIND=IKIND), PARAMETER :: prL = 10

      INTEGER(KIND=IKIND) iM, iFa, fid, i, iDmn
      REAL(KIND=RKIND) tmp

      fid = 1
      DO i=1, 2
         IF (.NOT.wtn(i)) CYCLE

         IF (cm%mas()) OPEN(fid, FILE=TRIM(fName(i)), STATUS='OLD',
     2      POSITION='APPEND')

         IF (i .EQ. 1) THEN
            DO iM=1, ib%nMsh
               DO iFa=1, ib%msh(iM)%nFa
                  IF (m .EQ. 1) THEN
                     IF (div) THEN
                        tmp = ib%msh(iM)%fa(iFa)%area
                        tmp = Integ(ib%msh(iM)%fa(iFa),tmpV,1)/tmp
                     ELSE
                        tmp = Integ(ib%msh(iM)%fa(iFa),tmpV,1)
                     END IF
                  ELSE IF (m .EQ. nsd) THEN
                     tmp = Integ(ib%msh(iM)%fa(iFa),tmpV,1,m)
                  ELSE
                     err = "WTXT only accepts 1 and nsd"
                  END IF
                  IF (cm%mas())
     2               WRITE(fid,'(A)',ADVANCE='NO') STR(tmp,prL)//" "
               END DO
            END DO
         ELSE
            DO iDmn=1, lEq%nDmnIB
               IF (div) THEN
                  tmp = lEq%dmnIB(iDmn)%v
                  tmp = Integ(lEq%dmnIB(iDmn)%Id, tmpV, 1, m)/tmp
               ELSE
                  tmp = Integ(lEq%dmnIB(iDmn)%Id, tmpV, 1, m)
               END IF
               IF (cm%mas())
     2            WRITE(fid,'(A)', ADVANCE='NO') STR(tmp,prL)//" "
            END DO
         END IF
         IF (cm%mas()) THEN
            WRITE(fid,*)
            CLOSE(fid)
         END IF
      END DO

      RETURN
      END SUBROUTINE IB_WTXT
!####################################################################
!     Apply Dirichlet boundary conditions on the immersed faces
      SUBROUTINE IB_SETBCDIR(Yb, Ub)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(INOUT) :: Yb(nsd+1,ib%tnNo),
     2   Ub(nsd,ib%tnNo)

      LOGICAL :: eDir(maxnsd)
      INTEGER(KIND=IKIND) :: iFa, iM, iEq, iBc, a, Ac, nNo, lDof, i

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: nV(:,:), tmpA(:,:), tmpY(:,:)

      iEq = ib%cEq
      IF (iEq .EQ. 0) RETURN

      DO iBc=1, eq(iEq)%nBcIB
         IF (.NOT.BTEST(eq(iEq)%bcIB(iBc)%bType,bType_Dir)) CYCLE
         eDir = .FALSE.
         lDof = 0
         DO i=1, nsd
            IF (eq(iEq)%bcIB(iBc)%eDrn(i) .NE. 0) THEN
               eDir(i) = .TRUE.
               lDof = lDof + 1
            END IF
         END DO
         IF (lDof .EQ. 0) lDof = nsd

         iFa = eq(iEq)%bcIB(iBc)%iFa
         iM  = eq(iEq)%bcIB(iBc)%iM

!     Prepare a pointer list and normals for a face or a shell
         nNo = ib%msh(iM)%fa(iFa)%nNo
         ALLOCATE(ptr(nNo), nV(nsd,nNo))
         ptr = 0
         DO a=1, nNo
            ptr(a)  = ib%msh(iM)%fa(iFa)%gN(a)
            nV(:,a) = ib%msh(iM)%fa(iFa)%nV(:,a)
         END DO

         ALLOCATE(tmpA(lDof,nNo), tmpY(lDof,nNo))
         CALL IB_SETBCDIRL(eq(iEq)%bcIB(iBc), nNo, lDof, nV, tmpA, tmpY)

         IF (ANY(eDir)) THEN
            DO a=1, nNo
               Ac = ptr(a)
               IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                  DO i=1, nsd
                     lDof = 0
                     IF (eDir(i)) THEN
                        lDof = lDof + 1
                        Yb(i,Ac) = tmpA(lDof,a)
                        Ub(i,Ac) = tmpY(lDof,a)
                     END IF
                  END DO
               ELSE
                  DO i=1, nsd
                     lDof = 0
                     IF (eDir(i)) THEN
                        lDof = lDof + 1
                        Yb(i,Ac) = tmpY(lDof,a)
                        Ub(i,Ac) = Ub(i,Ac) + Yb(i,Ac)*dt
                     END IF
                  END DO
               END IF
            END DO
         ELSE
            DO a=1, nNo
               Ac = ptr(a)
               IF (BTEST(eq(iEq)%bcIB(iBc)%bType,bType_impD)) THEN
                  Yb(1:nsd,Ac) = tmpA(:,a)
                  Ub(1:nsd,Ac) = tmpY(:,a)
               ELSE
                  Yb(1:nsd,Ac) = tmpY(:,a)
                  Ub(1:nsd,Ac) = Ub(1:nsd,Ac) + Yb(1:nsd,Ac)*dt
               END IF
            END DO
         END IF

         DEALLOCATE(ptr, nV, tmpA, tmpY)
      END DO

      RETURN
      END SUBROUTINE IB_SETBCDIR
!--------------------------------------------------------------------
      SUBROUTINE IB_SETBCDIRL(lBc, nNo, lDof, nvL, lA, lY)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bcType), INTENT(IN) :: lBc
      INTEGER(KIND=IKIND), INTENT(IN) :: lDof, nNo
      REAL(KIND=RKIND), INTENT(IN) :: nvL(nsd,nNo)
      REAL(KIND=RKIND), INTENT(INOUT) :: lA(lDof,nNo), lY(lDof,nNo)

      INTEGER(KIND=IKIND) :: a, i
      REAL(KIND=RKIND) :: dirY, dirA, nV(nsd)

      IF (BTEST(lBc%bType,bType_gen)) THEN
         IF (lDof .NE. lBc%gm%dof) err = "Inconsistent DOF to apply "//
     2      "Gen BC"
         IF (nNo .NE. SIZE(lBc%gm%d,2)) err = "Inconsistent nodes "//
     2      "to apply Gen BC"
         CALL IGBC(lBc%gm, lY, lA)
         RETURN
      ELSE IF (BTEST(lBc%bType,bType_ustd)) THEN
         CALL IFFT(lBc%gt, dirY, dirA)
      ELSE ! std / cpl
         dirA = 0._RKIND
         dirY = lBc%g
      END IF

      IF (lDof .EQ. nsd) THEN
         DO a=1, nNo
            nV      = nvL(:,a)
            lA(:,a) = dirA*lBc%gx(a)*nV
            lY(:,a) = dirY*lBc%gx(a)*nV
         END DO
      ELSE
         DO a=1, nNo
            DO i=1, lDof
               lA(i,a) = dirA*lBc%gx(a)
               lY(i,a) = dirY*lBc%gx(a)
            END DO
         END DO
      END IF

      RETURN
      END SUBROUTINE IB_SETBCDIRL
!####################################################################
!     Compute FSI force on the immersed bodies using Immersed Finite
!     Element Method (IFEM)
      SUBROUTINE IB_CALCFFSI(Ag, Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Ag(tDof,tnNo), Yg(tDof,tnNo),
     2   Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, b, e, g, i, j, is, ie, Ac, Bc, Ec, iM, jM,
     2   eNoNb, eNoN, nFn
      REAL(KIND=RKIND) w, Jac, rt, xp(nsd), xi(nsd), Gmat(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: Nb(:), Nbx(:,:), N(:), Nxi(:,:),
     2   Nx(:,:), xbl(:,:), xl(:,:), al(:,:), yl(:,:), ul(:,:), fN(:,:),
     3   lR(:,:)

!     TODO: This is temporary. Later we get domain based on each IB node
!     and communicate across IB process boundaries
      cDmn = DOMAIN(msh(1), ib%cEq, 1)

      is = eq(ib%cEq)%s
      ie = eq(ib%cEq)%e
      IF (ie .EQ. nsd+1) ie = ie - 1

!     We compute the fluid-structure interaction force on the IB domain
!     and assemble the residue on the background mesh. The test or the
!     weighting function spaces used in the weak form for residue is
!     based on the background fluid mesh, while the integration is
!     performed on the IB mesh in its reference configuration.
!     Loop over all IB mesh
      ib%R = 0._RKIND
      DO iM=1, ib%nMsh
         eNoNb = ib%msh(iM)%eNoN
         nFn   = MAX(ib%msh(iM)%nFn, 1)
         ALLOCATE(Nb(eNoNb), Nbx(nsd,eNoNb), xbl(nsd,eNoNb),
     2      ul(nsd,eNoNb), fN(nsd,nFn))
!        Loop over each trace of IB integration points
         DO i=1, ib%msh(iM)%trc%n
            e  = ib%msh(iM)%trc%gE(1,i)
            g  = ib%msh(iM)%trc%gE(2,i)
            Ec = ib%msh(iM)%trc%ptr(1,i)
            jM = ib%msh(iM)%trc%ptr(2,i)

!           Get IB domain
            ib%cDmn = IB_DOMAIN(ib%msh(iM), ib%cEq, e)

!           Get IB shape functions at the integration point
            Nb = ib%msh(iM)%N(:,g)
            IF (ib%msh(iM)%eType .EQ. eType_NRB)
     2         CALL NRBNNX(ib%msh(iM), e)

!           Transfer to element-level local arrays
            fN = 0._RKIND
            DO b=1, eNoNb
               Bc = ib%msh(iM)%IEN(b,e)
               xbl(:,b) = ib%x(:,Bc)
               ul(:,b)  = ib%Un(:,Bc)
               IF (ALLOCATED(ib%msh(iM)%fN)) THEN
                  DO j=1, nFn
                     fN(:,j) = ib%msh(iM)%fN((j-1)*nsd+1:j*nsd,e)
                  END DO
               END IF
            END DO

!           Compute shapefunction gradients and element Jacobian in the
!           reference configuration. The shapefunction gradients will be
!           used to compute deformation gradient tensor
            CALL GNN(eNoNb, nsd, ib%msh(iM)%Nx(:,:,g), xbl, Nbx, Jac,
     2         Gmat)
            IF (ISZERO(Jac)) err = "Jac < 0 @ element "//e

!           Scaled Gauss weight
            w = ib%msh(iM)%w(g) * Jac

!           Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nxi(nsd,eNoN), Nx(nsd,eNoN), xl(nsd,eNoN),
     2         al(nsd,eNoN), yl(nsd+1,eNoN), lR(nsd+1,eNoN))
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               al(:,a) = Ag(is:ie,Ac)
               yl(:,a) = Yg(is:ie+1,Ac)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Dg(ie+2:ie+nsd+1,Ac)
            END DO

!           To compute test/weighting function and its gradient, get the
!           trace of the IB integration point and the shapefunctions
!           with respect to the background mesh
!           Get IB integration point
            xp = 0._RKIND
            DO b=1, eNoNb
               Bc = ib%msh(iM)%IEN(b,e)
               xp(:) = xp(:) + Nb(b)*(xbl(:,b) + ib%Un(:,Bc))
            END DO

!           Initialize parametric coordinate
            xi = 0._RKIND
            DO a=1, msh(jM)%nG
               xi = xi + msh(jM)%xi(:,a)
            END DO
            xi = xi / REAL(msh(jM)%nG, KIND=RKIND)

!           Find shape functions and derivatives on the background mesh
!           at the integration point.
            CALL GETNNX(msh(jM)%eType, eNoN, xl, msh(jM)%xib,
     2         msh(jM)%Nb, xp, xi, N, Nxi)

!           Get the shapefunction derivatives in physical space
            CALL GNN(eNoN, nsd, Nxi, xl, Nx, rt, Gmat)

!           Compute the local residue due to IB-FSI forcing
            CALL IB_CALCFFSIL(eNoNb, eNoN, nFn, w, Jac, Nb, Nbx, N, Nx,
     2         Gmat, al, yl, ul, fN, lR)

!           Assemble to global residue
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               ib%R(:,Ac) = ib%R(:,Ac) + lR(:,a)
            END DO

            DEALLOCATE(N, Nxi, Nx, xl, al, yl, lR)
         END DO
         DEALLOCATE(Nb, Nbx, xbl, ul, fN)
      END DO

      CALL COMMU(ib%R)

      RETURN
      END SUBROUTINE IB_CALCFFSI
!--------------------------------------------------------------------
!     Compute the 3D FSI force due to IB in reference configuration
      SUBROUTINE IB_CALCFFSIL(eNoNb, eNoN, nFn, w, Je, Nb, Nbx, N, Nx,
     2   Gm, al, yl, ul, fN, lR)
      USE COMMOD
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNb, eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, Je, Nb(eNoNb), Nbx(nsd,eNoNb),
     2   N(eNoN), Nx(nsd,eNoN), Gm(nsd,nsd), al(nsd,eNoN),
     3   yl(nsd+1,eNoN), ul(nsd,eNoNb), fN(nsd,nFn)
      REAL(KIND=RKIND), INTENT(OUT) :: lR(nsd+1,eNoN)

      lR = 0._RKIND
      IF (nsd .EQ. 3) THEN
         CALL IB_IFEM3D(eNoNb, eNoN, nFn, w, Nb, Nbx, N, Nx, al, yl,
     2      ul, fN, lR)
      ELSE
         CALL IB_IFEM2D(eNoNb, eNoN, nFn, w, Nb, Nbx, N, Nx, al, yl,
     2      ul, fN, lR)
      END IF

      RETURN
      END SUBROUTINE IB_CALCFFSIL
!--------------------------------------------------------------------
!     Compute the 3D FSI force due to IB in reference configuration
      SUBROUTINE IB_IFEM3D(eNoNb, eNoN, nFn, w, Nb, Nbx, N, Nx, al, yl,
     2   ul, fN, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNb, eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, Nb(eNoNb), Nbx(3,eNoNb),
     2   N(eNoN), Nx(3,eNoN), al(3,eNoN), yl(4,eNoN), ul(3,eNoNb),
     3   fN(3,nFn)
      REAL(KIND=RKIND), INTENT(OUT) :: lR(4,eNoN)

      INTEGER(KIND=IKIND) :: a, b, iEq, iDmn
      REAL(KIND=RKIND) :: ya, rt, fb(3), vd(3), v(3), vx(3,3), vVx(3),
     2   rV(3), rM(3,3)

!     Solid domain parameters
      REAL(KIND=RKIND) :: Jac, rho_s, mu_s, eM_s, nu_s, F(3,3), Fi(3,3),
     2   S(3,3), PFt(3,3), CC(3,3,3,3)
      TYPE(stModelType) :: stModel

!     Fluid domain parameters
      REAL(KIND=RKIND) :: rho_f, mu_f, gam, es(3,3), mu_x

!     Define IB struct parameters
      iEq     = ib%cEq
      iDmn    = ib%cDmn
      stModel = eq(iEq)%dmnIB(iDmn)%stM
      rho_s   = eq(iEq)%dmnIB(iDmn)%prop(solid_density)
      mu_s    = eq(iEq)%dmnIB(iDmn)%prop(solid_viscosity)
      eM_s    = eq(iEq)%dmnIB(iDmn)%prop(elasticity_modulus)
      nu_s    = eq(iEq)%dmnIB(iDmn)%prop(poisson_ratio)

!     Define fluid parameters
      rho_f   = eq(iEq)%dmn(cDmn)%prop(fluid_density)

!     Body force
      fb(1)   = eq(iEq)%dmnIB(iDmn)%prop(f_x)
      fb(2)   = eq(iEq)%dmnIB(iDmn)%prop(f_y)
      fb(3)   = eq(iEq)%dmnIB(iDmn)%prop(f_z)

!     Active strain/stress parameter (dummy value)
      ya      = 0._RKIND

!     Inertia and body force
      v       = 0._RKIND
      vd      = -fb
      vx      = 0._RKIND
      DO a=1, eNoN
!        Velocity, v
         v(1) = v(1)  + N(a)*yl(1,a)
         v(2) = v(2)  + N(a)*yl(2,a)
         v(3) = v(3)  + N(a)*yl(3,a)

!        Acceleration (dv_i/dt)
         vd(1) = vd(1) + N(a)*al(1,a)
         vd(2) = vd(2) + N(a)*al(2,a)
         vd(3) = vd(3) + N(a)*al(3,a)

!        Grad v (v_i,j) w.r.t. current coordinates
         vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a)
         vx(1,2) = vx(1,2) + Nx(2,a)*yl(1,a)
         vx(1,3) = vx(1,3) + Nx(3,a)*yl(1,a)

         vx(2,1) = vx(2,1) + Nx(1,a)*yl(2,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(2,a)
         vx(2,3) = vx(2,3) + Nx(3,a)*yl(2,a)

         vx(3,1) = vx(3,1) + Nx(1,a)*yl(3,a)
         vx(3,2) = vx(3,2) + Nx(2,a)*yl(3,a)
         vx(3,3) = vx(3,3) + Nx(3,a)*yl(3,a)
      END DO

!     Deformation gradient tensor
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      F(3,3) = 1._RKIND
      DO b=1, eNoNb
!        Deformation tensor, F_{iI}. Here the shape function derivatives
!        are w.r.t. the reference coordinates
         F(1,1) = F(1,1) + Nbx(1,b)*ul(1,b)
         F(1,2) = F(1,2) + Nbx(2,b)*ul(1,b)
         F(1,3) = F(1,3) + Nbx(3,b)*ul(1,b)

         F(2,1) = F(2,1) + Nbx(1,b)*ul(2,b)
         F(2,2) = F(2,2) + Nbx(2,b)*ul(2,b)
         F(2,3) = F(2,3) + Nbx(3,b)*ul(2,b)

         F(3,1) = F(3,1) + Nbx(1,b)*ul(3,b)
         F(3,2) = F(3,2) + Nbx(2,b)*ul(3,b)
         F(3,3) = F(3,3) + Nbx(3,b)*ul(3,b)
      END DO
      Jac = MAT_DET(F, nsd)
      Fi  = MAT_INV(F, nsd)

!     V. grad V
      vVx(1) = v(1)*vx(1,1) + v(2)*vx(1,2) + v(3)*vx(1,3)
      vVx(2) = v(1)*vx(2,1) + v(2)*vx(2,2) + v(3)*vx(2,3)
      vVx(3) = v(1)*vx(3,1) + v(2)*vx(3,2) + v(3)*vx(3,3)

!     Strain rate tensor 2*e_ij := (v_ij + v_ji)
      es(1,1) = vx(1,1) + vx(1,1)
      es(2,2) = vx(2,2) + vx(2,2)
      es(3,3) = vx(3,3) + vx(3,3)

      es(1,2) = vx(1,2) + vx(2,1)
      es(1,3) = vx(1,3) + vx(3,1)
      es(2,3) = vx(2,3) + vx(3,2)

      es(2,1) = es(1,2)
      es(3,1) = es(1,3)
      es(3,2) = es(2,3)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2) + es(1,3)*es(1,3) +
     2      es(2,1)*es(2,1) + es(2,2)*es(2,2) + es(2,3)*es(2,3) +
     3      es(3,1)*es(3,1) + es(3,2)*es(3,2) + es(3,3)*es(3,3)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(iEq)%dmn(cDmn), gam, mu_f, rt, mu_x)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
!     Note that S = Siso + Svol. For incompressible solids, Svol = 0
!     For compressible solids, Svol = J*p*Cinv where p is computed
!     based on the dilational penalty model
      CALL GETPK2CC(stModel, F, nFn, fN, ya, S, CC)

!     1st Piola-Kirchhoff stress, P = F * S
      rM(1,1)  = F(1,1)*S(1,1)  + F(1,2)*S(2,1)  + F(1,3)*S(3,1)
      rM(1,2)  = F(1,1)*S(1,2)  + F(1,2)*S(2,2)  + F(1,3)*S(3,2)
      rM(1,3)  = F(1,1)*S(1,3)  + F(1,2)*S(2,3)  + F(1,3)*S(3,3)

      rM(2,1)  = F(2,1)*S(1,1)  + F(2,2)*S(2,1)  + F(2,3)*S(3,1)
      rM(2,2)  = F(2,1)*S(1,2)  + F(2,2)*S(2,2)  + F(2,3)*S(3,2)
      rM(2,3)  = F(2,1)*S(1,3)  + F(2,2)*S(2,3)  + F(2,3)*S(3,3)

      rM(3,1)  = F(3,1)*S(1,1)  + F(3,2)*S(2,1)  + F(3,3)*S(3,1)
      rM(3,2)  = F(3,1)*S(1,2)  + F(3,2)*S(2,2)  + F(3,3)*S(3,2)
      rM(3,3)  = F(3,1)*S(1,3)  + F(3,2)*S(2,3)  + F(3,3)*S(3,3)

!     P * F_transpose
      PFt(1,1) = rM(1,1)*F(1,1) + rM(1,2)*F(1,2) + rM(1,3)*F(1,3)
      PFt(1,2) = rM(1,1)*F(2,1) + rM(1,2)*F(2,2) + rM(1,3)*F(2,3)
      PFt(1,3) = rM(1,1)*F(3,1) + rM(1,2)*F(3,2) + rM(1,3)*F(3,3)

      PFt(2,1) = rM(2,1)*F(1,1) + rM(2,2)*F(1,2) + rM(2,3)*F(1,3)
      PFt(2,2) = rM(2,1)*F(2,1) + rM(2,2)*F(2,2) + rM(2,3)*F(2,3)
      PFt(2,3) = rM(2,1)*F(3,1) + rM(2,2)*F(3,2) + rM(2,3)*F(3,3)

      PFt(3,1) = rM(3,1)*F(1,1) + rM(3,2)*F(1,2) + rM(3,3)*F(1,3)
      PFt(3,2) = rM(3,1)*F(2,1) + rM(3,2)*F(2,2) + rM(3,3)*F(2,3)
      PFt(3,3) = rM(3,1)*F(3,1) + rM(3,2)*F(3,2) + rM(3,3)*F(3,3)

!     Stress contribution
      rM(1,1) = (mu_f-mu_s)*Jac*es(1,1) - PFt(1,1)
      rM(1,2) = (mu_f-mu_s)*Jac*es(1,2) - PFt(1,2)
      rM(1,3) = (mu_f-mu_s)*Jac*es(1,3) - PFt(1,3)

      rM(2,1) = (mu_f-mu_s)*Jac*es(2,1) - PFt(2,1)
      rM(2,2) = (mu_f-mu_s)*Jac*es(2,2) - PFt(2,2)
      rM(2,3) = (mu_f-mu_s)*Jac*es(2,3) - PFt(2,3)

      rM(3,1) = (mu_f-mu_s)*Jac*es(3,1) - PFt(3,1)
      rM(3,2) = (mu_f-mu_s)*Jac*es(3,2) - PFt(3,2)
      rM(3,3) = (mu_f-mu_s)*Jac*es(3,3) - PFt(3,3)

!     Inertia contribution
      rV(1) = (rho_f*Jac-rho_s)*(vd(1) + vVx(1))
      rV(2) = (rho_f*Jac-rho_s)*(vd(2) + vVx(2))
      rV(3) = (rho_f*Jac-rho_s)*(vd(3) + vVx(3))

!     Residue
      DO a=1, eNoN
         lR(1,a) = w*(N(a)*rV(1) + Nx(1,a)*rM(1,1) + Nx(2,a)*rM(1,2)
     2      + Nx(3,a)*rM(1,3))
         lR(2,a) = w*(N(a)*rV(2) + Nx(1,a)*rM(2,1) + Nx(2,a)*rM(2,2)
     2      + Nx(3,a)*rM(2,3))
         lR(3,a) = w*(N(a)*rV(3) + Nx(1,a)*rM(3,1) + Nx(2,a)*rM(3,2)
     2      + Nx(3,a)*rM(3,3))
      END DO

      RETURN
      END SUBROUTINE IB_IFEM3D
!--------------------------------------------------------------------
!     Compute the 2D FSI force due to IB in reference configuration
      SUBROUTINE IB_IFEM2D(eNoNb, eNoN, nFn, w, Nb, Nbx, N, Nx, al, yl,
     2   ul, fN, lR)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: eNoNb, eNoN, nFn
      REAL(KIND=RKIND), INTENT(IN) :: w, Nb(eNoNb), Nbx(2,eNoNb),
     2   N(eNoN), Nx(2,eNoN), al(2,eNoN), yl(3,eNoN), ul(2,eNoNb),
     3   fN(2,nFn)
      REAL(KIND=RKIND), INTENT(OUT) :: lR(3,eNoN)

      INTEGER(KIND=IKIND) :: a, b, iEq, iDmn
      REAL(KIND=RKIND) :: ya, rt, fb(2), vd(2), v(2), vx(2,2), vVx(2),
     2   rV(2), rM(2,2)

!     Solid domain parameters
      REAL(KIND=RKIND) :: Jac, rho_s, mu_s, eM_s, nu_s, F(2,2), Fi(2,2),
     2   S(2,2), PFt(2,2), CC(2,2,2,2)
      TYPE(stModelType) :: stModel

!     Fluid domain parameters
      REAL(KIND=RKIND) :: rho_f, mu_f, gam, es(2,2), mu_x

!     Define IB struct parameters
      iEq     = ib%cEq
      iDmn    = ib%cDmn
      stModel = eq(iEq)%dmnIB(iDmn)%stM
      rho_s   = eq(iEq)%dmnIB(iDmn)%prop(solid_density)
      mu_s    = eq(iEq)%dmnIB(iDmn)%prop(solid_viscosity)
      eM_s    = eq(iEq)%dmnIB(iDmn)%prop(elasticity_modulus)
      nu_s    = eq(iEq)%dmnIB(iDmn)%prop(poisson_ratio)

!     Define fluid parameters
      rho_f   = eq(iEq)%dmn(cDmn)%prop(fluid_density)

!     Body force
      fb(1)   = eq(iEq)%dmnIB(iDmn)%prop(f_x)
      fb(2)   = eq(iEq)%dmnIB(iDmn)%prop(f_y)

!     Active strain/stress parameter (dummy value)
      ya      = 0._RKIND

!     Inertia and body force
      v       = 0._RKIND
      vd      = -fb
      vx      = 0._RKIND
      DO a=1, eNoN
!        Velocity, v
         v(1) = v(1)  + N(a)*yl(1,a)
         v(2) = v(2)  + N(a)*yl(2,a)

!        Acceleration (dv_i/dt)
         vd(1) = vd(1) + N(a)*al(1,a)
         vd(2) = vd(2) + N(a)*al(2,a)

!        Grad v (v_i,j) w.r.t. current coordinates
         vx(1,1) = vx(1,1) + Nx(1,a)*yl(1,a)
         vx(1,2) = vx(1,2) + Nx(2,a)*yl(1,a)
         vx(2,1) = vx(2,1) + Nx(1,a)*yl(2,a)
         vx(2,2) = vx(2,2) + Nx(2,a)*yl(2,a)
      END DO

!     Deformation gradient tensor
      F      = 0._RKIND
      F(1,1) = 1._RKIND
      F(2,2) = 1._RKIND
      DO b=1, eNoNb
!        Deformation tensor, F_{iI}. Here the shape function derivatives
!        are w.r.t. the reference coordinates
         F(1,1) = F(1,1) + Nbx(1,b)*ul(1,b)
         F(1,2) = F(1,2) + Nbx(2,b)*ul(1,b)
         F(2,1) = F(2,1) + Nbx(1,b)*ul(2,b)
         F(2,2) = F(2,2) + Nbx(2,b)*ul(2,b)
      END DO
      Jac = MAT_DET(F, nsd)
      Fi  = MAT_INV(F, nsd)

!     V. grad V
      vVx(1) = v(1)*vx(1,1) + v(2)*vx(1,2)
      vVx(2) = v(1)*vx(2,1) + v(2)*vx(2,2)

!     Strain rate tensor 2*e_ij := (v_ij + v_ji)
      es(1,1) = vx(1,1) + vx(1,1)
      es(2,2) = vx(2,2) + vx(2,2)
      es(1,2) = vx(1,2) + vx(2,1)
      es(2,1) = es(1,2)

!     Shear-rate := (2*e_ij*e_ij)^.5
      gam = es(1,1)*es(1,1) + es(1,2)*es(1,2)
     2    + es(2,1)*es(2,1) + es(2,2)*es(2,2)
      gam = SQRT(0.5_RKIND*gam)

!     Compute viscosity based on shear-rate and chosen viscosity model
      CALL GETVISCOSITY(eq(iEq)%dmn(cDmn), gam, mu_f, rt, mu_x)

!     2nd Piola-Kirchhoff tensor (S) and material stiffness tensor (CC)
!     Note that S = Siso + Svol. For incompressible solids, Svol = 0
!     For compressible solids, Svol = J*p*Cinv where p is computed
!     based on the dilational penalty model
      CALL GETPK2CC(stModel, F, nFn, fN, ya, S, CC)

!     1st Piola-Kirchhoff stress, P = F * S
      rM(1,1)  = F(1,1)*S(1,1)  + F(1,2)*S(2,1)
      rM(1,2)  = F(1,1)*S(1,2)  + F(1,2)*S(2,2)
      rM(2,1)  = F(2,1)*S(1,1)  + F(2,2)*S(2,1)
      rM(2,2)  = F(2,1)*S(1,2)  + F(2,2)*S(2,2)

!     P * F_transpose
      PFt(1,1) = rM(1,1)*F(1,1) + rM(1,2)*F(1,2)
      PFt(1,2) = rM(1,1)*F(2,1) + rM(1,2)*F(2,2)
      PFt(2,1) = rM(2,1)*F(1,1) + rM(2,2)*F(1,2)
      PFt(2,2) = rM(2,1)*F(2,1) + rM(2,2)*F(2,2)

!     Stress contribution
      rM(1,1) = (mu_f-mu_s)*Jac*es(1,1) - PFt(1,1)
      rM(1,2) = (mu_f-mu_s)*Jac*es(1,2) - PFt(1,2)
      rM(2,1) = (mu_f-mu_s)*Jac*es(2,1) - PFt(2,1)
      rM(2,2) = (mu_f-mu_s)*Jac*es(2,2) - PFt(2,2)

!     Inertia contribution
      rV(1) = (rho_f*Jac-rho_s)*(vd(1) + vVx(1))
      rV(2) = (rho_f*Jac-rho_s)*(vd(2) + vVx(2))

!     Residue
      DO a=1, eNoN
         lR(1,a) = w*(N(a)*rV(1) + Nx(1,a)*rM(1,1) + Nx(2,a)*rM(1,2))
         lR(2,a) = w*(N(a)*rV(2) + Nx(1,a)*rM(2,1) + Nx(2,a)*rM(2,2))
      END DO

      RETURN
      END SUBROUTINE IB_IFEM2D
!####################################################################
!     Add contribution from IB to the residue (RHS)
      SUBROUTINE IB_CONSTRUCT()
      USE COMMOD
      IMPLICIT NONE

      INTEGER(KIND=IKIND) a

      DO a=1, tnNo
         R(1:nsd+1,a) = R(1:nsd+1,a) - ib%R(:,a)
      END DO

      RETURN
      END SUBROUTINE IB_CONSTRUCT
!####################################################################
!     Poject fluid velocity, pressure and acceleration onto IB
      SUBROUTINE IB_PROJFVAR(Yg, Dg)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Yg(tDof,tnNo), Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) :: a, b, e, g, i, is, ie, Ac, Bc, Ec, iM, jM,
     2   eNoN, eNoNb
      REAL(KIND=RKIND) :: w, Jac, xi(nsd), xp(nsd), yp(nsd+1),
     2   Gmat(nsd,nsd)

      REAL(KIND=RKIND), ALLOCATABLE :: Nb(:), Nbx(:,:), xbl(:,:), N(:),
     2   Nx(:,:), xl(:,:), yl(:,:), sA(:), Uo(:,:)

      is = eq(ib%cEq)%s
      ie = eq(ib%cEq)%e
      IF (ie .EQ. nsd+1) ie = ie - 1

!     Copy ib%Un
      ALLOCATE(sA(ib%tnNo), Uo(nsd,ib%tnNo))
      sA = 0._RKIND
      Uo = ib%Un

!     Initialize to 0
      ib%Yn = 0._RKIND
      ib%Un = 0._RKIND

!     Use L2 projection with mass lumping to project flow variables from
!     background fluid mesh to IB
      DO iM=1, ib%nMsh
         eNoNb = ib%msh(iM)%eNoN
         ALLOCATE(Nb(eNoNb), Nbx(nsd,eNoNb), xbl(nsd,eNoNb))
!        Loop over each trace, as we need to first interpolate flow var
!        at the IB integration points based on its trace
         DO i=1, ib%msh(iM)%trc%n
            e  = ib%msh(iM)%trc%gE(1,i)
            g  = ib%msh(iM)%trc%gE(2,i)
            Ec = ib%msh(iM)%trc%ptr(1,i)
            jM = ib%msh(iM)%trc%ptr(2,i)

!           Transfer to local arrays: IB mesh variables
            Nb = ib%msh(iM)%N(:,g)
            IF (ib%msh(iM)%eType .EQ. eType_NRB)
     2         CALL NRBNNX(ib%msh(iM), e)
            DO b=1, eNoNb
               Bc = ib%msh(iM)%IEN(b,e)
               xbl(:,b) = ib%x(:,Bc) + Uo(:,Bc)
            END DO
            CALL GNN(eNoNb, nsd, ib%msh(iM)%Nx(:,:,g), xbl, Nbx, Jac,
     2         Gmat)
            IF (ISZERO(Jac)) err = " Jac < 0 @ element "//e
            w = ib%msh(iM)%w(g) * Jac

!           Transfer to local arrays: background mesh variables
            eNoN = msh(jM)%eNoN
            ALLOCATE(N(eNoN), Nx(nsd,eNoN), xl(nsd,eNoN),
     2         yl(nsd+1,eNoN))
            yl = 0._RKIND
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               yl(:,a) = Yg(is:ie+1,Ac)
               xl(:,a) = x(:,Ac)
               IF (mvMsh) xl(:,a) = xl(:,a) + Dg(ie+2:ie+nsd+1,Ac)
            END DO

!           Coordinates of the integration point
            xp = 0._RKIND
            DO b=1, eNoNb
               xp = xp + Nb(b)*xbl(:,b)
            END DO

!           Initialize parametric coordinate
            xi = 0._RKIND
            DO a=1, msh(jM)%nG
               xi = xi + msh(jM)%xi(:,a)
            END DO
            xi = xi / REAL(msh(jM)%nG, KIND=RKIND)

!           Find shape functions and derivatives on the background mesh
!           at the integration point.
            CALL GETNNX(msh(jM)%eType, eNoN, xl, msh(jM)%xib,
     2         msh(jM)%Nb, xp, xi, N, Nx)

!           Use the computed shape functions to interpolate flow var at
!           the IB integration point
            yp = 0._RKIND
            DO a=1, eNoN
               Ac = msh(jM)%IEN(a,Ec)
               yp = yp + N(a)*yl(:,a)
            END DO

!           Project flow variables to IB nodes
            DO b=1, eNoNb
               Bc = ib%msh(iM)%IEN(b,e)
               ib%Yn(:,Bc) = ib%Yn(:,Bc) + w*Nb(b)*yp(:)
               sA(Bc) = sA(Bc) + w*Nb(b)
            END DO
            DEALLOCATE(N, Nx, xl, yl)
         END DO
         DEALLOCATE(Nb, Nbx, xbl)
      END DO

!     Synchronize Yb across all the processes
      ib%callD(3) = CPUT()
      CALL IB_SYNC(ib%Yn)
      CALL IB_SYNC(sA)
      ib%callD(3) = CPUT() - ib%callD(3)

      DO a=1, ib%tnNo
         IF (.NOT.ISZERO(sA(a))) THEN
            ib%Yn(:,a) = ib%Yn(:,a) / sA(a)
         END IF
         DO i=1, nsd
            ib%Un(i,a) = Uo(i,a) + ib%Yn(i,a)*dt
         END DO
      END DO

      DEALLOCATE(sA, Uo)

      RETURN
      END SUBROUTINE IB_PROJFVAR
!####################################################################
!     Write IB call duration
      SUBROUTINE IB_OUTR()
      USE COMMOD
      IMPLICIT NONE

      REAL(KIND=RKIND) rtmp
      CHARACTER(LEN=stdL) sOut

      std = REPEAT("-",55)
      WRITE(sOut,'(F6.2)') ib%callD(1)
      WRITE(sOut,'(A)') " IB call duration: "//TRIM(sOut)//' sec'
      rtmp = 100._RKIND*ib%callD(3)/ib%callD(1)
      WRITE(sOut,'(A)') TRIM(sOut)//" (comm."//
     2   STR(NINT(rtmp, KIND=IKIND),3)//"%)"
      rtmp = 100._RKIND*ib%callD(2)/ib%callD(1)
      WRITE(sOut,'(A)') TRIM(sOut)//", (updt."//
     2   STR(NINT(rtmp, KIND=IKIND),3)//"%)"
      std = sOut
      std = REPEAT("-",55)

      RETURN
      END SUBROUTINE IB_OUTR
!####################################################################
!     Debugs FSI force on the IB
      SUBROUTINE DEBUGIBR(incNd)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: incNd(ib%tnNo)

      INTEGER(KIND=IKIND) :: a, i, Ac, iM, fid
      REAL(KIND=RKIND) :: s, lo, hi, av
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: lI(:), gI(:)
      REAL(KIND=RKIND), ALLOCATABLE :: lR(:,:), gR(:,:)

!     DEBUG ib%R
      fid = 1289
      IF (cm%mas()) THEN
         WRITE(fName,'(A)') TRIM(appPath)//"dbg_ibR_hist.dat"
         IF (cTS .EQ. 1) THEN
            OPEN(fid,FILE=TRIM(fName))
            CLOSE(fid,STATUS='DELETE')
         END IF
         OPEN(fid,FILE=TRIM(fName),POSITION='APPEND')
         WRITE(fid,'(A)',ADVANCE='NO') STR(cTS)
      END IF

      DO iM=1, nMsh
         ALLOCATE(lR(nsd+1,msh(iM)%nNo), lI(msh(iM)%nNo))
         IF (cm%mas()) THEN
            ALLOCATE(gR(nsd+1,msh(iM)%gnNo), gI(msh(iM)%gnNo))
         ELSE
            ALLOCATE(gR(0,0), gI(0))
         END IF
         lR = 0._RKIND
         lI = 0
         DO a=1, msh(iM)%nNo
            Ac = msh(iM)%gN(a)
            lR(:,a) = ib%R(:,Ac)
            lI(a)   = incNd(Ac)
         END DO
         gR = GLOBAL(msh(iM), lR)
         gI = GLOBAL(msh(iM), lI)
         DEALLOCATE(lR, lI)

         IF (cm%mas()) THEN
!           first, momentum residue
            av =  0._RKIND
            lo =  1.E+6_RKIND
            hi = -1.E+6_RKIND
            DO a=1, msh(iM)%gnNo
               IF (gI(a) .EQ. 0) CYCLE
               s = 0._RKIND
               DO i=1, nsd
                  s = s + gR(i,a)**2._RKIND
               END DO
               s = SQRT(s)
               av = av + s
               IF (s .LT. lo) lo = s
               IF (s .GT. hi) hi = s
            END DO
            av = av / REAL(SUM(gI(:)), KIND=RKIND)
            WRITE(fid,'(A)',ADVANCE='NO') " "//STR(av)//" "//STR(hi-lo)

!           Pressure/continuity
            av =  0._RKIND
            lo =  1.E+6_RKIND
            hi = -1.E+6_RKIND
            DO a=1, msh(iM)%gnNo
               IF (gI(a) .EQ. 0) CYCLE
               s = gR(nsd+1,a)
               av = av + s
               IF (s .LT. lo) lo = s
               IF (s .GT. hi) hi = s
            END DO
            av = av / REAL(SUM(gI(:)), KIND=RKIND)
            WRITE(fid,'(A)') " "//STR(av)//" "//STR(hi-lo)
            CLOSE(fid)
         END IF
         DEALLOCATE(gR, gI)
      END DO

      iM = cm%reduce(iM)

      RETURN
      END SUBROUTINE DEBUGIBR
!--------------------------------------------------------------------
!     Debugs IB mesh traces
      SUBROUTINE DEBUGIBMSHTRC(lM)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM

      INTEGER(KIND=IKIND) :: i, a, e, g, n, Ac, Ec, iM, fid, ierr
      CHARACTER(LEN=stdL) :: fName
      REAL(KIND=RKIND) xp(nsd)

      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disps(:), lE(:),
     2   gE(:), eptr(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:)

      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lM%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      n = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(lE(4*lM%trc%n), gE(4*n))
      DO i=1, lM%trc%n
         lE(4*i-3) = lM%trc%gE(1,i)
         lE(4*i-2) = lM%trc%gE(2,i)
         lE(4*i-1) = lM%trc%ptr(1,i)
         lE(4*i)   = lM%trc%ptr(2,i)
      END DO

      i = 4*lM%trc%n
      CALL MPI_GATHERV(lE, i, mpint, gE, sCount, disps, mpint, master,
     2   cm%com(), ierr)

      DEALLOCATE(lE, disps, sCount)

      WRITE(fName,'(A)') TRIM(appPath)//"dbg_ib_trc_"//STR(cTS)//".dat"
      IF (cm%mas()) THEN
         ALLOCATE(eptr(2,lM%nG,lM%nEl))
         eptr = 0
         DO i=1, n
            e  = gE(4*i-3)
            g  = gE(4*i-2)
            Ec = gE(4*i-1)
            iM = gE(4*i)
            ePtr(1,g,e) = Ec
            ePtr(2,g,e) = iM
         END DO

         fid = 1289
         OPEN(fid,FILE=TRIM(fName))
         WRITE(fid,'(A)') "List of failed traces on mesh: "//
     2      TRIM(lM%name)
         ALLOCATE(xl(nsd,lM%eNoN))
         DO e=1, lM%nEl
            DO a=1, lM%eNoN
               Ac = lM%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lM%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
                  xp = 0._RKIND
                  DO a=1, lM%eNoN
                     xp(:) = xp(:) + lM%N(a,g)*xl(:,a)
                  END DO
                  WRITE(fid,'(2X,A)',ADVANCE='NO') STR(g)//" "//STR(e)
                  DO i=1, nsd
                     WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
                  END DO
                  WRITE(fid,'(A)')
               END IF
            END DO
         END DO
         DEALLOCATE(xl)
         CLOSE(fid)
      ELSE
         ALLOCATE(ePtr(0,0,0))
      END IF
      DEALLOCATE(gE, ePtr)

      i = cm%reduce(i)
      err = "ERROR: Failed to detect all the traces on "//
     2   TRIM(lM%name)//". See "//TRIM(fName)//" for more information."

      RETURN
      END SUBROUTINE DEBUGIBMSHTRC
!--------------------------------------------------------------------
!     Debugs IB face traces
      SUBROUTINE DEBUGIBFATRC(lFa)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(faceType), INTENT(IN) :: lFa

      INTEGER(KIND=IKIND) :: i, a, e, g, n, Ac, Ec, iM, fid, ierr
      REAL(KIND=RKIND) xp(nsd)
      CHARACTER(LEN=stdL) :: fName

      INTEGER(KIND=IKIND), ALLOCATABLE :: sCount(:), disps(:), lE(:),
     2   gE(:), eptr(:,:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:)

      ALLOCATE(sCount(cm%np()), disps(cm%np()))
      sCount = 0
      disps  = 0
      i = lFa%trc%n
      CALL MPI_GATHER(i, 1, mpint, sCount, 1, mpint, master, cm%com(),
     2   ierr)

      n = SUM(sCount(:))
      sCount = 4*sCount(:)
      DO i=2, cm%np()
         disps(i) = disps(i-1) + sCount(i-1)
      END DO

      ALLOCATE(lE(4*lFa%trc%n), gE(4*n))
      DO i=1, lFa%trc%n
         lE(4*i-3) = lFa%trc%gE(1,i)
         lE(4*i-2) = lFa%trc%gE(2,i)
         lE(4*i-1) = lFa%trc%ptr(1,i)
         lE(4*i)   = lFa%trc%ptr(2,i)
      END DO

      i = 4*lFa%trc%n
      CALL MPI_GATHERV(lE, i, mpint, gE, sCount, disps, mpint, master,
     2   cm%com(), ierr)

      DEALLOCATE(lE, disps, sCount)

      WRITE(fName,'(A)') TRIM(appPath)//"dbg_ib_trc_"//STR(cTS)//".dat"
      IF (cm%mas()) THEN
         ALLOCATE(eptr(2,lFa%nG,lFa%nEl))
         eptr = 0
         DO i=1, n
            e  = gE(4*i-3)
            g  = gE(4*i-2)
            Ec = gE(4*i-1)
            iM = gE(4*i)
            ePtr(1,g,e) = Ec
            ePtr(2,g,e) = iM
         END DO

         fid = 1289
         OPEN(fid,FILE=TRIM(fName))
         WRITE(fid,'(A)') "List of failed traces on mesh: "//
     2      TRIM(lFa%name)
         ALLOCATE(xl(nsd,lFa%eNoN))
         DO e=1, lFa%nEl
            DO a=1, lFa%eNoN
               Ac = lFa%IEN(a,e)
               xl(:,a) = ib%x(:,Ac) + ib%Un(:,Ac)
            END DO
            DO g=1, lFa%nG
               Ec = ePtr(1,g,e)
               iM = ePtr(2,g,e)
               IF (Ec.EQ.0 .OR. iM.EQ.0) THEN
                  xp = 0._RKIND
                  DO a=1, lFa%eNoN
                     xp(:) = xp(:) + lFa%N(a,g)*xl(:,a)
                  END DO
                  WRITE(fid,'(2X,A)',ADVANCE='NO') STR(g)//" "//STR(e)
                  DO i=1, nsd
                     WRITE(fid,'(A)',ADVANCE='NO') " "//STR(xp(i))
                  END DO
                  WRITE(fid,'(A)')
               END IF
            END DO
         END DO
         DEALLOCATE(xl)
         CLOSE(fid)
      ELSE
         ALLOCATE(ePtr(0,0,0))
      END IF
      DEALLOCATE(gE, ePtr)

      i = cm%reduce(i)
      err = "ERROR: Failed to detect all the traces on "//
     2   TRIM(lFa%name)//". See "//TRIM(fName)//" for more information."

      RETURN
      END SUBROUTINE DEBUGIBFATRC
!####################################################################

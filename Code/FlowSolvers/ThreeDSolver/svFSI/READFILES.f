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
!     This subroutine is intended for reading the input information
!     and allocating required space for program
!
!--------------------------------------------------------------------

      SUBROUTINE READFILES

      USE COMMOD
      USE ALLFUN
      USE LISTMOD

      IMPLICIT NONE

      LOGICAL :: ltmp
      INTEGER :: i, iEq
      INTEGER :: tArray(8)
      REAL(KIND=8) :: roInf
      CHARACTER(LEN=8) :: date
      CHARACTER(LEN=stdL) :: ctmp
      CHARACTER(LEN=stdL) :: mfsIn
      TYPE(listType) :: list
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) :: fTmp
      SAVE mfsIn, roInf

      IF (.NOT.resetSim) THEN
         ctmp = " P#"//STR(cm%id())
         IF (cm%mas()) ctmp = ""
         CALL io%o%new(CHNL_O,tag=ctmp,oTS=cm%mas(),oTF=cm%mas())
         CALL io%e%new(CHNL_E,tag=ctmp)
         CALL io%w%new(CHNL_W,tag=ctmp)
         CALL io%d%new(CHNL_D,tag=ctmp)

         std => io%o
         err => io%e
         wrn => io%w
         dbg => io%d
      END IF

      IF (cm%slv()) RETURN

!     These are the default values
      IF (.NOT.resetSim) THEN
         pClr         = .TRUE.
         mvMsh        = .FALSE.
         stFileFlag   = .FALSE.
         stFileRepl   = .FALSE.
         saveAve      = .FALSE.
         sepOutput    = .FALSE.
         legacyFmt    = .FALSE.
         saveATS      = 1
         saveIncr     = 10
         nITs         = 0
         roInf        = 2D-1
         stFileName   = "stFile"
         iniFilePath  = ""
         stopTrigName = "STOP_SIM"
         rmsh%isReqd  = .FALSE.
         ichckIEN     = .TRUE.
         zeroAve      = .FALSE.

         useTrilinosLS         = .FALSE.
         useTrilinosAssemAndLS = .FALSE.

         i = IARGC()
         IF (i .EQ. 0) THEN
            err = "Configuration file name is required as an argument"
         ELSEIF (i .GT. 1) THEN
            err = "Too many arguments"
         END IF
         CALL GETARG(1,ctmp)
         mfsIn = ctmp
      END IF

      list = listType(mfsIn,io)

      IF (.NOT.resetSim) THEN
         ltmp = .TRUE.
         lPtr => list%get(ltmp,"Save results in a folder")
         IF (ltmp) THEN
             saveName = ""
             lPtr => list%get(saveName,"Save results in folder")
             IF (saveName == "") THEN
                appPath = STR(cm%np())//"-procs"//delimiter
             ELSE
                appPath = TRIM(saveName)//delimiter
             END IF
         END IF
         lPtr => list%get(std%oTS,"Verbose")
         lPtr => list%get(wrn%oTS,"Warning")
         lPtr => list%get(dbg%oTS,"Debug")
         std%oTF = std%oTS
         wrn%oTF = wrn%oTS
         dbg%oTF = dbg%oTS

         lPtr => list%get(pClr,"Colorful terminal output")
         lPtr => list%get(sepOutput,"Use separator in the history file")
         lPtr => list%get(stFileFlag,"Continue previous simulation",1)

         IF (.NOT.stFileFlag) REWIND(std%fId)

         CALL DATE_AND_TIME(DATE=date, VALUES=tArray)
         std = "                                               "
         std = " ----------------------------------------------"
         std = "                                               "
         std = "    SimVascular Fluid-Structure Interaction    "
         std = "                  (svFSI v"//version//")       "
         std = "                                               "
         std = " ----------------------------------------------"
         std = "                                               "
         std = " This program comes with ABSOLUTELY NO WARRANTY"
         std = " Compiled on "//date(5:6)//"-"//date(7:8)//
     2      "-"//date(1:4)
         std = "                                               "

         IF (tArray(1).GE.expDate(1) .AND. tArray(2).GE.expDate(2)
     2      .AND. tArray(3).GE.expDate(3)) err = "License is expired"

         IF (cm%seq() .AND. cm%nT().EQ.1) THEN
            std = " Running sequentially"
         ELSE
            IF (.NOT.cm%seq()) std = " Running in parallel with "//
     2         cm%np()//" tasks"
            IF (cm%nT() .NE. 1) std = " Multi-threading with "//
     2         cm%nT()//" proccesors"
         END IF

         IF (.NOT.stFileFlag) THEN
            lPtr => list%get(fTmp,"Simulation initialization file path")
            IF (ASSOCIATED(lPtr))
     2         iniFilePath = fTmp%fname
         END IF

         lPtr => list%get(nsd,"Number of spatial dimensions",
     2      1,ll=2,ul=3)
         lPtr => list%get(nTs,"Number of time steps",1,ll=1)
         lPtr => list%get(dt,"Time step size",1,lb=0D0)
         lPtr => list%get(nITs,"Number of initialization time steps",
     2      ll=0)
         lPtr => list%get(roInf,"Spectral radius of infinite time step",
     2      ll=0D0,ul=1D0)

         lPtr =>list%get(stopTrigName,
     2      "Searched file name to trigger stop")
         stopTrigName = TRIM(appPath)//stopTrigName
         lPtr => list%get(ichckIEN, "Check IEN order")

         lPtr => list%get(legacyFmt,"Use legacy file format")

         IF (legacyFmt) THEN
            wrn = "Legacy file format is NOT fully supported. "//
     2         " May lead to errors"
            ctmp = 'VTKB'
            lPtr => list%get(ctmp,"Format of saved files")
            SELECT CASE (TRIM(ctmp))
            CASE('none')
               saveFormat = saveF_none
            CASE('VTK')
               saveFormat = saveF_VTK
            CASE('VTKB')
               saveFormat = saveF_VTKB
            CASE DEFAULT
               saveFormat = saveF_VTKB
            END SELECT
         END IF

         lPtr => list%get(saveIncr,"Increment in saving files",ll=1)
         lPtr => list%get(saveATS,"Start saving after time step",
     2         ll=1)
         lPtr => list%get(saveAve,"Save averaged results")
         lPtr => list%get(zeroAve,"Start averaging from zero")
         lPtr => list%get(saveName,"Name prefix of saved files",1)
         saveName = TRIM(appPath)//saveName

         lPtr => list%get(stFileRepl,"Overwrite restart file")
         lPtr => list%get(stFileName,"Restart file name")
         stFileName = TRIM(appPath)//stFileName

         stFileIncr = saveIncr
         lPtr => list%get(stFileIncr,
     2      "Increment in saving restart files",ll=0)
         lPtr => list%get(rmsh%isReqd, "Simulation requires remeshing")

      END IF ! resetSim

!--------------------------------------------------------------------
!     Reading the mesh
      CALL READMSH(list)

!--------------------------------------------------------------------
!     Reading equations
      nEq = list%srch("Add equation",ll=1)
      std = " Number of equations: "//nEq
      ALLOCATE(eq(nEq))
      eq%roInf = roInf
      DO iEq=1, nEq
!     This is pointing to an equation list
         lPtr => list%get(ctmp,"Add equation",iEq)
         CALL READEQ(eq(iEq), lPtr, ctmp)

         IF (eq(iEq)%phys .EQ. phys_fluid) THEN
            IF (iEq .NE. 1) err = "fluid equation must come first"
         END IF
         IF (eq(iEq)%phys .EQ. phys_FSI) THEN
            IF (iEq .NE. 1) err = "FSI equation must come first"
            lPtr => list%get(ctmp,"Add equation",2)
            IF (ctmp .NE. "mesh") err = "mesh equation"//
     2         " has to be specified after FSI equation"
         END IF
         IF (rmsh%isReqd .AND. eq(1)%phys .NE. phys_FSI) THEN
            err = "Remeshing is applicable only for FSI equation"
         END IF
         lPtr => list%get(ctmp,"Add equation",1)
         IF ((eq(iEq)%phys .EQ. phys_BBO) .OR.
     2       (eq(iEq)%phys .EQ. phys_heatF)) THEN
            IF (ctmp.NE."fluid" .AND. ctmp.NE."FSI") THEN
               err = "BBO/heatF equation has to be specified after"//
     2            " fluid/FSI equation"
            END IF
         END IF
         IF (eq(iEq)%phys .EQ. phys_mesh) THEN
            IF (.NOT.mvMsh) err = "mesh equation can only"//
     2         " be specified after FSI equation"
         END IF
      END DO
      IF (.NOT.ALLOCATED(cplBC%xo)) THEN
         cplBC%nX = 0
         ALLOCATE(cplBC%xo(cplBC%nX))
      END IF

      IF (.NOT.resetSim) THEN
         CALL list%check()
         std = CLR(" Configuration is completed successfully",3)
      ELSE
         std =
     2   CLR(" Configuration is completed successfully after remesh",3)
      END IF

      CALL DESTROYLIST(list)

      RETURN
      END SUBROUTINE READFILES

!####################################################################
!     This routine reads a Eq
      SUBROUTINE READEQ(lEq, list, eqName)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(listType), INTENT(INOUT) :: list
      CHARACTER(LEN=stdL), INTENT(IN) :: eqName

      INTEGER, PARAMETER :: maxOutput = 10
      INTEGER fid, iBc, phys(2), propL(maxNProp,10),
     2   outPuts(maxOutput), nDOP(4)
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPBC
      TYPE(fileType) fTmp

      lPtr => list%get(lEq%coupled,"Coupled")
      lPtr => list%get(lEq%minItr,"Min iterations",ll=1)
      lPtr => list%get(lEq%maxItr,"Max iterations",ll=1)
      lPtr => list%get(lEq%dBr,"Residual dB reduction",ul=0D0)
      lPtr => list%get(lEq%tol,"Tolerance",ll=0D0)

!     Coupled BC stuff
!     Initializing coupled BC if neccessary
      IF (.NOT.ALLOCATED(cplBC%xo)) THEN
         lPBC => list%get(ctmp,"Couple to cplBC")
         IF (ASSOCIATED(lPBC)) THEN
            SELECT CASE(ctmp)
            CASE('N')
               cplBC%schm = cplBC_NA
            CASE('I')
               cplBC%schm = cplBC_I
            CASE('SI')
               cplBC%schm = cplBC_SI
            CASE('E')
               cplBC%schm = cplBC_E
            CASE DEFAULT
               err = "Undefined cplBC%schm: "//ctmp
            END SELECT
         END IF
         IF (cplBC%schm .NE. cplBC_NA) THEN
!     This pointing to the section containing BC info under fluid eqn
            lPtr => lPBC%get(cplBC%nX,"Number of unknowns",1,ll=0)
            ALLOCATE(cplBC%xo(cplBC%nX))

            cplBC%xo = 0D0
            lPtr => lPBC%get(fTmp,"Unknowns initialization file path")
            IF (ASSOCIATED(lPtr)) THEN
               fid = fTmp%open()
               READ (fid,*) cplBC%xo
               CLOSE (fid)
            END IF

            lPtr => lPBC%get(cplBC%commuName,
     2         "File name for 0D-3D communication")
            cplBC%commuName = TRIM(appPath)//cplBC%commuName

            lPtr => lPBC%get(cplBC%saveName,
     2         "File name for saving unknowns")
            cplBC%saveName = TRIM(appPath)//cplBC%saveName

            lPtr => lPBC%get(fTmp,"0D code file path",1)
            cplBC%binPath = fTmp%fname
         END IF
      END IF

      propL = prop_NA
      SELECT CASE (eqName)
!     FLUID Navier-Stokes solver ------------------------------------
      CASE ('fluid')
         lEq%phys = phys_fluid

         propL(1,1) = fluid_density
         propL(2,1) = viscosity
         propL(3,1) = permeability
         propL(4,1) = backflow_stab
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/8,2,3,0/)
         outPuts(1) = out_velocity
         outPuts(2) = out_pressure
         outPuts(3) = out_energyFlux
         outPuts(4) = out_acceleration
         outPuts(5) = out_WSS
         outPuts(6) = out_vorticity
         outPuts(7) = out_strainInv
         outPuts(8) = out_vortex

         CALL READLS(lSolver_NS, lEq, list)
!     HEAT FLUID advection diffusion solver -------------------------
      CASE ('heatF')
         lEq%phys = phys_heatF

         propL(1,1) = conductivity
         propL(2,1) = source_term
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/2,1,2,0/)
         outPuts(1) = out_temperature
         outPuts(2) = out_heatFlux

         CALL READLS(lSolver_GMRES, lEq, list)
!     HEAT SOLID Laplac equation solver------------------------------
      CASE ('heatS')
         lEq%phys = phys_heatS

         propL(1,1) = conductivity
         propL(2,1) = source_term
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/2,1,2,0/)
         outPuts(1) = out_temperature
         outPuts(2) = out_heatFlux

         CALL READLS(lSolver_GMRES, lEq, list)
!     LINEAR ELASTICITY equation solver------------------------------
      CASE ('lElas')
         lEq%phys = phys_lElas

         propL(1,1) = solid_density
         propL(2,1) = elasticity_modulus
         propL(3,1) = poisson_ratio
         propL(4,1) = f_x
         propL(5,1) = f_y
         IF (nsd .EQ. 3) propL(6,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/3,1,0,0/)
         outPuts(1) = out_displacement
         outPuts(2) = out_velocity
         outPuts(3) = out_acceleration

         CALL READLS(lSolver_CG, lEq, list)
!     STRUCTURAL with nonlinear displacement equation solver---------
      CASE ('struct')
         lEq%phys = phys_struct

         propL(1,1) = solid_density
         propL(2,1) = elasticity_modulus
         propL(3,1) = poisson_ratio
         propL(4,1) = f_x
         propL(5,1) = f_y
         IF (nsd .EQ. 3) propL(6,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/4,2,0,0/)
         outPuts(1) = out_displacement
         outputs(2) = out_strainInv
         outPuts(3) = out_velocity
         outPuts(4) = out_acceleration

         CALL READLS(lSolver_CG, lEq, list)
!     FLUID STRUCTUR INTERACTION equation solver---------------------
      CASE ('FSI')
         mvMsh = .TRUE.
         lEq%phys = phys_FSI

         propL(1,1) = fluid_density
         propL(2,1) = viscosity
         propL(3,1) = permeability
         propL(4,1) = backflow_stab
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z
         propL(1,2) = solid_density
         propL(2,2) = elasticity_modulus
         propL(3,2) = poisson_ratio
         propL(4,2) = damping
         propL(5,2) = f_x
         propL(6,2) = f_y
         IF (nsd .EQ. 3) propL(7,2) = f_z
         phys(1) = phys_fluid
         phys(2) = phys_struct
         CALL READDOMAIN(lEq, propL, list, phys)

         nDOP = (/10,2,4,0/)
         outPuts(1)  = out_velocity
         outPuts(2)  = out_pressure
         outPuts(3)  = out_energyFlux
         outPuts(4)  = out_absVelocity
         outPuts(5)  = out_acceleration
         outPuts(6)  = out_WSS
         outPuts(7)  = out_vorticity
         outPuts(8)  = out_vortex
         outPuts(9)  = out_displacement
         outPuts(10) = out_strainInv

         CALL READLS(lSolver_GMRES, lEq, list)

         IF (rmsh%isReqd .AND. .NOT.resetSim)
     2      CALL READRMSH(list)

!     MESH motion equation solver------------------------------------
      CASE ('mesh')
         lEq%phys = phys_mesh

         propL(1,1) = poisson_ratio
         CALL READDOMAIN(lEq, propL, list)
         lEq%dmn%prop(solid_density) = 0D0
         lEq%dmn%prop(elasticity_modulus) = 1D0

         nDOP = (/3,1,0,0/)
         outPuts(1) = out_displacement
         outPuts(2) = out_velocity
         outPuts(3) = out_acceleration

         lEq%ls%reltol = 2D-1
         CALL READLS(LS_TYPE_CG, lEq, list)

!     Basset-Boussinesq-Oseen equation solver------------------------
      CASE ('BBO')
         wrn = "you need to double check the equations/TAG"
         lEq%phys = phys_BBO

         propL(1,1) = fluid_density
         propL(2,1) = viscosity
         propL(3,1) = particle_density
         propL(4,1) = particle_diameter
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/2,1,1,0/)
         outPuts(1) = out_velocity
         outPuts(2) = out_acceleration

         CALL READLS(lSolver_GMRES, lEq, list)
      CASE DEFAULT
         err = "Equation type "//TRIM(eqName)//" is not defined"
      END SELECT
      CALL READOUTPUTS(lEq, nDOP, outPuts, list)
!--------------------------------------------------------------------
!     Searching for BCs
      lEq%nBc = list%srch("Add BC")
      ALLOCATE(lEq%bc(lEq%nBc))
      std = " Number of imposed BC for equation <"//TRIM(eqName)//
     2   ">: "//lEq%nBc
      DO iBc=1, lEq%nBc
         lPBC => list%get(ctmp,"Add BC",iBc)
         CALL FINDFACE(ctmp, lEq%bc(iBc)%iM, lEq%bc(iBc)%iFa)
         CALL READBC(lEq%bc(iBc), lPBC, lEq%phys)
      END DO

      RETURN
      CONTAINS
!--------------------------------------------------------------------
!     This routine is intended to read properties of domains
      SUBROUTINE READDOMAIN(lEq, propList, list, physI)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER, INTENT(IN) :: propList(maxNProp,10)
      TYPE (listType), TARGET, INTENT(INOUT) :: list
      INTEGER, INTENT(IN), OPTIONAL :: physI(:)

      LOGICAL flag
      INTEGER i, iDmn, iPhys, iProp, prop, nPhys
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPD
      INTEGER, ALLOCATABLE :: phys(:)

      IF (PRESENT(physI)) THEN
         nPhys = SIZE(physI)
         ALLOCATE(phys(nPhys))
         phys = physI
      ELSE
         nPhys = 1
         ALLOCATE(phys(nPhys))
         phys = lEq%phys
      END IF

      flag = .FALSE.
      lEq%nDmn = list%srch("Domain")
      IF (lEq%nDmn .EQ. 0) THEN
         lEq%nDmn = 1
         flag = .TRUE.
      END IF

      ALLOCATE(lEq%dmn(lEq%nDmn))
      DO iDmn=1, lEq%nDmn
         IF (flag) THEN
!     If no domain keywork found, we search upper level for data
            lPD => list
            lEq%dmn(iDmn)%Id = -1
         ELSE
            lPD => list%get(lEq%dmn(iDmn)%Id,"Domain",iDmn,
     2         ll=0,ul=(BIT_SIZE(dmnId)-1))
!     If following is the case, then there is only one domain that must
!     be default (Id=-1)
            IF (.NOT.ALLOCATED(msh)) err = "No mesh is read yet"
            IF (.NOT.ALLOCATED(dmnId)) err = TRIM(list%ping("Domain",
     2         lPD))//" No domain file is read yet"

            DO i=1, iDmn-1
               IF (lEq%dmn(iDmn)%Id .EQ. lEq%dmn(i)%Id) err =
     2            TRIM(list%ping("Domain",lPD))//" Repeated domain ID"
            END DO
         END IF
         IF (nPhys .GT. 1) THEN
            IF (lEq%phys .NE. phys_FSI) THEN
               err = "Correct THIS, you have an equation that"//
     2            " requires more than one block and is not FSI"
            END IF
            lPtr => lPD%get(ctmp,"Equation",1)
            SELECT CASE(TRIM(ctmp))
            CASE("fluid")
               lEq%dmn(iDmn)%phys = phys_fluid
            CASE("struct")
               lEq%dmn(iDmn)%phys = phys_struct
            CASE("lElas")
               lEq%dmn(iDmn)%phys = phys_lElas
            CASE DEFAULT
               err = TRIM(lPD%ping("Equation",lPtr))//
     2            "Equation must be fluid/struct/lElas"
            END SELECT
         ELSE
            lEq%dmn(iDmn)%phys = lEq%phys
         END IF
         DO iPhys=1, nPhys
            IF (lEq%dmn(iDmn)%phys .EQ. phys(iPhys)) EXIT
         END DO
         IF (iPhys .GT. nPhys) err = "Undefined phys is used"

         lEq%dmn(iDmn)%cModel = cModel_NA
         IF (lEq%dmn(iDmn)%phys .EQ. phys_struct) THEN
            lPtr => lPD%get(ctmp,"Constitutive model")
            SELECT CASE (TRIM(ctmp))
            CASE ("stVK","stVenantKirchhoff")
               lEq%dmn(iDmn)%cModel = cModel_stVK
            CASE ("stVK85","stVenantKirchhoffSimo85")
               lEq%dmn(iDmn)%cModel = cModel_mStVK
            CASE ("nHK","nHK91","neoHookean","neoHookeanSimo91")
               lEq%dmn(iDmn)%cModel = cModel_nHook
            CASE DEFAULT
               err = "Undefined constitutive model used"
            END SELECT
         END IF

         DO iProp=1, maxNProp
            rtmp = 0D0
            prop = propList(iProp,iPhys)
            SELECT CASE (prop)
            CASE (prop_NA)
               EXIT
            CASE (fluid_density,solid_density)
               lPtr => lPD%get(rtmp,"Density",1,lb=0D0)
            CASE (viscosity)
               lPtr => lPD%get(rtmp,"Viscosity",1,lb=0D0)
            CASE (elasticity_modulus)
               lPtr => lPD%get(rtmp,"Elasticity modulus",1,lb=0D0)
            CASE (poisson_ratio)
               lPtr => lPD%get(rtmp,"Poisson ratio",1,ll=0D0,ub=5D-1)
            CASE (conductivity)
               lPtr => lPD%get(rtmp,"Conductivity",1,ll=0D0)
            CASE (f_x)
               lPtr => lPD%get(rtmp,"Force_X")
            CASE (f_y)
               lPtr => lPD%get(rtmp,"Force_Y")
            CASE (f_z)
               lPtr => lPD%get(rtmp,"Force_Z")
            CASE (particle_density)
               lPtr => lPD%get(rtmp,"Particle density",1,lb=0D0)
            CASE (particle_diameter)
               lPtr => lPD%get(rtmp,"Particle diameter",1,ll=0D0)
            CASE (permeability)
               rtmp = HUGE(lEq%dmn(iDmn)%prop)
               lPtr => lPD%get(rtmp,"Permeability",lb=0D0)
            CASE (backflow_stab)
               rtmp = 2D-1
               lPtr => lPD%get(rtmp,"Backflow stabilization coefficient"
     2            ,ll=0D0)
            CASE (source_term)
               lPtr => lPD%get(rtmp,"Source term")
            CASE (damping)
               lPtr => lPD%get(rtmp,"Damping")
            CASE DEFAULT
               err = "Undefined properties"
            END SELECT
            lEq%dmn(iDmn)%prop(prop) = rtmp
         END DO
      END DO

      RETURN
      END SUBROUTINE READDOMAIN
!--------------------------------------------------------------------
      END SUBROUTINE READEQ
!####################################################################
!     This routine reads a BC
      SUBROUTINE READBC(lBc, list, phys)

      USE COMMOD
      USE ALLFUN
      USE LISTMOD

      IMPLICIT NONE

      TYPE(bcType), INTENT(INOUT) :: lBc
      TYPE(listType), INTENT(INOUT) :: list
      INTEGER, INTENT(IN) :: phys

      LOGICAL ltmp
      INTEGER iFa, iM, a, b, fid, i, j, Ac
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) fTmp

      INTEGER, ALLOCATABLE :: ptr(:)

!     Reading the type: Dir/Neu/Per
      lPtr => list%get(ctmp,"Type")
      SELECT CASE (ctmp)
      CASE ("Dirichlet","Dir")
         lBc%bType = IBSET(lBc%bType,bType_Dir)
      CASE ("Neumann","Neu")
         lBc%bType = IBSET(lBc%bType,bType_Neu)
         IF (phys.EQ.phys_fluid .OR. phys.EQ.phys_FSI)
     2      lBc%bType = IBSET(lBc%bType,bType_bfs)
      CASE ("Periodic","Per")
         lBc%bType = IBSET(lBc%bType,bType_per)
         err = "Periodic BC hasn't been implemented yet"
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected BC type"
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
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lBc%gt%n = j
            ALLOCATE(lBc%gt%r(j))
            ALLOCATE(lBc%gt%i(j))
            CALL FFT(fid, i, lBc%gt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
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
         lBc%bType = IBSET(lBc%bType,bType_cpl)
         cplBC%nFa = cplBC%nFa + 1
         lBc%cplBcPtr = cplBC%nFa
         IF (cplBC%schm .EQ. cplBC_NA) THEN
            err = "'Couple to cplBC' must be specified before"//
     2         " using Coupled BC"
         END IF
      CASE ('Resistance')
         lBc%bType = IBSET(lBc%bType,bType_res)
         IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "Resistance"//
     2      " is only defined for Neu BC"
         IF (phys.NE.phys_fluid .AND. phys.NE.phys_FSI) err =
     2      "Resistance is only defined for fluid/FSI equations"

         lPtr => list%get(lBc%r,"Value",1)
      CASE ('General')
         lBc%bType = IBSET(lBc%bType,bType_gen)
         lPtr =>list%get(fTmp,"Temporal and spatial values file path",1)
         fid = fTmp%open()
         READ (fid,*) i, j, a
         iM  = lBc%iM
         iFa = lBc%iFa
         IF (a .NE. msh(iM)%fa(iFa)%nNo) THEN
            err = "Number of nodes does not match between "//
     2         TRIM(msh(iM)%fa(iFa)%name)//" and "//TRIM(fTmp%fname)
         END IF
         IF (i.LT.1 .OR. i.GT.nsd) err = "0 < dof <= "//nsd//
     2      " is violated in "//TRIM(fTmp%fname)

         ALLOCATE(lBc%gm)
         ALLOCATE(lBc%gm%t(j), lBc%gm%d(i,a,j), ptr(msh(iM)%gnNo))
!     I am seting all the nodes to zero just in case a node is not set
         lBc%gm%d   = 0D0
         lBc%gm%dof = i
         lBc%gm%nTP = j
         ptr        = 0
!     Preparing the pointer array
         DO a=1, msh(iM)%fa(iFa)%nNo
            Ac = msh(iM)%fa(iFa)%gN(a)
            ptr(Ac) = a
         END DO
         DO i=1, j
            READ (fid,*) rtmp
            lBc%gm%t(i) = rtmp
            IF (i .EQ. 1) THEN
               IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2            " should be zero in <"//TRIM(ftmp%fname)//">"
            ELSE
               rtmp = rtmp - lBc%gm%t(i-1)
               IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err = "Non-"//
     2            "increasing time trend is found in <"//
     3            TRIM(ftmp%fname)//">"
            END IF
         END DO
         lBc%gm%period = lBc%gm%t(j)
         DO b=1, msh(iM)%fa(iFa)%nNo
            READ (fid,*) Ac
            IF (Ac.GT.msh(iM)%gnNo .OR. Ac.LE.0) THEN
               err = "Entry "//Ac//" "//b//" is out of bound in "//
     2            TRIM(ftmp%fname)
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry "//Ac//" "//b//" not found in face "//
     2            TRIM(ftmp%fname)
            END IF
            DO i=1, j
               READ (fid,*) lBc%gm%d(:,a,i)
            END DO
         END DO
         CLOSE(fid)
      CASE DEFAULT
         err=TRIM(list%ping("Time dependence",lPtr))//" Unexpected type"
      END SELECT

!     To impose value or flux
      ltmp = .FALSE.
      lPtr => list%get(ltmp,"Impose flux")
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_flx)

!     To zero-out perimeter or not. Default is .true. for Dir
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) ltmp = .TRUE.
      lPtr => list%get(ltmp,"Zero out perimeter")
      lBc%bType = IBCLR(lBc%bType,bType_zp)
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_zp)

!     Impose BC on the state variable or its integral
      ltmp = .FALSE.
      IF (phys.EQ.phys_lElas .OR. phys.EQ.phys_mesh .OR.
     2    phys.EQ.phys_struct) ltmp = .TRUE.
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
      CASE ('D_dependent')
         lBc%bType = IBSET(lBc%bType,bType_ddep)
      CASE ('User_defined')
         lBc%bType = IBSET(lBc%bType,bType_ud)
         lPtr => list%get(fTmp,"Spatial profile file path",1)
         fid = fTmp%open()
         iM  = lBc%iM
         iFa = lBc%iFa
         ALLOCATE(lBc%gx(msh(iM)%fa(iFa)%nNo), ptr(msh(iM)%gnNo))
!     I am seting all the nodes to zero just in case a node is not set
         ptr    = 0
!     Preparing the pointer array
         DO a=1, msh(iM)%fa(iFa)%nNo
            lBc%gx(a) = 0D0
            Ac        = msh(iM)%fa(iFa)%gN(a)
            ptr(Ac)   = a
         END DO
         DO b=1, msh(iM)%fa(iFa)%nNo
            READ (fid,*) Ac, rtmp
            IF (Ac.GT.msh(iM)%gnNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//
     2            TRIM(fTmp%fname)
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry <"//b//" not found in face "//
     2            TRIM(fTmp%fname)
            END IF
            lBc%gx(a) = rtmp
         END DO
         CLOSE(fid)
      CASE DEFAULT
         err = TRIM(list%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

      RETURN
      END SUBROUTINE READBC
!####################################################################
!     This subroutine is to read inputs of VTK files or boundaries.
!     nDOP(1) is the total number of outputs, nDOP(2) is the default
!     number of outputs for VTK files, nDOP(3) is for boundaries, and
!     nDOP(4) is for volume
      SUBROUTINE READOUTPUTS(lEq, nDOP, outputs, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER, INTENT(IN) :: nDOP(4), outputs(nDOP(1))
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER nOut, iOut, i, j
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPO

      lEq%nOutput = nDOP(1)
      ALLOCATE(lEq%output(nDOP(1)))

      DO iOut=1, nDOP(1)
         SELECT CASE (outputs(iOut))
         CASE (out_velocity)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Velocity"
         CASE (out_pressure)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = nsd
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Pressure"
         CASE (out_acceleration)
            lEq%output(iOut)%grp  = outGrp_A
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Acceleration"
         CASE (out_temperature)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Temperature"
         CASE (out_displacement)
            lEq%output(iOut)%grp  = outGrp_D
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Displacement"
         CASE (out_WSS)
            lEq%output(iOut)%grp  = outGrp_WSS
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = maxnsd
            lEq%output(iOut)%name = "WSS"
         CASE (out_vorticity)
            lEq%output(iOut)%grp  = outGrp_vort
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = maxnsd
            lEq%output(iOut)%name = "Vorticity"
         CASE (out_energyFlux)
            lEq%output(iOut)%grp  = outGrp_eFlx
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Energy_flux"
         CASE (out_heatFlux)
            lEq%output(iOut)%grp  = outGrp_hFlx
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Heat_flux"
         CASE (out_absVelocity)
            lEq%output(iOut)%grp  = outGrp_absV
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Absolute_velocity"
         CASE (out_strainInv)
            lEq%output(iOut)%grp  = outGrp_stInv
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Strain_invariants"
         CASE (out_vortex)
            lEq%output(iOut)%grp  = outGrp_vortex
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Vortex"
         CASE DEFAULT
            err = "Internal output undefined"
         END SELECT
      END DO

!     These are the default values, we use the first nDef/nBDef outputs
      DO j=1, 3
         lEq%output(1:nDOP(j+1))%wtn(j) = .TRUE.
      END DO
!     First reading the outputs for VTK files and then for boundaries
!     and last for the volume
      nOut = list%srch("Output")
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output",iOut)
         SELECT CASE(TRIM(ctmp))
         CASE("Spatial")
            j = 1
         CASE("B_INT")
            j = 2
         CASE("V_INT")
            j = 3
         CASE DEFAULT
            j = -1
            err = TRIM(list%ping("Output",lPO))//" Undefined keyword"
         END SELECT
         DO i=1, lEq%nOutput
            lPtr => lPO%get(lEq%output(i)%wtn(j),lEq%output(i)%name)
            IF (TRIM(lEq%output(i)%name) .EQ. "Vortex") THEN
               IF (nsd .ne. maxNSD) THEN
                  lEq%output(i)%wtn(j) = .FALSE.
               END IF
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE READOUTPUTS
!####################################################################
!     This subroutine reads LS type, tolerance, ...
      SUBROUTINE READLS(ilsType, lEq, list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ilsType
      TYPE(eqType), INTENT(INOUT) :: lEq
      TYPE(listType), INTENT(INOUT) :: list

      INTEGER lSolverType, FSILSType
      LOGICAL flag
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr, lPL

!     Default LS
      lSolverType = ilsType
      FSILSType   = ilsType
      lPL => list%get(ctmp,"LS type")
      IF (ASSOCIATED(lPL)) THEN
         SELECT CASE(TRIM(ctmp))
         CASE('NS')
            lSolverType = lSolver_NS
            FSILSType   = LS_TYPE_NS
         CASE('GMRES')
            lSolverType = lSolver_GMRES
            FSILSType   = LS_TYPE_GMRES
         CASE('CG')
            lSolverType = lSolver_CG
            FSILSType   = LS_TYPE_CG
         CASE('BICG')
            lSolverType = lSolver_BICGS
            FSILSType   = LS_TYPE_BICGS
         CASE DEFAULT
            err = TRIM(list%ping("LS type",lPL))//" Undefined type"
         END SELECT
      END IF

      lEq%ls%LS_Type = lSolverType
         CALL FSILS_LS_CREATE(lEq%FSILS, FSILSType)

!     Default Preconditioners
      IF (useTrilinosLS) THEN
         lEq%ls%PREC_Type = PREC_TRILINOS_DIAGONAL
      ELSE
         lEq%ls%PREC_Type = PREC_NONE
      END IF

      IF (ASSOCIATED(lPL)) THEN
         lPtr => lPL%get(ctmp, "Preconditioner")
         IF (ASSOCIATED(lPtr)) THEN
            SELECT CASE(TRIM(ctmp))
            CASE('fsiLS', 'FSILS', 'svFSI')
               lEq%ls%PREC_Type = PREC_FSILS
            CASE('Trilinos-Diagonal')
               lEq%ls%PREC_Type = PREC_TRILINOS_DIAGONAL
               useTrilinosLS = .TRUE.
            CASE('Trilinos-BlockJacobi', 'BlockJacobi')
               lEq%ls%PREC_Type = PREC_TRILINOS_BLOCK_JACOBI
               useTrilinosLS = .TRUE.
            CASE('Trilinos-ILU')
               lEq%ls%PREC_Type = PREC_TRILINOS_ILU
               useTrilinosLS = .TRUE.
            CASE('Trilinos-ILUT')
               lEq%ls%PREC_Type = PREC_TRILINOS_ILUT
               useTrilinosLS = .TRUE.
            CASE('Trilinos-IC')
               lEq%ls%PREC_Type = PREC_TRILINOS_IC
               useTrilinosLS = .TRUE.
            CASE('Trilinos-ICT')
               lEq%ls%PREC_Type = PREC_TRILINOS_ICT
               useTrilinosLS = .TRUE.
            CASE ('Trilinos-ML')
               lEq%ls%PREC_Type = PREC_TRILINOS_ML
               useTrilinosLS = .TRUE.
            CASE DEFAULT
               err = TRIM(list%ping("Preconditioner",lPtr))
     2          //" Undefined type"
               lEq%ls%PREC_Type = PREC_NONE
            END SELECT
            print *, "Using Preconditioner: ", TRIM(ctmp)
         END IF

!        Set solver options file
         lPtr => lPL%get(lEq%ls%optionsFile,
     2            "Solver options file")
         IF (.NOT.useTrilinosAssemAndLS) THEN
            lPtr => lPL%get(flag, "Use Trilinos for assembly")
            IF (ASSOCIATED(lPtr)) useTrilinosAssemAndLS = flag
         END IF

         lPtr => lPL%get(lEq%ls%mItr,"Max iterations",ll=1)
         IF (ASSOCIATED(lPtr)) lEq%FSILS%RI%mItr = lEq%ls%mItr

         lPtr => lPL%get(lEq%ls%relTol,"Tolerance",lb=0D0,ub=1D0)
         IF (ASSOCIATED(lPtr)) lEq%FSILS%RI%relTol = lEq%ls%relTol

         lPtr => lPL%get(lEq%ls%absTol,"Absolute tolerance",
     2      lb=0D0,ub=1D0)
         IF (ASSOCIATED(lPtr)) lEq%FSILS%RI%absTol = lEq%ls%absTol

         lPtr => lPL%get(lEq%ls%sD,"Krylov space dimension",ll=1)
         IF (ASSOCIATED(lPtr)) lEq%FSILS%RI%sD = lEq%ls%sD

         IF (.NOT.useTrilinosLS) THEN
            lPtr => lPL%get(lEq%FSILS%RI%mItr,"Max iterations",ll=1)
            lPtr => lPL%get(lEq%FSILS%GM%mItr,"NS-GM max iterations",
     2         ll=1)
            lPtr => lPL%get(lEq%FSILS%CG%mItr,"NS-CG max iterations",
     2         ll=1)

            lPtr => lPL%get(lEq%FSILS%RI%relTol,"Tolerance",
     2         lb=0D0,ub=1D0)
            lPtr => lPL%get(lEq%FSILS%GM%relTol,"NS-GM tolerance",
     2         lb=0D0,ub=1D0)
            lPtr => lPL%get(lEq%FSILS%CG%relTol,"NS-CG tolerance",
     2         lb=0D0,ub=1D0)

            lPtr =>lPL%get(lEq%FSILS%RI%absTol,"Absolute tolerance",
     2         lb=0D0,ub=1D0)
            lEq%FSILS%GM%absTol = lEq%FSILS%RI%absTol
            lEq%FSILS%CG%absTol = lEq%FSILS%RI%absTol

            lEq%FSILS%GM%sD = lEq%FSILS%RI%sD
         END IF
      END IF

!     Check LS inputs
      IF (useTrilinosAssemAndLS .AND. .NOT.useTrilinosLS) err =
     2   "Error in LS inputs. Use Trilinos based LS"

      IF (useTrilinosLS) THEN
         IF (lSolverType .EQ. lSolver_NS) err =
     2   "NS solver is not supported with Trilinos or PETSc. Use FSILS"
      END IF

      IF (useTrilinosLS .AND. lEq%ls%PREC_Type.EQ.PREC_FSILS) err =
     2   "Cannot combine FSILS preconditioner with Trilinos"

      RETURN
      END SUBROUTINE READLS
!####################################################################
!     This subroutine reads properties of remesher (only for FSI)
      SUBROUTINE READRMSH(list)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list

      CHARACTER(LEN=stdL) ctmp
      INTEGER :: iM, jM
      TYPE(listType), POINTER :: lPtr, lPR, llP

      lPR => list%get(ctmp,"Remesher")
      IF (ASSOCIATED(lPR)) THEN
         SELECT CASE(TRIM(ctmp))
         CASE('Tetgen')
            rmsh%method = RMSH_TETGEN
         CASE('Meshsim')
            rmsh%method = RMSH_MESHSIM
         CASE DEFAULT
            err = TRIM(list%ping("Remesher",lPR))//" Undefined"
         END SELECT
      ELSE
         err = "Remesher properties not provided"
         RETURN
      END IF

      IF (ASSOCIATED(lPR)) THEN
         ALLOCATE(rmsh%maxEdgeSize(nMsh))
         rmsh%maxEdgeSize(:) = 1.5D-1
         rmsh%minDihedAng = 1.0D1
         rmsh%maxRadRatio = 1.15D0
         rmsh%freq = 100
         rmsh%cpVar = 10

         DO iM=1, nMsh
            lPtr => lPR%get(ctmp,"Max edge size",iM)
            IF (ASSOCIATED(lPtr)) THEN
               DO jM=1, nMsh
                  IF (ctmp .EQ. msh(jM)%name)
     2                llP => lPtr%get(rmsh%maxEdgeSize(jM),"val")
               END DO
            END IF
         END DO
         lPtr => lPR%get(rmsh%minDihedAng,"Min dihedral angle")
         lPtr => lPR%get(rmsh%maxRadRatio,"Max radius ratio",ll=1D0)
         lPtr => lPR%get(rmsh%cpVar,"Frequency for copying data")
         lPtr => lPR%get(rmsh%freq,"Remesh frequency")
      END IF

      RETURN
      END SUBROUTINE READRMSH
!####################################################################

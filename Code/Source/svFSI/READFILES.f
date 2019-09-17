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
!     and allocating required space for program.
!
!--------------------------------------------------------------------

      SUBROUTINE READFILES
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE

      LOGICAL :: flag
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
         saveVTK      = .FALSE.
         bin2VTK      = .FALSE.
         saveAve      = .FALSE.
         sepOutput    = .FALSE.
         saveATS      = 1
         saveIncr     = 10
         nITs         = 0
         startTS      = 0
         roInf        = 2D-1
         stFileName   = "stFile"
         iniFilePath  = ""
         stopTrigName = "STOP_SIM"
         rmsh%isReqd  = .FALSE.
         ichckIEN     = .TRUE.
         zeroAve      = .FALSE.
         useTrilinosLS         = .FALSE.
         useTrilinosAssemAndLS = .FALSE.
         cmmInit      = .FALSE.
         cmmVarWall   = .FALSE.
         shlEq        = .FALSE.
         pstEq        = .FALSE.
         sstEq        = .FALSE.
         ibFlag       = .FALSE.

         flag  = .FALSE.
         mfsIn = "svFSI.inp"
         INQUIRE(FILE=mfsIn, EXIST=flag)
         IF (.NOT.flag) THEN
            i = IARGC()
            IF (i .EQ. 0) THEN
               err = "Configuration file is required as an argument"
            ELSEIF (i .GT. 1) THEN
               err = "Too many arguments"
            END IF
            CALL GETARG(1,ctmp)
            mfsIn = ctmp
         END IF
      END IF

      list = listType(mfsIn,io)

      IF (.NOT.resetSim) THEN
         lPtr => list%get(ctmp, "Save results in folder")
         IF (ASSOCIATED(lPtr)) THEN
            saveName = TRIM(ctmp)
            appPath = TRIM(saveName)//delimiter
         ELSE
            saveName = ""
            appPath = STR(cm%np())//"-procs"//delimiter
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
         std = "                   (svFSI)                     "
         std = "                                               "
         std = " ----------------------------------------------"
         std = "                                               "
         std = " This program comes with ABSOLUTELY NO WARRANTY"
         std = " Compiled on "//date(5:6)//"-"//date(7:8)//
     2      "-"//date(1:4)
         std = "                                               "

         IF (cm%seq() .AND. cm%nT().EQ.1) THEN
            std = " Running sequentially"
         ELSE
            IF (.NOT.cm%seq()) std = " Running in parallel with "//
     2         cm%np()//" tasks"
            IF (cm%nT() .NE. 1) std = " Multi-threading with "//
     2         cm%nT()//" proccesors"
         END IF

         lPtr => list%get(fTmp,"Simulation initialization file path")
         IF (ASSOCIATED(lPtr)) iniFilePath = fTmp%fname

         lPtr => list%get(nsd,"Number of spatial dimensions",
     2      1,ll=2,ul=3)
         nstd = 6
         IF (nsd .EQ. 2) nstd = 3

         lPtr => list%get(nTs,"Number of time steps",1,ll=1)
         lPtr => list%get(startTS,"Starting time step",ll=0)
         lPtr => list%get(dt,"Time step size",1,lb=0D0)
         lPtr => list%get(nITs,"Number of initialization time steps",
     2      ll=0)
         lPtr => list%get(roInf,"Spectral radius of infinite time step",
     2      ll=0D0,ul=1D0)

         lPtr =>list%get(stopTrigName,
     2      "Searched file name to trigger stop")
         stopTrigName = TRIM(appPath)//stopTrigName
         lPtr => list%get(ichckIEN, "Check IEN order")

         saveName = "result"
         lPtr => list%get(saveVTK, "Save results to VTK format")
         lPtr => list%get(saveName,"Name prefix of saved VTK files")
         lPtr => list%get(saveIncr,"Increment in saving VTK files",ll=1)
         lPtr => list%get(saveATS,"Start saving after time step",ll=1)
         saveName = TRIM(appPath)//saveName

         lPtr => list%get(saveAve,"Save averaged results")
         lPtr => list%get(zeroAve,"Start averaging from zero")

         lPtr => list%get(stFileRepl,"Overwrite restart file")
         IF (.NOT.saveVTK .AND. stFileRepl) wrn = " Overwriting "//
     2      "restart files is not a good idea when not saving to VTK"
         lPtr => list%get(stFileName,"Restart file name")
         stFileName = TRIM(appPath)//stFileName

         stFileIncr = saveIncr
         lPtr => list%get(stFileIncr,
     2      "Increment in saving restart files",ll=0)
         lPtr => list%get(bin2VTK,"Convert BIN to VTK format")

         lPtr => list%get(rmsh%isReqd, "Simulation requires remeshing")
         IF (rmsh%isReqd) THEN
            saveVTK = .TRUE.
            IF (bin2VTK) err = "BIN to VTK conversion is not allowed"//
     2         " with dynamic remeshing"
         END IF
      END IF ! resetSim

!--------------------------------------------------------------------
!     Reading the mesh
      CALL READMSH(list)

!     Reading immersed boundary mesh data
      i = list%srch("Add IB")
      IF (i .GT. 0) THEN
         ibFlag = .TRUE.
         ALLOCATE(ib)
         CALL IB_READMSH(list)
         CALL IB_READOPTS(list)
      END IF

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
         IF (ibFlag) THEN
            IF (eq(iEq)%phys .EQ. phys_fluid .OR.
     2          eq(iEq)%phys .EQ. phys_FSI) THEN
               CALL IB_READEQ(eq(iEq), lPtr, ctmp)
            END IF
         END IF

         IF (eq(iEq)%phys .EQ. phys_fluid) THEN
            IF (iEq .NE. 1) err = "fluid equation must come first"
         END IF
         IF (eq(iEq)%phys .EQ. phys_CMM) THEN
            IF (iEq .NE. 1) err = "CMM equation must come first"
            IF (mvMsh) err = "Both CMM and ALE cannot be solved "//
     2         "together"
         END IF
         IF (eq(iEq)%phys .EQ. phys_FSI) THEN
            IF (iEq .NE. 1) err = "FSI equation must come first"
            lPtr => list%get(ctmp,"Add equation",2)
            IF (ctmp .NE. "mesh") err = "mesh equation"//
     2         " has to be specified after FSI equation"
         END IF
         IF (rmsh%isReqd .AND. eq(1)%phys.NE.phys_FSI) THEN
            err = "Remeshing is applicable only for FSI equation"
         END IF
         IF (iCntct) THEN
            IF (eq(iEq)%phys .NE. phys_shell) err =
     2         "Contact model is applicable for shell problems only"
            IF (nMsh .EQ. 1) err =
     2         "More than one mesh is needed to apply contact model"
         END IF
         lPtr => list%get(ctmp,"Add equation",1)
         IF (eq(iEq)%phys .EQ. phys_heatF) THEN
            IF (ctmp.NE."fluid" .AND. ctmp.NE."FSI") THEN
               err = "heatF equation has to be specified after"//
     2            " fluid/FSI equation"
            END IF
         END IF
         IF (eq(iEq)%phys .EQ. phys_mesh) THEN
            IF (.NOT.mvMsh) err = "mesh equation can only"//
     2         " be specified after FSI equation"
         END IF
      END DO

      IF (useTrilinosLS .AND. cplBC%schm.NE.cplBC_NA) THEN
         IF (cplBC%schm .NE. cplBC_E) err = "Using Trilinos with "//
     2      "semi-implicit or implicit CPLBC can lead to NaNs. "//
     3      "Try Explicit coupling or NS solver instead."
      END IF

      IF (cem%cpld) THEN
         IF (nEq .EQ. 1) err = "Min equations (2) not solved for"//
     2      " electro-mechanics coupling"
         i = 0

         DO iEq=1, nEq
            IF (eq(iEq)%phys .EQ. phys_CEP .OR.
     2          eq(iEq)%phys .EQ. phys_struct .OR.
     3          eq(iEq)%phys .EQ. phys_vms_struct) i = i + 1
         END DO
         IF (i .NE. 2) err = "Both electrophysiology and struct have"//
     2      " to be solved for electro-mechanics"

         IF (cem%aStress .AND. cem%aStrain) err = "Cannot set "//
     2      "both active strain and active stress coupling"

         IF (cem%aStrain) THEN
            IF (nsd .NE. 3) err = "Active strain coupling is allowed"//
     2         " only for 3D bodies"
            DO iEq=1, nEq
               DO i=1, eq(iEq)%nDmn
                  IF (eq(iEq)%dmn(i)%phys .NE. phys_vms_struct .AND.
     2                eq(iEq)%dmn(i)%phys .NE. phys_struct) CYCLE
                  IF (eq(iEq)%dmn(i)%stM%isoType .NE. stIso_HO) err =
     2               "Active strain is allowed with Holzapfel-Ogden "//
     3               "passive constitutive model only"
               END DO
            END DO
         END IF
      END IF

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

      INTEGER, PARAMETER :: maxOutput = 18

      INTEGER fid, iBc, iBf, iM, iFa, phys(3), propL(maxNProp,10),
     2   outPuts(maxOutput), nDOP(4)
      CHARACTER(LEN=stdL) ctmp

      TYPE(listType), POINTER :: lPtr, lPBC, lPBF
      TYPE(fileType) fTmp

      lPtr => list%get(lEq%coupled,"Coupled")
      lPtr => list%get(lEq%minItr,"Min iterations",ll=1)
      lPtr => list%get(lEq%maxItr,"Max iterations",ll=1)
      lPtr => list%get(lEq%dBr,"Residual dB reduction",ul=0D0)
      lPtr => list%get(lEq%tol,"Tolerance",ll=0D0)

!     Coupled BC stuff
!     Initializing coupled BC if neccessary
      IF (.NOT.ALLOCATED(cplBC%xo)) THEN
         lPBC => list%get(ctmp,"Couple to genBC")
         IF (ASSOCIATED(lPBC)) THEN
            cplBC%useGenBC = .TRUE.
         ELSE
            lPBC => list%get(ctmp,"Couple to cplBC")
         END IF
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
            IF (cplBC%useGenBC) THEN
               lPtr => lPBC%get(fTmp,"0D code file path",1)
               cplBC%binPath = fTmp%fname
               cplBC%commuName = "GenBC.int"
               cplBC%nX = 0
               ALLOCATE(cplBC%xo(cplBC%nX))
            ELSE
               lPtr => lPBC%get(cplBC%nX,"Number of unknowns",1,ll=0)
               ALLOCATE(cplBC%xo(cplBC%nX))
               cplBC%xo = 0D0

               lPtr => lPBC%get(fTmp,"0D code file path",1)
               cplBC%binPath = fTmp%fname

               lPtr => lPBC%get(fTmp,
     2            "Unknowns initialization file path")
               IF (ASSOCIATED(lPtr)) THEN
                  fid = fTmp%open()
                  READ (fid,*) cplBC%xo
                  CLOSE (fid)
               END IF

               lPtr => lPBC%get(cplBC%commuName,
     2            "File name for 0D-3D communication")
               cplBC%commuName = TRIM(appPath)//cplBC%commuName

               lPtr => lPBC%get(cplBC%saveName,
     2            "File name for saving unknowns")
               cplBC%saveName = TRIM(appPath)//cplBC%saveName

               lPtr => lPBC%get(cplBC%nXp,
     2            "Number of user-defined outputs")
               ALLOCATE(cplBC%xp(cplBC%nXp))
            END IF
         END IF
      END IF

      propL = prop_NA
      SELECT CASE (eqName)
!     FLUID Navier-Stokes solver ------------------------------------
      CASE ('fluid')
         lEq%phys = phys_fluid

         propL(1,1) = fluid_density
         propL(2,1) = permeability
         propL(3,1) = backflow_stab
         propL(4,1) = f_x
         propL(5,1) = f_y
         IF (nsd .EQ. 3) propL(6,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/10,2,3,0/)
         outPuts(1)  = out_velocity
         outPuts(2)  = out_pressure
         outPuts(3)  = out_energyFlux
         outPuts(4)  = out_acceleration
         outPuts(5)  = out_WSS
         outPuts(6)  = out_vorticity
         outPuts(7)  = out_strainInv
         outPuts(8)  = out_vortex
         outPuts(9)  = out_traction
         outPuts(10) = out_viscosity

         CALL READLS(lSolver_NS, lEq, list)

!     HEAT FLUID advection diffusion solver -------------------------
      CASE ('heatF', 'dyeTransport', 'scalarTransport', 'AD')
         lEq%phys = phys_heatF

         propL(1,1) = conductivity
         propL(2,1) = source_term
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/2,1,2,0/)
         outPuts(1) = out_temperature
         outPuts(2) = out_heatFlux

         CALL READLS(lSolver_GMRES, lEq, list)

!     HEAT SOLID Laplac equation solver------------------------------
      CASE ('heatS', 'laplace', 'poisson')
         lEq%phys = phys_heatS

         propL(1,1) = conductivity
         propL(2,1) = source_term
         propL(3,1) = solid_density
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

         nDOP = (/4,1,0,0/)
         outPuts(1) = out_displacement
         outPuts(2) = out_velocity
         outPuts(3) = out_acceleration
         outPuts(4) = out_integ

         CALL READLS(lSolver_CG, lEq, list)

!     STRUCTURAL with nonlinear displacement equation solver---------
      CASE ('struct')
         lEq%phys = phys_struct

         propL(1,1) = solid_density
         propL(2,1) = damping
         propL(3,1) = elasticity_modulus
         propL(4,1) = poisson_ratio
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/10,1,0,0/)
         outPuts(1)  = out_displacement
         outPuts(2)  = out_velocity
         outPuts(3)  = out_acceleration
         outPuts(4)  = out_stress
         outPuts(5)  = out_fibDir
         outPuts(6)  = out_fibAlign
         outputs(7)  = out_strainInv
         outPuts(8)  = out_integ
         outPuts(9)  = out_jacobian
         outPuts(10) = out_defGrad

         CALL READLS(lSolver_CG, lEq, list)

!     VMS-stabilized nonlinear STRUCTURAL mechanics solver ---
      CASE ('vms_struct')
         lEq%phys = phys_vms_struct
         sstEq    = .TRUE.

         propL(1,1) = solid_density
         propL(2,1) = elasticity_modulus
         propL(3,1) = poisson_ratio
         propL(4,1) = ctau_M
         propL(5,1) = ctau_C
         propL(6,1) = f_x
         propL(7,1) = f_y
         IF (nsd .EQ. 3) propL(8,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/11,1,0,0/)
         outPuts(1)  = out_displacement
         outPuts(2)  = out_velocity
         outPuts(3)  = out_pressure
         outPuts(4)  = out_acceleration
         outPuts(5)  = out_stress
         outPuts(6)  = out_fibDir
         outPuts(7)  = out_fibAlign
         outPuts(8)  = out_strainInv
         outPuts(9)  = out_integ
         outPuts(10) = out_jacobian
         outPuts(11) = out_defGrad

         CALL READLS(lSolver_CG, lEq, list)

!     PRESTRESS equation solver---------------------
      CASE ('prestress')
         lEq%phys = phys_preSt
         pstEq    = .TRUE.

         propL(1,1) = solid_density
         propL(2,1) = damping
         propL(3,1) = elasticity_modulus
         propL(4,1) = poisson_ratio
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/2,2,0,0/)
         outPuts(1) = out_displacement
         outPuts(2) = out_stress

         CALL READLS(lSolver_CG, lEq, list)

!     Nonlinear SHELL mechanics solver --------
      CASE ('shell')
         IF (nsd .NE. 3) err = "Shell mechanics can be solved only"//
     2      " in 3-dimensions"
         lEq%phys = phys_shell
         shlEq    = .TRUE.

         propL(1,1) = solid_density
         propL(2,1) = damping
         propL(3,1) = elasticity_modulus
         propL(4,1) = poisson_ratio
         propL(5,1) = shell_thickness
         propL(6,1) = f_x
         propL(7,1) = f_y
         propL(8,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/3,1,0,0/)
         outPuts(1) = out_displacement
         outPuts(2) = out_velocity
         outPuts(3) = out_integ

         CALL READLS(lSolver_CG, lEq, list)

!     COUPLED MOMENTUM FLUID STRUCTURE INTERACTION equation solver---
      CASE ('CMM')
         IF (nsd .NE. 3) err = "CMM eq. is not implemented for 2D "//
     2      "domains"
         lEq%phys = phys_CMM

         ALLOCATE(cmmBdry(gtnNo))
         cmmBdry = 0
         lPtr => list%get(ctmp, "Initialize")
         IF (ASSOCIATED(lPtr)) THEN
            cmmInit = .TRUE.
            IF (nEq .GT. 1) err = "More than one eqn. is not allowed"//
     2         "while initializing CMM"
            CALL TO_LOWER(ctmp)
            SELECT CASE (TRIM(ctmp))
            CASE ('inflate', 'inf')
               pstEq = .FALSE.
            CASE ('prestress', 'prest')
               pstEq = .TRUE.
            CASE DEFAULT
               err = " Unknown CMM initialization option"
            END SELECT

!           Set cmmBdry vector to be edge nodes of the wall
            DO iM=1, nMsh
               CALL SETCMMBDRY(msh(iM), cmmBdry)
            END DO
         END IF

         lPBC => list%get(ctmp, "Variable wall properties")
         IF (ASSOCIATED(lPBC)) THEN
            cmmVarWall = .TRUE.
            IF (.NOT.ALLOCATED(varWallProps)) THEN
               ALLOCATE(varWallProps(2,gtnNo))
               varWallProps = 0D0
            END IF

            iM  = 0
            iFa = 0
            IF (cmmInit) THEN
               CALL FINDMSH(ctmp, iM)
            ELSE
               CALL FINDFACE(ctmp, iM, iFa)
            END IF
            lPtr => lPBC%get(ctmp, "Wall properties file path", 1)
            CALL READWALLPROPSFF(ctmp, iM, iFa)
            NULLIFY(lPBC)
         END IF

         IF (.NOT.cmmInit) THEN
            propL(1,1) = fluid_density
            propL(2,1) = permeability
            propL(3,1) = backflow_stab
            propL(4,1) = solid_density
            propL(5,1) = poisson_ratio
            propL(6,1) = damping
            IF (.NOT.cmmVarWall) THEN
               propL(7,1) = shell_thickness
               propL(8,1) = elasticity_modulus
            END IF

            nDOP = (/10,2,4,0/)
            outPuts(1)  = out_velocity
            outPuts(2)  = out_pressure
            outPuts(3)  = out_energyFlux
            outPuts(4)  = out_absVelocity
            outPuts(5)  = out_acceleration
            outPuts(6)  = out_WSS
            outPuts(7)  = out_vorticity
            outPuts(8)  = out_displacement
            outPuts(9)  = out_strainInv
            outPuts(10) = out_viscosity
         ELSE
            propL(1,1) = poisson_ratio
            IF (.NOT.cmmVarWall) THEN
               propL(2,1) = shell_thickness
               propL(3,1) = elasticity_modulus
            END IF

            IF (pstEq) THEN
               nDOP = (/2,2,0,0/)
               outPuts(1) = out_displacement
               outPuts(2) = out_stress
            ELSE
               nDOP = (/1,1,0,0/)
               outPuts(1) = out_displacement
            END IF
         END IF

         CALL READDOMAIN(lEq, propL, list)
         IF (cmmInit) lEq%dmn(:)%prop(solid_density) = 0D0

         CALL READLS(lSolver_GMRES, lEq, list)

!     FLUID STRUCTURE INTERACTION equation solver---------------------
      CASE ('FSI')
         mvMsh = .TRUE.
         lEq%phys = phys_FSI

!        3 possible equations: fluid (must), struct/vms_struct
         phys(1) = phys_fluid
         phys(2) = phys_struct
         phys(3) = phys_vms_struct

!        fluid properties
         propL(1,1) = fluid_density
         propL(2,1) = permeability
         propL(3,1) = backflow_stab
         propL(4,1) = f_x
         propL(5,1) = f_y
         IF (nsd .EQ. 3) propL(6,1) = f_z

!        struct properties
         propL(1,2) = solid_density
         propL(2,2) = elasticity_modulus
         propL(3,2) = poisson_ratio
         propL(4,2) = damping
         propL(5,2) = f_x
         propL(6,2) = f_y
         IF (nsd .EQ. 3) propL(7,2) = f_z

!        vms_struct properties
         propL(1,3) = solid_density
         propL(2,3) = elasticity_modulus
         propL(3,3) = poisson_ratio
         propL(4,3) = ctau_M
         propL(5,3) = ctau_C
         propL(6,3) = f_x
         propL(7,3) = f_y
         IF (nsd .EQ. 3) propL(8,3) = f_z

         CALL READDOMAIN(lEq, propL, list, phys)

         nDOP = (/18,3,2,0/)
         outPuts(1)  = out_velocity
         outPuts(2)  = out_pressure
         outPuts(3)  = out_displacement
         outPuts(4)  = out_energyFlux
         outPuts(5)  = out_absVelocity
         outPuts(6)  = out_acceleration
         outPuts(7)  = out_WSS
         outPuts(8)  = out_vorticity
         outPuts(9)  = out_vortex
         outPuts(10) = out_traction
         outPuts(11) = out_stress
         outPuts(12) = out_fibDir
         outPuts(13) = out_fibAlign
         outPuts(14) = out_strainInv
         outPuts(15) = out_integ
         outPuts(16) = out_jacobian
         outPuts(17) = out_defGrad
         outPuts(18) = out_viscosity

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
!     Cardiac Electro-Physiology equation----------------------------
      CASE ('CEP')
         lEq%phys = phys_CEP
         cepEq    = .TRUE.

         lPtr => list%get(ctmp, "Coupling with mechanics")
         IF (ASSOCIATED(lPtr)) THEN
            cem%cpld = .TRUE.
            CALL TO_LOWER(ctmp)
            SELECT CASE (TRIM(ctmp))
            CASE ("active stress")
               cem%aStress = .TRUE.
            CASE ("active strain")
               cem%aStrain = .TRUE.
            CASE DEFAULT
               err = "Undefined coupling for cardiac electromechanics"
            END SELECT
         END IF

         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/1,1,0,0/)
         outPuts(1) = out_actionPotential

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

!     Initialize cplBC for RCR-type BC
      IF (ANY(BTEST(lEq%bc(:)%bType,bType_RCR))) THEN
         cplBC%schm = cplBC_SI
         cplBC%nX   = cplBC%nFa
         cplBC%nXp  = cplBC%nFa + 1
         IF (ALLOCATED(cplBC%xo)) err = "ERROR: cplBC structure is "//
     2      "already initialized. Unexpected behavior."
         ALLOCATE(cplBc%xo(cplBc%nX), cplBC%xp(cplBC%nXp))
         cplBC%xo = 0D0
         cplBC%xp = 0D0
         cplBC%saveName = TRIM(appPath)//"RCR.dat"
      END IF
!--------------------------------------------------------------------
!     Searching for body forces
      lEq%nBf = list%srch("Add BF")
      ALLOCATE(lEq%bf(lEq%nBf))
      DO iBf=1, lEq%nBf
         lPBF => list%get(ctmp,"Add BF",iBf)
         CALL FINDMSH(ctmp, lEq%bf(iBf)%iM)
         CALL READBF(lEq%bf(iBf), lPBF)
         IF (BTEST(lEq%bf(iBf)%bType,bType_Neu) .OR.
     2      BTEST(lEq%bf(iBf)%bType,bType_trac)) THEN
            IF (.NOT.shlEq .AND. .NOT.cmmInit) THEN
               err = "Pressure or traction load can be applied only"//
     2            " for shells or for initializing CMM"
            END IF
         END IF
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
            CASE("vms_struct")
               lEq%dmn(iDmn)%phys = phys_vms_struct
               IF (.NOT.sstEq) sstEq = .TRUE.
            CASE DEFAULT
               err = TRIM(lPD%ping("Equation",lPtr))//
     2            "Equation must be fluid/struct/vms_struct"
            END SELECT
         ELSE
            lEq%dmn(iDmn)%phys = lEq%phys
         END IF
         DO iPhys=1, nPhys
            IF (lEq%dmn(iDmn)%phys .EQ. phys(iPhys)) EXIT
         END DO
         IF (iPhys .GT. nPhys) err = "Undefined phys is used"

         DO iProp=1, maxNProp
            rtmp = 0D0
            prop = propList(iProp,iPhys)
            SELECT CASE (prop)
            CASE (prop_NA)
               EXIT
            CASE (fluid_density)
               IF(lEq%phys .EQ. phys_CMM) THEN
                  lPtr => lPD%get(rtmp,"Fluid density",1,lb=0D0)
               ELSE
                  lPtr => lPD%get(rtmp,"Density",1,lb=0D0)
               END IF
            CASE (solid_density)
               IF(lEq%phys .EQ. phys_CMM) THEN
                  lPtr => lPD%get(rtmp,"Solid density",1,ll=0D0)
               ELSE
                  lPtr => lPD%get(rtmp,"Density",ll=0D0)
               END IF
            CASE (elasticity_modulus)
               lPtr => lPD%get(rtmp,"Elasticity modulus",1,lb=0D0)
            CASE (poisson_ratio)
               lPtr => lPD%get(rtmp,"Poisson ratio",1,ll=0D0,ul=5D-1)
            CASE (conductivity)
               lPtr => lPD%get(rtmp,"Conductivity",1,ll=0D0)
            CASE (f_x)
               lPtr => lPD%get(rtmp,"Force_X")
            CASE (f_y)
               lPtr => lPD%get(rtmp,"Force_Y")
            CASE (f_z)
               lPtr => lPD%get(rtmp,"Force_Z")
            CASE (permeability)
               rtmp = 100000000000.0D0
               lPtr => lPD%get(rtmp,"Permeability",lb=0D0)
            CASE (backflow_stab)
               rtmp = 2D-1
               lPtr => lPD%get(rtmp,"Backflow stabilization coefficient"
     2            ,ll=0D0)
            CASE (source_term)
               lPtr => lPD%get(rtmp,"Source term")
            CASE (damping)
               lPtr => lPD%get(rtmp,"Mass damping")
            CASE (shell_thickness)
               lPtr => lPD%get(rtmp,"Shell thickness",1,lb=0D0)
            CASE (ctau_M)
               lPtr => lPD%get(rtmp,
     2            "Momentum stabilization coefficient")
               IF (.NOT.ASSOCIATED(lPtr)) rtmp = 0.001D0
            CASE (ctau_C)
               lPtr => lPD%get(rtmp,
     2            "Continuity stabilization coefficient")
               IF (.NOT.ASSOCIATED(lPtr)) rtmp = 0.0D0
            CASE DEFAULT
               err = "Undefined properties"
            END SELECT
            lEq%dmn(iDmn)%prop(prop) = rtmp
         END DO

         IF (lEq%dmn(iDmn)%phys .EQ. phys_CEP) THEN
            CALL READCEP(lEq%dmn(iDmn), lPD)
         END IF

         IF (lEq%dmn(iDmn)%phys.EQ.phys_struct  .OR.
     2       lEq%dmn(iDmn)%phys.EQ.phys_vms_struct .OR.
     3       lEq%dmn(iDmn)%phys.EQ.phys_preSt) THEN
            CALL READMATMODEL(lEq%dmn(iDmn), lPD)
         END IF

         IF (lEq%dmn(iDmn)%phys .EQ. phys_fluid .OR.
     2       (lEq%dmn(iDmn)%phys.EQ.phys_CMM .AND. .NOT.cmmInit)) THEN
            CALL READVISCMODEL(lEq%dmn(iDmn), lPD)
         END IF
      END DO

      RETURN
      END SUBROUTINE READDOMAIN
!--------------------------------------------------------------------
      END SUBROUTINE READEQ
!--------------------------------------------------------------------
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
         CALL TO_LOWER(ctmp)
         SELECT CASE(TRIM(ctmp))
         CASE('ns')
            lSolverType = lSolver_NS
            FSILSType   = LS_TYPE_NS
         CASE('gmres')
            lSolverType = lSolver_GMRES
            FSILSType   = LS_TYPE_GMRES
         CASE('cg')
            lSolverType = lSolver_CG
            FSILSType   = LS_TYPE_CG
         CASE('bicg')
            lSolverType = lSolver_BICGS
            FSILSType   = LS_TYPE_BICGS
         CASE DEFAULT
            err = TRIM(list%ping("LS type",lPL))//" Undefined type"
         END SELECT
      END IF

      lEq%ls%LS_Type = lSolverType
      CALL FSILS_LS_CREATE(lEq%FSILS, FSILSType)

!     Default Preconditioners
      lEq%ls%PREC_Type = PREC_FSILS
#ifdef WITH_TRILINOS
      IF (FSILSType .EQ. LS_TYPE_NS) THEN
         lEq%ls%PREC_Type = PREC_FSILS
      ELSE
         useTrilinosLS = .TRUE.
         lEq%ls%PREC_Type = PREC_TRILINOS_DIAGONAL
      END IF
#endif

      IF (ASSOCIATED(lPL)) THEN
         lPtr => lPL%get(ctmp, "Preconditioner")
         IF (ASSOCIATED(lPtr)) THEN
            CALL TO_LOWER(ctmp)
            SELECT CASE(TRIM(ctmp))
            CASE('fsils', 'svfsi')
               lEq%ls%PREC_Type = PREC_FSILS
               useTrilinosLS = .FALSE.
#ifdef WITH_TRILINOS
            CASE('trilinos-diagonal')
               lEq%ls%PREC_Type = PREC_TRILINOS_DIAGONAL
               useTrilinosLS = .TRUE.
            CASE('trilinos-blockjacobi', 'blockjacobi')
               lEq%ls%PREC_Type = PREC_TRILINOS_BLOCK_JACOBI
               useTrilinosLS = .TRUE.
            CASE('trilinos-ilu')
               lEq%ls%PREC_Type = PREC_TRILINOS_ILU
               useTrilinosLS = .TRUE.
            CASE('trilinos-ilut')
               lEq%ls%PREC_Type = PREC_TRILINOS_ILUT
               useTrilinosLS = .TRUE.
            CASE('trilinos-ic')
               lEq%ls%PREC_Type = PREC_TRILINOS_IC
               useTrilinosLS = .TRUE.
            CASE('trilinos-ict')
               lEq%ls%PREC_Type = PREC_TRILINOS_ICT
               useTrilinosLS = .TRUE.
            CASE ('trilinos-ml')
               lEq%ls%PREC_Type = PREC_TRILINOS_ML
               useTrilinosLS = .TRUE.
#endif
            CASE DEFAULT
               err = TRIM(list%ping("Preconditioner",lPtr))
     2          //" Undefined type"
            END SELECT
            std = " Using preconditioner: "//TRIM(ctmp)
         ELSE
            SELECT CASE (lEq%ls%PREC_Type)
            CASE (PREC_FSILS)
               std = " Using preconditioner: FSILS"
            CASE (PREC_TRILINOS_DIAGONAL)
               std = " Using preconditioner: Trilinos-Diagonal"
            CASE DEFAULT
               err = " Undefined preconditioner"
            END SELECT
         END IF

         IF (.NOT.useTrilinosAssemAndLS) THEN
            lPtr => lPL%get(flag, "Use Trilinos for assembly")
            IF (ASSOCIATED(lPtr)) useTrilinosAssemAndLS = flag
            IF (useTrilinosAssemAndLS .AND. ibFlag) err =
     2         "Cannnot assemble immersed boundaries using Trilinos."//
     3         " Use Trilinos for linear solver only."
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
     2   "NS solver is not supported with Trilinos"//
     3   "Use GMRES or BICG instead."
      END IF

      RETURN
      END SUBROUTINE READLS
!--------------------------------------------------------------------
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
      CHARACTER(LEN=stdL) ctmp, stmp
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
         CASE (out_integ)
            lEq%output(iOut)%grp  = outGrp_I
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            IF (nsd .EQ. 2) THEN
               lEq%output(iOut)%name = "Area"
            ELSE
               lEq%output(iOut)%name = "Volume"
            END IF
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
         CASE (out_strainInv)
            lEq%output(iOut)%grp  = outGrp_stInv
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Strain_invariants"
         CASE (out_absVelocity)
            lEq%output(iOut)%grp  = outGrp_absV
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Absolute_velocity"
         CASE (out_vortex)
            lEq%output(iOut)%grp  = outGrp_vortex
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Vortex"
         CASE (out_traction)
            lEq%output(iOut)%grp  = outGrp_trac
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Traction"
         CASE (out_stress)
            lEq%output(iOut)%grp  = outGrp_stress
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nstd
            lEq%output(iOut)%name = "Stress"
         CASE (out_fibDir)
            lEq%output(iOut)%grp  = outGrp_fN
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd
            lEq%output(iOut)%name = "Fiber_direction"
         CASE (out_fibAlign)
            lEq%output(iOut)%grp  = outGrp_fA
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Fiber_alignment"
         CASE (out_actionPotential)
            lEq%output(iOut)%grp  = outGrp_Y
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Action_potential"
         CASE (out_jacobian)
            lEq%output(iOut)%grp  = outGrp_J
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Jacobian"
         CASE (out_defGrad)
            lEq%output(iOut)%grp  = outGrp_F
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = nsd*nsd
            lEq%output(iOut)%name = "Deformation_gradient"
         CASE (out_viscosity)
            lEq%output(iOut)%grp  = outGrp_Visc
            lEq%output(iOut)%o    = 0
            lEq%output(iOut)%l    = 1
            lEq%output(iOut)%name = "Viscosity"
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
         CASE("B_INT", "Boundary_integral")
            j = 2
         CASE("V_INT", "Volume_integral")
            j = 3
         CASE ("Alias")
            CYCLE
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

!     Read any alias names for outputs
      DO iOut=1, nOut
         lPO => list%get(ctmp,"Output",iOut)
         SELECT CASE (TRIM(ctmp))
         CASE ("Alias")
            DO i=1, lEq%nOutput
               lPtr => lPO%get(stmp, TRIM(lEq%output(i)%name))
               IF (ASSOCIATED(lPtr)) lEq%output(i)%name = TRIM(stmp)
            END DO
            EXIT
         CASE DEFAULT
            CYCLE
         END SELECT
      END DO

      RETURN
      END SUBROUTINE READOUTPUTS
!--------------------------------------------------------------------
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
      REAL(KIND=8) rtmp, RCR(3)
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
      CASE ("Traction","Trac")
         lBc%bType = IBSET(lBc%bType,bType_trac)

         lPtr => list%get(fTmp, "Traction values file path (vtp)")
         IF (ASSOCIATED(lPtr)) THEN
            iM  = lBc%iM
            iFa = lBc%iFa
            i = nsd
            j = 2
            a = msh(iM)%fa(iFa)%nNo
            ALLOCATE(lBc%gm)
            lBc%gm%nTP = j
            lBc%gm%dof = i
            ALLOCATE(lBc%gm%t(j), lBc%gm%d(i,a,j))
            lBc%gm%t(1) = 0D0
            lBc%gm%t(2) = HUGE(rtmp)

            CALL READTRACBCFF(lBc%gm, msh(iM)%fa(iFa), fTmp%fname)
            lBc%bType = IBSET(lBc%bType,bType_gen)
            lBc%bType = IBSET(lBc%bType,bType_flat)

            lPtr => list%get(rtmp, "Traction multiplier",1)
            IF (.NOT.ASSOCIATED(lPtr)) rtmp = 1D0
            lBc%gm%d(:,:,:) = lBc%gm%d(:,:,:) * rtmp

            ALLOCATE(lBc%eDrn(nsd), lBc%h(nsd))
            lBc%eDrn    = 0
            lBc%h       = 0D0
            lBc%weakDir = .FALSE.
            RETURN
         END IF
      CASE ("Robin", "Rbn")
         lBc%bType = IBSET(lBc%bType,bType_Robin)
         lBc%bType = IBSET(lBc%bType,bType_Neu)
      CASE ("Coupled Momentum","CMM")
         lBc%bType = IBSET(lBc%bType,bType_CMM)
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected BC type"
      END SELECT

!     Allocate traction
      ALLOCATE(lBc%h(nsd))
      lBc%h = 0D0

      ALLOCATE(lBc%eDrn(nsd))
      lBc%eDrn = 0
      lPtr => list%get(lBc%eDrn,"Effective direction")

      ctmp = "Steady"
      lPtr => list%get(ctmp,"Time dependence")
      SELECT CASE (ctmp)
      CASE ('Steady')
         lBc%bType = IBSET(lBc%bType,bType_std)
         IF (BTEST(lBc%bType, bType_trac)) THEN
            lBc%h = 0D0
            lPtr => list%get(lBc%h,"Value")
         ELSE
            lBc%g = 0D0
            lPtr => list%get(lBc%g,"Value")
         END IF
      CASE ('Unsteady')
         lBc%bType = IBSET(lBc%bType,bType_ustd)
         IF (BTEST(lBc%bType,bType_trac)) err =
     2      "For unsteady traction, use general time dependence"
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
            ALLOCATE(lBc%gt%r(lBc%gt%n))
            ALLOCATE(lBc%gt%i(lBc%gt%n))
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
         IF (BTEST(lBc%bType,bType_Robin)) err = "Resistance"//
     2      " is not defined for Robin BC"
         IF (phys.NE.phys_fluid .AND.
     2       phys.NE.phys_FSI   .AND.
     3       phys.NE.phys_CMM) THEN
            err = "Resistance is only defined for fluid/FSI equations"
         END IF

         lPtr => list%get(lBc%r,"Value",1)
      CASE ('RCR', 'Windkessel')
         lBc%bType = IBSET(lBc%bType, bType_RCR)
         IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "RCR BC is only"//
     2      " defined for Neu BC"
         IF (phys.NE.phys_fluid .AND.
     2       phys.NE.phys_FSI   .AND.
     3       phys.NE.phys_CMM) THEN
            err = "RCR BC is only defined for fluid/CMM/FSI equations"
         END IF

         lPtr => list%get(RCR,"RCR values",1)
         lBc%RCR%Rp = RCR(1)
         lBc%RCR%C  = RCR(2)
         lBc%RCR%Rd = RCR(3)
         lPtr => list%get(lBc%RCR%Pd, "Distal pressure")

         IF (cplBC%schm.NE.cplBC_NA .OR. ALLOCATED(cplBC%xo)) err =
     2      "RCR cannot be used in conjunction with cplBC."
         cplBC%nFa = cplBC%nFa + 1
         lBc%cplBcPtr = cplBC%nFa
      CASE ('General')
         iM  = lBc%iM
         iFa = lBc%iFa
         lBc%bType = IBSET(lBc%bType,bType_gen)
         lPtr => list%get(ftmp,"BCT file path (vtp)")
         IF (ASSOCIATED(lPtr)) THEN
            ALLOCATE(lBc%gm)
            CALL READBCT(lBc%gm, msh(iM)%fa(iFa), ftmp%fname)
         ELSE
            lPtr =>list%get(fTmp,
     2         "Temporal and spatial values file path")
            IF (.NOT.ASSOCIATED(lPtr)) err = "Input file for General"//
     2         " BC (bct.vtp) is not provided"
            fid = fTmp%open()
            READ (fid,*) i, j, a
            IF (a .NE. msh(iM)%fa(iFa)%nNo) THEN
               err = "Number of nodes does not match between "//
     2            TRIM(msh(iM)%fa(iFa)%name)//" and "//TRIM(fTmp%fname)
            END IF
            IF (i.LT.1 .OR. i.GT.nsd) err = "0 < dof <= "//nsd//
     2         " is violated in "//fTmp%fname

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
               Ac = msh(iM)%lN(Ac)
               IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2            "detected for BC. Mesh: "//TRIM(msh(iM)%name)//
     3            ", Face: "//TRIM(msh(iM)%fa(iFa)%name)//", Node: "//
     4            STR(a)//" gN: "//STR(msh(iM)%fa(iFa)%gN(a))
               ptr(Ac) = a
            END DO
            DO i=1, j
               READ (fid,*) rtmp
               lBc%gm%t(i) = rtmp
               IF (i .EQ. 1) THEN
                  IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2               " should be zero in <"//TRIM(ctmp)//">"
               ELSE
                  rtmp = rtmp - lBc%gm%t(i-1)
                  IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err = "Non-incre"//
     2               "asing time trend is found in <"//TRIM(ctmp)//">"
               END IF
            END DO
            lBc%gm%period = lBc%gm%t(j)
            DO b=1, msh(iM)%fa(iFa)%nNo
               READ (fid,*) Ac
               IF (Ac.GT.msh(iM)%gnNo .OR. Ac.LE.0) THEN
                  err = "Entry "//b//" is out of bound in "//ctmp
               END IF
               a = ptr(Ac)
               IF (a .EQ. 0) THEN
                  err = "Entry "//b//" not found in face "//ctmp
               END IF
               DO i=1, j
                  READ (fid,*) lBc%gm%d(:,a,i)
               END DO
            END DO
            CLOSE(fid)
         END IF
      CASE DEFAULT
         err=TRIM(list%ping("Time dependence",lPtr))//" Unexpected type"
      END SELECT

!     Stiffness and damping parameters for Robin BC
      IF (BTEST(lBc%bType, bType_Robin)) THEN
         lPtr => list%get(lBc%k, "Stiffness", 1)
         lPtr => list%get(lBc%c, "Damping", 1)
      END IF

!     To impose value or flux
      ltmp = .FALSE.
      lPtr => list%get(ltmp,"Impose flux")
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_flx)

!     To zero-out perimeter or not. Default is .true. for Dir/CMM
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) ltmp = .TRUE.
      lPtr => list%get(ltmp, "Zero out perimeter")
      lBc%bType = IBCLR(lBc%bType,bType_zp)
      IF (ltmp .OR. BTEST(lBc%bType,bType_CMM))
     2   lBc%bType = IBSET(lBc%bType,bType_zp)

!     Impose BC on the state variable or its integral
      ltmp = .FALSE.
      IF (phys.EQ.phys_lElas .OR. phys.EQ.phys_mesh .OR.
     2    phys.EQ.phys_struct .OR. phys.EQ.phys_shell) ltmp = .TRUE.
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
         iM  = lBc%iM
         iFa = lBc%iFa
         ALLOCATE(lBc%gx(msh(iM)%fa(iFa)%nNo), ptr(msh(iM)%gnNo))
!     I am seting all the nodes to zero just in case a node is not set
         ptr    = 0
!     Preparing the pointer array
         DO a=1, msh(iM)%fa(iFa)%nNo
            lBc%gx(a) = 0D0
            Ac        = msh(iM)%fa(iFa)%gN(a)
            Ac = msh(iM)%lN(Ac)
            IF (Ac .EQ. 0) err = "Incorrect global node number "//
     2         "detected for BC. Mesh: "//TRIM(msh(iM)%name)//
     3         ", Face: "//TRIM(msh(iM)%fa(iFa)%name)//
     4         ", Node: "//STR(a)//" gN: "//STR(msh(iM)%fa(iFa)%gN(a))
            ptr(Ac)   = a
         END DO
         DO b=1, msh(iM)%fa(iFa)%nNo
            READ (fid,*) Ac, rtmp
            IF (Ac.GT.msh(iM)%gnNo .OR. Ac.LE.0) THEN
               err = "Entry "//b//" is out of bound in "//fTmp%fname
            END IF
            a = ptr(Ac)
            IF (a .EQ. 0) THEN
               err = "Entry <"//b//" not found in face "//fTmp%fname
            END IF
            lBc%gx(a) = rtmp
         END DO
         CLOSE(fid)
      CASE DEFAULT
         err = TRIM(list%ping("Profile",lPtr))//" Unexpected profile"
      END SELECT

!     Weak Dirichlet BC for fluid/FSI equations
      lBc%weakDir = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (phys.EQ.phys_fluid .OR. phys.EQ.phys_FSI)
     2      lPtr => list%get(lBc%weakDir, "Weakly applied")
      END IF

!     Read penalty values for weakly applied Dir BC
      IF (lBc%weakDir) THEN
         lBc%bType = IBCLR(lBc%bType,bType_zp)
         lPtr => list%get(rtmp, "Penalty parameter")
         IF (ASSOCIATED(lPtr)) lBc%tauB = rtmp
         lPtr => list%get(rtmp, "Penalty parameter (tangential)")
         IF (ASSOCIATED(lPtr)) lBc%tauB(1) = rtmp
         lPtr => list%get(rtmp, "Penalty parameter (normal)")
         IF (ASSOCIATED(lPtr)) lBc%tauB(2) = rtmp
      END IF

!     Read BCs for shells with triangular elements. Not necessary for
!     NURBS elements
      lPtr => list%get(ctmp,"Shell BC type")
      IF (ASSOCIATED(lPtr)) THEN
         SELECT CASE (ctmp)
         CASE ("Fixed", "fixed", "Clamped", "clamped")
            lBc%bType = IBSET(lBc%bType,bType_fix)
            IF (.NOT.BTEST(lBc%bType,bType_Dir)) err = "Fixed BC "//
     2         "can be applied for Dirichlet boundaries only"
         CASE ("Hinged", "hinged")
            lBc%bType = IBSET(lBc%bType,bType_hing)
            IF (.NOT.BTEST(lBc%bType,bType_Dir)) err = "Hinged BC "//
     2         "can be applied for Dirichlet boundaries only"
         CASE ("Free", "free")
            lBc%bType = IBSET(lBc%bType,bType_free)
            IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "Free BC "//
     2         "can be applied for Neumann boundaries only"
!        Symm BC needs to be verified
c         CASE ("Symm", "symm", "Symmetric", "symmetric")
c            lBc%bType = IBSET(lBc%bType,bType_symm)
c            IF (.NOT.BTEST(lBc%bType,bType_Neu)) err = "Symm BC "//
c     2         "can be applied for Neumann boundaries only"
         CASE DEFAULT
            err = TRIM(list%ping("Shell BC type",lPtr))//
     2         " Unexpected Shell BC type"
         END SELECT
      END IF

!     If a Neumann BC face is undeforming
      lBc%masN = 0
      IF (BTEST(lBc%bType,bType_Neu)) THEN
         ltmp = .FALSE.
         lBc%bType = IBCLR(lBc%bType,bType_undefNeu)
         lPtr => list%get(ltmp, "Undeforming Neu face")
         IF (ltmp) THEN
            IF (phys .NE. phys_vms_struct) err = "Undeforming Neu "//
     2         "face is currently formulated for VMS_STRUCT only"

            IF (BTEST(lBc%bType,bType_cpl) .OR.
     2          BTEST(lBc%bType,bType_res)) err = "Undeforming Neu "//
     2         "BC cannot be used with resistance or couple BC yet"

!           Clear any BC profile
            lBc%bType = IBCLR(lBc%bType,bType_flat)
            lBc%bType = IBCLR(lBc%bType,bType_para)
            IF (BTEST(lBc%bType,bType_ud) .OR.
     2          BTEST(lBc%bType,bType_gen)) THEN
               err = "General BC or user defined spatial profile "//
     2            "cannot be imposed on an undeforming Neu face"
            END IF

!           Clear zero perimeter flag
            lBc%bType = IBCLR(lBc%bType,bType_zp)

!           Reset profile to flat and set undeforming Neumann BC flag
            lBc%bType = IBSET(lBc%bType,bType_flat)
            lBc%bType = IBSET(lBc%bType,bType_undefNeu)

!           Set master-slave node parameters. Set a master node that
!           is not part of any other face (not on perimeter)
            iM  = lBc%iM
            iFa = lBc%iFa
            IF (ALLOCATED(ptr)) DEALLOCATE(ptr)
            ALLOCATE(ptr(gtnNo))
            ptr = 0
            DO a=1, msh(iM)%fa(iFa)%nNo
               Ac = msh(iM)%fa(iFa)%gN(a)
               ptr(Ac) = 1
            END DO

            DO j=1, msh(iM)%nFa
               IF (j .EQ. iFa) CYCLE
               DO a=1, msh(iM)%fa(j)%nNo
                  Ac = msh(iM)%fa(j)%gN(a)
                  ptr(Ac) = 0
               END DO
            END DO

            DO a=1, msh(iM)%fa(iFa)%nNo
               Ac = msh(iM)%fa(iFa)%gN(a)
               IF (ptr(Ac) .EQ. 1) THEN
                  lBc%masN = Ac
                  EXIT
               END IF
            END DO
            DEALLOCATE(ptr)
         END IF
      END IF

!     For CMM BC, load wall displacements
      IF (BTEST(lBC%bType, bType_CMM)) THEN
         lPtr => list%get(cTmp, "Initial displacements file path (vtp)")
         IF (.NOT.ASSOCIATED(lPtr) .AND. .NOT.ALLOCATED(Dinit)) THEN
            lPtr => list%get(cTmp, "Prestress file path (vtp)")
            IF (.NOT.ASSOCIATED(lPtr) .AND. .NOT.ALLOCATED(pS0)) THEN
               err = "Either wall displacement field or prestress "//
     2           "is required for CMM eqn"
            END IF

!           Read prestress tensor here
            IF (ASSOCIATED(lPtr)) THEN
               iM  = lBc%iM
               iFa = lBc%iFa
               IF (.NOT.ALLOCATED(msh(iM)%fa(iFa)%x)) THEN
                  ALLOCATE(msh(iM)%fa(iFa)%x(nstd,msh(iM)%fa(iFa)%nNo))
                  msh(iM)%fa(iFa)%x = 0D0
               END IF
               CALL READVTPPDATA(msh(iM)%fa(iFa), cTmp, "Stress",nstd,1)
               IF (.NOT.ALLOCATED(pS0)) THEN
                  ALLOCATE(pS0(nstd,gtnNo))
                  pS0 = 0D0
               END IF
               DO a=1, msh(iM)%fa(iFa)%nNo
                  Ac = msh(iM)%fa(iFa)%gN(a)
                  Ac = msh(iM)%gN(Ac)
                  pS0(:,Ac) = msh(iM)%fa(iFa)%x(:,a)
               END DO
               DEALLOCATE(msh(iM)%fa(iFa)%x)
               NULLIFY(lPtr)
            END IF
         END IF

!        Read displacement field here
         IF (ASSOCIATED(lPtr)) THEN
            iM  = lBc%iM
            iFa = lBc%iFa
            IF (.NOT.ALLOCATED(msh(iM)%fa(iFa)%x)) THEN
               ALLOCATE(msh(iM)%fa(iFa)%x(nsd,msh(iM)%fa(iFa)%nNo))
               msh(iM)%fa(iFa)%x = 0D0
            END IF
            CALL READVTPPDATA(msh(iM)%fa(iFa), cTmp, "Displacement",
     2         nsd, 1)
            IF (.NOT.ALLOCATED(Dinit)) THEN
               ALLOCATE(Dinit(nsd,gtnNo))
               Dinit = 0D0
            END IF
            DO a=1, msh(iM)%fa(iFa)%nNo
               Ac = msh(iM)%fa(iFa)%gN(a)
               Ac = msh(iM)%gN(Ac)
               Dinit(:,Ac) = msh(iM)%fa(iFa)%x(:,a)
            END DO
            DEALLOCATE(msh(iM)%fa(iFa)%x)
         END IF

!        Set cmmBdry vector for wall nodes
         iM  = lBc%iM
         iFa = lBc%iFa
         DO a=1, msh(iM)%fa(iFa)%nNo
            Ac = msh(iM)%fa(iFa)%gN(a)
            Ac = msh(iM)%gN(Ac)
            cmmBdry(Ac) = 1
         END DO
      END IF

      RETURN
      END SUBROUTINE READBC
!--------------------------------------------------------------------
!     This routine reads a body force
      SUBROUTINE READBF(lBf, list)
      USE COMMOD
      USE ALLFUN
      USE LISTMOD
      IMPLICIT NONE

      TYPE(bfType), INTENT(INOUT) :: lBf
      TYPE(listType), INTENT(INOUT) :: list

      LOGICAL flag
      INTEGER iM, a, i, j, Ac, fid
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp
      TYPE(listType), POINTER :: lPtr
      TYPE(fileType) fTmp

      iM  = lBf%iM
      lBf%dof = nsd

!     Reading the type: Vol/Neu/Trac
      lPtr => list%get(ctmp,"Type")
      CALL TO_LOWER(ctmp)
      SELECT CASE (ctmp)
      CASE ("volumetric","vol","internal","int")
         lBf%bType = IBSET(lBf%bType,bfType_vol)
      CASE ("traction","trac")
         lBf%bType = IBSET(lBf%bType,bfType_trac)
      CASE ("neumann","neu","pressure")
         lBf%bType = IBSET(lBf%bType,bfType_Neu)
         lBf%dof   = 1
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected BF type"
      END SELECT

!     Time dependence
      ctmp = "steady"
      lPtr => list%get(ctmp,"Time dependence")
      CALL TO_LOWER(ctmp)
      SELECT CASE (ctmp)
      CASE ('steady')
         lBf%bType = IBSET(lBf%bType,bfType_std)
         ALLOCATE(lBf%b(lBf%dof))
         lBf%b = 0D0
         IF (lBf%dof .EQ. 1) THEN
            lPtr => list%get(lBf%b(1),"Value",1)
         ELSE
            lPtr => list%get(lBf%b,"Value",1)
         END IF

      CASE ('unsteady')
         lBf%bType = IBSET(lBf%bType,bfType_ustd)

         ALLOCATE(lBf%bt(lBf%dof))
         flag = .FALSE.
         lPtr => list%get(flag,"Ramp function")
         DO a=1, lBf%dof
            lBf%bt(a)%lrmp = flag
            lPtr => list%get(ftmp, "Temporal values file path", a)
            IF (ASSOCIATED(lPtr)) THEN
               fid = fTmp%open()
               READ(fid,*) i, j
               IF (i .LT. 2) THEN
                  std = "Enter nPnts nFCoef; nPts*(t Q)"
                  err = "Wrong format in: "//fTmp%fname
               END IF
               lBf%bt(a)%n = j
               ALLOCATE(lBf%bt(a)%r(j))
               ALLOCATE(lBf%bt(a)%i(j))
               CALL FFT(fid, i, lBf%bt(1))
               CLOSE(fid)
            ELSE
               lPtr => list%get(fTmp,"Fourier coefficients file path",a)
               fid = fTmp%open()
               READ (fid,*) lBf%bt(a)%ti
               READ (fid,*) lBf%bt(a)%T
               READ (fid,*) lBf%bt(a)%qi
               READ (fid,*) lBf%bt(a)%qs
               READ (fid,*) j
               lBf%bt(a)%n = j
               ALLOCATE(lBf%bt(a)%r(j))
               ALLOCATE(lBf%bt(a)%i(j))
               DO i=1, j
                  READ (fid,*) lBf%bt(a)%r(i), lBf%bt(a)%i(i)
               END DO
               CLOSE(fid)
            END IF
         END DO

      CASE ('spatial')
         lBf%bType = IBSET(lBf%bType,bfType_spl)
         ALLOCATE(lBf%bx(lBf%dof,gtnNo))
         lBf%bx = 0D0
         lPtr => list%get(cTmp,"Spatial values file path")

         ALLOCATE(msh(iM)%x(lBf%dof,msh(iM)%gnNo))
         msh(iM)%x = 0D0
         IF (BTEST(lBf%bType,bfType_vol)) THEN
            CALL READVTUPDATA(msh(iM), ctmp, "Body_force", lBf%dof, 1)
         ELSE IF (BTEST(lBf%bType,bfType_trac)) THEN
            CALL READVTUPDATA(msh(iM), ctmp, "Traction", lBf%dof, 1)
         ELSE IF (BTEST(lBf%bType,bfType_Neu)) THEN
            CALL READVTUPDATA(msh(iM), ctmp, "Pressure", lBf%dof, 1)
         END  IF

         DO a=1, msh(iM)%gnNo
            Ac = msh(iM)%gN(a)
            lBf%bx(:,Ac) = msh(iM)%x(:,a)
         END DO
         DEALLOCATE(msh(iM)%x)

      CASE ('general')
         lBf%bType = IBSET(lBf%bType,bfType_gen)

         lPtr =>list%get(fTmp,"Temporal and spatial values file path",1)
         fid = fTmp%open()

         READ (fid,*) i, j, a
         IF (a .GT. msh(iM)%gnNo) err = "No. of nodes out of bounds"//
     2      " (body force for <"//TRIM(msh(iM)%name)//">)"
         IF (i .NE. lBf%dof) THEN
            err = "Mismatch in DOF for general type body force on <"
     2            //TRIM(msh(iM)%name)//">"
         END IF

         ALLOCATE(lBf%bm)
         lBf%bm%dof = lBf%dof
         lBf%bm%nTP = j

         ALLOCATE(lBf%bm%t(j), lBf%bm%d(lBf%dof,gtnNo,j))
         lBf%bm%t = 0D0
         lBf%bm%d = 0D0

         DO j=1, lBf%bm%nTP
            READ(fid,*) rtmp
            lBf%bm%t(j) = rtmp
            IF (j .EQ. 1) THEN
               IF (.NOT.ISZERO(rtmp)) err = "First time step should "//
     2            "be 0 (body force for <"//TRIM(msh(iM)%name)//">)"
            ELSE
               rtmp = rtmp - lBf%bm%t(j-1)
               IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err =
     2            "Non-increasing time trend found (body force for <"//
     3            TRIM(msh(iM)%name)//">)"
            END IF
         END DO
         lBf%bm%period = lBf%bm%t(lBf%bm%nTP)

         DO a=1, msh(iM)%gnNo
            READ (fid,*) Ac
            IF (Ac .GT. msh(iM)%gnNo) err =
     2       "Entry "//Ac//" is out of bound in "//TRIM(fTmp%fname)
            Ac = msh(iM)%gN(Ac)
            DO j=1, lBf%bm%nTP
               READ(fid,*) (lBf%bm%d(i,Ac,j), i=1, lBf%bm%dof)
            END DO
         END DO
         CLOSE(fid)

      CASE DEFAULT
         err = TRIM(list%ping("Time dependence",lPtr))//
     2      " Unexpected type (body force)"
      END SELECT

      RETURN
      END SUBROUTINE READBF
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
!#######################################################################
!     Reads properties of cardiac-electrophysiology model
      SUBROUTINE READCEP(lDmn, lPD)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(dmnType), INTENT(INOUT) :: lDmn
      TYPE(listType), INTENT(INOUT) :: lPD

      TYPE(listType), POINTER :: lPtr, list

      INTEGER i
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp

      lPtr => lPD%get(ctmp, "Electrophysiology model")
      lDmn%cep%nG = 0
      CALL TO_LOWER(ctmp)
      SELECT CASE(TRIM(ctmp))
      CASE ("ap", "aliev-panfilov")
         lDmn%cep%cepType = cepModel_AP
         lDmn%cep%nX = 2
         IF (cem%aStrain) err = " Active strain is not formulated "//
     2      "Aliev-Panfilov model"

      CASE ("bo", "bueno-orovio")
         lDmn%cep%cepType = cepModel_BO
         lDmn%cep%nX = 4

      CASE ("fn", "fitzhugh-nagumo")
         lDmn%cep%cepType = cepModel_FN
         lDmn%cep%nX = 2
         IF (cem%cpld) err = " Electromechanics is not formulated "//
     2      "Fitzhugh-Nagumo model"

      CASE ("ttp", "tentusscher-panfilov")
         lDmn%cep%cepType = cepModel_TTP
         lDmn%cep%nX = 7
         lDmn%cep%nG = 12

      CASE DEFAULT
         err = "Unknown electrophysiology model"
      END SELECT

      IF (nXion .LT. lDmn%cep%nX + lDmn%cep%nG) THEN
         nXion = lDmn%cep%nX + lDmn%cep%nG
      END IF

      lPtr => lPD%get(lDmn%cep%Diso,"Conductivity (iso)",ll=0D0)
      i = lPD%srch("Conductivity (ani)")
      lDmn%cep%nFn = i
      IF (i .EQ. 0) lDmn%cep%nFn = 1

      ALLOCATE(lDmn%cep%Dani(lDmn%cep%nFn))
      lDmn%cep%Dani = 0D0
      IF (i .NE. 0) THEN
         DO i=1, lDmn%cep%nFn
            lPtr => lPD%get(lDmn%cep%Dani(i), "Conductivity (ani)", i)
         END DO
      END IF

      lDmn%cep%imyo = 1
      lPtr => lPD%get(ctmp, "Myocardial zone")
      IF (ASSOCIATED(lPtr)) THEN
         CALL TO_LOWER(ctmp)
         SELECT CASE (TRIM(ctmp))
         CASE ("epi", "epicardium")
            lDmn%cep%imyo = 1

         CASE ("endo", "endocardium", "pfib", "purkinje")
            lDmn%cep%imyo = 2

         CASE ("myo", "mid-myo", "myocardium")
            lDmn%cep%imyo = 3

         CASE DEFAULT
            err = "Undefined myocardium zone"
         END SELECT
      END IF

      lDmn%cep%Istim%A  = 0.0d0
      lDmn%cep%Istim%Ts = 9999999.0d0
      lDmn%cep%Istim%CL = 9999999.0d0
      lDmn%cep%Istim%Td = 0.0d0

      list => lPD%get(ctmp, "Stimulus")
      IF (ASSOCIATED(list)) THEN
         lPtr => list%get(lDmn%cep%Istim%A, "Amplitude")
         IF (.NOT.ISZERO(lDmn%cep%Istim%A)) THEN
            lPtr => list%get(lDmn%cep%Istim%Ts, "Start time")
            lPtr => list%get(lDmn%cep%Istim%Td, "Duration")
            lPtr => list%get(lDmn%cep%Istim%CL, "Cycle length")
            IF (.NOT.ASSOCIATED(lPtr)) THEN
               lDmn%cep%Istim%CL = REAL(nTS, KIND=8) * dt
            END IF
         END IF
      END IF

!     Dual time stepping for cellular activation model
      lPtr => lPD%get(rtmp, "Time step for integration")
      IF (ASSOCIATED(lPtr)) THEN
         lDmn%cep%dt = rtmp
      ELSE
         lDmn%cep%dt = dt
      END IF

      lDmn%cep%odes%tIntType = tIntType_FE
      lPtr => lPD%get(ctmp, "ODE solver")
      IF (ASSOCIATED(lPtr)) THEN
         CALL TO_LOWER(ctmp)
         SELECT CASE (ctmp)
         CASE ("fe", "euler", "explicit")
            lDmn%cep%odes%tIntType = tIntType_FE

         CASE ("rk", "rk4", "runge")
            lDmn%cep%odes%tIntType = tIntType_RK4

         CASE ("cn", "cn2", "implicit")
            lDmn%cep%odes%tIntType = tIntType_CN2
            IF (lDmn%cep%cepType .EQ. cepModel_TTP) THEN
               err = "Implicit time integration for tenTusscher-"//
     2            "Panfilov model can give unexpected results. "//
     3            "Use FE or RK4 instead"
            END IF

         CASE DEFAULT
            err = " Unknown ODE time integrator"
         END SELECT
      END IF

      IF (lDmn%cep%odes%tIntType .EQ. tIntType_CN2) THEN
         list => lPtr
         lDmn%cep%odes%maxItr = 5
         lDmn%cep%odes%absTol = 1D-8
         lDmn%cep%odes%relTol = 1D-4
         lPtr => list%get(lDmn%cep%odes%maxItr, "Maximum iterations")
         lPtr => list%get(lDmn%cep%odes%absTol, "Absolute tolerance")
         lPtr => list%get(lDmn%cep%odes%relTol, "Relative tolerance")
      END IF

      lDmn%cep%Ksac = 0D0
      lPtr => lPD%get(rtmp, "Feedback parameter for "//
     2   "stretch-activated-currents")
      IF (ASSOCIATED(lPtr)) lDmn%cep%Ksac = rtmp
      IF (.NOT.cem%cpld) lDmn%cep%Ksac = 0D0

      RETURN
      END SUBROUTINE READCEP
!####################################################################
!     This subroutine reads properties of structure constitutive model
      SUBROUTINE READMATMODEL(lDmn, lPD)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(dmnType), INTENT(INOUT) :: lDmn
      TYPE(listType), INTENT(INOUT) :: lPD

      TYPE(listType), POINTER :: lPtr, lSt
      LOGICAL incompFlag
      CHARACTER(LEN=stdL) ctmp
      REAL(KIND=8) :: E, nu, lam, mu, kap, rtmp

!     Domain properties: elasticity modulus, poisson ratio
      E   = lDmn%prop(elasticity_modulus)
      nu  = lDmn%prop(poisson_ratio)

!     Shear modulus
      mu  = E/(1D0+nu)/2D0

!     Incompressible material
      incompFlag = .FALSE.
      IF (ISZERO(nu-0.5D0)) incompFlag = .TRUE.

!     Bulk modulus for compressible case
      IF (.NOT.incompFlag) THEN
         kap = E/(1D0-2D0*nu)/3D0
         lam = E*nu/(1D0+nu)/(1D0-2D0*nu)
      END IF

      lSt => lPD%get(ctmp, "Constitutive model")

!     Default: NeoHookean model
      IF (.NOT.ASSOCIATED(lSt)) THEN
         lDmn%stM%isoType = stIso_nHook
         lDmn%stM%C10 = mu/2.0D0
         RETURN
      END IF

      SELECT CASE (TRIM(ctmp))
      CASE ("lin", "linear")
         lDmn%stM%isoType = stIso_lin
         lDmn%stM%C10 = mu
         RETURN

      CASE ("stVK", "stVenantKirchhoff")
         IF (incompFlag) err = "Cannot choose stVK model for Poisson "//
     2      "ratio 0.5"
         lDmn%stM%isoType = stIso_stVK
         lDmn%stM%C10 = lam
         lDmn%stM%C01 = mu
         RETURN

      CASE ("m-stVK", "modified-stVK",  "modified-stVenantKirchhoff")
         IF (incompFlag) err = "Cannot choose stVK model for Poisson "//
     2      "ratio 0.5"
         lDmn%stM%isoType = stIso_mStVK
         lDmn%stM%C10 = kap
         lDmn%stM%C01 = mu
         RETURN

      CASE ("nHK", "nHK91", "neoHookean", "neoHookeanSimo91")
         lDmn%stM%isoType = stIso_nHook
         lDmn%stM%C10 = mu/2D0

      CASE ("MR", "Mooney-Rivlin")
         lDmn%stM%isoType = stIso_MR
         lPtr => lSt%get(lDmn%stM%C10, "c1")
         lPtr => lSt%get(lDmn%stM%C01, "c2")

      CASE ("HGO")
      ! Neo-Hookean ground matrix + quad penalty + anistropic fibers !
         lDmn%stM%isoType = stIso_HGO
         lDmn%stM%C10 = mu/2D0
         lPtr => lSt%get(lDmn%stM%aff, "a4")
         lPtr => lSt%get(lDmn%stM%bff, "b4")
         lPtr => lSt%get(lDmn%stM%ass, "a6")
         lPtr => lSt%get(lDmn%stM%bss, "b6")
         lPtr => lSt%get(lDmn%stM%kap, "kappa")

      CASE ("Guccione", "Gucci")
         lDmn%stM%isoType = stIso_Gucci
         lPtr => lSt%get(lDmn%stM%C10, "C")
         lPtr => lSt%get(lDmn%stM%bff, "bf")
         lPtr => lSt%get(lDmn%stM%bss, "bt")
         lPtr => lSt%get(lDmn%stM%bfs, "bfs")
         IF (nsd .NE. 3) THEN
            err = "Guccione material model is used for 3D problems "//
     2         "with 2 family of directions"
         END IF

      CASE ("HO", "Holzapfel")
      ! Holzapefel and Ogden model for myocardium !
         lDmn%stM%isoType = stIso_HO
         lPtr => lSt%get(lDmn%stM%aff, "a")
         lPtr => lSt%get(lDmn%stM%bff, "b")
         lPtr => lSt%get(lDmn%stM%aff, "a4f")
         lPtr => lSt%get(lDmn%stM%bff, "b4f")
         lPtr => lSt%get(lDmn%stM%ass, "a4s")
         lPtr => lSt%get(lDmn%stM%bss, "b4s")
         lPtr => lSt%get(lDmn%stM%afs, "afs")
         lPtr => lSt%get(lDmn%stM%bfs, "bfs")

      CASE DEFAULT
         err = "Undefined constitutive model used"
      END SELECT

      IF (incompFlag) RETURN

!     Look for dilational penalty model. HGO uses quadratic penalty model
      lPtr => lPD%get(ctmp, "Dilational penalty model")
      IF (.NOT.ASSOCIATED(lPtr)) wrn =
     2   "Couldn't find any penalty model"
      SELECT CASE(TRIM(ctmp))
      CASE ("quad", "Quad", "quadratic", "Quadratic")
         lDmn%stM%volType = stVol_Quad

      CASE ("ST91", "Simo-Taylor91")
         lDmn%stM%volType = stVol_ST91

      CASE ("M94", "Miehe94")
         lDmn%stM%volType = stVol_M94

      CASE DEFAULT
         err = "Undefined dilational penalty model"
      END SELECT

!     Default penalty parameter is equal to bulk modulus
      lDmn%stM%Kpen = kap
      lPtr => lPD%get(rtmp, "Penalty parameter")
      IF (ASSOCIATED(lPtr)) lDmn%stM%Kpen = rtmp

      RETURN
      END SUBROUTINE READMATMODEL
!####################################################################
!     This subroutine reads parameters of non-Newtonian viscosity model
      SUBROUTINE READVISCMODEL(lDmn, lPD)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE

      TYPE(dmnType), INTENT(INOUT) :: lDmn
      TYPE(listType), INTENT(INOUT) :: lPD

      TYPE(listType), POINTER :: lPtr, lVis
      CHARACTER(LEN=stdL) ctmp
      REAL(KIND=8) :: rtmp

      lVis => lPD%get(ctmp,"Viscosity",1)

      CALL TO_LOWER(ctmp)
      SELECT CASE (TRIM(ctmp))
      CASE ("constant", "const", "newtonian")
         lDmn%visc%viscType = viscType_Const
         lPtr => lVis%get(lDmn%visc%mu_i,"Value",1,lb=0D0)
         RETURN

      CASE ("carreau-yasuda", "cy")
         lDmn%visc%viscType = viscType_CY
         lPtr => lVis%get(lDmn%visc%mu_i,
     2      "Limiting high shear-rate viscosity",1,lb=0D0)
         lPtr => lVis%get(lDmn%visc%mu_o,
     2      "Limiting low shear-rate viscosity",1,lb=0D0)
         lPtr => lVis%get(lDmn%visc%lam,
     2      "Shear-rate tensor multiplier (lamda)",1,lb=0D0)
         lPtr => lVis%get(lDmn%visc%a,
     2      "Shear-rate tensor exponent (a)",1,lb=0D0)
         lPtr => lVis%get(lDmn%visc%n,"Power-law index (n)",1,lb=0D0)
         IF (lDmn%visc%mu_i .GT. lDmn%visc%mu_o) THEN
            err = "Unexpected inputs for Carreau-Yasuda model. "//
     2         "High shear-rate viscosity value should be higher than"//
     3         " low shear-rate value"
         END IF
         RETURN

      CASE ("cassons", "cass")
         lDmn%visc%viscType = viscType_Cass
         lPtr => lVis%get(lDmn%visc%mu_i,
     2      "Asymptotic viscosity parameter",1,lb=0D0)
         lPtr => lVis%get(lDmn%visc%mu_o,
     2      "Yield stress parameter",1,lb=0D0)
         lDmn%visc%lam = 0.5D0
         lPtr => lVis%get(rtmp,"Low shear-rate threshold")
         IF (ASSOCIATED(lPtr)) lDmn%visc%lam = rtmp
         RETURN

      CASE DEFAULT
         err = "Undefined constitutive model for viscosity used"
      END SELECT

      RETURN
      END SUBROUTINE READVISCMODEL
!####################################################################
!     This subroutine reads general velocity data from bct.vtp
      SUBROUTINE READBCT(lMB, lFa, fName)
      USE COMMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE

      TYPE(MBType), INTENT(INOUT) :: lMB
      TYPE(faceType), INTENT(IN) :: lFa
      CHARACTER(LEN=stdL), INTENT(IN) :: fName

      CHARACTER(LEN=*), PARAMETER :: shdr = "velocity_"
      TYPE(vtkXMLType) :: vtp
      INTEGER :: i, j, a, Ac, n, nNo, ntime, iM, istat
      REAL(KIND=8) :: t
      CHARACTER(LEN=stdL) :: stmp

      INTEGER, ALLOCATABLE :: ptr(:), gN(:)
      REAL(KIND=8), ALLOCATABLE :: tmpR(:,:)
      CHARACTER(LEN=stdL), ALLOCATABLE :: namesL(:)

      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtp, fName, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (init)"

      CALL getVTK_numPoints(vtp, nNo, iStat)
      IF (nNo .NE. lFa%nNo) err = "Mismatch in num points for face <"//
     2   TRIM(lFa%name)//">"

!     Get all the point data starting with "velocity_"
      CALL getVTK_numPointData(vtp, n, istat)
      ALLOCATE(namesL(n))
      namesL(:) = ""
      CALL getVTK_pointDataNames(vtp, namesL, istat)
      IF (istat .LT. 0) err = "VTP file read error (point names)"

      ntime = 0
      DO i=1, n
         stmp = namesL(i)
         DO j=1, 9
            IF (shdr(j:j) .NE. stmp(j:j)) EXIT
         END DO
         IF (j .LT. 9) THEN
            namesL(i) = ""
            CYCLE
         END IF
         ntime = ntime + 1
      END DO

!     Initialize lMB data structure
      lMB%dof = nsd
      lMB%nTP = ntime
      iM = lFa%iM
      ALLOCATE(lMB%t(ntime), lMB%d(nsd,nNo,ntime), ptr(msh(iM)%gnNo))
      lMB%t = 0D0
      lMB%d = 0D0
      ptr   = 0

!     Load time values from point data variable names
      ntime = 0
      DO i=1, n
         j = LEN(TRIM(namesL(i)))
         IF (j .EQ. 0) CYCLE
         stmp = namesL(i)
         stmp = stmp(10:j)
         ntime = ntime + 1
         READ(stmp,*) t
         lMB%t(ntime) = t
         IF (ntime .EQ. 1) THEN
            IF (.NOT.ISZERO(t)) err = "First time step should be zero"//
     2         " in <bct.vtp>"
         ELSE
            t = t - lMB%t(ntime-1)
            IF (ISZERO(t) .OR. t.LT.0D0) err = "Non-increasing time "//
     2         "trend is found in <bct.vtp>"
         END IF
      END DO
      lMB%period = lMB%t(ntime)

!     Prepare pointer array
      DO a=1, lFa%nNo
         Ac = lFa%gN(a)
         Ac = msh(iM)%lN(Ac)
         IF (Ac .EQ. 0) err = "Incorrect global node number detected "//
     2      "for BC. Mesh: "//TRIM(msh(iM)%name)//", Face: "//
     3      TRIM(lFa%name)//", Node: "//STR(a)//" gN: "//STR(lFa%gN(a))
         ptr(Ac) = a
      END DO

!     Get GlobalNodeID from vtp file and make sure it is consistent
!     with mesh structure
      ALLOCATE(gN(nNo))
      gN = 0
      CALL getVTK_pointData(vtp, "GlobalNodeID", gN, istat)
      DO a=1, nNo
         Ac = gN(a)
         IF (Ac.GT.msh(iM)%gnNo .OR. Ac.LE.0) THEN
            err = "Entry "//a//" is out of bound in <bct.vtp>"
         END IF
         Ac = ptr(Ac)
         IF (Ac .EQ. 0) THEN
            err = "Entry "//a//" not found in face <bct.vtp>"
         END IF
      END DO

!     Load spatial data for each time point from vtp file
      ALLOCATE(tmpR(nsd,nNo))
      ntime = 0
      DO i=1, n
         j = LEN(TRIM(namesL(i)))
         IF (j .EQ. 0) CYCLE
         ntime = ntime + 1
         tmpR  = 0D0
         CALL getVTK_pointData(vtp, namesL(i), tmpR, istat)
         DO a=1, nNo
            Ac = gN(a)
            Ac = ptr(Ac)
            lMB%d(:,Ac,ntime) = tmpR(:,a)
         END DO
      END DO

      CALL flushVTK(vtp)
      DEALLOCATE(namesL, gN, tmpR)

      RETURN
      END SUBROUTINE READBCT
!####################################################################
!     This subroutine reads traction data from a vtp file and stores
!     in moving BC data structure
      SUBROUTINE READTRACBCFF(lMB, lFa, fName)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      USE vtkXMLMod
      IMPLICIT NONE

      TYPE(MBType), INTENT(INOUT) :: lMB
      TYPE(faceType), INTENT(INOUT) :: lFa
      CHARACTER(LEN=stdL) :: fName

      TYPE(vtkXMLType) :: vtp
      TYPE(faceType) :: gFa
      INTEGER :: iStat, a, e, Ac

      INTEGER, ALLOCATABLE :: ptr(:), lIEN(:,:)
      REAL(KIND=8), ALLOCATABLE :: tmpX(:,:)

!     Read Traction data from VTP file
      iStat = 0
      std = " <VTK XML Parser> Loading file <"//TRIM(fName)//">"
      CALL loadVTK(vtp, fName, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (init)"

      CALL getVTK_numPoints(vtp, gFa%nNo, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (num points)"

      CALL getVTK_numElems(vtp, gFa%nEl, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (num cells)"

      CALL getVTK_nodesPerElem(vtp, gFa%eNoN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (nodes per cell)"

      ALLOCATE(gFa%x(nsd,gFa%nNo), tmpX(maxNSD,gFa%nNo))
      CALL getVTK_pointCoords(vtp, tmpX, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (coords)"
      gFa%x(:,:) = tmpX(1:nsd,:)

      ALLOCATE(gFa%IEN(gFa%eNoN,gFa%nEl))
      CALL getVTK_elemIEN(vtp, gFa%IEN, iStat)
      IF (iStat .LT. 0) err = "VTP file read error (ien)"
      gFa%IEN = gFa%IEN + 1

!     Use gFa%N as temporary array to load traction data
      ALLOCATE(gFa%N(nsd,gFa%nNo))
      CALL getVTK_pointData(vtp, "NS_Traction", tmpX, iStat)
      gFa%N(:,:) = tmpX(1:nsd,:)
      IF (iStat .LT. 0) err = "VTP file read error (point data)"

      DEALLOCATE(tmpX)
      CALL flushVTK(vtp)

!     Project traction from gFa to lFa. First prepare lFa%x, lFa%IEN
      ALLOCATE(ptr(gtnNo), lFa%x(nsd,lFa%nNo), lIEN(lFa%eNoN,lFa%nEl))

      ptr = 0
      DO a=1, lFa%nNo
         Ac = lFa%gN(a)
         ptr(Ac) = a
         lFa%x(:,a) = x(:,Ac)
      END DO

      lIEN = lFa%IEN
      DO e=1, lFa%nEl
         DO a=1, lFa%eNoN
            Ac = lFa%IEN(a,e)
            lFa%IEN(a,e) = ptr(Ac)
         END DO
      END DO
      DEALLOCATE(ptr)

      ALLOCATE(ptr(lFa%nNo))
      ptr = 0
      CALL FACEMATCH(lFa, gFa, ptr)
      lFa%IEN = lIEN

!     Copy traction data to MB data structure
      DO a=1, lFa%nNo
         Ac = ptr(a)
         lMB%d(:,a,1) = gFa%N(:,Ac)
         lMB%d(:,a,2) = gFa%N(:,Ac)
      END DO

      CALL DESTROY(gFa)
      DEALLOCATE(lFa%x, lIEN, ptr)

      RETURN
      END SUBROUTINE READTRACBCFF
!--------------------------------------------------------------------
      SUBROUTINE FACEMATCH(lFa, gFa, ptr)
      USE COMMOD
      USE ALLFUN
      USE UTILMOD
      IMPLICIT NONE

      TYPE(faceType), INTENT(IN) :: lFa, gFa
      INTEGER, INTENT(INOUT) :: ptr(lFa%nNo)

      INTEGER :: a, b, i
      REAL(KIND=8) :: ds
      TYPE(adjType) :: lAdj, gAdj
      TYPE(queueType) :: lnQ, gnQ

      INTEGER, ALLOCATABLE :: seed(:)
      LOGICAL, ALLOCATABLE :: flag(:)

      CALL SETNADJ(lFa, lAdj)

      CALL SETNADJ(gFa, gAdj)

      ALLOCATE(flag(lFa%nNo), seed(lFa%nNo))
      flag = .FALSE.
      seed = 0
      DO a=1, lFa%nNo
         DO b=1, gFa%nNo
            ds = SQRT( SUM( (lFa%x(:,a) - gFa%x(:,b))**2 ) )
            IF (ds .LT. 1D3*eps) THEN
               ptr(a)  = b
               flag(a) = .TRUE.
               seed(a) = b
               DO i=lAdj%prow(a), lAdj%prow(a+1)-1
                  CALL ENQUEUE(lnQ, lAdj%pcol(i))
                  seed(lAdj%pcol(i)) = b
               END DO
               EXIT
            END IF
         END DO
         IF (flag(a)) EXIT
      END DO
      IF (lnQ%n .EQ. 0) err = " Failed to find a matching node"

      DO WHILE (DEQUEUE(lnQ, a))
         IF (ALL(flag(:))) EXIT
         IF (flag(a)) CYCLE

         CALL DESTROY(gnQ)
         b = seed(a)
         DO i=gAdj%prow(b), gAdj%prow(b+1)-1
            CALL ENQUEUE(gnQ, gAdj%pcol(i))
         END DO

         DO WHILE (DEQUEUE(gnQ, b))
            ds = SQRT( SUM( (lFa%x(:,a) - gFa%x(:,b))**2 ) )
            IF (ds .LT. 1D3*eps) THEN
               ptr(a)  = b
               flag(a) = .TRUE.
               IF (seed(a) .EQ. 0) seed(a) = b
               DO i=lAdj%prow(a), lAdj%prow(a+1)-1
                  CALL ENQUEUE(lnQ, lAdj%pcol(i))
                  IF (seed(lAdj%pcol(i)) .EQ. 0) THEN
                     seed(lAdj%pcol(i)) = b
                  END IF
               END DO
               EXIT
            END IF
         END DO

         IF (.NOT.flag(a)) THEN
            DO b=1, gFa%nNo
               ds = SQRT( SUM( (lFa%x(:,a) - gFa%x(:,b))**2 ) )
               IF (ds .LT. 1D3*eps) THEN
                  ptr(a)  = b
                  flag(a) = .TRUE.
                  seed(a) = b
                  DO i=lAdj%prow(a), lAdj%prow(a+1)-1
                     CALL ENQUEUE(lnQ, lAdj%pcol(i))
                     seed(lAdj%pcol(i)) = b
                  END DO
                  EXIT
               END IF
            END DO
            IF (.NOT.flag(a)) err = " Failed to map node "//STR(a)
         END IF
      END DO

      DEALLOCATE(flag, seed, lAdj%prow, lAdj%pcol, gAdj%prow, gAdj%pcol)
      RETURN
      CONTAINS
!--------------------------------------------------------------------
         SUBROUTINE SETNADJ(pFa, pAdj)
         USE COMMOD
         IMPLICIT NONE

         TYPE(faceType), INTENT(IN) :: pFa
         TYPE(adjType), INTENT(INOUT) :: pAdj

         INTEGER :: a, b, e, i, Ac, Bc, maxN
         LOGICAL :: flag

         INTEGER, ALLOCATABLE :: incNd(:), adjL(:,:)

         ALLOCATE(incNd(pFa%nNo))
         incNd = 0
         DO e=1, pFa%nEl
            DO a=1, pFa%eNoN
               Ac = pFa%IEN(a,e)
               DO b=1, pFa%eNoN
                  IF (b .Eq. a) CYCLE
                  Bc = pFa%IEN(b,e)
                  incNd(Ac) = incNd(Ac) + 1
               END DO
            END DO
         END DO
         maxN = MAXVAL(incNd)

         ALLOCATE(adjL(maxN, pFa%nNo))
         adjL = 0
         incNd = 0
         DO e=1, pFa%nEl
            DO a=1, pFa%eNoN
               Ac = pFa%IEN(a,e)
               DO b=1, pFa%eNoN
                  IF (a .EQ. b) CYCLE
                  Bc = pFa%IEN(b,e)
                  flag = .TRUE.
                  DO i=1, incNd(Ac)
                     IF (adjL(i,Ac) .EQ. Bc) THEN
                        flag = .FALSE.
                        EXIT
                     END IF
                  END DO
                  IF (flag) THEN
                     incNd(Ac) = incNd(Ac) + 1
                     adjL(incNd(Ac),Ac) = Bc
                  END IF
               END DO
            END DO
         END DO
         maxN = MAXVAL(incNd)

         pAdj%nnz = 0
         DO a=1, pFa%nNo
            DO i=1, maxN
               Ac = adjL(i,a)
               IF (Ac .EQ. 0) EXIT
               pAdj%nnz = pAdj%nnz + 1
            END DO
         END DO

         ALLOCATE(pAdj%prow(pFa%nNo+1), pAdj%pcol(pAdj%nnz))
         b = 0
         pAdj%prow(1) = 1
         DO a=1, pFa%nNo
            DO i=1, maxN
               Ac = adjL(i,a)
               IF (Ac .EQ. 0) EXIT
               b = b + 1
               pAdj%pcol(b) = Ac
            END DO
            pAdj%prow(a+1) = b + 1
         END DO

         DEALLOCATE(adjL, incNd)
         RETURN
         END SUBROUTINE SETNADJ
!--------------------------------------------------------------------
      END SUBROUTINE FACEMATCH
!####################################################################
!     Find boundary edge nodes for CMM initialization
      SUBROUTINE SETCMMBDRY(lM, bNds)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(INOUT) :: lM
      INTEGER, INTENT(INOUT) :: bNds(gtnNo)

      INTEGER i, j, a, a1, b, b1, e, e1, Ac, Ac1, Bc, Bc1, nAdj
      LOGICAL flag

      INTEGER, ALLOCATABLE :: incL(:), adjL(:,:), tmpI(:,:)

!--------------------------------------------------------------------
!     First, get mesh adjacency
      ALLOCATE(incL(lM%gnNo))
      incL = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            incL(Ac) = incL(Ac) + 1
         END DO
      END DO

      nAdj = MAXVAL(incL)
      ALLOCATE(tmpI(nAdj, lM%gnNo))
      incL = 0
      tmpI = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            incL(Ac) = incL(Ac) + 1
            tmpI(incL(Ac), Ac) = e
         END DO
      END DO
      b = 2*nAdj

 001  b = b + nAdj
      DEALLOCATE(incL)
      ALLOCATE(incL(lM%gnEl))
      IF (ALLOCATED(adjL)) DEALLOCATE(adjL)
      ALLOCATE(adjL(b, lM%gnEl))
      adjL = 0
      incL = 0
      DO e=1, lM%gnEl
         DO a=1, lM%eNoN
            Ac = lM%gIEN(a,e)
            DO i=1, nAdj
               IF (tmpI(i,Ac) .EQ. 0) EXIT
               flag = .TRUE.
               DO j=1, incL(e)
                  IF (adjL(j,e) .EQ. tmpI(i,Ac)) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO
               IF (flag) THEN
                  incL(e) = incL(e) + 1
                  IF (incL(e) .GE. b) GOTO 001
                  adjL(incL(e),e) = tmpI(i,Ac)
               END IF
            END DO
         END DO
      END DO
      nAdj = MAXVAL(incL)
      DEALLOCATE(tmpI, incL)

      lM%eAdj%nnz = 0
      DO e=1, lM%gnEl
         DO i=1, nAdj
            IF (adjL(i,e) .NE. 0) THEN
               lM%eAdj%nnz = lM%eAdj%nnz + 1
            ELSE
               EXIT
            END IF
         END DO
      END DO

      ALLOCATE(lM%eAdj%prow(lM%gnEl+1), lM%eAdj%pcol(lM%eAdj%nnz))
      j = 0
      lM%eAdj%prow(1) = j + 1
      DO e=1, lM%gnEl
         DO i=1, nAdj
            IF (adjL(i,e) .NE. 0) THEN
               j = j + 1
               lM%eAdj%pcol(j) = adjL(i,e)
            ELSE
               EXIT
            END IF
         END DO
         lM%eAdj%prow(e+1) = j + 1
      END DO
      DEALLOCATE(adjL)

!--------------------------------------------------------------------
!     Loop over all elements to find edge nodes
      DO e=1, lM%gnEl
!        Create a list of neighboring elements
         nAdj = lM%eAdj%prow(e+1)-lM%eAdj%prow(e)
         ALLOCATE(incL(nAdj))
         i = 0
         DO a=lM%eAdj%prow(e), lM%eAdj%prow(e+1)-1
            i = i + 1
            incL(i) = lM%eAdj%pcol(a)
         END DO

!        Select an edge pair
         DO a=1, lM%eNoN
            b = a + 1
            IF (a .EQ. lM%eNoN) b = 1
            Ac = lM%gIEN(a,e)
            Bc = lM%gIEN(b,e)
            IF (Ac .GT. Bc) CALL SWAP(Ac, Bc)

!        Loop over all neighbors and check if this edge is found
            flag = .true.
            DO i=1, nAdj
               e1 = incL(i)
               IF (e1 .EQ. e) CYCLE
               DO a1=1, lM%eNoN
                  b1 = a1 + 1
                  IF (a1 .EQ. lM%eNoN) b1 = 1
                  Ac1 = lM%gIEN(a1,e1)
                  Bc1 = lM%gIEN(b1,e1)
                  IF (Ac1 .GT. Bc1) CALL SWAP(Ac1, Bc1)
                  IF (Ac .EQ. Ac1 .AND. Bc .EQ. Bc1) THEN
                     flag = .FALSE.
                     EXIT
                  END IF
               END DO ! b
               IF (.NOT.flag) EXIT
            END DO ! i
            IF (flag) THEN
               Ac = lM%gN(Ac)
               Bc = lM%gN(Bc)
               bNds(Ac) = 1
               bNds(Bc) = 1
            END IF
         END DO ! a
         DEALLOCATE(incL)
      END DO

      CALL DESTROY(lM%eAdj)

      RETURN
      END SUBROUTINE SETCMMBDRY
!####################################################################
!     Read CMM variable wall properties from file
      SUBROUTINE READWALLPROPSFF(fname, iM, iFa)
      USE COMMOD
      IMPLICIT NONE

      CHARACTER(LEN=stdL), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: iM, iFa

      INTEGER a, Ac

      IF (cmmInit) THEN
         IF (ALLOCATED(msh(iM)%x)) DEALLOCATE(msh(iM)%x)
         ALLOCATE(msh(iM)%x(1,msh(iM)%gnNo))
!        Read thickness
         msh(iM)%x = 0D0
         CALL READVTUPDATA(msh(iM), fname, "Thickness", 1, 1)
         DO a=1, msh(iM)%gnNo
            Ac = msh(iM)%gN(a)
            varWallProps(1,Ac) = msh(iM)%x(1,a)
         END DO

!        Read elasticity modulus
         msh(iM)%x = 0D0
         CALL READVTUPDATA(msh(iM), fname, "Elasticity_modulus", 1, 1)
         DO a=1, msh(iM)%gnNo
            Ac = msh(iM)%gN(a)
            varWallProps(2,Ac) = msh(iM)%x(1,a)
         END DO
         DEALLOCATE(msh(iM)%x)
      ELSE
         IF (ALLOCATED(msh(iM)%fa(iFa)%x)) DEALLOCATE(msh(iM)%fa(iFa)%x)
         ALLOCATE(msh(iM)%fa(iFa)%x(1,msh(iM)%fa(iFa)%nNo))
!        Read thickness
         msh(iM)%fa(iFa)%x = 0D0
         CALL READVTPPDATA(msh(iM)%fa(iFa), fname, "Thickness", 1, 1)
         DO a=1, msh(iM)%fa(iFa)%nNo
            Ac = msh(iM)%fa(iFa)%gN(a)
            Ac = msh(iM)%gN(Ac)
            varWallProps(1,Ac) = msh(iM)%fa(iFa)%x(1,a)
         END DO

!        Read elasticity modulus
         msh(iM)%fa(iFa)%x = 0D0
         CALL READVTPPDATA(msh(iM)%fa(iFa), fname, "Elasticity_modulus",
     2      1, 1)
         DO a=1, msh(iM)%fa(iFa)%nNo
            Ac = msh(iM)%fa(iFa)%gN(a)
            Ac = msh(iM)%gN(Ac)
            varWallProps(2,Ac) = msh(iM)%fa(iFa)%x(1,a)
         END DO
         DEALLOCATE(msh(iM)%fa(iFa)%x)
      END IF

      RETURN
      END SUBROUTINE READWALLPROPSFF
!####################################################################

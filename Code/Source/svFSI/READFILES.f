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
         ibFlag       = .FALSE.
         useTrilinosLS         = .FALSE.
         useTrilinosAssemAndLS = .FALSE.
         cplEM        = .FALSE.

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
         lPtr => list%get(dt,"Time step size",1,lb=0D0)
         lPtr => list%get(nITs,"Number of initialization time steps",
     2      ll=0)
         lPtr => list%get(roInf,"Spectral radius of infinite time step",
     2      ll=0D0,ul=1D0)

         lPtr =>list%get(stopTrigName,
     2      "Searched file name to trigger stop")
         stopTrigName = TRIM(appPath)//stopTrigName
         lPtr => list%get(ichckIEN, "Check IEN order")

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
      IF (cplEM) THEN
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

      INTEGER, PARAMETER :: maxOutput = 17
      INTEGER fid, iBc, phys(3), propL(maxNProp,10),
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

         nDOP = (/9,2,3,0/)
         outPuts(1) = out_velocity
         outPuts(2) = out_pressure
         outPuts(3) = out_energyFlux
         outPuts(4) = out_acceleration
         outPuts(5) = out_WSS
         outPuts(6) = out_vorticity
         outPuts(7) = out_strainInv
         outPuts(8) = out_vortex
         outPuts(9) = out_traction

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
         lEq%phys = phys_CMM

         propL(1,1) = fluid_density
         propL(2,1) = viscosity
         propL(3,1) = permeability
         propL(4,1) = backflow_stab
         propL(5,1) = solid_density
         propL(6,1) = elasticity_modulus
         propL(7,1) = poisson_ratio
         propL(8,1) = damping
         propL(9,1) = shell_thickness
         propL(10,1) = initialization_pressure
         propL(11,1) = f_x
         propL(12,1) = f_y
         IF (nsd .EQ. 3) propL(13,1) = f_z
         CALL READDOMAIN(lEq, propL, list)

         nDOP = (/9,2,4,0/)
         outPuts(1) = out_velocity
         outPuts(2) = out_pressure
         outPuts(3) = out_energyFlux
         outPuts(4) = out_absVelocity
         outPuts(5) = out_acceleration
         outPuts(6) = out_WSS
         outPuts(7) = out_vorticity
         outPuts(8) = out_displacement
         outPuts(9) = out_strainInv

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
         propL(2,1) = viscosity
         propL(3,1) = permeability
         propL(4,1) = backflow_stab
         propL(5,1) = f_x
         propL(6,1) = f_y
         IF (nsd .EQ. 3) propL(7,1) = f_z

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

         nDOP = (/17,3,2,0/)
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

!     Cardiac Electro-Physiology equation----------------------------
      CASE ('CEP')
         lEq%phys = phys_CEP
         lPtr => list%get(cplEM, "Coupled to mechanics")
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
      TYPE(listType), POINTER :: lPtr, lPD, lSFp
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
                  lPtr => lPD%get(rtmp,"Fluid Density",1,lb=0D0)
               ELSE
                  lPtr => lPD%get(rtmp,"Density",1,lb=0D0)
               END IF
            CASE (solid_density)
               IF(lEq%phys .EQ. phys_CMM) THEN
                  lPtr => lPD%get(rtmp,"Solid Density",1,lb=0D0)
               ELSE
                  lPtr => lPD%get(rtmp,"Density",1,lb=0D0)
               END IF
            CASE (viscosity)
               lPtr => lPD%get(rtmp,"Viscosity",1,lb=0D0)
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
            CASE (initialization_pressure)
               lPtr => lPD%get(rtmp,"Initialization Pressure",1,lb=0D0)
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

         IF (lEq%dmn(iDmn)%phys .EQ. phys_shell) THEN
            lSFp => lPD%get(ctmp,"Follower pressure load")
            IF (ASSOCIATED(lSFp)) THEN
               ALLOCATE(lEq%dmn(iDmn)%shlFp)
               CALL READShlFp(lSFp, lEq%dmn(iDmn)%shlFp, ctmp)
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE READDOMAIN
!--------------------------------------------------------------------
      SUBROUTINE READShlFp(list, lFp, ctmp)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(listType), INTENT(INOUT) :: list
      TYPE(shlFpType), INTENT(INOUT) :: lFp
      CHARACTER(LEN=stdL), INTENT(IN) :: ctmp

      INTEGER i, j
      TYPE(fileType) fTmp

      SELECT CASE (TRIM(ctmp))
      CASE ("steady")
         lFp%bType = IBSET(lFp%bType, bType_std)
         lPtr => list%get(lFp%p,"Value",1)
      CASE ("unsteady")
         lFp%bType = IBSET(lFp%bType, bType_ustd)
         ALLOCATE(lFp%pt)
         lPtr => list%get(ftmp, "Temporal values file path")
         IF (ASSOCIATED(lPtr)) THEN
            fid = fTmp%open()
            READ(fid,*) i, j
            IF (i .LT. 2) THEN
               std = "Enter nPnts nFCoef; nPts*(t Q)"
               err = "Wrong format in: "//fTmp%fname
            END IF
            lFp%pt%n = j
            ALLOCATE(lFp%pt%r(j))
            ALLOCATE(lFp%pt%i(j))
            CALL FFT(fid, i, lFp%pt)
            CLOSE(fid)
         ELSE
            lPtr => list%get(fTmp,"Fourier coefficients file path",1)
            fid = fTmp%open()
            READ (fid,*) lFp%pt%ti
            READ (fid,*) lFp%pt%T
            READ (fid,*) lFp%pt%qi
            READ (fid,*) lFp%pt%qs
            READ (fid,*) j
            lFp%pt%n = j
            ALLOCATE(lFp%pt%r(j))
            ALLOCATE(lFp%pt%i(j))
            DO i=1, j
               READ (fid,*) lFp%pt%r(i), lFp%pt%i(i)
            END DO
            CLOSE(fid)
         END IF
      CASE DEFAULT
         err = "Undefined follower load type"
      END SELECT

      RETURN
      END SUBROUTINE READShlFp
!--------------------------------------------------------------------
      END SUBROUTINE READEQ
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
      CASE ("Coupled Momentum","CMM")
         lBc%bType = IBSET(lBc%bType,bType_CMM)
      CASE ("Periodic","Per")
         lBc%bType = IBSET(lBc%bType,bType_per)
         err = "Periodic BC hasn't been implemented yet"
      CASE DEFAULT
         err = TRIM(list%ping("Type",lPtr))//" Unexpected BC type"
      END SELECT

!     Allocate traction
      ALLOCATE(lBc%h(nsd))
      lBc%h = 0D0

!     Weak Dirichlet BC for fluid/FSI equations
      lBc%weakDir = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir)) THEN
         IF (phys.EQ.phys_fluid .OR. phys.EQ.phys_FSI)
     2      lPtr => list%get(lBc%weakDir, "Weakly applied")
      END IF

!     Read penalty values for weakly applied Dir BC
      IF (lBc%weakDir) THEN
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

      ALLOCATE(lBc%eDrn(nsd))
      lBc%eDrn = 0
      lPtr => list%get(lBc%eDrn,"Effective direction")

      ctmp = "Steady"
      lPtr => list%get(ctmp,"Time dependence")
      SELECT CASE (ctmp)
      CASE ('Steady')
         lBc%bType = IBSET(lBc%bType,bType_std)
         IF (BTEST(lBc%bType, bType_trac)) THEN
            lPtr => list%get(lBc%h,"Value",1)
         ELSE
            lPtr => list%get(lBc%g,"Value",1)
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
         IF (phys.NE.phys_fluid .AND.
     2       phys.NE.phys_FSI   .AND.
     3       phys.NE.phys_CMM) THEN
            err = "Resistance is only defined for fluid/FSI equations"
         END IF

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
     2      " is violated in "//fTmp%fname

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
     2         "detected for BC. Mesh: "//TRIM(msh(iM)%name)//
     3         ", Face: "//TRIM(msh(iM)%fa(iFa)%name)//
     4         ", Node: "//STR(a)//" gN: "//STR(msh(iM)%fa(iFa)%gN(a))
            ptr(Ac) = a
         END DO
         DO i=1, j
            READ (fid,*) rtmp
            lBc%gm%t(i) = rtmp
            IF (i .EQ. 1) THEN
               IF (.NOT.ISZERO(rtmp)) err = "First time step"//
     2            " should be zero in <"//TRIM(ctmp)//">"
            ELSE
               rtmp = rtmp - lBc%gm%t(i-1)
               IF (ISZERO(rtmp) .OR. rtmp.LT.0D0) err = "Non-in"//
     2            "creasing time trend is found in <"//TRIM(ctmp)//">"
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
      CASE DEFAULT
         err=TRIM(list%ping("Time dependence",lPtr))//" Unexpected type"
      END SELECT

!     To impose value or flux
      ltmp = .FALSE.
      lPtr => list%get(ltmp,"Impose flux")
      IF (ltmp) lBc%bType = IBSET(lBc%bType,bType_flx)

!     To zero-out perimeter or not. Default is .true. for Dir/CMM
      ltmp = .FALSE.
      IF (BTEST(lBc%bType,bType_Dir) .OR.
     2    BTEST(lBc%bType,bType_CMM)) ltmp = .TRUE.
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
            lBc%bType = IBCLR(lBc%bType,bType_ddep)
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

      RETURN
      END SUBROUTINE READBC
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
            lEq%output(iOut)%name = "Fiber"
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
            SELECT CASE(TRIM(ctmp))
            CASE('fsiLS', 'FSILS', 'svFSI')
               lEq%ls%PREC_Type = PREC_FSILS
               useTrilinosLS = .FALSE.
#ifdef WITH_TRILINOS
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
      LOGICAL flag
      REAL(KIND=8) rtmp
      CHARACTER(LEN=stdL) ctmp

      lPtr => lPD%get(ctmp, "Electrophysiology model")
      SELECT CASE(TRIM(ctmp))
      CASE ("AP", "ap", "Aliev-Panfilov", "aliev-panfilov")
         lDmn%cep%cepType = cepModel_AP
      CASE ("FN", "fn", "Fitzhugh-Nagumo", "fitzhugh-nagumo")
         lDmn%cep%cepType = cepModel_FN
      CASE ("TTP", "tTP", "ttp", "tenTusscher-Panfilov")
         lDmn%cep%cepType = cepModel_TTP
      CASE DEFAULT
         err = "Unknown electrophysiology model"
      END SELECT

      lPtr => lPD%get(lDmn%cep%nX, "State variables", ll=1)

      lPtr => lPD%get(lDmn%cep%Diso,"Conductivity (iso)",ll=0D0)
      lDmn%cep%nFn = lPD%srch("Conductivity (ani)")
      IF (lDmn%cep%nFn .EQ. 0) lDmn%cep%nFn = 1

      ALLOCATE(lDmn%cep%Dani(lDmn%cep%nFn))
      lDmn%cep%Dani = 0D0
      IF (lDmn%cep%nFn .NE. 0) THEN
         DO i=1, lDmn%cep%nFn
            lPtr => lPD%get(lDmn%cep%Dani(i), "Conductivity (ani)", i)
         END DO
      END IF

      lDmn%cep%Istim%A  = 0.0d0
      lDmn%cep%Istim%Ts = HUGE(rtmp)
      lDmn%cep%Istim%Tp = HUGE(rtmp)
      lDmn%cep%Istim%Td = 0.0d0

      list => lPD%get(ctmp, "Stimulus")
      IF (ASSOCIATED(list)) THEN
         lPtr => list%get(lDmn%cep%Istim%A, "Amplitude")
         IF (.NOT.ISZERO(lDmn%cep%Istim%A)) THEN
            lPtr => list%get(lDmn%cep%Istim%Ts, "Start time")
            lPtr => list%get(lDmn%cep%Istim%Td, "Duration")
            lPtr => list%get(lDmn%cep%Istim%Tp, "Period")
         END IF
      END IF

      lDmn%cep%odes%tIntType = tIntType_FE
      lPtr => lPD%get(ctmp, "ODE solver")
      IF (ASSOCIATED(lPtr)) THEN
         SELECT CASE (ctmp)
         CASE ("FE", "Euler", "Explicit")
            lDmn%cep%odes%tIntType = tIntType_FE
         CASE ("RK4", "Runge")
            lDmn%cep%odes%tIntType = tIntType_RK4
         CASE ("CN2", "Implicit")
            lDmn%cep%odes%tIntType = tIntType_CN2
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

      lDmn%cep%Vrst = -80.0D0
      lPtr => lPD%get(rtmp, "Resting potential")
      IF (ASSOCIATED(lPtr)) lDmn%cep%Vrst = rtmp

      lDmn%cep%Kmef = 0D0
      lPtr => lPD%get(rtmp, "Parameter for mechano-electric feedback")
      IF (ASSOCIATED(lPtr)) lDmn%cep%Kmef = rtmp
      IF (.NOT.cplEM) lDmn%cep%Kmef = 0D0

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

!     Use gFa%w as temporary array to load traction data
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

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
!     All the data structures are defined in this module.
!
!--------------------------------------------------------------------

      MODULE COMMOD
      USE CMMOD
      USE CHNLMOD

      INCLUDE "FSILS.h"
      INCLUDE "cplBC.h"

!--------------------------------------------------------------------
!     Constants and upperbounds

!     maximum possible nsd, standard length for strings,
!     random for history file handel, maximum possible number
!     of output variables, size of blocks for openMP communications,
!     master is assumed to have zero ID, version, maximum number of
!     properties, License expiration date
      INTEGER, PARAMETER :: maxnsd = 3, version = 8, maxNProp = 30,
     3   expDate(3)=(/2020,1,1/)

!     Gauss points and their corresponding weights, upto 5 points
      REAL(KIND=8), PARAMETER :: gW(5,5)=RESHAPE((/2D0,0D0,0D0,0D0,0D0,
     2   1D0,1D0,0D0,0D0,0D0, 0.5555555555555556D0,0.8888888888888889D0,
     3   0.5555555555555556D0,0D0,0D0, 0.3478548451374538D0,
     4   0.6521451548625462D0,0.6521451548625462D0,0.3478548451374538D0,
     5   0D0, 0.236926885056189D0,0.4786286704993665D0,
     6   0.5688888888888889D0,0.4786286704993665D0,0.236926885056189D0/)
     7   ,(/5,5/))
      REAL(KIND=8), PARAMETER :: gXi(5,5)=RESHAPE((/0D0,0D0,0D0,0D0,0D0,
     2   -0.57735026918962584D0,0.57735026918962584D0,0D0,0D0,0D0,
     3   -0.7745966692414834D0,0D0,0.7745966692414834D0,0D0,0D0,
     4   -0.86113631159405257D0,-0.33998104358485631D0,
     5   0.33998104358485631D0,0.86113631159405257D0,0D0,
     6   -0.90617984593866396D0,-0.53846931010568311D0,0D0,
     7   0.53846931010568311D0,0.90617984593866396D0/),(/5,5/))

      CHARACTER, PARAMETER :: delimiter = "/"

!     Possible physical properties. Current maxNPror is 20.
!     When adding more properties, remember to increase maxNProp
!     Density of fluid, viscosity of fluid, density of solid, elsticity
!     modulus, Poisson's ratio, conductivity, internal force(X,Y,Z),
!     particles diameter, particle density, stabilization coefficient
!     for backflow divergence
      INTEGER, PARAMETER :: prop_NA = 0, fluid_density = 1,
     2   viscosity = 2, solid_density = 3, elasticity_modulus = 4,
     3   poisson_ratio = 5, conductivity = 6, f_x = 7, f_y = 8, f_z = 9,
     4   particle_diameter = 10, particle_density = 11,
     5   permeability = 12, backflow_stab = 13, source_term = 14,
     6   damping = 15

!     Types of accepted elements
!     Linear (1D), triangle (2D), tetrahedral (3D), bilinear (2D), quad
!     (1D), biquad (2D), brick (3D), general NURBS (1-3D)
      INTEGER, PARAMETER :: eType_NA = 100, eType_LIN = 101,
     2   eType_TRI = 102, eType_TET = 103, eType_BIL = 104,
     3   eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107,
     4   eType_NRB = 108, eType_WDG = 109

!     Types of equations that are included in this solver
!     Fluid equation (Navier-Stokes), structure (non-linear), heat
!     equation, linear elasticity, heat in a fluid
!     (advection-diffusion), fluid-structure-interaction, elector
!     magnetic, mesh motion, and Basset-Boussinesq Oseen equation
      INTEGER, PARAMETER :: phys_NA = 200, phys_fluid = 201,
     2   phys_struct = 202, phys_heatS = 203, phys_lElas = 204,
     3   phys_heatF = 205, phys_FSI = 206, phys_elcMag = 207,
     4   phys_mesh = 208, phys_BBO = 209

!     Saving formats
!     Don't save, VTK ASCII format, VTK binary format
      INTEGER, PARAMETER :: saveF_NA = 300, saveF_none = 301,
     2   saveF_VTK = 302, saveF_VTKB = 303

!     Differenty type of coupling for cplBC
!     Not-available, implicit, semi-implicit, and explicit
      INTEGER, PARAMETER :: cplBC_NA = 400, cplBC_I = 401,
     2   cplBC_SI = 402, cplBC_E = 403

!     boundary conditions types. Items of this list can be combined
!     BCs from imposing perspective can be Neu/Dir/per
!     BCs time dependence can be std/ustd/cpl/gen/res
!     BCs spatial distribution can be flat/para/ud
!     Beside these nodes at the boundary perimeter can be set to
!     zero and flux through surface can be assigned instead of nodal
!     values.
!     For non-matching meshes, variables can be projected at the
!     interface.
!     Drichlet, Neumann, periodic, steady, unsteady, coupled,
!     general (combination of ud/ustd), resistance, flat profile,
!     parabolic profile, user defined profile, zero out perimeter,
!     impose flux instead of value, L2 projection, backflow
!     stabilization, impose BC on the integral of state variable or D
!     (instead of Y), diplacement dependent
      INTEGER, PARAMETER :: bType_Dir = 0, bType_Neu = 1,
     2   bType_per = 2, bType_std = 3, bType_ustd = 4, bType_cpl = 5,
     3   bType_gen = 6, bType_res = 7, bType_flx = 8, bType_flat = 9,
     4   bType_para = 10, bType_ud = 11, bType_zp = 12, bType_bfs = 13,
     5   bType_impD = 14, bType_ddep = 15

!     Possible senarios for the output, followed by the possible outputs
!     Undefined output, extract it from A, extract it from Y, extract it
!     from D, calculate WSS, calculate vorticity, energy flux, heat
!     flux, absolute velocity (for FSI)
      INTEGER, PARAMETER :: outGrp_NA = 500, outGrp_A = 501,
     2   outGrp_Y = 502, outGrp_D = 503, outGrp_WSS = 504,
     3   outGrp_vort = 505, outGrp_eFlx = 506, outGrp_hFlx = 507,
     4   outGrp_absV = 508, outGrp_stInv = 509, outGrp_vortex = 510

      INTEGER, PARAMETER :: out_velocity = 599, out_pressure = 598,
     2   out_acceleration = 597, out_temperature = 596, out_WSS = 595,
     3   out_vorticity = 594, out_displacement = 593,
     4   out_energyFlux = 592, out_heatFlux = 591,
     5   out_absVelocity = 590, out_strainInv = 589, out_vortex = 588
      INTEGER, PARAMETER :: RMSH_TETGEN = 1, RMSH_MESHSIM = 2

!     Type of constitutive model for structure equation
      INTEGER, PARAMETER :: cModel_NA = 600, cModel_stVK = 601,
     2   cModel_mStVK = 602, cModel_nHook = 603, cModel_Guccione = 604,
     3   cModel_reducedHolzapfel = 605

!     Preconditioner definitions
      INTEGER, PARAMETER :: PREC_NONE = 700, PREC_FSILS = 701,
     2   PREC_TRILINOS_DIAGONAL = 702, PREC_TRILINOS_BLOCK_JACOBI = 703,
     3   PREC_TRILINOS_ILU = 704, PREC_TRILINOS_ILUT = 705,
     4   PREC_TRILINOS_IC = 706, PREC_TRILINOS_ICT = 707,
     5   PREC_TRILINOS_ML = 708

!     Solver definitions
      INTEGER, PARAMETER :: lSolver_NA = 799, lSolver_CG=798,
     3   lSolver_GMRES=797, lSolver_NS=796, lSolver_BICGS = 795,
     4   lSolver_MINRES = 794, lSolver_MUMPS = 793,
     5   lSolver_UMFPACK = 792, lSolver_PASTIX = 791,
     6   lSolver_SUPERLU = 790

!--------------------------------------------------------------------
!     Here comes subTypes definitions later used in other derived types
!     This is the container for B-Splines
      TYPE bsType
!        Number of knots (p + nNo + 1)
         INTEGER :: n = 0
!        Number of Gauss points for integration
         INTEGER nG
!        Number of knot spans (element)
         INTEGER :: nEl = 0
!        Number of control points (nodes)
         INTEGER :: nNo = 0
!        Number of sample points in each element (for output)
         INTEGER nSl
!        The order
         INTEGER p
!        Knot vector
         REAL(KIND=8), ALLOCATABLE :: xi(:)
      END TYPE bsType

!     Domain type is to keep track with element belong to which domain
!     and also different hysical quantities
      TYPE dmnType
!        The domain ID. Default includes entire domain
         INTEGER :: Id = -1
!        which physics must be solved in this domain
         INTEGER :: phys
!        Type of constitutive model (only for struct/FSI)
         INTEGER :: cModel = cModel_NA
!        The volume of this domain
         REAL(KIND=8) :: v = 0D0
!        physical properties, such as viscosity, density, ...
         REAL(KIND=8) :: prop(maxNProp) = 0D0
      END TYPE dmnType

!     The face type containing mesh at boundary
      TYPE faceType
!        Parametric direction normal to this face (NURBS)
         INTEGER d
!        Number of nodes (control points) in a single element
         INTEGER eNoN
!        Element type
         INTEGER :: eType = eType_NA
!        The mesh index that this face belongs to
         INTEGER :: iM
!        Number of elements
         INTEGER :: nEl = 0
!        Global number of elements
         INTEGER :: gnEl = 0
!        Number of Gauss points for integration
         INTEGER nG
!        Number of nodes
         INTEGER :: nNo = 0
!        Global element Ids
         INTEGER, ALLOCATABLE :: gE(:)
!        Global node Ids
         INTEGER, ALLOCATABLE :: gN(:)
!        Connectivity array
         INTEGER, ALLOCATABLE :: IEN(:,:)
!        EBC array (gE + gIEN)
         INTEGER, ALLOCATABLE :: gebc(:,:)
!        Surface area
         REAL(KIND=8) area
!        Gauss point weights
         REAL(KIND=8), ALLOCATABLE :: w(:)
!        Position coordinates
         REAL(KIND=8), ALLOCATABLE :: x(:,:)
!        Shape functions at Gauss points
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
!        Normal vector to each nodal point
         REAL(KIND=8), ALLOCATABLE :: nV(:,:)
!        Shape functions derivative at Gauss points
         REAL(KIND=8), ALLOCATABLE :: Nx(:,:,:)
!        Face name for flux files
         CHARACTER(LEN=stdL) name
      END TYPE faceType

!     Declared type for outputed variables
      TYPE outputType
!        Is this output suppose to be written into VTK, boundary, vol
         LOGICAL :: wtn(3) = .FALSE.
!        The group that this belong to (one of outType_*)
         INTEGER :: grp = outGrp_NA
!        Length of the outputed variable
         INTEGER l
!        Offset from the first index
         INTEGER o
!        The name to be used for the output and also in input file
         CHARACTER(LEN=stdL) name
      END TYPE outputType

!     Moving boundary data structure (used for general BC)
      TYPE MBType
!     Degrees of freedom of d(:,.,.)
         INTEGER dof
!     Number of time points to be read
         INTEGER :: nTP = 0
!     The period of data
         REAL(KIND=8) period
!     Time points
         REAL(KIND=8), ALLOCATABLE :: t(:)
!     Displacements at each direction, location, and time point
         REAL(KIND=8), ALLOCATABLE :: d(:,:,:)
      END TYPE MBType

!     Fourier coefficients that are used to specify unsteady BCs
      TYPE fcType
!        Number of Fourier coefficient
         INTEGER :: n = 0
!        Initial value
         REAL(KIND=8) qi
!        Time derivative of linear part
         REAL(KIND=8) qs
!        Period
         REAL(KIND=8) T
!        Initial time
         REAL(KIND=8) ti
!        Imaginary part of coefficint
         REAL(KIND=8), ALLOCATABLE :: i(:)
!        Real part of coefficint
         REAL(KIND=8), ALLOCATABLE :: r(:)
      END TYPE fcType

!     Boundary condition data type
      TYPE bcType
!        Pre/Res/Flat/Para... boundary types
         INTEGER :: bType = 0
!        Pointer to coupledBC%face
         INTEGER :: cplBCptr = 0
!        The face index that corresponds to this BC
         INTEGER iFa
!        The mesh index that corresponds to this BC
         INTEGER iM
!        Pointer to FSILS%bc
         INTEGER lsPtr
!        Defined steady value
         REAL(KIND=8) :: g = 0D0
!        Neu: defined resistance
         REAL(KIND=8) :: r = 0D0
!        Direction vector for imposing the BC
         INTEGER, ALLOCATABLE :: eDrn(:)
!        Spatial depanadant BC (profile data)
         REAL(KIND=8), ALLOCATABLE :: gx(:)
!        General BC (unsteady and UD combination)
         TYPE(MBType), ALLOCATABLE :: gm
!        Time depandant BC (Unsteady imposed value)
         TYPE(fcType), ALLOCATABLE :: gt
      END TYPE bcType

!     Linear system of equations solver type
      TYPE lsType
!        LS solver                     (IN)
         INTEGER LS_type
!        Preconditioner                (IN)
         INTEGER PREC_Type
!        Successful solving            (OUT)
         LOGICAL :: suc = .FALSE.
!        Maximum iterations            (IN)
         INTEGER :: mItr = 1000
!        Space dimension               (IN)
         INTEGER sD
!        Number of iteration           (OUT)
         INTEGER itr
!        Number of Ax multiple         (OUT)
         INTEGER cM
!        Number of |x| norms           (OUT)
         INTEGER cN
!        Number of <x.y> dot products  (OUT)
         INTEGER cD
!        Only for data alignment       (-)
         INTEGER reserve
!        Absolute tolerance            (IN)
         REAL(KIND=8) :: absTol = 1D-8
!        Relative tolerance            (IN)
         REAL(KIND=8) :: relTol = 1D-8
!        Initial norm of residual      (OUT)
         REAL(KIND=8) iNorm
!        Final norm of residual        (OUT)
         REAL(KIND=8) fNorm
!        Res. rduction in last itr.    (OUT)
         REAL(KIND=8) dB
!        Calling duration              (OUT)
         REAL(KIND=8) callD
!        Solver options file           (IN)
         TYPE(fileType) :: optionsFile
      END TYPE lsType

!--------------------------------------------------------------------
!     All the subTypes are defined, now defining the major types that
!     will be directly allocated

!     For coupled 0D-3D problems
      TYPE cplBCType
!        Is multi-domain active
         LOGICAL :: coupled = .FALSE.
!        Number of coupled faces
         INTEGER :: nFa = 0
!        Number of unknowns in the 0D domain
         INTEGER :: nX = 0
!        Implicit/Explicit/Semi-implicit schemes
         INTEGER :: schm = cplBC_NA
!        Path to the 0D code binary file
         CHARACTER(LEN=stdL) :: binPath
!        File name for communication between 0D and 3D
         CHARACTER(LEN=stdL) :: commuName = ".CPLBC_0D_3D.tmp"
!        The name of history file containing "X"
         CHARACTER(LEN=stdL) :: saveName = "LPN.dat"
!        New time step unknowns in the 0D domain
         REAL(KIND=8), ALLOCATABLE :: xn(:)
!        Old time step unknowns in the 0D domain
         REAL(KIND=8), ALLOCATABLE :: xo(:)
!        Data structure used for communicating with 0D code
         TYPE(cplFaceType), ALLOCATABLE :: fa(:)
      END TYPE cplBCType

!     This is the container for a mesh or NURBS patch, those specific
!     to NURBS are noted
      TYPE mshType
!        Whether the shape function is linear
         LOGICAL lShpF
!        Element type
         INTEGER :: eType = eType_NA
!        Number of nodes (control points) in a single element
         INTEGER eNoN
!        Global number of elements (knot spanes)
         INTEGER :: gnEl = 0
!        Global number of nodes (control points)
         INTEGER :: gnNo = 0
!        Number of element face
         INTEGER nEf
!        Number of elements (knot spanes)
         INTEGER :: nEl = 0
!        Number of faces
         INTEGER :: nFa = 0
!        Number of Gauss points for integration
         INTEGER nG
!        Number of nodes (control points)
         INTEGER :: nNo = 0
!        Number of elements sample points to be outputs (NURBS)
         INTEGER nSl
!        the element type recognized by VTK format
         INTEGER vtkType
!        Element distribution between processors
         INTEGER, ALLOCATABLE :: eDist(:)
!        Element domain ID number
         INTEGER, ALLOCATABLE :: eId(:)
!        Global nodes maping nNo --> tnNo
         INTEGER, ALLOCATABLE :: gN(:)
!        GLobal projected nodes mapping
!        projected -> unprojected mapping
         INTEGER, ALLOCATABLE :: gpN(:)
!        Global connectivity array mappig eNoN,nEl --> gnNo
         INTEGER, ALLOCATABLE :: gIEN(:,:)
!        The connectivity array mapping eNoN,nEl --> nNo
         INTEGER, ALLOCATABLE :: IEN(:,:)
!        gIEN mapper from old to new
         INTEGER, ALLOCATABLE :: otnIEN(:)
!        Local knot pointer (NURBS)
         INTEGER, ALLOCATABLE :: INN(:,:)
!        Global to local maping tnNo --> nNo
         INTEGER, ALLOCATABLE :: lN(:)
!        Control points weights (NURBS)
         REAL(KIND=8), ALLOCATABLE :: nW(:)
!        Gauss weights
         REAL(KIND=8), ALLOCATABLE :: w(:)
!        Position coordinates
         REAL(KIND=8), ALLOCATABLE :: x(:,:)
!        Parent shape function
         REAL(KIND=8), ALLOCATABLE :: N(:,:)
!        Parent shape functions gradient
         REAL(KIND=8), ALLOCATABLE :: Nx(:,:,:)
!        Mesh Name
         CHARACTER(LEN=stdL) :: name
!        BSpline in different directions (NURBS)
         TYPE(bsType), ALLOCATABLE :: bs(:)
!        Faces are stored in this variable
         TYPE(faceType), ALLOCATABLE :: fa(:)
      END TYPE mshType

!     Equation type
      TYPE eqType
!        Should be satisfied in a coupled/uncoupled fashion
         LOGICAL :: coupled = .TRUE.
!        Satisfied/not satisfied
         LOGICAL ok
!        Degrees of freedom
         INTEGER :: dof = 0
!        Pointer to end of unknown Yo(:,s:e)
         INTEGER e
!        Number of performed iterations
         INTEGER itr
!        Maximum iteration for this eq.
         INTEGER :: maxItr = 5
!        Minimum iteration for this eq.
         INTEGER :: minItr = 1
!        Number of possible outputs
         INTEGER :: nOutput = 0
!        Number of domains
         INTEGER :: nDmn = 0
!        Number of BCs
         INTEGER :: nBc = 0
!        Type of equation fluid/heatF/heatS/lElas/FSI
         INTEGER phys
!        Pointer to start of unknown Yo(:,s:e)
         INTEGER s
!        \alpha_f
         REAL(KIND=8) af
!        \alpha_m
         REAL(KIND=8) am
!        \beta
         REAL(KIND=8) beta
!        dB reduction in residual
         REAL(KIND=8) :: dBr = -4D1
!        \gamma
         REAL(KIND=8) gam
!        Initial norm of residual
         REAL(KIND=8) iNorm
!        First iteration norm
         REAL(KIND=8) pNorm
!        \rho_{infinity}
         REAL(KIND=8) roInf
!        Accepted relative tolerance
         REAL(KIND=8) :: tol = 1D64
!        Equation symbol
         CHARACTER(LEN=2) :: sym = "NA"
!        type of linear solver
         TYPE(lsType) ls
!        FSILS type of linear solver
         TYPE(FSILS_lsType) FSILS
!        BCs associated with this equation
         TYPE(bcType), ALLOCATABLE :: bc(:)
!        domains that this equation must be solved
         TYPE(dmnType), ALLOCATABLE :: dmn(:)
!        Outputs
         TYPE(outputType), ALLOCATABLE :: output(:)
      END TYPE eqType

!     This type will be used to write data in the VTK files.
      TYPE dataType
!        Element number of nodes
         INTEGER eNoN
!        Number of elements
         INTEGER nEl
!        Number of nodes
         INTEGER nNo
!        vtk type
         INTEGER vtkType
!        Connectivity array
         INTEGER, ALLOCATABLE :: IEN(:,:)
!        Element based variables to be written
         INTEGER, ALLOCATABLE :: xe(:,:)
!        All the variables after transformation to global format
         REAL(KIND=8), ALLOCATABLE :: gx(:,:)
!        All the variables to be written (including position)
         REAL(KIND=8), ALLOCATABLE :: x(:,:)
      END TYPE dataType

      TYPE rmshType
!     Whether remesh is required for problem or not
         LOGICAL :: isReqd
!     Method for remeshing: 1-TetGen, 2-MeshSim
         INTEGER :: method
!     Counter to track number of remesh done
         INTEGER :: cntr
!     Time step from which remeshing is done
         INTEGER :: rTS
!     Time step freq for saving data
         INTEGER :: cpVar
!     Time step at which forced remeshing is done
         INTEGER :: fTS
!     Time step frequency for forced remeshing
         INTEGER :: freq
!     Time where remeshing starts
         REAL(KIND=8) :: time
!     Mesh quality parameters
         REAL(KIND=8) :: minDihedAng
         REAL(KIND=8) :: maxRadRatio
!     Edge size of mesh
         REAL(KIND=8), ALLOCATABLE :: maxEdgeSize(:)
!     Initial norm of an equation
         REAL(KIND=8), ALLOCATABLE :: iNorm(:)
!     Copy of solution variables where remeshing starts
         REAL(KIND=8), ALLOCATABLE :: A0(:,:)
         REAL(KIND=8), ALLOCATABLE :: Y0(:,:)
         REAL(KIND=8), ALLOCATABLE :: D0(:,:)
!     Solution variables used for averaging
         REAL(KIND=8), ALLOCATABLE :: Aav(:,:)
         REAL(KIND=8), ALLOCATABLE :: Yav(:,:)
         REAL(KIND=8), ALLOCATABLE :: Dav(:,:)
!     Flag is set if remeshing is required for each mesh
         LOGICAL, ALLOCATABLE :: flag(:)
      END TYPE rmshType

!--------------------------------------------------------------------
!     All the types are defined, time to use them

!     LOGICAL VARIABLES
!     Whether there is a requirement to update mesh and Dn-Do variables
      LOGICAL dFlag
!     Whether mesh is moving
      LOGICAL mvMsh
!     Whether to averaged results
      LOGICAL saveAve
!     Whether any file being saved
      LOGICAL savedOnce
!     Whether to use separator in output
      LOGICAL sepOutput
!     Whether start from beginning or from simulations
      LOGICAL stFileFlag
!     Whether to overwrite restart file or not
      LOGICAL stFileRepl
!     Use Legacy mesh input format
      LOGICAL legacyFmt
!     Restart simulation after remeshing
      LOGICAL resetSim
!     Check IEN array for initial mesh
      LOGICAL ichckIEN
!     Reset averaging variables from zero
      LOGICAL zeroAve
!     Use C++ Trilinos framework for the linear solvers
      LOGICAL useTrilinosLS
!     Use C++ Trilinos framework for assembly and for linear solvers
      LOGICAL useTrilinosAssemAndLS

!     INTEGER VARIABLES
!     Current domain
      INTEGER cDmn
!     Current equation
      INTEGER cEq
!     Current time step
      INTEGER cTS
!     Current equation degrees of freedom
      INTEGER dof
!     Global total number of nodes
      INTEGER gtnNo
!     Number of equations
      INTEGER nEq
!     Number of faces in the LHS passed to FSILS
      INTEGER nFacesLS
!     Number of meshes
      INTEGER nMsh
!     Number of spatial dimensions
      INTEGER nsd
!     Number of time steps
      INTEGER nTS
!     Number of initialization time steps
      INTEGER nITS
!     stFiles record length
      INTEGER recLn
!     Start saving after this number of time step
      INTEGER saveATS
!     Format of output file
      INTEGER saveFormat
!     Increment in saving solutions
      INTEGER saveIncr
!     Stamp ID to make sure simulation is compatible with stFiles
      INTEGER stamp(8)
!     Increment in saving restart file
      INTEGER stFileIncr
!     Total number of degrees of freedom per node
      INTEGER tDof
!     Total number of nodes
      INTEGER tnNo
!     Restart Time Step
      INTEGER rsTS

!     REAL VARIABLES
!     Time step size
      REAL(KIND=8) dt
!     Time
      REAL(KIND=8) time

!     CHARACTER VARIABLES
!     Initialization file path
      CHARACTER(LEN=stdL) iniFilePath
!     Saved output file name
      CHARACTER(LEN=stdL) saveName
!     Restart file name
      CHARACTER(LEN=stdL) stFileName
!     Stop_trigger file name
      CHARACTER(LEN=stdL) stopTrigName

!     ALLOCATABLE DATA
!     Column pointer (for sparse LHS matrix structure)
      INTEGER, ALLOCATABLE :: colPtr(:)
!     Domain ID
      INTEGER, ALLOCATABLE :: dmnId(:)
!     Local to global pointer tnNo --> gtnNo
      INTEGER, ALLOCATABLE :: ltg(:)
!     Row pointer (for sparse LHS matrix structure)
      INTEGER, ALLOCATABLE :: rowPtr(:)

!     Old time derivative of variables (acceleration)
      REAL(KIND=8), ALLOCATABLE :: Ao(:,:)
!     New time derivative of variables
      REAL(KIND=8), ALLOCATABLE :: An(:,:)
!     Old integrated variables (dissplacement)
      REAL(KIND=8), ALLOCATABLE :: Do(:,:)
!     New integrated variables
      REAL(KIND=8), ALLOCATABLE :: Dn(:,:)
!     Residual vector
      REAL(KIND=8), ALLOCATABLE :: R(:,:)
!     LHS matrix
      REAL(KIND=8), ALLOCATABLE :: Val(:,:)
!     Position vector
      REAL(KIND=8), ALLOCATABLE :: x(:,:)
!     Old variables (velocity)
      REAL(KIND=8), ALLOCATABLE :: Yo(:,:)
!     New variables
      REAL(KIND=8), ALLOCATABLE :: Yn(:,:)
!     Fiber direction (for electrophysiology / structure mechanics)
      REAL(KIND=8), ALLOCATABLE :: fN(:,:)

!     DERIVED TYPE VARIABLES
!     Coupled BCs structures used for multidomain simulations
      TYPE(cplBCType), SAVE :: cplBC
!     All data related to equations are stored in this container
      TYPE(eqType), ALLOCATABLE :: eq(:)
!     FSILS data structure to produce LHS sparse matrix
      TYPE(FSILS_lhsType) lhs
!     All the meshes are stored in this variable
      TYPE(mshType), ALLOCATABLE :: msh(:)
!     Input/output to the screen is handled by this structure
      TYPE(chnlType), POINTER :: std, err, wrn, dbg
!     To group above channels
      TYPE(ioType), TARGET :: io
!     The general communicator
      TYPE(cmType) cm
!     Remesher type
      TYPE(rmshType) :: rmsh

      END MODULE COMMOD

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
      USE CEPMOD

      INCLUDE "FSILS.h"

!--------------------------------------------------------------------
!     Constants and upperbounds
      INCLUDE "CONSTS.f"

!--------------------------------------------------------------------
!     Here comes subTypes definitions later used in other derived types
!     Function spaces (basis) type
      TYPE fsType
!        Whether the basis function is linear
         LOGICAL lShpF
!        Element type
         INTEGER(KIND=IKIND) :: eType = eType_NA
!        Number of basis functions, typically equals msh%eNoN
         INTEGER(KIND=IKIND) eNoN
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        Gauss weights
         REAL(KIND=RKIND), ALLOCATABLE :: w(:)
!        Gauss integration points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Bounds on Gauss integration points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xib(:,:)
!        Parent shape function
         REAL(KIND=RKIND), ALLOCATABLE :: N(:,:)
!        Bounds on shape functions
         REAL(KIND=RKIND), ALLOCATABLE :: Nb(:,:)
!        Parent shape functions gradient
         REAL(KIND=RKIND), ALLOCATABLE :: Nx(:,:,:)
!        Second derivatives of shape functions - used for shells & IGA
         REAL(KIND=RKIND), ALLOCATABLE :: Nxx(:,:,:)
      END TYPE fsType

!     This is the container for B-Splines
      TYPE bsType
!        Number of knots (p + nNo + 1)
         INTEGER(KIND=IKIND) :: n = 0
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        Number of knot spans (element)
         INTEGER(KIND=IKIND) :: nEl = 0
!        Number of control points (nodes)
         INTEGER(KIND=IKIND) :: nNo = 0
!        Number of sample points in each element (for output)
         INTEGER(KIND=IKIND) nSl
!        The order
         INTEGER(KIND=IKIND) p
!        Knot vector
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:)
      END TYPE bsType

!     Fourier coefficients that are used to specify unsteady BCs
      TYPE fcType
!        If this is a ramp function
         LOGICAL lrmp
!        Number of Fourier coefficient
         INTEGER(KIND=IKIND) :: n = 0
!        No. of dimensions (scalar or vector)
         INTEGER(KIND=IKIND) :: d
!        Initial value
         REAL(KIND=RKIND), ALLOCATABLE :: qi(:)
!        Time derivative of linear part
         REAL(KIND=RKIND), ALLOCATABLE :: qs(:)
!        Period
         REAL(KIND=RKIND) T
!        Initial time
         REAL(KIND=RKIND) ti
!        Imaginary part of coefficint
         REAL(KIND=RKIND), ALLOCATABLE :: i(:,:)
!        Real part of coefficint
         REAL(KIND=RKIND), ALLOCATABLE :: r(:,:)
      END TYPE fcType

!     Moving boundary data structure (used for general BC)
      TYPE MBType
!     Degrees of freedom of d(:,.,.)
         INTEGER(KIND=IKIND) dof
!     Number of time points to be read
         INTEGER(KIND=IKIND) :: nTP = 0
!     The period of data
         REAL(KIND=RKIND) period
!     Time points
         REAL(KIND=RKIND), ALLOCATABLE :: t(:)
!     Displacements at each direction, location, and time point
         REAL(KIND=RKIND), ALLOCATABLE :: d(:,:,:)
      END TYPE MBType

      TYPE rcrType
!        Proximal resistance
         REAL(KIND=RKIND) :: Rp = 0._RKIND
!        Capacitance
         REAL(KIND=RKIND) :: C  = 0._RKIND
!        Distance resistance
         REAL(KIND=RKIND) :: Rd = 0._RKIND
!        Distal pressure
         REAL(KIND=RKIND) :: Pd = 0._RKIND
!        Initial value
         REAL(KIND=RKIND) :: Xo = 0._RKIND
      END TYPE rcrType

!     Boundary condition data type
      TYPE bcType
!        Strong/Weak application of Dirichlet BC
         LOGICAL :: weakDir
!        Whether load vector changes with deformation
!        (Neu - struct/ustruct only)
         LOGICAL :: flwP = .FALSE.
!        Robin: apply only in normal direction
         LOGICAL :: rbnN = .FALSE.
!        Pre/Res/Flat/Para... boundary types
         INTEGER(KIND=IKIND) :: bType = 0
!        Pointer to coupledBC%face
         INTEGER(KIND=IKIND) :: cplBCptr = 0
!        The face index that corresponds to this BC
         INTEGER(KIND=IKIND) iFa
!        The mesh index that corresponds to this BC
         INTEGER(KIND=IKIND) iM
!        Pointer to FSILS%bc
         INTEGER(KIND=IKIND) lsPtr
!        Clamped Neu BC master node parameter
         INTEGER(KIND=IKIND) masN
!        Defined steady value
         REAL(KIND=RKIND) :: g = 0._RKIND
!        Neu: defined resistance
         REAL(KIND=RKIND) :: r = 0._RKIND
!        Robin: stiffness
         REAL(KIND=RKIND) :: k = 0._RKIND
!        Robin: damping
         REAL(KIND=RKIND) :: c = 0._RKIND
!        Penalty parameters for weakly applied Dir BC
         REAL(KIND=RKIND) :: tauB(2) = 0._RKIND
!        Direction vector for imposing the BC
         INTEGER(KIND=IKIND), ALLOCATABLE :: eDrn(:)
!        Defined steady vector (traction)
         REAL(KIND=RKIND), ALLOCATABLE :: h(:)
!        Spatial dependant BC (profile data)
         REAL(KIND=RKIND), ALLOCATABLE :: gx(:)
!        General BC (unsteady and UD combination)
         TYPE(MBType), ALLOCATABLE :: gm
!        Time dependant BC (Unsteady imposed value)
         TYPE(fcType), ALLOCATABLE :: gt
!        Neu: RCR
         TYPE(rcrType) :: RCR
      END TYPE bcType

!     Body force data structure type
      TYPE bfType
!        Type of body force applied
         INTEGER(KIND=IKIND) :: bType = 0
!        No. of dimensions (1 or nsd)
         INTEGER(KIND=IKIND) :: dof
!        Mesh index corresponding to this body force
         INTEGER(KIND=IKIND) :: iM
!        Steady value
         REAL(KIND=RKIND), ALLOCATABLE :: b(:)
!        Steady but spatially dependant
         REAL(KIND=RKIND), ALLOCATABLE :: bx(:,:)
!        Time dependant (unsteady imposed value)
         TYPE(fcType), ALLOCATABLE :: bt
!        General (unsteady and spatially dependent combination)
         TYPE(MBType), ALLOCATABLE :: bm
      END TYPE bfType

!     Imposed fiber stress type
      TYPE fibStrsType
!        Time dependence of fiber stress (steady/unsteady)
         INTEGER(KIND=IKIND) :: fType = 0
!        Constant steady value
         REAL(KIND=8) :: g = 0._RKIND
!        Unsteady time-dependent values
         TYPE(fcType) :: gt
!        Cross-fiber (sheet) stress parameter
         REAL(KIND=RKIND) :: eta_s = 0._RKIND
      END TYPE fibStrsType

!     Structural domain type
      TYPE stModelType
!        Type of constitutive model (volumetric) for struct/FSI
         INTEGER(KIND=IKIND) :: volType = stVol_NA
!        Penalty parameter
         REAL(KIND=RKIND) :: Kpen = 0._RKIND
!        Type of constitutive model (isochoric) for struct/FSI
         INTEGER(KIND=IKIND) :: isoType = stIso_NA
!        Parameters specific to the constitutive model (isochoric)
!        NeoHookean model (C10 = mu/2)
         REAL(KIND=RKIND) :: C10 = 0._RKIND
!        Mooney-Rivlin model (C10, C01)
         REAL(KIND=RKIND) :: C01 = 0._RKIND
!        Holzapfel-Ogden / HGO model
!        (a, b, aff, bff, ass, bss, afs, bfs, kap, khs)
         REAL(KIND=RKIND) :: a   = 0._RKIND
         REAL(KIND=RKIND) :: b   = 0._RKIND
         REAL(KIND=RKIND) :: aff = 0._RKIND
         REAL(KIND=RKIND) :: bff = 0._RKIND
         REAL(KIND=RKIND) :: ass = 0._RKIND
         REAL(KIND=RKIND) :: bss = 0._RKIND
         REAL(KIND=RKIND) :: afs = 0._RKIND
         REAL(KIND=RKIND) :: bfs = 0._RKIND
!        Collagen fiber dispersion parameter (HGO model)
         REAL(KIND=RKIND) :: kap = 0._RKIND
!        Heaviside function parameter (Holzapfel-Ogden model)
         REAL(KIND=RKIND) :: khs = 100._RKIND
!        Lee-Sacks model
         REAL(KIND=RKIND) :: a0   = 0._RKIND
         REAL(KIND=RKIND) :: b1   = 0._RKIND
         REAL(KIND=RKIND) :: b2   = 0._RKIND
         REAL(KIND=RKIND) :: mu0  = 0._RKIND
!        Fiber reinforcement stress
         TYPE(fibStrsType) :: Tf
      END TYPE stModelType

!     Excitation-contraction model type for electromechanics
      TYPE eccModelType
!        Active stress coupling
         LOGICAL :: astress = .FALSE.
!        Active strain coupling
         LOGICAL :: astrain = .FALSE.

!        If excitation is coupled with cellular activation model or
!        imposed using an analytical function
         LOGICAL :: caCpld = .TRUE.
!        Type of active strain coupling
         INTEGER(KIND=IKIND) :: asnType = asnType_NA
!        Orthotropy parameter for active strain
         REAL(KIND=RKIND) :: k = 1._RKIND

!        Below variables are for decoupled excitation-contraction
!        Type of decoupling: analytically function or prescribed
         INTEGER :: dType = 0
!        Input parameters file path
         CHARACTER(LEN=stdL) :: fpar_in
!        Time step for integration
         REAL(KIND=RKIND) :: dt
!        Time integration options
         TYPE(odeType) :: odeS
!        State variable for excitation-contraction coupling
!          := activation force for active stress model
!          := fiber contraction parameter for active strain model
         REAL(KIND=RKIND) :: Ya = 0._RKIND
!        Cross-fiber stress component for active stress coupling
         REAL(KIND=RKIND) :: eta_s = 0._RKIND
!        Unsteady time-dependent values for prescribed fiber-shortening
         TYPE(fcType) :: Yat
      END TYPE eccModelType

!     Fluid viscosity model type
      TYPE viscModelType
!        Type of constitutive model for fluid viscosity
         INTEGER(KIND=IKIND) :: viscType = viscType_NA
!        Limiting zero shear-rate viscosity value
         REAL(KIND=RKIND) :: mu_o = 0._RKIND
!        Limiting high shear-rate viscosity (asymptotic) value
         REAL(KIND=RKIND) :: mu_i = 0._RKIND
!        Strain-rate tensor multiplier
         REAL(KIND=RKIND) :: lam = 0._RKIND
!        Strain-rate tensor exponent
         REAL(KIND=RKIND) :: a = 0._RKIND
!        Power-law exponent
         REAL(KIND=RKIND) :: n = 0._RKIND
      END TYPE viscModelType

!     Domain type is to keep track with element belong to which domain
!     and also different physical quantities
      TYPE dmnType
!        The domain ID. Default includes entire domain
         INTEGER(KIND=IKIND) :: Id = -1
!        Which physics must be solved in this domain
         INTEGER(KIND=IKIND) :: phys
!        The volume of this domain
         REAL(KIND=RKIND) :: v = 0._RKIND
!        General physical properties such as density, elastic modulus...
         REAL(KIND=RKIND) :: prop(maxNProp) = 0._RKIND
!        Electrophysiology model
         TYPE(cepModelType) :: cep
!        Structure material model
         TYPE(stModelType) :: stM
!        Excitation-contraction coupling
         TYPE(eccModelType) :: ec
!        Viscosity model for fluids
         TYPE(viscModelType) :: visc
      END TYPE dmnType

!     Mesh adjacency (neighboring element for each element)
      TYPE adjType
!        No of non-zeros
         INTEGER(KIND=IKIND) :: nnz = 0
!        Column pointer
         INTEGER(KIND=IKIND), ALLOCATABLE :: pcol(:)
!        Row pointer
         INTEGER(KIND=IKIND), ALLOCATABLE :: prow(:)
      END TYPE adjType

!     Tracer type used for immersed boundaries. Identifies traces of
!     nodes and integration points on background mesh elements
      TYPE traceType
!        No. of non-zero nodal traces
         INTEGER(KIND=IKIND) :: n = 0
!        No. of non-zero integration point traces
         INTEGER(KIND=IKIND) :: nG = 0
!        Self pointer of each trace to the IB global node
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
!        Self pointer of each trace to the IB integration point and
!        element ID
         INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:,:)
!        Nodal trace pointer array stores two values for each trace.
!        (1) background mesh element to which the trace points to,
!        (2) mesh ID
         INTEGER(KIND=IKIND), ALLOCATABLE :: nptr(:,:)
!        Integration point tracer array stores two values for each trace
!        (1) background mesh element to which the trace points to,
!        (2) mesh ID
         INTEGER(KIND=IKIND), ALLOCATABLE :: gptr(:,:)
!        Parametric coordinate for each nodal trace
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Parametric coordinate for each Gauss point trace
         REAL(KIND=RKIND), ALLOCATABLE :: xiG(:,:)
      END TYPE traceType

!     The face type containing mesh at boundary
      TYPE faceType
!        Parametric direction normal to this face (NURBS)
         INTEGER(KIND=IKIND) d
!        Number of nodes (control points) in a single element
         INTEGER(KIND=IKIND) eNoN
!        Element type
         INTEGER(KIND=IKIND) :: eType = eType_NA
!        The mesh index that this face belongs to
         INTEGER(KIND=IKIND) :: iM
!        Number of elements
         INTEGER(KIND=IKIND) :: nEl = 0
!        Global number of elements
         INTEGER(KIND=IKIND) :: gnEl = 0
!        Number of function spaces
         INTEGER(KIND=IKIND) nFs
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        Number of nodes
         INTEGER(KIND=IKIND) :: nNo = 0
!        Global element Ids
         INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:)
!        Global node Ids
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
!        Global to local maping tnNo --> nNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: lN(:)
!        Connectivity array
         INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:)
!        EBC array (gE + gIEN)
         INTEGER(KIND=IKIND), ALLOCATABLE :: gebc(:,:)
!        Surface area
         REAL(KIND=RKIND) area
!        Gauss point weights
         REAL(KIND=RKIND), ALLOCATABLE :: w(:)
!        Position coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!        Gauss points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Shape functions at Gauss points
         REAL(KIND=RKIND), ALLOCATABLE :: N(:,:)
!        Normal vector to each nodal point
         REAL(KIND=RKIND), ALLOCATABLE :: nV(:,:)
!        Shape functions derivative at Gauss points
         REAL(KIND=RKIND), ALLOCATABLE :: Nx(:,:,:)
!        Second derivatives of shape functions - for shells & IGA
         REAL(KIND=RKIND), ALLOCATABLE :: Nxx(:,:,:)
!        Face name for flux files
         CHARACTER(LEN=stdL) name
!        Face nodal adjacency
         TYPE(adjType) :: nAdj
!        Face element adjacency
         TYPE(adjType) :: eAdj
!        Function spaces (basis)
         TYPE(fsType), ALLOCATABLE :: fs(:)
      END TYPE faceType

!     Declared type for outputed variables
      TYPE outputType
!        Is this output suppose to be written into VTK, boundary, vol
         LOGICAL :: wtn(3) = .FALSE.
!        The group that this belong to (one of outType_*)
         INTEGER(KIND=IKIND) :: grp = outGrp_NA
!        Length of the outputed variable
         INTEGER(KIND=IKIND) l
!        Offset from the first index
         INTEGER(KIND=IKIND) o
!        The name to be used for the output and also in input file
         CHARACTER(LEN=stdL) name
      END TYPE outputType

!     Linear system of equations solver type
      TYPE lsType
!        LS solver                     (IN)
         INTEGER(KIND=IKIND) LS_type
!        Preconditioner                (IN)
         INTEGER(KIND=IKIND) PREC_Type
!        Successful solving            (OUT)
         LOGICAL :: suc = .FALSE.
!        Maximum iterations            (IN)
         INTEGER(KIND=IKIND) :: mItr = 1000
!        Space dimension               (IN)
         INTEGER(KIND=IKIND) sD
!        Number of iteration           (OUT)
         INTEGER(KIND=IKIND) itr
!        Number of Ax multiple         (OUT)
         INTEGER(KIND=IKIND) cM
!        Number of |x| norms           (OUT)
         INTEGER(KIND=IKIND) cN
!        Number of <x.y> dot products  (OUT)
         INTEGER(KIND=IKIND) cD
!        Only for data alignment       (-)
         INTEGER(KIND=IKIND) reserve
!        Absolute tolerance            (IN)
         REAL(KIND=RKIND) :: absTol = 1.E-12_RKIND
!        Relative tolerance            (IN)
         REAL(KIND=RKIND) :: relTol = 1.E-8_RKIND
!        Initial norm of residual      (OUT)
         REAL(KIND=RKIND) iNorm
!        Final norm of residual        (OUT)
         REAL(KIND=RKIND) fNorm
!        Res. rduction in last itr.    (OUT)
         REAL(KIND=RKIND) dB
!        Calling duration              (OUT)
         REAL(KIND=RKIND) callD
      END TYPE lsType

!     Contact model type
      TYPE cntctModelType
!        Contact model
         INTEGER(KIND=IKIND) :: cType = cntctM_NA
!        Penalty parameter
         REAL(KIND=RKIND) k
!        Min depth of penetration
         REAL(KIND=RKIND) h
!        Max depth of penetration
         REAL(KIND=RKIND) c
!        Min norm of face normals in contact
         REAL(KIND=RKIND) al
!        Potentail exponential index
         REAL(KIND=RKIND) p
!        Rin
         REAL(KIND=RKIND) Rin
!        Rout
         REAL(KIND=RKIND) Rout
!        Gap
         REAL(KIND=RKIND) gap
!        Tolerance
         REAL(KIND=RKIND) :: tol = 1.E-6_RKIND
      END TYPE cntctModelType

!--------------------------------------------------------------------
!     All the subTypes are defined, now defining the major types that
!     will be directly allocated

      TYPE cplFaceType
!        GenBC_Dir/GenBC_Neu
         INTEGER(KIND=IKIND) :: bGrp
!        Pointer to X
         INTEGER(KIND=IKIND) :: Xptr
!        Internal genBC use
         INTEGER(KIND=IKIND) :: eqv = 0
!        Flow rates at t
         REAL(KIND=RKIND) :: Qo = 0._RKIND
!        Flow rates at t+dt
         REAL(KIND=RKIND) :: Qn = 0._RKIND
!        Pressures at t
         REAL(KIND=RKIND) :: Po = 0._RKIND
!        Pressures at t+dt
         REAL(KIND=RKIND) :: Pn = 0._RKIND
!        Imposed flow/pressure
         REAL(KIND=RKIND) :: y = 0._RKIND
!        Name of the face
         CHARACTER(LEN=128) name
!        RCR type BC
         TYPE(rcrType) :: RCR
      END TYPE cplFaceType

!     For coupled 0D-3D problems
      TYPE cplBCType
!        Is multi-domain active
         LOGICAL :: coupled = .FALSE.
!        Whether to use genBC
         LOGICAL :: useGenBC = .FALSE.
!        Whether to initialize RCR from flow data
         LOGICAL :: initRCR = .FALSE.
!        Number of coupled faces
         INTEGER(KIND=IKIND) :: nFa = 0
!        Number of unknowns in the 0D domain
         INTEGER(KIND=IKIND) :: nX = 0
!        Number of output variables addition to nX
         INTEGER(KIND=IKIND) :: nXp = 0
!        Implicit/Explicit/Semi-implicit schemes
         INTEGER(KIND=IKIND) :: schm = cplBC_NA
!        Path to the 0D code binary file
         CHARACTER(LEN=stdL) :: binPath
!        File name for communication between 0D and 3D
         CHARACTER(LEN=stdL) :: commuName = ".CPLBC_0D_3D.tmp"
!        The name of history file containing "X"
         CHARACTER(LEN=stdL) :: saveName = "LPN.dat"
!        New time step unknowns in the 0D domain
         REAL(KIND=RKIND), ALLOCATABLE :: xn(:)
!        Old time step unknowns in the 0D domain
         REAL(KIND=RKIND), ALLOCATABLE :: xo(:)
!        Output variables to be printed
         REAL(KIND=RKIND), ALLOCATABLE :: xp(:)
!        Data structure used for communicating with 0D code
         TYPE(cplFaceType), ALLOCATABLE :: fa(:)
      END TYPE cplBCType

!     This is the container for a mesh or NURBS patch, those specific
!     to NURBS are noted
      TYPE mshType
!        Whether the shape function is linear
         LOGICAL lShpF
!        Whether the mesh is shell
         LOGICAL :: lShl = .FALSE.
!        Whether the mesh is fibers (Purkinje)
         LOGICAL :: lFib = .FALSE.
!        Element type
         INTEGER(KIND=IKIND) :: eType = eType_NA
!        Number of nodes (control points) in a single element
         INTEGER(KIND=IKIND) eNoN
!        Global number of elements (knot spans)
         INTEGER(KIND=IKIND) :: gnEl = 0
!        Global number of nodes (control points)
         INTEGER(KIND=IKIND) :: gnNo = 0
!        Number of element face. Used for reading Gambit mesh files
         INTEGER(KIND=IKIND) nEf
!        Number of elements (knot spans)
         INTEGER(KIND=IKIND) :: nEl = 0
!        Number of faces
         INTEGER(KIND=IKIND) :: nFa = 0
!        Number of function spaces
         INTEGER(KIND=IKIND) :: nFs
!        Number of Gauss points for integration
         INTEGER(KIND=IKIND) nG
!        Number of nodes (control points)
         INTEGER(KIND=IKIND) :: nNo = 0
!        Number of elements sample points to be outputs (NURBS)
         INTEGER(KIND=IKIND) nSl
!        The element type recognized by VTK format
         INTEGER(KIND=IKIND) vtkType
!        Number of fiber directions
         INTEGER(KIND=IKIND) nFn
!        Mesh scale factor
         REAL(KIND=RKIND) scF
!        IB: Mesh size parameter
         REAL(KIND=RKIND) dx
!        Element distribution between processors
         INTEGER(KIND=IKIND), ALLOCATABLE :: eDist(:)
!        Element domain ID number
         INTEGER(KIND=IKIND), ALLOCATABLE :: eId(:)
!        Global nodes maping nNo --> tnNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
!        GLobal projected nodes mapping
!        projected -> unprojected mapping
         INTEGER(KIND=IKIND), ALLOCATABLE :: gpN(:)
!        Global connectivity array mappig eNoN,nEl --> gnNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: gIEN(:,:)
!        The connectivity array mapping eNoN,nEl --> nNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:)
!        gIEN mapper from old to new
         INTEGER(KIND=IKIND), ALLOCATABLE :: otnIEN(:)
!        Local knot pointer (NURBS)
         INTEGER(KIND=IKIND), ALLOCATABLE :: INN(:,:)
!        Global to local maping tnNo --> nNo
         INTEGER(KIND=IKIND), ALLOCATABLE :: lN(:)
!        Shells: extended IEN array with neighboring nodes
         INTEGER(KIND=IKIND), ALLOCATABLE :: eIEN(:,:)
!        Shells: boundary condition variable
         INTEGER(KIND=IKIND), ALLOCATABLE :: sbc(:,:)
!        IB: Whether a cell is a ghost cell or not
         INTEGER(KIND=IKIND), ALLOCATABLE :: iGC(:)
!        Control points weights (NURBS)
         REAL(KIND=RKIND), ALLOCATABLE :: nW(:)
!        Gauss weights
         REAL(KIND=RKIND), ALLOCATABLE :: w(:)
!        Gauss integration points in parametric space
         REAL(KIND=RKIND), ALLOCATABLE :: xi(:,:)
!        Bounds on parameteric coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: xib(:,:)
!        Position coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!        Parent shape function
         REAL(KIND=RKIND), ALLOCATABLE :: N(:,:)
!        Shape function bounds
         REAL(KIND=RKIND), ALLOCATABLE :: Nb(:,:)
!        Normal vector to each nodal point (for Shells)
         REAL(KIND=RKIND), ALLOCATABLE :: nV(:,:)
!        Fiber orientations stored at the element level - used for
!        electrophysiology and solid mechanics
         REAL(KIND=RKIND), ALLOCATABLE :: fN(:,:)
!        Transmural coordinate for active strain type coupling
         REAL(KIND=RKIND), ALLOCATABLE :: tmX(:)
!        Parent shape functions gradient
         REAL(KIND=RKIND), ALLOCATABLE :: Nx(:,:,:)
!        Second derivatives of shape functions - used for shells & IGA
         REAL(KIND=RKIND), ALLOCATABLE :: Nxx(:,:,:)
!        Mesh Name
         CHARACTER(LEN=stdL) :: name
!        Mesh nodal adjacency
         TYPE(adjType) :: nAdj
!        Mesh element adjacency
         TYPE(adjType) :: eAdj
!        Function spaces (basis)
         TYPE(fsType), ALLOCATABLE :: fs(:)
!        BSpline in different directions (NURBS)
         TYPE(bsType), ALLOCATABLE :: bs(:)
!        Faces are stored in this variable
         TYPE(faceType), ALLOCATABLE :: fa(:)
!        IB: tracers
         TYPE(traceType) :: trc
      END TYPE mshType

!     Equation type
      TYPE eqType
!        Should be satisfied in a coupled/uncoupled fashion
         LOGICAL :: coupled = .TRUE.
!        Satisfied/not satisfied
         LOGICAL ok
!        Use C++ Trilinos framework for the linear solvers
         LOGICAL useTLS
!        Use C++ Trilinos framework for assembly and for linear solvers
         LOGICAL assmTLS
!        Degrees of freedom
         INTEGER(KIND=IKIND) :: dof = 0
!        Pointer to end of unknown Yo(:,s:e)
         INTEGER(KIND=IKIND) e
!        Number of performed iterations
         INTEGER(KIND=IKIND) itr
!        Maximum iteration for this eq.
         INTEGER(KIND=IKIND) :: maxItr = 5
!        Minimum iteration for this eq.
         INTEGER(KIND=IKIND) :: minItr = 1
!        Number of possible outputs
         INTEGER(KIND=IKIND) :: nOutput = 0
!        IB: Number of possible outputs
         INTEGER(KIND=IKIND) :: nOutIB = 0
!        Number of domains
         INTEGER(KIND=IKIND) :: nDmn = 0
!        IB: Number of immersed domains
         INTEGER(KIND=IKIND) :: nDmnIB = 0
!        Number of BCs
         INTEGER(KIND=IKIND) :: nBc = 0
!        IB: Number of BCs on immersed surfaces
         INTEGER(KIND=IKIND) :: nBcIB = 0
!        Number of BFs
         INTEGER(KIND=IKIND) :: nBf = 0
!        Type of equation fluid/heatF/heatS/lElas/FSI
         INTEGER(KIND=IKIND) phys
!        Pointer to start of unknown Yo(:,s:e)
         INTEGER(KIND=IKIND) s
!        \alpha_f
         REAL(KIND=RKIND) af
!        \alpha_m
         REAL(KIND=RKIND) am
!        \beta
         REAL(KIND=RKIND) beta
!        \gamma
         REAL(KIND=RKIND) gam
!        Initial norm of residual
         REAL(KIND=RKIND) iNorm
!        First iteration preconditioned relative residual norm
         REAL(KIND=RKIND) pNorm
!        \rho_{infinity}
         REAL(KIND=RKIND) roInf
!        Accepted relative tolerance
         REAL(KIND=RKIND) :: tol
!        Accepted absolute tolerance
         REAL(KIND=RKIND) :: absTol = 1.E-15_RKIND
!        Equation symbol
         CHARACTER(LEN=2) :: sym = "NA"
!        type of linear solver
         TYPE(lsType) ls
!        FSILS type of linear solver
         TYPE(FSILS_lsType) FSILS
!        BCs associated with this equation
         TYPE(bcType), ALLOCATABLE :: bc(:)
!        IB: BCs associated with this equation on immersed surfaces
         TYPE(bcType), ALLOCATABLE :: bcIB(:)
!        domains that this equation must be solved
         TYPE(dmnType), ALLOCATABLE :: dmn(:)
!        IB: immersed domains that this equation must be solved
         TYPE(dmnType), ALLOCATABLE :: dmnIB(:)
!        Outputs
         TYPE(outputType), ALLOCATABLE :: output(:)
!        IB: Outputs
         TYPE(outputType), ALLOCATABLE :: outIB(:)
!        Body force associated with this equation
         TYPE(bfType), ALLOCATABLE :: bf(:)
      END TYPE eqType

!     This type will be used to write data in the VTK files.
      TYPE dataType
!        Element number of nodes
         INTEGER(KIND=IKIND) eNoN
!        Number of elements
         INTEGER(KIND=IKIND) nEl
!        Number of nodes
         INTEGER(KIND=IKIND) nNo
!        vtk type
         INTEGER(KIND=IKIND) vtkType
!        Connectivity array
         INTEGER(KIND=IKIND), ALLOCATABLE :: IEN(:,:)
!        Element based variables to be written
         REAL(KIND=RKIND), ALLOCATABLE :: xe(:,:)
!        All the variables after transformation to global format
         REAL(KIND=RKIND), ALLOCATABLE :: gx(:,:)
!        All the variables to be written (including position)
         REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
      END TYPE dataType

      TYPE rmshType
!     Whether remesh is required for problem or not
         LOGICAL :: isReqd
!     Method for remeshing: 1-TetGen, 2-MeshSim
         INTEGER(KIND=IKIND) :: method
!     Counter to track number of remesh done
         INTEGER(KIND=IKIND) :: cntr
!     Time step from which remeshing is done
         INTEGER(KIND=IKIND) :: rTS
!     Time step freq for saving data
         INTEGER(KIND=IKIND) :: cpVar
!     Time step at which forced remeshing is done
         INTEGER(KIND=IKIND) :: fTS
!     Time step frequency for forced remeshing
         INTEGER(KIND=IKIND) :: freq
!     Time where remeshing starts
         REAL(KIND=RKIND) :: time
!     Mesh quality parameters
         REAL(KIND=RKIND) :: minDihedAng
         REAL(KIND=RKIND) :: maxRadRatio
!     Edge size of mesh
         REAL(KIND=RKIND), ALLOCATABLE :: maxEdgeSize(:)
!     Initial norm of an equation
         REAL(KIND=RKIND), ALLOCATABLE :: iNorm(:)
!     Copy of solution variables where remeshing starts
         REAL(KIND=RKIND), ALLOCATABLE :: A0(:,:)
         REAL(KIND=RKIND), ALLOCATABLE :: Y0(:,:)
         REAL(KIND=RKIND), ALLOCATABLE :: D0(:,:)
!     Flag is set if remeshing is required for each mesh
         LOGICAL, ALLOCATABLE :: flag(:)
      END TYPE rmshType

      TYPE ibCommType
!        Num traces (nodes) local to each process
         INTEGER(KIND=IKIND), ALLOCATABLE :: n(:)
!        Pointer to global trace (node num) stacked contiguously
         INTEGER(KIND=IKIND), ALLOCATABLE :: gN(:)
!        Num traces (Gauss points) local to each process
         INTEGER(KIND=IKIND), ALLOCATABLE :: nG(:)
!        Pointer to global trace (Gauss point) stacked contiguously
         INTEGER(KIND=IKIND), ALLOCATABLE :: gE(:)
      END TYPE ibCommType

!     Immersed Boundary (IB) data type
      TYPE ibType
!        Whether any file being saved
         LOGICAL :: savedOnce = .FALSE.
!        IB method
         INTEGER(KIND=IKIND) :: mthd = ibMthd_NA
!        IB coupling
         INTEGER(KIND=IKIND) :: cpld = ibCpld_NA
!        IB interpolation method
         INTEGER(KIND=IKIND) :: intrp = ibIntrp_NA
!        Current IB domain ID
         INTEGER(KIND=IKIND) :: cDmn
!        Current equation
         INTEGER(KIND=IKIND) :: cEq = 0
!        Total number of IB nodes
         INTEGER(KIND=IKIND) :: tnNo
!        Number of IB meshes
         INTEGER(KIND=IKIND) :: nMsh
!        IB call duration (1: total time; 2: update; 3,4: communication)
         REAL(KIND=RKIND) :: callD(4)
!        IB Domain ID
         INTEGER(KIND=IKIND), ALLOCATABLE :: dmnID(:)
!        Row pointer (for sparse LHS matrix storage)
         INTEGER(KIND=IKIND), ALLOCATABLE :: rowPtr(:)
!        Column pointer (for sparse LHS matrix storage)
         INTEGER(KIND=IKIND), ALLOCATABLE :: colPtr(:)
!        IB position coordinates
         REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!        Velocity (new)
         REAL(KIND=RKIND), ALLOCATABLE :: Yb(:,:)
!        Time derivative of displacement (old)
         REAL(KIND=RKIND), ALLOCATABLE :: Auo(:,:)
!        Time derivative of displacement (new)
         REAL(KIND=RKIND), ALLOCATABLE :: Aun(:,:)
!        Time derivative of displacement (n+am)
         REAL(KIND=RKIND), ALLOCATABLE :: Auk(:,:)
!        Displacement (old)
         REAL(KIND=RKIND), ALLOCATABLE :: Ubo(:,:)
!        Displacement (new)
         REAL(KIND=RKIND), ALLOCATABLE :: Ubn(:,:)
!        Displacement (n+af)
         REAL(KIND=RKIND), ALLOCATABLE :: Ubk(:,:)
!        Displacement (projected on background mesh, old)
         REAL(KIND=RKIND), ALLOCATABLE :: Uo(:,:)
!        Displacement (projected on background mesh, new, n+af)
         REAL(KIND=RKIND), ALLOCATABLE :: Un(:,:)
!        Residue (FSI force)
         REAL(KIND=RKIND), ALLOCATABLE :: R(:,:)
!        Residue (displacement, background mesh)
         REAL(KIND=RKIND), ALLOCATABLE :: Ru(:,:)
!        Residue (displacement, IB mesh)
         REAL(KIND=RKIND), ALLOCATABLE :: Rub(:,:)
!        LHS tangent matrix for displacement
         REAL(KIND=RKIND), ALLOCATABLE :: Ku(:,:)

!        DERIVED TYPE VARIABLES
!        IB meshes
         TYPE(mshType), ALLOCATABLE :: msh(:)
!        IB communicator
         TYPE(ibCommType) :: cm
      END TYPE ibType

!     Data type for Trilinos Linear Solver related arrays
      TYPE tlsType
!        Local to global mapping
         INTEGER(KIND=IKIND), ALLOCATABLE :: ltg(:)
!        Factor for Dirichlet BCs
         REAL(KIND=RKIND), ALLOCATABLE :: W(:,:)
!        Residue
         REAL(KIND=RKIND), ALLOCATABLE :: R(:,:)
      END TYPE tlsType
!--------------------------------------------------------------------
!     All the types are defined, time to use them

!     LOGICAL VARIABLES
!     Whether there is a requirement to update mesh and Dn-Do variables
      LOGICAL dFlag
!     Whether mesh is moving
      LOGICAL mvMsh
!     Whether to averaged results
      LOGICAL saveAve
!     Whether to save to VTK files
      LOGICAL saveVTK
!     Whether any file being saved
      LOGICAL savedOnce
!     Whether to use separator in output
      LOGICAL sepOutput
!     Whether start from beginning or from simulations
      LOGICAL stFileFlag
!     Whether to overwrite restart file or not
      LOGICAL stFileRepl
!     Restart simulation after remeshing
      LOGICAL resetSim
!     Check IEN array for initial mesh
      LOGICAL ichckIEN
!     Reset averaging variables from zero
      LOGICAL zeroAve
!     Whether CMM equation is initialized
      LOGICAL cmmInit
!     Whether variable wall properties are used for CMM
      LOGICAL cmmVarWall
!     Whether shell equation is being solved
      LOGICAL shlEq
!     Whether PRESTRESS is being solved
      LOGICAL pstEq
!     Whether velocity-pressure based structural dynamics solver is used
      LOGICAL sstEq
!     Whether excitation-contraction is coupled
      LOGICAL ecCpld
!     Whether to detect and apply any contact model
      LOGICAL iCntct
!     Whether any Immersed Boundary (IB) treatment is required
      LOGICAL ibFlag
!     Postprocess step - convert bin to vtk
      LOGICAL bin2VTK

!     INTEGER(KIND=IKIND) VARIABLES
!     Current domain
      INTEGER(KIND=IKIND) cDmn
!     Current equation
      INTEGER(KIND=IKIND) cEq
!     Current time step
      INTEGER(KIND=IKIND) cTS
!     Starting time step
      INTEGER(KIND=IKIND) startTS
!     Current equation degrees of freedom
      INTEGER(KIND=IKIND) dof
!     Global total number of nodes
      INTEGER(KIND=IKIND) gtnNo
!     Number of equations
      INTEGER(KIND=IKIND) nEq
!     Number of faces in the LHS passed to FSILS
      INTEGER(KIND=IKIND) nFacesLS
!     Number of meshes
      INTEGER(KIND=IKIND) nMsh
!     Number of spatial dimensions
      INTEGER(KIND=IKIND) nsd
!     Number of time steps
      INTEGER(KIND=IKIND) nTS
!     Number of initialization time steps
      INTEGER(KIND=IKIND) nITS
!     stFiles record length
      INTEGER(KIND=IKIND) recLn
!     Start saving after this number of time step
      INTEGER(KIND=IKIND) saveATS
!     Increment in saving solutions
      INTEGER(KIND=IKIND) saveIncr
!     Stamp ID to make sure simulation is compatible with stFiles
      INTEGER(KIND=IKIND) stamp(7)
!     Increment in saving restart file
      INTEGER(KIND=IKIND) stFileIncr
!     Total number of degrees of freedom per node
      INTEGER(KIND=IKIND) tDof
!     Total number of nodes
      INTEGER(KIND=IKIND) tnNo
!     Restart Time Step
      INTEGER(KIND=IKIND) rsTS
!     Number of stress values to be stored
      INTEGER(KIND=IKIND) nsymd

!     REAL VARIABLES
!     Time step size
      REAL(KIND=RKIND) dt
!     Time
      REAL(KIND=RKIND) time
!     Simulation starting time
      REAL(KIND=RKIND) :: start_time = 0._RKIND

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
      INTEGER(KIND=IKIND), ALLOCATABLE :: colPtr(:)
!     Domain ID
      INTEGER(KIND=IKIND), ALLOCATABLE :: dmnId(:)
!     Local to global pointer tnNo --> gtnNo
      INTEGER(KIND=IKIND), ALLOCATABLE :: ltg(:)
!     Row pointer (for sparse LHS matrix structure)
      INTEGER(KIND=IKIND), ALLOCATABLE :: rowPtr(:)
!     Array that maps global node id to rowN in the matrix
      INTEGER(KIND=IKIND), ALLOCATABLE :: idMap(:)

!     Boundary nodes set for CMM initialization and for zeroing-out
!     non-wall nodal displacements
      INTEGER(KIND=IKIND), ALLOCATABLE :: cmmBdry(:)

!     IB: iblank used for immersed boundaries (1 => solid, 0 => fluid)
      INTEGER, ALLOCATABLE :: iblank(:)

!     Old time derivative of variables (acceleration)
      REAL(KIND=RKIND), ALLOCATABLE :: Ao(:,:)
!     New time derivative of variables
      REAL(KIND=RKIND), ALLOCATABLE :: An(:,:)
!     Old integrated variables (dissplacement)
      REAL(KIND=RKIND), ALLOCATABLE :: Do(:,:)
!     New integrated variables
      REAL(KIND=RKIND), ALLOCATABLE :: Dn(:,:)
!     Residual vector
      REAL(KIND=RKIND), ALLOCATABLE :: R(:,:)
!     LHS matrix
      REAL(KIND=RKIND), ALLOCATABLE :: Val(:,:)
!     Position vector
      REAL(KIND=RKIND), ALLOCATABLE :: x(:,:)
!     Old variables (velocity)
      REAL(KIND=RKIND), ALLOCATABLE :: Yo(:,:)
!     New variables
      REAL(KIND=RKIND), ALLOCATABLE :: Yn(:,:)
!     Body force
      REAL(KIND=RKIND), ALLOCATABLE :: Bf(:,:)

!     Additional arrays for velocity-based formulation of nonlinear
!     solid mechanics
!     Time derivative of displacement
      REAL(KIND=RKIND), ALLOCATABLE :: Ad(:,:)
!     Residue of the displacement equation
      REAL(KIND=RKIND), ALLOCATABLE :: Rd(:,:)
!     LHS matrix for displacement equation
      REAL(KIND=RKIND), ALLOCATABLE :: Kd(:,:)

!     Variables for prestress calculations
      REAL(KIND=RKIND), ALLOCATABLE :: pS0(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: pSn(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: pSa(:)

!     Temporary storage for initializing state variables
      REAL(KIND=RKIND), ALLOCATABLE :: Pinit(:)
      REAL(KIND=RKIND), ALLOCATABLE :: Vinit(:,:)
      REAL(KIND=RKIND), ALLOCATABLE :: Dinit(:,:)

!     State variable for excitation-contraction coupling
!       := activation force for active stress model
!       := fiber contraction parameter for active strain model
      REAL(KIND=RKIND), ALLOCATABLE :: ec_Ya(:)

!     CMM-variable wall properties: 1-thickness, 2-Elasticity modulus
      REAL(KIND=RKIND), ALLOCATABLE :: varWallProps(:,:)

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
      TYPE(rmshType) rmsh
!     Contact model type
      TYPE(cntctModelType) cntctM
!     IB: Immersed boundary data structure
      TYPE(ibType), ALLOCATABLE :: ib
!     Trilinos Linear Solver data type
      TYPE(tlsType), ALLOCATABLE :: tls

      END MODULE COMMOD

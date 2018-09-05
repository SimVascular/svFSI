!#######################################################################
!     Parameters used for Aliev-Panfilov Ventricular Myocyte Model
!#######################################################################
!     Scaling factors
!     Voltage scaling
      REAL(KIND=8) :: Vscale  = 100.0D0
!     Time scaling
      REAL(KIND=8) :: Tscale  = 12.9D0
!     Voltage offset parameter
      REAL(KIND=8) :: Voffset = -80.0D0
!-----------------------------------------------------------------------
!     Model parameters
      REAL(KIND=8) :: alpha = 0.01D0
      REAL(KIND=8) :: a     = 0.002D0
      REAL(KIND=8) :: b     = 0.15D0
      REAL(KIND=8) :: c     = 8.0D0
      REAL(KIND=8) :: mu1   = 0.2D0
      REAL(KIND=8) :: mu2   = 0.3D0

!     Cm: Cell capacitance per unit surface area
      REAL(KIND=8) :: Cm  = 1.0d0
!     sV: Surface to volume ratio
      REAL(KIND=8) :: sV  = 1.0D0
!     rho: Cellular resistivity
      REAL(KIND=8) :: rho = 1.0D0
!#######################################################################

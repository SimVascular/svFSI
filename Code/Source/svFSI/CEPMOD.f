!#######################################################################
      MODULE CEPMOD
      USE APMOD
      USE FNMOD
      USE TTPMOD
      IMPLICIT NONE

!     Type of cardiac electrophysiology models: Aliev-Panfilov model,
!     tenTusscher-Panfilov model
      INTEGER, PARAMETER :: cepModel_NA = 100, cepModel_AP = 101,
     2   cepModel_FN = 102, cepModel_TTP = 103

!     Time integration scheme: Forward-Euler, Runge-Kutta 4th order,
!     Crank-Nicholson
      INTEGER, PARAMETER :: tIntType_NA  = 200, tIntType_FE = 201,
     2   tIntType_RK4 = 202, tIntType_CN2 = 203

!     Time integration scheme and related parameters
      TYPE odeType
!        Time integration method type
         INTEGER :: tIntType = tIntType_NA
!        Max. iterations for Newton-Raphson method
         INTEGER :: maxItr
!        Absolute tolerance
         REAL(KIND=8) :: absTol
!        Relative tolerance
         REAL(KIND=8) :: relTol
      END TYPE odeType

!     External stimulus type
      TYPE stimType
!        start time
         REAL(KIND=8) :: Ts
!        duration of stimulus
         REAL(KIND=8) :: Td
!        time period
         REAL(KIND=8) :: Tp
!        stimulus amplitude
         REAL(KIND=8) :: A
      END TYPE stimType

!     Cardiac electrophysiology model type
      TYPE cepModelType
!        Type of cardiac electrophysiology model
         INTEGER :: cepType = cepModel_NA
!        Number of unknowns
         INTEGER :: nX
!        Isotropic conductivity
         REAL(KIND=8) :: Diso = 0D0
!        Anisotropic conductivity
         REAL(KIND=8) :: Dani = 0D0
!        External stimulus
         TYPE(stimType) :: Istim
!        Time integration options
         TYPE(odeType) :: odeS
      END TYPE cepModelType

      END MODULE CEPMOD
!#######################################################################

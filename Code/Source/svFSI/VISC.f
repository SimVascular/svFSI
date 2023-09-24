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
!     Set of subroutines to read viscosity models for fluid and solid,
!     compute fluid's viscosity, and compute solid's viscous stresses
!     and stiffness matrices
!
!--------------------------------------------------------------------

!####################################################################
!     This subroutine reads parameters of a fluid's viscosity model
      SUBROUTINE READ_VISC_FLUID(lDmn, lPD)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(dmnType), INTENT(INOUT) :: lDmn
      TYPE(listType), INTENT(INOUT) :: lPD

      TYPE(listType), POINTER :: lPtr, lVis
      REAL(KIND=RKIND) :: rtmp
      CHARACTER(LEN=stdL) ctmp

      lVis => lPD%get(ctmp,"Viscosity",1)

      CALL TO_LOWER(ctmp)
      SELECT CASE (TRIM(ctmp))
      CASE ("constant", "const", "newtonian")
         lDmn%viscF%viscTypeF = viscTypeF_Const
         lPtr => lVis%get(lDmn%viscF%mu_i,"Value",1,lb=0._RKIND)

      CASE ("carreau-yasuda", "cy")
         lDmn%viscF%viscTypeF = viscTypeF_CY
         lPtr => lVis%get(lDmn%viscF%mu_i,
     2      "Limiting high shear-rate viscosity",1,lb=0._RKIND)
         lPtr => lVis%get(lDmn%viscF%mu_o,
     2      "Limiting low shear-rate viscosity",1,lb=0._RKIND)
         lPtr => lVis%get(lDmn%viscF%lam,
     2      "Shear-rate tensor multiplier (lamda)",1,lb=0._RKIND)
         lPtr => lVis%get(lDmn%viscF%a,
     2      "Shear-rate tensor exponent (a)",1,lb=0._RKIND)
         lPtr => lVis%get(lDmn%viscF%n,"Power-law index (n)",1,
     2      lb=0._RKIND)
         IF (lDmn%viscF%mu_i .GT. lDmn%viscF%mu_o) THEN
            err = "Unexpected inputs for Carreau-Yasuda model. "//
     2         "High shear-rate viscosity value should be higher than"//
     3         " low shear-rate value"
         END IF

      CASE ("cassons", "cass")
         lDmn%viscF%viscTypeF = viscTypeF_Cass
         lPtr => lVis%get(lDmn%viscF%mu_i,
     2      "Asymptotic viscosity parameter",1,lb=0._RKIND)
         lPtr => lVis%get(lDmn%viscF%mu_o,
     2      "Yield stress parameter",1,lb=0._RKIND)
         lDmn%viscF%lam = 0.5_RKIND
         lPtr => lVis%get(rtmp,"Low shear-rate threshold")
         IF (ASSOCIATED(lPtr)) lDmn%viscF%lam = rtmp

      CASE DEFAULT
         err = "Undefined constitutive model for viscosity used"
      END SELECT

      IF ((lDmn%phys .EQ. phys_stokes) .AND.
     2    (lDmn%viscF%viscTypeF .NE. viscTypeF_Const)) THEN
         err = "Only constant viscosity is allowed for Stokes flow"
      END IF

      RETURN
      END SUBROUTINE READ_VISC_FLUID
!####################################################################
!     This subroutine reads parameters of a solid's viscosity model
      SUBROUTINE READ_VISC_SOLID(lDmn, lPD)
      USE COMMOD
      USE LISTMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(dmnType), INTENT(INOUT) :: lDmn
      TYPE(listType), INTENT(INOUT) :: lPD

      TYPE(listType), POINTER :: lPtr, lVis
      CHARACTER(LEN=stdL) ctmp

      lVis => lPD%get(ctmp, "Viscosity", 1)

      CALL TO_LOWER(ctmp)
      SELECT CASE (TRIM(ctmp))
      CASE ("potential", "pot")
         lDmn%viscS%viscTypeS = viscTypeS_Pot
         lPtr => lVis%get(lDmn%viscS%mu, "Value", 1, ll=0._RKIND)

      CASE ("newtonian", "newt", "fluid")
         lDmn%viscS%viscTypeS = viscTypeS_Newt
         lPtr => lVis%get(lDmn%viscS%mu, "Value", 1, ll=0._RKIND)

      CASE DEFAULT
         err = "Undefined constitutive model for viscosity used"
      END SELECT

      RETURN
      END SUBROUTINE READ_VISC_SOLID
!####################################################################
!     Computes fluid's viscosity depending on the chosen model:
!     Newtonian or non-Newtonian
      SUBROUTINE GET_FLUID_VISC(lDmn, gamma, mu, mu_s, mu_x)
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(INOUT)  :: gamma
      REAL(KIND=RKIND), INTENT(OUT) :: mu, mu_s, mu_x

      REAL(KIND=RKIND) :: mu_i, mu_o, lam, a, n, T1, T2

      SELECT CASE (lDmn%viscF%viscTypeF)
!     Constant viscosity (Newtonian)
      CASE (viscTypeF_Const)
         mu   = lDmn%viscF%mu_i
         mu_s = mu
         mu_x = 0._RKIND

!     Carreau-Yasuda model
      CASE (viscTypeF_CY)
         mu_i = lDmn%viscF%mu_i
         mu_o = lDmn%viscF%mu_o
         lam  = lDmn%viscF%lam
         a    = lDmn%viscF%a
         n    = lDmn%viscF%n

         T1   = 1._RKIND + (lam*gamma)**a
         T2   = T1**((n-1._RKIND)/a)
         mu   = mu_i + (mu_o-mu_i)*T2
         mu_s = mu_i

         T1   = T2/T1
         T2   = lam**a * gamma**(a-1._RKIND) * T1
         mu_x = (mu_o-mu_i)*(n-1._RKIND)*T2

!     Cassons model
      CASE (viscTypeF_Cass)
         mu_i = lDmn%viscF%mu_i
         mu_o = lDmn%viscF%mu_o
         lam  = lDmn%viscF%lam

         IF (gamma .LT. lam) THEN
            mu_o  = mu_o/SQRT(lam)
            gamma = lam
         ELSE
            mu_o = mu_o/SQRT(gamma)
         END IF
         mu   = (mu_i + mu_o) * (mu_i + mu_o)
         mu_s = mu_i*mu_i
         mu_x = 2._RKIND*mu_o*(mu_o + mu_i)/gamma

      END SELECT

      RETURN
      END SUBROUTINE GET_FLUID_VISC
!####################################################################
!     Computes solid's viscous stress and stifness matrices
      SUBROUTINE GET_VISC_STRS(lDmn, eNoN, Nx, vx, F, Svis, Kvis_u,
     2   Kvis_v)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: Nx(nsd,eNoN), vx(nsd,nsd),
     2   F(nsd,nsd)
      REAL(KIND=RKIND), INTENT(OUT) :: Svis(nsd,nsd),
     2   Kvis_u(nsd*nsd,eNoN,eNoN), Kvis_v(nsd*nsd,eNoN,eNoN)

      IF (lDmn%viscS%viscTypeS .EQ. viscTypeS_Pot) THEN
         CALL GET_VISC_STRS_POT(lDmn%viscS%mu, eNoN, Nx, vx, F, Svis,
     2      Kvis_u, Kvis_v)

      ELSE IF (lDmn%viscS%viscTypeS .EQ. viscTypeS_Newt) THEN
         CALL GET_VISC_STRS_NEWT(lDmn%viscS%mu, eNoN, Nx, vx, F, Svis,
     2      Kvis_u, Kvis_v)
      END IF

      RETURN
      END SUBROUTINE GET_VISC_STRS
!--------------------------------------------------------------------
!     Computes solid's viscous stress and stifness matrices for
!     strain energy potential-based viscous model
      SUBROUTINE GET_VISC_STRS_POT(mu, eNoN, Nx, vx, F, Svis, Kvis_u,
     2   Kvis_v)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: mu, Nx(nsd,eNoN), vx(nsd,nsd),
     2   F(nsd,nsd)
      REAL(KIND=RKIND), INTENT(OUT) :: Svis(nsd,nsd),
     2   Kvis_u(nsd*nsd,eNoN,eNoN), Kvis_v(nsd*nsd,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j, ii
      REAL(KIND=RKIND) :: NxNx, FFt(nsd,nsd), FtVx(nsd,nsd),
     2   FVxt(nsd,nsd), FNx(nsd,eNoN), VxNx(nsd,eNoN)

      Svis   = 0._RKIND
      Kvis_u = 0._RKIND
      Kvis_v = 0._RKIND

!     Additional terms for computing stiffness due to viscous stress
      FFt  = MATMUL(F, TRANSPOSE(F))
      FtVx = MATMUL(TRANSPOSE(F), vx)
      FVxt = MATMUL(F, TRANSPOSE(vx))
      FNx  = 0._RKIND
      VxNx = 0._RKIND
      DO a=1, eNoN
         DO i=1, nsd
            DO j=1, nsd
               FNx(i,a)  = FNx(i,a)  + F(i,j)*Nx(j,a)
               VxNx(i,a) = VxNx(i,a) + vx(i,j)*Nx(j,a)
            END DO
         END DO
      END DO

!     2nd Piola-Kirchhoff stress due to viscosity
      Svis = mu*MAT_SYMM(FtVx, nsd)

!     Stiffness matrix due to solid viscosity
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx = 0._RKIND
            DO i=1, nsd
               NxNx = NxNx + Nx(i,a)*Nx(i,b)
            END DO

            DO i=1, nsd
               DO j=1, nsd
                  ii = (i-1)*nsd + j
                  Kvis_u(ii,a,b) = 0.5_RKIND*mu*
     2               (FNx(i,b)*VxNx(j,a) + NxNx*FVxt(i,j))
                  Kvis_v(ii,a,b) = 0.5_RKIND*mu*
     2               (NxNx*FFt(i,j) + FNx(i,b)*FNx(j,a))
               END DO
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE GET_VISC_STRS_POT
!--------------------------------------------------------------------
!     Computes solid's viscous stress and stifness matrices for
!     a Newtonian-fluid like viscosity model
      SUBROUTINE GET_VISC_STRS_NEWT(mu, eNoN, Nx, vx, F, Svis, Kvis_u,
     2   Kvis_v)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: eNoN
      REAL(KIND=RKIND), INTENT(IN) :: mu, Nx(nsd,eNoN), vx(nsd,nsd),
     2   F(nsd,nsd)
      REAL(KIND=RKIND), INTENT(OUT) :: Svis(nsd,nsd),
     2   Kvis_u(nsd*nsd,eNoN,eNoN), Kvis_v(nsd*nsd,eNoN,eNoN)

      INTEGER(KIND=IKIND) :: a, b, i, j, ii
      REAL(KIND=RKIND) :: Jac, NxNx, r2d, IDm(nsd,nsd), Fi(nsd,nsd),
     2   VxFi(nsd,nsd), ddev(nsd,nsd), NxFi(nsd,eNoN), DdNx(nsd,eNoN),
     3   VxNx(nsd,eNoN)

      Svis   = 0._RKIND
      Kvis_u = 0._RKIND
      Kvis_v = 0._RKIND

      IDm = MAT_ID(nsd)
      Jac = MAT_DET(F, nsd)
      Fi  = MAT_INV(F, nsd)

!     Velocity gradient in current configuration
      VxFi = MATMUL(vx, Fi)

!     Deviatoric strain tensor
      ddev = MAT_DEV(MAT_SYMM(VxFi,nsd), nsd)

!     Additional quantities required for computing stiffness
      NxFi = 0._RKIND
      DdNx = 0._RKIND
      VxNx = 0._RKIND
      DO a=1, eNoN
         DO i=1, nsd
            DO j=1, nsd
               NxFi(i,a) = NxFi(i,a) + Nx(j,a)*Fi(j,i)
            END DO
         END DO

         DO i=1, nsd
            DO j=1, nsd
               DdNx(i,a) = DdNx(i,a) + ddev(i,j)*NxFi(j,a)
               VxNx(i,a) = VxNx(i,a) + VxFi(j,i)*NxFi(j,a)
            END DO
         END DO
      END DO

!     2nd Piola-Kirchhoff stress due to viscosity
      Svis = MATMUL(ddev, TRANSPOSE(Fi))
      Svis = 2._RKIND*mu*Jac*MATMUL(Fi, Svis)

!     Stiffness matrix due to solid viscosity
      r2d = 2._RKIND / REAL(nsd, KIND=RKIND)
      DO b=1, eNoN
         DO a=1, eNoN
            NxNx = 0._RKIND
            DO i=1, nsd
               NxNx = NxNx + NxFi(i,a)*NxFi(i,b)
            END DO

            DO i=1, nsd
               DO j=1, nsd
                  ii = (i-1)*nsd + j

!              Derivative of the residue w.r.t. displacement
                  Kvis_u(ii,a,b) = mu*Jac*(2._RKIND*
     2               (DdNx(i,a)*NxFi(j,b) - DdNx(i,b)*NxFi(j,a)) -
     3               (NxNx*VxFi(i,j) + NxFi(i,b)*VxNx(j,a) -
     4                r2d*NxFi(i,a)*VxNx(j,b)))

!              Derivative of the residue w.r.t. velocity
                  Kvis_v(ii,a,b) = mu*Jac*(NxNx*IDm(i,j) +
     2               NxFi(i,b)*NxFi(j,a) - r2d*NxFi(i,a)*NxFi(j,b))
               END DO
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE GET_VISC_STRS_NEWT
!####################################################################

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
!     Here, second Piola-Kirchhoff stress tensor and the material
!     stiffness tensors are computed for material constitutive models.
!
!--------------------------------------------------------------------

!     Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
!     including both dilational and isochoric components
      SUBROUTINE GETPK2CC(lDmn, F, nfd, fl, tmX, ya, S, Dm)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      INTEGER(KIND=IKIND), INTENT(IN) :: nfd
      REAL(KIND=RKIND), INTENT(IN) :: F(nsd,nsd), fl(nsd,nfd), tmX, ya
      REAL(KIND=RKIND), INTENT(OUT) :: S(nsd,nsd), Dm(nsymd,nsymd)

      TYPE(stModelType) :: stM
      REAL(KIND=RKIND) :: nd, Kp, J, Ja, J2d, J4d, trE, p, pl, Inv1,
     2   Inv2, Inv4, Inv6, Inv8, Tfa, Tsa, IDm(nsd,nsd), C(nsd,nsd),
     3   E(nsd,nsd), Ci(nsd,nsd), Sb(nsd,nsd), CCb(nsd,nsd,nsd,nsd),
     4   PP(nsd,nsd,nsd,nsd), CC(nsd,nsd,nsd,nsd)
      REAL(KIND=RKIND) :: r1, r2, g1, g2, g3, rexp
!     Guccione
      REAL(KIND=RKIND) :: QQ, Rm(nsd,nsd), Es(nsd,nsd), RmRm(nsd,nsd,6)
!     HGO/HO model
      REAL(KIND=RKIND) :: Eff, Ess, Efs, kap, c4f, c4s, dc4f, dc4s,
     2   Hff(nsd,nsd), Hss(nsd,nsd), Hfs(nsd,nsd)
!     Active strain for electromechanics
      REAL(KIND=RKIND) :: Fe(nsd,nsd), Fa(nsd,nsd), Fai(nsd,nsd)

      S    = 0._RKIND
      Dm   = 0._RKIND

!     Some preliminaries
      stM  = lDmn%stM
      nd   = REAL(nsd, KIND=RKIND)
      Kp   = stM%Kpen

!     Fiber-reinforced stress
      CALL GET_FIB_STRESS(stM%Tf, Tfa)
      Tsa = Tfa*stM%Tf%eta_s

!     Electromechanics coupling - active stress
      IF (lDmn%ec%astress) THEN
         Tfa = Tfa + ya
         Tsa = Tsa + ya*lDmn%ec%eta_s
      END IF

!     Electromechanics coupling - active strain
      Fe   = F
      Fa   = MAT_ID(nsd)
      Fai  = Fa
      IF (lDmn%ec%astrain) THEN
         CALL GET_ACTV_DEFGRAD(lDmn%ec, nfd, fl, tmX, ya, Fa)
         Fai  = MAT_INV(Fa, nsd)
         Fe   = MATMUL(F, Fai)
      END IF

      Ja   = MAT_DET(Fa, nsd)
      J    = MAT_DET(Fe, nsd)
      J2d  = J**(-2._RKIND/nd)
      J4d  = J2d*J2d

      IDm  = MAT_ID(nsd)
      C    = MATMUL(TRANSPOSE(Fe), Fe)
      E    = 0.5_RKIND * (C - IDm)
      Ci   = MAT_INV(C, nsd)

      trE  = MAT_TRACE(E, nsd)
      Inv1 = J2d*MAT_TRACE(C,nsd)
      Inv2 = 0.5_RKIND*( Inv1*Inv1 - J4d*MAT_TRACE(MATMUL(C,C), nsd) )

!     Contribution of dilational penalty terms to S and CC
      p  = 0._RKIND
      pl = 0._RKIND
      IF (.NOT.ISZERO(Kp)) CALL GETSVOLP(stM, J, p, pl)

!     Now, compute isochoric and total stress, elasticity tensors
      SELECT CASE (stM%isoType)
      CASE (stIso_lin)
         g1 = stM%C10            ! mu
         S  = g1*Idm
         RETURN

!     St.Venant-Kirchhoff
      CASE (stIso_stVK)
         g1 = stM%C10            ! lambda
         g2 = stM%C01 * 2._RKIND      ! 2*mu

         S  = g1*trE*IDm + g2*E
         CC = g1*TEN_DYADPROD(IDm, IDm, nsd) + g2*TEN_IDs(nsd)

!     modified St.Venant-Kirchhoff
      CASE (stIso_mStVK)
         g1 = stM%C10 ! kappa
         g2 = stM%C01 ! mu

         S  = g1*LOG(J)*Ci + g2*(C-IDm)
         CC = g1 * ( -2._RKIND*LOG(J)*TEN_SYMMPROD(Ci, Ci, nsd) +
     2      TEN_DYADPROD(Ci, Ci, nsd) ) + 2._RKIND*g2*TEN_IDs(nsd)

!     NeoHookean model
      CASE (stIso_nHook)
         g1 = 2._RKIND * stM%C10
         Sb = g1*IDm

!        Fiber reinforcement/active stress
         Sb = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)

         r1 = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S  = J2d*Sb - r1*Ci
         CC = (-2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                      TEN_DYADPROD(S, Ci, nsd))

         S  = S + p*J*Ci
         CC = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd) +
     2         (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!     Mooney-Rivlin model
      CASE (stIso_MR)
         g1  = 2._RKIND * (stM%C10 + Inv1*stM%C01)
         g2  = -2._RKIND * stM%C01
         Sb  = g1*IDm + g2*J2d*C

!        Fiber reinforcement/active stress
         Sb  = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)

         g1  = 4._RKIND*J4d* stM%C01
         CCb = g1 * (TEN_DYADPROD(IDm, IDm, nsd) - TEN_IDs(nsd))

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC - (2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                           TEN_DYADPROD(S, Ci, nsd) )

         S   = S + p*J*Ci
         CC  = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd) +
     2          (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!     HGO (Holzapfel-Gasser-Ogden) model for arteries with isochoric
!     invariants for the anisotropy terms (decoupled)
      CASE (stIso_HGO_d)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "HGO material model (2)"
         kap  = stM%kap
         Inv4 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = J2d*NORM(fl(:,2), MATMUL(C, fl(:,2)))

         Eff  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv4 - 1._RKIND
         Ess  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv6 - 1._RKIND

         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Hff  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hff
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Hss  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hss

         g1   = stM%C10
         g2   = stM%aff * Eff * EXP(stM%bff*Eff*Eff)
         g3   = stM%ass * Ess * EXP(stM%bss*Ess*Ess)
         Sb   = 2._RKIND*(g1*IDm + g2*Hff + g3*Hss)

!        Fiber reinforcement/active + additional cross-fiber stresses
         Sb   = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
     2             + Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd)

         g1   = stM%aff*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)*
     2      EXP(stM%bff*Eff*Eff)
         g2   = stM%ass*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)*
     2      EXP(stM%bss*Ess*Ess)
         g1   = 4._RKIND*J4d*g1
         g2   = 4._RKIND*J4d*g2

         CCb  = g1 * TEN_DYADPROD(Hff, Hff, nsd) +
     2          g2 * TEN_DYADPROD(Hss, Hss, nsd)

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC - (2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                           TEN_DYADPROD(S, Ci, nsd) )

         S   = S + p*J*Ci
         CC  = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd) +
     2          (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!     HGO (Holzapfel-Gasser-Ogden) model for arteries with full
!     invariants for the anisotropy terms (modified-anisotropy)
      CASE (stIso_HGO_ma)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "HGO material model (2)"
         kap  = stM%kap
         Inv1 = MAT_TRACE(C,nsd)
         Inv4 = NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = NORM(fl(:,2), MATMUL(C, fl(:,2)))

         Eff  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv4 - 1._RKIND
         Ess  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv6 - 1._RKIND

!        Isochoric and volumetric contribution to stress and stiffness
!        tensors
         Sb   = 2._RKIND*stM%C10*IDm
         r1   = J2d*MAT_DDOT(C, Sb, nsd) / nd

         S    = J2d*Sb - r1*Ci
         CC   = -2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     2                           TEN_DYADPROD(S, Ci, nsd) )

         S    = S + p*J*Ci
         CC   = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd)
     2        +  (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!        Anisotropic contribution to stress and stiffness tensors
!        Fiber-Fiber interaction + additional fiber reinforcement
         rexp = EXP(stM%bff*Eff*Eff)
         g1   = 2._RKIND*stM%aff*Eff*rexp
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Hff  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hff
         S    = S + (g1*Hff) + (Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd))

         g1   = (1._RKIND + (2._RKIND*stM%bff*Eff*Eff))
         g1   = 4._RKIND*stM%aff*g1*rexp
         CC   = CC + (g1*TEN_DYADPROD(Hff, Hff, nsd))

!        Sheet-Sheet interaction + additional cross-fiber stress
         rexp = EXP(stM%bss*Ess*Ess)
         g2   = 2._RKIND*stM%ass*Ess*rexp
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Hss  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hss
         S    = S + (g2*Hss) + (Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd))

         g2   = (1._RKIND + (2._RKIND*stM%bss*Ess*Ess))
         g2   = 4._RKIND*stM%ass*g2*rexp
         CC   = CC + (g2*TEN_DYADPROD(Hss, Hss, nsd))

!     Guccione (1995) transversely isotropic model
      CASE (stIso_Gucci)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Guccione material model (2)"
!        Compute isochoric component of E
         E = 0.5_RKIND * (J2d*C - Idm)

!        Transform into local orthogonal coordinate system
         Rm(:,1) = fl(:,1)
         Rm(:,2) = fl(:,2)
         Rm(:,3) = CROSS(fl)

!        Project E to local orthogocal coordinate system
         Es = MATMUL(E, Rm)
         Es = MATMUL(TRANSPOSE(Rm), Es)

         g1 = stM%bff
         g2 = stM%bss
         g3 = stM%bfs

         QQ = g1 * Es(1,1)*Es(1,1) +
     2        g2 *(Es(2,2)*Es(2,2) + Es(3,3)*Es(3,3)  +
     3             Es(2,3)*Es(2,3) + Es(3,2)*Es(3,2)) +
     4        g3 *(Es(1,2)*Es(1,2) + Es(2,1)*Es(2,1)  +
     5             Es(1,3)*Es(1,3) + Es(3,1)*Es(3,1))

         r2 = stM%C10 * EXP(QQ)

!        Fiber stiffness contribution := (dE*_ab / dE_IJ)
         RmRm(:,:,1) = MAT_DYADPROD(Rm(:,1), Rm(:,1), nsd)
         RmRm(:,:,2) = MAT_DYADPROD(Rm(:,2), Rm(:,2), nsd)
         RmRm(:,:,3) = MAT_DYADPROD(Rm(:,3), Rm(:,3), nsd)

         RmRm(:,:,4) = MAT_SYMMPROD(Rm(:,1), Rm(:,2), nsd)
         RmRm(:,:,5) = MAT_SYMMPROD(Rm(:,2), Rm(:,3), nsd)
         RmRm(:,:,6) = MAT_SYMMPROD(Rm(:,3), Rm(:,1), nsd)

         Sb = g1*Es(1,1)*RmRm(:,:,1) + g2*(Es(2,2)*RmRm(:,:,2) +
     2      Es(3,3)*RmRm(:,:,3) + 2._RKIND*Es(2,3)*RmRm(:,:,5)) +
     4      2._RKIND*g3*(Es(1,2)*RmRm(:,:,4) + Es(1,3)*RmRm(:,:,6))

         CCb = 2._RKIND*TEN_DYADPROD(Sb, Sb, nsd)
         Sb  = Sb * r2

!        Fiber reinforcement/active + additional cross-fiber stresses
         Sb  = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
     2            + Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd)

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         r2  = r2*J4d
         CCb = r2*(CCb + g1*TEN_DYADPROD(RmRm(:,:,1), RmRm(:,:,1), nsd)
     2       + g2*(TEN_DYADPROD(RmRm(:,:,2), RmRm(:,:,2), nsd)
     3       + TEN_DYADPROD(RmRm(:,:,3), RmRm(:,:,3), nsd)
     4       + TEN_DYADPROD(RmRm(:,:,5), RmRm(:,:,5), nsd)*2._RKIND)
     5       + 2._RKIND*g3*(TEN_DYADPROD(RmRm(:,:,4), RmRm(:,:,4), nsd)
     6       + TEN_DYADPROD(RmRm(:,:,6), RmRm(:,:,6), nsd)))


         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC - (2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                           TEN_DYADPROD(S, Ci, nsd) )

         S   = S + p*J*Ci
         CC  = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd) +
     2          (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!     HO (Holzapfel-Ogden) model for myocardium with isoschoric
!     invariants for the anisotropy terms (decoupled)
      CASE (stIso_HO_d)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Holzapfel material model (2)"

!        Compute fiber-based invariants
         Inv4 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = J2d*NORM(fl(:,2), MATMUL(C, fl(:,2)))
         Inv8 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,2)))

         Eff  = Inv4 - 1._RKIND
         Ess  = Inv6 - 1._RKIND
         Efs  = Inv8

!        Smoothed heaviside function
         c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
         c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

!        Approx. derivative of smoothed heaviside function
         dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
         dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))
c         dc4f = stM%khs /
c     2      (EXP(stM%khs*Eff) + EXP(-stM%khs*Eff) + 2.0_RKIND)
c         dc4s = stM%khs /
c     2      (EXP(stM%khs*Ess) + EXP(-stM%khs*Ess) + 2.0_RKIND)

!        Isotropic + fiber-sheet interaction stress
         g1   = stM%a * EXP(stM%b*(Inv1-3._RKIND))
         g2   = 2._RKIND * stM%afs * EXP(stM%bfs*Efs*Efs)
         Hfs  = MAT_SYMMPROD(fl(:,1), fl(:,2), nsd)
         Sb   = g1*IDm + g2*Efs*Hfs

!        Isotropic + fiber-sheet interaction stiffness
         g1   = g1 * 2._RKIND*J4d*stM%b
         g2   = g2 * 2._RKIND*J4d*(1._RKIND + 2._RKIND*stM%bfs*Efs*Efs)
         CCb  = g1 * TEN_DYADPROD(IDm, IDm, nsd) +
     2          g2 * TEN_DYADPROD(Hfs, Hfs, nsd)

!        Fiber-fiber interaction stress + additional reinforcement (Tfa)
         rexp = EXP(stM%bff * Eff * Eff)
         g1   = c4f*Eff*rexp
         g1   = g1 + (0.5_RKIND*dc4f/stM%bff)*(rexp - 1._RKIND)
         g1   = 2._RKIND*stM%aff*g1 + Tfa
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Sb   = Sb + g1*Hff

!        Fiber-fiber interaction stiffness
         g1   = c4f*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)
         g1   = (g1 + 2._RKIND*dc4f*Eff)*rexp
         g1   = 4._RKIND*J4d*stM%aff*g1
         CCb  = CCb + g1*TEN_DYADPROD(Hff, Hff, nsd)

!        Sheet-sheet interaction stress + additional cross-fiber stress
         rexp = EXP(stM%bss * Ess * Ess)
         g2   = c4s*Ess*rexp
         g2   = g2 + (0.5_RKIND*dc4s/stM%bss)*(rexp - 1._RKIND)
         g2   = 2._RKIND*stM%ass*g2 + Tsa
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Sb   = Sb + g2*Hss

!        Sheet-sheet interaction stiffness
         g2   = c4s*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)
         g2   = (g2 + 2._RKIND*dc4s*Ess)*rexp
         g2   = 4._RKIND*J4d*stM%ass*g2
         CCb  = CCb + g2*TEN_DYADPROD(Hss, Hss, nsd)

!        Isochoric 2nd-Piola-Kirchhoff stress and stiffness tensors
         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC - (2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                           TEN_DYADPROD(S, Ci, nsd) )

!        Add pressure contribution
         S   = S + p*J*Ci
         CC  = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd) +
     2          (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!     HO (Holzapfel-Ogden) model for myocardium with full invariants
!     for the anisotropy terms (modified-anisotropy)
      CASE (stIso_HO_ma)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Holzapfel material model (2)"

!        Compute fiber-based full invariants (not isochoric)
         Inv4 = NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = NORM(fl(:,2), MATMUL(C, fl(:,2)))
         Inv8 = NORM(fl(:,1), MATMUL(C, fl(:,2)))

         Eff  = Inv4 - 1._RKIND
         Ess  = Inv6 - 1._RKIND
         Efs  = Inv8

!        Smoothed heaviside function
         c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
         c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

!        Approx. derivative of smoothed heaviside function
         dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
         dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))
c         dc4f = stM%khs /
c     2      (EXP(stM%khs*Eff) + EXP(-stM%khs*Eff) + 2.0_RKIND)
c         dc4s = stM%khs /
c     2      (EXP(stM%khs*Ess) + EXP(-stM%khs*Ess) + 2.0_RKIND)

!        Isochoric stress and stiffness
         g1   = stM%a * EXP(stM%b*(Inv1-3._RKIND))
         Sb   = g1*IDm
         r1   = J2d*MAT_DDOT(C, Sb, nsd) / nd

         g1   = g1 * 2._RKIND*J4d*stM%b
         CCb  = g1 * TEN_DYADPROD(IDm, IDm, nsd)

!        Add isochoric stress and stiffness contribution
         S    = J2d*Sb - r1*Ci

         PP   = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC   = TEN_DDOT(CCb, PP, nsd)
         CC   = TEN_TRANSPOSE(CC, nsd)
         CC   = TEN_DDOT(PP, CC, nsd)
         CC   = CC - (2._RKIND/nd) * ( TEN_DYADPROD(Ci, S, nsd) +
     2                                 TEN_DYADPROD(S, Ci, nsd) )

!        Add pressure contribution to stress and stiffness
         S    = S + p*J*Ci
         CC   = CC + 2._RKIND*(r1 - p*J) * TEN_SYMMPROD(Ci, Ci, nsd)
     2         + (pl*J - 2._RKIND*r1/nd) * TEN_DYADPROD(Ci, Ci, nsd)

!        Now that both isochoric and volumetric components were added,
!        anisotropic components need to be added

!        Fiber-sheet interaction terms
         g1   = 2._RKIND * stM%afs * EXP(stM%bfs*Efs*Efs)
         Hfs  = MAT_SYMMPROD(fl(:,1), fl(:,2), nsd)
         S    = S + (g1*Efs*Hfs)

         g1   = g1 * 2._RKIND*(1._RKIND + 2._RKIND*stM%bfs*Efs*Efs)
         CC   = CC + (g1*TEN_DYADPROD(Hfs, Hfs, nsd))

!        Fiber-fiber interaction stress + additional reinforcement (Tfa)
         rexp = EXP(stM%bff * Eff * Eff)
         g1   = c4f*Eff*rexp
         g1   = g1 + (0.5_RKIND*dc4f/stM%bff)*(rexp - 1._RKIND)
         g1   = (2._RKIND*stM%aff*g1) + Tfa
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         S    = S + (g1*Hff)

!        Fiber-fiber interaction stiffness
         g1   = c4f*(1._RKIND + (2._RKIND*stM%bff*Eff*Eff))
         g1   = (g1 + (2._RKIND*dc4f*Eff))*rexp
         g1   = 4._RKIND*stM%aff*g1
         CC   = CC + (g1*TEN_DYADPROD(Hff, Hff, nsd))

!        Sheet-sheet interaction stress + additional cross-fiber stress
         rexp = EXP(stM%bss * Ess * Ess)
         g2   = c4s*Ess*rexp
         g2   = g2 + (0.5_RKIND*dc4s/stM%bss)*(rexp - 1._RKIND)
         g2   = 2._RKIND*stM%ass*g2 + Tsa
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         S    = S + (g2*Hss)

!        Sheet-sheet interaction stiffness
         g2   = c4s*(1._RKIND + (2._RKIND*stM%bss*Ess*Ess))
         g2   = (g2 + (2._RKIND*dc4s*Ess))*rexp
         g2   = 4._RKIND*stM%ass*g2
         CC   = CC + (g2*TEN_DYADPROD(Hss, Hss, nsd))

      CASE DEFAULT
         err = "Undefined material constitutive model"
      END SELECT

!     Contribution from active strain
      IF (lDmn%ec%astrain) THEN
         S   = Ja * MATMUL(Fai, S)
         S   = MATMUL(S, TRANSPOSE(Fai))
         CCb = TEN_DYADPROD(Fai, Fai, nsd)
         CC  = TEN_DDOT_3424(CC, CCb, nsd)
         CC  = Ja * TEN_DDOT_2412(CCb, CC, nsd)
      END IF

!     Convert to Voigt Notation
      CALL CCTOVOIGT(CC, Dm)

      RETURN
      END SUBROUTINE GETPK2CC
!--------------------------------------------------------------------
      SUBROUTINE GETSVOLP(stM, J, p, pl)
      USE COMMOD
      IMPLICIT NONE
      TYPE(stModelType), INTENT(IN) :: stM
      REAL(KIND=RKIND), INTENT(IN) :: J
      REAL(KIND=RKIND), INTENT(INOUT) :: p, pl

      REAL(KIND=RKIND) Kp

      Kp = stM%Kpen
      SELECT CASE (stM%volType)
      CASE (stVol_Quad)
         p  = Kp*(J-1._RKIND)
         pl = Kp*(2._RKIND*J-1._RKIND)

      CASE (stVol_ST91)
         p  = 0.5_RKIND*Kp*(J-1._RKIND/J)
         pl = Kp*J

      CASE (stVol_M94)
         p  = Kp*(1._RKIND-1._RKIND/J)
         pl = Kp

      END SELECT

      END SUBROUTINE GETSVOLP
!####################################################################
!     Compute isochoric (deviatoric) component of 2nd Piola-Kirchhoff
!     stress and material stiffness tensors
      SUBROUTINE GETPK2CCdev(lDmn, F, nfd, fl, tmX, ya, S, Dm, Ja)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      INTEGER(KIND=IKIND), INTENT(IN) :: nfd
      REAL(KIND=RKIND), INTENT(IN) :: F(nsd,nsd), fl(nsd,nfd), tmX, ya
      REAL(KIND=RKIND), INTENT(OUT) :: S(nsd,nsd), Dm(nsymd,nsymd), Ja

      TYPE(stModelType) :: stM
      REAL(KIND=RKIND) :: nd, J, J2d, J4d, trE, Inv1, Inv2, Inv4, Inv6,
     2   Inv8, Tfa, IDm(nsd,nsd), C(nsd,nsd), E(nsd,nsd), Ci(nsd,nsd),
     3   Tsa, Sb(nsd,nsd), CCb(nsd,nsd,nsd,nsd), PP(nsd,nsd,nsd,nsd),
     4   CC(nsd,nsd,nsd,nsd)
      REAL(KIND=RKIND) :: r1, r2, g1, g2, g3, rexp
      ! Guccione !
      REAL(KIND=RKIND) :: QQ, Rm(nsd,nsd), Es(nsd,nsd), RmRm(nsd,nsd,6)
      ! HGO, HO !
      REAL(KIND=RKIND) :: Eff, Ess, Efs, kap, c4f, c4s, dc4f, dc4s,
     2   Hff(nsd,nsd), Hss(nsd,nsd), Hfs(nsd,nsd)
!     Active strain for electromechanics
      REAL(KIND=RKIND) :: Fe(nsd,nsd), Fa(nsd,nsd), Fai(nsd,nsd)

      S    = 0._RKIND
      Dm   = 0._RKIND

!     Some preliminaries
      stM  = lDmn%stM
      nd   = REAL(nsd, KIND=RKIND)

!     Fiber-reinforced stress
      CALL GET_FIB_STRESS(stM%Tf, Tfa)
      Tsa = Tfa*stM%Tf%eta_s

!     Electromechanics coupling - active stress
      IF (lDmn%ec%astress) THEN
         Tfa = Tfa + ya
         Tsa = Tsa + ya*lDmn%ec%eta_s
      END IF

!     Electromechanics coupling - active strain
      Fe   = F
      Fa   = MAT_ID(nsd)
      Fai  = Fa
      IF (lDmn%ec%astrain) THEN
         CALL GET_ACTV_DEFGRAD(lDmn%ec, nfd, fl, tmX, ya, Fa)
         Fai  = MAT_INV(Fa, nsd)
         Fe   = MATMUL(F, Fai)
      END IF

      Ja   = MAT_DET(Fa, nsd)
      J    = MAT_DET(Fe, nsd)
      J2d  = J**(-2._RKIND/nd)
      J4d  = J2d*J2d

      IDm  = MAT_ID(nsd)
      C    = MATMUL(TRANSPOSE(Fe), Fe)
      E    = 0.5_RKIND * (C - IDm)
      Ci   = MAT_INV(C, nsd)

      trE  = MAT_TRACE(E, nsd)
      Inv1 = J2d*MAT_TRACE(C,nsd)
      Inv2 = 0.5_RKIND*( Inv1*Inv1 - J4d*MAT_TRACE(MATMUL(C,C), nsd) )

!     Isochoric part of 2nd Piola-Kirchhoff and elasticity tensors
      SELECT CASE (stM%isoType)
!     NeoHookean model
      CASE (stIso_nHook)
         g1 = 2._RKIND * stM%C10
         Sb = g1*IDm

!        Fiber reinforcement/active stress
         Sb = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)

         r1 = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S  = J2d*Sb - r1*Ci
         CC = 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2        1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3      - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                   TEN_DYADPROD(S, Ci, nsd) )

!     Mooney-Rivlin model
      CASE (stIso_MR)
         g1  = 2._RKIND * (stM%C10 + Inv1*stM%C01)
         g2  = -2._RKIND * stM%C01
         Sb  = g1*IDm + g2*J2d*C

!        Fiber reinforcement/active stress
         Sb  = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)

         g1  = 4._RKIND*J4d* stM%C01
         CCb = g1 * (TEN_DYADPROD(IDm, IDm, nsd) - TEN_IDs(nsd))

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC + 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2              1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3            - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                         TEN_DYADPROD(S, Ci, nsd) )

!     HGO (Holzapfel-Gasser-Ogden) model for arteries with isochoric
!     invariants for the anisotropy terms (decoupled)
      CASE (stIso_HGO_d)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "HGO material model (2)"
         kap  = stM%kap
         Inv4 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = J2d*NORM(fl(:,2), MATMUL(C, fl(:,2)))

         Eff  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv4 - 1._RKIND
         Ess  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv6 - 1._RKIND

         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Hff  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hff
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Hss  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hss

         g1   = stM%C10
         g2   = stM%aff * Eff * EXP(stM%bff*Eff*Eff)
         g3   = stM%ass * Ess * EXP(stM%bss*Ess*Ess)
         Sb   = 2._RKIND*(g1*IDm + g2*Hff + g3*Hss)

!        Fiber reinforcement/active + additional cross-fiber stresses
         Sb   = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
     2             + Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd)

         g1   = stM%aff*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)*
     2      EXP(stM%bff*Eff*Eff)
         g2   = stM%ass*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)*
     2      EXP(stM%bss*Ess*Ess)
         g1   = 4._RKIND*J4d * g1
         g2   = 4._RKIND*J4d * g2

         CCb  = g1 * TEN_DYADPROD(Hff, Hff, nsd) +
     2          g2 * TEN_DYADPROD(Hss, Hss, nsd)

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC + 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2              1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3            - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                         TEN_DYADPROD(S, Ci, nsd) )

!     HGO (Holzapfel-Gasser-Ogden) model for arteries with full
!     invariants for the anisotropy terms (modified-anisotropy)
      CASE (stIso_HGO_ma)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "HGO material model (2)"
         kap  = stM%kap
         Inv1 = MAT_TRACE(C,nsd)
         Inv4 = NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = NORM(fl(:,2), MATMUL(C, fl(:,2)))

         Eff  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv4 - 1._RKIND
         Ess  = kap*Inv1 + (1._RKIND-3._RKIND*kap)*Inv6 - 1._RKIND

!        Isochoric contribution to stress and stiffness tensors
         Sb   = 2._RKIND*stM%C10*IDm
         r1   = J2d*MAT_DDOT(C, Sb, nsd) / nd

         S    = J2d*Sb - r1*Ci
         CC   = 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2          1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3        - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                          TEN_DYADPROD(S, Ci, nsd) )

!        Anisotropic contribution to stress and stiffness tensors
!        Fiber-Fiber interaction + additional fiber reinforcement
         rexp = EXP(stM%bff*Eff*Eff)
         g1   = 2._RKIND*stM%aff*Eff*rexp
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Hff  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hff
         S    = S + (g1*Hff) + (Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd))

         g1   = (1._RKIND + (2._RKIND*stM%bff*Eff*Eff))
         g1   = 4._RKIND*stM%aff*g1*rexp
         CC   = CC + (g1*TEN_DYADPROD(Hff, Hff, nsd))

!        Sheet-Sheet interaction + additional cross-fiber stress
         rexp = EXP(stM%bss*Ess*Ess)
         g2   = 2._RKIND*stM%ass*Ess*rexp
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Hss  = kap*IDm + (1._RKIND-3._RKIND*kap)*Hss
         S    = S + (g2*Hss) + (Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd))

         g2   = (1._RKIND + (2._RKIND*stM%bss*Ess*Ess))
         g2   = 4._RKIND*stM%ass*g2*rexp
         CC   = CC + (g2*TEN_DYADPROD(Hss, Hss, nsd))

!     Guccione (1995) transversely isotropic model
      CASE (stIso_Gucci)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Guccione material model (2)"
!        Compute isochoric component of E
         E = 0.5_RKIND * (J2d*C - Idm)

!        Transform into local orthogonal coordinate system
         Rm(:,1) = fl(:,1)
         Rm(:,2) = fl(:,2)
         Rm(:,3) = CROSS(fl)

!        Project E to local orthogocal coordinate system
         Es = MATMUL(E, Rm)
         Es = MATMUL(TRANSPOSE(Rm), Es)

         g1 = stM%bff
         g2 = stM%bss
         g3 = stM%bfs

         QQ = g1 * Es(1,1)*Es(1,1) +
     2        g2 *(Es(2,2)*Es(2,2) + Es(3,3)*Es(3,3)  +
     3             Es(2,3)*Es(2,3) + Es(3,2)*Es(3,2)) +
     4        g3 *(Es(1,2)*Es(1,2) + Es(2,1)*Es(2,1)  +
     5             Es(1,3)*Es(1,3) + Es(3,1)*Es(3,1))

         r2 = stM%C10 * EXP(QQ)

!        Fiber stiffness contribution := (dE*_ab / dE_IJ)
         RmRm(:,:,1) = MAT_DYADPROD(Rm(:,1), Rm(:,1), nsd)
         RmRm(:,:,2) = MAT_DYADPROD(Rm(:,2), Rm(:,2), nsd)
         RmRm(:,:,3) = MAT_DYADPROD(Rm(:,3), Rm(:,3), nsd)

         RmRm(:,:,4) = MAT_SYMMPROD(Rm(:,1), Rm(:,2), nsd)
         RmRm(:,:,5) = MAT_SYMMPROD(Rm(:,2), Rm(:,3), nsd)
         RmRm(:,:,6) = MAT_SYMMPROD(Rm(:,3), Rm(:,1), nsd)

         Sb = g1*Es(1,1)*RmRm(:,:,1) + g2*(Es(2,2)*RmRm(:,:,2) +
     2      Es(3,3)*RmRm(:,:,3) + 2._RKIND*Es(2,3)*RmRm(:,:,5)) +
     4      2._RKIND*g3*(Es(1,2)*RmRm(:,:,4) + Es(1,3)*RmRm(:,:,6))

         CCb = 2._RKIND*TEN_DYADPROD(Sb, Sb, nsd)
         Sb  = Sb * r2

!        Fiber reinforcement/active + addtional cross-fiber stresses
         Sb  = Sb + Tfa*MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
     2            + Tsa*MAT_DYADPROD(fl(:,2), fl(:,2), nsd)

         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         r2  = r2*J4d
         CCb = r2*(CCb + g1*TEN_DYADPROD(RmRm(:,:,1), RmRm(:,:,1), nsd)
     2       + g2*(TEN_DYADPROD(RmRm(:,:,2), RmRm(:,:,2), nsd)
     3       + TEN_DYADPROD(RmRm(:,:,3), RmRm(:,:,3), nsd)
     4       + TEN_DYADPROD(RmRm(:,:,5), RmRm(:,:,5), nsd)*2._RKIND)
     5       + 2._RKIND*g3*(TEN_DYADPROD(RmRm(:,:,4), RmRm(:,:,4), nsd)
     6       + TEN_DYADPROD(RmRm(:,:,6), RmRm(:,:,6), nsd)))

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC + 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2              1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3            - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                         TEN_DYADPROD(S, Ci, nsd) )

!     HO (Holzapfel-Ogden) model for myocardium with isoschoric
!     invariants for the anisotropy terms (decoupled)
      CASE (stIso_HO_d)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Holzapfel material model (2)"

!        Compute fiber-based invariants
         Inv4 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = J2d*NORM(fl(:,2), MATMUL(C, fl(:,2)))
         Inv8 = J2d*NORM(fl(:,1), MATMUL(C, fl(:,2)))

         Eff  = Inv4 - 1._RKIND
         Ess  = Inv6 - 1._RKIND
         Efs  = Inv8

!        Smoothed heaviside function
         c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
         c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

!        Approx. derivative of smoothed heaviside function
         dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
         dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))
c         dc4f = stM%khs /
c     2      (EXP(stM%khs*Eff) + EXP(-stM%khs*Eff) + 2.0_RKIND)
c         dc4s = stM%khs /
c     2      (EXP(stM%khs*Ess) + EXP(-stM%khs*Ess) + 2.0_RKIND)

!        Isotropic + fiber-sheet interaction stress
         g1   = stM%a * EXP(stM%b*(Inv1-3._RKIND))
         g2   = 2._RKIND * stM%afs * EXP(stM%bfs*Efs*Efs)
         Hfs  = MAT_SYMMPROD(fl(:,1), fl(:,2), nsd)
         Sb   = g1*IDm + g2*Efs*Hfs

!        Isotropic + fiber-sheet interaction stiffness
         g1   = g1 * 2._RKIND*J4d*stM%b
         g2   = g2 * 2._RKIND*J4d*(1._RKIND + 2._RKIND*stM%bfs*Efs*Efs)
         CCb  = g1 * TEN_DYADPROD(IDm, IDm, nsd) +
     2          g2 * TEN_DYADPROD(Hfs, Hfs, nsd)

!        Fiber-fiber interaction stress + additional reinforcement (Tfa)
         rexp = EXP(stM%bff * Eff * Eff)
         g1   = c4f*Eff*rexp
         g1   = g1 + (0.5_RKIND*dc4f/stM%bff)*(rexp - 1._RKIND)
         g1   = 2._RKIND*stM%aff*g1 + Tfa
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         Sb   = Sb + g1*Hff

!        Fiber-fiber interaction stiffness
         g1   = c4f*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)
         g1   = (g1 + 2._RKIND*dc4f*Eff)*rexp
         g1   = 4._RKIND*J4d*stM%aff*g1
         CCb  = CCb + g1*TEN_DYADPROD(Hff, Hff, nsd)

!        Sheet-sheet interaction stress + additional cross-fiber stress
         rexp = EXP(stM%bss * Ess * Ess)
         g2   = c4s*Ess*rexp
         g2   = g2 + (0.5_RKIND*dc4s/stM%bss)*(rexp - 1._RKIND)
         g2   = 2._RKIND*stM%ass*g2 + Tsa
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         Sb   = Sb + g2*Hss

!        Sheet-sheet interaction stiffness
         g2   = c4s*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)
         g2   = (g2 + 2._RKIND*dc4s*Ess)*rexp
         g2   = 4._RKIND*J4d*stM%ass*g2
         CCb  = CCb + g2*TEN_DYADPROD(Hss, Hss, nsd)

!        Isochoric 2nd-Piola-Kirchhoff stress and stiffness tensors
         r1  = J2d*MAT_DDOT(C, Sb, nsd) / nd
         S   = J2d*Sb - r1*Ci

         PP  = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC  = TEN_DDOT(CCb, PP, nsd)
         CC  = TEN_TRANSPOSE(CC, nsd)
         CC  = TEN_DDOT(PP, CC, nsd)
         CC  = CC + 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2              1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3            - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                         TEN_DYADPROD(S, Ci, nsd) )

!     HO (Holzapfel-Ogden) model for myocardium with full invariants
!     for the anisotropy terms (modified-anisotropy)
      CASE (stIso_HO_ma)
         IF (nfd .LT. 2) err = "Min fiber directions not defined for "//
     2      "Holzapfel material model (2)"

!        Compute fiber-based invariants
         Inv4 = NORM(fl(:,1), MATMUL(C, fl(:,1)))
         Inv6 = NORM(fl(:,2), MATMUL(C, fl(:,2)))
         Inv8 = NORM(fl(:,1), MATMUL(C, fl(:,2)))

         Eff  = Inv4 - 1._RKIND
         Ess  = Inv6 - 1._RKIND
         Efs  = Inv8

!        Smoothed heaviside function
         c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
         c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

!        Approx. derivative of smoothed heaviside function
         dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
         dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))
c         dc4f = stM%khs /
c     2      (EXP(stM%khs*Eff) + EXP(-stM%khs*Eff) + 2.0_RKIND)
c         dc4s = stM%khs /
c     2      (EXP(stM%khs*Ess) + EXP(-stM%khs*Ess) + 2.0_RKIND)

!        Isochoric stress and stiffness
         g1   = stM%a * EXP(stM%b*(Inv1-3._RKIND))
         Sb   = g1*IDm
         r1   = J2d*MAT_DDOT(C, Sb, nsd) / nd

         g1   = g1 * 2._RKIND*J4d*stM%b
         CCb  = g1 * TEN_DYADPROD(IDm, IDm, nsd)

!        Add isochoric stress and stiffness contribution
         S    = J2d*Sb - r1*Ci

         PP   = TEN_IDs(nsd) - (1._RKIND/nd) * TEN_DYADPROD(Ci, C, nsd)
         CC   = TEN_DDOT(CCb, PP, nsd)
         CC   = TEN_TRANSPOSE(CC, nsd)
         CC   = TEN_DDOT(PP, CC, nsd)
         CC   = CC + 2._RKIND*r1 * ( TEN_SYMMPROD(Ci, Ci, nsd) -
     2               1._RKIND/nd *   TEN_DYADPROD(Ci, Ci, nsd) )
     3             - 2._RKIND/nd * ( TEN_DYADPROD(Ci, S, nsd) +
     4                               TEN_DYADPROD(S, Ci, nsd) )

!        Now add anisotropic components
!        Fiber-sheet interaction terms
         g1   = 2._RKIND * stM%afs * EXP(stM%bfs*Efs*Efs)
         Hfs  = MAT_SYMMPROD(fl(:,1), fl(:,2), nsd)
         S    = S + (g1*Efs*Hfs)

         g1   = g1 * 2._RKIND*(1._RKIND + 2._RKIND*stM%bfs*Efs*Efs)
         CC   = CC + (g1*TEN_DYADPROD(Hfs, Hfs, nsd))

!        Fiber-fiber interaction stress + additional reinforcement (Tfa)
         rexp = EXP(stM%bff * Eff * Eff)
         g1   = c4f*Eff*rexp
         g1   = g1 + (0.5_RKIND*dc4f/stM%bff)*(rexp - 1._RKIND)
         g1   = (2._RKIND*stM%aff*g1) + Tfa
         Hff  = MAT_DYADPROD(fl(:,1), fl(:,1), nsd)
         S    = S + (g1*Hff)

!        Fiber-fiber interaction stiffness
         g1   = c4f*(1._RKIND + (2._RKIND*stM%bff*Eff*Eff))
         g1   = (g1 + (2._RKIND*dc4f*Eff))*rexp
         g1   = 4._RKIND*stM%aff*g1
         CC   = CC + (g1*TEN_DYADPROD(Hff, Hff, nsd))

!        Sheet-sheet interaction stress + additional cross-fiber stress
         rexp = EXP(stM%bss * Ess * Ess)
         g2   = c4s*Ess*rexp
         g2   = g2 + (0.5_RKIND*dc4s/stM%bss)*(rexp - 1._RKIND)
         g2   = 2._RKIND*stM%ass*g2 + Tsa
         Hss  = MAT_DYADPROD(fl(:,2), fl(:,2), nsd)
         S    = S + (g2*Hss)

!        Sheet-sheet interaction stiffness
         g2   = c4s*(1._RKIND + (2._RKIND*stM%bss*Ess*Ess))
         g2   = (g2 + (2._RKIND*dc4s*Ess))*rexp
         g2   = 4._RKIND*stM%ass*g2
         CC   = CC + (g2*TEN_DYADPROD(Hss, Hss, nsd))

      CASE DEFAULT
         err = "Undefined isochoric material constitutive model"
      END SELECT

!     Contribution from active strain
      IF (lDmn%ec%astrain) THEN
         S   = Ja * MATMUL(Fai, S)
         S   = MATMUL(S, TRANSPOSE(Fai))
         CCb = TEN_DYADPROD(Fai, Fai, nsd)
         CC  = TEN_DDOT_3424(CC, CCb, nsd)
         CC  = Ja * TEN_DDOT_2412(CCb, CC, nsd)
      END IF

!     Convert to Voigt Notation
      CALL CCTOVOIGT(CC, Dm)

      RETURN
      END SUBROUTINE GETPK2CCdev
!--------------------------------------------------------------------
!     Compute rho and beta depending on the Gibb's free-energy based
!     volumetric penalty model
      SUBROUTINE GVOLPEN(lDmn, p, ro, bt, dro, dbt, Ja)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(IN) :: p, Ja
      REAL(KIND=RKIND), INTENT(OUT) :: ro, bt, dro, dbt

      REAL(KIND=RKIND) :: Kp, r1, r2

      ro  = lDmn%prop(solid_density)/Ja
      bt  = 0._RKIND
      dbt = 0._RKIND
      dro = 0._RKIND

      Kp  = lDmn%stM%Kpen
      IF (ISZERO(Kp)) RETURN

      SELECT CASE (lDmn%stM%volType)
      CASE (stVol_Quad)
         r1  = 1._RKIND/(Kp - p)

         ro  = ro*Kp*r1
         bt  = r1
         dro = ro*r1
         dbt = r1*r1

      CASE (stVol_ST91)
         r1  = ro/Kp
         r2  = SQRT(p*p + Kp*Kp)

         ro  = r1*(p + r2)
         bt  = 1._RKIND/r2
         dro = ro*bt
         dbt = -bt*p/(p*p + Kp*Kp)

      CASE (stVol_M94)
         r1  = ro/Kp
         r2  = Kp + p

         ro  = r1*r2
         bt  = 1._RKIND/r2
         dro = r1
         dbt = -bt*bt

      CASE DEFAULT
         err = "Undefined volumetric material constitutive model"
      END SELECT

      RETURN
      END SUBROUTINE GVOLPEN
!--------------------------------------------------------------------
!     Compute stabilization parameters tauM and tauC
      SUBROUTINE GETTAU(lDmn, detF, Je, tauM, tauC)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      REAL(KIND=RKIND), INTENT(IN)  :: detF, Je
      REAL(KIND=RKIND), INTENT(OUT) :: tauM, tauC

      REAL(KIND=RKIND) :: ctM, ctC, he, rho0, Em, nu, mu, lam, c

      he   = 0.5_RKIND * Je**(1._RKIND/REAL(nsd, KIND=RKIND))
      rho0 = lDmn%prop(solid_density)
      Em   = lDmn%prop(elasticity_modulus)
      nu   = lDmn%prop(poisson_ratio)
      ctM  = lDmn%prop(ctau_M)
      ctC  = lDmn%prop(ctau_C)

      mu   = 0.5_RKIND*Em / (1._RKIND + nu)
      IF (ISZERO(nu-0.5_RKIND)) THEN
         c = SQRT(mu / rho0)
      ELSE
         lam = 2._RKIND*mu*nu / (1._RKIND-2._RKIND*nu)
         c = SQRT((lam + 2._RKIND*mu)/rho0)
      END IF

      tauM = ctM * (he/c) * (detF/rho0)
      tauC = ctC * (he*c) * (rho0/detF)

      RETURN
      END SUBROUTINE GETTAU
!####################################################################
!     Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
!     for incompressible shell elements
      SUBROUTINE GETPK2CC_SHLi(lDmn, nfd, fNa0, gg_0, gg_x, g33, Sml,
     2  Dml)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      INTEGER(KIND=IKIND), INTENT(IN) :: nfd
      REAL(KIND=RKIND), INTENT(IN) :: gg_0(2,2), gg_x(2,2), fNa0(2,nfd)
      REAL(KIND=RKIND), INTENT(OUT) :: g33, Sml(3), Dml(3,3)

      INTEGER a, b, iFn
      REAL(KIND=RKIND) :: Jg2i, I1, mu, gi_0(2,2), gi_x(2,2), S(2,2),
     2   CC(2,2,2,2), SN(2,2), CCN(2,2,2,2), Inv4, Inv6, Inv8, Eff, Ess,
     3   Efs, flM(2,2), c4f, c4s, dc4f, dc4s, d1, fl(2,2,nfd), Hfs(2,2),
     4   EI1, g1, g2
      TYPE(stModelType) :: stM

      Sml  = 0._RKIND
      Dml  = 0._RKIND

!     Some preliminaries
      stM  = lDmn%stM

!     Inverse of metric coefficients in shell continuum
      gi_0 = MAT_INV(gg_0, 2)
      gi_x = MAT_INV(gg_x, 2)

!     Ratio of inplane Jacobian determinants
      Jg2i = MAT_DET(gg_x, 2)
      IF (ISZERO(Jg2i)) err = " Divide by zero in-plane Jacobian"//
     2   " determinant"
      Jg2i = MAT_DET(gg_0, 2) / Jg2i

      I1 = 0._RKIND
      DO a=1, 2
        DO b=1, 2
            I1 = I1 + gi_0(a,b)*gg_x(a,b)
        END DO
      END DO

      SELECT CASE(stM%isoType)
      CASE (stIso_nHook)
         mu  = 2._RKIND * stM%C10
         S   = mu*(gi_0 - Jg2i*gi_x)

         CC  = 2._RKIND*mu*Jg2i*(TEN_DYADPROD(gi_x, gi_x, 2) +
     2                     TEN_SYMMPROD(gi_x, gi_x, 2))

      CASE (stIso_MR)
         SN  = (gi_0 - Jg2i*gi_x)
         CCN = 2._RKIND*Jg2i*(TEN_DYADPROD(gi_x, gi_x, 2) +
     2                     TEN_SYMMPROD(gi_x, gi_x, 2))

         S   = stM%C10*SN + stM%C01*Jg2i* (gi_0 - I1*gi_x)
     2       + stM%C01/Jg2i*gi_x
         CC  = (stM%C10 + stM%C01*I1) * CCN - 2._RKIND*stM%C01 *
     2         Jg2i * (TEN_DYADPROD(gi_0, gi_x, 2) +
     3         TEN_DYADPROD(gi_x, gi_0, 2)) + 2._RKIND*stM%C01 /
     4         Jg2i *(TEN_DYADPROD(gi_x, gi_0, 2) -
     5         TEN_SYMMPROD(gi_x, gi_x, 2))

!     HO (Holzapfel-Ogden) model for myocardium with full invariants
!     for the anisotropy terms (modified-anisotropy)
      CASE (stIso_HO_ma)
         IF (nfd .NE. 2) err = "Min fiber directions not defined for "//
     2      "Holzapfel material model (2)"

         DO iFn=1, nfd
            fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn)
            fl(1,2,iFn) = fNa0(1,iFn)*fNa0(2,iFn)
            fl(2,1,iFn) = fNa0(2,iFn)*fNa0(1,iFn)
            fl(2,2,iFn) = fNa0(2,iFn)*fNa0(2,iFn)
         END DO

!        Compute fiber-based invariants
         Inv4 = gg_x(1,1)*fl(1,1,1) + gg_x(1,2)*fl(1,2,1)
     2          + gg_x(2,1)*fl(2,1,1) + gg_x(2,2)*fl(2,2,1)
         Inv6 = gg_x(1,1)*fl(1,1,2) + gg_x(1,2)*fl(1,2,2)
     2          + gg_x(2,1)*fl(2,1,2) + gg_x(2,2)*fl(2,2,2)
         Inv8 = gg_x(1,1)*fNa0(1,1)*fNa0(1,2)
     2            + gg_x(1,2)*fNa0(1,1)*fNa0(2,2)
     3            + gg_x(2,1)*fNa0(2,1)*fNa0(1,2)
     4            + gg_x(2,2)*fNa0(2,1)*fNa0(2,2)

         Eff  = Inv4 - 1._RKIND
         Ess  = Inv6 - 1._RKIND
         Efs  = Inv8

!        Smoothed heaviside function
         c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
         c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

!        Approx. derivative of smoothed heaviside function
         dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
         dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))

!        Add isochoric stress and stiffness contribution
         EI1  = I1 + Jg2i - 3._RKIND
         SN   = (gi_0 - Jg2i*gi_x)
         CCN  = 2._RKIND*Jg2i*(TEN_DYADPROD(gi_x, gi_x, 2) +
     2                         TEN_SYMMPROD(gi_x, gi_x, 2))

         d1   = stM%a*EXP(stM%b*EI1)
         S    = d1*SN
         CC   = d1*(CCN + 2._RKIND*stM%b*TEN_DYADPROD(SN, SN, 2))

!        Anisotropic part
!        Fiber sheet
         Hfs  = MAT_SYMMPROD(fNa0(:,1), fNa0(:,2), 2)
         g1   = 2._RKIND*stM%afs*EXP(stM%bfs*Efs*Efs)
         S    = S + g1*Efs*Hfs

         g2   = g1*2._RKIND*(1._RKIND + 2._RKIND*Efs*Efs)
         CC   = CC + g2*TEN_DYADPROD(Hfs, Hfs, 2)

!        Fiber
         IF (Eff > 0._RKIND) THEN
            flM = fl(:,:,1)
            ! S  = S + 2._RKIND*stM%aff*Eff*flM
            ! CC = CC + 4._RKIND*stM%aff*TEN_DYADPROD(flM, flM, 2)
            g1 = 2._RKIND*stM%aff*EXP(stM%bff*Eff*Eff)
            S  = S  + g1*Eff*flM
            g2 = g1*2._RKIND*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)
            CC = CC + g2*TEN_DYADPROD(flM, flM, 2)
         END IF

!        Sheet
         IF (Ess > 0._RKIND) THEN
            flM = fl(:,:,2)
            g1  = 2._RKIND*stM%ass*EXP(stM%bss*Ess*Ess)
            S   = S + g1*Ess*flM
            g2  = g1*2._RKIND*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)
            CC  = CC + g2*TEN_DYADPROD(flM, flM, 2)
         END IF

!     Lee Sacks model for aorta with full invariants
!     for the anisotropy terms (modified-anisotropy)
      CASE (stIso_LS)
         SN = (gi_0 - Jg2i*gi_x)
         CCN= 2._RKIND*Jg2i*(TEN_DYADPROD(gi_x, gi_x, 2) +
     2                       TEN_SYMMPROD(gi_x, gi_x, 2))

         DO iFn=1, nfd
            fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn)
            fl(1,2,iFn) = fNa0(1,iFn)*fNa0(2,iFn)
            fl(2,1,iFn) = fNa0(2,iFn)*fNa0(1,iFn)
            fl(2,2,iFn) = fNa0(2,iFn)*fNa0(2,iFn)
         END DO

!        Compute fiber-based invariants
         Inv4 = gg_x(1,1)*fl(1,1,1) + gg_x(1,2)*fl(1,2,1)
     2        + gg_x(2,1)*fl(2,1,1) + gg_x(2,2)*fl(2,2,1)
         Eff  = Inv4 - 1._RKIND

!        Isotropic contribution
         EI1 = (I1 + Jg2i - 3._RKIND)
         S   = stM%a*SN + 2._RKIND*stM%a0*stM%mu0*stM%b1*EI1
     2       * EXP(stM%b1*EI1*EI1) * SN
         CC  = stM%a*CCN + 4._RKIND*stM%a0*stM%mu0*stM%b1
     2       * EXP(stM%b1*EI1*EI1)
     3       * ((1._RKIND + 2._RKIND*stM%b1*EI1*EI1)
     4       * TEN_DYADPROD(SN, SN, 2) + 0.5_RKIND*EI1*CCN)

!        Anisotropic contribution
         IF (Eff > 0._RKIND) THEN
            flM = fl(:,:,1)
            S   = S + 2._RKIND*stM%a0*(1._RKIND - stM%mu0)*stM%b2*Eff
     2          * EXP(stM%b2*Eff*Eff)*flM
            CC  = CC + 4._RKIND*stM%a0*(1._RKIND - stM%mu0)*stM%b2
     2          * (1._RKIND + 2._RKIND*stM%b2*Eff*Eff)
     3          * EXP(stM%b2*Eff*Eff) * TEN_DYADPROD(flM, flM, 2)

         END IF

      CASE DEFAULT
         err = "Undefined material constitutive model"
      END SELECT

      g33 = Jg2i

!     Convert to Voigt notation
      Sml(1) = S(1,1)
      Sml(2) = S(2,2)
      Sml(3) = S(1,2)

      Dml(1,1) = CC(1,1,1,1)
      Dml(1,2) = CC(1,1,2,2)
      Dml(1,3) = CC(1,1,1,2)

      Dml(2,2) = CC(2,2,2,2)
      Dml(2,3) = CC(2,2,1,2)

      Dml(3,3) = CC(1,2,1,2)

      Dml(2,1) = Dml(1,2)
      Dml(3,1) = Dml(1,3)
      Dml(3,2) = Dml(2,3)

      RETURN
      END SUBROUTINE GETPK2CC_SHLi
!--------------------------------------------------------------------
!     Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
!     for compressible shell elements
      SUBROUTINE GETPK2CC_SHLc(lDmn, nfd, fNa0, gg_0, gg_x, g33, Sml, 
     2  Dml)
      USE MATFUN
      USE COMMOD
      IMPLICIT NONE
      TYPE(dmnType), INTENT(IN) :: lDmn
      INTEGER(KIND=IKIND), INTENT(IN) :: nfd
      REAL(KIND=RKIND), INTENT(IN) :: gg_0(2,2), gg_x(2,2), fNa0(2,nfd)
      REAL(KIND=RKIND), INTENT(OUT) :: g33, Sml(3), Dml(3,3)

      INTEGER(KIND=IKIND), PARAMETER :: MAXITR = 20
      REAL(KIND=RKIND), PARAMETER :: ATOL = 1E-10

      INTEGER(KIND=IKIND) :: itr, i, j, k, l, iFn
      REAL(KIND=RKIND) :: Jg2, J2, J23, f13, f23, trC3, C33, kap, mu,
     2   pJ, plJ, gi_x(2,2), gi_0(3,3), Ci(3,3), S(3,3), CC(3,3,3,3),
     3   fl(2,2,nfd), C1, C2, J43, Gi4AS(3,3,3,3), I2, I2ij(3,3),
     4   I2ijkl(3,3,3,3), Cikl(3,3,3,3), Inv4, Inv6, Inv8, Eff, Ess,
     5   Efs, c4f, c4s, dc4f, dc4s, EI1, d1, g1, g2, Hfs(3,3), 
     6   SN(3,3), CCN(3,3,3,3), flM(3,3)
      TYPE(stModelType) :: stM

      Sml  = 0._RKIND
      Dml  = 0._RKIND

!     Initialize tensor operations
      CALL TEN_INIT(3)

!     Some preliminaries
      stM  = lDmn%stM
      kap  = stM%Kpen
      mu   = 2._RKIND * stM%C10
      f13  = 1._RKIND / 3._RKIND
      f23  = 2._RKIND / 3._RKIND

!     Inverse of metric coefficients in shell continuum
      gi_x = MAT_INV(gg_x, 2)

      gi_0 = 0._RKIND
      gi_0(1:2,1:2) = MAT_INV(gg_0, 2)
      gi_0(3,3) = 1._RKIND

!     Ratio of inplane Jacobian determinant squared
      Jg2 = MAT_DET(gg_x, 2) / MAT_DET(gg_0, 2)

!     Begin Newton iterations to satisfy plane-stress condition.
!     The objective is to find C33 that satisfies S33 = 0.
      itr = 0
      C33 = 1._RKIND
      DO
         itr  = itr + 1

!        Trace (C)
         trC3 = (gg_x(1,1)*gi_0(1,1) + gg_x(1,2)*gi_0(1,2)
     2        +  gg_x(2,1)*gi_0(2,1) + gg_x(2,2)*gi_0(2,2) + C33)*f13

!        Jacobian-related quantities
         J2  = Jg2*C33
         J23 = J2**(-f13)

!        Inverse of curvilinear Cauchy-Green deformation tensor
         Ci(:,:) = 0._RKIND
         Ci(3,3) = 1._RKIND/C33
         Ci(1:2,1:2) = gi_x(:,:)

!        Contribution from dilational penalty terms to S and CC
         pJ  = 0.5_RKIND*kap*(J2 - 1._RKIND)
         plJ = kap*J2

         SELECT CASE (stM%isoType)
         CASE (stIso_nHook)
!           2nd Piola Kirchhoff stress
            S  = mu*J23*(gi_0 - trC3*Ci) + pJ*Ci

!           Elasticity tensor
            CC = (mu*J23*f23*trC3 + plJ)*TEN_DYADPROD(Ci, Ci, 3)
     2         + (mu*J23*trC3 - pJ)*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)
     3         - f23*mu*J23*(TEN_DYADPROD(gi_0, Ci, 3) +
     4                       TEN_DYADPROD(Ci, gi_0, 3))

         CASE (stIso_MR)
!           2nd Piola Kirchhoff stress
            C1 = stM%C10
            C2 = stM%C01
            J43 = J2**(-f23)
            I2ijkl = TEN_DYADPROD(gi_0, gi_0, 3) - 
     2                              TEN_SYMMPROD(gi_0, gi_0, 3)          
            I2ij = TEN_MDDOT(I2ijkl, Ci, 3)
            Gi4AS = TEN_ASYMPROD12(gi_0, gi_0, 3)
            I2  = MAT_DDOT(Ci, TEN_MDDOT(Gi4AS, Ci, 3), 3)
            ! I2  = MAT_DDOT(Ci, I2ij, 3)
            Cikl = - TEN_SYMMPROD(Ci, Ci, 3)

            S  = C1*J23*(gi_0 - trC3*Ci) + pJ*Ci
            S  = S + C2*J43*(-f23*I2*Ci + I2ij)

!           Elasticity tensor
            CC = (C1*J23*f23*trC3 + plJ)*TEN_DYADPROD(Ci, Ci, 3)
     2         + (C1*J23*trC3 - pJ)*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)
     3         - f23*C1*J23*(TEN_DYADPROD(gi_0, Ci, 3) +
     4                       TEN_DYADPROD(Ci, gi_0, 3))
            CC = CC + 2._RKIND*f23*C2*J43*( 
     2               TEN_DYADPROD((f23*I2*Ci-I2ij), Ci, 3) - I2*Cikl
     3                - TEN_DYADPROD(Ci, I2ij, 3) + I2ijkl)
                     
         CASE (stIso_HO_ma)
            IF (nfd .NE. 2) err = "Min fiber directions not defined"//  
     2          "for Holzapfel material model (2)"

            DO iFn=1, nfd
                fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn)
                fl(1,2,iFn) = fNa0(1,iFn)*fNa0(2,iFn)
                fl(2,1,iFn) = fNa0(2,iFn)*fNa0(1,iFn)
                fl(2,2,iFn) = fNa0(2,iFn)*fNa0(2,iFn)
            END DO

    !        Compute fiber-based invariants
            Inv4 = gg_x(1,1)*fl(1,1,1) + gg_x(1,2)*fl(1,2,1)
     2              + gg_x(2,1)*fl(2,1,1) + gg_x(2,2)*fl(2,2,1)
            Inv6 = gg_x(1,1)*fl(1,1,2) + gg_x(1,2)*fl(1,2,2)
     2            + gg_x(2,1)*fl(2,1,2) + gg_x(2,2)*fl(2,2,2)
            Inv8 = gg_x(1,1)*fNa0(1,1)*fNa0(1,2)
     2              + gg_x(1,2)*fNa0(1,1)*fNa0(2,2)
     3              + gg_x(2,1)*fNa0(2,1)*fNa0(1,2)
     4              + gg_x(2,2)*fNa0(2,1)*fNa0(2,2)

            Eff  = Inv4 - 1._RKIND
            Ess  = Inv6 - 1._RKIND
            Efs  = Inv8

    !        Smoothed heaviside function
            c4f  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Eff))
            c4s  = 1._RKIND / (1._RKIND + EXP(-stM%khs*Ess))

    !        Approx. derivative of smoothed heaviside function
            dc4f = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Eff))
            dc4s = 0.25_RKIND*stM%khs*EXP(-stM%khs*ABS(Ess))

    !        Add isochoric stress and stiffness contribution
            ! EI1  = I1 + Jg2i - 3._RKIND
            d1 = stM%a*J23*EXP(3._RKIND*stM%b*(trC3*J23 - 1._RKIND))
            SN = (gi_0 - trC3*Ci)
            S  = d1*SN + pJ*Ci
            CC = f23*trC3*TEN_DYADPROD(Ci, Ci, 3)
     2         + trC3*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)
     3         - f23*(TEN_DYADPROD(gi_0, Ci, 3)
     4         + TEN_DYADPROD(Ci, gi_0, 3))
     5         + 2._RKIND*stM%b*J23*TEN_DYADPROD(SN, SN, 3)
            CC = d1*CC + plJ*TEN_DYADPROD(Ci, Ci, 3)
     2         - pJ*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)     

    !        Anisotropic part
    !        Fiber sheet
            Hfs = 0._RKIND
            Hfs(1:2,1:2)  = MAT_SYMMPROD(fNa0(:,1), fNa0(:,2), 2)
            g1   = 2._RKIND*stM%afs*EXP(stM%bfs*Efs*Efs)
            S    = S + g1*Efs*Hfs

            g2   = g1*2._RKIND*(1._RKIND + 2._RKIND*Efs*Efs)
            CC   = CC + g2*TEN_DYADPROD(Hfs, Hfs, 3)

    !        Fiber
            flM = 0._RKIND
            IF (Eff > 0._RKIND) THEN
                flM(1:2,1:2) = fl(:,:,1)
                ! S  = S + 2._RKIND*stM%aff*Eff*flM
                ! CC = CC + 4._RKIND*stM%aff*TEN_DYADPROD(flM, flM, 2)
                g1 = 2._RKIND*stM%aff*EXP(stM%bff*Eff*Eff)
                S  = S  + g1*Eff*flM
                g2 = g1*2._RKIND*(1._RKIND + 2._RKIND*stM%bff*Eff*Eff)
                CC = CC + g2*TEN_DYADPROD(flM, flM, 3)
            END IF

    !        Sheet
            IF (Ess > 0._RKIND) THEN
                flM(1:2,1:2) = fl(:,:,2)
                g1  = 2._RKIND*stM%ass*EXP(stM%bss*Ess*Ess)
                S   = S + g1*Ess*flM
                g2  = g1*2._RKIND*(1._RKIND + 2._RKIND*stM%bss*Ess*Ess)
                CC  = CC + g2*TEN_DYADPROD(flM, flM, 3)
            END IF

         CASE (stIso_LS)
            DO iFn=1, nfd
                fl(1,1,iFn) = fNa0(1,iFn)*fNa0(1,iFn)
                fl(1,2,iFn) = fNa0(1,iFn)*fNa0(2,iFn)
                fl(2,1,iFn) = fNa0(2,iFn)*fNa0(1,iFn)
                fl(2,2,iFn) = fNa0(2,iFn)*fNa0(2,iFn)
            END DO

!           Compute fiber-based invariants
            Inv4 = gg_x(1,1)*fl(1,1,1) + gg_x(1,2)*fl(1,2,1)
     2          + gg_x(2,1)*fl(2,1,1) + gg_x(2,2)*fl(2,2,1)
            Eff  = Inv4 - 1._RKIND

!           Isotropic contribution
            d1 = 2._RKIND*stM%a0*stM%mu0*stM%b1*J23
     2          * EXP(3._RKIND*stM%b1*(trC3*J23 - 1._RKIND))
            SN = (gi_0 - trC3*Ci)
            CCN = f23*trC3*TEN_DYADPROD(Ci, Ci, 3)
     2         + trC3*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)
     3         - f23*(TEN_DYADPROD(gi_0, Ci, 3)
     4         + TEN_DYADPROD(Ci, gi_0, 3))
            S  = (stM%a + d1*(3._RKIND*(trC3*J23 - 1._RKIND)))*SN
     2         + pJ*Ci
            CC = (1._RKIND + 18._RKIND*stM%b1*(trC3*J23 - 1._RKIND)
     2         * (trC3*J23 - 1._RKIND))*J23*TEN_DYADPROD(SN, SN, 3)
     3         + 3._RKIND*(trC3*J23 - 1._RKIND)*CCN
            CC = stM%a*CCN + 2._RKIND*d1*CC 
     2         + plJ*TEN_DYADPROD(Ci, Ci, 3)
     3         - pJ*2._RKIND*TEN_SYMMPROD(Ci, Ci, 3)     

!           Anisotropic contribution
            flM = 0._RKIND
            IF (Eff > 0._RKIND) THEN
                flM(1:2,1:2) = fl(:,:,1)
                S   = S + 2._RKIND*stM%a0*(1._RKIND - stM%mu0)
     2              *stM%b2*Eff * EXP(stM%b2*Eff*Eff)*flM
                CC  = CC + 4._RKIND*stM%a0*(1._RKIND - stM%mu0)
     2              * stM%b2 * (1._RKIND + 2._RKIND*stM%b2*Eff*Eff)
     3              * EXP(stM%b2*Eff*Eff) * TEN_DYADPROD(flM, flM, 3)
            END IF
            

         CASE DEFAULT
            err = "Undefined material constitutive model"

         END SELECT

         IF (ABS(S(3,3)) .LE. ATOL) EXIT
         IF (itr .GT. MAXITR) THEN
            wrn = " Failed to converge plane-stress condition"
            EXIT
         END IF

         C33 = C33 - (2._RKIND*S(3,3)/CC(3,3,3,3))
      END DO

      g33 = C33

!     Statically condense CC
      DO i=1, 2
         DO j=1, 2
            DO k=1, 2
               DO l=1, 2
                  C33 = CC(i,j,3,3)*CC(3,3,k,l)/CC(3,3,3,3)
                  CC(i,j,k,l) = CC(i,j,k,l) - C33
               END DO
            END DO
         END DO
      END DO

      g33 = C33

!     Convert the in-plane components to Voigt notation
      Sml(1) = S(1,1)
      Sml(2) = S(2,2)
      Sml(3) = S(1,2)

      Dml(1,1) = CC(1,1,1,1)
      Dml(1,2) = CC(1,1,2,2)
      Dml(1,3) = CC(1,1,1,2)

      Dml(2,2) = CC(2,2,2,2)
      Dml(2,3) = CC(2,2,1,2)

      Dml(3,3) = CC(1,2,1,2)

      Dml(2,1) = Dml(1,2)
      Dml(3,1) = Dml(1,3)
      Dml(3,2) = Dml(2,3)

      RETURN
      END SUBROUTINE GETPK2CC_SHLc
!####################################################################
!     Convert elasticity tensor to Voigt notation
      SUBROUTINE CCTOVOIGT(CC, Dm)
      USE COMMOD, ONLY : IKIND, RKIND, nsd, nsymd
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: CC(nsd,nsd,nsd,nsd)
      REAL(KIND=RKIND), INTENT(INOUT) :: Dm(nsymd,nsymd)

      INTEGER(KIND=IKIND) :: i, j

      IF (nsd .EQ. 3) THEN
         Dm(1,1) = CC(1,1,1,1)
         Dm(1,2) = CC(1,1,2,2)
         Dm(1,3) = CC(1,1,3,3)
         Dm(1,4) = CC(1,1,1,2)
         Dm(1,5) = CC(1,1,2,3)
         Dm(1,6) = CC(1,1,3,1)

         Dm(2,2) = CC(2,2,2,2)
         Dm(2,3) = CC(2,2,3,3)
         Dm(2,4) = CC(2,2,1,2)
         Dm(2,5) = CC(2,2,2,3)
         Dm(2,6) = CC(2,2,3,1)

         Dm(3,3) = CC(3,3,3,3)
         Dm(3,4) = CC(3,3,1,2)
         Dm(3,5) = CC(3,3,2,3)
         Dm(3,6) = CC(3,3,3,1)

         Dm(4,4) = CC(1,2,1,2)
         Dm(4,5) = CC(1,2,2,3)
         Dm(4,6) = CC(1,2,3,1)

         Dm(5,5) = CC(2,3,2,3)
         Dm(5,6) = CC(2,3,3,1)

         Dm(6,6) = CC(3,1,3,1)

         DO i=2, 6
            DO j=1, i-1
               Dm(i,j) = Dm(j,i)
            END DO
         END DO

      ELSE IF (nsd .EQ. 2) THEN
         Dm(1,1) = CC(1,1,1,1)
         Dm(1,2) = CC(1,1,2,2)
         Dm(1,3) = CC(1,1,1,2)

         Dm(2,2) = CC(2,2,2,2)
         Dm(2,3) = CC(2,2,1,2)

         Dm(3,3) = CC(1,2,1,2)

         Dm(2,1) = Dm(1,2)
         Dm(3,1) = Dm(1,3)
         Dm(3,2) = Dm(2,3)

      END IF

      RETURN
      END SUBROUTINE CCTOVOIGT
!####################################################################
!     Compute additional fiber-reinforcement stress
      SUBROUTINE GET_FIB_STRESS(Tfl, g)
      USE COMMOD
      IMPLICIT NONE
      TYPE(fibStrsType), INTENT(IN) :: Tfl
      REAL(KIND=RKIND), INTENT(OUT) :: g

      REAL(KIND=RKIND) rtmp

      g = 0._RKIND
      IF (BTEST(Tfl%fType, bType_std)) THEN
         g = Tfl%g
      ELSE IF (BTEST(Tfl%fType, bType_ustd)) THEN
         CALL IFFT(Tfl%gt, g, rtmp)
      END IF

      RETURN
      END SUBROUTINE GET_FIB_STRESS
!####################################################################
!     Compute active component of deformation gradient tensor for
!     electromechanics coupling based on active strain formulation
      SUBROUTINE GET_ACTV_DEFGRAD(ec, nfd, fl, lam, yf, Fa)
      USE MATFUN
      USE UTILMOD
      USE COMMOD, ONLY : asnType_tiso, asnType_ortho, asnType_hetortho,
     2   nsd, eccModelType
      IMPLICIT NONE
      TYPE(eccModelType), INTENT(IN) :: ec
      INTEGER(KIND=IKIND), INTENT(IN) :: nfd
      REAL(KIND=RKIND), INTENT(IN) :: fl(nsd,nfd), lam, yf
      REAL(KIND=RKIND), INTENT(OUT) :: Fa(nsd,nsd)

      REAL(KIND=RKIND) :: ys, yn, f(nsd), s(nsd), n(nsd), Hf(nsd,nsd),
     2   Hs(nsd,nsd), Hn(nsd,nsd)

      f  = fl(:,1)
      s  = fl(:,2)
      n  = CROSS(fl)

      Hf = MAT_DYADPROD(f, f, nsd)

!     Transversely isotropic activation
      IF (ec%asnType .EQ. asnType_tiso) THEN
         yn = (1._RKIND/SQRT(1._RKIND+yf)) - 1._RKIND

!     Orthotropic activation
      ELSE IF (ec%asnType .EQ. asnType_ortho) THEN
         yn = ec%k*yf

!     Transmurally heteregenous orthotropic activation
!     lam: transmural coordinate (=0, endo; =1, epi)
      ELSE IF (ec%asnType .EQ. asnType_hetortho) THEN
         yn = (1._RKIND - lam)*ec%k*yf
     2      + lam*(1._RKIND/SQRT(1._RKIND + yf) - 1._RKIND)

      END IF

      ys = 1._RKIND/((1._RKIND+yf)*(1._RKIND+yn)) - 1._RKIND
      Hs = MAT_DYADPROD(s, s, nsd)
      Hn = MAT_DYADPROD(n, n, nsd)

      Fa = (1._RKIND+yf)*Hf + (1._RKIND+ys)*Hs + (1._RKIND+yn)*Hn
!      Fa = MAT_ID(nsd) + (yf*Hf) + (ys*Hs) + (yn*Hn)

      RETURN
      END SUBROUTINE GET_ACTV_DEFGRAD
!####################################################################

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
!     This routine computes body force for the current equation and
!     assembles it to the residue
!
!--------------------------------------------------------------------

      SUBROUTINE SETBF(Dg)
      USE COMMOD
      IMPLICIT NONE
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) iBf, iM

      Bf(:,:) = 0._RKIND
      DO iBf=1, eq(cEq)%nBf
         iM = eq(cEq)%bf(iBf)%iM
         CALL SETBFL(eq(cEq)%bf(iBf), msh(iM), Dg)
      END DO

      RETURN
      END SUBROUTINE SETBF
!--------------------------------------------------------------------
      SUBROUTINE SETBFL(lBf, lM, Dg)
      USE COMMOD
      IMPLICIT NONE
      TYPE(bfType), INTENT(IN) :: lBf
      TYPE(mshType), INTENT(IN) :: lM
      REAL(KIND=RKIND), INTENT(IN) :: Dg(tDof,tnNo)

      INTEGER(KIND=IKIND) a, e, Ac, idof, nNo, eNoN
      REAL(KIND=RKIND) rtmp

      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)
      REAL(KIND=RKIND), ALLOCATABLE :: f(:), bfl(:,:), bfg(:,:),
     2   xl(:,:), dl(:,:)

      nNo  = lM%nNo
      idof = lBf%dof
      eNoN = lM%eNoN

      ALLOCATE(f(idof))

      IF (BTEST(lBf%bType, bfType_std)) THEN
         f(:) = lBf%b(:)

      ELSE IF (BTEST(lBf%bType,bfType_ustd)) THEN
         CALL IFFT(lBf%bt, f, rtmp)

      ELSE IF (BTEST(lBf%bType, bfType_gen)) THEN
         ALLOCATE(bfl(idof,nNo), xl(idof,nNo))
         CALL IGBC(lBf%bm, bfl, xl)
         DEALLOCATE(xl)

      END IF

      ALLOCATE(bfg(idof,tnNo))
      bfg = 0._RKIND

      IF (BTEST(lBf%bType,bfType_gen)) THEN
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            bfg(:,Ac) = bfl(:,a)
         END DO
         DEALLOCATE(bfl)
      ELSE IF (BTEST(lBf%bType,bfType_spl)) THEN
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            bfg(:,Ac) = lBf%bx(:,a)
         END DO
      ELSE
         DO a=1, lM%nNo
            Ac = lM%gN(a)
            bfg(:,Ac) = f(:)
         END DO
      END IF

      DEALLOCATE(f)

!     Assemble pressure/traction load (shells/CMM initialization) to
!     residue. For general body force (vector), assemble later with
!     other volumetric forces
      IF (BTEST(lBf%bType, bfType_vol)) THEN
         DO a=1, nNo
            Ac = lM%gN(a)
            Bf(:,Ac) = bfg(:,Ac)
         END DO

      ELSE
         ALLOCATE(bfl(idof,eNoN), xl(nsd,eNoN), dl(tDof,eNoN),
     2      ptr(eNoN))
         DO e=1, lM%nEl
            DO a=1, eNoN
               Ac = lM%IEN(a,e)
               ptr(a)   = Ac
               xl(:,a)  = x(:,Ac)
               dl(:,a)  = Dg(:,Ac)
               bfl(:,a) = bfg(:,Ac)
            END DO
            CALL BFCONSTRUCT(lM, e, eNoN, idof, xl, dl, bfl, ptr)
         END DO
         DEALLOCATE(bfl, xl, dl, ptr)
      END IF

      DEALLOCATE(bfg)

      RETURN
      END SUBROUTINE SETBFL
!####################################################################
!     This subroutine is reached only for shell follower pressure loads
!     or applying initialization pressure for CMM method. nsd must be
!     equal to 3
      SUBROUTINE BFCONSTRUCT(lM, e, eNoN, idof, xl, dl, bfl, ptr)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      TYPE(mshType), INTENT(IN) :: lM
      INTEGER(KIND=IKIND), INTENT(IN) :: e, eNoN, idof, ptr(eNoN)
      REAL(KIND=RKIND), INTENT(IN) :: dl(tDof,eNoN), bfl(idof,eNoN)
      REAL(KIND=RKIND), INTENT(INOUT) :: xl(nsd,eNoN)

      INTEGER(KIND=IKIND) :: g, cPhys
      REAL(KIND=RKIND) :: w

      REAL(KIND=RKIND), ALLOCATABLE :: N(:), Nx(:,:), lR(:,:), lK(:,:,:)

      IF (nsd .NE. 3) err = "Unexpected encounter. Correction needed "//
     2   "in BFCONSTRUCT"

      ALLOCATE(N(eNoN), Nx(2,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

      cDmn  = DOMAIN(lM, cEq, e)
      cPhys = eq(cEq)%dmn(cDmn)%phys

!     Updating the shape functions, if neccessary
      IF (lM%eType .EQ. eType_NRB) CALL NRBNNX(lM, e)

!     Setting intial values
      lR = 0._RKIND
      lK = 0._RKIND

      DO g=1, lM%nG
         w  = lM%w(g)
         N  = lM%N(:,g)
         Nx = lM%Nx(:,:,g)
         SELECT CASE (cPhys)
         CASE (phys_shell)
            CALL SHELLFP(eNoN, w, N, Nx, dl, xl, bfl, lR, lK)

         CASE (phys_CMM)
            CALL BCMMi(eNoN, idof, w, N, Nx, xl, bfl, lR)

         CASE DEFAULT
            err = "Undefined phys in BFCONSTRUCT"
         END SELECT
      END DO

!     Now doing the assembly part
#ifdef WITH_TRILINOS
      IF (eq(cEq)%assmTLS) THEN
         CALL TRILINOS_DOASSEM(eNoN, ptr, lK, lR)
      ELSE
#endif
         CALL DOASSEM(eNoN, ptr, lK, lR)
#ifdef WITH_TRILINOS
      END IF
#endif

      RETURN
      END SUBROUTINE BFCONSTRUCT
!####################################################################

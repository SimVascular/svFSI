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
!     Subroutines related to initializing linear solver arrays and
!     function calls to svFSILS and Trilinos solver library
!
!--------------------------------------------------------------------

      SUBROUTINE LSALLOC(lEq)
      USE COMMOD
      IMPLICIT NONE
      TYPE(eqType), INTENT(IN) :: lEq

      IF (ALLOCATED(R)) DEALLOCATE(R)
      ALLOCATE (R(dof,tnNo))
      R(:,:) = 0._RKIND

      IF (.NOT.lEq%assmTLS) THEN
         IF (ALLOCATED(Val)) DEALLOCATE(Val)
         ALLOCATE (Val(dof*dof,lhs%nnz))
         Val(:,:) = 0._RKIND
      END IF

#ifdef WITH_TRILINOS
      IF (lEq%ls%LS_Packg .EQ. lSPackg_TRILINOS) THEN
         IF (ALLOCATED(tls%W)) THEN
            DEALLOCATE(tls%W, tls%R)
            CALL TRILINOS_LHS_FREE()
         END IF
         ALLOCATE(tls%W(dof,tnNo), tls%R(dof,tnNo))
         CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo, lhs%nnz,
     2      tls%ltg, ltg, rowPtr, colPtr, dof)
      END IF
#endif

      RETURN
      END SUBROUTINE LSALLOC
!####################################################################
      SUBROUTINE LSSOLVE(lEq, incL, res)
      USE COMMOD
#ifdef WITH_HYPRE
      USE HYPREMOD
#endif
      IMPLICIT NONE
      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER(KIND=IKIND), INTENT(IN) :: incL(nFacesLS)
      REAL(KIND=RKIND), INTENT(IN) :: res(nFacesLS)
      REAL(KIND=RKIND), ALLOCATABLE :: Wr(:,:), Wc(:,:)

#ifdef WITH_TRILINOS
      INTEGER(KIND=IKIND) a

      IF (lEq%assmTLS) THEN
         lEq%FSILS%RI%suc = .FALSE.
         CALL TRILINOS_SOLVE(tls%R, tls%W, lEq%FSILS%RI%fNorm,
     2      lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr, lEq%FSILS%RI%callD,
     3      lEq%FSILS%RI%dB, lEq%FSILS%RI%suc, lEq%ls%LS_type,
     4      lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr, lEq%FSILS%RI%sD,
     5      lEq%ls%PREC_Type, lEq%assmTLS)

      ELSE IF(lEq%useTLS) THEN
         CALL TRILINOS_GLOBAL_SOLVE(Val, R, tls%R, tls%W,
     2      lEq%FSILS%RI%fNorm, lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr,
     3      lEq%FSILS%RI%callD, lEq%FSILS%RI%dB, lEq%FSILS%RI%suc,
     4      lEq%ls%LS_type, lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr,
     5      lEq%FSILS%RI%sD, lEq%ls%PREC_Type)

      ELSE
#endif
         CALL FSILS_SOLVE(lhs, lEq%FSILS, dof, R, Val,
     2      lEq%ls%PREC_Type, incL=incL, res=res)
#ifdef WITH_TRILINOS
      END IF

      IF (lEq%ls%LS_Packg .EQ. lSPackg_TRILINOS) THEN
         DO a=1, tnNo
            R(:,a) = tls%R(:,lhs%map(a))
         END DO
         R = Wc*R
         DEALLOCATE(Wr,Wc)
      END IF
#endif

      RETURN
      END SUBROUTINE LSSOLVE
!--------------------------------------------------------------------
      SUBROUTINE INIT_DIR_AND_COUPNEU_BC(incL, res)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE
      INTEGER(KIND=IKIND), INTENT(IN) :: incL(lhs%nFaces)
      REAL(KIND=RKIND), INTENT(IN) :: res(lhs%nFaces)

      REAL(KIND=RKIND), ALLOCATABLE :: v(:,:)

      INTEGER(KIND=IKIND) a, i, Ac, faIn, faDof
      LOGICAL flag, isCoupledBC

      IF (lhs%nFaces .NE. 0) THEN
         lhs%face%incFlag = .TRUE.
         DO faIn=1, lhs%nFaces
            IF (incL(faIn) .EQ. 0)  lhs%face(faIn)%incFlag = .FALSE.
         END DO
         DO faIn=1, lhs%nFaces
            lhs%face(faIn)%coupledFlag = .FALSE.
            IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
            flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
            IF (flag .AND. res(faIn).NE.0._RKIND) THEN
               lhs%face(faIn)%res = res(faIn)
               lhs%face(faIn)%coupledFlag = .TRUE.
            END IF
         END DO
      END IF

      IF (eq(cEq)%ls%LS_Packg .EQ. lSPackg_TRILINOS) THEN
            tls%W = 1._RKIND
      ! DO faIn=1, lhs%nFaces
      !    IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
      !    faDof = MIN(lhs%face(faIn)%dof,dof)
      !    IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
      !       DO a=1, lhs%face(faIn)%nNo
      !          Ac = lhs%face(faIn)%glob(a)
      !          DO i=1, faDof
      !             tls%W(i,Ac) = tls%W(i,Ac) * lhs%face(faIn)%val(i,a)
      !          END DO
      !       END DO
      !    END IF
      ! END DO
      END IF

      ALLOCATE(v(dof,tnNo))
      v = 0._RKIND
      isCoupledBC = .FALSE.
      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            isCoupledBC = .TRUE.
            faDof = MIN(lhs%face(faIn)%dof,dof)
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, faDof
                  v(i,Ac) = v(i,Ac) +
     2               SQRT(ABS(res(faIn)))*lhs%face(faIn)%val(i,a)
               END DO
            END DO
         END IF
      END DO

#ifdef WITH_TRILINOS
      IF (eq(cEq)%ls%LS_Packg .EQ. lSPackg_TRILINOS) 
     2   CALL TRILINOS_BC_CREATE(v, isCoupledBC)
#endif
      DEALLOCATE(v)

      END SUBROUTINE INIT_DIR_AND_COUPNEU_BC
!####################################################################

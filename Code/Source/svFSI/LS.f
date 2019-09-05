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

      SUBROUTINE LSALLOC()
      USE COMMOD
      IMPLICIT NONE

      IF (ALLOCATED(R)) THEN
         IF (SIZE(R,1) .NE. dof) THEN
            DEALLOCATE(R)
            ALLOCATE (R(dof,tnNo))
            IF (.NOT.useTrilinosAssemAndLS) THEN
               DEALLOCATE(Val)
               ALLOCATE (Val(dof*dof,lhs%nnz))
            END IF

#ifdef WITH_TRILINOS
            IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
               DEALLOCATE(tls%W, tls%R)
               ALLOCATE(tls%W(dof,tnNo), tls%R(dof,tnNo))
               CALL TRILINOS_LHS_FREE()
               CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo,
     2            lhs%nnz, tls%ltg, ltg, rowPtr, colPtr, dof)
            END IF
#endif
         END IF
      ELSE
         ALLOCATE (R(dof,tnNo))
         IF (.NOT.useTrilinosAssemAndLS) THEN
            ALLOCATE(Val(dof*dof,lhs%nnz))
         END IF
#ifdef WITH_TRILINOS
         IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
            CALL TRILINOS_LHS_CREATE(gtnNo, lhs%mynNo, tnNo, lhs%nnz,
     2         tls%ltg, ltg, rowPtr, colPtr, dof)
            ALLOCATE(tls%W(dof,tnNo), tls%R(dof,tnNo))
         END IF
#endif
      END IF

      IF (.NOT. useTrilinosAssemAndLS) Val(:,:) = 0D0
      R(:,:) = 0D0

      RETURN
      END SUBROUTINE LSALLOC
!####################################################################
      SUBROUTINE LSSOLVE(lEq, incL, res)
      USE COMMOD
      IMPLICIT NONE

      TYPE(eqType), INTENT(INOUT) :: lEq
      INTEGER, INTENT(IN) :: incL(nFacesLS)
      REAL(KIND=8), INTENT(IN) :: res(nFacesLS)

#ifdef WITH_TRILINOS
      INTEGER a

      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
         CALL INIT_DIR_AND_COUPNEU_BC(incL, res)
      END IF

      IF (useTrilinosAssemAndLS) THEN
         lEq%FSILS%RI%suc = .FALSE.
         CALL TRILINOS_SOLVE(tls%R, tls%W, lEq%FSILS%RI%fNorm,
     2      lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr, lEq%FSILS%RI%callD,
     3      lEq%FSILS%RI%dB, lEq%FSILS%RI%suc, lEq%ls%LS_type,
     4      lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr, lEq%FSILS%RI%sD,
     5      lEq%ls%PREC_Type, useTrilinosAssemAndLS)

      ELSE IF(useTrilinosLS) THEN
         CALL TRILINOS_GLOBAL_SOLVE(Val, R, tls%R, tls%W,
     2      lEq%FSILS%RI%fNorm, lEq%FSILS%RI%iNorm, lEq%FSILS%RI%itr,
     3      lEq%FSILS%RI%callD, lEq%FSILS%RI%dB, lEq%FSILS%RI%suc,
     4      lEq%ls%LS_type, lEq%FSILS%RI%reltol, lEq%FSILS%RI%mItr,
     5      lEq%FSILS%RI%sD, lEq%ls%PREC_Type)

      ELSE
#endif
         CALL FSILS_SOLVE(lhs, lEq%FSILS, dof, R, Val, incL=incL,
     2      res=res)
#ifdef WITH_TRILINOS
      END IF

      IF (useTrilinosLS .OR. useTrilinosAssemAndLS) THEN
         DO a=1, tnNo
            R(:,a) = tls%R(:,lhs%map(a))
         END DO
      END IF
#endif

      RETURN
      END SUBROUTINE LSSOLVE
!--------------------------------------------------------------------
      SUBROUTINE INIT_DIR_AND_COUPNEU_BC(incL, res)
      USE COMMOD
      USE ALLFUN
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: incL(lhs%nFaces)
      REAL(KIND=8), INTENT(IN) :: res(lhs%nFaces)

      REAL(KIND=8), ALLOCATABLE :: v(:,:)

      INTEGER a, i, Ac, faIn, faDof
      LOGICAL flag, isCoupledBC

      IF (lhs%nFaces .NE. 0) THEN
         lhs%face%incFlag = .TRUE.
         DO faIn=1, lhs%nFaces
            IF (incL(faIn) .EQ. 0) lhs%face(faIn)%incFlag = .FALSE.
         END DO
         DO faIn=1, lhs%nFaces
            lhs%face(faIn)%coupledFlag = .FALSE.
            IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
            flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
            IF (flag .AND. res(faIn).NE.0D0) THEN
               lhs%face(faIn)%res = res(faIn)
               lhs%face(faIn)%coupledFlag = .TRUE.
            END IF
         END DO
      END IF

      tls%W = 1D0
      DO faIn=1, lhs%nFaces
         IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
         faDof = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, faDof
                  tls%W(i,Ac) = tls%W(i,Ac) * lhs%face(faIn)%val(i,a)
               END DO
            END DO
         END IF
      END DO

      ALLOCATE(v(dof,tnNo))
      v = 0D0
      isCoupledBC = .FALSE.
      DO faIn=1, lhs%nFaces
         IF (lhs%face(faIn)%coupledFlag) THEN
            isCoupledBC = .TRUE.
            faDof = MIN(lhs%face(faIn)%dof,dof)
            DO a=1, lhs%face(faIn)%nNo
               Ac = lhs%face(faIn)%glob(a)
               DO i=1, faDof
                  v(i,Ac) = v(i,Ac) +
     2               SQRT(res(faIn))*lhs%face(faIn)%val(i,a)
               END DO
            END DO
         END IF
      END DO

#ifdef WITH_TRILINOS
      CALL TRILINOS_BC_CREATE(v, isCoupledBC)
#endif
      DEALLOCATE(v)

      END SUBROUTINE INIT_DIR_AND_COUPNEU_BC
!####################################################################

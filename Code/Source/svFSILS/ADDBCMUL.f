!--------------------------------------------------------------------
!     Created by Mahdi Esmaily Moghadam
!     contact memt63@gmail.com for reporting the bugs.
!--------------------------------------------------------------------
!
!     UC Copyright Notice
!     This software is Copyright ©2012 The Regents of the University of
!     California. All Rights Reserved.
!
!     Permission to copy and modify this software and its documentation
!     for educational, research and non-profit purposes, without fee,
!     and without a written agreement is hereby granted, provided that
!     the above copyright notice, this paragraph and the following three
!     paragraphs appear in all copies.
!
!     Permission to make commercial use of this software may be obtained
!     by contacting:
!     Technology Transfer Office
!     9500 Gilman Drive, Mail Code 0910
!     University of California
!     La Jolla, CA 92093-0910
!     (858) 534-5815
!     invent@ucsd.edu
!
!     This software program and documentation are copyrighted by The
!     Regents of the University of California. The software program and
!     documentation are supplied "as is", without any accompanying
!     services from The Regents. The Regents does not warrant that the
!     operation of the program will be uninterrupted or error-free. The
!     end-user understands that the program was developed for research
!     purposes and is advised not to rely exclusively on the program for
!     any reason.
!
!     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
!     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
!     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
!     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
!     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
!     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
!     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
!     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
!     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
!
!--------------------------------------------------------------------
!     The contribution of coupled BCs is added to the matrix-vector
!     product operation. Depending on the type of operation (adding the
!     contribution or compution the PC contribution) different
!     coefficients are used.
!     The resistance is stored in lhs%face(faIn)%res
!     The current matrix vector product is in Y. The vector to be 
!     multiplied and added is in X. The matrix is represented by res.
!        
!--------------------------------------------------------------------

      SUBROUTINE ADDBCMUL(lhs, op_Type, dof, X, Y)
      INCLUDE "FSILS_STD.h"
      TYPE(FSILS_lhsType), INTENT(INOUT) :: lhs
      INTEGER(KIND=LSIP), INTENT(IN) :: op_type, dof
      REAL(KIND=LSRP), INTENT(IN) :: X(dof, lhs%nNo)
      REAL(KIND=LSRP), INTENT(INOUT) :: Y(dof, lhs%nNo)

      INTEGER(KIND=LSIP) faIn, i, a, Ac, nsd
      REAL(KIND=LSRP) S, FSILS_DOTV

      REAL(KIND=LSRP), ALLOCATABLE :: v(:,:), coef(:)

      ALLOCATE(coef(lhs%nFaces), v(dof,lhs%nNo))

!     Here is where the res(istance) value is "added" to the stiffness
!     matrix (Moghadam et al. 2013 eq. 27).
!     res is transfered to a variable coef
!     See FSILS_STRUCT.h for lhs%face variable descriptions
      IF (op_Type .EQ. BCOP_TYPE_ADD) THEN
         coef = lhs%face%res
      ELSE IF(op_Type .EQ. BCOP_TYPE_PRE) THEN
         coef = -lhs%face%res/(1._LSRP + (lhs%face%res*lhs%face%nS))
      ELSE
         PRINT *, "FSILS: op_Type is not defined"
         STOP "FSILS: FATAL ERROR"
      END IF

      DO faIn=1, lhs%nFaces
         nsd = MIN(lhs%face(faIn)%dof,dof)
         IF (lhs%face(faIn)%coupledFlag) THEN
            IF (lhs%face(faIn)%sharedFlag) THEN ! if face is shared between procs
               v = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     v(i,Ac) = lhs%face(faIn)%valM(i,a)
                  END DO
               END DO
               S = coef(faIn)*FSILS_DOTV(dof,lhs%mynNo, lhs%commu, v, X)
               DO a=1, lhs%face(faIn)%nNo
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S
                  END DO
               END DO
            ELSE
!              What is valM(i,a)? See PRECOND.F for where is it set
!              lhs%face(faIn)%valM(i,a) =                            &
!     &               lhs%face(faIn)%val(i,a)*W(i,Ac)
!              val(i,a) contains the integrals of Na * n_i over a surface
!              which is precisely the integral in Moghadam et al. 2013, eq. 27
!              Note, val should contain integrals over deformed config surface
!              if follower pressure load is used (BNEUFOLWP)
!              I believe valM is val, scaled by some factor to account
!              for preconditioning
               S = 0._LSRP
               DO a=1, lhs%face(faIn)%nNo ! Loop over nodes on face
                  Ac = lhs%face(faIn)%glob(a) ! Get global node number
                  DO i=1, nsd
                     S = S + lhs%face(faIn)%valM(i,a)*X(i,Ac)
                  END DO
               END DO
!              Multiply S by the resistance (equal to coef, see above)
               S = coef(faIn)*S
!              Add S times second integral to the current matrix-vector product Y
               DO a=1, lhs%face(faIn)%nNo
!                 If the coupled surface is virtually capped to compute flow rate
!                 then the right integral should be over the capped surface, while
!                 the left integral should be over the uncapped surface. We
!                 can satisfy this by skipping this addition if the face is virtual
                  IF (lhs%face(faIn)%virtual) CYCLE
                  Ac = lhs%face(faIn)%glob(a)
                  DO i=1, nsd
                     Y(i,Ac) = Y(i,Ac) + lhs%face(faIn)%valM(i,a)*S
                  END DO
               END DO
            END IF
         END IF
      END DO

      RETURN
      END SUBROUTINE ADDBCMUL
!####################################################################

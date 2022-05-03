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
!     The data structures used in FSILS.
!--------------------------------------------------------------------

      INCLUDE "FSILS_TYPEDEF.h"

!     Some definitions
      INTEGER(KIND=LSIP), PARAMETER :: LS_TYPE_CG = 798,
     &   LS_TYPE_GMRES = 797, LS_TYPE_NS = 796, LS_TYPE_BICGS = 795

      INTEGER(KIND=LSIP), PARAMETER :: PRECOND_FSILS = 701,
     &   PRECOND_RCS = 709

      INTEGER(KIND=LSIP), PARAMETER :: BC_TYPE_Dir = 0, BC_TYPE_Neu = 1

      INTEGER(KIND=LSIP), PARAMETER :: BCOP_TYPE_ADD = 0,
     &   BCOP_TYPE_PRE = 1

!     Communication structure
      TYPE FSILS_commuType
         SEQUENCE
!        Free of created          (USE)
         LOGICAL :: foC = .FALSE.
!        If this the master       (USE)
         LOGICAL masF
!        Master ID                (USE)
         INTEGER(KIND=LSIP) master
!        ID of this proc.         (USE)
         INTEGER(KIND=LSIP) task
!        Task in FORTRAN indexing (USE)
         INTEGER(KIND=LSIP) tF
!        Total number of tasks    (USE)
         INTEGER(KIND=LSIP) nTasks
!        MPI communicator         (IN)
         INTEGER(KIND=LSIP) comm
!        Only for data alignment
         INTEGER(KIND=LSIP) reserved
      END TYPE FSILS_commuType

!     LHS matrix related data
      TYPE FSILS_faceType
         SEQUENCE
!        Free or created                (USE)
         LOGICAL :: foC = .FALSE.
!        Neu: P/Q coupling              (USE)
         LOGICAL coupledFlag
!        Neu: shared between proces     (USE)
         LOGICAL :: sharedFlag = .FALSE.
!        Included in the computations   (IN)
         LOGICAL incFlag
!        Number of nodes                (IN)
         INTEGER(KIND=LSIP) :: nNo = 0
!        Degrees of freedom for val     (IN)
         INTEGER(KIND=LSIP) dof
!        Dir/Neu                        (IN)
         INTEGER(KIND=LSIP) :: bGrp = BC_TYPE_Dir
!        Only for data alignment
         INTEGER(KIND=LSIP) reserved
!        Global node number             (IN)
         INTEGER(KIND=LSIP), ALLOCATABLE :: glob(:)
!        ||Sai||**2._LSRP                   (USE)
         REAL(KIND=LSRP) nS
!        Neu: P = res*Q                 (IN)
         REAL(KIND=LSRP) :: res = 0._LSRP    ! Resistance value for Neumann face
!        nodal Sai for Neu              (IN)
         REAL(KIND=LSRP), ALLOCATABLE :: val(:,:)
!        Neu W*Sai                      (TMP)
         REAL(KIND=LSRP), ALLOCATABLE :: valM(:,:)
      END TYPE FSILS_faceType

!     All following are in (USE)
      TYPE FSILS_cSType
         SEQUENCE
!        The processor to communicate with
         INTEGER(KIND=LSIP) iP
!        Number of data to be commu
         INTEGER(KIND=LSIP) n
!        Pointer to the data for commu
         INTEGER(KIND=LSIP), ALLOCATABLE :: ptr(:)
      END TYPE FSILS_cSType

      TYPE FSILS_lhsType
         SEQUENCE
!        Free of created                     (USE)
         LOGICAL :: foC = .FALSE.
!        Global number of nodes              (IN)
         INTEGER(KIND=LSIP) :: gnNo = 0
!        Number of nodes                     (IN)
         INTEGER(KIND=LSIP) :: nNo = 0
!        Number of non-zero in lhs           (IN)
         INTEGER(KIND=LSIP) :: nnz = 0
!        Number of faces                     (IN)
         INTEGER(KIND=LSIP) :: nFaces = 0
!        nNo of this proc                    (USE)
         INTEGER(KIND=LSIP) mynNo
!        nNo of shared with lower proc       (USE)
         INTEGER(KIND=LSIP) shnNo
!        Number of communication requests    (USE)
         INTEGER(KIND=LSIP) :: nReq = 0
!        Column pointer                      (USE)
         INTEGER(KIND=LSIP), ALLOCATABLE :: colPtr(:)
!        Row pointer                         (USE)
         INTEGER(KIND=LSIP), ALLOCATABLE :: rowPtr(:,:)
!        Diagonal pointer                    (USE)
         INTEGER(KIND=LSIP), ALLOCATABLE :: diagPtr(:)
!        Mapping of nodes                    (USE)
         INTEGER(KIND=LSIP), ALLOCATABLE :: map(:)
         TYPE(FSILS_commuType) commu
         TYPE(FSILS_cSType), ALLOCATABLE :: cS(:)
         TYPE(FSILS_faceType), ALLOCATABLE :: face(:)
      END TYPE FSILS_lhsType

!     LS related structures
      TYPE FSILS_subLsType
         SEQUENCE
!        Successful solving            (OUT)
         LOGICAL suc
!        Maximum iteration             (IN)
         INTEGER(KIND=LSIP) mItr
!        Space dimension               (IN)
         INTEGER(KIND=LSIP) sD
!        Number of iteration           (OUT)
         INTEGER(KIND=LSIP) itr
!        Number of Ax multiply         (OUT)
         INTEGER(KIND=LSIP) cM
!        Number of |x| norms           (OUT)
         INTEGER(KIND=LSIP) cN
!        Number of <x.y> dot products  (OUT)
         INTEGER(KIND=LSIP) cD
!        Only for data alignment       (-)
         INTEGER(KIND=LSIP) reserve
!        Absolute tolerance            (IN)
         REAL(KIND=LSRP) absTol
!        Relative tolerance            (IN)
         REAL(KIND=LSRP) relTol
!        Initial norm of residual      (OUT)
         REAL(KIND=LSRP) iNorm
!        Final norm of residual        (OUT)
         REAL(KIND=LSRP) fNorm
!        Res. rduction in last itr.    (OUT)
         REAL(KIND=LSRP) dB
!        Calling duration              (OUT)
         REAL(KIND=LSRP) callD
      END TYPE FSILS_subLsType

      TYPE FSILS_lsType
         SEQUENCE
!        Free of created             (USE)
         LOGICAL :: foC = .FALSE.
!        Which one of LS             (IN)
         INTEGER(KIND=LSIP) LS_type
!        Contribution of mom. res.   (OUT)
         INTEGER(KIND=LSIP) Resm
!        Contribution of cont. res.  (OUT)
         INTEGER(KIND=LSIP) Resc
         TYPE(FSILS_subLsType) GM
         TYPE(FSILS_subLsType) CG
         TYPE(FSILS_subLsType) RI
      END TYPE FSILS_lsType

!####################################################################

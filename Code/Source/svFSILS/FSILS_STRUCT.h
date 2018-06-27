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

!     Some definitions
      INTEGER, PARAMETER :: LS_TYPE_CG=798, LS_TYPE_GMRES=797,
     &   LS_TYPE_NS=796, LS_TYPE_BICGS=795

      INTEGER, PARAMETER :: BC_TYPE_Dir = 0, BC_TYPE_Neu = 1
      INTEGER, PARAMETER :: BCOP_TYPE_ADD = 0, BCOP_TYPE_PRE = 1


!     Communication structure
      TYPE FSILS_commuType
         SEQUENCE
!        Free of created          (USE)
         LOGICAL :: foC = .FALSE.
!        If this the master       (USE)
         LOGICAL masF
!        Master ID                (USE)
         INTEGER master
!        ID of this proc.         (USE)
         INTEGER task
!        Task in FORTRAN indexing (USE)
         INTEGER tF
!        Total number of tasks    (USE)
         INTEGER nTasks
!        MPI communicator         (IN)
         INTEGER comm
!        Only for data alignment
         INTEGER reserved
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
         INTEGER :: nNo = 0
!        Degrees of freedom for val     (IN)
         INTEGER dof
!        Dir/Neu                        (IN)
         INTEGER :: bGrp = BC_TYPE_Dir
!        Only for data alignment
         INTEGER reserved
!        Global node number             (IN)
         INTEGER, ALLOCATABLE :: glob(:)
!        ||Sai||**2D0                   (USE)
         REAL(KIND=8) nS
!        Neu: P = res*Q                 (IN)
         REAL(KIND=8) :: res = 0D0
!        nodal Sai for Neu              (IN)
         REAL(KIND=8), ALLOCATABLE :: val(:,:)
!        Neu W*Sai                      (TMP)
         REAL(KIND=8), ALLOCATABLE :: valM(:,:)
      END TYPE FSILS_faceType

!     All following are in (USE)
      TYPE FSILS_cSType
         SEQUENCE
!        The processor to communicate with
         INTEGER iP
!        Number of data to be commu
         INTEGER n
!        Pointer to the data for commu
         INTEGER, ALLOCATABLE :: ptr(:)
      END TYPE FSILS_cSType

      TYPE FSILS_lhsType
         SEQUENCE
!        Free of created                     (USE)
         LOGICAL :: foC = .FALSE.
!        Global number of nodes              (IN)
         INTEGER :: gnNo = 0
!        Number of nodes                     (IN)
         INTEGER :: nNo = 0
!        Number of non-zero in lhs           (IN)
         INTEGER :: nnz = 0
!        Number of faces                     (IN)
         INTEGER :: nFaces = 0
!        nNo of this proc                    (USE)
         INTEGER mynNo
!        nNo of shared with lower proc       (USE)
         INTEGER shnNo
!        Number of communication requests    (USE)
         INTEGER :: nReq = 0
!        Column pointer                      (USE)
         INTEGER, ALLOCATABLE :: colPtr(:)
!        Row pointer                         (USE)
         INTEGER, ALLOCATABLE :: rowPtr(:,:)
!        Diagonal pointer                    (USE)
         INTEGER, ALLOCATABLE :: diagPtr(:)
!        Mapping of nodes                    (USE)
         INTEGER, ALLOCATABLE :: map(:)
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
         INTEGER mItr
!        Space dimension               (IN)
         INTEGER sD
!        Number of iteration           (OUT)
         INTEGER itr
!        Number of Ax multiply         (OUT)
         INTEGER cM
!        Number of |x| norms           (OUT)
         INTEGER cN
!        Number of <x.y> dot products  (OUT)
         INTEGER cD
!        Only for data alignment       (-)
         INTEGER reserve
!        Absolute tolerance            (IN)
         REAL(KIND=8) absTol
!        Relative tolerance            (IN)
         REAL(KIND=8) relTol
!        Initial norm of residual      (OUT)
         REAL(KIND=8) iNorm
!        Final norm of residual        (OUT)
         REAL(KIND=8) fNorm
!        Res. rduction in last itr.    (OUT)
         REAL(KIND=8) dB
!        Calling duration              (OUT)
         REAL(KIND=8) callD
      END TYPE FSILS_subLsType

      TYPE FSILS_lsType
         SEQUENCE
!        Free of created             (USE)
         LOGICAL :: foC = .FALSE.
!        Which one of LS             (IN)
         INTEGER LS_type
!        Contribution of mom. res.   (OUT)
         INTEGER Resm
!        Contribution of cont. res.  (OUT)
         INTEGER Resc
         TYPE(FSILS_subLsType) GM
         TYPE(FSILS_subLsType) CG
         TYPE(FSILS_subLsType) RI
      END TYPE FSILS_lsType

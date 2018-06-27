/*
 *--------------------------------------------------------------------
 *  Created by Mahdi Esmaily Moghadam
 *  contact memt63@gmail.com for reporting the bugs.
 *--------------------------------------------------------------------
 *  UC Copyright Notice
 *  This software is Copyright ©2012 The Regents of the University of
 *  California. All Rights Reserved.

 *  Permission to copy and modify this software and its documentation
 *  for educational, research and non-profit purposes, without fee,
 *  and without a written agreement is hereby granted, provided that
 *  the above copyright notice, this paragraph and the following three
 *  paragraphs appear in all copies.
 *
 *  Permission to make commercial use of this software may be obtained
 *  by contacting:
 *  Technology Transfer Office
 *  9500 Gilman Drive, Mail Code 0910
 *  University of California
 *  La Jolla, CA 92093-0910
 *  (858) 534-5815
 *  invent@ucsd.edu
 *
 *  This software program and documentation are copyrighted by The
 *  Regents of the University of California. The software program and
 *  documentation are supplied "as is", without any accompanying
 *  services from The Regents. The Regents does not warrant that the
 *  operation of the program will be uninterrupted or error-free. The
 *  end-user understands that the program was developed for research
 *  purposes and is advised not to rely exclusively on the program for
 *  any reason.
 *
 *  IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
 *  PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
 *  DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
 *  SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
 *  CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *  THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
 *  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 *  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
 *  SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
 *  UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
 *  MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 *
 *--------------------------------------------------------------------
 *  The data structures used in FSILS.
 *--------------------------------------------------------------------
 */
//     Some defenitions
#define LS_TYPE_CG 798
#define LS_TYPE_GMRES 797
#define LS_TYPE_NS 796
#define BC_TYPE_Dir 0
#define BC_TYPE_Neu 1
#define BCOP_TYPE_ADD 0
#define BCOP_TYPE_PRE 1

//     Communication structure
typedef struct {
   int foC;// = false;              // Free of created          (USE)
   int masF;                        // If this the master       (USE)
   int master;                      // Master ID                (USE)
   int task;                        // ID of this proc.         (USE)
   int tF;                          // Task in FORTRAN indexing (USE)
   int nTasks;                      // Total number of tasks    (USE)
   MPI_Aint comm;                   // MPI communicator         (IN)
} FSILS_commuType;
/*
//     LHS matrix related data
typedef struct {
   int foC;// = false;        // Free or created                (USE)
   int coupledFlag;        // Neu: P/Q coupling              (USE)
   int sharedFlag;// = false; // Neu: shared between proces     (USE)
   int incFlag;                  // Included in the computations   (IN)
   int nNo;// = 0;            // Number of nodes                (IN)
   int dof;                // Degrees of freedom for val     (IN)
   int bGrp;// = BC_TYPE_Dir; // Dir/Neu                        (IN)
   int reserved;                 // Only for data alignment
   int *glob;              // Global node number             (IN)
   double nS;              // ||Sai||**2D0                   (USE)
   double res;// = 0D0;       // Neu: P = res*Q                 (IN)
   double **val;           // nodal Sai for Neu              (IN)
   double **valM;          // Neu W*Sai                      (TMP)
} FSILS_faceType;

typedef struct {
   int ptr;       // Pointer to start of data for commu (only 2 proc shared points)
   int n;         // Number of data to be commu  (only 2 proc shared points)
   int tag;       // Commu tag
   int req;       // Commu req
   int nBl;       // Number of blocks for commu  (for 3 < proc shared points)
   int reserved;  // Only for data alignment
   int *blPtr;    // Pointer to beggining of each block (for 3 < proc shared points)
   int *blN;      // Length of each block (for 3 < proc shared points)
} FSILS_cSType;

typedef struct {
   int foC;// = false;        // Free or created                (USE)
   int gnNo;// = 0;           // Global number of nodes      (IN)
   int nNo;// = 0;            // Number of nodes             (IN)
   int nnz;// = 0;            // Number of non-zero in lhs   (IN)
   int nFaces;// = 0;         // Number of faces             (IN)
   int mynNo;              // nNo of this proc            (USE)
   int *colPtr;            // Column pointer              (USE)
   int **rowPtr;           // Row pointer                 (USE)
   int *diagPt;            // Diagonal pointer            (USE)
   int *map;               // Mapping of nodes            (USE)
   FSILS_commuType commu;
   FSILS_cSType *cS;
   FSILS_faceType *face;
} FSILS_lhsType;
*/
//     LS related structures
typedef struct {
   int suc;       // Successful solving          (OUT)
   int mItr;      // Maximum iteration           (IN)
   int sD;        // Space dimension             (IN)
   int itr;       // Number of iteration         (OUT)
   int cM;        // Number of Ax multiply       (OUT)
   int cN;        // Number of |x| norms         (OUT)
   int cD;        // Number of <x.y> dot products(OUT)
   int reserve;   // Only for data alignment     (-)
   double absTol; // Absolute tolerance          (IN)
   double relTol; // Relative tolerance          (IN)
   double iNorm;  // Initial norm of residual    (OUT)
   double fNorm;  // Final norm of residual      (OUT)
   double dB;     // Res. rduction in last itr.  (OUT)
   double callD;  // Calling duration            (OUT)
} FSILS_subLsType;

typedef struct {
   int foC;// = false;   // Free of created             (USE)
   int LS_type;            // Which one of LS             (IN)
   int Resm;               // Contribution of mom. res.   (OUT)
   int Resc;               // Contribution of cont. res.  (OUT)
   FSILS_subLsType GM;
   FSILS_subLsType CG;
   FSILS_subLsType RI;
} FSILS_lsType;

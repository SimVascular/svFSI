/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*!
  \file    trilinos_linear_solver.cpp
  \brief   wrap Trilinos solver functions
*/

#include "trilinos_linear_solver.h"
#define NOOUTPUT

// --- Define global Trilinos variables to be used in below functions ---------

/// Unique block map consisting of nodes owned by each processor
Epetra_BlockMap *Trilinos::blockMap;

/// Global block force vector
Epetra_FEVector *Trilinos::F;

/// Global block stiffness matrix
Epetra_FEVbrMatrix *Trilinos::K;

/// Solution vector consisting of unique nodes owned by the processor
Epetra_Vector *Trilinos::X;

/// Solution vector with nodes owned by processor followed by its ghost nodes
Epetra_Vector *Trilinos::ghostX;

/// Import ghostMap into blockMap to create ghost map
Epetra_Import *Trilinos::Importer;

/// Contribution from coupled neumann boundary conditions
Epetra_FEVector *Trilinos::bdryVec;

Epetra_FECrsGraph *Trilinos::K_graph;

ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;

Ifpack_Preconditioner* ifpackPrec = NULL;

// --- Define global variables to be shared amongst the below functions -------

/// Nodal degrees of freedom
int dof;

/// Total number of nodes including the ghost nodes
int ghostAndLocalNodes;

/// Nodes owned by processor
int localNodes;

/// Converts local proc column indices to global indices to be inserted
std::vector<int> globalColInd;

/// Converts local indices to global indices in unsorted ghost node order
std::vector<int> localToGlobalUnsorted;

/// Stores number of nonzeros per row for the topology
std::vector<int> nnzPerRow;

std::vector<int> localToGlobalSorted;

int timecount = 0;

bool coupledBC;

// ----------------------------------------------------------------------------
/**
 * Define the matrix vector multiplication operation to do at each iteration
 * (K + v*v') *x = K*x + v*(v'*x), v = bdryVec
 *
 * uses efficient vectorized operation rather than explicitly forming
 * the rank 1 outer product matrix
 *
 * \param x vector to be applied on the operator
 * \param y result of sprase matrix vector multiplication
 */
int TrilinosMatVec::Apply(const Epetra_MultiVector &x,
        Epetra_MultiVector &y) const
{
  //store initial matrix vector product result in Y
  Trilinos::K->Apply(x, y); //K*x
  if (coupledBC)
  {
    //now need to add on bdry term v*(v'*x)
    double *dot = new double[1]; //only 1 multivector to store result in it

    //compute dot product (v'*x)
    Trilinos::bdryVec->Dot(x,dot);

    //Y = 1*Y + dot*v daxpy operation
    y.Update(*dot, //scalar for v
             *Trilinos::bdryVec, //FE_Vector
             1.0); ///scalar for Y

    delete[] dot;
  }
  return 0;
}

// ----------------------------------------------------------------------------
/**
 * make the graph global and optimize storage on it
 * this function creates the LHS structure topology and RHS vector based on the
 * block map
 * \param numGlobalNodes total/global number of nodes each node can have dof
 *                       for the spatial dim
 * \param numLocalNodes  number of nodes owned by this proc in block coordiantes
 * \param numGhostAndLocalNodes number of nodes owned and shared by this
 *                              processor includes ghost nodes
 * \param nnz            number of nonzeros in LHS matrix of calling proc
 * \param ltgSorted      integer pointer of size numLocalNodes returns global
 *                       node number of local node a to give blockrow
 * \param ltgUnsorted    unsorted/not-reordered version
 * \param rowPtr         CSR row ptr of size numLocalNodes + 1 to block rows
 * \param colInd         CSR column indices ptr (size nnz points) to block rows
 * \param Dof            size of each block element to give dim of each block

 */
void trilinos_lhs_create_(unsigned &numGlobalNodes, unsigned &numLocalNodes,
        unsigned &numGhostAndLocalNodes, unsigned &nnz, const int *ltgSorted,
        const int *ltgUnsorted, const int *rowPtr, const int *colInd, int &Dof)
{
  int indexBase = 1; //0 based indexing for C/C++ / 1 for Fortran
  dof = Dof; //constant size dof blocks
  ghostAndLocalNodes = numGhostAndLocalNodes;
  localNodes = numLocalNodes;
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  bool new_mapping_pattern = false; //could change with remeshing

  // check if mesh was changed and vectors need to be resized
  if (localToGlobalSorted.size() != numGhostAndLocalNodes) {
    localToGlobalSorted.clear();   // sets size to zero
    localToGlobalUnsorted.clear();
    nnzPerRow.clear();
    globalColInd.clear();
  }

  // allocate memory for vectors
  if (localToGlobalSorted.size() == 0) { // vectors not allocated or cleared
    localToGlobalSorted.reserve(numGhostAndLocalNodes);
    localToGlobalUnsorted.reserve(numGhostAndLocalNodes);
    nnzPerRow.reserve(numGhostAndLocalNodes);
    globalColInd.reserve(nnz);
    new_mapping_pattern = true;
  }

  //define localtoglobal to be used for unqiue partition map
  //only take ltgSorted(1:numLocalNodes) since those are owned by the processor
  if (new_mapping_pattern) {
    for (unsigned i = 0; i < numLocalNodes; ++i) {
      //any nodes following are ghost nodes so that subset
      localToGlobalSorted.emplace_back(ltgSorted[i]);
    }
  }

  // Create blockmap for arbitrary distributon of constant size, dof, element
  // amongst the procs. This is a non-overlapping map defining a unique partition
  // localToGlobalSorted[0] contains global index value of ith owned element on
  // the calling processor
  Trilinos::blockMap = new Epetra_BlockMap(numGlobalNodes, numLocalNodes,
                       &localToGlobalSorted[0], dof, indexBase, comm);

  //Create blockmap which includes the ghost nodes to be imported into the
  //solution vector at the end
  int inc_ghost = -1; // since including ghost nodes sum will not equal total
  Epetra_BlockMap ghostMap(inc_ghost, numGhostAndLocalNodes,  ltgSorted, dof,
          indexBase, comm);

  // Calculate nnzPerRow to pass into graph constructor
  if (new_mapping_pattern) {
    nnzPerRow.clear (); //destroy vector
    for (unsigned i = 0; i < numLocalNodes; ++i)
      nnzPerRow.emplace_back(rowPtr[i+1] - rowPtr[i]);
  } //fi

  // Construct graph based on nnz per row
  Trilinos::K_graph = new Epetra_FECrsGraph(Copy, *Trilinos::blockMap,
              &nnzPerRow[0]);

  unsigned nnzCount = 0; //cumulate count of block nnz per rows

  // Use unsortedltg to map local col ind to global indices
  if (globalColInd.size() != nnz) { // only if nnz changed
    globalColInd.clear(); // destroy
    for (unsigned i = 0; i < nnz; ++i) {
      // Convert to global indexing subtract 1 for fortran vs C array indexing
      globalColInd.emplace_back(ltgUnsorted[colInd[i] - 1]);
    }
  } //if

  //loop over block rows owned by current proc using localToGlobal index pointer
  for (unsigned i = 0; i < numGhostAndLocalNodes; ++i)
  {
    if (new_mapping_pattern && i >= numLocalNodes)
      nnzPerRow.emplace_back(rowPtr[i+1] - rowPtr[i]);
    int numEntries = nnzPerRow[i]; //nnz per row
    int error = 0;
    int num_rows_inserting = 1; //number of rows-inserting one row at a time
    error = Trilinos::K_graph->InsertGlobalIndices(num_rows_inserting,
            &ltgUnsorted[i], numEntries, &globalColInd[nnzCount]);

    //colInd is indexed by nnz per rows processed so far
    if (error != 0)
    {
      std::cout << "ERROR: Inserting Global Indices into Map" << std::endl;
      exit(1);
    }

    // final result should be nnz and so will not go out of bounds on colInd
    // up to nnz -1
    nnzCount += numEntries;
    //store into global vector
    if (new_mapping_pattern)
      localToGlobalUnsorted.emplace_back(ltgUnsorted[i]);
  } // for

  //by end of iterations nnzCount should equal nnz-otherwise there is an error
  if (nnzCount != nnz)
  {
    std::cout << "Number of Entries Entered in Graph does not equal nnz " << std::endl;
    exit(1);
  }

  //check if trilinos methods are successful
  if (Trilinos::K_graph->GlobalAssemble() != 0) //Calls FillComplete
  {
    std::cout << "ERROR: Calling Fill Complete on Graph" << std::endl;
    exit(1);
  }

  //Trilinos::K_graph->OptimizeStorage(); //TODO: that causes it to fail

  // --- Create block finite element matrix from graph with fillcomplete ------
  // construct matrix from filled graph
  Trilinos::K = new Epetra_FEVbrMatrix(Copy, *Trilinos::K_graph);
  //construct RHS force vector F topology
  Trilinos::F = new Epetra_FEVector(*Trilinos::blockMap);
  Trilinos::bdryVec = new Epetra_FEVector(*Trilinos::blockMap);

  // Initialize solution vector which is unique and does not include the ghost
  // indices using the unique map
  Trilinos::X = new Epetra_Vector(*Trilinos::blockMap);
  //initialize vector which will import the ghost nodes using the ghost map
  Trilinos::ghostX = new Epetra_Vector(ghostMap);
  //Create importer of the two maps
  Trilinos::Importer = new Epetra_Import(ghostMap, *Trilinos::blockMap);
  //need to give v a map-same as F
  if (new_mapping_pattern)
    for (unsigned i = numLocalNodes; i < numGhostAndLocalNodes; ++i)
      localToGlobalSorted.emplace_back(ltgSorted[i]);
} // trilinos_lhs_create_

// ----------------------------------------------------------------------------
/**
 * This function assembles the dense element block stiffness matrix and element
 * block force vector into the global block stiffness K and global block force
 * vector F.
 * Happens on a single element and will be called repeatedly as we loop over
 * all the elements.
 *
 * \param numNodesPerElement number of nodes (d) for element e
 * \param eqN                converts local element node number to local proc
 *                           equation number-size is numNodesPerElement
 * \param lK                 element local stiffness of size
 *                           (dof*dof, numNodesPerElement, numNodesPerElement)
 * \param lR                 element local force vector of size
 *                           (dof, numNodesPerElement)
 */
void trilinos_doassem_(int &numNodesPerElement, const int *eqN,
        const double *lK, double *lR)
{
  //dof values per global ID in the force vector
  int numValuesPerID = dof;

  //converts eqN in local proc values to global values
  std::vector<int> localToGlobal(numNodesPerElement);

  //subtract 1 since eqN is 1 based and C is 0 based indeing
  for (int i = 0; i < numNodesPerElement; ++i)
    localToGlobal[i] = localToGlobalUnsorted[eqN[i] - 1];

  //loop over local nodes on the element
  for (int a = 0; a < numNodesPerElement; ++a)
  {
    // Sum into contributions from element node-global assemble will take those
    // from shared nodes on other processors since FE routine.
    int error = Trilinos::K->BeginSumIntoGlobalValues(localToGlobal[a],
            numNodesPerElement, &localToGlobal[0]);

    if (error != 0)
    {
      std::cout << "ERROR: Setting block row and block column of summing into "
                << "global values!"
                << std::endl;
      exit(1);
    }

    //submit global F values
    int num_block_rows = 1; //number of global block rows put in 1 at a time
    Trilinos::F->SumIntoGlobalValues (num_block_rows, &localToGlobal[a],
            &numValuesPerID, &lR[a*dof]);

    if (error != 0) //1 or more indices not associated with calling processor!
    {
      std::cout << error << std::endl;
      std::cout << "ERROR: Summing into global vector values!" << std::endl;
      exit(1);
    }

    //loop over local nodes for columns
     std::vector<double> values(dof*dof);
    for (int b = 0; b < numNodesPerElement; ++b)
    {
      //transpose block since Trilinos takes in SerialMAtrix in column major
      for (int i = 0; i < dof; ++i)
      {
        //premult by diagonal W so W(a) multiplies row a
        for (int j = 0; j < dof; ++j)
        {
          //taking transpose of block so flip i & j
          values[i*dof + j]
              = lK[b*dof*dof*numNodesPerElement + a*dof*dof + j*dof + i];
        }
      }

      error = Trilinos::K->SubmitBlockEntry(&values[0], dof, dof, dof);

      if (error != 0)
      {
        std::cout << "ERROR: Inputting values of summing into matrix "
                  << "global values!"
                  << std::endl;
        exit(1);
      }
    }

    //finish submitting the block entry for the current row
    error = Trilinos::K->EndSubmitEntries();
    if (error != 0)
    {
      std::cout << "ERROR: End submitting block entries!" << std::endl;
      exit(1);
    }
  }
} // trilinos_doassem_

// ----------------------------------------------------------------------------
/**
 * Tthis function creates the LHS structure topology and RHS vector based on
 * the block map using the global stiffness matrix and global force vector
 * assuming the assembly part has been done in svFSI
 *
 * \param Val         nonzero values of assembled global stiffness
 * \param RHS         global force vector
 * \param x           solution vector
 * \param dirW        gives information about the Dirichlet BC
 * \param resNorm     norm of solution
 * \param initNorm    initial norm of RHS
 * \param numIters    linear solver number of iterations
 * \param solverTime  time in linear solver
 * \param dB          log ratio for output
 * \param converged   neither or not converged in the max number of iterations
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 */
void trilinos_global_solve_(const double *Val, const double *RHS, double *x,
        const double *dirW, double &resNorm, double &initNorm, int &numIters,
        double &solverTime, double &dB, bool &converged, int &lsType,
        double &relTol, int &maxIters, int &kspace, int &precondType)
{
  int nnzCount = 0; //cumulate count of block nnz per rows
  int count = 0;
  int numValuesPerID = dof; //dof values per id pointer to dof
  std::vector<double> values(dof*dof); // holds local matrix entries

  //loop over block rows owned by current proc using localToGlobal index pointer
  for (int i = 0; i < ghostAndLocalNodes; ++i)
  {
    int numEntries = nnzPerRow[i]; //block per of entries per row
    //copy global stiffness values
    int error = Trilinos::K->BeginReplaceGlobalValues(localToGlobalUnsorted[i],
            numEntries, &globalColInd[nnzCount]);

    //need to check if globalColInd and localToGlobal equal for block diag
    if (error != 0)
    {
      std::cout << error << std::endl;
      std::cout << "ERROR: Setting block row and block column of setting "
                << "global values!"
                << std::endl;
      exit(1);
    }

    // Copy F values-loop over dof for bool
    int num_block_rows = 1; //number of global block rows put in 1 at a time
    error = Trilinos::F->ReplaceGlobalValues (num_block_rows,
            &localToGlobalUnsorted[i], &numValuesPerID, &RHS[i*dof]);

    //check is bool true or false whether to give 0 if on diagonal 1
    if (error != 0)
    {
      std::cout << "ERROR: Setting global vector values!" << std::endl;
      exit(1);
    }
    for (int j = 0; j < numEntries; ++j)
    {
      for (int l = 0; l < dof; ++l) //loop over dof for bool to contruct
      {
        for (int m = 0; m < dof; ++m)
          values[l*dof + m] = Val[count*dof*dof + m*dof + l]; //transpose it
      }
      //submit square dof*dof blocks
      error = Trilinos::K->SubmitBlockEntry(&values[0], dof, dof, dof);

      if (error != 0)
      {
        std::cout << "ERROR: Inputting values of setting matrix global "
                  << "values!" << std::endl;
        exit(1);
      }
      count++;
    }
    error = Trilinos::K->EndSubmitEntries(); //for current block row
    if (error != 0)
    {
      std::cout << "ERROR: End submitting block entries!" << std::endl;
      exit(1);
    }
    nnzCount += numEntries;
  }

  //call solver code which assembles K and F for shared processors
  trilinos_solve_(x, dirW, resNorm, initNorm, numIters,
          solverTime, dB, converged, lsType,
          relTol, maxIters, kspace, precondType);

} // trilinos_global_solve_

// ----------------------------------------------------------------------------
/**
 * This function uses the established data structure to solve the block linear
 * system and passes the solution vector back to fortran with the ghost nodes
 * filled in
 *
 * \param x           solution vector
 * \param dirW        gives information about Dirichlet BC
 * \param resNorm     norm of the final residual
 * \param initNorm    norm of initial resiual x_init = 0 ||b||
 * \param numIters    number of iterations in the linear solver
 * \param solverTime  time in the linear solver
 * \param dB          log ratio for output
 * \param converged   can pass in solver and preconditioner type too
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 */
void trilinos_solve_(double *x, const double *dirW, double &resNorm,
        double &initNorm, int &numIters, double &solverTime, double &dB,
        bool &converged, int &lsType, double &relTol, int &maxIters,
        int &kspace, int &precondType)
{
  // Already filled from graph so does not need to call fillcomplete
  // routine will sum in contributions from elements on shared nodes amongst
  // processors
  int error = Trilinos::K->GlobalAssemble(false);
  if (error != 0)
  {
    std::cout << "ERROR: Global Assembling stiffness matrix" << std::endl;
    exit(1);
  }

  //very important for performance-makes memory contiguous
  Trilinos::K->OptimizeStorage();

  //sum in values from shared nodes amongst the processors
  error = Trilinos::F->GlobalAssemble();

  if (error != 0)
  {
    std::cout << "ERROR: Global Assembling force vector" << std::endl;
    exit(1);
  }

  // Construct Jacobi scaling vector which uses dirW to take the Dirichlet BC
  // into account
  Epetra_Vector diagonal(*Trilinos::blockMap);
  constructJacobiScaling(dirW, diagonal);

  // Compute norm of preconditioned multivector F that we will be solving
  // problem with
  Trilinos::F->Norm2(&initNorm); //pass preconditioned norm W*F

  // Define Epetra_Operator which is global stiffness with coupled boundary
  // conditions included
  TrilinosMatVec K_bdry;

  //Define linear problem if v is 0 does standard matvec product with K
  Epetra_LinearProblem Problem(&K_bdry, Trilinos::X, Trilinos::F);
  AztecOO Solver(Problem);

  //Can set output solver parameter options below
#ifdef NOOUTPUT
  Solver.SetAztecOption(AZ_diagnostics, AZ_none);
  Solver.SetAztecOption(AZ_output, AZ_none);
#endif

  setPreconditioner(precondType, Solver);

  // Set convergence type as relative ||r|| <= relTol||b||
  Solver.SetAztecOption(AZ_conv, AZ_rhs);

  //Set solver options
  int numRestarts = 1; //also changes for gmres
  int maxItersPerRestart = maxIters;

  // Solver is GMRES by default
  if (lsType ==TRILINOS_GMRES_SOLVER)
  {
    //special parameters to set orthog and kspace
    numRestarts = maxIters; //different definition for gmres
    Solver.SetAztecOption(AZ_orthog, AZ_modified); //modified GS
    Solver.SetAztecOption(AZ_kspace, kspace);
    Solver.SetAztecOption(AZ_solver, AZ_gmres);
    maxItersPerRestart = kspace; //total maxIters is kspace * numRestarts
  }
  else if (lsType == TRILINOS_BICGSTAB_SOLVER)
    Solver.SetAztecOption(AZ_solver, AZ_bicgstab);

  else if (lsType == TRILINOS_CG_SOLVER)
    Solver.SetAztecOption(AZ_solver, AZ_cg);

  //checkStatus to calculate residual norm
  AztecOO_StatusTestResNorm restartResNorm(K_bdry, *Trilinos::X,
                      (Epetra_Vector&) Trilinos::F[0], relTol);
  restartResNorm.DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                 AztecOO_StatusTestResNorm::TwoNorm);
  Solver.SetStatusTest(&restartResNorm);

  int status;
  numIters = 0;
  solverTime = 0.0;
  for (int iter = 0; iter < numRestarts; ++iter)
  {
    dB = (iter > 0) ? resNorm : initNorm;
    if (numRestarts > 1) numIters += 1; //gmres case
    status = Solver.Iterate(maxItersPerRestart, relTol);
    numIters += Solver.NumIters();
    resNorm = restartResNorm.GetResNormValue();
    solverTime += Solver.SolveTime();
    if (resNorm < relTol * initNorm) break;
  }
  converged = (status == 0) ? true : false;
  dB = 10 * log(restartResNorm.GetResNormValue()/dB); //fits with gmres def

  //Right scaling so need to multiply x by diagonal
  Trilinos::X->Multiply(1.0, *Trilinos::X, diagonal, 0.0);

  //Fill ghost X with x communicating ghost nodes amongst processors
  error = Trilinos::ghostX->Import(*Trilinos::X, *Trilinos::Importer, Insert);
  //check imported correctly
  if (error != 0)
  {
    std::cout << "ERROR: Map ghost node importer!" << std::endl;
     exit(1);
  }

    error = Trilinos::ghostX->ExtractCopy(x);
   if (error != 0)
   {
     std::cout << "ERROR: Extracting copy of solution vector!" << std::endl;
     exit(1);
  }
  //set to 0 for the next time iteration
  Trilinos::K->PutScalar(0.0);
  Trilinos::F->PutScalar(0.0);
  if (coupledBC) Trilinos::bdryVec->PutScalar(0.0);
  //0 out initial guess for iteration
  Trilinos::X->PutScalar(0.0);
  // Free memory if MLPrec is invoked
  if (ifpackPrec) {
      delete ifpackPrec;
      ifpackPrec = NULL;
  }
  if (MLPrec) {
      delete MLPrec;
      MLPrec = NULL;
  }
} // trilinos_solve_

// ----------------------------------------------------------------------------
void setPreconditioner(int precondType, AztecOO &Solver)
{
  //initialize reordering for ILU/ILUT preconditioners
  Solver.SetAztecOption(AZ_reorder, 1);
  //solve precond structure into separate function
  if (precondType == TRILINOS_DIAGONAL_PRECONDITIONER ||
      precondType == NO_PRECONDITIONER )
    Solver.SetAztecOption(AZ_precond, AZ_none);
  else if (precondType == TRILINOS_BLOCK_JACOBI_PRECONDITIONER)
  {
    Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    checkDiagonalIsZero();
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if(precondType == TRILINOS_ILU_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    Solver.SetAztecOption(AZ_overlap,1);
    Solver.SetAztecOption(AZ_graph_fill,0);
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if (precondType == TRILINOS_ILUT_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    Solver.SetAztecOption(AZ_overlap,1);
    Solver.SetAztecOption(AZ_graph_fill,4);
    Solver.SetAztecParam(AZ_ilut_fill, 2);
    Solver.SetAztecParam(AZ_drop, 1e-2);

    // Sets the global stiffness prior to rank 1 update as the matrix off of
    // which the preconditioner to calculated from to utilize the native
    // preconditioners from trilinos
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if (precondType == TRILINOS_IC_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    setIFPACKPrec(Solver);
  }
  else if (precondType == TRILINOS_ICT_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    setIFPACKPrec(Solver); //add in parameter for string for different types
  }
  else if (precondType == TRILINOS_ML_PRECONDITIONER)
    setMLPrec(Solver);
  else
  {
    std::cout << "ERROR: Preconditioner Type is undefined" << std::endl;
    exit(1);
  }
} // setPreconditioner

// ----------------------------------------------------------------------------
/**
 * Tune parameters for htis and IFPACK
 * Ref: https://trilinos.org/oldsite/packages/ml/mlguide5.pdf
 */
void setMLPrec(AztecOO &Solver)
{
  //break up into initializer
  Teuchos::ParameterList MLList;
  int *options = new int[AZ_OPTIONS_SIZE];
  double *params = new double[AZ_PARAMS_SIZE];
  ML_Epetra::SetDefaults("SA",MLList, options, params);

  // ML general options
  // output level, 0 being silent and 10 verbose
  bool verbose = false;
  if (verbose)
    MLList.set("ML output", 10);
  else
    MLList.set("ML output", 0);

  // ML Cycle options
  // maximum number of levels possible
  MLList.set("max levels",4);
  //If set to increasing, level 0 will correspond to the finest level.
  //If set to decreasing,max levels - 1 will correspond to the finest level.
  //Default: increasing.
  MLList.set("increasing or decreasing","increasing");

  // Aggregation & prolongator parameters
  // coarsening options:  Uncoupled, MIS, Uncoupled-MIS (uncoupled on the finer grids, then switch to MIS)
  //MLList.set("aggregation: type", "Uncoupled");
  MLList.set("aggregation: type", "MIS");
  MLList.set("aggregation: threshold", 0.0);

  // Smoother parameters
  //common smoother options: Chebyshev, Gauss-Seidel, symmetric Gauss-Seidel, Jacobi, ILU, IC
  MLList.set("smoother: ifpack type","ILU");
  MLList.set("smoother: ifpack overlap",1);
  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill",1);
  MLList.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");
  MLList.set("subsmoother: type", "symmetric Gauss-Seidel");

  //Coarse grid parameters
  MLList.set("coarse: max size", 15);
  // Coarse solver Default: AmesosKLU.
  MLList.set("coarse: type","symmetric Gauss-Seidel");

  // Load balancing options
  MLList.set("repartition: enable",1);
  MLList.set("repartition: max min ratio",1.3);
  MLList.set("repartition: min per proc",500);
  MLList.set("repartition: partitioner","Zoltan");
  MLList.set("repartition: Zoltan dimensions",2);

  // create the preconditioner object based on options in MLList and compute hierarchy
  if (MLPrec == NULL)
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*Trilinos::K, MLList, false);

  timecount = 0;
  if (timecount == 0)
    MLPrec->ComputePreconditioner();
  else {//switch to recompute
    MLPrec->ReComputePreconditioner();
  }
  Solver.SetPrecOperator(MLPrec);

  //solver to separte function only recompute at separate time iter
  timecount += 1;

  delete[] options;
  delete[] params;
}// setMLPrec

// ----------------------------------------------------------------------------
/**
 * pass in IC, ICT
 * pass in string for which to turn on right now set to IC
 */
void setIFPACKPrec(AztecOO &Solver)
{
  //Ifpack Factory;
  //std::string PrecType = "ILUT"; // exact solve on each subdomain
  //int OverlapLevel = 0; // one row of overlap among the processes
  //ifpackPrec = Factory.Create(PrecType, Trilinos::K, OverlapLevel);
  //Teuchos::ParameterList List;
  //List.set("fact: level-of-fill", 2);
  //List.set("fact: drop tolerance", 1e-2);
  //ifpackPrec->SetParameters(List);
  //ifpackPrec->Initialize();
  //ifpackPrec->Compute();
  //Solver.SetPrecOperator(&*ifpackPrec);

  Teuchos::ParameterList List;
  int OverlapLevel = 0;
  ifpackPrec = new Ifpack_AdditiveSchwarz<Ifpack_ILUT> (Trilinos::K, OverlapLevel);
  List.set("fact: ict level-of-fill", 2.0);
  List.set("fact: drop tolerance", 1e-2);
  ifpackPrec->SetParameters(List);
  ifpackPrec->Initialize();
  ifpackPrec->Compute();
  Solver.SetPrecOperator(&*ifpackPrec);

} // setIFPACKPrec


// ----------------------------------------------------------------------------
/**
 * This routine is to be used with preconditioners such as ILUT which require
 * 1s on the diagonal
 */
void checkDiagonalIsZero()
{
  Epetra_Vector diagonal(*Trilinos::blockMap);
  Trilinos::K->ExtractDiagonalCopy(diagonal);
  bool isZeroDiag = false; //initialize to false
  for (int i = 0; i < diagonal.MyLength(); ++i)
  {
    //if diagonal is 0 change it to be 1
    if (diagonal[i] == 0.0)
    {
      diagonal[i] = 1.0;
      isZeroDiag = true;
    }
  }
  if (isZeroDiag) Trilinos::K->ReplaceDiagonalValues(diagonal);
} // void checkDiagonalIsZero()

// ----------------------------------------------------------------------------
/**
 * To be called within solve-output diagonal from here
 *
 * \param dirW    pass in array with Dirichlet boundary face nodes marked
 * \paramdiagonal diagonal scaling vector need to output to multiply solution by
 */
void constructJacobiScaling(const double *dirW, Epetra_Vector &diagonal)
{
  //loop over nodes owned by that processor
  for (int i = 0; i < localNodes; ++i)
  {
    for (int j = 0; j < dof; ++j)
    {
      int error = diagonal.ReplaceGlobalValue(localToGlobalSorted[i],
                          j, //block offset
                          0, //multivector of 1 vector
                          dirW[i*dof+j]); //value to insert
      if (error != 0)
      {
        std::cout << "ERROR: Setting Dirichlet diagonal scaling values!" << std::endl;
        exit(1);
      }
    }
  }

  //Extract diagonal of K
  Epetra_Vector Kdiag(*Trilinos::blockMap);
  Trilinos::K->ExtractDiagonalCopy(Kdiag);
  for (int i = 0; i < Kdiag.MyLength(); ++i)
  {
    if (Kdiag[i] == 0.0) Kdiag[i] = 1.0;
    Kdiag[i] = 1 / sqrt(abs(Kdiag[i])); //Jacobi scaling 1 / sqrt|aii|
  }

  //multiply the two vectors together
  diagonal.Multiply(1.0, diagonal, Kdiag, 0.0);

  //Multiply K and Fon the left by diagonal
  //Let diag = W, solving W*K*W*y = W*F, where X = W*y
  Trilinos::K->LeftScale(diagonal);

  // Elem by elem multiplication to support diagonal matrix multiply of vector
  // this -> 0.0*this + 1.0*F*diag
  Trilinos::F->Multiply(1.0, *Trilinos::F, diagonal, 0.0);

  //Right scaling of K need to multiply solution vector as well
  Trilinos::K->RightScale(diagonal);

  //Need to multiply v in vv' by diagonal
  if (coupledBC)
      Trilinos::bdryVec->Multiply(1.0, *Trilinos::bdryVec, diagonal, 0.0);
} // void constructJacobiScaling()

// ----------------------------------------------------------------------------
/**
 * \param  v            coupled boundary vector
 * \param  isCoupledBC  determines if coupled resistance BC is turned on
 */
void trilinos_bc_create_(const double *v, bool &isCoupledBC)
{
  //store as global to determine which matvec multiply to use in solver
  coupledBC = isCoupledBC;

  int error = 0;

  if (isCoupledBC)
  {
    //loop over block rows owned by current proc using localToGlobal index ptr
    for (int i = 0; i < ghostAndLocalNodes; ++i)
    {
      //sum values into v fe case
      int num_global_rows = 1; //number of global block rows put in 1 at a time
      error = Trilinos::bdryVec->ReplaceGlobalValues (num_global_rows,
               &localToGlobalSorted[i],  //global idx id-inputting sorted array
               &dof, //dof values per id pointer to dof
               &v[i*dof]); //values of size dof
      if (error != 0)
      {
        std::cout << "ERROR: Setting boundary vector values!" << std::endl;
        exit(1);
      }
    }
  }
} // void trilinos_bc_create_()

// ----------------------------------------------------------------------------
/**
 * free the memory of the global data structures
 */
void trilinos_lhs_free_()
{
  if (MLPrec) {
      MLPrec->DestroyPreconditioner();
      delete MLPrec;
      MLPrec = NULL;
  }
  if (Trilinos::blockMap) {
      delete Trilinos::blockMap;
      Trilinos::blockMap = NULL;
  }
  if (Trilinos::F) {
      delete Trilinos::F;
      Trilinos::F = NULL;
  }
  if (Trilinos::K) {
      delete Trilinos::K;
      Trilinos::K = NULL;
  }
  if (Trilinos::X) {
      delete Trilinos::X;
      Trilinos::X = NULL;
  }
  if (Trilinos::ghostX) {
      delete Trilinos::ghostX;
      Trilinos::ghostX = NULL;
  }
  if (Trilinos::Importer) {
      delete Trilinos::Importer;
      Trilinos::Importer = NULL;
  }
  if (Trilinos::bdryVec) {
      delete Trilinos::bdryVec;
      Trilinos::bdryVec = NULL;
  }
  if (Trilinos::K_graph) {
      delete Trilinos::K_graph;
      Trilinos::K_graph = NULL;
  }

}

// ----------------------------------------------------------------------------
/**
 * for debugging purposes here are routines to print the matrix and RHS vector
 */
void printMatrixToFile()
{
  std::ofstream Kfile("K.txt");
  Kfile << std::scientific;
  Kfile.precision(17);

  //loop over block rows owned by current proc using localToGlobal index pointer
  for (int i = 0; i < ghostAndLocalNodes; ++i)
  {
    int numEntries = nnzPerRow[i]; //block per of entries per row
    //copy global stiffness values
    int rowDim, numBlockEntries;
    std::vector<int> blockIndices(numEntries);
    //int *blockIndices = new int[numEntries];
    std::vector<int> colDims(numEntries);
    //int * colDim = new int[numEntries];
    Trilinos::K->BeginExtractGlobalBlockRowCopy(localToGlobalUnsorted[i],
                        numEntries, rowDim, numBlockEntries, &blockIndices[0],
                        &colDims[0]);

    for (int j = 0; j < numEntries; ++j)
    {
      std::vector<double> values(dof*dof);
      int sizeofValues = dof*dof;
      int LDA = dof;
      Trilinos::K->ExtractEntryCopy(sizeofValues, &values[0], LDA, false);
      //print returned block
      for (int k = 0; k < dof; ++k)
      {
        for (int l = 0; l < dof; ++l)
        {
          Kfile << values[l*dof + k] << " ";
        }
      }
    }
  }
  Kfile.close();
} // printMatrixToFile()


void printRHSToFile()
{
  std::ofstream Ffile("F.txt");
  Ffile.precision(17);
  std::vector<double> F(Trilinos::F->MyLength());
  //extract copy to print values
  Trilinos::F->ExtractCopy(&F[0], 0);
  //print values to file
  for (int i = 0; i < Trilinos::F->MyLength(); ++i)
    Ffile << F[i] << std::endl; //Jacobi preconditioning on the lefts;
  Ffile.close();
}

// ----------------------------------------------------------------------------
/**
 */
void printSolutionToFile()
{
  std::ofstream Xfile("X.txt");
  Xfile.precision(17);
  std::vector<double> X(Trilinos::X->MyLength());
  //extract copy to print values
  Trilinos::X->ExtractCopy(&X[0], 0);
  //print values to file
  for (int i = 0; i < Trilinos::X->MyLength(); ++i)
    Xfile << X[i] << std::endl; //Jacobi preconditioning on the lefts;
  Xfile.close();
}

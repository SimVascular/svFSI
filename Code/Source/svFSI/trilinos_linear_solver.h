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

#ifndef TRILINOS_LINEAR_SOLVER_H
#define TRILINOS_LINEAR_SOLVER_H
/*!
  \file    trilinos_linear_solver.h
  \brief   wrap Trilinos solver functions
*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <stdio.h>
#include <vector>
#include <iostream>
#include "mpi.h"
#include <time.h>
#include <numeric>

// Epetra includes
#include "Epetra_MpiComm.h" //include MPI communication
#include "Epetra_Map.h" //need to create block map
#include "Epetra_FEVbrMatrix.h" //sparse matrix for FE
#include "Epetra_FEVector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"

// AztecOO includes
#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"

// ML includes
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

// Ifpack includes
#include "Ifpack.h"
#include "Ifpack_ConfigDefs.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"

/**************************************************************/
/*                      Macro Definitions                     */
/**************************************************************/

// Define linear solver as following naming in FSILS_struct
#define TRILINOS_CG_SOLVER 798
#define TRILINOS_GMRES_SOLVER 797
#define TRILINOS_BICGSTAB_SOLVER 795

// Define preconditioners as following naming in FSILS_struct
#define NO_PRECONDITIONER 700
#define TRILINOS_DIAGONAL_PRECONDITIONER 702
#define TRILINOS_BLOCK_JACOBI_PRECONDITIONER 703
#define TRILINOS_ILU_PRECONDITIONER 704
#define TRILINOS_ILUT_PRECONDITIONER 705
#define TRILINOS_IC_PRECONDITIONER 706
#define TRILINOS_ICT_PRECONDITIONER 707
#define TRILINOS_ML_PRECONDITIONER 708

// Initialize all Epetra types we need separate from Fortran
struct Trilinos
{
  static Epetra_BlockMap *blockMap;
  static Epetra_FEVector *F;
  static Epetra_FEVbrMatrix *K;
  static Epetra_Vector *X;
  static Epetra_Vector *ghostX;
  static Epetra_Import *Importer;
  static Epetra_FEVector *bdryVec;
  static Epetra_MpiComm *comm;
  static Epetra_FECrsGraph *K_graph;
};

/**
 * \class TrilinosMatVec
 * \brief This class implements the pure virtual class Epetra_Operator for the
 *        AztecOO iterative solve which only uses the Apply() method to compute
 *        the matrix vector product
 */
class TrilinosMatVec: public virtual Epetra_Operator
{
public:

  /** Define matrix vector operation at each iteration of the linear solver
   *  adds on the coupled neuman boundary contribution to the matrix
   *
   *  \param x vector to be applied on the operator
   *  \param y result of sprase matrix vector multiplication
   */
  int Apply(const Epetra_MultiVector &x, Epetra_MultiVector &y) const;

  /** Tells whether to use the transpose of the matrix in each matrix
   * vector product */
  int SetUseTranspose(bool use_transpose)
  {
    return Trilinos::K->SetUseTranspose(use_transpose);
  }

  /// Computes A_inv*x
  int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    return Trilinos::K->ApplyInverse(X,Y);
  }

  /// Infinity norm for global stiffness does not add in the boundary term
  double NormInf() const
  {
    return Trilinos::K->NormInf();
  }

  /// Returns a character string describing the operator
  const char * Label() const
  {
    return Trilinos::K->Label();
  }

  /// Returns current UseTranspose setting
  bool UseTranspose() const
  {
    return Trilinos::K->UseTranspose();
  }

  /// Returns true if this object can provide an approx Inf-norm false otherwise
  bool HasNormInf() const
  {
    return Trilinos::K->HasNormInf();
  }

  /// Returns pointer to Epetra_Comm communicator associated with this operator
  const Epetra_Comm &Comm() const
  {
    return Trilinos::K->Comm();
  }

  /// Returns Epetra_Map object assoicated with domain of this operator
  const Epetra_Map &OperatorDomainMap() const
  {
    return Trilinos::K->OperatorDomainMap();
  }

  /// Returns the Epetra_Map object associated with teh range of this operator
  const Epetra_Map &OperatorRangeMap() const
  {
    return Trilinos::K->OperatorRangeMap();
  }

};// class TrilinosMatVec

//  --- Functions to be called in fortran -------------------------------------

#ifdef __cplusplus
  extern "C"
  {
#endif
  /// Give function definitions which will be called through fortran
  void trilinos_lhs_create_(unsigned &numGlobalNodes, unsigned &numLocalNodes,
          unsigned &numGhostAndLocalNodes, unsigned &nnz, const int *ltgSorted,
          const int *ltgUnsorted, const int *rowPtr, const int *colInd,
          int &dof);

  /**
   * \param v           coeff in the scalar product
   * \param isCoupledBC determines if coupled resistance BC is turned on
   */
  void trilinos_bc_create_(const double *v, bool &isCoupledBC);

  void trilinos_doassem_(int &numNodesPerElement, const int *eqN,
          const double *lK, double *lR);

  void trilinos_global_solve_(const double *Val, const double *RHS,
          double *x, const double *dirW, double &resNorm, double &initNorm,
          int &numIters, double &solverTime, double &dB, bool &converged,
          int &lsType, double &relTol, int &maxIters, int &kspace,
          int &precondType);

  void trilinos_solve_(double *x, const double *dirW, double &resNorm,
          double &initNorm, int &numIters, double &solverTime,
          double &dB, bool &converged, int &lsType, double &relTol,
          int &maxIters, int &kspace, int &precondType);

  void trilinos_lhs_free_();

#ifdef __cplusplus  /* this brace matches the one on the extern "C" line */
  }
#endif

// --- Define functions to only be called in C++ ------------------------------
void setPreconditioner(int precondType, AztecOO &Solver);

void setMLPrec(AztecOO &Solver);

void setIFPACKPrec(AztecOO &Solver);

void checkDiagonalIsZero();

void constructJacobiScaling(const double *dirW,
              Epetra_Vector &diagonal);

// --- Debugging functions ----------------------------------------------------
void printMatrixToFile();

void printRHSToFile();

void printSolutionToFile();

#endif //TRILINOS_LINEAR_SOLVER_H

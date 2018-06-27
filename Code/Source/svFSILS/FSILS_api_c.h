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
 *  Routines prototypes are declared here. This provides the details 
 *  of how to call these functions
 *--------------------------------------------------------------------
 */
   
void fsils_commu_create_(FSILS_commuType* commu, MPI_Aint* commi);
void fsils_commu_free_(FSILS_commuType* commu);
         
void fsils_lhs_create_c_(void** lhs, FSILS_commuType* commu, int* gnNo, int* nNo, int* nnz, int* gNodes, int* rowPtr, int* colPtr, int* nFaces);
void fsils_lhs_free_(void* lhs);

void fsils_ls_create_(FSILS_lsType* ls, int* LS_type, double* relTol, double* absTol, int* maxItr, int* dimKry, double* relTolIn, double* absTolIn, int* maxItrIn);
void fsils_ls_free_(FSILS_lsType* ls);
         
void fsils_bc_create_(void* lhs, int* faIn, int* nNo, int* dof, int* BC_type, int* gNodes, double* Val);
void fsils_bc_free_(void* lhs, int* faIn);
         
void fsils_solve_(void* lhs, FSILS_lsType* ls, int* dof, double* Ri, double* Val, int* incL, double* res);

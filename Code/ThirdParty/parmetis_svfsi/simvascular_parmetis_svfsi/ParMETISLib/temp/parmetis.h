/*
 * Copyright 1997-2003, Regents of the University of Minnesota
 *
 * parmetis.h
 *
 * This file contains function prototypes and constrant definitions for 
 * ParMETIS
 *
 * Started 7/21/03
 * George
 *
 */

#ifndef __parmetis_h__
#define __parmetis_h__

#include <mpi.h>

#ifndef _MSC_VER
#define __cdecl
#endif


/*************************************************************************
* Data-structures
**************************************************************************/
/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype
#else
typedef short idxtype
#endif


/*************************************************************************
* Constants 
**************************************************************************/
#define PARMETIS_MAJOR_VERSION        3
#define PARMETIS_MINOR_VERSION        2
#define PARMETIS_SUBMINOR_VERSION     0


/*************************************************************************
* Function prototypes
**************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

/*-------------------------------------------------------------------
* API Introduced with Release 3.0 (current API) 
*--------------------------------------------------------------------*/
void __cdecl ParMETIS_V3_AdaptiveRepart(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *vsize, idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, 
	     int *nparts, float *tpwgts, float *ubvec, float *ipc2redist, 
	     int *options, int *edgecut, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_V3_PartGeomKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, 
	     int *ncon, int *nparts, float *tpwgts, float *ubvec, int *options, 
	     int *edgecut, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_V3_PartGeom(
             idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_V3_PartKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
	     float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
	     MPI_Comm *comm)

void __cdecl ParMETIS_V3_Mesh2Dual(
             idxtype *elmdist, idxtype *eptr, idxtype *eind, int *numflag, 
	     int *ncommonnodes, idxtype **xadj, idxtype **adjncy, MPI_Comm *comm)

void __cdecl ParMETIS_V3_PartMeshKway(
             idxtype *elmdist, idxtype *eptr, idxtype *eind, idxtype *elmwgt, 
	     int *wgtflag, int *numflag, int *ncon, int *ncommonnodes, int *nparts, 
	     float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
	     MPI_Comm *comm)

void __cdecl ParMETIS_V3_NodeND(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, 
	     int *options, idxtype *order, idxtype *sizes, MPI_Comm *comm)

void __cdecl ParMETIS_V3_RefineKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
	     float *tpwgts, float *ubvec, int *options, int *edgecut, 
	     idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_V32_NodeND(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
             idxtype *numflag, int *mtype, int *rtype, int *p_nseps, int *s_nseps,
             float *ubfrac, int *seed, int *dbglvl, idxtype *order, 
             idxtype *sizes, MPI_Comm *comm)



/*------------------------------------------------------------------
* Backward compatibility routines with Release 2.0
*-------------------------------------------------------------------*/
void __cdecl ParMETIS_PartKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
             int *edgecut, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_PartGeomKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
	     int *wgtflag, int *numflag, int *ndims, float *xyz, int *nparts, int *options, 
	     int *edgecut, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_PartGeom(
             idxtype *vtxdist, int *ndims, float *xyz, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_PartGeomRefine(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *ndims, float *xyz, 
	     int *options, int *edgecut, idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_RefineKway(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
	     idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_RepartLDiffusion(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
	     idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_RepartGDiffusion(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
	     idxtype *adjwgt, int *wgtflag, int *numflag, int *options, int *edgecut, 
	     idxtype *part, MPI_Comm *comm)

void __cdecl ParMETIS_RepartRemap(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
	     int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, 
	     MPI_Comm *comm)

void __cdecl ParMETIS_RepartMLRemap(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, 
	     int *wgtflag, int *numflag, int *options, int *edgecut, idxtype *part, 
	     MPI_Comm *comm)

void __cdecl ParMETIS_NodeND(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
	     idxtype *order, idxtype *sizes, MPI_Comm *comm)

void __cdecl ParMETIS_SerialNodeND(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int *numflag, int *options, 
	     idxtype *order, idxtype *sizes, MPI_Comm *comm)




/*-------------------------------------------------------------------
* Backward compatibility routines with Release 1.0 
*--------------------------------------------------------------------*/
void __cdecl PARKMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, 
	     idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARGKMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
             int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARGRMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
             int ndims, float *xyz, idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARGMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int ndims, float *xyz,
             idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARRMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
	     idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARUAMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
	     idxtype *adjwgt, idxtype *part, int *options, MPI_Comm comm)

void __cdecl PARDAMETIS(
             idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt,
             idxtype *part, int *options, MPI_Comm comm)

#ifdef __cplusplus
}
#endif


/*************************************************************************
* Various constants used for the different parameters
**************************************************************************/
/* Matching types */
#define PARMETIS_MTYPE_LOCAL     1    /* Restrict matching to within processor vertices */
#define PARMETIS_MTYPE_GLOBAL    2    /* Remote vertices can be matched */

/* Separator refinement types */
#define PARMETIS_SRTYPE_GREEDY    1    /* Vertices are visted from highest to lowest gain */
#define PARMETIS_SRTYPE_2PHASE    2    /* Separators are refined in a two-phase fashion using
                                          PARMETIS_SRTYPE_GREEDY for the 2nd phase */

/* Coupling types for ParMETIS_V3_RefineKway & ParMETIS_V3_AdaptiveRepart */
#define PARMETIS_PSR_COUPLED    1    /* # of partitions == # of processors */
#define PARMETIS_PSR_UNCOUPLED  2    /* # of partitions != # of processors */


/* Debug levels (fields should be ORed) */
#define PARMETIS_DBGLVL_TIME        1      /* Perform timing analysis */
#define PARMETIS_DBGLVL_INFO        2      /* Perform timing analysis */
#define PARMETIS_DBGLVL_PROGRESS    4      /* Show the coarsening progress */
#define PARMETIS_DBGLVL_REFINEINFO  8      /* Show info on communication during folding */
#define PARMETIS_DBGLVL_MATCHINFO   16     /* Show info on matching */
#define PARMETIS_DBGLVL_RMOVEINFO   32     /* Show info on communication during folding */
#define PARMETIS_DBGLVL_REMAP       64     /* Determines if remapping will take place */

#endif 

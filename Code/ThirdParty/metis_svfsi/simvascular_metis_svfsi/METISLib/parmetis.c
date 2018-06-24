/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * parmetis.c
 *
 * This file contains top level routines that are used by ParMETIS
 *
 * Started 10/14/97
 * George
 *
 * $Id: parmetis.c,v 1.2 2003/07/24 18:39:11 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function is the entry point for KMETIS with seed specification
* in options[7] 
**************************************************************************/
void METIS_PartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                         idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                         int *options, int *edgecut, idxtype *part)
{
  int i;
  float *tpwgts;

  tpwgts = fmalloc(*nparts, "KMETIS: tpwgts");
  for (i=0; i<*nparts; i++) 
    tpwgts[i] = 1.0/(1.0*(*nparts));

  METIS_WPartGraphKway2(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, 
                       tpwgts, options, edgecut, part);
  GKfree((void **)&tpwgts, LTERM);
}


/*************************************************************************
* This function is the entry point for KWMETIS with seed specification
* in options[7] 
**************************************************************************/
void METIS_WPartGraphKway2(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
                          idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
                          float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  int i, j;
  GraphType graph;
  CtrlType ctrl;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_KMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, *wgtflag);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = KMETIS_CTYPE;
    ctrl.IType = KMETIS_ITYPE;
    ctrl.RType = KMETIS_RTYPE;
    ctrl.dbglvl = KMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_KMETIS;
  ctrl.CoarsenTo = 20*(*nparts);
  ctrl.maxvwgt = 1.5*((graph.vwgt ? idxsum(*nvtxs, graph.vwgt) : (*nvtxs))/ctrl.CoarsenTo);

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  *edgecut = MlevelKWayPartitioning(&ctrl, &graph, *nparts, part, tpwgts, 1.000);

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function is the entry point for the node ND code for ParMETIS
**************************************************************************/
void METIS_NodeNDP(int nvtxs, idxtype *xadj, idxtype *adjncy, int npes, 
                   int *options, idxtype *perm, idxtype *iperm, idxtype *sizes) 
{
  int i, ii, j, l, wflag, nflag;
  GraphType graph;
  CtrlType ctrl;
  idxtype *cptr, *cind;

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType   = ONMETIS_CTYPE;
    ctrl.IType   = ONMETIS_ITYPE;
    ctrl.RType   = ONMETIS_RTYPE;
    ctrl.dbglvl  = ONMETIS_DBGLVL;
    ctrl.oflags  = ONMETIS_OFLAGS;
    ctrl.pfactor = ONMETIS_PFACTOR;
    ctrl.nseps   = ONMETIS_NSEPS;
  }
  else {
    ctrl.CType   = options[OPTION_CTYPE];
    ctrl.IType   = options[OPTION_ITYPE];
    ctrl.RType   = options[OPTION_RTYPE];
    ctrl.dbglvl  = options[OPTION_DBGLVL];
    ctrl.oflags  = options[OPTION_OFLAGS];
    ctrl.pfactor = options[OPTION_PFACTOR];
    ctrl.nseps   = options[OPTION_NSEPS];
  }
  if (ctrl.nseps < 1)
    ctrl.nseps = 1;

  ctrl.optype = OP_ONMETIS;
  ctrl.CoarsenTo = 100;

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  InitRandom(-1);

  if (ctrl.oflags&OFLAG_COMPRESS) {
    /*============================================================
    * Compress the graph 
    ==============================================================*/
    cptr = idxmalloc(nvtxs+1, "ONMETIS: cptr");
    cind = idxmalloc(nvtxs, "ONMETIS: cind");

    CompressGraph(&ctrl, &graph, nvtxs, xadj, adjncy, cptr, cind);

    if (graph.nvtxs >= COMPRESSION_FRACTION*(nvtxs)) {
      ctrl.oflags--; /* We actually performed no compression */
      GKfree((void **)&cptr, &cind, LTERM);
    }
    else if (2*graph.nvtxs < nvtxs && ctrl.nseps == 1)
      ctrl.nseps = 2;
  }
  else {
    SetUpGraph(&graph, OP_ONMETIS, nvtxs, 1, xadj, adjncy, NULL, NULL, 0);
  }


  /*=============================================================
  * Do the nested dissection ordering 
  --=============================================================*/
  ctrl.maxvwgt = 1.5*(idxsum(graph.nvtxs, graph.vwgt)/ctrl.CoarsenTo);
  AllocateWorkSpace(&ctrl, &graph, 2);

  idxset(2*npes-1, 0, sizes);
  MlevelNestedDissectionP(&ctrl, &graph, iperm, graph.nvtxs, npes, 0, sizes);

  FreeWorkSpace(&ctrl, &graph);

  if (ctrl.oflags&OFLAG_COMPRESS) { /* Uncompress the ordering */
    if (graph.nvtxs < COMPRESSION_FRACTION*(nvtxs)) { 
      /* construct perm from iperm */
      for (i=0; i<graph.nvtxs; i++)
        perm[iperm[i]] = i; 
      for (l=ii=0; ii<graph.nvtxs; ii++) {
        i = perm[ii];
        for (j=cptr[i]; j<cptr[i+1]; j++)
          iperm[cind[j]] = l++;
      }
    }

    GKfree((void **)&cptr, &cind, LTERM);
  }


  for (i=0; i<nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
void MlevelNestedDissectionP(CtrlType *ctrl, GraphType *graph, idxtype *order, int lastvtx, 
                             int npes, int cpos, idxtype *sizes)
{
  int i, j, nvtxs, nbnd, tvwgt, tpwgts2[2];
  idxtype *label, *bndind;
  GraphType lgraph, rgraph;
  float ubfactor;

  nvtxs = graph->nvtxs;

  if (nvtxs == 0) {
    GKfree((void **)&graph->gdata, &graph->rdata, &graph->label, LTERM);
    return;
  }

  /* Determine the weights of the partitions */
  tvwgt = idxsum(nvtxs, graph->vwgt);
  tpwgts2[0] = tvwgt/2;
  tpwgts2[1] = tvwgt-tpwgts2[0];

  if (cpos >= npes-1) 
    ubfactor = ORDER_UNBALANCE_FRACTION;
  else 
    ubfactor = 1.05;


  MlevelNodeBisectionMultiple(ctrl, graph, tpwgts2, ubfactor);

  IFSET(ctrl->dbglvl, DBG_SEPINFO, printf("Nvtxs: %6d, [%6d %6d %6d]\n", graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));

  if (cpos < npes-1) {
    sizes[2*npes-2-cpos] = graph->pwgts[2];
    sizes[2*npes-2-(2*cpos+1)] = graph->pwgts[1];
    sizes[2*npes-2-(2*cpos+2)] = graph->pwgts[0];
  }

  /* Order the nodes in the separator */
  nbnd = graph->nbnd;
  bndind = graph->bndind;
  label = graph->label;
  for (i=0; i<nbnd; i++) 
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree((void **)&graph->gdata, &graph->rdata, &graph->label, LTERM);

  if (rgraph.nvtxs > MMDSWITCH || 2*cpos+1 < npes-1) 
    MlevelNestedDissectionP(ctrl, &rgraph, order, lastvtx, npes, 2*cpos+1, sizes);
  else {
    MMDOrder(ctrl, &rgraph, order, lastvtx); 
    GKfree((void **)&rgraph.gdata, &rgraph.rdata, &rgraph.label, LTERM);
  }
  if (lgraph.nvtxs > MMDSWITCH || 2*cpos+2 < npes-1) 
    MlevelNestedDissectionP(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs, npes, 2*cpos+2, sizes);
  else {
    MMDOrder(ctrl, &lgraph, order, lastvtx-rgraph.nvtxs); 
    GKfree((void **)&lgraph.gdata, &lgraph.rdata, &lgraph.label, LTERM);
  }
}


/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_NodeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
           idxtype *adjwgt, float *ubfactor, int *options, int *sepsize, idxtype *part) 
{
  int i, j, tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType  = ONMETIS_CTYPE;
    ctrl.IType  = ONMETIS_ITYPE;
    ctrl.RType  = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType  = options[OPTION_CTYPE];
    ctrl.IType  = options[OPTION_ITYPE];
    ctrl.RType  = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags    = OFLAG_COMPRESS; /* For by-passing the pre-coarsening for multiple runs */
  ctrl.RType     = 2;  /* Standard 1-sided node refinement code */
  ctrl.pfactor   = 0;
  ctrl.nseps     = 5;  /* This should match NUM_INIT_MSECTIONS in ParMETISLib/defs.h */
  ctrl.optype    = OP_ONMETIS;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/ 
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelNodeBisectionMultiple(&ctrl, &graph, tpwgts, *ubfactor*.95);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  GKfree((void **)&graph.gdata, &graph.rdata, &graph.label, LTERM);


  FreeWorkSpace(&ctrl, &graph);

}


/*************************************************************************
* This function is the entry point for ONWMETIS. It requires weights on the
* vertices. It is for the case that the matrix has been pre-compressed.
**************************************************************************/
void METIS_EdgeComputeSeparator(int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
           idxtype *adjwgt, int *options, int *sepsize, idxtype *part) 
{
  int i, j, tvwgt, tpwgts[2];
  GraphType graph;
  CtrlType ctrl;

  SetUpGraph(&graph, OP_ONMETIS, *nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);
  tvwgt = idxsum(*nvtxs, graph.vwgt);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType = ONMETIS_CTYPE;
    ctrl.IType = ONMETIS_ITYPE;
    ctrl.RType = ONMETIS_RTYPE;
    ctrl.dbglvl = ONMETIS_DBGLVL;
  }
  else {
    ctrl.CType = options[OPTION_CTYPE];
    ctrl.IType = options[OPTION_ITYPE];
    ctrl.RType = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }

  ctrl.oflags    = 0;
  ctrl.pfactor   = 0;
  ctrl.nseps     = 5;
  ctrl.optype    = OP_OEMETIS;
  ctrl.CoarsenTo = amin(100, *nvtxs-1);
  ctrl.maxvwgt   = 1.5*tvwgt/ctrl.CoarsenTo;

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, 2);

  /*============================================================
   * Perform the bisection
   *============================================================*/ 
  tpwgts[0] = tvwgt/2;
  tpwgts[1] = tvwgt-tpwgts[0];

  MlevelEdgeBisection(&ctrl, &graph, tpwgts, 1.05);
  ConstructMinCoverSeparator(&ctrl, &graph, 1.05);

  *sepsize = graph.pwgts[2];
  idxcopy(*nvtxs, graph.where, part);

  GKfree((void **)&graph.gdata, &graph.rdata, &graph.label, LTERM);


  FreeWorkSpace(&ctrl, &graph);

}


/*************************************************************************
* This function is the entry point for PWMETIS that accepts exact weights
* for the target partitions
**************************************************************************/
void METIS_mCPartGraphRecursive2(int *nvtxs, int *ncon, idxtype *xadj, idxtype *adjncy, 
       idxtype *vwgt, idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, 
       float *tpwgts, int *options, int *edgecut, idxtype *part)
{
  int i, j;
  GraphType graph;
  CtrlType ctrl;
  float *mytpwgts;
  float avgwgt;

  if (*numflag == 1)
    Change2CNumbering(*nvtxs, xadj, adjncy);

  SetUpGraph(&graph, OP_PMETIS, *nvtxs, *ncon, xadj, adjncy, vwgt, adjwgt, *wgtflag);
  graph.npwgts = NULL;
  mytpwgts = fmalloc(*nparts, "mytpwgts");
  scopy(*nparts, tpwgts, mytpwgts);

  if (options[0] == 0) {  /* Use the default parameters */
    ctrl.CType  = McPMETIS_CTYPE;
    ctrl.IType  = McPMETIS_ITYPE;
    ctrl.RType  = McPMETIS_RTYPE;
    ctrl.dbglvl = McPMETIS_DBGLVL;
  }
  else {
    ctrl.CType  = options[OPTION_CTYPE];
    ctrl.IType  = options[OPTION_ITYPE];
    ctrl.RType  = options[OPTION_RTYPE];
    ctrl.dbglvl = options[OPTION_DBGLVL];
  }
  ctrl.optype = OP_PMETIS;
  ctrl.CoarsenTo = 100;

  ctrl.nmaxvwgt = 1.5/(1.0*ctrl.CoarsenTo);

  InitRandom(options[7]);

  AllocateWorkSpace(&ctrl, &graph, *nparts);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  ASSERT(CheckGraph(&graph));
  *edgecut = MCMlevelRecursiveBisection2(&ctrl, &graph, *nparts, mytpwgts, part, 1.000, 0);

/* 
{
idxtype wgt[2048], minwgt, maxwgt, sumwgt;

printf("nvtxs: %d, nparts: %d, ncon: %d\n", graph.nvtxs, *nparts, *ncon);
for (i=0; i<(*nparts)*(*ncon); i++)
  wgt[i] = 0;
for (i=0; i<graph.nvtxs; i++)
  for (j=0; j<*ncon; j++)
    wgt[part[i]*(*ncon)+j] += vwgt[i*(*ncon)+j];

for (j=0; j<*ncon; j++) {
 minwgt = maxwgt = sumwgt = 0;
 for (i=0; i<(*nparts); i++) {
   minwgt = (wgt[i*(*ncon)+j] < wgt[minwgt*(*ncon)+j]) ? i : minwgt;
   maxwgt = (wgt[i*(*ncon)+j] > wgt[maxwgt*(*ncon)+j]) ? i : maxwgt;
   sumwgt += wgt[i*(*ncon)+j];
 }
 avgwgt = (float)sumwgt / (float)*nparts;
 printf("min: %5d, max: %5d, avg: %5.2f, balance: %6.3f\n", wgt[minwgt*(*ncon)+j], wgt[maxwgt*(*ncon)+j], avgwgt, (float)wgt[maxwgt*(*ncon)+j] / avgwgt);
}
printf("\n");
}
*/

  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(ctrl.dbglvl, DBG_TIME, PrintTimers(&ctrl));

  FreeWorkSpace(&ctrl, &graph);
  GKfree((void **)&mytpwgts, LTERM);

  if (*numflag == 1)
    Change2FNumbering(*nvtxs, xadj, adjncy, part);
}


/*************************************************************************
* This function takes a graph and produces a bisection of it
**************************************************************************/
int MCMlevelRecursiveBisection2(CtrlType *ctrl, GraphType *graph, int nparts,
    float *tpwgts, idxtype *part, float ubfactor, int fpart)
{
  int i, nvtxs, cut;
  float wsum, tpwgts2[2];
  idxtype *label, *where;
  GraphType lgraph, rgraph;

  nvtxs = graph->nvtxs;
  if (nvtxs == 0) 
    return 0;

  /* Determine the weights of the partitions */
  tpwgts2[0] = ssum(nparts/2, tpwgts);
  tpwgts2[1] = 1.0-tpwgts2[0];

  MCMlevelEdgeBisection(ctrl, graph, tpwgts2, ubfactor);
  cut = graph->mincut;

  label = graph->label;
  where = graph->where;
  for (i=0; i<nvtxs; i++)
    part[label[i]] = where[i] + fpart;

  if (nparts > 2) 
    SplitGraphPart(ctrl, graph, &lgraph, &rgraph);

  /* Free the memory of the top level graph */
  GKfree((void **)&graph->gdata, &graph->nvwgt, &graph->rdata, &graph->label, &graph->npwgts, LTERM);

  /* Scale the fractions in the tpwgts according to the true weight */
  wsum = ssum(nparts/2, tpwgts);
  sscale(nparts/2, 1.0/wsum, tpwgts);
  sscale(nparts-nparts/2, 1.0/(1.0-wsum), tpwgts+nparts/2);

  /* Do the recursive call */
  if (nparts > 3) {
    cut += MCMlevelRecursiveBisection2(ctrl, &lgraph, nparts/2, tpwgts, part, ubfactor, fpart);
    cut += MCMlevelRecursiveBisection2(ctrl, &rgraph, nparts-nparts/2, tpwgts+nparts/2, part, ubfactor, fpart+nparts/2);
  }
  else if (nparts == 3) {
    cut += MCMlevelRecursiveBisection2(ctrl, &rgraph, nparts-nparts/2, tpwgts+nparts/2, part, ubfactor, fpart+nparts/2);
    GKfree((void **)&lgraph.gdata, &lgraph.nvwgt, &lgraph.label, LTERM);
  }

  return cut;

}



/*************************************************************************/
/*! This function is the entry point of a node-based separator refinement
    of the nodes with an hmarker[] of 0. */
/*************************************************************************/
void METIS_NodeRefine(int nvtxs, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
           idxtype *adjwgt, idxtype *where, idxtype *hmarker, float ubfactor)
{
  GraphType *graph;
  CtrlType ctrl;

  ctrl.dbglvl    = ONMETIS_DBGLVL;
  ctrl.optype    = OP_ONMETIS;

  graph = CreateGraph();
  SetUpGraph(graph, OP_ONMETIS, nvtxs, 1, xadj, adjncy, vwgt, adjwgt, 3);

  AllocateWorkSpace(&ctrl, graph, 2);

  Allocate2WayNodePartitionMemory(&ctrl, graph);
  idxcopy(nvtxs, where, graph->where);

  Compute2WayNodePartitionParams(&ctrl, graph);

  FM_2WayNodeRefine_OneSidedP(&ctrl, graph, hmarker, ubfactor, 10); 
  /* FM_2WayNodeRefine_TwoSidedP(&ctrl, graph, hmarker, ubfactor, 10); */

  FreeWorkSpace(&ctrl, graph);

  idxcopy(nvtxs, graph->where, where);

  FreeGraph(graph);

}


/*************************************************************************/
/*! This function performs a node-based 1-sided FM refinement that moves
    only nodes whose hmarker[] == -1. It is used by Parmetis. */
/*************************************************************************/
void FM_2WayNodeRefine_OneSidedP(CtrlType *ctrl, GraphType *graph, 
          idxtype *hmarker, float ubfactor, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind, nbad, qsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *swaps, *perm, *inqueue;
  PQueueType parts; 
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;	
  int pass, from, to, limit;
  int badmaxpwgt, mindiff, newdiff;

  ASSERT(graph->mincut == graph->pwgts[2]);

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  vwgt   = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where  = graph->where;
  pwgts  = graph->pwgts;
  rinfo  = graph->nrinfo;

  PQueueInit(ctrl, &parts, nvtxs, ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt));
      
  perm    = idxwspacemalloc(ctrl, nvtxs);
  swaps   = idxwspacemalloc(ctrl, nvtxs);
  mptr    = idxwspacemalloc(ctrl, nvtxs+1);
  mind    = idxwspacemalloc(ctrl, nvtxs);
  inqueue = idxwspacemalloc(ctrl, nvtxs);

  idxset(nvtxs, -1, inqueue);

  badmaxpwgt = (int)(ubfactor*amax(pwgts[0], pwgts[1]));

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions-N1: [%6d %6d] Nv-Nb[%6d %6d] MaxPwgt[%6d]. ISep: %6d\n", 
        pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, badmaxpwgt, graph->mincut));

  to = (pwgts[0] < pwgts[1] ? 1 : 0);
  for (pass=0; pass<npasses; pass++) {
    from = to; 
    to   = (from+1)%2;

    PQueueReset(&parts);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      if (hmarker[i] == -1 || hmarker[i] == to) {
        PQueueInsert(&parts, i, vwgt[i]-rinfo[i].edegrees[from]);
        inqueue[i] = pass;
      }
    }
    qsize = parts.nnodes;

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = nbnd;

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = nbad = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      if ((higain = PQueueGetMax(&parts)) == -1) 
        break;

      inqueue[higain] = -1;

      ASSERT(bndptr[higain] != -1);

      if (pwgts[to]+vwgt[higain] > badmaxpwgt) { /* Skip this vertex */
        if (nbad++ > limit) 
          break; 
        else {
          nswaps--;
          continue;  
        }
      }

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[from]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[from]-rinfo[higain].edegrees[from]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut      = pwgts[2];
        mincutorder = nswaps;
        mindiff     = newdiff;
        nbad        = 0;
      }
      else {
        if (nbad++ > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[from]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      swaps[nswaps] = higain;  


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          rinfo[k].edegrees[to] += vwgt[higain];
        }
        else if (where[k] == from) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k]      = 2;
          pwgts[from]  -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2) 
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[from];
              rinfo[kk].edegrees[from] -= vwgt[k];

              /* Update the gain of this node if it was skipped */
              if (inqueue[kk] == pass)
                PQueueUpdateUp(&parts, kk, oldgain, oldgain+vwgt[k]); 
            }
          }

          /* Insert the new vertex into the priority queue. Safe due to one-sided moves */
          if (hmarker[k] == -1 || hmarker[k] == to) {
            PQueueInsert(&parts, k, vwgt[k]-edegrees[from]);
            inqueue[k] = pass;
          }
        }
      }
      mptr[nswaps+1] = nmind;


      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] \t[%5d %5d %5d] [%3d %2d]\n", 
                       higain, to, (vwgt[higain]-rinfo[higain].edegrees[from]), vwgt[higain], pwgts[0], pwgts[1], pwgts[2], nswaps, limit));

    }


    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));
      ASSERT(where[higain] == to);

      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) 
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = from;
        INC_DEC(pwgts[from], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2) 
            rinfo[kk].edegrees[from] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d, QSIZE: %6d\n", 
          mincut, mincutorder, pwgts[0], pwgts[1], nbnd, qsize));

    graph->mincut = mincut;
    graph->nbnd   = nbnd;

    if (pass%2 == 1 && (mincutorder == -1 || mincut >= initcut))
      break;
  }

  PQueueFree(ctrl, &parts);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************/
/*! This function performs a node-based (two-sided) FM refinement that 
    moves only nodes whose hmarker[] == -1. It is used by Parmetis. */
/*************************************************************************/
void FM_2WayNodeRefine_TwoSidedP(CtrlType *ctrl, GraphType *graph, 
          idxtype *hmarker, float ubfactor, int npasses)
{
  int i, ii, j, k, jj, kk, nvtxs, nbnd, nswaps, nmind;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *edegrees, *bndind, *bndptr;
  idxtype *mptr, *mind, *moved, *swaps, *perm;
  PQueueType parts[2]; 
  NRInfoType *rinfo;
  int higain, oldgain, mincut, initcut, mincutorder;	
  int pass, to, other, limit;
  int badmaxpwgt, mindiff, newdiff;
  int u[2], g[2];

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  vwgt   = graph->vwgt;

  bndind = graph->bndind;
  bndptr = graph->bndptr;
  where  = graph->where;
  pwgts  = graph->pwgts;
  rinfo  = graph->nrinfo;


  i = ComputeMaxNodeGain(nvtxs, xadj, adjncy, vwgt);
  PQueueInit(ctrl, &parts[0], nvtxs, i);
  PQueueInit(ctrl, &parts[1], nvtxs, i);

  moved    = idxwspacemalloc(ctrl, nvtxs);
  swaps    = idxwspacemalloc(ctrl, nvtxs);
  mptr     = idxwspacemalloc(ctrl, nvtxs+1);
  mind     = idxwspacemalloc(ctrl, nvtxs);
  perm     = idxwspacemalloc(ctrl, nvtxs);

  IFSET(ctrl->dbglvl, DBG_REFINE,
    printf("Partitions: [%6d %6d] Nv-Nb[%6d %6d]. ISep: %6d\n", pwgts[0], pwgts[1], graph->nvtxs, graph->nbnd, graph->mincut));

  badmaxpwgt = (int)(ubfactor*amax(pwgts[0], pwgts[1]));

  for (pass=0; pass<npasses; pass++) {
    idxset(nvtxs, -1, moved);
    PQueueReset(&parts[0]);
    PQueueReset(&parts[1]);

    mincutorder = -1;
    initcut = mincut = graph->mincut;
    nbnd = graph->nbnd;

    RandomPermute(nbnd, perm, 1);
    for (ii=0; ii<nbnd; ii++) {
      i = bndind[perm[ii]];
      ASSERT(where[i] == 2);
      if (hmarker[i] == -1) {
        PQueueInsert(&parts[0], i, vwgt[i]-rinfo[i].edegrees[1]);
        PQueueInsert(&parts[1], i, vwgt[i]-rinfo[i].edegrees[0]);
        moved[i] = -5;
      }
      else if (hmarker[i] != 2) {
        PQueueInsert(&parts[hmarker[i]], i, vwgt[i]-rinfo[i].edegrees[(hmarker[i]+1)%2]);
        moved[i] = -(10+hmarker[i]);
      }
    }

    ASSERT(CheckNodeBnd(graph, nbnd));
    ASSERT(CheckNodePartitionParams(graph));

    limit = nbnd;

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    mptr[0] = nmind = 0;
    mindiff = abs(pwgts[0]-pwgts[1]);
    to = (pwgts[0] < pwgts[1] ? 0 : 1);
    for (nswaps=0; nswaps<nvtxs; nswaps++) {
      u[0] = PQueueSeeMax(&parts[0]);  
      u[1] = PQueueSeeMax(&parts[1]);
      if (u[0] != -1 && u[1] != -1) {
        g[0] = vwgt[u[0]]-rinfo[u[0]].edegrees[1];
        g[1] = vwgt[u[1]]-rinfo[u[1]].edegrees[0];

        to = (g[0] > g[1] ? 0 : (g[0] < g[1] ? 1 : pass%2)); 

        if (pwgts[to]+vwgt[u[to]] > badmaxpwgt) 
          to = (to+1)%2;
      }
      else if (u[0] == -1 && u[1] == -1) {
        break;
      }
      else if (u[0] != -1 && pwgts[0]+vwgt[u[0]] <= badmaxpwgt) {
        to = 0;
      }
      else if (u[1] != -1 && pwgts[1]+vwgt[u[1]] <= badmaxpwgt) {
        to = 1;
      }
      else
        break;

      other = (to+1)%2;

      higain = PQueueGetMax(&parts[to]);

      /* Delete its matching entry in the other queue */
      if (moved[higain] == -5) 
        PQueueDelete(&parts[other], higain, vwgt[higain]-rinfo[higain].edegrees[to]);

      ASSERT(bndptr[higain] != -1);

      pwgts[2] -= (vwgt[higain]-rinfo[higain].edegrees[other]);

      newdiff = abs(pwgts[to]+vwgt[higain] - (pwgts[other]-rinfo[higain].edegrees[other]));
      if (pwgts[2] < mincut || (pwgts[2] == mincut && newdiff < mindiff)) {
        mincut      = pwgts[2];
        mincutorder = nswaps;
        mindiff     = newdiff;
      }
      else {
        if (nswaps - mincutorder > limit) {
          pwgts[2] += (vwgt[higain]-rinfo[higain].edegrees[other]);
          break; /* No further improvement, break out */
        }
      }

      BNDDelete(nbnd, bndind, bndptr, higain);
      pwgts[to] += vwgt[higain];
      where[higain] = to;
      moved[higain] = nswaps;
      swaps[nswaps] = higain;  


      /**********************************************************
      * Update the degrees of the affected nodes
      ***********************************************************/
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) { /* For the in-separator vertices modify their edegree[to] */
          oldgain = vwgt[k]-rinfo[k].edegrees[to];
          rinfo[k].edegrees[to] += vwgt[higain];
          if (moved[k] == -5 || moved[k] == -(10+other)) 
            PQueueUpdate(&parts[other], k, oldgain, oldgain-vwgt[higain]);
        }
        else if (where[k] == other) { /* This vertex is pulled into the separator */
          ASSERTP(bndptr[k] == -1, ("%d %d %d\n", k, bndptr[k], where[k]));
          BNDInsert(nbnd, bndind, bndptr, k);

          mind[nmind++] = k;  /* Keep track for rollback */
          where[k] = 2;
          pwgts[other] -= vwgt[k];

          edegrees = rinfo[k].edegrees;
          edegrees[0] = edegrees[1] = 0;
          for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
            kk = adjncy[jj];
            if (where[kk] != 2) 
              edegrees[where[kk]] += vwgt[kk];
            else {
              oldgain = vwgt[kk]-rinfo[kk].edegrees[other];
              rinfo[kk].edegrees[other] -= vwgt[k];
              if (moved[kk] == -5 || moved[kk] == -(10+to))
                PQueueUpdate(&parts[to], kk, oldgain, oldgain+vwgt[k]);
            }
          }

          /* Insert the new vertex into the priority queue (if it has not been moved). */
          if (moved[k] == -1 && (hmarker[k] == -1 || hmarker[k] == to)) {
            PQueueInsert(&parts[to], k, vwgt[k]-edegrees[other]);
            moved[k] = -(10+to);
          }
#ifdef FULLMOVES  /* this does not work as well as the above partial one */
          if (moved[k] == -1) {
            if (hmarker[k] == -1) {
              PQueueInsert(&parts[0], k, vwgt[k]-edegrees[1]);
              PQueueInsert(&parts[1], k, vwgt[k]-edegrees[0]);
              moved[k] = -5;
            }
            else if (hmarker[k] != 2) {
              PQueueInsert(&parts[hmarker[k]], k, vwgt[k]-edegrees[(hmarker[k]+1)%2]);
              moved[k] = -(10+hmarker[k]);
            }
          }
#endif
        }
      }
      mptr[nswaps+1] = nmind;

      IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moved %6d to %3d, Gain: %5d [%5d] [%4d %4d] \t[%5d %5d %5d]\n", higain, to, g[to], g[other], vwgt[u[to]], vwgt[u[other]], pwgts[0], pwgts[1], pwgts[2]));

    }


    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    for (nswaps--; nswaps>mincutorder; nswaps--) {
      higain = swaps[nswaps];

      ASSERT(CheckNodePartitionParams(graph));

      to = where[higain];
      other = (to+1)%2;
      INC_DEC(pwgts[2], pwgts[to], vwgt[higain]);
      where[higain] = 2;
      BNDInsert(nbnd, bndind, bndptr, higain);

      edegrees = rinfo[higain].edegrees;
      edegrees[0] = edegrees[1] = 0;
      for (j=xadj[higain]; j<xadj[higain+1]; j++) {
        k = adjncy[j];
        if (where[k] == 2) 
          rinfo[k].edegrees[to] -= vwgt[higain];
        else
          edegrees[where[k]] += vwgt[k];
      }

      /* Push nodes out of the separator */
      for (j=mptr[nswaps]; j<mptr[nswaps+1]; j++) {
        k = mind[j];
        ASSERT(where[k] == 2);
        where[k] = other;
        INC_DEC(pwgts[other], pwgts[2], vwgt[k]);
        BNDDelete(nbnd, bndind, bndptr, k);
        for (jj=xadj[k]; jj<xadj[k+1]; jj++) {
          kk = adjncy[jj];
          if (where[kk] == 2) 
            rinfo[kk].edegrees[other] += vwgt[k];
        }
      }
    }

    ASSERT(mincut == pwgts[2]);

    IFSET(ctrl->dbglvl, DBG_REFINE,
      printf("\tMinimum sep: %6d at %5d, PWGTS: [%6d %6d], NBND: %6d\n", mincut, mincutorder, pwgts[0], pwgts[1], nbnd));

    graph->mincut = mincut;
    graph->nbnd = nbnd;

    if (mincutorder == -1 || mincut >= initcut)
      break;
  }

  PQueueFree(ctrl, &parts[0]);
  PQueueFree(ctrl, &parts[1]);

  idxwspacefree(ctrl, nvtxs+1);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}

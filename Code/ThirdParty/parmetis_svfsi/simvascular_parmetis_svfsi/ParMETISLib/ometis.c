/*!
 * Copyright 1997, Regents of the University of Minnesota
 *
 * \file
 * \brief This is the entry point of parallel ordering routines
 *
 * \date Started 8/1/2008
 * \author George Karypis
 * \version\verbatim $Id: ometis.c 9716 2011-04-06 15:50:54Z karypis $ \endverbatime
 *
 */

#include <parmetislib.h>


/***********************************************************************************/
/*! This function is the entry point of the parallel ordering algorithm. 
    It simply translates the arguments to the tunable version. */
/***********************************************************************************/
void ParMETIS_V3_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
         idxtype *numflag, int *options, idxtype *order, idxtype *sizes,
         MPI_Comm *comm)
{
  int seed   = (options != NULL && options[0] != 0 ? options[PMV3_OPTION_SEED] : -1);
  int dbglvl = (options != NULL && options[0] != 0 ? options[PMV3_OPTION_DBGLVL] : -1);

  ParMETIS_V32_NodeND(vtxdist, xadj, adjncy, 
      /*vwgt=*/NULL, 
      numflag, 
      /*mtype=*/NULL, 
      /*rtype=*/NULL, 
      /*p_nseps=*/NULL, 
      /*s_nseps=*/NULL, 
      /*ubfrac=*/NULL,
      /*seed=*/(options==NULL | options[0] == 0 ? NULL : &seed),
      /*dbglvl=*/(options==NULL | options[0] == 0 ? NULL : &dbglvl),
      order, sizes, comm);
}


/***********************************************************************************/
/*! This function is the entry point of the parallel ordering algorithm using the
    old API */
/***********************************************************************************/
void PAROMETIS(idxtype *vtxdist, idxtype *xadj, idxtype *vwgt, idxtype *adjncy, 
         idxtype *adjwgt, idxtype *order, idxtype *sizes, int *options, MPI_Comm comm)
{
  int numflag, newoptions[5];

  newoptions[0] = 1;
  newoptions[PMV3_OPTION_DBGLVL] = options[4];
  newoptions[PMV3_OPTION_SEED] = GLOBAL_SEED;

  numflag = options[3];

  ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, &numflag, newoptions, order, sizes, &comm);

  options[0] = -1;

}


/***********************************************************************************/
/*! This function is the entry point of the tunable parallel ordering algorithm. 
    This is the main ordering algorithm and implements a multilevel nested 
    dissection ordering approach. 
*/
/***********************************************************************************/
void ParMETIS_V32_NodeND(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
              idxtype *numflag, int *mtype, int *rtype, int *p_nseps, int *s_nseps, 
              float *ubfrac, int *seed, int *dbglvl, idxtype *order, idxtype *sizes, 
              MPI_Comm *comm)
{
  int i, j;
  int ltvwgts[MAXNCON];
  int npes, mype, wgtflag;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  idxtype *morder;
  int minnvtxs, dbglvl_original;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  /* Deal with poor vertex distributions */
  ctrl.comm = *comm;
  if (GlobalSEMin(&ctrl, vtxdist[mype+1]-vtxdist[mype]) < 1) {
    if (mype == 0)
      printf("Error: Poor vertex distribution (processor with no vertices).\n");
    return;
  }

#ifdef XXX
  /* Increase all weights by one to eliminate potentially zero weight vertices */
  if (vwgt) {
    for (i=0; i<vtxdist[mype+1]-vtxdist[mype]; i++)
      vwgt[i]++;
  }
#endif


  /* Renumber the adjacency list if appropriate */ 
  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 1);


  SetUpComm(&ctrl, *comm);

  dbglvl_original = (dbglvl == NULL ? 0 : *dbglvl);

  /*=======================================================================*/
  /*! Compute the initial k-way partitioning */
  /*=======================================================================*/
  ctrl.nparts      = 5*npes;
  ctrl.partType    = STATIC_PARTITION;
  ctrl.tpwgts      = fsmalloc(ctrl.nparts, 1.0/(float)(ctrl.nparts), "tpwgts");
  ctrl.ubvec[0]    = 1.03;
  ctrl.CoarsenTo   = amin(vtxdist[npes]+1, 200*amax(npes, ctrl.nparts));
  ctrl.ps_relation = -1;
  ctrl.dbglvl      = 0;

  ctrl.seed        = (seed == NULL ? GLOBAL_SEED : *seed);
  ctrl.seed        = (ctrl.seed == 0 ? mype : ctrl.seed*mype);
  ctrl.sync        = GlobalSEMax(&ctrl, ctrl.seed);

  wgtflag = 0;
  graph = Mc_SetUpGraph(&ctrl, 1, vtxdist, xadj, NULL, adjncy, NULL, &wgtflag);

  AllocateWSpace(&ctrl, graph, &wspace);

  IFSET(dbglvl_original, DBG_TIME, InitTimers(&ctrl));
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, starttimer(ctrl.TotalTmr));

  Mc_Global_Partition(&ctrl, graph, &wspace);

  /* Collapse the number of partitions to be from 0..npes-1 */
  for (i=0; i<graph->nvtxs; i++)
    graph->where[i] = graph->where[i]%npes;
  ctrl.nparts = npes;

  /* Put back the real vertex weights */
  if (vwgt) {
    GKfree((void **)&graph->vwgt, LTERM);
    graph->vwgt = vwgt;
    wgtflag = 2;
  }


  /*=======================================================================*/
  /*! Move the graph according to the partitioning */
  /*=======================================================================*/
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, starttimer(ctrl.MoveTmr));

  mgraph = Mc_MoveGraph(&ctrl, graph, &wspace);

  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, stoptimer(ctrl.MoveTmr));


  /*=======================================================================*/
  /*! Now compute an ordering of the moved graph */
  /*=======================================================================*/
  AdjustWSpace(&ctrl, mgraph, &wspace);

  ctrl.partType  = ORDER_PARTITION;
  ctrl.mtype     = (mtype  == NULL ? PARMETIS_MTYPE_GLOBAL  : *mtype);
  ctrl.rtype     = (rtype  == NULL ? PARMETIS_SRTYPE_2PHASE : *rtype);
  ctrl.p_nseps   = (p_nseps  == NULL ? 1 : *p_nseps);
  ctrl.s_nseps   = (s_nseps  == NULL ? 1 : *s_nseps);
  ctrl.ubfrac    = (ubfrac == NULL ? ORDER_UNBALANCE_FRACTION : *ubfrac);
  ctrl.dbglvl    = dbglvl_original;
  ctrl.ipart     = ISEP_NODE;
  ctrl.CoarsenTo = amin(graph->gnvtxs-1,
                        amax(1500*npes, graph->gnvtxs/(5*NUM_INIT_MSECTIONS*npes)));

  /* compute tvwgts */
  for (j=0; j<mgraph->ncon; j++)
    ltvwgts[j] = 0;

  for (i=0; i<mgraph->nvtxs; i++)
    for (j=0; j<mgraph->ncon; j++)
      ltvwgts[j] += mgraph->vwgt[i*mgraph->ncon+j];

  for (j=0; j<mgraph->ncon; j++)
    ctrl.tvwgts[j] = GlobalSESum(&ctrl, ltvwgts[j]);

  mgraph->nvwgt = fmalloc(mgraph->nvtxs*mgraph->ncon, "mgraph->nvwgt");
  for (i=0; i<mgraph->nvtxs; i++)
    for (j=0; j<mgraph->ncon; j++)
      mgraph->nvwgt[i*mgraph->ncon+j] = (float)(mgraph->vwgt[i*mgraph->ncon+j]) / (float)(ctrl.tvwgts[j]);


  morder = idxmalloc(mgraph->nvtxs, "PAROMETIS: morder");
  MultilevelOrder(&ctrl, mgraph, morder, sizes, &wspace);

  /* Invert the ordering back to the original graph */
  ProjectInfoBack(&ctrl, graph, order, morder, &wspace);

  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(dbglvl_original, DBG_TIME, stoptimer(ctrl.TotalTmr));
  IFSET(dbglvl_original, DBG_TIME, PrintTimingInfo(&ctrl));
  IFSET(dbglvl_original, DBG_TIME, MPI_Barrier(ctrl.gcomm));

  GKfree((void **)&ctrl.tpwgts, &morder, LTERM);
  FreeGraph(mgraph);
  FreeInitialGraphAndRemap(graph, wgtflag, 1);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  /* If required, restore the graph numbering */
  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, order, npes, mype, 0);

#ifdef XXX
  /* Decrease the earlier increased weights */
  if (vwgt) {
    for (i=0; i<vtxdist[mype+1]-vtxdist[mype]; i++)
      vwgt[i]--;
  }
#endif

}


/*********************************************************************************/
/*!
  This is the top level ordering routine. 
  \param order is the computed ordering.
  \param sizes is the 2*nparts array that will store the sizes of each subdomains 
               and the sizes of the separators at each level. Note that the 
               top-level separator is stores at \c sizes[2*nparts-2].
*/
/*********************************************************************************/
void MultilevelOrder(CtrlType *ctrl, GraphType *graph, idxtype *order, idxtype *sizes, 
         WorkSpaceType *wspace)
{
  int i, nparts, nvtxs, npes;
  idxtype *perm, *lastnode, *morder, *porder;
  GraphType *mgraph;

  nvtxs = graph->nvtxs;

  npes = 1<<log2Int(ctrl->npes); /* # of nested dissection levels = floor(log_2(npes)) */

  perm     = idxmalloc(nvtxs, "MultilevelOrder: perm");
  lastnode = idxsmalloc(4*npes, -1, "MultilevelOrder: lastnode");

  for (i=0; i<nvtxs; i++) 
    perm[i] = i;
  lastnode[2] = graph->gnvtxs;

  idxset(nvtxs, -1, order);

  /* This is used as a pointer to the end of the sizes[] array (i.e., >=nparts)
     that has not yet been filled in so that the separator sizes of the succesive
     levels will be stored correctly. It is used in LabelSeparatos() */
  sizes[0] = 2*npes-1;

  graph->where = idxsmalloc(nvtxs, 0, "MultilevelOrder: graph->where");

  for (nparts=2; nparts<=npes; nparts*=2) {
    ctrl->nparts = nparts;

    Order_Partition_Multiple(ctrl, graph, wspace);

    LabelSeparators(ctrl, graph, lastnode, perm, order, sizes, wspace);

    CompactGraph(ctrl, graph, perm, wspace);

    if (ctrl->CoarsenTo < 100*nparts) {
      ctrl->CoarsenTo = 1.5*ctrl->CoarsenTo;
    }
    ctrl->CoarsenTo = amin(ctrl->CoarsenTo, graph->gnvtxs-1);
  }


  /*-----------------------------------------------------------------
   / Move the graph so that each processor gets its partition 
   -----------------------------------------------------------------*/
  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MoveTmr));

  SetUp(ctrl, graph, wspace);
  graph->ncon = 1; /* needed for Mc_MoveGraph */
  mgraph = Mc_MoveGraph(ctrl, graph, wspace);

  /* Fill in the sizes[] array for the local part. Just the vtxdist of the mgraph */
  for (i=0; i<npes; i++)
    sizes[i] = mgraph->vtxdist[i+1]-mgraph->vtxdist[i];

  porder = idxmalloc(graph->nvtxs, "MultilevelOrder: porder");
  morder = idxmalloc(mgraph->nvtxs, "MultilevelOrder: morder");

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MoveTmr));

  /* Find the local ordering */
  if (ctrl->mype < npes)
    LocalNDOrder(ctrl, mgraph, morder, lastnode[2*(npes+ctrl->mype)]-mgraph->nvtxs, wspace);

  /* Project the ordering back to the before-move graph */
  ProjectInfoBack(ctrl, graph, porder, morder, wspace);

  /* Copy the ordering from porder to order using perm */
  for (i=0; i<graph->nvtxs; i++) {
    ASSERT(ctrl, order[perm[i]] == -1);
    order[perm[i]] = porder[i];
  }


  FreeGraph(mgraph);
  GKfree((void **)&perm, (void **)&lastnode, (void **)&porder, (void **)&morder, LTERM);

  /* PrintVector(ctrl, 2*npes-1, 0, sizes, "SIZES"); */
}


/**************************************************************************/
/*! This is the top-level driver of the multiple multisection ordering
    code. */
/***************************************************************************/
void Order_Partition_Multiple(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, sid, iter, nvtxs, nparts, nlevels;
  idxtype *xadj, *adjncy, *where, *gpwgts, *imap;
  idxtype *bestseps, *bestwhere, *origwhere;

  SetUp(ctrl, graph, wspace);

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  bestseps  = idxsmalloc(2*nparts, -1, "Order_Partition_Multiple: bestseps");
  bestwhere = idxmalloc(nvtxs+graph->nrecv, "Order_Partition_Multiple: bestwhere");

  origwhere = graph->where;

  for (nlevels=-1, iter=0; iter<ctrl->p_nseps; iter++) {
    graph->where = idxmalloc(nvtxs, "Order_Partition_Multiple: where");
    idxcopy(nvtxs, origwhere, graph->where);

    Order_Partition(ctrl, graph, wspace, &nlevels, 0);

    where  = graph->where;
    gpwgts = graph->gpwgts;
    /* Update the where[] vectors of the subdomains that improved */
    for (i=0; i<nvtxs; i++) {
      sid = (where[i] < nparts ? nparts + where[i] - (where[i]%2) : where[i]);
      if (iter == 0 || bestseps[sid] > gpwgts[sid])
        bestwhere[i] = where[i];
    }
    /* Update the size of the separators for those improved subdomains */
    for (i=0; i<nparts; i+=2) {
      sid = nparts+i;
      if (iter == 0 || bestseps[sid] > gpwgts[sid]) 
        bestseps[sid] = gpwgts[sid];
    }

    /* free all the memory allocated for coarsening/refinement, but keep
       the setup fields so that they will not be re-computed */
    FreeNonGraphNonSetupFields(graph);
  }

  graph->where = bestwhere;
  AllocateNodePartitionParams(ctrl, graph, wspace);
  ComputeNodePartitionParams(ctrl, graph, wspace);

  for (i=0; i<nparts; i+=2) 
    ASSERT(ctrl, bestseps[nparts+i] == graph->gpwgts[nparts+i]);

  GKfree((void **)&bestseps, &origwhere, LTERM);

  /* PrintVector(ctrl, 2*nparts-1, 0, bestseps, "bestseps"); */

}


/**************************************************************************/
/*! The driver of the multilvelel separator finding algorithm */
/**************************************************************************/
void Order_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace,
         int *nlevels, int clevel)
{

  SetUp(ctrl, graph, wspace);
  graph->ncon = 1;

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d]\n",
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));

  if ((*nlevels != -1 && *nlevels == clevel) ||
      (*nlevels == -1 && 
       ((graph->gnvtxs < 1.66*ctrl->CoarsenTo) || 
        (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)))) {
    /* set the nlevels to where coarsening stopped */
    *nlevels = clevel;

    /* Compute the initial npart-way multisection */
    InitMultisection(ctrl, graph, wspace);

    if (graph->finer == NULL) { /* Do that only if no-coarsening took place */
      AllocateNodePartitionParams(ctrl, graph, wspace);
      ComputeNodePartitionParams(ctrl, graph, wspace);
      switch (ctrl->rtype) {
        case PARMETIS_SRTYPE_GREEDY:
          KWayNodeRefine_Greedy(ctrl, graph, wspace, NGR_PASSES, ctrl->ubfrac);
          break;
        case PARMETIS_SRTYPE_2PHASE:
          KWayNodeRefine2Phase(ctrl, graph, wspace, NGR_PASSES, ctrl->ubfrac);
          break;
        default:
          errexit("Unknown rtype of %d\n", ctrl->rtype);
      }
    }
  }
  else { /* Coarsen it and then partition it */
    switch (ctrl->mtype) {
      case PARMETIS_MTYPE_LOCAL:
        Match_Local(ctrl, graph, wspace);
        break;
      case PARMETIS_MTYPE_GLOBAL:
        Match_Global(ctrl, graph, wspace);
        break;
      default:
        errexit("Unknown mtype of %d\n", ctrl->mtype);
    }

    Order_Partition(ctrl, graph->coarser, wspace, nlevels, clevel+1);

    Mc_ProjectPartition(ctrl, graph, wspace);
    AllocateNodePartitionParams(ctrl, graph, wspace);
    ComputeNodePartitionParams(ctrl, graph, wspace);

    switch (ctrl->rtype) {
      case PARMETIS_SRTYPE_GREEDY:
        KWayNodeRefine_Greedy(ctrl, graph, wspace, NGR_PASSES, ctrl->ubfrac);
        break;
      case PARMETIS_SRTYPE_2PHASE:
        KWayNodeRefine2Phase(ctrl, graph, wspace, NGR_PASSES, ctrl->ubfrac);
        break;
      default:
        errexit("Unknown rtype of %d\n", ctrl->rtype);
    }
  }
}



/*********************************************************************************/
/*! This function is used to assign labels to the nodes in the separators. 
    It uses the appropriate entry in the lastnode array to select label boundaries
    and adjusts it for the next level. */
/*********************************************************************************/
void LabelSeparators(CtrlType *ctrl, GraphType *graph, idxtype *lastnode, idxtype
    *perm, idxtype *order, idxtype *sizes, WorkSpaceType *wspace) 
{ 
  int i, nvtxs, nparts, sid; idxtype *where, *lpwgts, *gpwgts, *sizescan;

  nparts = ctrl->nparts;

  nvtxs  = graph->nvtxs;
  where  = graph->where;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  if (ctrl->dbglvl&DBG_INFO) { 
    if (ctrl->mype == 0) {
      printf("SepWgts:  ");
      for (i=0; i<nparts; i+=2)
        printf(" %d [%d %d]", gpwgts[nparts+i], gpwgts[i], gpwgts[i+1]);
      printf("\n");
    }
    MPI_Barrier(ctrl->comm);
  }

  /* Compute the local size of the separator. This is required in case the 
     graph has vertex weights */
  idxset(2*nparts, 0, lpwgts);
  for (i=0; i<nvtxs; i++) 
    lpwgts[where[i]]++;

  sizescan = idxmalloc(2*nparts, "LabelSeparators: sizescan");

  /* Perform a Prefix scan of the separator sizes to determine the boundaries */
  MPI_Scan((void *)lpwgts, (void *)sizescan, 2*nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
  MPI_Allreduce((void *)lpwgts, (void *)gpwgts, 2*nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

#ifdef DEBUG_ORDER
  PrintVector(ctrl, 2*nparts, 0, lpwgts, "Lpwgts");
  PrintVector(ctrl, 2*nparts, 0, sizescan, "SizeScan");
  PrintVector(ctrl, 2*nparts, 0, lastnode, "LastNode");
#endif

  /* Fillin the sizes[] array. See the comment on MultilevelOrder() on the 
     purpose of the sizes[0] value. */
  for (i=nparts-2; i>=0; i-=2) 
    sizes[--sizes[0]] = gpwgts[nparts+i];

  if (ctrl->dbglvl&DBG_INFO) { 
    if (ctrl->mype == 0) {
      printf("SepSizes: ");
      for (i=0; i<nparts; i+=2)
        printf(" %d [%d %d]", gpwgts[nparts+i], gpwgts[i], gpwgts[i+1]);
      printf("\n");
    }
    MPI_Barrier(ctrl->comm);
  }

  for (i=0; i<2*nparts; i++)
    sizescan[i] -= lpwgts[i];

  /* Assign the order[] values to the separator nodes */
  for (i=0; i<nvtxs; i++) {
    if (where[i] >= nparts) {
      sid = where[i];
      sizescan[sid]++;
      ASSERT(ctrl, order[perm[i]] == -1);
      order[perm[i]] = lastnode[sid] - sizescan[sid];
      /*myprintf(ctrl, "order[%d] = %d, %d\n", perm[i], order[perm[i]], sid); */
    }
  }

  /* Update lastnode array */
  idxcopy(2*nparts, lastnode, sizescan);
  for (i=0; i<nparts; i+=2) {
    lastnode[2*nparts+2*i]     = sizescan[nparts+i]-gpwgts[nparts+i]-gpwgts[i+1];
    lastnode[2*nparts+2*(i+1)] = sizescan[nparts+i]-gpwgts[nparts+i];
    /*myprintf(ctrl, "lastnode: %d %d\n", lastnode[2*nparts+2*i], * lastnode[2*nparts+2*(i+1)]);*/
  }

  GKfree((void **)&sizescan, LTERM);

}



/*************************************************************************
* This function compacts a graph by removing the vertex separator
**************************************************************************/
void CompactGraph(CtrlType *ctrl, GraphType *graph, idxtype *perm, 
         WorkSpaceType *wspace)
{
  int i, j, l, nvtxs, cnvtxs, cfirstvtx, nparts, npes; 
  idxtype *xadj, *adjncy, *adjwgt, *vtxdist, *where;
  idxtype *cmap, *cvtxdist, *newwhere;

  nparts = ctrl->nparts;
  npes   = ctrl->npes;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  if (graph->cmap == NULL)
    graph->cmap = idxmalloc(nvtxs+graph->nrecv, "CompactGraph: cmap");
  cmap = graph->cmap;

  vtxdist = graph->vtxdist;

  /*************************************************************
  * Construct the cvtxdist of the contracted graph. Uses the fact
  * that lpwgts stores the local non separator vertices.
  **************************************************************/
  cvtxdist = wspace->pv1;
  cnvtxs   = idxsum(nparts, graph->lpwgts);

  MPI_Allgather((void *)&cnvtxs, 1, IDX_DATATYPE, (void *)cvtxdist, 1, IDX_DATATYPE, 
       ctrl->comm);
  MAKECSR(i, npes, cvtxdist);

#ifdef DEBUG_ORDER
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[ctrl->mype];

  /* Create the cmap of what you know so far locally */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (where[i] < nparts) {
      perm[cnvtxs] = perm[i];
      cmap[i] = cfirstvtx + cnvtxs++;
    }
  }

  CommInterfaceData(ctrl, graph, cmap, wspace->indices, cmap+nvtxs);


  /*************************************************************
  * Finally, compact the graph
  **************************************************************/
  newwhere = idxmalloc(cnvtxs, "CompactGraph: newwhere");
  cnvtxs = l = 0;
  for (i=0; i<nvtxs; i++) {
    if (where[i] < nparts) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        ASSERT(ctrl, where[i] == where[adjncy[j]] || where[adjncy[j]] >= nparts);
        if (where[i] == where[adjncy[j]]) {
          adjncy[l]   = cmap[adjncy[j]];
          adjwgt[l++] = adjwgt[j];
        }
      }

      xadj[cnvtxs] = l;
      graph->vwgt[cnvtxs] = graph->vwgt[i];
      newwhere[cnvtxs]    = where[i];
      cnvtxs++;
    }
  }
  SHIFTCSR(i, cnvtxs, xadj);

  GKfree((void **)&graph->match, (void **)&graph->cmap, (void **)&graph->lperm, 
         (void **)&graph->where, (void **)&graph->label, (void **)&graph->rinfo,
         (void **)&graph->nrinfo, (void **)&graph->lpwgts, (void **)&graph->gpwgts, 
         (void **)&graph->sepind, (void **)&graph->peind,
         (void **)&graph->sendptr, (void **)&graph->sendind, 
         (void **)&graph->recvptr, (void **)&graph->recvind, 
         (void **)&graph->imap, (void **)&graph->rlens, (void **)&graph->slens, 
         (void **)&graph->rcand, (void **)&graph->pexadj, 
         (void **)&graph->peadjncy, (void **)&graph->peadjloc, LTERM);
 
  graph->nvtxs  = cnvtxs;
  graph->nedges = l;
  graph->gnvtxs = cvtxdist[npes];
  graph->where  = newwhere;
  idxcopy(npes+1, cvtxdist, graph->vtxdist);

  /*
  {
    int i, j, k;
    int *mylpwgts;

    mylpwgts = idxsmalloc(nparts, 0, "mylpwgts");
    for (i=0; i<cnvtxs; i++) 
      mylpwgts[newwhere[i]]++;

    PrintVector(ctrl, nparts, 0, mylpwgts, "mylpwgts");

    GKfree((void **)&mylpwgts, LTERM);
  }
  */
}


/*************************************************************************/
/*! This function orders the locally stored graph using MMD. The vertices 
    will be ordered from firstnode onwards. */
/*************************************************************************/
void LocalNDOrder(CtrlType *ctrl, GraphType *graph, idxtype *order, int firstnode, WorkSpaceType *wspace)
{
  int i, j, nvtxs, firstvtx, lastvtx;
  idxtype *xadj, *adjncy;
  idxtype *perm, *iperm;
  int numflag=0, options[10];

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->SerialTmr));

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;

  firstvtx = graph->vtxdist[ctrl->mype];
  lastvtx  = graph->vtxdist[ctrl->mype+1];

  /* Relabel the vertices so that they are in local index space */
  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ASSERT(ctrl, adjncy[j]>=firstvtx && adjncy[j]<lastvtx);
      adjncy[j] -= firstvtx;
    }
  }

  ASSERT(ctrl, 2*(nvtxs+5) < wspace->maxcore);

  perm = wspace->core;
  iperm = perm + nvtxs + 5;

  options[0] = 1;
  options[1] = 3;              /* ctype = SHEM */
  options[2] = 1;              /* itype = GGPKL */
  options[3] = 2;              /* rtype = sep1sided */
  options[4] = 0;              /* dbglvl */
  options[5] = 1;              /* oflags = compress */
  options[6] = -1;             /* pfactor */
  options[7] = ctrl->s_nseps;  /* nseps */

  if (graph->vwgt) 
    METIS_NodeWND(&nvtxs, xadj, adjncy, graph->vwgt, &numflag, options, perm, iperm);
  else
    METIS_NodeND(&nvtxs, xadj, adjncy, &numflag, options, perm, iperm);

  for (i=0; i<nvtxs; i++) {
    ASSERT(ctrl, iperm[i]>=0 && iperm[i]<nvtxs);
    order[i] = firstnode+iperm[i];
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->SerialTmr));
}




/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: metis.h,v 1.3 2003/07/25 13:52:00 karypis Exp $
 */

/*
#define	DEBUG		1
#define	DMALLOC		1
*/

#include <metis_svfsi_stdheaders.h>

#ifdef DMALLOC
#include <dmalloc.h>
#endif

// moved this to metis_svfsi directory 
// #include "../metis_svfsi_parmetis.h"  /* Get the idxtype definition */
#include <metis_svfsi_parmetis.h>  /* Get the idxtype definition */
#include <metis_svfsi_defs.h>
#include <metis_svfsi_struct.h>
#include <metis_svfsi_macros.h>
#include <metis_svfsi_rename.h>
#include <metis_svfsi_proto.h>


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * par_metis.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: parmetislib.h,v 1.2 2003/07/21 17:50:22 karypis Exp $
 */

/*
#define DEBUG			1
#define DMALLOC			1
*/

#include <parmetis_svfsi_stdheaders.h>

// updated to place in current directory
#include <parmetis_svfsi_parmetis.h>
// #include "../parmetis_svfsi_parmetis.h"

#ifdef DMALLOC
#include <dmalloc.h>
#endif

#include <parmetis_svfsi_rename.h>
#include <parmetis_svfsi_defs.h>
#include <parmetis_svfsi_struct.h>
#include <parmetis_svfsi_macros.h>
#include <parmetis_svfsi_proto.h>


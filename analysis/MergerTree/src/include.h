//==================================================================
// will include all
//   - DEFINEFLAGS,
//   - relevant prototypes (general C and MergerTree specific),
//   - typedefs
//     -> but *not* the common variables
//==================================================================

#ifndef INCLUDE_INCLUDE_H
#define INCLUDE_INCLUDE_H

// all the DEFINEFLAGS defining the mode of operation of MergerTree
#include "define.h"

// C standard includes
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>
#include <libgen.h>
#include <time.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

// all typedefs for the structures etc.
#include "tdef.h"

// own prototypes
#include "libutil.h"
#include "libio.h"
#include "libmtree.h"
#include "libmerit.h"
#ifdef USE_PIDMAP
#include "libpidmap.h"
#endif
#ifdef SNAPSKIPPING
#include "libsnapskipping.h"
#endif

#endif

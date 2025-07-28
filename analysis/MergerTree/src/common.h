//=================================================
// all those (inelegant) common variables
// (yes, this is bad coding, but c'est la vie.)
//=================================================

#ifndef INCLUDE_COMMON_H
#define INCLUDE_COMMON_H

#include "stdint.h"
#include "inttypes.h"

#include "tdef.h"

extern HALOptr     halos[2];
extern PARTptr     parts[2];
extern uint64_t    nHalos[2];
extern uint64_t    PidMax[2], PidMax_global;
extern uint64_t    PidMin;

extern uint64_t   *PidMap[2]; // only used with "#define USE_PIDMAP"
extern uint64_t    NPids[2];
#endif

#include "../param.h"
#include "../tdef.h"

#ifndef GENERATE_GRIDS_INCLUDED
#define GENERATE_GRIDS_INCLUDED

gridls  *gen_domgrids     (int *no_grids);
boolean  gen_refgrid      (gridls **grid_list, int *no_grids);
boolean  gen_AMRhierarchy (gridls **grid_list, int *no_grids);

#endif


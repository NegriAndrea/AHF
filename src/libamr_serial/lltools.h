#include "../param.h"
#include "../tdef.h"

#ifndef LLTOOLS_INCLUDED
#define LLTOOLS_INCLUDED

void     ll                (long unsigned npart, partptr fst_part, gridls *grid);
void     NULL_ll           (gridls *cur_grid);
void     NULL_newll        (gridls *cur_grid);
void     extended_llsearch (partptr cur_part, gridls *cur_grid);
void     coagrid_llsearch  (partptr cur_part, gridls *cur_grid);
void     add_part_to_ll    (nptr new_node, partptr cur_part);

#endif


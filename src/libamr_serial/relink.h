#include "../param.h"
#include "../tdef.h"

#ifndef RELINK_INCLUDED
#define RELINK_INCLUDED

boolean  relink            (gridls *coa_grid, gridls *fin_grid);
void     relink_back       (gridls *coa_grid, gridls *fin_grid);
void     relink_bk_edge    (gridls *coa_grid, gridls *fin_grid);

#endif


#include "../param.h"
#include "../tdef.h"

#ifndef COARSE_TO_FINE_INCLUDED
#define COARSE_TO_FINE_INCLUDED

void     c2f_pot          (gridls *coa_grid, gridls *fin_grid);
void     c2f_temp2        (gridls *coa_grid, gridls *fin_grid);
void     cf_Laplace       (gridls *coa_grid, gridls *fin_grid);

#endif


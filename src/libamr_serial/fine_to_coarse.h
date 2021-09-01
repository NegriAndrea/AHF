#include "../param.h"
#include "../tdef.h"

#ifndef FINE_TO_COARSE_INCLUDED
#define FINE_TO_COARSE_INCLUDED

void     f2c_dens         (gridls *fin_grid);
void     f2c_pot          (gridls *fin_grid);
void     f2c_volume       (gridls *fin_grid);
void     fc_overlap       (gridls *fin_grid);

#endif


#include "../common.h"
#include "../param.h"
#include "../tdef.h"

#ifndef SOLVE_GRAVITY_INCLUDED
#define SOLVE_GRAVITY_INCLUDED

/* the major wrapper for all things gravity */
void    solve_gravity      (gridls *grid_list, int curgrid_no);

/* the actual routines doing the work... */
void    solve_dom_gravity  (gridls *grid_list);
void    solve_ref_gravity  (gridls *grid_list, int grid_no);
void    solve_cg           (gridls *cur_grid);
boolean converged          (gridls *cur_grid);
boolean slow_conv          (gridls *cur_grid);

#endif


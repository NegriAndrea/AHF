#include "../param.h"
#include "../tdef.h"

#ifndef DENSITY_INCLUDED
#define DENSITY_INCLUDED

void    zero_dens        (gridls *cur_grid);
void    zero_pot         (gridls *cur_grid);
void    zero_temp1       (gridls *cur_grid);
double  sum_dens         (gridls *cur_grid);
void    assign_part      (gridls *cur_grid, pqptr fst_pquad, cqptr fst_cquad, 
			      nqptr fst_nquad, nptr cur_node, 
			      long z, long y, long x, partptr cur_part);
void    unassign_part    (gridls *cur_grid, pqptr fst_pquad, cqptr fst_cquad, 
			      nqptr fst_nquad, nptr cur_node, 
			      long z, long y, long x, partptr cur_part);
boolean assign_dens      (gridls *cur_grid);
boolean assign_npart     (gridls *cur_grid);

void    overlap_dens     (gridls *cur_grid);
void    restore_dens     (gridls *cur_grid);
void    refill_dens      (gridls *coa_grid);
void    reflect_dens     (gridls *cur_grid, boolean reflect);
void    adjust_dens      (gridls *cur_grid, boolean reflect);
void    stack_dens       (gridls *first_grid, gridls *last_grid);
void    zero_temp1       (gridls *cur_grid);

#endif


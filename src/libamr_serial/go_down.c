/* the important definitions have to be included first */
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "amr_serial.h"

/*================================================================================
 * interpolate potential values from coarse to fine grid
 *
 * this ensures two things:
 *   a) we get the correct boundary values
 *   b) we start iterating with a solution not too far from the truth
 *================================================================================*/
gridls *go_down(gridls *cur_grid)
{
  gridls *coa_grid;       /* coarse grid  */
  gridls *fin_grid;       /* fine grid    */
  pqptr coa_pquad;        /* coarse pquad */
  pqptr fin_pquad;        /* fine pquad   */
  
  /* initialize fine and coarse pointers */
  coa_grid  = cur_grid;
  fin_grid  = cur_grid + 1;
  fin_pquad = fin_grid->pquad;
  coa_pquad = coa_grid->pquad;
  
#ifndef AHFlean
  c2f_pot(coa_grid, fin_grid);
#endif
  
  return(fin_grid);
}


#include "include.h"
#include "common.h"


/*==================================================================================================
 * max_merit
 *==================================================================================================*/
uint64_t max_merit(uint64_t jhalo, int isimu)
{
  uint64_t ihalo;
  
  /* mtree[] is ordered by merit and hence we only need to check the first entry */
  if(halos[isimu][jhalo].ncroco > 0) {
    return(halos[isimu][jhalo].mtree[0].id[1]);
  }
  else {
#ifdef DEBUG
    fprintf(stderr,"jhalo=%ld in isimu=%d does not point to anywhere!?\n",jhalo,isimu);
#endif
    return(0);
  }
}

/*==================================================================================================
 * merit_sort
 *==================================================================================================*/
int merit_sort(const void *mtree1, const void *mtree2)
{
  double merit1, merit2;
  
  merit1 = ((MTREEptr)mtree1)->merit;
  merit2 = ((MTREEptr)mtree2)->merit;
  
  if(merit1 > merit2)
    return(-1);
  else if(merit1 < merit2)
    return(+1);
  else
    return(0);
}


#include "include.h"
#include "common.h"

/*==================================================================================================
 * pid_cmp
 *==================================================================================================*/
int pid_cmp(const void *id1, const void *id2)
{
  uint64_t *pid1, *pid2;
  
  pid1 = (uint64_t *)id1;
  pid2 = (uint64_t *)id2;
  
  if(*pid1 < *pid2)
    return(-1);
  else if(*pid1 > *pid2)
    return(+1);
  else
    return(0);
}

/*==================================================================================================
 * create_PidMap
 *
 * PidMap[isimu] shall be used like this:
 *
 *     PidMap[isimu][ipart] = Pid;
 *     ipart                =  (uint64_t) bsearch(&(Pid), PidMap[isimu], NPids[isimu], sizeof(uint64_t), pid_cmp);
 *
 *       -> PidMap[isimu] maps the unique 'Pid' onto a unique 'ipart = [0, NPids[simu]-1]'
 *==================================================================================================*/
void create_PidMap(int isimu)
{
  uint64_t  ihalo, ipart;
  uint64_t *tmp_PidMap, tmp_NPids;
  clock_t   elapsed;
  elapsed = clock();
  
  fprintf(stderr,"  o creating PidMap for simu=%d ... ",isimu);
  
  tmp_NPids   = (uint64_t)0;
  tmp_PidMap  = NULL;
  
  for(ihalo=(uint64_t)0; ihalo<nHalos[isimu]; ihalo++) {
    for(ipart=(uint64_t)0; ipart<halos[isimu][ihalo].npart; ipart++) {
      tmp_PidMap            = (uint64_t *)realloc(tmp_PidMap, (tmp_NPids+1)*sizeof(uint64_t));
      tmp_PidMap[tmp_NPids] = halos[isimu][ihalo].Pid[ipart];
      tmp_NPids++;
    }
  }
  
  // make PidMap[] searchable via bsearch()
  qsort((void *)tmp_PidMap, tmp_NPids, sizeof(uint64_t), pid_cmp);
  
  // remove duplicates from sorted list
  PidMap[isimu] = NULL;
  NPids[isimu]  = (uint64_t)0;
  for(ipart=0; ipart<tmp_NPids-1; ipart++) {
    while(tmp_PidMap[ipart] == tmp_PidMap[ipart+1]) {
      ipart++;
    }
    PidMap[isimu]               = (uint64_t *)realloc(PidMap[isimu] , (NPids[isimu]+1)*sizeof(uint64_t));
    PidMap[isimu][NPids[isimu]] = tmp_PidMap[ipart];
    NPids[isimu]++;
  }
  PidMap[isimu]               = (uint64_t *)realloc(PidMap[isimu] , (NPids[isimu]+1)*sizeof(uint64_t));
  PidMap[isimu][NPids[isimu]] = tmp_PidMap[ipart];
  NPids[isimu]++;
  free(tmp_PidMap);
  
  PidMax[isimu] = NPids[isimu];
  
  elapsed = clock()-elapsed;
  
  fprintf(stderr," done in %4.2f sec. (removed %"PRIu64" duplicates, keeping %"PRIu64" unique Pids)\n", (float)elapsed/CLOCKS_PER_SEC,(tmp_NPids-NPids[isimu]),NPids[isimu]);
}


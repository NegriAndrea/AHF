#ifdef DARK_ENERGY
#include <math.h>
#include <stdio.h>
#include <string.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "darkenergy.h"
#include "utility.h"

struct dark_energy DarkEnergy;

void read_dark_energy_table(char *dark_energy_url)
{
  FILE   *fdark_energy=NULL;
  char    dummyline[MAXSTRING];
  long   *idx, i;
  double *a, *Hubble, *OmegaM;
  
  fdark_energy = fopen(dark_energy_url,"r");
  
  if(fdark_energy==NULL)
    fprintf(stderr, "\nError. File %s not found, required by read_dark_energy_table().\n", dark_energy_url);
  
  DarkEnergy.n_entries = -1;
  DarkEnergy.a = calloc(1,sizeof(double));
  DarkEnergy.Hubble = calloc(1,sizeof(double));
  DarkEnergy.OmegaM = calloc(1,sizeof(double));
  
  while(fgets(dummyline,MAXSTRING,fdark_energy) != NULL)
   {
    DarkEnergy.n_entries++;
    DarkEnergy.a = realloc(DarkEnergy.a, (DarkEnergy.n_entries+1)*sizeof(double));
    DarkEnergy.Hubble = realloc(DarkEnergy.Hubble, (DarkEnergy.n_entries+1)*sizeof(double));
    DarkEnergy.OmegaM = realloc(DarkEnergy.OmegaM, (DarkEnergy.n_entries+1)*sizeof(double));
    
    sscanf(dummyline,"%lf %lf %lf", &(DarkEnergy.a[DarkEnergy.n_entries]),
           &(DarkEnergy.Hubble[DarkEnergy.n_entries]),
           &(DarkEnergy.OmegaM[DarkEnergy.n_entries]));
   }
  
  // sort DarkEnergy. entries by a[] (interpolate() expects this ascending ordering!)
  idx = (long *) calloc(DarkEnergy.n_entries, sizeof(long));
  indexx(DarkEnergy.n_entries, (DarkEnergy.a)-1, idx-1);

  a      = (double *) calloc(DarkEnergy.n_entries, sizeof(double));
  Hubble = (double *) calloc(DarkEnergy.n_entries, sizeof(double));
  OmegaM = (double *) calloc(DarkEnergy.n_entries, sizeof(double));
  for(i=0; i<DarkEnergy.n_entries; i++) {
    a[i]      = DarkEnergy.a[idx[i]-1];
    Hubble[i] = DarkEnergy.Hubble[idx[i]-1];
    OmegaM[i] = DarkEnergy.OmegaM[idx[i]-1];
  }
  for(i=0; i<DarkEnergy.n_entries; i++) {
    DarkEnergy.a[i]      = a[i];
    DarkEnergy.Hubble[i] = Hubble[i];
    DarkEnergy.OmegaM[i] = OmegaM[i];
    
    fprintf(stderr,"%lf %lf %f\n",DarkEnergy.a[i],DarkEnergy.Hubble[i],DarkEnergy.OmegaM[i]);
  }
  
  free(idx);
  free(a);
  free(Hubble);
  free(OmegaM);
  
  
  fclose(fdark_energy);
}

double Hubble_DE(double a)
{
	return interpolate(DarkEnergy.a, DarkEnergy.Hubble, DarkEnergy.n_entries, a);
}

double OmegaM_DE(double a)
{
  	return interpolate(DarkEnergy.a, DarkEnergy.OmegaM, DarkEnergy.n_entries, a);
}

#endif

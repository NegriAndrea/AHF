#include <math.h>
#include <stdio.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"


/*
 * we need to cut-and-paste the create_timeline() routine from ../src/cosmology.c
 * in order to be able to bypass the #MAXTIME definition...
 */
void create_timeline2(double a_init, double a_fin, tlptr timeline);
#define MAXTIME2 200


typedef struct {
  int    isnap;
  double zred;
} zred_t;

zred_t *read_zredlist(char *, int *);

int main(argc,argv)
int argc;
char **argv;
{
  int     no_outputs, i, nzred;
  double  z_init, z_final, a_init, a_final, a, t, t0, hubble;
  char    zredlist[MAXSTRING];
  zred_t *data;
  
  if(argc != 5) {
    fprintf(stderr,"usage: agelist omega0 lambda0 hubble0 zredlist\n");
    exit(0);
  }
  
  simu.omega0  = atof(argv[1]);
  simu.lambda0 = atof(argv[2]);
  hubble       = atof(argv[3]);
  strcpy(zredlist,argv[4]);
  
  data = read_zredlist(zredlist, &nzred);
  
  z_init  = data[nzred-1].zred*2;
  z_final = 0.;
  a_init  = 1.0/(1.0+z_init);
  a_final = 1.0/(1.0+z_final);
  create_timeline2(a_init, a_final, &simu.timeline);

  t0 = calc_t(1./(1.+data[0].zred))*Mpc/H0/Gyr/hubble;
  
  fprintf(stderr,"#    snapnum       a             z           t(t0)      t(year)\n");
  for(i=0; i<nzred; i++) {
    a = 1./(1.+data[i].zred);
    t = calc_t(a)*Mpc/H0/Gyr/hubble;
    fprintf(stderr,"%d %lf %lf %lf %g\n",
            data[i].isnap,
            1./(1.+data[i].zred),
            data[i].zred,
            t/t0,
            t*1e9);
  }
  
  
  
  
  free(data);
}

/*===========================================================================
 * read_zredlist()
 *===========================================================================*/
zred_t *read_zredlist(char *zredlist, int *nzred)
{
  zred_t *data;
  double  z;
  int     isnap;
  FILE   *fp;
  char    line[MAXSTRING];
  
  
  fp = fopen(zredlist,"r");
  fgets(line,MAXSTRING,fp);
  *nzred = 0;
  data   = NULL;
  while(!feof(fp)) {
    sscanf(line,"%d %lf",&isnap,&z);
    
    (*nzred)++;
    data = (zred_t *) realloc(data, (*nzred)*sizeof(zred_t));
    data[*nzred-1].isnap = isnap;
    data[*nzred-1].zred  = z;

    fgets(line,MAXSTRING,fp);
  }
  fclose(fp);
  
  return(data);
}

/*===========================================================================
 * create a timeline for current cosmological model
 *===========================================================================*/
void create_timeline2(double a_init, double a_fin, tlptr timeline)
{
  int    iloop;
  double a,t,omega,lambda,rhoc,hubble;
  
  for(iloop=0; iloop<MAXTIME2; iloop++)
   {
    a      = ((double)iloop+1.0)/(double)MAXTIME2 * (a_fin-a_init) + a_init;
    t      = calc_t(a);
    omega  = calc_omega(a);
    lambda = calc_lambda(a);
    hubble = H0    * sqrt(
                          simu.lambda0 * (1.-1./pow2(a))         +
                          simu.omega0  * (1./pow3(a)-1./pow2(a)) + 1./pow2(a));
    rhoc   = rhoc0 *     (
                          simu.lambda0*(1.-1./pow2(a))           +
                          simu.omega0*(1./pow3(a)-1./pow2(a))    + 1./pow2(a));
    
    timeline->a[iloop]      = a;
    timeline->t[iloop]      = t;
    timeline->omega[iloop]  = omega;
    timeline->lambda[iloop] = lambda;
    timeline->hubble[iloop] = hubble;
    timeline->age[iloop]    = t*Mpc/H0/Gyr;
    timeline->virial[iloop] = calc_virial(a);
   }
}
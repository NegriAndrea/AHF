#include <math.h>
#include <stdio.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

int main()
{
  int     no_outputs, i, UseRhoBack;
  double  zred, z_init, z_final, a_init, a_final, a;
  
  printf("===============================\n");
  printf(" calculate virial overdensity\n");
  printf("===============================\n");
  printf("please give omega0:                     ");
  scanf("%lf",&simu.omega0);
  printf("please give lambda0:                    ");
  scanf("%lf",&simu.lambda0);
  printf("please give redshift z:                 ");
  scanf("%lf",&zred);
  printf("in units of RhoBack(1) or RhoCrit(0):   ");
  scanf("%d",&UseRhoBack);
  printf("\n");
  
  simu.UserDvir   = -1;
  simu.UseRhoBack = UseRhoBack;
  
  z_init  = 100.;
  z_final = 0.;
  a_init  = 1.0/(1.0+z_init);
  a_final = 1.0/(1.0+z_final);
  
  create_timeline(a_init, a_final, &simu.timeline);
  
  
  a = 1./(1.+zred);
  printf("virial overdensity = %lf\n",calc_virial(a));
}

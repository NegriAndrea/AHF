#include <math.h>
#include <stdio.h>
#include "../param.h"
#include "../tdef.h"

#ifndef COSMOLOGY_INCLUDED
#define COSMOLOGY_INCLUDED

double calc_rho_crit  (double a);
double calc_rho_vir   (double a);
double calc_omega     (double a);
double calc_lambda    (double a);
double calc_Hubble    (double a);
double calc_Hubble_VDE(double a);
double calc_a         (double t);
double calc_t         (double a);
double calc_super_t   (double a);
double calc_super_a   (double super_t);
double calc_a_dot     (double a);
double calc_a_ddot    (double a);
double calc_q         (double a);
double calc_virial    (double a);
double calc_growth    (double a);
double collapse       (double age, double acc);

void   create_timeline(double a_init, double a_fin, tlptr timeline);

#endif


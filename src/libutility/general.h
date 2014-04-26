#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef GENERAL_INCLUDED
#define GENERAL_INCLUDED

long    modulo             (long a, long  b);
boolean set_opposite       (boolean tf);
boolean is_even            (long a);
void    freen              (void **ptr);
void    itoa_              (int n, char s[]);
void    find_max           (double *x, double *y, int num, int smooth, double *maxx, double *maxy);
void    find_min           (double *x, double *y, int num, int smooth, double *minx, double *miny);
void    find_max_spline    (double *x, double *y, int num, int smooth, double *maxx, double *maxy, int steps);
void    find_min_spline    (double *x, double *y, int num, int smooth, double *minx, double *miny, int steps);
int     descending         (void *p1, void *p2);
int     ascending          (void *p1, void *p2);
void    sexchange          (void *,size_t);
int     test_endian        ();
int     ReadString         (FILE *fptr, char         *s, int swap);
int     ReadChars          (FILE *fptr,char *s,int n);
int     ReadLong           (FILE *fptr, long         *s, int swap);
int     ReadULong          (FILE *fptr,unsigned long *n,int swap);
int     ReadLongLong       (FILE *fptr, long long    *s, int swap);
int     ReadDouble         (FILE *fptr, double       *s, int swap);
int     ReadFloat          (FILE *fptr, float        *s, int swap);
int     ReadInt            (FILE *fptr, int          *s, int swap);
int     ReadUInt           (FILE *fptr, unsigned int *s, int swap);
void    smooth3            (double *data, int nbins, int nsmooth);
void    smooth5            (double *data, int nbins, int nsmooth);
int     sign               (double a);
double  cubic              (double a1, double a2, double a3, double a4);
double  INTEGRATE          (double (*func)(double), double a, double b, double step, double eps);
void    RUNGE              (double y, double dydx, double x, double h, double(*func)(double),
                            double *yout, double *yerr);
void    RUNGE5VAR          (double *y, double dydx, double *x, double htry, double eps, 
                            double yscale, double *hnext, double (*func)(double));
int     zbrac              (double (*func)(double, double, double), double *x1, double *x2,
                            double param1, double param2);  // additinal variables required by mond.c
double  rtbis              (double (*func)(double,double,double), double x1, double x2, double xacc, 
                            double param1, double param2);  // additinal variables required by mond.c
void    spline             (double x[], double y[], int n, double yp1, double ypn, double y2[]);
void    splint             (double xa[], double ya[], double y2a[], int n, double x, double *y);
void    polint             (double xa[], double ya[], int n, double x, double *y, double *dy);
void    jacobi             (double a[4][4], int n, double d[], double v[4][4], int *nrot);
float   ran3               (int *idum);
void    indexx             (unsigned long n, double arr[], unsigned long indx[]);
int    *ivector            (long nl, long nh);
void    free_ivector       (int *v, long nl, long nh);
double  interpolate        (double*, double*, int, double);

void    fourn              (flouble data[], unsigned long nn[], int ndim, int isign);

#endif

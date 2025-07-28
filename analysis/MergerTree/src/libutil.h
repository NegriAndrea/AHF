
#ifndef INCLUDE_LIBUTIL_H
#define INCLUDE_LIBUTIL_H

#include "include.h"

/*=============================================================================
 * some convenient shortcuts/abbreviations
 *=============================================================================*/
#define MAXSTRING    2048   /* used for char statement, i.e. filenames etc.     */
#define PI           3.14159265358979323846264338
#define TRUE         1
#define FALSE        0
#define ZERO         (1e-6)
#define MIN(A,B)        ((A)<(B)?(A):(B))
#define MAX(A,B)        ((A)>(B)?(A):(B))
#define FUNC(x)         ((*func)(x))
#define SWAP(a,b,temp)  {temp=(a);(a)=(b);(b)=temp;}
#define Re(ix,iy,iz,L)  (2*((iz)*(L)*(L) + (iy)*(L) + (ix))    )
#define Im(ix,iy,iz,L)  (2*((iz)*(L)*(L) + (iy)*(L) + (ix)) + 1)
#define FRAC(x)         (((x)>1.)?(1./(x)):(x))
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))

/*=============================================================================
 * some physical constants...
 *=============================================================================*/
#define Gyr       3.1558e16         /* [sec]                   */
#define Mpc       3.08567782e19     /* [km]                    */
#define H0        100.              /* [h*km]/[sec*Mpc]        */
#define rhoc0     2.7755397e11      /* [h^2*Msun]/[Mpc^3]      */
#define Grav      4.3006485e-9      /* [Mpc*km^2]/[Msun*sec^2] */
#define cH0	      2998.0		        /* c/H0 (in h^-1 Mpc)      */
#define kB_per_mp 0.825481286614E-2 /* [(km/sec)^2/K]          */
#define kBoltzman 6.9416792         /* [(km/sec)^2 Msun/K]     */
#define Msun      1.9891e30         /* [kg]                    */

/*=============================================================================
 * PROTOYPES OF MergerTree specific ROUTINES
 *=============================================================================*/
char    *xdirname               (const char *path);
char    *xbasename              (const char *path);
void     construct_filename     (char filename[MAXSTRING], int32_t, char infile[MAXSTRING]);
uint64_t count_halos            (char filename[MAXSTRING], int32_t);
int32_t  count_particles_files  (char filename[MAXSTRING]);
void     dump_defines           ();
int      check_Ptype            (uint64_t Ptype);

/*=============================================================================
 * PROTOTYPES OF NUMERICAL RECIPES ET AL. ROUTINES
 *=============================================================================*/
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

void    fourn              (double data[], unsigned long nn[], int ndim, int isign);

#endif

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utility.h"

/* set_opposite return opposite boolean value */
boolean set_opposite(boolean tf)
{
   if(tf == TRUE)
      return FALSE;
   else
      return TRUE;
}

/* is_even: test if int argument is even */
boolean is_even(long a)
{
   int b = 2;
   
   if((fmod(a,b)) > 0.0)
      return FALSE;
   else
      return TRUE;
}

/* freen: free memory nullify pointer */
void freen(void **ptr)
{
   free(*ptr);
   *ptr = NULL;
}

/* reverse: reverses  string s in place */
void reverse(char s[])
{
   int i, j;
   char c;
   
   for(i = 0, j = strlen(s) - 1; i < j; i++, j--){
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

/* itoa: convert n to charcters in s[] */
void itoa_(int n, char s[])
{
   int i;
   
   i = 0;
   do
     {     
        s[i++] = n % 10 + '0';
     }
   while ((n /= 10) > 0);
   
   s[i] = '\0';
   reverse(s);
}

/* modulo: a modulo b for integer numbers */
long modulo(long a, long  b)
{
   if(a < 0)
      return a + b;
   if(a >= b)
      return a - b;
   else
      return a;
}

/* emulates FORTRAN pause statement */
void pause(void)
{
   char input[100]; 
   
   printf("\nPAUSE - type 'g' to continue\n");
   scanf("%s",input);
   
   if(input[0] != 'g') exit(0);
   
   printf("\n");
}

/* these two functions descending and ascending are for qsort() */
int descending(void *p1, void *p2)
{
   
   long unsigned g1, g2;
   
   g1 = *(long unsigned *) p1;
   g2 = *(long unsigned *) p2;
   
   if (g1 > g2)
      return (-1);
   else if (g1 < g2)
      return (1);
   else
      return (0);
}

int ascending(void *p1, void *p2)
{
   
   long unsigned g1, g2;
   
   g1 = *(long unsigned *) p1;
   g2 = *(long unsigned *) p2;
   
   if (g1 > g2)
      return (-1);
   else if (g1 < g2)
      return (1);
   else
      return (0);
}

void sexchange(void *p,size_t s)
{
  int n;
  unsigned char ptmp,*pc;
  
  pc = (unsigned char *)p;
  
  for(n=0; n < s/2; n++) {
    ptmp = pc[n];
    pc[n] = pc[s - n - 1];
    pc[s - n - 1] = ptmp;
  }
}

int test_endian()
{
  /* Are we little or big endian?  From Harbison & Steele.  */
  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  if (u.c[sizeof (long) - 1] == 1){
    return(1); /* Big Endian */
  }
  else if (u.c[0] == 1) {
    return(0); /* Little Endian */
  }
  else {
    printf("We should not be here.  VERY bad things are happening!\n");
  }
  return (-1);
}

double interpolate(double *x_array, double *y_array, int n_array, double x)
{
  /* interpolate using IORD points around x_array[j] and y_array[j] */
  int     jl,ju,jm,j,IORD=6;
  double  h, dh;
  
  jl = 0;
  ju = n_array+1;
  while(ju-jl > 1)
   {
    jm = (ju+jl)/2;
    if(x < x_array[jm])
      jl=jm;
    else
      ju=jm;
   }
  j=jl;
  if(j >= n_array-IORD)
    j = n_array-IORD-1;
  polint(&x_array[j],&y_array[j],IORD,x,&h,&dh);
  
  return(h);
}

/*==========================================================================
 * cubic:
 c
 c     Function for evaluating a cubic equation.  In the case where
 c     there is one real root, this routine returns it.  If there are three
 c     roots, it returns the smallest positive one.
 c
 c     The cubic equation is a1*x**3 + a2*x**2 + a3*x + a4 = 0.
 c
 c     Ref: Tallarida, R. J., "Pocket Book of Integrals and Mathematical
 c     Formulas," 2nd ed., Boca Raton: CRC Press 1992, pp. 8-9.
 c
 *==========================================================================*/
double cubic(double a1, double a2, double a3, double a4)
{
  double a,b,d,r,theta,pi,ar,ai,y,y1,y2,y3,p,q;
  
  if (fabs(a1) < ZERO)
    {
      fprintf(io.logfile,"ERROR: Quadratic/linear equation passed to CUBIC.");
      return(0.0);
    }
  
  
  p = a3/a1 - pow2(a2/a1)/3.;
  q = a4/a1 - a2*a3/pow2(a1)/3. + 2.*pow3(a2/a1)/27.;
  
  d = pow3(p)/27. + pow2(q)/4.;
  
  if (d > 0.)
    {
      a = -0.5*q + sqrt(d);
      b = -0.5*q - sqrt(d);
      
      if (a > 0.)
        a =   pow( a, 1./3.);
      else
        a = - pow(-a, 1./3.);
      
      if (b > 0.)
        b =   pow( b, 1./3.);
      else
        b = - pow(-b, 1./3.);
      
      y = a + b;
    }
  else
    {
      ar    = -0.5*q;
      ai    = sqrt(-d);
      r     = pow( pow2(ar)+pow2(ai), 1./6.);
      theta = atan2(ai,ar);
      y1    = 2. * r * cos(theta/3.) - a2/a1/3.;
      y     = y1;
      pi    = 4.*atan(1.);
      y2    = 2. * r * cos(theta/3.+2.*pi/3.) - a2/a1/3.;
      
      if (y < 0. || (y2 > 0. && y2 < y))
        y = y2;
      
      y3 = 2. * r * cos(theta/3.-2.*pi/3.) - a2/a1/3.;
      
      if (y < 0. || (y3 > 0. && y3 < y))
        y = y3;
    }
  
  return(y);
}

/*
 Read a string of n characters
 */
int ReadString(FILE *fptr,char *s,int n)
{
  int i,c;
  
  if(sizeof(char) != 1)
    {
      fprintf(stderr,"ReadString: sizeof(char)=%ld and not 1\n",sizeof(char));
      exit(0);
    }
  
  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(FALSE);
    s[i] = c;
    s[i+1] = '\0';
  }
  return(TRUE);
}

/*
 Read an array of n characters
 NOTE: the difference to ReadString() is that we do not '\0'-terminate the array
 */
int ReadChars(FILE *fptr,char *s,int n)
{
   int i,c;
   
   if(sizeof(char) != 1)
     {
      fprintf(stderr,"ReadChars: sizeof(char)=%ld and not 1\n",sizeof(char));
      exit(0);
     }
   
   s[0] = '\0';
   for (i=0;i<n;i++) {
      c = fgetc(fptr);
      if (c == EOF)
         return(FALSE);
      s[i] = c;
   }
   return(TRUE);
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
    {
      fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
      exit(0);
    }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
    {
      fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
      exit(0);
    }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped long integer
 */
int ReadLong(FILE *fptr,long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(long) == 4)
    {
      if (fread(n,4,1,fptr) != 1)
        return(FALSE);
      if (swap) {
        cptr = (unsigned char *)n;
        tmp     = cptr[0];
        cptr[0] = cptr[3];
        cptr[3] = tmp;
        tmp     = cptr[1];
        cptr[1] = cptr[2];
        cptr[2] = tmp;
      }
    }
  else if(sizeof(long) == 8)
    {
      if (fread(n,8,1,fptr) != 1)
        return(FALSE);
      if (swap) {
        cptr = (unsigned char *)n;
        tmp     = cptr[0];
        cptr[0] = cptr[7];
        cptr[7] = tmp;
        tmp     = cptr[1];
        cptr[1] = cptr[6];
        cptr[6] = tmp;
        tmp     = cptr[2];
        cptr[2] = cptr[5];
        cptr[5] = tmp;
        tmp     = cptr[3];
        cptr[3] = cptr[4];
        cptr[4] = tmp;
      }
    }
  else
    {
      fprintf(stderr,"ReadLong: something wrong...cannot read long\n");
      exit(0);
    }
  
  
  
  return(TRUE);
}

/*
 Read a possibly byte swapped unsigned long integer
 */
int ReadULong(FILE *fptr,unsigned long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(unsigned long) == 4)
   {
    if (fread(n,4,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
   }
  else if(sizeof(unsigned long) == 8)
   {
    if (fread(n,8,1,fptr) != 1)
      return(FALSE);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
    }
   }
  else
   {
    fprintf(stderr,"ReadULong: something wrong...cannot read long\n");
    exit(0);
   }
  
  return(TRUE);
}

/*
 Read a possibly byte swapped long long integer
 */
int ReadLongLong(FILE *fptr,long long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if (fread(n,8,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  return(TRUE);
}

/*
 Read a possibly byte swapped double precision number
 Assume IEEE
 */
int ReadDouble(FILE *fptr,double *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(double) != 8)
    {
      fprintf(stderr,"ReadDouble: sizeof(double)=%ld and not 8\n",sizeof(double));
      exit(0);
    }
  
  if (fread(n,8,1,fptr) != 1)
    return(FALSE);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  
  return(TRUE);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(float) != 4)
    {
      fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
      exit(0);
    }
  
  if (fread(n,4,1,fptr) != 1)
    return(FALSE);
  if (swap) 
    {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
  return(TRUE);
}

/*==============================================================================
 * smooth data array
 *==============================================================================*/
void smooth3(double *data, int nbins, int nsmooth)
{
  int    i,j;
  double *tmpdata;
  
  if (nsmooth == 0 || nbins < 3)
    return;
  
  tmpdata = (double *)calloc(nbins, sizeof(double));
  
  for(j=0; j<nsmooth; j++)
    { 
      tmpdata[0] = (data[0]+data[1])/2.;
      
      for(i=1; i<nbins-1; i++)
        tmpdata[i] = (data[i-1]+data[i]+data[i+1])/3.;
      
      tmpdata[nbins-1] = (data[nbins-1]+data[nbins-2])/2.;
      
      for(i=0; i<nbins; i++)
        data[i] = tmpdata[i];
    }
  
  free(tmpdata);
  
}

/*==============================================================================
 * smooth data array
 *==============================================================================*/
void smooth5(double *data, int nbins, int nsmooth)
{
  int    i,j;
  double *tmpdata;
  
  tmpdata = (double *)calloc(nbins, sizeof(double));
  
  for(j=0; j<nsmooth; j++)
    { 
      tmpdata[0] = (data[0]+data[1]+data[2])/3.;
      tmpdata[1] = (data[0]+data[1]+data[2]+data[3])/4.;
      
      for(i=2; i<nbins-2; i++)
        tmpdata[i] = (data[i-2]+data[i-1]+data[i]+data[i+1]+data[i+2])/5.;
      
      tmpdata[nbins-2] = ( data[nbins-1]+data[nbins-2]
                          +data[nbins-3]+data[nbins-4])/4.;
      tmpdata[nbins-1] = (data[nbins-1]+data[nbins-2]+data[nbins-3])/3.;
      
      for(i=0; i<nbins-1; i++)
        data[i] = tmpdata[i];
    }
  
  free(tmpdata);
  
}

/*==============================================================================
 *  find maximum within array y[]
 *==============================================================================*/
void find_max(double *x, double *y_raw, int nbins, int ismooth, double *x_max, double *y_max)
{
  int    ibin, ibin_max[2];
  double ymax;
  double *y;
   
  /* do not smooth the actual array as this will obscure the peak value! */
  y = calloc(nbins,sizeof(double));
  for(ibin=0; ibin<nbins; ibin++)
     y[ibin]=y_raw[ibin];
  
  smooth3(y, nbins, ismooth);
  
  ymax = -10.0;
  
  /* approach from left */
  ibin_max[0] = nbins-1;
  for(ibin=0; ibin<nbins-1; ibin++)
    {
      if(y[ibin] > ymax)
        {
          ymax        = y[ibin];
          ibin_max[0] = ibin;
        }
    }
  
  /* approch from right */
  ibin_max[1] = ibin_max[0];
  for(ibin=(nbins-1+ibin_max[0])/2; ibin>ibin_max[0]; ibin--)
    {
      if(y[ibin] > ymax)
        {
          ymax        = y[ibin];
          ibin_max[1] = ibin;
        }
    }
  
  *x_max = (     x[ibin_max[0]] +     x[ibin_max[1]] )/2;
  *y_max = ( y_raw[ibin_max[0]] + y_raw[ibin_max[1]] )/2;
   
  /* be nice */
  free(y);
}


/**
 * \brief  Finds the maximum in a given array.
 *
 * \param  *x      An array containing the x-values.
 * \param  *y      An array containing the y-values.
 * \param  num     The length of the both arrays.
 * \param  steps   Total number of interpolation steps.
 * \param  smooth  Total number of smoothing steps.
 * \param  *maxx   A pointer to a variable that will hold the position of
 *                 the maximum.
 * \param  *maxy   A pointer to a variable that will hold the value of
 *                 the function at the found maximum.
 *
 * \return  Returns nothing.
 */
void find_max_spline(double *x, double *y, int num, int smooth, double *maxx, double *maxy, int steps)
{
  double *y2, xi, yi, incr;
  
  /* Smooth the data */
  smooth3(y, num, smooth);
  
  /* obtain second derivatives to be used with splint(), i.e. "spline interpolation" */
  y2 = (double *) calloc(num, sizeof(double));
  spline(x-1, y-1, num, 2E33, 2E33, y2-1);
  
  /* scan the spline interpolation using splint() */
  *maxy = -1e10;
  incr = (x[num-1] - x[0])/(double)(steps-1);
  xi = x[0];
  while (xi <= x[num-1]) {
    splint(x-1, y-1, y2-1, num, xi, &yi);
    if (yi > *maxy) {
      *maxx = xi;
      *maxy = yi;
    }
    xi += incr;
  }
  
  free(y2);
  
  return;
}


/*==============================================================================
 *  find minimum within array y[]
 *==============================================================================*/
void find_min(double *x, double *y, int nbins, int ismooth, double *x_min, double *y_min)
{
  int    ibin, ibin_min[2];
  double ymax;
  
  smooth3(y, nbins, ismooth);
  
  ymax = -10.0;
  
  /* approach from left */
  ibin_min[0] = nbins-1;
  for(ibin=0; ibin<nbins-1; ibin++)
    {
      if(y[ibin] < ymax)
        {
          ymax        = y[ibin];
          ibin_min[0] = ibin;
        }
    }
  
  /* approch from right */
  ibin_min[1] = ibin_min[0];
  for(ibin=(nbins-1+ibin_min[0])/2; ibin>ibin_min[0]; ibin--)
    {
      if(y[ibin] < ymax)
        {
          ymax        = y[ibin];
          ibin_min[1] = ibin;
        }
    }
  
  *x_min = ( x[ibin_min[0]] + x[ibin_min[1]] )/2;
  *y_min = ( y[ibin_min[0]] + y[ibin_min[1]] )/2;   
}


/**
 * \brief  Finds the minimum in a given array.
 *
 * \param  *x      An array containing the x-values.
 * \param  *y      An array containing the y-values.
 * \param  num     The length of the both arrays.
 * \param  steps   Total number of interpolation steps.
 * \param  smooth  Total number of smoothing steps.
 * \param  *maxx   A pointer to a variable that will hold the position of
 *                 the maximum.
 * \param  *maxy   A pointer to a variable that will hold the value of
 *                 the function at the found maximum.
 *
 * \return  Returns nothing.
 */
void find_min_spline(double *x, double *y, int num, int smooth, double *minx, double *miny, int steps)
{
  double *y2, xi, yi, incr;
  
  /* Smooth the data */
  smooth3(y, num, smooth);
  
  /* obtain second derivatives to be used with splint(), i.e. "spline interpolation" */
  y2 = (double *) calloc(num, sizeof(double));
  spline(x-1, y-1, num, 2E33, 2E33, y2-1);
  
  /* scan the spline interpolation using splint() */
  *miny = 1e10;
  incr = (x[num-1] - x[0])/(double)(steps-1);
  xi = x[0];
  while (xi <= x[num-1]) {
    splint(x-1, y-1, y2-1, num, xi, &yi);
    if (yi < *miny) {
      *minx = xi;
      *miny = yi;
    }
    xi += incr;
  }
  
  free(y2);
  
  return;
}


/*===================================================================================
 * sign:   is there a C function for this?
 *===================================================================================*/
int sign(double a)
{
  if(a>0)
    return(+1);
  else
    return(-1);
}


/*===============================================================================
 * INTEGRATE:  integrates a supplied function func(x) from a to b accurately !
 *
 *             uses RUNGE5VAR() and RUNGE()
 *===============================================================================*/
#define safety  (double) 0.9e0
#define pgrow   (double) -0.2e0
#define pshrink (double) -0.25e0
#define errcon  (double)  1.89e-4

void RUNGE5VAR(double *y, double dydx, double *x, double htry, double eps, 
               double yscale, double *hnext, double (*func)(double))
{
  double errmax,h,hold,htemp,xnew,yerr,ytemp,temp;
  
  h      = htry;
  errmax = (double) 10.0;
  while(errmax > (double) 1.0)
    {
      RUNGE(*y,dydx,*x,h,func,&ytemp,&yerr);
      
      errmax = fabs(yerr/yscale)/eps;
      if(errmax > (double) 1.0)
        { 
          htemp = safety*h * pow(errmax,pshrink);
          hold  = h;
          temp  = MAX(fabs(htemp),(double)0.1*fabs(h));
          
          if(h  >= 0)
            h = fabs(temp);
          else
            h = -fabs(temp);
          
          xnew = *x + h;
          
          if (fabs(xnew-*x) < (double)1e-5)
            {
              /*    	      printf("WARNING: Stepsize underflow in RUNGE5VAR()\n"); */
              h      = hold;
              errmax = (double)0.0;
            }
        }
    }
  if (errmax > errcon)
    {
      *hnext = safety*h * pow(errmax,pgrow);
    }
  else
    {
      *hnext = (double) 5.0 * h;
    }
  
  *x = *x + h;
  *y = ytemp;
  
  return;
}

#define a2      ((double) 0.2e0)
#define a3      ((double) 0.3e0)
#define a4      ((double) 0.6e0)
#define a5      ((double) 1.e0)
#define a6      ((double) 0.875e0)
#define c1      ((double) 37.e0/378.e0)
#define c3      ((double) 250.e0/621.e0)
#define c4      ((double) 125.e0/594.e0)
#define c6      ((double) 512.e0/1771.e0)
#define dc1     (c1 -  (double) 2825.e0/27648.e0)
#define dc3     (c3 -  (double) 18575.e0/48384.e0)
#define dc4     (c4 -  (double) 13525.e0/55296.e0)
#define dc5     ((double) -277.e0/14336.e0)
#define dc6     (c6 -  (double) 0.25e0)

void RUNGE(double y, double dydx, double x, double h, double(*func)(double),
           double *yout, double *yerr)
{
  double ak3, ak4, ak5, ak6;
  
  ak3  = FUNC((double)(x+a3*h));
  ak4  = FUNC((double)(x+a4*h));
  ak5  = FUNC((double)(x+a5*h));
  ak6  = FUNC((double)(x+a6*h));
  
  *yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6);
  *yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);
}

#define MAXSTEPS_INTEGRATE  (1000000000)  /* needed for INTEGRATE() function         */

double INTEGRATE (double (*func)(double), double a, double b, double step,
                  double eps)
{
   double x, dx, y, dydx, yscale, dxnext;
   int    Nstep;
   
   
   x     = a;
   dx    = step;
   y     = (double)0.0;
   Nstep = 0;
   
   while( ((x-b)*(b-a) < (double)0.0) && (Nstep < MAXSTEPS_INTEGRATE))
     {
      Nstep  = Nstep + 1;
      dydx   = FUNC(x);
      yscale = MAX(fabs(y) + fabs(dx*dydx), (double)1.e-8);
      
      if ((x+dx-b)*(x+dx-a) > (double)0.0)
        {
         dx = b - x;
        }
      
      RUNGE5VAR(&y,dydx,&x,dx,eps,yscale,&dxnext,func);
      
      dx = dxnext;
     }
   
   return(y);
}


/*===========================================================================
 * NUMERICAL RECIPES
 *
 * most of the following routines have been borrowed from Numerica Recipes;
 * however, they may have one or the other adaption to better suit the
 * needs of AMIGA...
 *===========================================================================*/
#define NR_END 1
#define FREE_ARG char*
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}
double *vector(long nl, long nh)
/* allocate a vector with subscript range v[nl..nh] */
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in vector()");
   return v-nl+NR_END;
}
void free_vector(double *v, long nl, long nh)
/* free a vector allocated with vector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}

#define FACTOR 1.6
#define NTRY   50

int zbrac(double (*func)(double, double, double), double *x1, double *x2,
          double param1, double param2)
{
  int j;
  double f1,f2;
  
  f1=(*func)(*x1,param1,param2);
  f2=(*func)(*x2,param1,param2);
  for (j=1;j<=NTRY;j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2))
      f1=(*func)(*x1 += FACTOR*(*x1-*x2),param1,param2);
    else
      f2=(*func)(*x2 += FACTOR*(*x2-*x1),param1,param2);
  }
  return 0;
}
#undef FACTOR
#undef NTRY

/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */
#define JMAX 40
double rtbis(double (*func)(double,double,double), 
             double x1, double x2, double xacc, double param1, double param2)
{
  int j;
  double dx,f,fmid,xmid,rtb;
  
  f=(*func)(x1,param1,param2);
  fmid=(*func)(x2,param1,param2);
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5),param1,param2);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  fprintf(stderr,"Too many bisections in rtbis\n");
  return 0.0;
}
#undef JMAX

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c,*d;
   
   dif=fabs(x-xa[1]);
   c=vector(1,n);
   d=vector(1,n);
   for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
         ns=i;
         dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
         ho=xa[i]-x;
         hp=xa[i+m]-x;
         w=c[i+1]-d[i];
         if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
         den=w/den;
         d[i]=hp*den;
         c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   free_vector(d,1,n);
   free_vector(c,1,n);
}
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

/* ran3: random generator */
float ran3(int *idum)
{
   static int inext,inextp;
   static long ma[56];
   static int iff=0;
   long mj,mk;
   int i,ii,k;
   
   if (*idum < 0 || iff == 0) {
      iff=1;
      mj=MSEED-(*idum < 0 ? -*idum : *idum);
      mj %= MBIG;
      ma[55]=mj;
      mk=1;
      for (i=1;i<=54;i++) {
         ii=(21*i) % 55;
         ma[ii]=mk;
         mk=mj-mk;
         if (mk < MZ) mk += MBIG;
         mj=ma[ii];
      }
      for (k=1;k<=4;k++)
         for (i=1;i<=55;i++) {
            ma[i] -= ma[1+(i+30) % 55];
            if (ma[i] < MZ) ma[i] += MBIG;
         }
            inext=0;
      inextp=31;
      *idum=1;
   }
   if (++inext == 56) inext=1;
   if (++inextp == 56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if (mj < MZ) mj += MBIG;
   ma[inext]=mj;
   return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

#define NRANSI
#define M 7
#define NSTACK 50

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
   int *v;
   
   v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
   if (!v) nrerror("allocation failure in ivector()");
   return v-nl+NR_END;
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
   free((FREE_ARG) (v+nl-NR_END));
}


void indexx(unsigned long n, double arr[], unsigned long indx[])
{
   unsigned long i,indxt,ir=n,itemp,j,k,l=1;
   int jstack=0,*istack;
   double a;
   
   istack=ivector(1,NSTACK);
   for (j=1;j<=n;j++) indx[j]=j;
   for (;;) {
      if (ir-l < M) {
         for (j=l+1;j<=ir;j++) {
            indxt=indx[j];
            a=arr[indxt];
            for (i=j-1;i>=1;i--) {
               if (arr[indx[i]] <= a) break;
               indx[i+1]=indx[i];
            }
            indx[i+1]=indxt;
         }
         if (jstack == 0) break;
         ir=istack[jstack--];
         l=istack[jstack--];
      } else {
         k=(l+ir) >> 1;
         SWAP(indx[k],indx[l+1],itemp);
         if (arr[indx[l+1]] > arr[indx[ir]]) {
            SWAP(indx[l+1],indx[ir],itemp)
         }
         if (arr[indx[l]] > arr[indx[ir]]) {
            SWAP(indx[l],indx[ir],itemp)
         }
         if (arr[indx[l+1]] > arr[indx[l]]) {
            SWAP(indx[l+1],indx[l],itemp)
         }
         i=l+1;
         j=ir;
         indxt=indx[l];
         a=arr[indxt];
         for (;;) {
            do i++; while (arr[indx[i]] < a);
            do j--; while (arr[indx[j]] > a);
            if (j < i) break;
            SWAP(indx[i],indx[j],itemp)
         }
         indx[l]=indx[j];
         indx[j]=indxt;
         jstack += 2;
         if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
         if (ir-i+1 >= j-l) {
            istack[jstack]=ir;
            istack[jstack-1]=i;
            ir=j-1;
         } else {
            istack[jstack]=j-1;
            istack[jstack-1]=l;
            l=i;
         }
      }
   }
   free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */

#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);	\
a[k][l]=h+s*(g-h*tau);

void jacobi(double a[4][4], int n, double d[], double v[4][4], int *nrot)
{
   int    j,iq,ip,i;
   double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
   
   b=vector(1,n);
   z=vector(1,n);
   for (ip=1;ip<=n;ip++) {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
   }
   for (ip=1;ip<=n;ip++) {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
   }
   *nrot=0;
   for (i=1;i<=50;i++) {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) {
         for (iq=ip+1;iq<=n;iq++)
            sm += fabs(a[ip][iq]);
      }
      if (sm == 0.0) {
         free_vector(z,1,n);
         free_vector(b,1,n);
         return;
      }
      if (i < 4)
         tresh=0.2*sm/(n*n);
      else
         tresh=0.0;
      for (ip=1;ip<=n-1;ip++) {
         for (iq=ip+1;iq<=n;iq++) {
            g=100.0*fabs(a[ip][iq]);
            if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
                && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
               a[ip][iq]=0.0;
            else if (fabs(a[ip][iq]) > tresh) {
               h=d[iq]-d[ip];
               if ((double)(fabs(h)+g) == (double)fabs(h))
                  t=(a[ip][iq])/h;
               else {
                  theta=0.5*h/(a[ip][iq]);
                  t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                  if (theta < 0.0) t = -t;
               }
               c=1.0/sqrt(1+t*t);
               s=t*c;
               tau=s/(1.0+c);
               h=t*a[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               a[ip][iq]=0.0;
               for (j=1;j<=ip-1;j++) {
                  ROTATE(a,j,ip,j,iq)
               }
               for (j=ip+1;j<=iq-1;j++) {
                  ROTATE(a,ip,j,j,iq)
               }
               for (j=iq+1;j<=n;j++) {
                  ROTATE(a,ip,j,iq,j)
               }
               for (j=1;j<=n;j++) {
                  ROTATE(v,j,ip,j,iq)
               }
               ++(*nrot);
            }
         }
      }
      for (ip=1;ip<=n;ip++) {
         b[ip] += z[ip];
         d[ip]=b[ip];
         z[ip]=0.0;
      }
   }
   nrerror("Too many iterations in routine jacobi");
}

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
 int i,k;
 double p,qn,sig,un,*u;
 
 u=vector((long)1,(long)(n-1));
 if (yp1 > 0.99e30)
  y2[1]=u[1]=0.0;
 else {
  y2[1] = -0.5;
  u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
 }
 for (i=2;i<=n-1;i++) {
  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p=sig*y2[i-1]+2.0;
  y2[i]=(sig-1.0)/p;
  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
 }
 if (ypn > 0.99e30)
  qn=un=0.0;
 else {
  qn=0.5;
  un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
 }
 y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
 for (k=n-1;k>=1;k--)
  y2[k]=y2[k]*y2[k+1]+u[k];
 free_vector(u,(long)1,(long)(n-1));
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
 void nrerror(char error_text[]);
 int klo,khi,k;
 double h,b,a;
 
 klo=1;
 khi=n;
 while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
 }
 h=xa[khi]-xa[klo];
 if (h == 0.0) nrerror("Bad xa input to routine splint");
 a=(xa[khi]-x)/h;
 b=(x-xa[klo])/h;
 *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#undef ROTATE
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */


void fourn(flouble data[], unsigned long nn[], int ndim, int isign)
{
  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  double tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;
  
  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
        for (i1=i2;i1<=i2+ip1-2;i1+=2) {
          for (i3=i1;i3<=ip3;i3+=ip2) {
            i3rev=i2rev+i3-i2;
            SWAP(data[i3],data[i3rev],tempr);
            SWAP(data[i3+1],data[i3rev+1],tempr);
          }
        }
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
        for (i1=i3;i1<=i3+ip1-2;i1+=2) {
          for (i2=i1;i2<=ip3;i2+=ifp2) {
            k1=i2;
            k2=k1+ifp1;
            tempr=(double)wr*data[k2]-(double)wi*data[k2+1];
            tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
            data[k2]=data[k1]-tempr;
            data[k2+1]=data[k1+1]-tempi;
            data[k1] += tempr;
            data[k1+1] += tempi;
          }
        }
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software ,2kB. */



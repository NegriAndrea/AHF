#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../param.h"
#include "../tdef.h"

#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

/* define some MIN and MAX funtions */
#define MIN(A,B)        ((A)<(B)?(A):(B))
#define MAX(A,B)        ((A)>(B)?(A):(B))
#define FUNC(x)         ((*func)(x))
#define SWAP(a,b,temp)  {temp=(a);(a)=(b);(b)=temp;}
#define Re(ix,iy,iz,L)  (2*((iz)*(L)*(L) + (iy)*(L) + (ix))    )
#define Im(ix,iy,iz,L)  (2*((iz)*(L)*(L) + (iy)*(L) + (ix)) + 1)
#define FRAC(x)         (((x)>1.)?(1./(x)):(x))


/* faster way to calculate power of 2, 3 and 4 */
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define pow5(x) ((x)*(x)*(x)*(x)*(x))
#define pow6(x) ((x)*(x)*(x)*(x)*(x)*(x))


#include "general.h"
#include "specific.h"
#include "alloc_struct.h"
#include "cosmology.h"
#include "loadbalance.h"

#endif

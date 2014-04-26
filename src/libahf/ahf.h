#ifndef AHF_INCLUDED
#define AHF_INCLUDED

#include "../param.h"
#include "../tdef.h"

#ifdef AHF 

/*-------------------
 * visible functions
 *-------------------*/
void ahf_gridinfo (gridls *grid, int curgrid_no);
void ahf_halos    (gridls *grid);
int  ahf_profiles (void);


/*------------------------
 * structures used by ahf
 *------------------------*/
typedef struct {
   int num;              /* Number of times this colour appears                      */
   int name;             /* Name of the colour                                       */
} COLGATH;

typedef struct {
   int numReplace;       /* Number of rules to change the colour to the lowest value */
   int holder[NDIM];     /* Holds the NDIM possible colour values                    */
   int numCol;           /* Number of *different* colours that the nodes point to    */
   COLGATH *info;			 /* The initial information abut the surrounding nodes       */
} COLOURHOLDER;

typedef struct spatialRef {
   int   name;           /* Name of the colour                                      */
   int   refLevel;       /* The refinement level                                    */
	intXYZ periodic;	    /* tells the periodic nature of this refinement            */
   pqptr	cur_pquad;      /* Pointer to the pquad holding the cquad                  */
   cqptr cur_cquad;      /* Pointer to the cquad holding the nquad                  */ 
   cqptr icur_cquad;
   nqptr cur_nquad;      /* Pointer to the nquad holding the node                   */
   nqptr icur_nquad;
   nptr  cur_node;       /* Pointer to the starting node for this colour            */
   float	x,y,z;
   struct spatialRef *next;  /* Pointer to the next SPATIALREF                      */
} SPATIALREF;


SPATIALREF *spatialRefTail;

COLOURHOLDER	col;


/* internally used functions */
int	      checkColourInfo  (nptr tsc_box[3][3][3]);
void	      colourInfo       (nptr tsc_box[3][3][3]);
intXYZ	   testBound        (nptr tsc_box[3][3][3], int x, int y, int z);
int	      colourCompare    (const void *, const void *);
SPATIALREF* insertColour     ( SPATIALREF* newColour, int firstCOLOUR );
SPATIALREF* deleteColour     ( SPATIALREF* spatialRefHead, int name);

#endif /* AHF */

#endif

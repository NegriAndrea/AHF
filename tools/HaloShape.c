#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <xlocale.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"
#include "../src/libutility/utility.h"

#define PALL   -10   // not defined in param.h (included via utility.h!)

#define NLIMIT 10
#define OWAIN

/***********************************************************************
 *                          STRUCTURES  ETC.
 ***********************************************************************/
typedef struct PART *PARTPTR;
typedef struct PART
{
   double  pos[3];
   double  vel[3];
   double  mass;
   int     ptype;
} PART;  

typedef struct SHAPE *SHAPEPTR;
typedef struct SHAPE
{
   int     nbins;
   
   double *r;
   
   long   *npart;
   double *mass;
   double *ovdens;
   double *dens;
   double *v2_circ;
   
   /* standard inertia tensor */
   double *a, *b, *c;
   double *Eax, *Eay, *Eaz;
   double *Ebx, *Eby, *Ebz;
   double *Ecx, *Ecy, *Ecz;
   
   /* reduced inertia tensor */
   double *a_red, *b_red, *c_red;
   double *Eax_red, *Eay_red, *Eaz_red;
   double *Ebx_red, *Eby_red, *Ebz_red;
   double *Ecx_red, *Ecy_red, *Ecz_red;
} SHAPE;

/***********************************************************************
 *                          COMMON VARIABLES
 ***********************************************************************/
double    Xc, Yc, Zc, BoxSize, BoxHalf;
long      Npart, Ndm, Ngas, Nstars;
PARTPTR   fst_part;

/***********************************************************************
 *                             FUNCTIONS
 ***********************************************************************/
void    read_particles       (char *infile);
void    calc_shapeprofile    (PARTPTR fst_ptype, long Ntype, SHAPEPTR shape);
void    write_shapeprofile   (FILE *fpout, SHAPE shape);
int     distance_sort        (const void *p0, const void *p1);
double  box_distance         (double U, double V);
void    create_profile       (SHAPEPTR shape, int nbins);
void    rotate               (double *x, double *y, double *z, double phi, double theta);
void    write_geomfile       (char *infile);
PARTPTR get_particles        (int ptype, int Ntype);

/***********************************************************************
 *                                MAIN
 ***********************************************************************/
int main(argc,argv)
int argc;
char **argv;
{  
   SHAPE    all, dm, gas, stars;
   PARTPTR  fst_pall, fst_pdm, fst_pgas, fst_pstars;
   char     outfile[2048];
   int      slen;
   FILE    *fpout;
   
   if(argc<6)
     {
      fprintf(stderr,"usage: %s particlefile Xhalo Yhalo Zhalo BoxSize\n",*argv);
      exit(1);
     }
   else
     {
      fprintf(stderr,"====================================================================\n");
      fprintf(stderr,"           calculate the shape profile of a given halo\n");
      fprintf(stderr,"====================================================================\n");
     }
   
   /*-------------------------------------------------------------------
    * init those COMMON variables using the command lines values
    *-------------------------------------------------------------------*/
   Xc      = (double) atof(argv[2]);
   Yc      = (double) atof(argv[3]);
   Zc      = (double) atof(argv[4]);
   BoxSize = (double) atof(argv[4]);
   BoxHalf = BoxSize/2.;
   
   /* come up with some name for the outfile (assuming "*.dat" for the particles input file!) */
   strcpy(outfile,argv[1]);
   slen=strlen(outfile);
   outfile[slen-4]='.';
   outfile[slen-3]='s';
   outfile[slen-2]='h';
   outfile[slen-1]='a';
   outfile[slen  ]='p';
   outfile[slen+1]='e';
   outfile[slen+2]='s';
   outfile[slen+3]='\0';
   fprintf(stderr,"o writing shape profiles to %s\n",outfile);
   if(!(fpout=fopen(outfile,"w")))
     {
      fprintf(stderr,"FATAL: could not open %s\n",outfile);
      exit(0);
     }
   else
      fclose(fpout); // we were just testing!
   
   
   /*-------------------------------------------------------------------
    * read particles into COMMON structure fst_part
    *-------------------------------------------------------------------*/
   read_particles(argv[1]);
   write_geomfile(argv[1]);
   
   /*-------------------------------------------------------------------
    * calculate shape profiles
    *-------------------------------------------------------------------*/
   if(Npart > NLIMIT)
     {
      fprintf(stderr,"o calculating shape profile for ptype=%4d ",PALL);
      fst_pall   = get_particles(PALL, Npart);
      calc_shapeprofile(fst_pall,   Npart,  &all);
      fprintf(stderr,"done\n");
     }
   else
      all.nbins = 0;
   
   if(Ndm > NLIMIT)
     {
      fprintf(stderr,"o calculating shape profile for ptype=%4d ",(int)PDM);
      fst_pdm    = get_particles(PDM, Ndm);
      calc_shapeprofile(fst_pdm,    Ndm,    &dm);
      fprintf(stderr,"done\n");
     }
   else
      dm.nbins = 0;
   
   if(Ngas > NLIMIT)
     {
      fprintf(stderr,"o calculating shape profile for ptype=%4d ",(int)PGAS);
      fst_pgas   = get_particles(PGAS, Ngas);
      calc_shapeprofile(fst_pgas,   Ngas,   &gas);
      fprintf(stderr,"done\n");
     }
   else
      gas.nbins = 0;
   
   if(Nstars > NLIMIT)
     {
      fprintf(stderr,"o calculating shape profile for ptype=%4d ",(int)PSTAR);
      fst_pstars = get_particles(PSTAR, Nstars);
      calc_shapeprofile(fst_pstars, Nstars, &stars);
      fprintf(stderr,"done\n");
     }
   else
      stars.nbins = 0;
   
   
   
   /*-------------------------------------------------------------------
    * write all shape profiles into the same file
    *-------------------------------------------------------------------*/
   fprintf(stderr,"o writing all profiles into single file %s ... ",outfile);
   fpout = fopen(outfile,"w");
   write_shapeprofile(fpout, all);
   fclose(fpout);
   
   // append the remaining profiles
   fpout = fopen(outfile,"a");
   write_shapeprofile(fpout, dm);
   write_shapeprofile(fpout, gas);
   write_shapeprofile(fpout, stars);
   fclose(fpout);
   fprintf(stderr,"done\n");
   
   
   
   
   /*-------------------------------------------------------------------
    * free all memory
    *-------------------------------------------------------------------*/
   // would need to write a routine that free's all arrays within each SHAPE structure,
   // but am too lazy at the moment ;-)
}


/***********************************************************************
 *                   calculate the shape profiles
 *      (NOTE: this routine assumes distance-ordered particles)
 ***********************************************************************/
void calc_shapeprofile(PARTPTR fst_ptype, long Ntype, SHAPEPTR shape)
{
   int     nbins, ibin;
   PARTPTR cur_part;
   double  dX, dY, dZ;
   double  dist_min,   dist_max;
   double  ldist_min, ldist_max, ldr;
   double  cur_dist, cur_rad, lcur_rad, rad_prev, rad_mid;
   double  weight, M_sphere, dM, M_prev;
   double  Volume, dV, V_prev;
   long    npart, n_prev;
	double  itensor[3][3];
	double  a11, a22, a33, a12, a13, a23, axis1, axis2, axis3;
	double  itensor_red[3][3];
	double  a11_red, a22_red, a33_red, a12_red, a13_red, a23_red, axis1_red, axis2_red, axis3_red;
   
   /* how many bins should be used for cur_halo */
   nbins   = (int) (6.2*(log10((double)Ntype))-3.5);
   nbins   = MAX(1,nbins);
   shape->nbins = nbins;
   
   /* allocate space for profile */
   create_profile(shape, nbins);
   
   /* minimum distance of current particle species */
   cur_part = fst_ptype+NLIMIT/2;
   dX       = box_distance(cur_part->pos[X], Xc);
   dY       = box_distance(cur_part->pos[Y], Yc);
   dZ       = box_distance(cur_part->pos[Z], Zc);
   dist_min = sqrt(pow2(dX)+pow2(dY)+pow2(dZ));
   
   /* maximum distance of current particle species */
   cur_part = fst_ptype+(Ntype-1);
   dX       = box_distance(cur_part->pos[X], Xc);
   dY       = box_distance(cur_part->pos[Y], Yc);
   dZ       = box_distance(cur_part->pos[Z], Zc);
   dist_max = sqrt(pow2(dX)+pow2(dY)+pow2(dZ));
   
   fprintf(stderr,"using %12ld particles and %4d bins from %e to %e... ",Ntype,nbins,dist_min*1000,dist_max*1000);
   
   /* logarithmical binning from dist_min to dist_max */
	ldist_min = log10(dist_min);
	ldist_max = log10(dist_max);
	ldr       = (ldist_max - ldist_min) / (double)nbins;
   
   
	/* reset cumulative counters */
   cur_part = fst_ptype;
	cur_dist = -1.0;
   npart    = 0;
   M_sphere = 0.;
   a11_red  = 0.;
   a22_red  = 0.;
   a33_red  = 0.;
   a12_red  = 0.;
   a13_red  = 0.;
   a23_red  = 0.;
   a11      = 0.;
   a22      = 0.;
   a33      = 0.;
   a12      = 0.;
   a13      = 0.;
   a23      = 0.;
   
   /* loop over all bins */
   for (ibin = 0; ibin < nbins; ibin++) 
     {
      /* get current outer radius using logarithmic radial bins */
      lcur_rad = ldist_min + ((double)ibin + 1) * ldr;
      cur_rad  = pow(10., lcur_rad);
      
      /* collect all particles up to cur_rad */
      while (cur_dist < cur_rad && cur_part < fst_ptype+Ntype)
        {
         /* increment number of particles counter */
         npart++;
         
         /*----------------------
          * calculate properties
          *----------------------*/
         weight = cur_part->mass;
         
         M_sphere += weight;
         
         /* put particle into halo rest frame */
         dX = box_distance(cur_part->pos[X], Xc);
         dY = box_distance(cur_part->pos[Y], Yc);
         dZ = box_distance(cur_part->pos[Z], Zc);
                  
         /* pow2(distance) of current particle */
         cur_dist = (pow2(dX) + pow2(dY) + pow2(dZ));
         
         /* reduced inertia tensor of all particles within sphere(!) */
         if(cur_dist > MACHINE_ZERO)
           {
            a11_red += weight * dX * dX / (cur_dist);
            a22_red += weight * dY * dY / (cur_dist);
            a33_red += weight * dZ * dZ / (cur_dist);
            a12_red += weight * dX * dY / (cur_dist);
            a13_red += weight * dX * dZ / (cur_dist);
            a23_red += weight * dY * dZ / (cur_dist);
           }
         
         /* inertia tensor of all particles within sphere(!) */
         a11 += weight * dX * dX;
         a22 += weight * dY * dY;
         a33 += weight * dZ * dZ;
         a12 += weight * dX * dY;
         a13 += weight * dX * dZ;
         a23 += weight * dY * dZ;
         
         /* do the sqrt() here... */
         cur_dist = sqrt(cur_dist);
         
         /* jump to next particle */
         cur_part++;
        } /* while(cur_dist < cur_rad && cur_part < fst_part+Npart) */
      
      /* get eigenavalues of inertia tensor (all part's within sphere!) */
      itensor[0][0] = a11;
      itensor[1][1] = a22;
      itensor[2][2] = a33;
      itensor[0][1] = a12;
      itensor[1][0] = a12;                   // standard inertia tensor
      itensor[0][2] = a13;
      itensor[2][0] = a13;
      itensor[1][2] = a23;
      itensor[2][1] = a23;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      itensor_red[0][0] = a11_red;
      itensor_red[1][1] = a22_red;
      itensor_red[2][2] = a33_red;
      itensor_red[0][1] = a12_red;
      itensor_red[1][0] = a12_red;          // reduced inertia tensor
      itensor_red[0][2] = a13_red;
      itensor_red[2][0] = a13_red;
      itensor_red[1][2] = a23_red;
      itensor_red[2][1] = a23_red;
      get_axes(itensor_red, &axis1_red, &axis2_red, &axis3_red);
      
      
      /* volume of sphere [0, cur_rad] */
      Volume = 4. * PI / 3. * pow3(cur_rad);
      
      /* differential values */
      dM = M_sphere - M_prev; /* mass   in shell [prev_rad, cur_rad] */
      dV = Volume   - V_prev; /* volume of shell [prev_rad, cur_rad] */
      
      /* store values in shape profile */
      rad_mid              = (cur_rad + rad_prev) / 2.;
      shape->r[ibin]       = cur_rad;
      shape->npart[ibin]   = npart;
      shape->mass[ibin]    = M_sphere;
      shape->ovdens[ibin]  = M_sphere / Volume;
      shape->dens[ibin]    = dM       / dV;
      shape->v2_circ[ibin] = M_sphere / cur_rad;
      shape->Eax[ibin]     = itensor[0][0];
      shape->Eay[ibin]     = itensor[1][0];
      shape->Eaz[ibin]     = itensor[2][0];
      shape->Ebx[ibin]     = itensor[0][1];
      shape->Eby[ibin]     = itensor[1][1];
      shape->Ebz[ibin]     = itensor[2][1];
      shape->Ecx[ibin]     = itensor[0][2];
      shape->Ecy[ibin]     = itensor[1][2];
      shape->Ecz[ibin]     = itensor[2][2];
      shape->a[ibin]       = 1.0;
      shape->b[ibin]       = sqrt(axis2 / axis1);
      shape->c[ibin]       = sqrt(axis3 / axis1);
      shape->Eax_red[ibin] = itensor_red[0][0];
      shape->Eay_red[ibin] = itensor_red[1][0];
      shape->Eaz_red[ibin] = itensor_red[2][0];
      shape->Ebx_red[ibin] = itensor_red[0][1];
      shape->Eby_red[ibin] = itensor_red[1][1];
      shape->Ebz_red[ibin] = itensor_red[2][1];
      shape->Ecx_red[ibin] = itensor_red[0][2];
      shape->Ecy_red[ibin] = itensor_red[1][2];
      shape->Ecz_red[ibin] = itensor_red[2][2];
      shape->a_red[ibin]   = 1.0;
      shape->b_red[ibin]   = sqrt(axis2_red / axis1_red);
      shape->c_red[ibin]   = sqrt(axis3_red / axis1_red);
      
      /* store old values */
      M_prev   = M_sphere;
      V_prev   = Volume;
      n_prev   = npart;
      rad_prev = cur_rad;
      
     } /* ibin loop */
   
}


/***********************************************************************
 *                  read the halo particles from file
 ***********************************************************************/
void write_shapeprofile(FILE *fpout, SHAPE shape)
{
   int ibin;
   
   fprintf(fpout,"%d\n",shape.nbins);
   for(ibin=0; ibin<shape.nbins; ibin++)
     {
      fprintf(fpout,"%e %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
              shape.r[ibin]*1000,                          // 1
              shape.npart[ibin],                           // 2
              shape.mass[ibin],                            // 3
              shape.ovdens[ibin] / (0.24*rhoc0),           // 4
              shape.dens[ibin] / (0.24*rhoc0),             // 5
              sqrt(shape.v2_circ[ibin]),                   // 6
              shape.a[ibin],                               // 7
              shape.b[ibin],                               // 8
              shape.c[ibin],                               // 9
              shape.Eax[ibin],                             // 10
              shape.Eay[ibin],                             // 11
              shape.Eaz[ibin],                             // 12
              shape.Ebx[ibin],                             // 13
              shape.Eby[ibin],                             // 14
              shape.Ebz[ibin],                             // 15
              shape.Ecx[ibin],                             // 16
              shape.Ecy[ibin],                             // 17
              shape.Ecz[ibin],                             // 18
              shape.a_red[ibin],                           // 19
              shape.b_red[ibin],                           // 20
              shape.c_red[ibin],                           // 21
              shape.Eax_red[ibin],                         // 22
              shape.Eay_red[ibin],                         // 23
              shape.Eaz_red[ibin],                         // 24
              shape.Ebx_red[ibin],                         // 25
              shape.Eby_red[ibin],                         // 26
              shape.Ebz_red[ibin],                         // 27
              shape.Ecx_red[ibin],                         // 28
              shape.Ecy_red[ibin],                         // 29
              shape.Ecz_red[ibin]                          // 30
              );
     }
}


/***********************************************************************
 *                  read the halo particles from file
 ***********************************************************************/
void read_particles(char *infile)
{
   FILE *fpin;
   char  indata[2048];
   long  ipart, ID;
   PARTPTR cur_part;
   
   if(!(fpin=fopen(infile,"r")))
     {
      fprintf(stderr,"FATAL: could not open %s\nAborting\n",infile);
      exit(0);
     }
   
   /* count number of particles in file */
   Npart = 0;
   while(fgets(indata,2048,fpin) != NULL)
      Npart++;
   fprintf(stderr,"o found %ld particles ... ",Npart);
   
   /* allocate memory for all those particles */
   fst_part = (PARTPTR) calloc(Npart, sizeof(PART));
   
   /* eventually read particles*/
   rewind(fpin);
   Ndm      = 0;
   Ngas     = 0;
   Nstars   = 0;
   cur_part = fst_part;
   for(ipart=0; ipart<Npart; ipart++)
     {
      fgets(indata,2048,fpin);
      
#ifdef OWAIN
      sscanf(indata,"%d %lf %lf %lf %lf", &(cur_part->ptype), &(cur_part->mass), &(cur_part->pos[X]), &(cur_part->pos[Y]), &(cur_part->pos[Z]));
      if (cur_part->ptype == 0)
         cur_part->ptype = PGAS;
      else if (cur_part->ptype == 1)
         cur_part->ptype = PDM;
      else if (cur_part->ptype == 2)
         cur_part->ptype = PSTAR;
#else
      sscanf(indata,"%ld %d %lf %lf %lf %lf",&ID, &(cur_part->ptype), &(cur_part->mass), &(cur_part->pos[X]), &(cur_part->pos[Y]), &(cur_part->pos[Z]));
      cur_part->ptype  = -cur_part->ptype;  // the PDM, PGAS, etc. identifiers are all negative
#endif
      
      cur_part->pos[X] /= 1000;
      cur_part->pos[Y] /= 1000;             // scale back to Mpc/h
      cur_part->pos[Z] /= 1000;
      
      /* the shape should be independent of any rotation of the object! */
      //rotate(&cur_part->pos[X], &cur_part->pos[Y], &cur_part->pos[Z], -63.*PI/180., -41.*PI/180.);
            
      /* count individual species */
      if(cur_part->ptype == PDM)
         Ndm++;
      else if(cur_part->ptype == PGAS)
         Ngas++;
      else if(cur_part->ptype == PSTAR)
         Nstars++;
      
      /* move to next particle */
      cur_part++;
     }
   
   /* double check particle content of halo */
   if(Npart != Ndm+Ngas+Nstars)
     {
      fprintf(stderr," FATAL: %ld != %ld+%ld+%ld\n",Npart,Ndm,Ngas,Nstars);
      //exit(0);
     }
   
   fprintf(stderr,"and read them successfully\n");
   fclose(fpin);
   
   
   /* qsort particles according to distance */
   qsort((void *)fst_part, Npart, sizeof(PART), distance_sort);
}

/***********************************************************************
 *    a little routine that calculates the distance between U and V 
 *           taking into account period boundary conditions
 ***********************************************************************/
double box_distance(double U, double V)
{
   double D;
   
   D = U-V;
   if (D >  BoxHalf) D -= BoxSize;
   if (D < -BoxHalf) D += BoxSize;
   
   return(D);
}

/***********************************************************************
 *       needed with qsort() in order to sort the particles
 *                    with respects to distance
 ***********************************************************************/
int distance_sort(const void *p0, const void *p1)
{
   double dx[2], dy[2], dz[2];
   double d2[2];
   
   dx[0] = box_distance(((PARTPTR)p0)->pos[X], Xc);
   dy[0] = box_distance(((PARTPTR)p0)->pos[Y], Yc);
   dz[0] = box_distance(((PARTPTR)p0)->pos[Z], Zc);
   dx[1] = box_distance(((PARTPTR)p1)->pos[X], Xc);
   dy[1] = box_distance(((PARTPTR)p1)->pos[Y], Yc);
   dz[1] = box_distance(((PARTPTR)p1)->pos[Z], Zc);
   
   d2[0] = pow2(dx[0])+pow2(dy[0])+pow2(dz[0]);
   d2[1] = pow2(dx[1])+pow2(dy[1])+pow2(dz[1]);
   
   if(d2[0]<d2[1])
      return(-1);
   else if(d2[0]>d2[1])
      return(+1);
   else
      return(0);
   
   
}

/***********************************************************************
 *         allocates arrays within shape structure for nbins
 ***********************************************************************/
void create_profile(SHAPEPTR shape, int nbins)
{
   shape->r       = (double *) calloc(nbins, sizeof(double));
   shape->npart   = (long *)   calloc(nbins, sizeof(long));
   shape->mass    = (double *) calloc(nbins, sizeof(double));
   shape->ovdens  = (double *) calloc(nbins, sizeof(double));
   shape->dens    = (double *) calloc(nbins, sizeof(double));
   shape->v2_circ = (double *) calloc(nbins, sizeof(double));
   shape->a       = (double *) calloc(nbins, sizeof(double));
   shape->b       = (double *) calloc(nbins, sizeof(double));
   shape->c       = (double *) calloc(nbins, sizeof(double));
   shape->Eax     = (double *) calloc(nbins, sizeof(double));
   shape->Eay     = (double *) calloc(nbins, sizeof(double));
   shape->Eaz     = (double *) calloc(nbins, sizeof(double));
   shape->Ebx     = (double *) calloc(nbins, sizeof(double));
   shape->Eby     = (double *) calloc(nbins, sizeof(double));
   shape->Ebz     = (double *) calloc(nbins, sizeof(double));
   shape->Ecx     = (double *) calloc(nbins, sizeof(double));
   shape->Ecy     = (double *) calloc(nbins, sizeof(double));
   shape->Ecz     = (double *) calloc(nbins, sizeof(double));
   shape->a_red   = (double *) calloc(nbins, sizeof(double));
   shape->b_red   = (double *) calloc(nbins, sizeof(double));
   shape->c_red   = (double *) calloc(nbins, sizeof(double));
   shape->Eax_red = (double *) calloc(nbins, sizeof(double));
   shape->Eay_red = (double *) calloc(nbins, sizeof(double));
   shape->Eaz_red = (double *) calloc(nbins, sizeof(double));
   shape->Ebx_red = (double *) calloc(nbins, sizeof(double));
   shape->Eby_red = (double *) calloc(nbins, sizeof(double));
   shape->Ebz_red = (double *) calloc(nbins, sizeof(double));
   shape->Ecx_red = (double *) calloc(nbins, sizeof(double));
   shape->Ecy_red = (double *) calloc(nbins, sizeof(double));
   shape->Ecz_red = (double *) calloc(nbins, sizeof(double));
}

/***********************************************************************
 *     copy particles of a given type over to a newly created array
 ***********************************************************************/
PARTPTR get_particles(int ptype, int Ntype)
{
   PARTPTR fst_ptype, cur_part;
   long    ipart;
   
   /* and now fill fst_ptype array with requested particle type */
   if(ptype == PALL)
     {
      /* as simple as this */
      fst_ptype = fst_part;
     }
   else
     {
      /* make room for particles */
      fst_ptype = (PARTPTR) calloc(Ntype, sizeof(PART));
      
      /* loop over all particles */
      ipart = 0;
      for(cur_part=fst_part; cur_part<fst_part+Npart; cur_part++)
        {
         /* we found a particle of the requested type*/
         if(cur_part->ptype == ptype)
           {
            /* copy all values over */
            (fst_ptype+ipart)->ptype  = cur_part->ptype;
            (fst_ptype+ipart)->mass   = cur_part->mass;
            (fst_ptype+ipart)->pos[X] = cur_part->pos[X];
            (fst_ptype+ipart)->pos[Y] = cur_part->pos[Y];
            (fst_ptype+ipart)->pos[Z] = cur_part->pos[Z];
            (fst_ptype+ipart)->vel[X] = cur_part->vel[X];
            (fst_ptype+ipart)->vel[Y] = cur_part->vel[Y];
            (fst_ptype+ipart)->vel[Z] = cur_part->vel[Z];
            
            /* next ptype particle */
            ipart++;
           }
        }
      
      /* ipart should be identical to the number of previously calculated Ntype */
      if(ipart != Ntype)
         fprintf(stderr,"o something strange happens: ipart=%ld != Ntype=%ld\n",ipart,Ntype);
     }
   
   return(fst_ptype);
}


/***********************************************************************
 *                     simply rotate a subhalo
 ***********************************************************************/
void rotate(double *x, double *y, double *z, double phi, double theta)
{
   double xx, yy, zz;
   
   xx = box_distance(*x,Xc) * cos(phi) - box_distance(*y,Yc) * sin(phi);
   yy = box_distance(*x,Xc) * sin(phi) + box_distance(*y,Yc) * cos(phi);
   zz = box_distance(*z,Zc);
   
   *x = fmod(xx                                  + Xc + BoxSize, BoxSize);
   *y = fmod(yy * cos(theta) - zz * sin(theta)   + Yc + BoxSize, BoxSize);
   *z = fmod(yy * sin(theta) + zz * cos(theta)   + Zc + BoxSize, BoxSize);
}

/***********************************************************************
 *                     for double-checking the halo
 ***********************************************************************/
void write_geomfile(char *infile)
{
   char outfile[2048];
   int  slen;
   FILE *fpout;
   PARTPTR cur_part;
   
   strcpy(outfile,infile);
   slen=strlen(outfile);
   outfile[slen-4]='.';
   outfile[slen-3]='g';
   outfile[slen-2]='e';
   outfile[slen-1]='o';
   outfile[slen  ]='m';
   outfile[slen+1]='\0';
   
   fpout = fopen(outfile,"w");
   
   for(cur_part=fst_part; cur_part<fst_part+Npart; cur_part++)
     {
      if(cur_part->ptype==PDM)
         fprintf(fpout,"P %e %e %e 0 0 0 2\n",cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
      if(cur_part->ptype==PGAS)
         fprintf(fpout,"P %e %e %e 1 0 0 3\n",cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
      if(cur_part->ptype==PSTAR)
         fprintf(fpout,"P %e %e %e 0 0 1 4\n",cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z]);
     }   
   fclose(fpout);
}
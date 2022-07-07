/* Determine the gravitational potential of the given halo, at every
 * particle position (maybe option: how many particles shall be used?)
 * 
 * K. Warnick, Jan. 2006, AIP
 *
 * STATUS, UPDATES: 24.01.2006: just created, nothing works yet
 *		    25.01.2006: works fine, it seems!
 *
 */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAXPART 2000000

#define G 4.3006485e-9  /* [Mpc/Msun] [km^2/s^2] */

int main()
{
   FILE *fpart, *fhalos;
   FILE *fpot;

   double Phi0, Phi, mass, I_now, I_mid, I_prev, d_prev, dr, M_r;
   double x, y, z, d[MAXPART];
   double xc, yc, zc;
   double tmpv, Mvir, Box, zred, a, phi_fac;
   long int i, j, npart, readnpart, readnvpart;
   int halonum;
   char line[500], infile[100], infile2[100], outfile[100];
   
   long unsigned *idx;
   long int index;
   double *dist;
   
   void indexx(unsigned long n, double arr[], unsigned long indx[]);
   
   /* open file with halo particles + velocities */
   
   fprintf(stderr,"File with extracted particle positions: \n");
   scanf("%s",&infile);
   fprintf(stderr,"_halos-file (for centre): \n");
   scanf("%s",&infile2);
   fprintf(stderr,"Use files %s, %s.\n",infile, infile2);
   fprintf(stderr,"Number of this halo (0=host): \n");
   scanf("%d",&halonum);
   fprintf(stderr,"Box size: \n");
   scanf("%lf",&Box);
   fprintf(stderr,"current redshift: \n");
   scanf("%lf",&zred);
   fprintf(stderr,"Output file: \n");
   scanf("%s",&outfile);
   
   a = 1./(1.+zred);
   phi_fac = G/a;

   /* read halo centre and determine average particle mass */
   if((fhalos=fopen(infile2,"r"))==NULL)
   {
      printf("Cannot open %s\n", infile2);
      exit(1);
   }
   fgets(line,sizeof(line),fhalos);	/* First line contains comments only! */
   
   i=0.;
   while(i <= halonum)
   {
      fscanf(fhalos, "%ld %ld  %lf %lf %lf  %lf %lf %lf  %lf", &readnpart, &readnvpart, &xc, &yc, &zc, &tmpv, &tmpv, &tmpv, &Mvir);
      if(feof(fhalos)) 
      {
	  fprintf(stderr,"Could not find halo number %d -- stop.\n", halonum); 
	  exit(1);
      }
      /* Read rest of line (there might be more data than I am reading here!) */
      fgets(line,sizeof(line),fhalos);
      i=i+1;
   }

   /* mean mass of all particles */
   mass = Mvir/readnvpart;
   fprintf(stderr,"Mean mass of all particles: %g\n", mass);      
   fprintf(stderr,"Halo centre was found at: %lf %lf %lf\n", xc, yc, zc);

   if((fpart=fopen(infile,"r"))==NULL)
   {
      printf("Cannot open %s\n", infile);
      exit(1);
   }

   /* read particle positions + velocities */
   i=0;
   while(1)
   {
      fscanf(fpart, "%lf %lf %lf", &x, &y, &z);
      if(feof(fpart)) 
         break;
      /* Read rest of line (there might be more data than I am reading here!) */
      fgets(line,sizeof(line),fpart);

      /* calculate positions in the halo rest frame: */
      x = x - xc;
      y = y - yc;
      z = z - zc;

      /* periodic boundary conditions */
      if (x >  0.5*Box) x -= Box;
      if (y >  0.5*Box) y -= Box;
      if (z >  0.5*Box) z -= Box;
      if (x < -0.5*Box) x += Box;
      if (y < -0.5*Box) y += Box;
      if (z < -0.5*Box) z += Box;

      d[i] = sqrt(x*x+y*y+z*z);

      i=i+1;
   }
   npart=i;	
   /* because at end of file: last line is read, i is incremented, end of file is found, break: now i is one more than there are values! */
   fprintf(stderr,"%ld particles read.\n", npart);
   fprintf(stderr,"%ld particles predicted by halos-file.\n", readnpart);
   fprintf(stderr,"d: %lf, %lf, %lf\n", d[0], d[1], d[2]);
   

   /* now I have all particles stored in the halo rest frame */
   /* just in case, they haven't been sorted properly (rounding errors?): sort them! */    
   
   /*fprintf(stderr,"first values: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14]);
   */
   
   /* sort according to distance */
   dist = (double *) calloc(npart+1, sizeof(double));
   idx = (long unsigned *) calloc(npart+1, sizeof(long unsigned));
   /* probably need npart+1 because it sorts values from 1 to npart! */
   
   for(i=0; i<npart; i++)
   {
      dist[i+1] = d[i];
   }
   
   indexx(npart,dist,idx);

   /* create ordered list */
   for(i=0; i<npart; i++)
   {
      index       = idx[i+1];
      d[i]	  = dist[index];
   }

   free(idx);
   free(dist);
   
   /*fprintf(stderr,"after sort: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14]);
   */
   
   /* determine Phi0, see HaloTracker.c */
   /* ==================================*/
   fprintf(stderr, "Determine Phi0 ...");

   M_r    = 0.0;
   Phi0   = 0.0;
   I_prev = 0.0;
   d_prev = 0.0;

   /* loop over all sorted particles from inside out */
   for(i=0; i<npart; i++)
   {
   
      /* cumulative mass */
      M_r   = M_r+mass;
      
      /* midpoint integration */
      I_now = M_r / (d[i]*d[i]);
      I_mid = 0.5 * (I_now + I_prev);
      dr    = d[i] - d_prev;
      Phi0  = Phi0 + I_mid*dr;

      /* store values for next integration step */
      d_prev    = d[i];
      I_prev    = I_now;   
      
      //fprintf(stderr,"Phi0: %g\n", Phi0);
   }
   
   /* finally: Phi0:*/
   Phi0 = M_r/d[npart-1] + Phi0; 
   fprintf(stderr, " done. Phi0 = %g \n", Phi0*phi_fac);
   fprintf(stderr, " phi_fac = %g \n", phi_fac);
   


   /* determine Phi */
   /*===============*/

   fprintf(stderr,"Determine Phi ...");
   M_r    = 0.0;
   Phi    = 0.0;
   I_prev = 0.0;
   d_prev = 0.0;


   /* open output file */
   if((fpot=fopen(outfile,"w"))==NULL)
   {
      printf("Cannot open %s\n", outfile);
      exit(1);
   }
   fprintf(fpot, "# d [Mpc/h]      Phi-Phi0 [km^2/s^2, physical units]\n");
   

   /* loop over all sorted particles from inside out */
   for(i=0; i<npart; i++)
   {
      /* cumulative mass */
      M_r   = M_r+mass;
      
      /* midpoint integration */
      I_now = M_r / (d[i]*d[i]);
      I_mid = 0.5 * (I_now + I_prev);
      dr    = d[i] - d_prev;
      
      Phi  = Phi + I_mid*dr;

      /* write current distance and Phi-Phi0 to file, in physical units! */
      /* only write every 20th value*/
      if (i*0.05 - (int)(i*0.05) <= 1.e-5) 
      {
         fprintf(fpot, "%14g  %14g\n", d[i], (Phi-Phi0)*phi_fac);
      }
      /* store values for next integration step */
      d_prev    = d[i];
      I_prev    = I_now;   
   }
   fprintf(stderr," done.\n");
   fprintf(stderr,"Data written to file %s\n", outfile);
   fprintf(stderr,"Potential in physical units, distance in co-moving coordinates! %s\n", outfile); 

}


/* sorting routine from utility.c (AMIGA), originally from Num. Rec.*/
#define NR_END 1
#define FREE_ARG char*

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define NSTACK 50
#define M 7
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}

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
         SWAP(indx[k],indx[l+1]);
         if (arr[indx[l+1]] > arr[indx[ir]]) {
            SWAP(indx[l+1],indx[ir])
         }
         if (arr[indx[l]] > arr[indx[ir]]) {
            SWAP(indx[l],indx[ir])
         }
         if (arr[indx[l+1]] > arr[indx[l]]) {
            SWAP(indx[l+1],indx[l])
         }
         i=l+1;
         j=ir;
         indxt=indx[l];
         a=arr[indxt];
         for (;;) {
            do i++; while (arr[indx[i]] < a);
            do j--; while (arr[indx[j]] > a);
            if (j < i) break;
            SWAP(indx[i],indx[j])
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



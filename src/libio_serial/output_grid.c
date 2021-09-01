#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"


/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

#define ASSIGN_DMDENS  // assigns density to grid prior to writing output

#ifdef ASSIGN_DMDENS
#include "../libamr_serial/amr_serial.h"
#endif

/*=============================================================================
 * dump the full quad information to file
 *=============================================================================*/
void output_grid(gridls *cur_grid, int dumpflag)
{
   char          filename[MAXSTRING], file_no[MAXSTRING], appendix[MAXSTRING];
   FILE          *outstream, *dbfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          x, y, z;
   long          no_quad;
   int           machine_sizeof_long;

   /* copy timestep dependent quantities over to io.header */
   io.header.no_timestep = global.no_timestep;
   io.header.a_current   = global.a;
   io.header.K_current   = energy.K_current;
   io.header.U_current   = energy.U_current;
   io.header.Eintegral   = energy.integral;
   io.header.Econst      = energy.econst;

#ifdef ASSIGN_DMDENS
   zero_dens(cur_grid);
   assign_dens(cur_grid);
#endif
   
   
   /* open file */
   if(dumpflag == 0)
     {
      //fprintf(io.logfile,"\nstarting to write grid file...");
      sprintf(file_no,"z%.3f",fabs(global.z));
#ifdef DEBUG_GRIDS
      strcpy(filename, "DEBUG_GRIDS-");
#else
      strcpy(filename, io.outfile_prefix);
#endif
      strcat(filename, file_no);
     }
   else
     {
      //fprintf(io.logfile,"\nstarting to dump grid file...");
      strcpy(filename, io.dumpfile_name);
     }
      
   fflush(io.logfile);
   sprintf(appendix,".grid-%09ld",cur_grid->l1dim);
   strcat(filename, appendix);

   if ((outstream = fopen(filename,"wb")) == NULL) 
     {
      fprintf(io.logfile,"output: could not open file %s\n", filename);
      exit(1);
     }


   /* write a simple "1" to file (for BYTESWAP testing) */
   machine_sizeof_long = sizeof(long);
   fwrite(&machine_sizeof_long, sizeof(int), 1, outstream);

   
   /* write io.header even into grid file... */
   if(fwrite(&(io.header), sizeof(io.header), 1, outstream) != 1)
     {
      fprintf(io.logfile,"\n\noutput_grid: could not write io.header\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }

   
   /* how many linked pquads ? */
   no_quad = 0;
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
      no_quad++;
   fwrite((void*)&no_quad, sizeof(long), 1, outstream);
#ifdef DEBUG_IO_GRIDS
   fprintf(stderr,"%ld -> no_pquad = %ld\n",cur_grid->l1dim, no_quad);
#endif
   
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
      fwrite((void*)&cur_pquad->z,      sizeof(int),  1,outstream);
      fwrite((void*)&cur_pquad->length, sizeof(int),  1,outstream);
#ifdef DEBUG_IO_GRIDS
      fprintf(stderr,"%ld -> z = %ld   length = %ld\n",cur_grid->l1dim, cur_pquad->z, cur_pquad->length);
#endif
      
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
        {
         /* how many linked cquads ? */
         no_quad = 0;
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
            no_quad++;
         fwrite((void*)&no_quad, sizeof(long), 1, outstream);
#ifdef DEBUG_IO_GRIDS
         fprintf(stderr,"%ld -> no_cquad = %ld\n",cur_grid->l1dim, no_quad);
#endif
         
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
           {
            fwrite((void*)&icur_cquad->y,      sizeof(int),  1,outstream);
            fwrite((void*)&icur_cquad->length, sizeof(int),  1,outstream);
#ifdef DEBUG_IO_GRIDS
            fprintf(stderr,"%ld -> y = %ld   length = %ld\n",cur_grid->l1dim, cur_cquad->y, cur_cquad->length);
#endif
            
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
              {
               /* how many linked nquads ? */
               no_quad = 0;
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                  no_quad++;
               fwrite((void*)&no_quad, sizeof(long), 1, outstream);
#ifdef DEBUG_IO_GRIDS
               fprintf(stderr,"%ld -> no_nquad = %ld\n",cur_grid->l1dim, no_quad);
#endif
               
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                 {
                  fwrite((void*)&icur_nquad->x,      sizeof(int),  1,outstream);
                  fwrite((void*)&icur_nquad->length, sizeof(int),  1,outstream);
#ifdef DEBUG_IO_GRIDS
                  fprintf(stderr,"%ld -> x = %ld   length = %ld\n",cur_grid->l1dim, cur_nquad->x, cur_nquad->length);
#endif
                  
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     fwrite((void*)&cur_node->dens, sizeof(flouble),      1, outstream);
                    }
                 }
              }
           }
        }
      
     }
      
   fclose(outstream);

   //fprintf(io.logfile,"done\n");
   //fflush(io.logfile);
}

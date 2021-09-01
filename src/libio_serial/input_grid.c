#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "io_serial.h"
#include "../libutility/utility.h"
#include "../libamr_serial/amr_serial.h"

/*=============================================================================
 * read the full quad information from file
 *=============================================================================*/
void input_grid(gridls *cur_grid)
{
   char          filename[MAXSTRING], file_no[MAXSTRING], appendix[MAXSTRING];
   FILE          *instream;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          x, y, z;
   long          no_pquad, no_cquad, no_nquad;
   long          ipquad, icquad, inquad;
   int           i, SWAPBYTES;
   info_io       io_grid;            

   /* open file */
   sprintf(file_no,"z%.3f",fabs(global.z));
   sprintf(appendix,".grid-%09ld",cur_grid->l1dim);
   strcpy(filename, io.outfile_prefix);
   //strcpy(filename, io.icfile_name);
   strcat(filename, file_no);
   strcat(filename, appendix);
   
   fprintf(stderr,"\n*******************************************************************************\n");
   fprintf(stderr,"                           starting input_grid()\n");
   fprintf(stderr,"*******************************************************************************\n");
   
   
   if ((instream = fopen(filename,"rb")) == NULL) 
     {
      fprintf(stderr,"input_grid: could not open file %s\n", filename);
      fprintf(io.logfile,"input_grid: could not open file %s\n", filename);
      exit(1);
     }

   /* read in io.header from file */
   read_amiga_header(instream, &io_grid, &SWAPBYTES);

   /* how many linked pquads ?  */
   ReadLong(instream, &no_pquad, SWAPBYTES);
#ifdef DEBUG_IO_GRIDS
   fprintf(stderr,"%ld -> no_pquad = %ld\n",cur_grid->l1dim, no_pquad);
#endif
   if(no_pquad == 0)
      cur_grid->pquad = NULL;
   else
      cur_grid->pquad = c_pquad(1);
   cur_pquad = cur_grid->pquad;
   
   /* loop over linked pquads */
   for(ipquad=0; ipquad<no_pquad; ipquad++)
     {
      ReadLong(instream, &cur_pquad->z,      SWAPBYTES);
      ReadLong(instream, &cur_pquad->length, SWAPBYTES);
      
#ifdef DEBUG_IO_GRIDS
      fprintf(stderr,"%ld -> z = %ld   length = %ld\n",cur_grid->l1dim, cur_pquad->z, cur_pquad->length);
#endif
      cur_pquad->loc = c_cquad(cur_pquad->length);
      
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
        {
         ReadLong(instream, &no_cquad, SWAPBYTES);
#ifdef DEBUG_IO_GRIDS
         fprintf(stderr,"%ld -> no_cquad = %ld\n",cur_grid->l1dim, no_cquad);
#endif
         /* the linked-list of cquad's starts at cur_cquad */
         icur_cquad = cur_cquad;
         
         /* loop over linked cquads */
         for(icquad=0; icquad<no_cquad; icquad++)
           {
            ReadLong(instream, &icur_cquad->y,      SWAPBYTES);
            ReadLong(instream, &icur_cquad->length, SWAPBYTES);
#ifdef DEBUG_IO_GRIDS
            fprintf(stderr,"%ld -> y = %ld   length = %ld\n",cur_grid->l1dim, icur_cquad->y, icur_cquad->length);
#endif
            
            icur_cquad->loc = c_nquad(icur_cquad->length);
            
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
              {
               ReadLong(instream, &no_nquad, SWAPBYTES);
#ifdef DEBUG_IO_GRIDS
               fprintf(stderr,"%ld -> no_nquad = %ld\n",cur_grid->l1dim, no_nquad);
#endif
               /* the linked-list of nquad's starts at cur_nquad */
               icur_nquad = cur_nquad;
               
               /* loop over linked nquads */
               for(inquad=0; inquad<no_nquad; inquad++)
                 {
                  ReadLong(instream, &icur_nquad->x,      SWAPBYTES);
                  ReadLong(instream, &icur_nquad->length, SWAPBYTES);
#ifdef DEBUG_IO_GRIDS
                  fprintf(stderr,"%ld -> x = %ld   length = %ld\n",cur_grid->l1dim, icur_nquad->x, icur_nquad->length);
#endif
                  
                  icur_nquad->loc = c_node(icur_nquad->length);
                  
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     if(io.header.double_precision)
                       {
                        ReadDouble(instream, (double *) &cur_node->dens, SWAPBYTES);
                       }
                     else
                       {
                        ReadFloat(instream, (float *) &cur_node->dens, SWAPBYTES);
                       }
                    }
                  
                  /* the last icur_nquad->next has to point to NULL */
                  if(inquad != no_nquad-1)
                    {
                     icur_nquad->next = c_nquad(1);
                     icur_nquad       = icur_nquad->next;
                    }
                  else
                    {
                     icur_nquad->next  = NULL;
                    }
                  
                 }
              }
            
            if(icquad != no_cquad-1)
              {
               icur_cquad->next = c_cquad(1);
               icur_cquad       = icur_cquad->next;
              }
            else
              {
               icur_cquad->next = NULL;
              }
            
           }
        }
      
      if(ipquad != no_pquad-1)
        {
         cur_pquad->next = c_pquad(1);
         cur_pquad       = cur_pquad->next;
        }
      else
        {
         cur_pquad->next = NULL;
        }
      
     }
      
   fclose(instream);

#ifdef VERBOSE
   fprintf(stderr,"*******************************************************************************\n");
   fprintf(stderr,"                           finished input_grid()\n");
   fprintf(stderr,"*******************************************************************************\n\n");
#endif
}

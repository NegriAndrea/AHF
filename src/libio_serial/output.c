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


/*================================================================================
 * write particles to a file (either dumpfile or output file)
 *================================================================================*/
void output(char *outfile_name, int dumpflag)
{
  partptr cur_part;                  /* current particle               */
  int     i,j;                       /* loop index                     */
  int     slen;                      /* string length                  */
  char    file_ncpy[MAXSTRING];      /* copy of file name              */
  char    file_no[10];               /* file number                    */
  char   *f_name;                    /* file name pointer              */
  FILE   *outstream;                 /* file pointer                   */
  int     machine_sizeof_long;       /* store sizeof(long) in outfile  */
  flouble pos, mom, fweight;
   
  
  /* copy timestep dependent quantities over to io.header */
  io.header.no_timestep = global.no_timestep;
  io.header.a_current   = global.a;
  io.header.K_current   = energy.K_current;
  io.header.U_current   = energy.U_current;
  io.header.Eintegral   = energy.integral;
  io.header.Econst      = energy.econst;
  
  
  if(dumpflag == 0)   /* normal output */
    {
     f_name = strcpy(file_ncpy, outfile_name);
	  
     /* OUTPUT CONVENTION: 3 digits*/
     sprintf(file_no,"z%.3f",fabs(global.z));
     strcat(f_name, file_no);
     
       /* open file */
       if ((outstream = fopen(f_name,"wb")) == NULL) 
       {
          fprintf(io.logfile,"output: could not open file %s\n", f_name);
          exit(1);
       }
    }
  else
    {
      if((outstream = fopen(outfile_name,"wb")) == NULL)
        {
         fprintf(io.logfile,"output: could not open dump file %s\n",
                 outfile_name);
         exit(1);
        }
    }
  
  
  /* write a simple "1" to file (for BYTESWAP testing) */
  machine_sizeof_long = sizeof(long);
  fwrite(&machine_sizeof_long, sizeof(int), 1, outstream);

  /* write IO header */
  if(fwrite(&(io.header), sizeof(io.header), 1, outstream) != 1)
    {
      fprintf(io.logfile,"\n\noutput: could not write io.header\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
    }

  cur_part = global.fst_part;
  for(i = 0; i < global.no_part; i++)
   {
     for(j = X; j <= Z; j++)
       {
        cur_part->pos[j] = fmod((double)cur_part->pos[j] + 1.0, 1.0);
        
        pos = cur_part->pos[j];
        mom = cur_part->mom[j];
        
        fwrite((void*)&pos, sizeof(flouble), 1, outstream);
        fwrite((void*)&mom, sizeof(flouble), 1, outstream);
       }
#ifdef MULTIMASS
     fweight  = cur_part->weight;
     fwrite((void*)&fweight, sizeof(flouble), 1, outstream);
#endif  // MULTIMASS

     /* move to next particle */
     cur_part++;
   }
  fflush(outstream);
  fclose(outstream);
}



/**
 * \file io_art_header.c
 *
 * Provides functions for reading and writing the header of ART
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
 \**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

#include "io_art_header.h"
#include "io_util.h"

/* taken from param.h */
#define rhoc0     2.7755397e11      /* [h^2*Msun]/[Mpc^3]      */

/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
 \**********************************************************************/
int dustbin;
#define SKIP(fp) {io_util_readint32(fp, &dustbin, f->swapped);}

/**********************************************************************\
 *    Prototypes of local functions                                   * 
 \**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
 \**********************************************************************/
extern io_art_header_t
io_art_header_get(io_logging_t log, io_art_t f)
{
  io_art_header_t dummy;
  long skipsize;
  float *wspecies;
  int iloop, *lspecies;
  FILE *fPMcrd, *fart;
  char artline[2048];
  char c;
  char PMcrd[2048];
  int  is, id, slen,i;
  id = is = slen = i = 0;

  
  /* Some sanity checks */
  if ((f == NULL) || (f->file == NULL))
    return NULL;
  
#ifndef MULTIMASS
  io_logging_fatal(log, "you need to switch on -DMULTIMASS even if you only have a one-species simulation!");
  return NULL;
#endif
  
  /* Check if there already is a header, do nothing then */
  if (f->header != NULL)
    return f->header;
  
  /* Create the header structure array */
  dummy = (io_art_header_t)malloc((size_t)sizeof(io_art_header_struct_t));
  if (dummy == NULL) {
    io_logging_memfatal(log, "ART header structure");
    return NULL;
  }
  
  /* open the file PMcrd.DAT */
  slen = strlen(f->fname);
  for(iloop=slen-1; iloop>0; iloop--)
   {
    if(f->fname[iloop] == 's')
     {
      is = iloop;
      break;
     }
   }
  if(is==0)
   {
    io_logging_fatal(log, "could not find the 's' in the input file name");
    return NULL;
   }
  for(iloop=is; iloop>0; iloop--)
   {
    if(f->fname[iloop] == '/')
     {
      id = iloop;
      break;
     }
   }
  
  if(id==0)
    id=-1;
  else
    for(iloop=0; iloop<=id; iloop++)
      PMcrd[iloop] = f->fname[iloop];
  
  PMcrd[id+1] = 'P';
  PMcrd[id+2] = 'M';
  PMcrd[id+3] = 'c';
  PMcrd[id+4] = 'r';
  PMcrd[id+5] = 'd';
  for(iloop=is+2; iloop<=slen; iloop++)
    PMcrd[id+6+iloop-(is+2)] = f->fname[iloop];
  PMcrd[iloop] = '\0';
  fprintf(stderr,"o reading header of size %lu from file %s (slen=%lu)\n",sizeof(io_art_header_struct_t),PMcrd,strlen(PMcrd));
  
  fPMcrd = fopen(PMcrd,"rb");
  if(fPMcrd == NULL)
   {
    io_logging_fatal(log, "could not open PMcrd.DAT with ART header information");
    return NULL;
   }
    
  /* Now start reading the header */
  SKIP(fPMcrd);
  for(iloop=0; iloop<45; iloop++)
   {
    fread(&c, sizeof(char), 1, fPMcrd);
    strcpy(&(dummy->header_string[iloop]), &c);
   }
  dummy->header_string[44] = '\0';
  fprintf(stderr,"o header string => %s\n",dummy->header_string);

  io_util_readfloat(fPMcrd, &(dummy->aexpn), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->aexp0), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->amplt), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->astep), f->swapped);
  io_util_readint32(fPMcrd, &(dummy->istep), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->partw), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->tintg), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->ekin), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->ekin1), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->ekin2), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->au0), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->aeu0), f->swapped);
  io_util_readint32(fPMcrd, &(dummy->nrowc), f->swapped);
  io_util_readint32(fPMcrd, &(dummy->ngridc), f->swapped);
  io_util_readint32(fPMcrd, &(dummy->nspecies), f->swapped);
  io_util_readint32(fPMcrd, &(dummy->nseed), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->Om0), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->Oml0), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->hubble), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->wp5), f->swapped);
  io_util_readfloat(fPMcrd, &(dummy->Ocurv), f->swapped);
  for(iloop=0; iloop<100; iloop++)
    io_util_readfloat(fPMcrd, &(dummy->extras[iloop]), f->swapped);
  
  
  wspecies =       &(dummy->extras[0]);
  lspecies = (int*)&(dummy->extras[10]);
  if(dummy->nspecies == 0)
   {
    dummy->N_particles = (dummy->nrowc)*(dummy->nrowc)*(dummy->nrowc);
    dummy->N_pages     = (dummy->nrowc);
    dummy->N_in_last   = (dummy->nrowc)*(dummy->nrowc);
   }
  else
   {
    dummy->N_particles =  lspecies[dummy->nspecies-1];
    dummy->N_pages     = (dummy->N_particles -1)/((dummy->nrowc)*(dummy->nrowc)) + 1;
    dummy->N_in_last   =  dummy->N_particles -((dummy->nrowc)*(dummy->nrowc))*(dummy->N_pages-1);
   }
  
  if(dummy->extras[100-1] > 0 && dummy->extras[59-1] > 0)
   {
    dummy->boxsize = dummy->extras[100-1];
    dummy->munit   = dummy->extras[59-1];
   }
  else
   {
    fprintf(stderr,"o reading ART units from file art.info\n");
    if((fart = fopen("art.info","r")) == NULL)
     {
      io_logging_fatal(log,"Could not open art.info containing the following information:");
      io_logging_fatal(log,"20.0                    ART_BOXSIZE           [Mpc/h]");
      io_logging_fatal(log,"2.2197e09               ART_MUNIT             [Msun/h]");
      io_logging_fatal(log,"Please generate this file. Aborting now!");
      free(dummy);
      return NULL;
     }
    fgets(artline, 2048, fart);
    sscanf(artline,"%lf",&(dummy->boxsize));
    fgets(artline, 2048, fart);
    sscanf(artline,"%lf",&(dummy->munit));  
    fclose(fart);    
   }
  
  
  
  /* make header information globally available */
  f->header = dummy;
  
  return dummy;
}

extern void
io_art_header_del(io_logging_t log, io_art_header_t *header)
{
  if ( (header == NULL) || (*header == NULL) )
    return;
  
  free(*header);
  
  *header = NULL;
  
  return;
}

extern void
io_art_header_write(io_logging_t log,
                    io_art_header_t header,
                    io_art_t f)
{
  if ( (header == NULL) || (f == NULL))
    return;
  
  if (f->mode != IO_FILE_WRITE)
    return;
  
  if (f->file == NULL)
    return;
  
  if (header != f->header) {
    io_logging_msg(log, INT32_C(1),
                   "Writing a different header than stored in "
                   "the file object to the file.");
  }
  
  /* TODO: Write the header */
  
  return;
}


extern void
io_art_header_log(io_logging_t log, io_art_header_t header)
{
  io_logging_msg(log, INT32_C(5),
                 "   aexpn:                        %e",
                 header->aexpn);
  io_logging_msg(log, INT32_C(5),
                 "   aexp0:                        %e",
                 header->aexp0);
  io_logging_msg(log, INT32_C(5),
                 "   astep:                        %e",
                 header->astep);
  io_logging_msg(log, INT32_C(5),
                 "   istep:                        %" PRIi32,
                 header->istep);
  io_logging_msg(log, INT32_C(5),
                 "   partw:                        %e",
                 header->partw);
  io_logging_msg(log, INT32_C(5),
                 "   nrowc:                        %" PRIi32,
                 header->nrowc);
  io_logging_msg(log, INT32_C(5),
                 "   ngridc:                       %" PRIi32,
                 header->ngridc);
  io_logging_msg(log, INT32_C(5),
                 "   nspecies:                     %" PRIi32,
                 header->nspecies);
  io_logging_msg(log, INT32_C(5),
                 "   Om0:                          %e",
                 header->Om0);
  io_logging_msg(log, INT32_C(5),
                 "   Oml0:                         %e",
                 header->Oml0);
  io_logging_msg(log, INT32_C(5),
                 "   extras[58] (pmass):           %e",
                 header->extras[58]);
  io_logging_msg(log, INT32_C(5),
                 "   extras[99] (boxsize):         %e",
                 header->extras[99]);
  
  
  return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
 \**********************************************************************/

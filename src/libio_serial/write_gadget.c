#if (defined MULTIMASS && defined GAS_PARTICLES)

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/
unsigned int blklen;
#define GADGET_WSKIP       fwrite(&blklen, sizeof(int), 1, fpout);

/*=============================================================================
 *                                STRUCTURES
 *=============================================================================*/
info_gadget gadget;

/*=============================================================================
 *                                PROTOTYPES
 *=============================================================================*/
void write_blockname(FILE *fpout, char *name);
int get_ngas(partptr fst_part, long npart);
int get_nstar(partptr fst_part, long npart);
int get_ndm_hr(partptr fst_part, long npart);


/*==============================================================================
 * write particles in Part[] as a GADGET2 file in float precision
 *==============================================================================*/
void write_gadget(char *fname)
{
  int ipart, itype, nwrite, id, ioffset;
  float f[3], m, u, rho, eps;
  int   i;
  int   ngas, nstar, ndm, ndm_hr, ndm_lr, ncheck_all, ncheck_gas, ncheck_dm, ncheck_dm_hr, ncheck_dm_lr, ncheck_star;
  FILE *fpout;
  double x_fac, v_fac, m_fac, rho_fac;
  partptr cur_part;
  
  x_fac   = simu.boxsize*1000;
  v_fac   = simu.boxsize / simu.t_unit / global.a / sqrt(global.a);
  m_fac   = simu.pmass/1e10;
  rho_fac = simu.pmass/pow3(simu.boxsize);
  
  fpout = fopen(fname,"wb");
  
  // number of particles in each category
  ngas   = get_ngas(global.fst_part, global.no_part);
  nstar  = get_nstar(global.fst_part, global.no_part);
  ndm_hr = get_ndm_hr(global.fst_part, global.no_part);
  ndm    = global.no_part-ngas-nstar;
  ndm_lr = ndm-ndm_hr;
  
  fprintf(stderr,"write_gadget: ndm=%d (ndmhr=%d ndmlr=%d) ngas=%d nstar=%d ntot=%ld\n",ndm,ndm_hr,ndm_lr,ngas,nstar,global.no_part);
  fprintf(stderr,"               -> writing to %s\n",fname);
  
  // fill header
  gadget.header.expansion    = global.a;
  gadget.header.redshift     = 1./global.a - 1.;
  if(nstar>0) {
    gadget.header.flagsfr      = 1;
    gadget.header.flagfeedback = 1;
  }
  else {
    gadget.header.flagsfr      = 0;
    gadget.header.flagfeedback = 0;
  }
  gadget.header.NumFiles     = 1;
  gadget.header.BoxSize      = x_fac;
  gadget.header.Omega0       = simu.omega0;
  gadget.header.OmegaLambda  = simu.lambda0;
  gadget.header.HubbleParam  = H0;
  gadget.header.np[0]        = ngas;
  gadget.header.np[1]        = ndm_hr;
  gadget.header.np[2]        = ndm_lr;
  gadget.header.np[3]        = 0;
  gadget.header.np[4]        = nstar;
  gadget.header.np[5]        = 0;
  gadget.header.nall[0]      = gadget.header.np[0];
  gadget.header.nall[1]      = gadget.header.np[1];
  gadget.header.nall[2]      = gadget.header.np[2];
  gadget.header.nall[3]      = gadget.header.np[3];
  gadget.header.nall[4]      = gadget.header.np[4];
  gadget.header.nall[5]      = gadget.header.np[5];
  
  // we will write a MASS block
  for(itype=0; itype<6; itype++) {
    gadget.header.massarr[itype] = 0.0;
  }
  
  // just to be sure
  rewind(fpout);
  
  //=============================================
  // HEADER
  //=============================================
  write_blockname(fpout, "HEAD");
  
  blklen = sizeof(gadget.header);
  GADGET_WSKIP;
  fwrite(&(gadget.header), sizeof(gadget.header), 1, fpout);
  GADGET_WSKIP;

  //=============================================
  // POS
  //=============================================
  write_blockname(fpout, "POS ");
  
  blklen = 0;
  for(itype=0; itype<6; itype++) {
    blklen += gadget.header.np[itype] * (3*sizeof(float));
  }
  GADGET_WSKIP;
  
  ncheck_gas   = 0;
  ncheck_star  = 0;
  ncheck_dm    = 0;
  ncheck_dm_hr = 0;
  ncheck_dm_lr = 0;
  ncheck_all   = 0;
  for(itype=0; itype<4; itype++) {
    for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
      
      f[X] = (float) cur_part->pos[X]*x_fac;
      f[Y] = (float) cur_part->pos[Y]*x_fac;
      f[Z] = (float) cur_part->pos[Z]*x_fac;
      
      switch(itype) {
        case 0: // gas
          if (isgreaterequal(cur_part->u, PGAS)) {
            fwrite(f, sizeof(float), 3, fpout);
            ncheck_gas++;
            ncheck_all++;
          }
          break;
        case 1: // dm_hr
          if ( fabs(cur_part->weight - 1) < ZERO ) {
            fwrite(f, sizeof(float), 3, fpout);
            ncheck_dm_hr++;
            ncheck_all++;
          }
          break;
        case 2: // dm_lr
          if ( !(fabs(cur_part->u - PSTAR) < ZERO || isgreaterequal(cur_part->u, PGAS) || fabs(cur_part->weight - 1) < ZERO) ) {
            fwrite(f, sizeof(float), 3, fpout);
            ncheck_dm_lr++;
            ncheck_all++;
          }
          break;
        case 3: // star
          if (fabs(cur_part->u - PSTAR) < ZERO) {
            fwrite(f, sizeof(float), 3, fpout);
            ncheck_star++;
            ncheck_all++;
          }
          break;
      } // switch
    } // cur_part
  } // itype
  if(ncheck_star != nstar) {
    fprintf(stderr,"ERROR: nstar = %d, ncheck_star = %d\n",nstar,ncheck_star);
  }
  if(ncheck_gas != ngas) {
    fprintf(stderr,"ERROR: ngas = %d, ncheck_gas  = %d\n",ngas,ncheck_gas);
  }
  if(ncheck_dm_hr != ndm_hr) {
    fprintf(stderr,"ERROR: ndm_hr = %d, ncheck_dm   = %d\n",ndm_hr,ncheck_dm_hr);
  }
  if(ncheck_dm_lr != ndm_lr) {
    fprintf(stderr,"ERROR: ndm_lr = %d, ncheck_dm   = %d\n",ndm_lr,ncheck_dm_lr);
  }
  if(ncheck_all != global.no_part) {
    fprintf(stderr,"ERROR: nall = %d, ncheck_all   = %d\n",global.no_part,ncheck_all);
  }
  
  GADGET_WSKIP;
  
  
  //=============================================
  // VEL
  //=============================================
  write_blockname(fpout, "VEL ");
  
  blklen = 0;
  for(itype=0; itype<6; itype++) {
    blklen += gadget.header.np[itype] * (3*sizeof(float));
  }
  
  GADGET_WSKIP;
  
  for(itype=0; itype<4; itype++) {
    for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
      
      f[X] = (float) cur_part->mom[X]*v_fac;
      f[Y] = (float) cur_part->mom[Y]*v_fac;
      f[Z] = (float) cur_part->mom[Z]*v_fac;
      
      switch(itype) {
        case 0: // gas
          if (isgreaterequal(cur_part->u, PGAS)) {
            fwrite(f, sizeof(float), 3, fpout);
          }
          break;
        case 1: // dm_hr
          if ( fabs(cur_part->weight - 1) < ZERO ) {
            fwrite(f, sizeof(float), 3, fpout);
          }
          break;
        case 2: // ndm_lr
          if ( !(fabs(cur_part->u - PSTAR) < ZERO || isgreaterequal(cur_part->u, PGAS) || fabs(cur_part->weight - 1) < ZERO) ) {
            fwrite(f, sizeof(float), 3, fpout);
          }
          break;
        case 3: // star
          if (fabs(cur_part->u - PSTAR) < ZERO) {
            fwrite(f, sizeof(float), 3, fpout);
          }
          break;
      } // switch
    } // cur_part
  } // itype

  GADGET_WSKIP;
  
  //=============================================
  // IDS
  //=============================================
  write_blockname(fpout, "ID  ");
  
  blklen = 0;
  for(itype=0; itype<6; itype++) {
    blklen += gadget.header.np[itype] * sizeof(int);
  }
  
  GADGET_WSKIP;
  
  for(itype=0; itype<4; itype++) {
    for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
      
      i = (int) cur_part->id;
      
      switch(itype) {
        case 0: // gas
          if (isgreaterequal(cur_part->u, PGAS)) {
            fwrite(&i, sizeof(int), 1, fpout);
          }
          break;
        case 1: // dm_hr
          if ( fabs(cur_part->weight - 1) < ZERO ) {
            fwrite(&i, sizeof(int), 1, fpout);
          }
          break;
        case 2: // dm_lr
          if ( !(fabs(cur_part->u - PSTAR) < ZERO || isgreaterequal(cur_part->u, PGAS) || fabs(cur_part->weight - 1) < ZERO) ) {
            fwrite(&i, sizeof(int), 1, fpout);
          }
          break;
        case 3: // star
          if (fabs(cur_part->u - PSTAR) < ZERO) {
            fwrite(&i, sizeof(int), 1, fpout);
          }
          break;
      } // switch
    } // cur_part
  } // itype
  
  GADGET_WSKIP;
  
  //=============================================
  // MASS
  //=============================================
  write_blockname(fpout, "MASS");
  
  blklen = 0;
  for(itype=0; itype<6; itype++) {
    blklen += gadget.header.np[itype] * sizeof(float);
  }
  
  GADGET_WSKIP;
  
  for(itype=0; itype<4; itype++) {
    for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
      
      m = (float) cur_part->weight*m_fac;
      
      switch(itype) {
        case 0: // gas
          if (isgreaterequal(cur_part->u, PGAS)) {
            fwrite(&m, sizeof(float), 1, fpout);
          }
          break;
        case 1: // dm_hr
          if ( fabs(cur_part->weight - 1) < ZERO ) {
            fwrite(&m, sizeof(float), 1, fpout);
          }
          break;
        case 2: // dm_lr
          if ( !(fabs(cur_part->u - PSTAR) < ZERO || isgreaterequal(cur_part->u, PGAS) || fabs(cur_part->weight - 1) < ZERO) ) {
            fwrite(&m, sizeof(float), 1, fpout);
          }
          break;
        case 3: // star
          if (fabs(cur_part->u - PSTAR) < ZERO) {
            fwrite(&m, sizeof(float), 1, fpout);
          }
          break;
      } // switch
    } // cur_part
  } // itype

  GADGET_WSKIP;
  
  //=============================================
  // U
  //=============================================
  write_blockname(fpout, "U   ");
  
  blklen = gadget.header.np[0] * sizeof(float);
  GADGET_WSKIP;
  
  for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
    
    u = (float) cur_part->u;  // always stored in physical units
    if (isgreaterequal(cur_part->u, PGAS)) {
      fwrite(&u, sizeof(float), 1, fpout);
    }
  } // cur_part
  
  GADGET_WSKIP;
  
#ifdef STORE_MORE
  //=============================================
  // RHO
  //=============================================
  write_blockname(fpout, "RHO ");
  
  blklen = gadget.header.np[0] * sizeof(float);
  GADGET_WSKIP;
  
  for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
    
    rho = (float) cur_part->rho * rho_fac;
    if (isgreaterequal(cur_part->u, PGAS)) {
      fwrite(&rho, sizeof(float), 1, fpout);
    }
  } // cur_part
  
  GADGET_WSKIP;

  //=============================================
  // HSL
  //=============================================
  write_blockname(fpout, "HSL ");
  
  blklen = 0;
  for(itype=0; itype<6; itype++) {
    blklen += gadget.header.np[itype] * sizeof(float);
  }
  GADGET_WSKIP;
  
  for(itype=0; itype<4; itype++) {
    for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++) {
      
      eps = (float) cur_part->eps*x_fac;
      
      switch(itype) {
        case 0: // gas
          if (isgreaterequal(cur_part->u, PGAS)) {
            fwrite(&eps, sizeof(float), 1, fpout);
          }
          break;
        case 1: // dm_hr
          if ( fabs(cur_part->weight - 1) < ZERO ) {
            fwrite(&eps, sizeof(float), 1, fpout);
          }
          break;
        case 2: // dm_lr
          if ( !(fabs(cur_part->u - PSTAR) < ZERO || isgreaterequal(cur_part->u, PGAS) || fabs(cur_part->weight - 1) < ZERO) ) {
            fwrite(&eps, sizeof(float), 1, fpout);
          }
          break;
        case 3: // star
          if (fabs(cur_part->u - PSTAR) < ZERO) {
            fwrite(&eps, sizeof(float), 1, fpout);
          }
          break;
      } // switch
    } // cur_part
  } // itype
  
  GADGET_WSKIP;
#endif
  
  fclose(fpout);
}

void write_blockname(FILE *fpout, char *name)
{
  int nextblock;
  
  blklen = 4*sizeof(char) + sizeof(int);
  GADGET_WSKIP;
  fwrite(name, sizeof(char), 4, fpout);
  nextblock = 0; // TODO: as none of my reading routines is using this number, just set it to zero
  fwrite(&nextblock, sizeof(int), 1, fpout);
  GADGET_WSKIP;

}

/*=============================================================================
 * get_ngas()
 *============================================================================*/
int get_ngas(partptr fst_part, long npart)
{
  long    ipart;
  int     ngas;
  partptr cur_part;
  
  ngas = 0;
  for(ipart=0; ipart<npart; ipart++) {
    cur_part = fst_part+ipart;
    if (isgreaterequal(cur_part->u, PGAS)) {
      ngas++;
    }
  }
  return(ngas);
}

/*=============================================================================
 * get_nstar()
 *============================================================================*/
int get_nstar(partptr fst_part, long npart)
{
  long    ipart;
  int     nstar;
  partptr cur_part;
  
  nstar = 0;
  for(ipart=0; ipart<npart; ipart++) {
    cur_part = fst_part+ipart;
    if (fabs(cur_part->u - PSTAR) < ZERO) {
      nstar++;
    }
  }
  return(nstar);
}

/*=============================================================================
 * get_ndm_hr()
 *============================================================================*/
int get_ndm_hr(partptr fst_part, long npart)
{
  long    ipart;
  int     ndm_hr;
  partptr cur_part;
  
  ndm_hr = 0;
  for(ipart=0; ipart<npart; ipart++) {
    cur_part = fst_part+ipart;
    if (fabs(cur_part->weight - 1) < ZERO) {
      ndm_hr++;
    }
  }
  return(ndm_hr);
}

#endif // MULTIMASS && GAS_PARTICLES


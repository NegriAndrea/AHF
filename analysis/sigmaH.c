//===============================================================================
// sigmaH:
// =======
//    calculate the dispersion of line-of-sight Hubble parameter measurements
//
//
// - the code is a copy of AHF's main.c and hence uses libio/ to read in the
//   actual data file
//
// - for that reason the sigmaH.input is a copy of AHF.input and hence contains
//   more information that needed for sigmaH.c ... c'est la vie!
//
// - in the present version the observer will be placed at a random subset
//   of the positions found in the input data file
//
// - the number of these positions is one additional input parameter (command line)
//
//
//
// logic of operation:
// -------------------
// - the points read in from file are referred to as "particles"
// - the neighbours to each observer (i.e. random particle) will be collected
//   into a HALO structure using existing routines from AHF
//===============================================================================


//========================================================
//  DEFINES
//========================================================
// where to place the observer: these options are mutually exclusive and only one must be defined!
//#define OBSERVER_CLOSEST_RANDOM_HALO   // finds halo closest to a random point
#define OBSERVER_RANDOM_HALO           // uses a random halo
//#define OBSERVER_RANDOM_POINT          // uses a random point (putting observer at rest)
//#define OBSERVER_FROM_FILE               // reads observer positions (and velocity) from file

// this will force the observer's velocity to be zero
//#define OBSERVER_AT_REST

//#define DEBUG
#define ISEED 123456         // seed for ran3() when picking random observers


#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "../src/common.h"

#include "../src/libsfc/sfc.h"
#include "../src/startrun.h"
#include "../src/libutility/utility.h"
#include "../src/libahf/ahf.h"
#include "../src/libahf/ahf_halos_sfc.h"
#include "../src/libahf/ahf_io.h"

#ifdef WITH_OPENMP
#include <omp.h>
#endif


//========================================================
//  prototypes
//========================================================
#define Hidx(i,j,L) ((i) + (j)*(L))


/*==============================================================================
 * MAIN: where everything starts ....
 *==============================================================================*/
int main(int argc, char **argv)
{
  FILE    *fpout;
  char     AMIGA_input[MAXSTRING], outfile[MAXSTRING];
  int      no_timestep;
  double   timecounter;
  double   timestep;
  uint64_t Nobserver;
  uint64_t iobserver, jobserver;
  partptr  cur_part;
  double   Rsphere,Rsphere_min,Rsphere_max,Rfrac;
  int      Nspheres;
  uint64_t isphere;
  int      iseed=ISEED;
  HALO     halo;
  double   x_fac, r_fac, v_fac;
  uint64_t ineighbour;
  double   dx,dy,dz,d;
  double   vx,vy,vz,vlos;
  double   Hubble;
  double  *Hlos, *Hsphere, *sigmaHsphere;
  uint64_t norm;
  int      omp_id;
  time_t   elapsed = (time_t)0;
  
#ifdef OBSERVER_FROM_FILE
  FILE    *fpin;
  char     observer_file[MAXSTRING], line[MAXSTRING];
  double  *xobs, *yobs, *zobs, *vxobs, *vyobs, *vzobs;
  int      idummy;
#endif

  //========================================================
  //               deal with command line
  //========================================================
#ifdef OBSERVER_FROM_FILE
  if(argc<6) {
    fprintf(stderr,"usage:    %s sigmaH.input observer_file Nspheres Rsphere_min Rsphere_max\n", argv[0]);
    fprintf(stderr,"       or %s --parameterfile\n", argv[0]);
    exit(1);
  }
#else
  if(argc<6) {
    fprintf(stderr,"usage:    %s sigmaH.input Nobserver Nspheres Rsphere_min Rsphere_max\n", argv[0]);
    fprintf(stderr,"       or %s --parameterfile\n", argv[0]);
    exit(1);
  }
#endif
  
  // maybe the user only wants the parameterfile?
  if(strcmp(argv[1],"--parameterfile") == 0)  {
    global_io.params                 = (io_parameter_t) calloc(1,sizeof(io_parameter_struct_t));
    global_io.params->outfile_prefix = (char *) calloc(MAXSTRING,sizeof(char));
    global.a                         = 1;
    strcpy(global_io.params->outfile_prefix,"AHF");
    write_parameterfile();
    exit(0);
  }
  else {
    strcpy(AMIGA_input, argv[1]);
    Nspheres     = (uint64_t) atoi(argv[3]);
    Rsphere_min  = (double) atof(argv[4]);
    Rsphere_max  = (double) atof(argv[5]);

#ifdef OBSERVER_FROM_FILE
    strcpy(observer_file, argv[2]);
#else
    Nobserver    = (uint64_t) atoi(argv[2]);
#endif
  }
  
#ifdef OBSERVER_FROM_FILE
  //========================================================
  //               get observers from file
  //========================================================
  fpin = fopen(observer_file,"r");
  if(fpin == NULL) {    
    fprintf(stderr,"FATAL: cannot open %s for reading\n",observer_file);
    exit(0);
  }
  
  Nobserver = 0;
  xobs      = NULL;
  yobs      = NULL;
  zobs      = NULL;
  vxobs     = NULL;
  vyobs     = NULL;
  vzobs     = NULL;
  // read first line
  fgets(line,MAXSTRING,fpin);
  
  while(!feof(fpin)) {
    if(strncmp(line,"#",1) != 0) {
      Nobserver++;
      xobs  = (double *) realloc(xobs, Nobserver*sizeof(double));
      yobs  = (double *) realloc(yobs, Nobserver*sizeof(double));
      zobs  = (double *) realloc(zobs, Nobserver*sizeof(double));
      vxobs = (double *) realloc(vxobs,Nobserver*sizeof(double));
      vyobs = (double *) realloc(vyobs,Nobserver*sizeof(double));
      vzobs = (double *) realloc(vzobs,Nobserver*sizeof(double));
      
      // scan line for positions and velocities (if there are no velocities, do not try to use them!)
      
      // generic input format
//      sscanf(line,"%lf %lf %lf %lf %lf %lf",
//             (xobs+(Nobserver-1)),
//             (yobs+(Nobserver-1)),
//             (zobs+(Nobserver-1)),
//             (vxobs+(Nobserver-1)),
//             (vyobs+(Nobserver-1)),
//             (vzobs+(Nobserver-1))  );
      
      // void files from Stefan
      sscanf(line,"%d %lf %lf %lf",
             &idummy,
             (xobs+(Nobserver-1)),
             (yobs+(Nobserver-1)),
             (zobs+(Nobserver-1)));
      vxobs[Nobserver-1] = 0.0;
      vyobs[Nobserver-1] = 0.0;
      vzobs[Nobserver-1] = 0.0;
    }
    
    // read next line
    fgets(line,MAXSTRING,fpin);
  }
  
  fclose(fpin);  
#endif
  
  //========================================================
  //                     startrun
  //========================================================
  elapsed  = (time_t)0;
  elapsed -= time(NULL);

  timing.io       -= time(NULL);
  timing.startrun -= time(NULL);
  
	startrun((argc > 1) ? argv[1] : NULL, &timecounter, &timestep, &no_timestep);
  
  timing.startrun += time(NULL);
  timing.io       += time(NULL);

  fprintf(stderr,"\nsigmaH parameters:\n");
  fprintf(stderr,"------------------\n");
  fprintf(stderr,"Nobserver     = %"PRIu64"\n",Nobserver);
  fprintf(stderr,"Nspheres      = %d\n",Nspheres);
  fprintf(stderr,"Rsphere_min   = %lf\n",Rsphere_min);
  fprintf(stderr,"Rsphere_max   = %lf\n\n",Rsphere_max);
  
  
  elapsed += time(NULL);
  fprintf(stderr,"o startrun done in %ld sec.\n",elapsed);

  // some relevant stuff
  global.fst_part = global_info.fst_part;
  
  //========================================================
	//                 generate SFC keys
  //========================================================
  elapsed  = (time_t)0;
  elapsed -= time(NULL);
  fprintf(stderr,"o generating sfc keys ... ");

  timing.sfckey -= time(NULL);
	for (uint64_t i=0; i<global_info.no_part; i++) {
		partptr part=global_info.fst_part+i;
		part->sfckey = sfc_curve_calcKey(global_info.ctype,
		                                 (double)(part->pos[0]),
		                                 (double)(part->pos[1]),
		                                 (double)(part->pos[2]),
		                                 BITS_PER_DIMENSION);
	}
	// sorting all particles to have fast access later on
	qsort(global_info.fst_part,
	      global_info.no_part,
	      sizeof(part),
	      &cmp_sfckey_part);
  timing.sfckey += time(NULL);

  elapsed += time(NULL);
  fprintf(stderr,"   done in %ld sec.\n",elapsed);
 
  
  
  
  
  
  
  
  //========================================================
	//                    analysis
  //========================================================
  elapsed  = (time_t)0;
  elapsed -= time(NULL);
  fprintf(stderr,"o placing observers:\n");
  
  ahf.time -= time(NULL);
    
  // derive mean Hlos for every observer and every sphere radius
  //=============================================================
  
  // dump observer positions to file
  sprintf(outfile,"%s-Nobserver%07"PRIu64,global_io.params->outfile_prefix,Nobserver);
  fpout = fopen(outfile,"w");
  if(fpout == NULL) {
    fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
    exit(0);
  }

  
  // prepare some conversion factors
  x_fac   = simu.boxsize;
  r_fac   = simu.boxsize * global.a;
  v_fac   = simu.boxsize / simu.t_unit / global.a;   //NOTE: AHF stores a^2 \dot{x} as velocity!!!!
	Hubble  = calc_Hubble(global.a);                   // [km/sec/Mpc]

  fprintf(stderr,"    Hubble  = %lf\n",Hubble);
  fprintf(stderr,"    boxsize = %lf\n",simu.boxsize);
  fprintf(stderr,"    a       = %lf\n",global.a);
  fprintf(stderr,"    r_fac   = %lf\n",r_fac);
  fprintf(stderr,"    v_fac   = %lf\n",v_fac);
  
  
  // array to hold all Hlos values
  Hlos = (double *) calloc(Nobserver*Nspheres, sizeof(double));
  
#ifdef WITH_OPENMP
  // obtain thread number to serve as seed
  omp_id = omp_get_thread_num();
  iseed *= omp_id;
#endif

  // loop over all observers
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic) private(jobserver,iobserver,cur_part,halo,isphere,Rsphere,ineighbour,dx,dy,dz,d,vx,vy,vz,vlos,iseed,norm) shared(Hlos,Nobserver,global_info,v_fac,r_fac)
#endif
  for(iobserver=0; iobserver<Nobserver; iobserver++) {
    
    
    
    /*===========================================================
     *                  OBSERVER_FROM_FILE
     *===========================================================*/
#ifdef OBSERVER_FROM_FILE
    halo.pos.x = xobs[iobserver]   /x_fac;
    halo.pos.y = yobs[iobserver]   /x_fac;
    halo.pos.z = zobs[iobserver]   /x_fac;
    halo.vel.x = vxobs[iobserver]  /v_fac;
    halo.vel.y = vyobs[iobserver]  /v_fac;
    halo.vel.z = vzobs[iobserver]  /v_fac;
#endif // OBSERVER_FROM_FILE
    
    
    /*===========================================================
     *              OBSERVER_CLOSEST_RANDOM_HALO
     *===========================================================*/
#ifdef OBSERVER_CLOSEST_RANDOM_HALO
    
    // pick random position throughout box
    halo.pos.x     = (float)(ran3(&iseed));
    halo.pos.y     = (float)(ran3(&iseed));
    halo.pos.z     = (float)(ran3(&iseed));
    halo.npart     = 0;
    Rfrac          = 0.25;
    
    while(halo.npart == 0) {
#ifdef DEBUG
      fprintf(stderr,"   observer at %f %f %f\n",halo.pos.x*simu.boxsize,halo.pos.y*simu.boxsize,halo.pos.z*simu.boxsize);
#endif
      
      // search for neighbours about that position
      halo.gatherRad = Rfrac*Rsphere_min/simu.boxsize;
      halo.npart     = 0;
      halo.ipart     = NULL;
      ahf_halos_sfc_gatherParts(&halo);
#ifdef DEBUG
      fprintf(stderr,"    -> found %ld nearest neighbours to observer\n",halo.npart);
#endif
      
      // pick closest particle
      if(halo.npart > 0) {
        sort_halo_particles(&halo);
        cur_part = global_info.fst_part + halo.ipart[0];
#ifdef DEBUG
        fprintf(stderr,"      -> using ipart=%ld\n",halo.ipart[0]);
#endif
      }
      else {
        Rfrac *= 2.;
#ifdef DEBUG
        fprintf(stderr,"      -> increasing neighbour search radius to %lf\n",Rfrac*Rsphere_min);
#endif
      }
    }

    // remove from halo structure again
    if(halo.npart > 0) free(halo.ipart);
    halo.ipart = NULL;
    halo.npart = 0;
    
    // prepare halo structure as we are re-using ahf_halos_sfc_gatherParts()
    halo.pos.x     = cur_part->pos[X];
    halo.pos.y     = cur_part->pos[Y];
    halo.pos.z     = cur_part->pos[Z];
    halo.vel.x     = cur_part->mom[X];
    halo.vel.y     = cur_part->mom[Y];
    halo.vel.z     = cur_part->mom[Z];
#endif // OBSERVER_CLOSEST_RANDOM_HALO
    
    

    /*===========================================================
     *                  OBSERVER_RANDOM_HALO
     *===========================================================*/
#ifdef OBSERVER_RANDOM_HALO
    // pick a random particle as the observer
    jobserver = (uint64_t) (ran3(&iseed) * global_info.no_part);
    while(jobserver >= global_info.no_part)
      jobserver = (uint64_t) (ran3(&iseed) * global_info.no_part);
    cur_part  = global_info.fst_part + jobserver;
    
    // prepare halo structure as we are re-using ahf_halos_sfc_gatherParts()
    halo.pos.x     = cur_part->pos[X];
    halo.pos.y     = cur_part->pos[Y];
    halo.pos.z     = cur_part->pos[Z];
    //fprintf(stderr,"cur_part=%ld jobserver=%ld global_info.no_part=%ld global_info.fst_part=%ld...",cur_part,jobserver,global_info.no_part,global_info.fst_part);
    halo.vel.x     = cur_part->mom[X];
    halo.vel.y     = cur_part->mom[Y];
    halo.vel.z     = cur_part->mom[Z];
    //fprintf(stderr,"set\n");
#endif
   

    /*===========================================================
     *                  OBSERVER_RANDOM_POINT
     *===========================================================*/
#ifdef OBSERVER_RANDOM_POINT
    halo.pos.x     = (float)(ran3(&iseed));
    halo.pos.y     = (float)(ran3(&iseed));
    halo.pos.z     = (float)(ran3(&iseed));
    halo.vel.x     = 0.0;
    halo.vel.y     = 0.0;
    halo.vel.z     = 0.0;
#endif
    
    
    /*===========================================================
     *                    OBSERVER_AT_REST
     *===========================================================*/
#ifdef OBSERVER_AT_REST
    halo.vel.x     = 0.0;
    halo.vel.y     = 0.0;
    halo.vel.z     = 0.0;
#endif
    
    fprintf(fpout,"%f %f %f %f %f %f\n",
            halo.pos.x*simu.boxsize,
            halo.pos.y*simu.boxsize,
            halo.pos.z*simu.boxsize,
            halo.vel.x*v_fac,
            halo.vel.y*v_fac,
            halo.vel.z*v_fac);
    fflush(fpout);

    
    
    
    
    
    
    
    
    
    
    
    
    /*===========================================================
     *                      SPHERE LOOP
     *===========================================================*/
    // loop over all spheres for this observer
    for(isphere=0; isphere<Nspheres; isphere++) {
      
      // sphere radius in internal units
      Rsphere = (Rsphere_min + (double)isphere*(Rsphere_max-Rsphere_min)/((double)(Nspheres-1)))/simu.boxsize;
      
      // add gathering radius and prepare array to hold neighbours
      halo.gatherRad = Rsphere;
      halo.npart     = 0;
      halo.ipart     = NULL;
      
      // collect all particles inside the gathering radius
      ahf_halos_sfc_gatherParts(&halo);


#ifdef DEBUG
      fprintf(stderr," (sfc gathered %ld neighbours)",halo.npart);
#endif

      // perform statistic for all neighbours
      Hlos[Hidx(isphere,iobserver,Nspheres)] = 0.0;
      norm                                   = 0;
      
      for(ineighbour=0; ineighbour<halo.npart; ineighbour++) {
        
        cur_part = global_info.fst_part+halo.ipart[ineighbour];
        
        // 3D distance
        dx  = (cur_part->pos[X] - halo.pos.x);
        dy  = (cur_part->pos[Y] - halo.pos.y);  // distance in internal units
        dz  = (cur_part->pos[Z] - halo.pos.z);
        
        if (dx >  0.5) dx -= 1.0;
        if (dy >  0.5) dy -= 1.0;
        if (dz >  0.5) dz -= 1.0;
        if (dx < -0.5) dx += 1.0;   // take care of periodic boundary conditions
        if (dy < -0.5) dy += 1.0;
        if (dz < -0.5) dz += 1.0;
        
        dx  *= r_fac;
        dy  *= r_fac;  // distance (eventually) in physical units
        dz  *= r_fac;

        d   = sqrt(pow2(dx)+pow2(dy)+pow2(dz));
        
        // 3D velocity
        vx  = (cur_part->mom[X] - halo.vel.x) * v_fac;
        vy  = (cur_part->mom[Y] - halo.vel.y) * v_fac; // peculiar velocity in physical units
        vz  = (cur_part->mom[Z] - halo.vel.z) * v_fac;
        vx += Hubble * (dx);
        vy += Hubble * (dy); // add correct Hubble flow
        vz += Hubble * (dz);

        // line-of-sight velocity as measured by present observer
        vlos = (vx*dx + vy*dy + vz*dz)/d;
        
        // accumulate Hlos
        if(d > MACHINE_ZERO) {
          // Hubble parameter as measured by present observer
          Hlos[Hidx(isphere,iobserver,Nspheres)] += vlos/d;
          norm++;
        }
      } // ineighbour
      
      // mean Hlos for this sphere
      if(norm > 0)
        Hlos[Hidx(isphere,iobserver,Nspheres)] /= (double)norm;
      else
        Hlos[Hidx(isphere,iobserver,Nspheres)] = -1.;
      
      // physically remove particles from memory
      if(halo.ipart) {
        free(halo.ipart);
        halo.npart = 0;
        halo.ipart = NULL;
      }
    } // isphere
   
  } // iobserver
  
  fclose(fpout);
  elapsed += time(NULL);
  fprintf(stderr,"   done in %ld sec.\n",elapsed);
 
  
  // collapse Hlos values obtaining mean and stddev for every sphere
  //=================================================================
  elapsed  = (time_t)0;
  elapsed -= time(NULL);
  fprintf(stderr,"o calculating means and stddevs ... ");

  Hsphere      = (double *) calloc(Nspheres, sizeof(double));
  sigmaHsphere = (double *) calloc(Nspheres, sizeof(double));
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static) private(isphere,iobserver) shared(Nobserver,Hlos,Hsphere)
#endif
  for(isphere=0; isphere<Nspheres; isphere++) {
    norm = 0;
    for(iobserver=0; iobserver<Nobserver; iobserver++) {
      if(Hlos[Hidx(isphere,iobserver,Nspheres)] > 0.) {
        Hsphere[isphere] += Hlos[Hidx(isphere,iobserver,Nspheres)];
        norm++;
      }
    }
    if(norm > 0)
      Hsphere[isphere] /= (double)norm;
    else
      Hsphere[isphere] = -1.;
  }
  
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static) private(isphere,iobserver) shared(Nobserver,Hlos,Hsphere,Hubble)
#endif
  for(isphere=0; isphere<Nspheres; isphere++) {
    norm = 0;
    for(iobserver=0; iobserver<Nobserver; iobserver++) {
      if(Hsphere[isphere] > 0.) {
        sigmaHsphere[isphere] += pow2((Hsphere[isphere]-Hubble)/Hubble);
        norm++;
      }
    }
    if(norm > 0)
      sigmaHsphere[isphere] /= (double)norm;
    else
      sigmaHsphere[isphere] -1.;
  }

  elapsed += time(NULL);
  fprintf(stderr,"   done in %ld sec.\n",elapsed);

  ahf.time += time(NULL);
  
  
  
  
 
  
  
  //========================================================
	//                       output
  //========================================================
  elapsed  = (time_t)0;
  elapsed -= time(NULL);
  
  // full information for each sphere
  //==================================
  for(isphere=0; isphere<Nspheres; isphere++) {
    
    // sphere radius
    Rsphere = (Rsphere_min + (double)isphere*(Rsphere_max-Rsphere_min)/((double)(Nspheres-1)));

    // construct outfile name
    sprintf(outfile,"%s-Nobserver%07"PRIu64"-Nspheres%03d-Rpshere%lf",global_io.params->outfile_prefix,Nobserver,Nspheres,Rsphere);
    fprintf(stderr,"o writing Rsphere=%lf information to: %s ... ",Rsphere,outfile);
    
    // open output file
    fpout = fopen(outfile,"w");
    if(fpout == NULL) {
      fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
      exit(0);
    }

    fprintf(fpout,"# Hlos(1)\n");

    for(iobserver=0; iobserver<Nobserver; iobserver++) {
      fprintf(fpout,"%lf\n",Hlos[Hidx(isphere,iobserver,Nspheres)]);
    }
  }
  fclose(fpout);

  
  // reduced information
  //=====================
  sprintf(outfile,"%s-Nobserver%07"PRIu64"-Nspheres%03d-Rsphere_min%lf-Rpshere_max%lf",global_io.params->outfile_prefix,Nobserver,Nspheres,Rsphere_min,Rsphere_max);
  fprintf(stderr,"o writing sphere information to: %s ... ",outfile);
  
  // open output file
  fpout = fopen(outfile,"w");
  if(fpout == NULL) {
    fprintf(stderr,"FATAL: cannot open %s for writing\n",outfile);
    exit(0);
  }

  fprintf(fpout,"# Rsphere(1) Hsphere(2) sigmaH(3)\n");
  for(isphere=0; isphere<Nspheres; isphere++) {
    Rsphere = (Rsphere_min + (double)isphere*(Rsphere_max-Rsphere_min)/((double)(Nspheres-1)));
    fprintf(fpout,"%lf %lf %lf\n",Rsphere,Hsphere[isphere],sigmaHsphere[isphere]);
  }
  fclose(fpout);
  
  elapsed += time(NULL);
  fprintf(stderr," done in %ld sec.\n",elapsed);
  
  
  
  //========================================================
	//                 update logfile
  //========================================================
  //write_logfile(timecounter, timestep, no_timestep);

  //========================================================
	//                      BYE BYE
  //========================================================
  free(Hlos);
  free(Hsphere);
  free(sigmaHsphere);
#ifdef OBSERVER_FROM_FILE
  free(xobs);
  free(yobs);
  free(zobs);
  free(vxobs);
  free(vyobs);
  free(vzobs);
#endif
  
  if(io.icfile_name)       free(io.icfile_name);
  if(io.dumpfile_name)     free(io.dumpfile_name);
  if(io.logfile_name)      free(io.logfile_name);
  if(io.outfile_prefix)    free(io.outfile_prefix);
  if(global.termfile_name) free(global.termfile_name);
  
  if(global.fst_part)      free(global.fst_part);
  if(global.fst_gas)       free(global.fst_gas);
  if(global.fst_star)      free(global.fst_star);
  
  fprintf(io.logfile, "==========================================================\n");
  fprintf(io.logfile, "                       FINISHED (v%3.1f/%03d)\n",VERSION,BUILD);
  fprintf(io.logfile, "==========================================================\n");
  fclose(io.logfile);
  
  
  return EXIT_SUCCESS;
}

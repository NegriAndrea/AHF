#include "include.h"
#include "common.h"


/*==================================================================================================
 * write_mtree:
 *==================================================================================================*/
int write_mtree(char OutFile[MAXSTRING])
{
  uint64_t  ihalo, nHalos0_good;
  int64_t   icroco;
  FILE *fpout, *fpout_idx, *fpout_croco;
  char outname[MAXSTRING], outname_idx[MAXSTRING], outname_croco[MAXSTRING];
  clock_t   elapsed;
  
  elapsed = clock();
  
  sprintf(outname,"%s_mtree",OutFile);
  sprintf(outname_croco,"%s_croco",OutFile);
  strcpy(outname_idx, outname);
  strcat(outname_idx, "_idx");
  
  fpout = fopen(outname,"w");
  if(fpout == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname);
    exit(0);
  }
  
  fpout_idx = fopen(outname_idx,"w");
  if(fpout_idx == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname_idx);
    exit(0);
  }
  
  fpout_croco = fopen(outname_croco,"w");
  if(fpout_croco == NULL) {
    fprintf(stderr,"could not open file %s\nexiting\n",outname_croco);
    exit(0);
  }
  
  // count the number of halos with a progenitor
  nHalos0_good = (uint64_t)0;
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    if(halos[0][ihalo].ncroco > 0) {
      nHalos0_good++;
    }
  }
  fprintf(stderr,"  o writing cross-correlation for %"PRIu64" haloes (out of %"PRIu64" in total) ...",nHalos0_good,nHalos[0]);
  fprintf(fpout,"%"PRIu64"\n",nHalos0_good); // we do not write the version number into the _mtree file as they are merged into one single SUSSING2013 file by the post-processing script
  fprintf(fpout_croco,"#MergerTree version: %f\n",mtree_version);
  fprintf(fpout_croco,"#   HaloID(1)   HaloPart(2)  NumProgenitors(3)\n");
  fprintf(fpout_croco,"#      SharedPart(1)    HaloID(2)   HaloPart(3) merit(4)\n");
  fprintf(fpout_idx,"# HaloID(1) HaloID(2), MergerTree version: %f\n",mtree_version);
  fflush(fpout);
  fflush(fpout_idx);
  fflush(fpout_croco);
  
  for(ihalo=0; ihalo<nHalos[0]; ihalo++) {
    
    if(halos[0][ihalo].ncroco > 0) {
      //      this is the old format where the haloid corresponds to the linenumber
      //      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].mtree[0].id[0], halos[0][ihalo].mtree[0].id[1]);
      
      // this is the case where we use the haloid as found in *_particles
      fprintf(fpout_idx,"%12"PRIu64" %12"PRIu64"\n", halos[0][ihalo].haloid, halos[0][ihalo].mtree[0].haloid[1]);
      fflush(fpout_idx);
      
      fprintf(fpout,"%"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].ncroco);
      
      fprintf(fpout_croco,"%"PRIu64"  %"PRIu64"  %"PRIu64"\n",
              halos[0][ihalo].haloid,
              halos[0][ihalo].npart,
              halos[0][ihalo].ncroco);
      fflush(fpout);
      fflush(fpout_croco);
      
      for(icroco=0; icroco<halos[0][ihalo].ncroco; icroco++) {
        fprintf(fpout,"%"PRIu64"\n",
                halos[0][ihalo].mtree[icroco].haloid[1]);
        
        fprintf(fpout_croco,"  %"PRIu64"  %"PRIu64"  %"PRIu64" %lf\n",
                halos[0][ihalo].mtree[icroco].common,
                halos[0][ihalo].mtree[icroco].haloid[1],
                halos[0][ihalo].mtree[icroco].npart[1],
                halos[0][ihalo].mtree[icroco].merit);
        fflush(fpout);
        fflush(fpout_croco);
      }
    } // ncroco > 0
  } // for(ihalo)
  
  /* close files */
  fclose(fpout);
  fclose(fpout_idx);
  fclose(fpout_croco);
  
  elapsed = clock()-elapsed;
  fprintf(stderr," done in %4.2f sec.\n",(float)elapsed/CLOCKS_PER_SEC);
  return(1);
}

/*==================================================================================================
 * read_particles_bin:
 *
 * read AHF_particles_bin files (only existent for CLUES-WMAP3 simulation!?)
 *
 *==================================================================================================*/
int read_particles_bin(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin, *fphids;
  char      line[MAXSTRING], infile[MAXSTRING], hidsname[MAXSTRING];
  int64_t   ihalo, jhalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, jpart, Pid, Ptype, haloid, haloid_from_file;
  uint64_t  PidMin_local=((uint64_t)1<<62);
  uint64_t  PidMax_local=0;
  clock_t   elapsed;
  
  char c;
  int  i32, nh, ih;
  long i64;
  
  elapsed = clock();
  
  fprintf(stderr,"  o reading file %s ...",filename);
  
  fpin = fopen(filename,"rb");
  if(fpin == NULL)
  {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
  }
  
  /* reset all variables */
  nHalos[isimu] = 0;
  halos[isimu]  = NULL;
  
  // rubbish ?
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadInt(fpin,&i32,SWAPBYTES);
  
  // nuber of haloes
  ReadInt(fpin,&i32,SWAPBYTES);
  nh = i32;
  
  // rubbish ?
  fread(&c,sizeof(char),1,fpin);
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadInt(fpin,&i32,SWAPBYTES);
  ReadLong(fpin,&i64,SWAPBYTES);
  
  // loop over halos and their particle ids
  for(ihalo=0; ihalo<nh; ihalo++) {
    
    // use halo counter as haloid
    haloid = ihalo;
    
    // nparticles in halo
    ReadLong(fpin,&i64,SWAPBYTES);
    nPartInHalo = i64;
    
    /* found yet another halo */
    nHalos[isimu] += 1;
    halos[isimu]   = (HALOptr) realloc(halos[isimu], (nHalos[isimu]+1)*sizeof(HALOS));
    
    /* store haloid */
    halos[isimu][ihalo].haloid = haloid;
    
    /* halos[][].Pid will be incrementally filled using realloc() */
    halos[isimu][ihalo].Pid   = NULL;
    halos[isimu][ihalo].mtree = NULL;
    
    /* read all their id's */
    nPartInUse = 0;
    for(ipart=0; ipart<nPartInHalo; ipart++)
    {
      // read particle id
      ReadInt(fpin,&i32,SWAPBYTES);
      Pid = i32;
      
      // simulation only contains DM particles
      Ptype = 1;
      
      // here we can restrict the cross-correlation to a ceratain sub-set of all particles
      if(check_Ptype(Ptype) == TRUE)
      {
        halos[isimu][ihalo].Pid             = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
        if(halos[isimu][ihalo].Pid == NULL) {
          fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
          exit(-1);
        }
        halos[isimu][ihalo].Pid[nPartInUse] = Pid;
        
        if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] = (halos[isimu][ihalo].Pid[nPartInUse]);
        if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin        = (halos[isimu][ihalo].Pid[nPartInUse]);
        if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
        if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
        
        nPartInUse++;
      }
    }
    
    /* store number of particles in halo */
    halos[isimu][ihalo].npart = nPartInUse;
  } // for()
  
  fclose(fpin);
  
  elapsed = clock()-elapsed;
  
  fprintf(stderr," done in %4.2f sec. (nhalos = %"PRIu64", full ID range = %"PRIu64" -> %"PRIu64", local ID range = %"PRIu64" -> %"PRIu64")\n",
          (float)elapsed/CLOCKS_PER_SEC,nHalos[isimu],PidMin,PidMax[isimu],PidMin_local,PidMax_local);
  
  
  return(1);
}

/*==================================================================================================
 * read_particles:
 *
 * read the file storing the particle IDs for each halo
 *
 *      nHalos = number of halos found in file
 *      Pid    = id's of all those particles
 *
 *==================================================================================================*/
int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE     *fpin, *fphids;
  char      line[MAXSTRING], infile[MAXSTRING], hidsname[MAXSTRING];
  int64_t   ihalo, jhalo;
  uint64_t  nPartInHalo, nPartInUse, ipart, jpart, Pid, Ptype, haloid, haloid_from_file;
  uint64_t  PidMin_local=((uint64_t)1<<62);
  uint64_t  PidMax_local=0;
  clock_t   elapsed;
  
  elapsed = clock();
  
#ifdef READ_MPARTICLES
  int32_t  nfiles, ifile;
  uint64_t nHalosInFile;
  
  fprintf(stderr,"  o reading multiple files %s ",filename);
  
  // count the number of particles files/snapshot
  nfiles = count_particles_files(filename);
  fprintf(stderr," (nfiles=%"PRIi32,nfiles);
  if(nfiles == 0) {
    fprintf(stderr,"Could not open multiple _particles file complying with the required filename convention\nABORTING\n");
    exit(0);
  }
  
  // accumulate the number of halos by summing the first line in all _particles files
  nHalos[isimu] = count_halos(filename, nfiles);
  
  fprintf(stderr," nhalos=%"PRIu64") file: ",nHalos[isimu]);
  
  // allocate memory for haloes as one block
  halos[isimu]  = (HALOptr) calloc(nHalos[isimu], sizeof(HALOS));
  if(halos[isimu] == NULL) {
    fprintf(stderr,"\nCould not allocate memory for halos[] array of isimu = %d (size = %f GB)\nAborting\n",
            isimu,(float)(nHalos[isimu]*sizeof(HALOS)/1024./1024./1024.));
    exit(0);
  }
  
  // reset halo and total particle counter
  ihalo = 0;
  for(ifile=0; ifile<nfiles; ifile++) {
    
    fprintf(stderr,"%"PRIi32,ifile);
    
    // open file
    construct_filename(filename, ifile, infile);
    fpin = fopen(infile,"r");
    
    // read total number of haloes in file
    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%"SCNu64, &nHalosInFile);
    
    fprintf(stderr," (%"PRIu64") ",nHalosInFile);
    
    // loop over all haloes in this file
    for(jhalo=0; jhalo<nHalosInFile; jhalo++) {
      
      // get actual information
      fgets(line,MAXSTRING,fpin);
      sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid);
      
      // store haloid
      halos[isimu][ihalo].haloid = haloid;
      
      // loop over all particles in halo
      nPartInUse              = 0;
      halos[isimu][ihalo].Pid = NULL;
      for(ipart=0; ipart<nPartInHalo; ipart++)
      {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
        }
        else if(Ptype < 0 || Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
        }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
        if(check_Ptype(Ptype) == TRUE)
        {
          halos[isimu][ihalo].Pid = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+2)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"\n read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64" particles (Pid=%ld)\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1),halos[isimu][ihalo].Pid);
            exit(-1);
          }
          
          halos[isimu][ihalo].Pid[nPartInUse] = Pid;
          
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin        = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          
          nPartInUse++;
        } // if(Ptype)
        
      } // for(ipart)
      
      // store number of particles in halo
      halos[isimu][ihalo].npart = nPartInUse;
      
      // move to next halo in overall list
      ihalo++;
      
    } // for(jhalo)
    
    // close file
    fclose(fpin);
    
  } // for(ifile)
  
#else // READ_MPARTICLES
  
  fprintf(stderr,"  o reading file %s ...",filename);
  
  fpin = fopen(filename,"r");
  if(fpin == NULL)
  {
    fprintf(stderr,"could not open file %s\nexiting!\n",filename);
    exit(0);
  }
  
#ifdef READ_HALOIDS_FROM_FILE
  sprintf(hidsname,"%s_hids",filename);
  fphids = fopen(hidsname,"r");
  if(fphids == NULL)
  {
    fprintf(stderr,"could not open file %s\nexiting!\n",hidsname);
    exit(0);
  }
#endif
  
  /* reset all variables */
  nHalos[isimu] = 0;
  ihalo         = -1;
  halos[isimu]  = NULL;
  
  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
  
#ifndef THERE_IS_NO_NHALOS_LINE
  /* for AHF_particles files the first line is numGoodHalos which we will happily ignore and count that number ourselves */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
    fgets(line,MAXSTRING,fpin);
#endif
  
  do {
    if(strncmp(line,"#",1) != 0)
    {
      /* has a haloid been written */
      if(sscanf(line,"%"SCNu64" %"SCNu64, &nPartInHalo, &haloid) == 1)
      {
        /* if not, just get the number of particles */
        sscanf(line,"%"SCNu64, &nPartInHalo);
        
        /* and use halo counter as id */
        haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
      }
#ifdef USE_LINENUMBER_AS_HALOID
      haloid = ihalo+1; // +1, because ihalo has not been incremented yet!
#endif
#ifdef READ_HALOIDS_FROM_FILE
      fscanf(fphids,"%"SCNi64,&haloid_from_file);
      haloid = haloid_from_file;
#endif
      
      /* found yet another halo */
      ihalo++;
      nHalos[isimu] += 1;
      halos[isimu]   = (HALOptr) realloc(halos[isimu], (nHalos[isimu]+1)*sizeof(HALOS));
      
      /* store haloid */
      halos[isimu][ihalo].haloid = haloid;
      
      /* halos[][].Pid will be incrementally filled using realloc() */
      halos[isimu][ihalo].Pid   = NULL;
      halos[isimu][ihalo].mtree = NULL;
      
      /* read all their id's */
      nPartInUse = 0;
      for(ipart=0; ipart<nPartInHalo; ipart++)
      {
        /* read line containing ID and possibly some more information */
        fgets(line,MAXSTRING,fpin);
        
        /* check whether we are able to read a meaningful particle type, too */
        if(sscanf(line,"%"SCNu64" %"SCNu64, &Pid, &Ptype) == 1) {
          /* if not, set Ptype to 1 as this is the type we will use below */
          Ptype = 1;
        }
        else if(Ptype < 0 || Ptype > abs(PDMbndry)) {
          /* not a meaningful type, maybe something else has been stored? */
          Ptype = 1;
        }
        
        // here we can restrict the cross-correlation to a ceratain sub-set of all particles
        if(check_Ptype(Ptype) == TRUE)
        {
          halos[isimu][ihalo].Pid = (uint64_t *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(uint64_t));
          if(halos[isimu][ihalo].Pid == NULL) {
            fprintf(stderr,"read_particles: could not realloc() halos[%d][%ld].Pid for %"PRIu64"particles\nABORTING\n",isimu,(long)ihalo,(nPartInUse+1));
            exit(-1);
          }
          halos[isimu][ihalo].Pid[nPartInUse] = Pid;
          
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax[isimu])       PidMax[isimu] = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin)              PidMin        = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] > PidMax_local)        PidMax_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          if(halos[isimu][ihalo].Pid[nPartInUse] < PidMin_local)        PidMin_local  = (halos[isimu][ihalo].Pid[nPartInUse]);
          
          nPartInUse++;
        }
      }
      
      /* store number of particles in halo */
      halos[isimu][ihalo].npart = nPartInUse;
    }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
  
  fclose(fpin);
#ifdef READ_HALOIDS_FROM_FILE
  fclose(fphids);
#endif
  
#endif // READ_MPARTICLES
  
  elapsed = clock()-elapsed;
  
  fprintf(stderr," done in %4.2f sec. (nhalos = %"PRIu64", full ID range = %"PRIu64" -> %"PRIu64", local ID range = %"PRIu64" -> %"PRIu64")\n",
          (float)elapsed/CLOCKS_PER_SEC,nHalos[isimu],PidMin,PidMax[isimu],PidMin_local,PidMax_local);
  
  
  return(1);
}


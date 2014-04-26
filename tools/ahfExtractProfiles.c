#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*==================================================================================================
 * main:
 *
 *      extract a given set of haloes from an AHF analysis writing for each a file containing
 *      the *.AHF_halos line as well as the corresponding *.AHF_profiles columns
 *
 *==================================================================================================*/
int main(argc,argv)
int argc;
char **argv;
{
    char  halolist[2048], AHFhalos[2048], AHFprofiles[2048], outfile[2048];
    char  line[2048], halosheader[2048], profilesheader[2048];
    long *haloID, nhalos, jhalo, ihalo, khalo, nbins, ibin;
    float dummy;
    
    FILE *fp_halolist, *fp_AHFhalos, *fp_AHFprofiles, *fp_outfile;
    
   /*===================================================================
    * check consistency
    *===================================================================*/
   if(argc<3)
     {
      fprintf(stderr,"usage: %s AHFprefix halolist\n",*argv);
      exit(1);
     }
   
    /*===================================================================
     * prepare filenames
     *===================================================================*/
    strcpy(AHFhalos,    argv[1]);
    strcpy(AHFprofiles, argv[1]);
    strcpy(halolist,    argv[2]);
    
    strcat(AHFhalos,    ".AHF_halos");
    strcat(AHFprofiles, ".AHF_profiles");
    
    fprintf(stderr,"=========================================================\n");
    fprintf(stderr,"                  ExtractProfiles\n");
    fprintf(stderr,"=========================================================\n");
    fprintf(stderr," file with halo IDs: \t%s\n",halolist);
    fprintf(stderr," reading  halos   from %s\n",AHFhalos);
    fprintf(stderr," reading profiles from %s\n",AHFprofiles);
    fprintf(stderr,"\n");
    
    
    /*===================================================================
     * get list of halo IDs from file
     *===================================================================*/
    if( !(fp_halolist = fopen(halolist,"r")) )
       {
          fprintf(stderr," could not open %s\nABORTING\n",halolist);
          exit(0);
       }
    
    fgets(line,2048,fp_halolist);
    sscanf(line,"%ld",&nhalos);
    
    haloID = (long *) calloc(nhalos,sizeof(long));
    
    for(ihalo=0; ihalo<nhalos; ihalo++)
        fscanf(fp_halolist,"%ld",&(haloID[ihalo]));

    fclose(fp_halolist);

    /*===================================================================
     * find halos in *.AHF_halos and *.AHF_profiles files
     *===================================================================*/
    for(jhalo=0; jhalo<nhalos; jhalo++)
     {
        /* get halo from list */
        ihalo = haloID[jhalo];
        
        /* be verbose */
        fprintf(stderr," extracting halo #%ld\n",ihalo);
        
        /* open files */
        if( !(fp_AHFhalos = fopen(AHFhalos,"r")) )
         {
            fprintf(stderr," could not open %s\nABORTING\n",AHFhalos);
            exit(0);
         }
        if( !(fp_AHFprofiles = fopen(AHFprofiles,"r")) )
         {
            fprintf(stderr," could not open %s\nABORTING\n",AHFprofiles);
            exit(0);
         }
        
        sprintf(outfile,"halo_%08ld.dat",ihalo);
        if(!(fp_outfile = fopen(outfile,"w")))
         {
            fprintf(stderr," could not open %s\nABORTING\n",outfile);
            exit(0);
        }
        
        /* read header lines */
        fgets(halosheader,2048,fp_AHFhalos);
        halosheader[strlen(halosheader)-1] = '\0';
        
        fgets(profilesheader,2048,fp_AHFprofiles);
        profilesheader[strlen(profilesheader)-1] = '\0';
        
        /* skip irrelevant halos */
        for(khalo=0; khalo<ihalo; khalo++)
         {
            /* be verbose */
          //fprintf(stderr,"      skipping halo #%ld\n",khalo);
            
            /* skip *.AHF_halos */
            fgets(line,2048,fp_AHFhalos);
            
            /* get nbins from line[] string */
            sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %ld",
                   &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
                   &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
                   &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
                   &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
                   &nbins);

            /* skip *.AHF_profiles */
            for(ibin=0; ibin<nbins; ibin++)
                fgets(line,2048,fp_AHFprofiles);
         }
        
        /* reading and writing actual halo */
        fgets(line,2048,fp_AHFhalos);
        line[strlen(line)-1] = '\0';
        fprintf(fp_outfile,"%s\n",halosheader);
        fprintf(fp_outfile,"%s\n",line);
        
        /* get nbins from line[] string */
        sscanf(line,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %ld",
               &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
               &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
               &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
               &dummy,&dummy,&dummy,&dummy,&dummy,&dummy,
               &nbins);

        fprintf(fp_outfile,"# %ld, %s\n",nbins,profilesheader);
        for(ibin=0; ibin<nbins; ibin++)
         {
            fgets(line,2048,fp_AHFprofiles);
            line[strlen(line)-1] = '\0';
            fprintf(fp_outfile,"%s\n",line);  
         }
       
        
        /* close files */
        fclose(fp_AHFhalos);
        fclose(fp_AHFprofiles);
        fclose(fp_outfile);
        
        fprintf(stderr,"\n");
     }
    
    free(haloID);

}
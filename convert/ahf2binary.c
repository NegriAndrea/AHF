#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

/*===================
 * main()
 *===================*/
int main(argc,argv)
int argc;
char **argv;
{   
  char   outfile[2048], indata[2048];
  FILE  *fpin, *fpout;
  float *Var;
  int    i, NCOLUMNS, NWRITE;

  fprintf(stderr,"==========================================================\n");
  fprintf(stderr," generate a binary file out of a _halos or _profiles file\n");
  fprintf(stderr," (note, the generated _bin file has a different format as\n");
  fprintf(stderr,"    the ones generated with AHF when using -DAHFbinary!)\n");
  fprintf(stderr,"==========================================================\n");

  
  if(argc<4)
    {
      fprintf(stderr,"usage: %s AHF_halos NCOLUMNS NWRITE\n",*argv);
      exit(1);
    }
  
  NCOLUMNS = atoi(argv[2]);
  NWRITE   = atoi(argv[3]);
  
  strcpy(outfile,argv[1]);
  strcat(outfile,"_bin");
  fprintf(stderr,"generating %s\n",outfile);
  fpout = fopen(outfile,"wb");


  Var = (float *) calloc(NCOLUMNS, sizeof(float));
  
  // read AHF file ... and simultaneously write binary file
  fpin  = fopen(argv[1],"r");
  fgets(indata, 2048, fpin);
  
  while(!feof(fpin))
    {
      for(i=0; i<NCOLUMNS; i++)
        {
          fscanf(fpin,"%f",&Var[i]);
          //fprintf(stderr,"%f ",Var[i]);
        }
      //fprintf(stderr,"\n");
      
      fwrite((void *)Var, sizeof(float), NWRITE, fpout);
    }

  fclose(fpin);
  fclose(fpout);
}

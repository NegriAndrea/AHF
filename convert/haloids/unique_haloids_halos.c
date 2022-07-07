#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>

#define MAXSTRING 2048

int main(int argc, char **argv)
{
  char     inname[MAXSTRING], outname[MAXSTRING], line[MAXSTRING], pattern[MAXSTRING], *rest;
  FILE    *fpin, *fpout;
  int32_t  nhalos, npart, ptype, ihalo, snapid, haloid, hostid, numSubStruct;
  uint64_t uniquehaloid, uniquehostid, pid, ipart;
  double   mass;
  
  snapid = atoi(argv[2]);
  sprintf(inname,  "%s",argv[1]);
  sprintf(outname, "%s_uniquehaloids",basename(inname));
  
  fpin  = fopen(inname, "r");
  fpout = fopen(outname,"w");
  
  // header line
  fgets(line,MAXSTRING,fpin);
  fprintf(fpout,"%s",line);
  
  // loop over all halos
  fgets(line,MAXSTRING,fpin);
  while(!feof(fpin)) {

    sscanf(line,"%"SCNi32" %"SCNi32" %"SCNi32" %lf", &haloid, &hostid, &numSubStruct, &mass);
    uniquehaloid = (uint64_t)snapid*1000000000000+(haloid+1);
    if(hostid<0) {
      uniquehostid = 0;
    }
    else {
      uniquehostid = (uint64_t)snapid*1000000000000+(hostid+1);
    }
    
    // the mass column is the best string to search for when trying to locate 'the rest'
    sprintf(pattern,"%12.6g",mass);
    rest = strstr(line,pattern);
    fprintf(fpout,"%"PRIu64"  %"PRIu64"  %"PRIi32" %s", uniquehaloid, uniquehostid, numSubStruct, rest);

    fgets(line,MAXSTRING,fpin);
  }
  
  fclose(fpin);
  fclose(fpout);
}

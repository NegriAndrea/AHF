#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#define MAXSTRING 2048

int main(int argc, char **argv)
{
    char     inname[MAXSTRING], outname[MAXSTRING], line[MAXSTRING];
    FILE    *fpin, *fpout;
    int32_t  nhalos, npart, ptype, ihalo, snapid;
    uint64_t haloid, uniquehaloid, pid, ipart;
    
    snapid = atoi(argv[2]);
    sprintf(inname,  "%s",argv[1]);
    sprintf(outname, "%s_uniquehaloids",inname);
    
    fpin  = fopen(inname, "r");
    fpout = fopen(outname,"w");

    fgets(line,MAXSTRING,fpin);
    sscanf(line,"%"SCNi32, &nhalos);
    fprintf(fpout,"%"PRIi32"\n",nhalos);

    for(ihalo=0; ihalo<nhalos; ihalo++){
        fgets(line,MAXSTRING,fpin);
        sscanf(line,"%"SCNi32" %"SCNi32, &npart, &haloid);
        uniquehaloid = (uint64_t)snapid*1000000000000+(haloid+1);
        fprintf(fpout,"%"PRIi32" %"PRIu64"\n",npart, uniquehaloid);
        for(ipart=0; ipart<npart; ipart++) {
            fgets(line,MAXSTRING,fpin);
            sscanf(line,"%"SCNu64" %"SCNi32, &pid, &ptype);
            fprintf(fpout,"%"PRIu64" %"PRIi32"\n",pid, ptype);
        }
    }
    
    fclose(fpin);
    fclose(fpout);
}

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

#include "../src/param.h"
#include "../src/tdef.h"
#include "../src/common.c"

#include "../src/libutility/utility.h"

// & = AND
// | = OR
// ^ = XOR
// <<
// >>

HALO *halos;

void uint64_to_binary(uint64_t x);

int main (int argc, char *argv[]) 
{
  uint64_t ID;
  uint64_t Npart, Nid, Xid, Yid, Zid;
  double Xc,Yc,Zc,BoxSize;

  if(argc != 6) {
    fprintf(stderr,"usage: ./ahfCalcHaloID BoxSize Npart X Y Z\n");
    exit(0);
  }
  
  BoxSize = (double)   atof(argv[1]);
  Npart   = (uint64_t) atoi(argv[2]);
  Xc      = (double)   atof(argv[3]);
  Yc      = (double)   atof(argv[4]);
  Zc      = (double)   atof(argv[5]);

  Xc  /= BoxSize;
  Yc  /= BoxSize;
  Zc  /= BoxSize;

  halos = (HALO *) calloc(1, sizeof(HALO));
  halos[0].npart  = Npart;
  halos[0].pos.x  = Xc;
  halos[0].pos.y  = Yc;
  halos[0].pos.z  = Zc;

  ID = getHaloID(halos, 0);
  fprintf(stderr,"ID = %22"PRIu64"\n",ID);  
}

void uint64_to_binary(uint64_t x)
{
  uint64_t cnt, mask;
  mask = 1LU << 63;
  
  for(cnt=1;cnt<=64;++cnt)
    {
      putchar(((x & mask) == 0) ? '0' : '1');
      x <<= 1;
      if(cnt % 8 == 0 && cnt != 64)
	putchar(' ');
      if(cnt == 64)
	putchar('\n');
    }
}

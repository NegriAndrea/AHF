#include <stdio.h>
#include <stdlib.h>
#include "cubekey.h"
#include <stdint.h>
#include <string.h>

int main(int argc, char* argv[]){
  int i,j,depth,flag_pos=0;
  uint8_t pos;
  cubekey_t ck, ck_aux, *pck=NULL;
  flouble x,y,z,edge;
  double edge_test=1.0;
  char string[256], string2[256];
  ck_adjacents_t ck_adj;
  ck_shift_t ckX,ckY,ckZ;
  //
  //depth=3;
  //ckX=9;
  //for(pos=(uint8_t)depth;;pos--){
    //fprintf(stdout,"%d ",get_bit64(ckX,pos));
    //if(pos==0) break;
  //}
  //fprintf(stdout,"\n");
  //return 0;
  
  
  //while(1){
    if(argc==2){
      strncpy(string,argv[1],256);
    }
    else{
      fprintf(stdout,"Enter cubekey: ");
      fgets(string, 256, stdin);
    }

#ifdef CUBEKEY_128
    ck=strtoull(string, &string2,10);
#else
    ck=atoll(string);
#endif
    ck=0;
    fprintf(stderr, "TODO CEROS:\n");
    fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
    fprintf(stderr, "lz=%3d, ck_get_flag_pos(ck)=%d, ck_get_depth=%d\n", clz(ck), ck_get_flag_pos(ck), ck_get_depth(ck));

    ck=~ck;
    fprintf(stderr, "NEGADO:\n");
    fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
    fprintf(stderr, "lz=%3d, ck_get_flag_pos(ck)=%d, ck_get_depth=%d\n", clz(ck), ck_get_flag_pos(ck), ck_get_depth(ck));


    ck=~ck;

    for(i=0;i<CUBEKEY_BITS-1;i+=2){
      fprintf(stderr, "Set pos %d:\n", i);
      set_bit_ck(&ck,i);
      fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
      fprintf(stderr, "lz=%3d, ck_get_flag_pos(ck)=%d, ck_get_depth=%d, get_bit_ck(%d)=%d, get_bit_ck(%d)=%d\n", clz(ck), ck_get_flag_pos(ck), ck_get_depth(ck), i, get_bit_ck(ck,i), i+1, get_bit_ck(ck,i+1));
      ck=~ck;
      fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
      fprintf(stderr, "lz=%3d, ck_get_flag_pos(ck)=%d, ck_get_depth=%d, get_bit_ck(%d)=%d, get_bit_ck(%d)=%d\n", clz(ck), ck_get_flag_pos(ck), ck_get_depth(ck), i, get_bit_ck(ck,i), i+1, get_bit_ck(ck,i+1));
      ck=~ck;
    }

    fprintf(stderr, "Ponemos a 1 solamente el bit-pos 127\n");
    ck = 1;
    ck <<= 127;
    fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
    if(ck==0) fprintf(stderr, " el bit 127 a 1 y ck==0 devuelve TRUE\n");
    else fprintf(stderr, " el bit 127 a 1 y ck==0 devuelve FALSE\n");


    exit(1);


    fprintf(stderr,"strtoull()-> ck=%llu\n", ck);
    fprintbits_cubekey(stderr,ck);fprintf(stderr, "\n");
    fprintf(stderr, "NEGADO:\n");
    fprintbits_cubekey(stderr,~ck);fprintf(stderr, "\n");
//    fprintbits_cubekey_3block(stderr,ck);


    /*
    for(i=1;i<MAX_DEPTH;i++){
      edge_test=edge_test/2.0;
      fprintf(stderr,"edge_test=%e, depth=%d\n", edge_test,i);
    }
    exit(1);
    */

    set_bit_ck(&ck,126);
    flag_pos=ck_get_flag_pos(ck);

    fprintf(stderr,"flag_pos(ck) returned %d\n", flag_pos);

    if(flag_pos==0 || flag_pos%3!=0){
      fprintf(stdout, "INVALID cubekey entered: %"PRIck"\n",ck);
      fprintbits_cubekey(stdout,ck);
      fprintf(stdout, "\nSizeof cubekey=           %lu\n",sizeof(cubekey_t));
      fprintf(stdout, "Sizeof int=                 %lu\n",sizeof(int));
      fprintf(stdout, "Sizeof long int=            %lu\n",sizeof(long int));
      fprintf(stdout, "Sizeof long long int=       %lu\n",sizeof(long long int));
      return -1;
    }
    
    


    ck2coor(&x,&y,&z,&edge,ck);
    ck_get_shifts(ck, &ckX, &ckY, &ckZ);
    depth=ck_get_depth(ck);
    fprintf(stdout, "==ENTERED== %15"PRIu64" (depth %2d) cubekey-> [ %e,   %e,   %e ] (%e edge)[SHIFTS",ck,depth,x,y,z,edge);
    fprintf(stdout," ");
    for(pos=(uint8_t)depth-1;;pos--){
      fprintf(stdout,"%d",get_bit64(ckX,pos));
      if(pos==0) break;
    }
    fprintf(stdout," ");
    for(pos=(uint8_t)depth-1;;pos--){
      fprintf(stdout,"%d",get_bit64(ckY,pos));
      if(pos==0) break;
    }
    fprintf(stdout," ");
    for(pos=(uint8_t)depth-1;;pos--){
      fprintf(stdout,"%d",get_bit64(ckZ,pos));
      if(pos==0) break;
    }
    fprintf(stdout,"]\n");
    
    //fprintf(stdout, "x = %10lf\n",x);
    //fprintf(stdout, "y = %10lf\n",y);
    //fprintf(stdout, "z = %10lf\n",z);
    //fprintf(stdout, "edge = %10f\n",edge);
    
    memset(&ck_adj,0,sizeof(ck_adjacents_t));
    
    ck_get_adjacents_side_edge(ck,&ck_adj);
    pck=&ck_adj.ckBBB;
    j=0;
    for(i=0; i<NUM_ADJACENT;i++){
      if(*pck!=0){
        j++;
        ck2coor(&x,&y,&z,&edge,*pck);
        ck_get_shifts(*pck, &ckX, &ckY, &ckZ);
        fprintf(stdout, " adj   [%2d/%d %2d] %15"PRIu64" cubekey->\t%4e,   %4e,   %4e   (%4e edge)   [SHIFTS ",i,NUM_ADJACENT,j,*pck,x,y,z,edge);
        for(pos=(uint8_t)depth-1;;pos--){
          fprintf(stdout,"%d",get_bit64(ckX,pos));
          if(pos==0) break;
        }
        fprintf(stdout," ");
        for(pos=(uint8_t)depth-1;;pos--){
          fprintf(stdout,"%d",get_bit64(ckY,pos));
          if(pos==0) break;
        }
        fprintf(stdout," ");
        for(pos=(uint8_t)depth-1;;pos--){
          fprintf(stdout,"%d",get_bit64(ckZ,pos));
          if(pos==0) break;
        }
        fprintf(stdout,"]\n");
      }
      pck++;
    }
    
  //}//while
  
}

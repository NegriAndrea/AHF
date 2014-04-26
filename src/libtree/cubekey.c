/*
 * cubekey.c
 * 
 * Copyright 2013 Fernando Campos <fernando.campos@uam.es>
 * 
 * This file will contain the functions to convert from (and to) a 
 * 64bits integer to the shifts on each dimension wich will
 * uniquely indentify a subcube in the volume. We will also need
 * the depth of the subcube. 
 * 
 * The uint64 variable will contain 21 bits per dimension (x,y,z)
 * and 1 unused bit (left most)
 * 
 */
#include "cubekey.h"
#include "stdio.h"
#include "stdlib.h"
#include <stdint.h>
#include <inttypes.h>
#include <string.h>

#include "../common.h"
#include "../tdef.h"

inline void set_bit64(uint64_t *x, uint8_t pos){
	*x |= (1L << pos);
}

inline void set_bit32(uint32_t *x, uint8_t pos){
	*x |= (1 << pos);
}

inline void clr_bit64(uint64_t *x, uint8_t pos){
	*x &= ~(1L << pos);
}

inline void clr_bit32(uint32_t *x, uint8_t pos){
	*x &= ~(1 << pos);
}

inline int get_bit64(uint64_t x, uint8_t pos) {
    return x & (1L << pos);
}

inline int get_bit32(uint32_t x, uint8_t pos) {
    return x & (1 << pos);
}

/*
 * TO DO: Use intrinsic clz for 64 bits
 */
#ifdef __GNUC__
	#include <x86intrin.h>
	#define clz(x) /*fprintf(stderr,"clzl builtin\n");*/ \
								__builtin_clzl(x)
//#define ctz(x) __builtin_ctz(x)
#else
inline int clz(uint64_t x){
	int n=0;
	if (x == 0) return 64;
	if ((x & 0xFFFFFFFF00000000) == 0){ n+=32; x<<=32;}
	if ((x & 0xFFFF000000000000) == 0){ n+=16; x<<=16;}
	if ((x & 0xFF00000000000000) == 0){ n+=8;  x<<= 8;}
	if ((x & 0xF000000000000000) == 0){ n+=4;  x<<= 4;}
	if ((x & 0xC000000000000000) == 0){ n+=2;  x<<= 2;}
	if ((x & 0x8000000000000000) == 0){ n+=1;  x<<= 1;}
	fprintf(stderr,"clz not builtin\n");
	return n;
}
#endif

inline int ck_get_flag_pos(cubekey_t ck){
	int lz;

	lz=clz(ck);
	if(lz==64)
		return 64;
	else
		return 63-lz;
}

inline int ck_get_depth(cubekey_t ck){
	return ck_get_flag_pos(ck)/3;
}

void printbits_cubekey(cubekey_t n){
	cubekey_t i;
	int count=0;
	i = 0x8000000000000000;

	while (i > 0){
		if (n & i)
			printf("1");
		else
			printf("0");
		i >>= 1;
		count++;
		if(count%8==0)printf(" | ");
		else if(count%4==0)printf(" ");
	}
}
void fprintbits_cubekey(FILE* f, cubekey_t n){
	cubekey_t i;
	int count=0;
	i = 0x8000000000000000;

	while (i > 0){
		if (n & i)
			fprintf(f,"1");
		else
			fprintf(f, "0");
		i >>= 1;
		count++;
		if(count%8==0) fprintf(f," | ");
		else if(count%4==0) fprintf(f," ");
	}
}

void printbits_cubekey_3block(cubekey_t n){
	cubekey_t i;
	int count=0,lz=clz(n);
	i = 0x8000000000000000;
	
	i>>=lz+1;
	count+=lz+1;
	printf("*|");
	while (i > 0){
		if (n & i)
			printf("1");
		else
			printf("0");
		i >>= 1;
		count++;
		if(count==1)printf("|");
		else if(count%3==1)printf("|");
	}
	printf("\b- ");
}


void fprintbits_cubekey_3block(FILE* f,cubekey_t n){
	cubekey_t i;
	int count=0,lz=clz(n);
	i = 0x8000000000000000;
	
	i>>=lz+1;
	count+=lz+1;
	fprintf(f,"*|");
	while (i > 0){
		if (n & i)
			fprintf(f,"1");
		else
			fprintf(f,"0");
		i >>= 1;
		count++;
		if(count==1)fprintf(f,"|");
		else if(count%3==1)fprintf(f,"|");
	}
	fprintf(f,"\b- ");
}

void printbits_ckcoor(uint32_t n){
	uint32_t i;
	int count=0;
	i = 0x80000000;

	while (i > 0){
		if (n & i)
			printf("1");
		else
			printf("0");
		i>>=1;
		count++;
		if(count%8==0)printf(" | ");
		else if(count%4==0)printf(" ");
	}
}

/*
 * Given a depth and 3 coordinates, return the corresponding key
 */
int_fast8_t coor2ck(uint64_t* ckey, const flouble xcoor, const flouble ycoor, const flouble zcoor, const uint_fast8_t depth){
	uint64_t x=0,y=0,z=0;
	uint64_t i=0;
	flouble cur_x=0, cur_y=0, cur_z=0;
	uint64_t local_ckey=0;
	//float edge_length=getBoxSize();
	flouble edge_length=1.0;
	flouble half;
	*ckey=0;

//Parameters check
	if(depth>MAX_DEPTH){
		fprintf(stderr,"depth (%"PRIuFAST8") bigger than MAX_DEPTH\n", depth);
		return ERRDEPTH;
	}
	if(depth<1){
		fprintf(stderr,"depth (%"PRIuFAST8") smaller than 1\n", depth);
		return ERRDEPTH;
	}
	if(ckey==NULL){
		fprintf (stderr, "path_key NULL\n");
		return ERRCKEY;
	}
	#ifdef LIBTREE_DEBUG
	fprintf(stderr, "BoxSize=%f, depth=%"PRIuFAST32", xcoor=%f, ycoor=%f, zcoor=%f\n", edge_length, depth, xcoor, ycoor, zcoor);
	#endif

	for (i=0; i<depth; i++){
		half=edge_length/2;
		//Clear previous depth-FLAG
		clr_bit64(&local_ckey,3*(i));
		//Set Depth-FLAG
		set_bit64(&local_ckey,3*(i+1));

		#ifdef LIBTREE_DEBUG
			fprintf(stderr, "Setting flat at pos %llu\t[i=%"PRIu64"]\n", 3*(i+1), i);
		#endif

		#ifdef LIBTREE_DEBUG
			fprintf(stderr,"X (depth=%"PRIu64")\t%f > %f (%f+%f):\n", i+1, xcoor, cur_x+half, cur_x, half);
		#endif
		if (xcoor >= cur_x + half){
			//x = x | 1<<i;
			set_bit64(&x,i);
			cur_x+=half;
			#ifdef LIBTREE_DEBUG
				fprintf(stderr,"\tcur_x=%F\tX = %#08llX MOVE\n", cur_x, x);
			#endif
			set_bit64(&local_ckey,3*i+2);
		}

		#ifdef LIBTREE_DEBUG
			fprintf(stdout,"Y (depth=%"PRIu64")\t%f > %f (%f+%f):\n", i+1, ycoor, cur_y+half, cur_y, half);
		#endif
		if (ycoor > cur_y + half){
			set_bit64(&y,i);
			cur_y+=half;
			#ifdef LIBTREE_DEBUG
				fprintf(stdout,"\tcur_y=%f\tY = %#08llX MOVE\n", cur_y, y);
			#endif
			set_bit64(&local_ckey,3*i+1);
		}

		#ifdef LIBTREE_DEBUG
			fprintf(stdout,"Z (depth=%"PRIu64")\t%f > %f (%f+%f):\n", i+1, zcoor, cur_z+half, cur_z, half);
		#endif
		if (zcoor > cur_z + half){
			set_bit64(&z,i);
			cur_z+=half;
			#ifdef LIBTREE_DEBUG
				fprintf(stdout,"\tcur_z=%f\tZ = %#08llX MOVE\n", cur_z, z);
			#endif
			set_bit64(&local_ckey,3*i);
		}
		//Update edge_length for the next iteration
		edge_length=half;
		#ifdef LIBTREE_DEBUG
		printf ("key=\t");printbits_uint64(local_ckey); printf("\n");
		#endif
	}
	#ifdef LIBTREE_DEBUG
	fprintf(stdout,"______________\n");
	printf ("key=\t");printbits_uint64(local_ckey); printf("\n");
	#endif
	*ckey=local_ckey;
	return 0;
}

void ck2coor(flouble* x, flouble* y, flouble* z, flouble* edge, uint64_t cubekey){
	uint_fast8_t flag_pos=0;
	flouble half;

	/*
	 * TO-DO: create a ck_is_valid function
	 */
	flag_pos=ck_get_flag_pos(cubekey);
	if(flag_pos%3!=0){
		fprintf(stderr, "[ck2coor] ERROR: flag_pos not valid. cubekey probably wrong.cubekey=%lu, flag_pos=%u\n",cubekey,flag_pos);
		*x=*y=*z=-1.0;
		*edge=1;
		return;
	}

	half=BOX_SIZE;
	*x=*y=*z=0.0;
	while(cubekey>1){
		half=half/2;
		if(cubekey& 0x01)
			*z+=half;
		if(cubekey& 0x02)
			*y+=half;
		if(cubekey& 0x04)
			*x+=half;
		cubekey>>=3;
	}
	*edge=half;
}

void ck2coor_center(flouble* x, flouble* y, flouble* z, cubekey_t ck){
	flouble edge;
	ck2coor(x,y,z,&edge,ck);
	edge=edge/2;
	*x+=edge;
	*y+=edge;
	*z+=edge;
	return;
}


inline void ck_get_daughters(cubekey_t ck, ck_daughters_t* ck_daughters){
  int depth=ck_get_depth(ck);
	int ckFlagPos=ck_get_flag_pos(ck);

	if(depth==MAX_DEPTH){
		fprintf(stderr,"[ck_get_daughters] cubekey has maximum depth (%d). Not possible to get daughters ID's\n",depth);
		memset(ck_daughters,0, sizeof(cubekey_t)*NUM_DAUGHTERS);
		//ck_daughters->ck000=0;
		//ck_daughters->ck001=0;
		//ck_daughters->ck010=0;
		//ck_daughters->ck011=0;
		//ck_daughters->ck100=0;
		//ck_daughters->ck101=0;
		//ck_daughters->ck110=0;
		//ck_daughters->ck111=0;
		return;
	}
	//printf("\nck   :     ");
	//printbits_cubekey_3block(ck);
	//ck2coor(&x,&y,&z,&edge,ck);
	//fprintf(stdout,"\t%4.8f %4.8f %4.8f (%4.8f)", x,y,z,edge);

	//Move Flag from original cubekey to deeper level. Will be the base for daughters ck.
	//ck is passed by value, no changes outside this function
	clr_bit64(&ck,ckFlagPos);
	set_bit64(&ck,ckFlagPos+3);
	
	//ck000
	memcpy(&(ck_daughters->ck000),&ck,SIZEOF_CUBEKEY);
	//ck001
	memcpy(&(ck_daughters->ck001),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck001),ckFlagPos); //Old Flag-pos is Z-pos now
	//ck010
	memcpy(&(ck_daughters->ck010),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck010),ckFlagPos+1);
	//ck011
	memcpy(&(ck_daughters->ck011),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck011),ckFlagPos+1);
	set_bit64(&(ck_daughters->ck011),ckFlagPos);
	
	//Set X to 1 in original ck to save set_bit calls
	set_bit64(&ck,ckFlagPos+2);

	//ck100
	memcpy(&(ck_daughters->ck100),&ck,SIZEOF_CUBEKEY);
	//ck101
	memcpy(&(ck_daughters->ck101),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck101),ckFlagPos); //Old Flag-pos is Z-pos now
	//ck110
	memcpy(&(ck_daughters->ck110),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck110),ckFlagPos+1);
	//ck111
	memcpy(&(ck_daughters->ck111),&ck,SIZEOF_CUBEKEY);
	set_bit64(&(ck_daughters->ck111),ckFlagPos+1);
	set_bit64(&(ck_daughters->ck111),ckFlagPos);
}

inline void ck_get_parent(cubekey_t ck, cubekey_t* ck_parent){
  int depth=ck_get_depth(ck);
	int ckFlagPos=ck_get_flag_pos(ck);

	if(depth==1){
		*ck_parent=0;
		return;
	}
	//printf("\nck   :     ");
	//printbits_cubekey_3block(ck);
	//ck2coor(&x,&y,&z,&edge,ck);
	//fprintf(stdout,"\t%4.8f %4.8f %4.8f (%4.8f)", x,y,z,edge);

	/*Move Flag from original cubekey to less deep level (ckFlagPos-3) 
	 * and clear the in between bits: ckFlagPos,ckFlagPos-1(X),ckFlagPos-2(Y).
	 * ckFlagPos-3, X bit in original cubekey, will be the Flag for the parent.
	 */
	memcpy(ck_parent,&ck,SIZEOF_CUBEKEY);
	clr_bit64(ck_parent,ckFlagPos);
	clr_bit64(ck_parent,ckFlagPos-1);
	clr_bit64(ck_parent,ckFlagPos-2);
	set_bit64(ck_parent,ckFlagPos-3);
	
	//fprintf(stderr,"\n[ck_get_parent] original: ");fprintbits_cubekey_3block(stderr,ck);fprintf(stderr,"\t parent: ");fprintbits_cubekey_3block(stderr,*ck_parent);fprintf(stderr,"\n");
}

inline void ck_get_shifts(cubekey_t ck, uint32_t* ckX,uint32_t* ckY,uint32_t* ckZ){
	uint32_t ckFlagPos=ck_get_flag_pos(ck);
	int i;
	*ckX=*ckY=*ckZ=0;
	
	for(i=0;i<ckFlagPos;i+=3){
		//Shift to the left the corresponding per-coordinate cubekey and copy the bit from cubekey
		//Z bit
		*ckZ<<=1;
		if(get_bit64(ck,i)!=0){
				set_bit32(ckZ,0);
		}
		//Y bit
		*ckY<<=1;
		if(get_bit64(ck,i+1)!=0){
			set_bit32(ckY,0);
		}
		//X bit
		*ckX<<=1;
		if(get_bit64(ck,i+2)!=0){
			set_bit32(ckX,0);
		}
	}//for i
	
}

inline void ck_get_adjacents(cubekey_t ck, ck_adjacents_t* ck_adjacents){
  int i=0;
  cubekey_t *p_ck_iter=NULL;
  int depth=ck_get_depth(ck);
  uint32_t ckX=0,ckY=0,ckZ=0,ckFlagPos=ck_get_flag_pos(ck);
  uint32_t ckX_B,ckX_S,ckX_F;
  uint32_t ckY_B,ckY_S,ckY_F;
  uint32_t ckZ_B,ckZ_S,ckZ_F;
  uint32_t module;
  
	//Point cubekey pointer iterator to the first adjacent cubekey
	p_ck_iter=&ck_adjacents->ckBBB;
	//Initialize all 26 adjacent subcubes to 0
	memset(p_ck_iter,0,SIZEOF_CUBEKEY*NUM_ADJACENT);
	
	
	//In the special case depth==1, any subcube has 7 adjacent subcubes. Let's calculate them
	if(depth==1){
		//fprintf(stderr,"\n[ck_get_adjacents] Original cubekey: ");fprintbits_cubekey_3block(stderr,ck);
		//Clear Flag-bit
		clr_bit64(&ck,ckFlagPos);		
		for(i=0;i<7;i++){ //7 is the number of adjacent subcubes
			*p_ck_iter=(ck+i+1)%8;
			set_bit64(p_ck_iter,ckFlagPos);
			//fprintf(stderr,"\n[ck_get_adjacents] adjacent %d: ",i+1);;fprintbits_cubekey_3block(stderr,*p_ck_iter);
			p_ck_iter++;
		}
		//fprintf(stderr,"\n");
		//fprintf(stderr,"[ck_get_adjacents] depth < 2 (%d). Make no sense to ask for adjacents cubekeys at this level\n",depth);
		return;
	}

	/*
	 * Extract shifts per dimension (X, Y and Z) from cubekey
	 * Since we insert the most significant bits from right to left in the cubekey,
	 * we will extract it from rigth to left and insert in the per-coordinate cubekeys (ckX,ckY,ckZ)
	 * shifting to the left.
	 */
	
	ck_get_shifts(ck,&ckX,&ckY,&ckZ);
	//
	//for(i=0;i<ckFlagPos;i+=3){
		//Shift to the left the corresponding per-coordinate cubekey and copy the bit from cubekey
		//Z bit
		//ckZ<<=1;
		//if(get_bit64(ck,i)!=0){
				//set_bit32(&ckZ,0);
		//}
		//Y bit
		//ckY<<=1;
		//if(get_bit64(ck,i+1)!=0){
			//set_bit32(&ckY,0);
		//}
		//X bit
		//ckX<<=1;
		//if(get_bit64(ck,i+2)!=0){
			//set_bit32(&ckX,0);
		//}
	//}//for i
	
	/*
	 * For periodic boundary condition, we must consider -1 (backward) and +1 (forward) shifts, even in the limits of the cube
	 * Applying module operation we will obtain, for example, 111 and 001 as adjacent to 000 (module 8).
	 * The module will depend on the refinement depth, this is the number of bits per dimension.
	 * i.e.: If we are using 4 bits per dimension (depth=4), module will be 2^4=16.
	 */
	module=1<<depth;
	ckX_S=ckX;
	ckX_B=(ckX-1)%(module);
	ckX_F=(ckX+1)%(module);
	
	ckY_S=ckY;
	ckY_B=(ckY-1)%(module);
	ckY_F=(ckY+1)%(module);
	
	ckZ_S=ckZ;
	ckZ_B=(ckZ-1)%(module);
	ckZ_F=(ckZ+1)%(module);
	
	//printf("\nckX: %d"); printbits_ckcoor(ckX);   printf("ckY:   "); printbits_ckcoor(ckY);   printf("ckZ:   "); printbits_ckcoor(ckZ);
	//printf("\nckX_B: "); printbits_ckcoor(ckX_B); printf("ckX_S: "); printbits_ckcoor(ckX_S); printf("ckX_F: "); printbits_ckcoor(ckX_F);
	//printf("\nckY_B: "); printbits_ckcoor(ckY_B); printf("ckY_S: "); printbits_ckcoor(ckY_S); printf("ckY_F: "); printbits_ckcoor(ckY_F);
	//printf("\nckZ_B: "); printbits_ckcoor(ckZ_B); printf("ckZ_S: "); printbits_ckcoor(ckZ_S); printf("ckZ_F: "); printbits_ckcoor(ckZ_F);

	//printf("\nckX:   %u ",ckX); printf("ckY:   %u ",ckY); printf("ckZ:   %u ",ckZ);
	//printf("\nckX_B: %u ",ckX_B);printf("ckX_S: %u ",ckX_S); printf("ckX_F: %u ",ckX_F);
	//printf("\nckY_B: %u ",ckY_B);printf("ckY_S: %u ",ckY_S); printf("ckY_F: %u ",ckY_F);
	//printf("\nckZ_B: %u ",ckZ_B);printf("ckZ_S: %u ",ckZ_S); printf("ckZ_F: %u ",ckZ_F);

//Compose cubekeys for 26 adjacent subcubes

	//Set Pos-Flag
	for(i=0;i<NUM_ADJACENT;i++){
		set_bit64(p_ck_iter,ckFlagPos);
		p_ck_iter++;
  }
	
	for(i=depth;i>0;i--){
		//Z Backward (at depth=i)
		if(get_bit32(ckZ_B,i-1)){
			set_bit64(&(ck_adjacents->ckBBB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSBB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFBB),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBFB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSFB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFFB),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBSB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSSB),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFSB),(depth-i)*3);
		}
		
		//Z Still (at depth=i)
		if(get_bit32(ckZ_S,i-1)){
			set_bit64(&(ck_adjacents->ckBBS),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSBS),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFBS),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBFS),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSFS),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFFS),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBSS),(depth-i)*3);
			//set_bit64(&(ck_adjacents->ckSSS),(i*3)-1); SSS is original cubekey-> No sense!
			set_bit64(&(ck_adjacents->ckFSS),(depth-i)*3);
		}
		
		//Z Forward (at depth=i)
		if(get_bit32(ckZ_F,i-1)){
			set_bit64(&(ck_adjacents->ckBBF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSBF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFBF),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBFF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSFF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFFF),(depth-i)*3);
			
			set_bit64(&(ck_adjacents->ckBSF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckSSF),(depth-i)*3);
			set_bit64(&(ck_adjacents->ckFSF),(depth-i)*3);
		}

/*************/

		//Y Backward (at depth=i)
		if(get_bit32(ckY_B,i-1)){
			set_bit64(&(ck_adjacents->ckBBB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSBB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFBB),(depth-i)*3+1);
			
			set_bit64(&(ck_adjacents->ckBBF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSBF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFBF),(depth-i)*3+1);
			
			set_bit64(&(ck_adjacents->ckBBS),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSBS),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFBS),(depth-i)*3+1);			
		}

		//Y Still (at depth=i)
		if(get_bit32(ckY_S,i-1)){
			set_bit64(&(ck_adjacents->ckBSB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSSB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFSB),(depth-i)*3+1);
			
			set_bit64(&(ck_adjacents->ckBSF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSSF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFSF),(depth-i)*3+1);
			
			set_bit64(&(ck_adjacents->ckBSS),(depth-i)*3+1);
			//set_bit64(&(ck_adjacents->ckSSS),(i*3)-2); SSS is original cubekey-> No sense!
			set_bit64(&(ck_adjacents->ckFSS),(depth-i)*3+1);			
		}

		//Y Forward (at depth=i)
		if(get_bit32(ckY_F,i-1)){
			set_bit64(&(ck_adjacents->ckBFB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSFB),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFFB),(depth-i)*3+1);
			
			set_bit64(&(ck_adjacents->ckBFS),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSFS),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFFS),(depth-i)*3+1);			
			
			set_bit64(&(ck_adjacents->ckBFF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckSFF),(depth-i)*3+1);
			set_bit64(&(ck_adjacents->ckFFF),(depth-i)*3+1);
		}
		
/*************/

		//X Backward (at depth=i)
		if(get_bit32(ckX_B,i-1)){
			set_bit64(&(ck_adjacents->ckBBB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBSB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBFB),(depth-i)*3+2);
			
			set_bit64(&(ck_adjacents->ckBBS),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBSS),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBFS),(depth-i)*3+2);			
			
			set_bit64(&(ck_adjacents->ckBBF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBSF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckBFF),(depth-i)*3+2);
		}

		//X Still (at depth=i)
		if(get_bit32(ckX_S,i-1)){
			set_bit64(&(ck_adjacents->ckSBB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckSSB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckSFB),(depth-i)*3+2);

			set_bit64(&(ck_adjacents->ckSBS),(depth-i)*3+2);
			//set_bit64(&(ck_adjacents->ckSSS),(depth-i)*3+2); SSS is original cubekey-> No sense!
			set_bit64(&(ck_adjacents->ckSFS),(depth-i)*3+2);			
			
			set_bit64(&(ck_adjacents->ckSBF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckSSF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckSFF),(depth-i)*3+2);			
		}

		//X Forward (at depth=i)
		if(get_bit32(ckX_F,i-1)){
			set_bit64(&(ck_adjacents->ckFBB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFSB),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFFB),(depth-i)*3+2);
			
			set_bit64(&(ck_adjacents->ckFBS),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFSS),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFFS),(depth-i)*3+2);
			
			set_bit64(&(ck_adjacents->ckFBF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFSF),(depth-i)*3+2);
			set_bit64(&(ck_adjacents->ckFFF),(depth-i)*3+2);			
		}
	}	
}

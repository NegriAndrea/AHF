/*
 * cubekey.h
 * 
 * Copyright 2014 F. Campos <fernandocdelpozo@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef _CUBEKEY_H_
#define _CUBEKEY_H_

#include <math.h>
#include <stdint.h>
#include "../common.h"

//Maximum depth in the subdivision is given by the 
//64bits for the position in key_t.pos divided
//by 3 dimensions

#ifdef CUBEKEY_128

#ifndef __GNUC__
  #error "__int128 only defined for GCC compiler"
#endif

  typedef unsigned __int128 cubekey_t;
  typedef uint64_t ck_shift_t;
  #define MAX_DEPTH 42
  #define CUBEKEY_BITS 128
  #define CUBEKEY_SHIFT_BITS 64
  #define PRIck "llu"

#else //Cubekey of 64 bits
  typedef uint64_t cubekey_t;
  typedef uint32_t ck_shift_t;
  #define MAX_DEPTH 21
  #define CUBEKEY_BITS 64
  #define CUBEKEY_SHIFT_BITS 32
  #define PRIck "lu"

#endif

#define SIZEOF_CUBEKEY (sizeof(cubekey_t))
#define SIZEOF_CK_SHIFT (sizeof(ck_shift_t))

#define NLEVELS MAX_DEPTH+1
#define NUM_DAUGHTERS 8
#define NUM_ADJACENT 26
#define BOX_SIZE 1.0L

#define ERRDEPTH -5
#define ERRCKEY -6

#ifdef CUBEKEY_128
  int  get_bit128(cubekey_t,uint8_t);
  void set_bit128(cubekey_t*,uint8_t);
  void clr_bit128(cubekey_t*,uint8_t);
  //64bits Cubekey_shift
  int  get_bit64(ck_shift_t, uint8_t);
  void set_bit64(ck_shift_t*,uint8_t);
  void clr_bit64(ck_shift_t*,uint8_t);

  /* 128bits Cubekeys, 64bits ck_shifts */
  #define clr_bit_ck clr_bit128
  #define get_bit_ck get_bit128
  #define set_bit_ck set_bit128

  #define clr_bit_ck_shift clr_bit64
  #define get_bit_ck_shift get_bit64
  #define set_bit_ck_shift set_bit64

#else
  int  get_bit64(cubekey_t, uint8_t);
  void set_bit64(cubekey_t*,uint8_t);
  void clr_bit64(cubekey_t*,uint8_t);

  //We only use 32bit words when not using 128bits cubekeys
  int  get_bit32(ck_shift_t, uint8_t);
  void set_bit32(ck_shift_t*,uint8_t);
  void clr_bit32(ck_shift_t*,uint8_t);

  /* 64bits Cubekeys, 32bits ck_shifts */
  #define clr_bit_ck clr_bit64
  #define get_bit_ck get_bit64
  #define set_bit_ck set_bit64

  #define clr_bit_ck_shift clr_bit32
  #define get_bit_ck_shift get_bit32
  #define set_bit_ck_shift set_bit32
#endif


typedef struct ck_adjacents_s{
//B for backward, S for Still, F for forward
//XYZ positions

//Combining Z -> Backward
  //Combining Y -> Backward
  cubekey_t ckBBB;
  cubekey_t ckSBB;
  cubekey_t ckFBB;
  //Combining Y -> Still
  cubekey_t ckBSB;
  cubekey_t ckSSB;
  cubekey_t ckFSB;
  //Combining Y -> Forward
  cubekey_t ckBFB;
  cubekey_t ckSFB;
  cubekey_t ckFFB;

//Combining Z -> Still
  //Combining Y -> Backward
  cubekey_t ckBBS;
  cubekey_t ckSBS;
  cubekey_t ckFBS;
  //Combining Y -> Still
  cubekey_t ckBSS;
//  cubekey_t ckSSS; Still for the 3 dimensions is the original subcube
  cubekey_t ckFSS;
  //Combining Y -> Forward
  cubekey_t ckBFS;
  cubekey_t ckSFS;
  cubekey_t ckFFS;

//Combining Z -> Forward 
  //Combining Y -> Backward
  cubekey_t ckBBF;
  cubekey_t ckSBF;
  cubekey_t ckFBF;
  //Combining Y -> Still
  cubekey_t ckBSF;
  cubekey_t ckSSF;
  cubekey_t ckFSF;
  //Combining Y -> Forward
  cubekey_t ckBFF;
  cubekey_t ckSFF;
  cubekey_t ckFFF;
  
} ck_adjacents_t;

typedef struct ck_daughters_s{
  cubekey_t ck000;
  cubekey_t ck001;
  cubekey_t ck010;
  cubekey_t ck011;
  cubekey_t ck100;
  cubekey_t ck101;
  cubekey_t ck110;
  cubekey_t ck111;
} ck_daughters_t;


int clz(cubekey_t x);
int ck_get_flag_pos(cubekey_t x);
int ck_get_depth(cubekey_t x);



//Generic functions to work with cubekeys, shifts and coordinates
int_fast8_t coor2ck(cubekey_t*, const flouble, const flouble, const flouble, const uint_fast8_t);
void ck2coor(flouble*, flouble*, flouble*, flouble*, cubekey_t);
void ck2coor_center(flouble*, flouble*, flouble*, cubekey_t);
void ck_get_daughters(cubekey_t, ck_daughters_t*);
void ck_get_adjacents(cubekey_t, ck_adjacents_t*);
void ck_get_adjacents_side_edge(cubekey_t, ck_adjacents_t*);
void ck_get_parent(cubekey_t, cubekey_t*);
void ck_get_shifts(cubekey_t, ck_shift_t*, ck_shift_t*, ck_shift_t*);

//Binary printing functions
void printbits_cubekey(cubekey_t);
void printbits_cubekey_3block(cubekey_t);
void fprintbits_cubekey_3block(FILE*,cubekey_t);
void fprintbits_cubekey(FILE*, cubekey_t);
void printbits_ckcoor(uint32_t);

#endif //_CUBEKEY_H_

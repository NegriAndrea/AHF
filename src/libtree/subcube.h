#ifndef _SUBCUBE_H_
#define _SUBCUBE_H_

//#include <math.h>
//#include <stdint.h>
//#include "../common.h"
#include "../tdef.h"
#include "cubekey.h"
#include "uthash.h"
#include "utarray.h"

#define BACKWARD	-1
#define STILL 		0
#define FORWARD		1

#define NO_PATCH         -1L
#define PATCH_REJECTED   -2L



typedef struct subcube_s{
  cubekey_t cubekey;
  UT_array *particles;    //particles contained in the subcube. NULL if the subcube have been refined
  uint64_t nparticles;
  UT_hash_handle hh;
  int64_t patch_id;      // allow for largest possible number (reserving highest bit for flagging things)
  //int level;
} subcube_t;

typedef subcube_t* psubcube_t;

#define SIZEOF_SUBCUBE (sizeof(subcube_t))
#define SIZEOF_PSUBCUBE (sizeof(psubcube_t))
#define SIZEOF_PARTID (sizeof(uint64_t))
#define SIZEOF_PARTPTR (sizeof(partptr))

//typedef struct subcube_s{
  //UT_array *particles;
  //partarray_ll* next;
//}partarray_ll;


//extern UT_icd cubekey_icd;

//subcubes Hash-table manipulation
void table_add_subcube(subcube_t**,subcube_t*);
void table_find_subcube(subcube_t*, cubekey_t*, subcube_t**);
void table_del_subcube(subcube_t*, subcube_t*);
uint64_t table_get_num_subcubes(subcube_t*);
uint64_t table_get_overhead(subcube_t*);
//void table_iterate_subcube(subcube_t* sc_table, subcube_t** sc_aux, subcube_t** sc_tmp){
//HASH_ITER macro implements a for statement waiting for the sentences block. 
//We MUST implement table_iterate_subcube as a macro in subcube.h
#define table_iterate_subcube(sc_table, sc_aux, sc_tmp)			\
HASH_ITER(hh, sc_table, sc_aux, sc_tmp)


//subcube manipulation
void subcube_create(subcube_t**,cubekey_t);
void subcube_add_particle(subcube_t*, partptr*);
partptr* subcube_next_particle(subcube_t*, partptr*);
uint64_t subcube_get_num_particles(subcube_t*);
uint64_t subcube_get_num_stored_particles(subcube_t*);
void subcube_free_particles_array(subcube_t*);
void subcube_free(subcube_t**);


// Depends on the particle data-structure access-> Moved to generate_tree.c (calling source)
//int subcube_refine(subcube_t**, subcube_t*);

//inline int insert_particle(subcube_t*, 

//TO-DO
int table_free(subcube_t*);

#endif


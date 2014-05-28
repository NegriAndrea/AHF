/*
 * subcube.c
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


#include "subcube.h"

//Structure needed by UTARRAY to create a particles array
UT_icd cubekey_icd = {SIZEOF_PARTPTR, NULL, NULL, NULL };

void subcube_create(subcube_t** sc,cubekey_t cubekey){
  if((*sc=(subcube_t*)malloc(sizeof(subcube_t)))==NULL){
    perror("malloc(subcube_t): ");
    fprintf(stderr,"[subcube_create] malloc for new subcube structure failed\n");
    return;
  }
  memcpy(&(*sc)->cubekey, &cubekey, SIZEOF_CUBEKEY);
  (*sc)->patch_id=NO_PATCH;
  (*sc)->nparticles=0;
  utarray_new((*sc)->particles, &cubekey_icd);
}

void table_add_subcube(subcube_t** sc_table, subcube_t* sc){
  HASH_ADD(hh,*sc_table,cubekey,sizeof(cubekey_t),sc);
}

void table_find_subcube(subcube_t* sc_table, cubekey_t* cubekey, subcube_t** sc){
  HASH_FIND(hh, sc_table, cubekey, sizeof(cubekey_t), *sc);
}

void table_del_subcube(subcube_t* sc_table, subcube_t* sc){
  HASH_DEL(sc_table,sc);
}

uint64_t table_get_num_subcubes(subcube_t* sc_table){
  return (uint64_t) HASH_COUNT(sc_table);
}

uint64_t table_get_overhead(subcube_t* sc_table){
  return (uint64_t) HASH_OVERHEAD(hh, sc_table);
}

//HASH_ITER macro implements a for statement waiting for the sentences block. 
//We MUST implement table_iterate_subcube as a macro in subcube.h
//void table_iterate_subcube(subcube_t* sc_table, subcube_t** sc_aux, subcube_t** sc_tmp){
  //HASH_ITER(hh, sc_table, *sc_aux, *sc_tmp);
//}

void subcube_add_particle(subcube_t* sc, partptr* ipart){
  utarray_push_back(sc->particles,ipart);
  sc->nparticles++;
}

partptr* subcube_next_particle(subcube_t* sc, partptr* part_it){
  if(sc->particles==NULL)
    return NULL;
  return (partptr*) utarray_next(sc->particles,part_it);
}

uint64_t subcube_get_num_particles(subcube_t* sc){
  return sc->nparticles;
}

uint64_t subcube_get_num_stored_particles(subcube_t* sc){
  if(sc->particles==NULL) return 0;
  return utarray_len(sc->particles);
}

void subcube_free_particles_array(subcube_t* sc){
  if(sc->particles){
    utarray_free(sc->particles);
    sc->particles=NULL;
  }
}
void subcube_free(subcube_t** sc){
  if(*sc!=NULL){
    subcube_free_particles_array(*sc);
    free(*sc);
    //We cannot set this pointer to NULL because it could be used as aux pointer in UThash iteration
    //*sc=NULL; 
  }
}

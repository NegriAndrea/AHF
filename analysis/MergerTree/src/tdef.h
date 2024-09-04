#ifndef INCLUDE_TDEF_H
#define INCLUDE_TDEF_H

typedef struct MTREE *MTREEptr;
typedef struct MTREE
{
  uint64_t haloid[2];
  uint64_t id[2];
  uint64_t npart[2];
  uint64_t common;
  double   merit;
} MTREE;

typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  uint64_t  haloid;
  uint64_t  npart;       // CAREFUL: we restrict the tree building to USE_PTYPE and hence this npart is not necessarily what is given in _halos as the total number of particles!
  uint64_t *Pid;
  
  uint64_t  ncroco;
  MTREEptr  mtree;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  uint64_t  nhalos;
  uint64_t *Hid;
}PARTS;

typedef char   boolean;          /* C has no boolean data type */
#endif

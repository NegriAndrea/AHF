#ifndef _HGN_RGADGET_H
#define _HGN_RGADGET_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint64_t PartId;	/* type of the particle id */

/*! The basica particle structure returned by routines */
typedef struct particle {
	PartId id;		/*!< the particles ID  */
	double pos[3];		/*!< the position x,y,z  */
	double vel[3];		/*!< velocity vx,vy,vz  */
	double mass;		/*!< mass of this particle */
} Particle;

/*! Descriptor for a reduced gadget file */
typedef struct RGfile {
	int fd;			/*!< file descriptor */
	const void *addr;	/*!< mapped address */
	PartId n_ids;		/*!< number of ids in the file */
	const PartId *idstart;	/*!< the list of IDs  */
	const Particle *pstart;	/*!< the list of particles */
} RGfile;

/*! Open a reduced gadget file
 * \return NULL on failure
 * \return structure on success
 */
RGfile *openrgadget(const char *filename);
/*! Fetch an individual particle
 * \return NULL on failure
 * \return structure on success
 * \note the return is read only memory
 */
const Particle *getparticle(RGfile *rg, PartId id);
/*! Fetch a group of particles
 * \param[in] rg - the file structure returned from openrgadget
 * \param[in] idlist - an arrary of particle IDs
 * \param[in] idcount - how many particles are in the array
 * \param[in] destination - an array of Particle pointers allocated by
 * the caller and filled in with the particle, or NULL if not in the file.
 * \return the number of found particles
 */

int getparticleset(RGfile *rg, PartId idlist[], 
		   unsigned int idcount, 
		   const Particle *destination[]);
/*! close a reduced gadget file and free off memory and things
 */
void closergadget(RGfile *rg);
#ifdef __cplusplus
}
#endif
#endif

//=================================================
// here the mode of operation of MergerTree is set
//=================================================

#ifndef INCLUDE_DEFINE_H
#define INCLUDE_DEFINE_H

//#define WITH_OPENMP                   // activate the (few) OpenMP parallel for-loops

#define USE_PTYPE                     // restrict analysis to certain particle species (define them yourself in read_particles()!)
                                      // note: when restricting to stars (or gas) it will only work with -DUSE_PIDMAP!
                                      // TODO: much more user-friendly implementation, please!
#define MINCOMMON      10             // we only cross-correlate haloes if they at least share MINCOMMON particles
#define MTREE_BOTH_WAYS               // make sure that every halo has only one descendant (this is relevant for SAM models)
//#define USE_PIDMAP                    // the Pids are not used as array indices anymore, e.g. one can now also use star particles for the tree building

// flags that only work together with SNAPSKIPPING
#define SNAPSKIPPING                         // whenever a connection [0]->[1] is not considered credible, the halo will be copied and considered in the connection [1]->[2] (recursively)
//#define SNAPSKIPPING_UNCREDIBLEMASSRATIO 2   // do not allow for mass jumps Mdesc = UNCREDIBLEMASSRATIO * Mprog
//#define SNAPSKIPPING_CONSIDERALLPROGENITORS  // consider all progenitors for UNCREDIBLEMASSRATIO test, in order of merit function
//#define DEBUG_SNAPSKIPPING


// support for AHF's MPI output
//#define READ_MPARTICLES               // support to read multiple _particle files, they must be of the latest AHF _particles file format!
#define NDIGITS                 4     // number of digits to be used for fileid


// support of ancient formats (without nhalos line, without haloids, whatever)
//#define THERE_IS_NO_NHALOS_LINE
//#define READ_HALOIDS_FROM_FILE
//#define USE_LINENUMBER_AS_HALOID      // overwrites(!) haloid as found in _particles
//#define EXCLUSIVE_PARTICLES           // each particle is only allowed to belong to one object (i.e. the lowest mass one): NOT fully tested yet!!!
//#define WITH_QSORT                    // uses qsort() instead of indexx() when ordering the progenitors according to merit function: NOT fully implemented yet!!!

// only relevant for read_particles_bin() (which is not being used anywhere at the moment)
#define SWAPBYTES 0

// GADGET convention for various particle types
#define PGAS                0.0   /* identifier for gas particles; has to be exactly 0.0!!!!                                     */
#define PDM                -1.0   /* identifier for dm particles; whatever negative value                                        */
#define PSTAR              -4.0   /* identifier for star particles; whatever negative value                                      */
#define PDMbndry           -5.0   /* identifier for bndry particles; whatever negative value                                     */

/*-------------------------------------------------------------------------------------
 *                                   DEPENDENCIES
 *-------------------------------------------------------------------------------------*/
#if (defined EXCLUSIVE_PARTICLES && defined WITH_OPENMP)
#define WITH_OPENMP_EXCLUSIVE
#endif

//#define MTREEDEVEL

// here we keep track of the version number manually
#define mtree_version 1.3

#endif


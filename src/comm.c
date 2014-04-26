/**
 * \file comm.c
 *
 * Provides functionality for communication between processes.
 */

/* Only needed in MPI mode */
#ifdef WITH_MPI

/**********************************************************************\
 *    Includes                                                        *
\**********************************************************************/
#include "comm.h"
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "libutility/alloc_struct.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               *
\**********************************************************************/

/**
 * Gives the largest number of sub chunks a large send or recieve is
 * splitted into.  The total number of bytes that can be send is thus
 * LOCAL_MAX_SUBSENDRECVS * LOCAL_MAX_MSGSIZE.
 */
#define LOCAL_MAX_SUBSENDRECVS 10

/**
 * This is the size of the largest MPI chunk to be send around (given in
 * bytes), if the message is larger, it will be broken up into smaller
 * sends/recv.
 * Note that 2147483647L corresponds to INT_MAX, the largest value you
 * can safely pass to the MPI functions as a count argument.
 */
#define LOCAL_MAX_MSGSIZE 2147483647L

/**
 * \brief Structure holding the information of what needs to be
 *        transfered.
 */
struct comm_struct {
	/** The number of elements to recieve */
	uint64_t recv;
	/** The offset from the beginning of the recieve array */
	uint64_t recv_displ;
	/** The number of elements to send */
	uint64_t send;
	/** The offset from the beginning of the send array */
	uint64_t send_displ;
	/** Request handler for the recieving */
	MPI_Request recv_req[LOCAL_MAX_SUBSENDRECVS];
	/** Request handler for the sending */
	MPI_Request send_req[LOCAL_MAX_SUBSENDRECVS];
};

/** Convenient typedef */
typedef struct comm_struct comm_struct_t;

/** Convenient typedef */
typedef comm_struct_t *comm_t;


/**
 * \brief Buffer element for the boundary update.
 */
struct comm_bound_struct {
	/** The key value */
	sfc_key_t key;
	/** The value that should be transfered */
	fpv_t val;
};

/** Convenient typedef */
typedef struct comm_bound_struct comm_bound_struct_t;

/** Convenient typedef */
typedef comm_bound_struct_t *comm_bound_t;


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/
#ifdef MPI_DEBUG
/**
 * \brief Prints debug information, only available in MPI_DEBUG
 *         mode.
 */
void
local_log_commarr(io_logging_t log, comm_t comm, int ncpu);
#endif


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern void
comm_dist_part(io_logging_t log,
               partptr *fst_part,
               uint64_t *no_part,
               loadbalance_t lb)
{
	comm_struct_t *comm;
	uint64_t totsend, totrecv;
	uint64_t tmp1, tmp2;
	partptr newpart;
	int i;
	uint64_t k;
	sfc_key_t target, j;
	MPI_Status status;

	/* Get memory for recv count and senddispl array */
	comm = (comm_struct_t *)malloc(sizeof(comm_struct_t)*lb->ncpu);
	if (comm == NULL) {
		io_logging_memfatal(log, "communication array");
		common_terminate(EXIT_FAILURE);
	}
	/* Initialize the requests (required later for usage checking) */
	for (i=0; i<lb->ncpu; i++) {
		int j;
		for (j=0; j<LOCAL_MAX_SUBSENDRECVS; j++) {
			comm[i].recv_req[j] = MPI_REQUEST_NULL;
			comm[i].send_req[j] = MPI_REQUEST_NULL;
		}
	}

	/* Telling everyone how many particles are going to be send */
	totrecv = UINT64_C(0);
	totsend = UINT64_C(0);
	for (i=0; i<lb->ncpu; i++) {
		MPI_Scatter((void *)(lb->no_parts_loc), 1, global_mpi.dt_uint64,
		            (void *)&(comm[i].recv), 1, global_mpi.dt_uint64,
		            i, MPI_COMM_WORLD);
		totrecv += comm[i].recv;
		if (i>0)
			comm[i].recv_displ = comm[i-1].recv_displ + comm[i-1].recv;
		else
			comm[i].recv_displ = UINT64_C(0);
		comm[i].send = lb->no_parts_loc[i];
		totsend += comm[i].send;
		if (i>0)
			comm[i].send_displ = comm[i-1].send_displ + comm[i-1].send;
		else
			comm[i].send_displ = UINT64_C(0);
	}

#	ifdef MPI_DEBUG
	local_log_commarr(global_io.log, comm, lb->ncpu);
#	endif

	/* Sanity check, this MUST be equal */
	if (totrecv != lb->no_parts[global_mpi.rank]) {
		io_logging_fatal(log,
		                 "Sanity check failed. The computed amount "
		                 "of particles I get from others plus the "
		                 "ones I keep do not match the amount of "
		                 "particles the loadbalancing scheme has "
		                 "determined.");
		io_logging_fatal(log,
		                 "totrecv = %" PRIu64 "\t"
		                 "lb->no_parts[%i] = %" PRIu64,
		                 totrecv, global_mpi.rank,
		                 lb->no_parts[global_mpi.rank]);
		free(comm);
		common_terminate(EXIT_FAILURE);
	}

	/* Get memory for the new particle array */
	newpart = c_part((long)totrecv);
	if (newpart == NULL) {
		io_logging_memfatal(log, "temporary particle storage (%fkB)",
		                    totrecv*sizeof(part));
		free(comm);
		common_terminate(EXIT_FAILURE);
	}

	/*
	 * Fire up the receiving and sending.
	 * Note that the receiving will FAIL if the message cannot be
	 * reveiced in at most LOCAL_MAX_SUBSENDRECVS chunks of size
	 * LOCAL_MAX_MSGSIZE (in bytes).  The offending receive will trigger
	 * the abort NOT then send.
	 */
	fflush(NULL);
	for (i=0; i<lb->ncpu; i++) {
		int j = 0;

		/* Chunked receiving */
		if (comm[i].recv > UINT64_C(0)) {
			uint64_t numPartRecvTodo = comm[i].recv;
			uint64_t numPartRecvDone = 0;
			uint64_t numPartRecvNow = 0;
			j=0;
			do {
				if (numPartRecvTodo * sizeof(part) > LOCAL_MAX_MSGSIZE) {
					numPartRecvNow = LOCAL_MAX_MSGSIZE / (sizeof(part));
				} else {
					numPartRecvNow = numPartRecvTodo;
				}
				MPI_Irecv((void *)(newpart+comm[i].recv_displ+numPartRecvDone),
				          numPartRecvNow * sizeof(part), MPI_BYTE,
				          i, j, MPI_COMM_WORLD, &(comm[i].recv_req[j]));
				j++;
				numPartRecvTodo -= numPartRecvNow;
				numPartRecvDone += numPartRecvNow;
			} while (j < LOCAL_MAX_SUBSENDRECVS && numPartRecvTodo > 0);
			/* Kill the program if not all particles could be received */
			if (j == LOCAL_MAX_SUBSENDRECVS && numPartRecvTodo > 0) {
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				io_logging_fatal(log,
				                 "Process %i could not receive %"PRIu64
				                 " particles from process %i in at most "
				                 "%i chunks. Increase LOCAL_MAX_MSGSIZE "
				                 "or LOCAL_MAX_SUBSENDRECVS in %s",
				                 rank, comm[i].recv, i,
				                 LOCAL_MAX_SUBSENDRECVS, __FILE__);
				common_terminate(EXIT_FAILURE);
			}
		}
		/* Chunked sending */
		if (comm[i].send > UINT64_C(0)) {
			uint64_t numPartSendTodo = comm[i].send;
			uint64_t numPartSendDone = 0;
			uint64_t numPartSendNow = 0;
			j = 0;
			do {
				if (numPartSendTodo * sizeof(part) > LOCAL_MAX_MSGSIZE) {
					numPartSendNow = LOCAL_MAX_MSGSIZE / (sizeof(part));
				} else {
					numPartSendNow = numPartSendTodo;
				}
				MPI_Isend((void *)((*fst_part)+comm[i].send_displ+numPartSendDone),
				          numPartSendNow * sizeof(part), MPI_BYTE,
				          i, j, MPI_COMM_WORLD, &(comm[i].send_req[j]));
				j++;
				numPartSendTodo -= numPartSendNow;
				numPartSendDone += numPartSendNow;
			} while (j < LOCAL_MAX_SUBSENDRECVS && numPartSendTodo > 0);
		}
	}

	/* Now sit there and wait for everything to finish */
	for (i=0; i<lb->ncpu; i++) {
		int j =0;
		if (comm[i].recv > UINT64_C(0)) {
			do {
				MPI_Wait(&(comm[i].recv_req[j]), &status);
				j++;
			} while (j<LOCAL_MAX_SUBSENDRECVS && comm[i].recv_req[j] != MPI_REQUEST_NULL);
		}
		j = 0;
		if (comm[i].send > UINT64_C(0)) {
			do {
				MPI_Wait(&(comm[i].send_req[j]), &status);
				j++;
			} while (j<LOCAL_MAX_SUBSENDRECVS && comm[i].send_req[j] != MPI_REQUEST_NULL);
		}
	}

	/* Throw away the old particle array */
	free(*fst_part);

	/* Get a new one */
	*fst_part = c_part((long)totrecv);
	if (*fst_part == NULL) {
		io_logging_memfatal(log, "new particle storage (%fkB)",
		                    totrecv*sizeof(part));
		free(newpart);
		free(comm);
		common_terminate(EXIT_FAILURE);
	}
	*no_part = totrecv;

	/* Prepare for sorting the particles back */
	tmp1 = tmp2 = UINT64_C(0);
	for (j=lb->fstkey[global_mpi.rank];
	     j<=lb->lstkey[global_mpi.rank];
	     j++) {
		tmp1 = lb->bf[j];
		lb->bf[j] = tmp2;
		tmp2 += tmp1;
	}

	/* Start sorting back */
	for (k=0; k<totrecv; k++) {
		target = sfc_curve_contract(lb->level, BITS_PER_DIMENSION,
		                            lb->ctype, newpart[k].sfckey);
		memcpy((void *)(*fst_part+lb->bf[target]),
		       (void *)(newpart+k),
		       sizeof(part));
		(newpart+k)->ll = NULL;
		lb->bf[target]++;
	}

	/* Clean up */
	free(newpart);
	free(comm);

	/* Refill bf properly */
	for (j=lb->lstkey[global_mpi.rank];
	     j>lb->fstkey[global_mpi.rank];
	     j--) {
		lb->bf[j] -= lb->bf[j-1];
	}

	return;
}

extern uint64_t
comm_dist_part_ahf(io_logging_t log,
                   partptr *fst_part,
                   uint64_t *no_part,
                   loadbalance_t lb)
{
	comm_struct_t *comm;
	uint64_t totsend, totrecv, new_no_part;
	partptr newpart, *parts, part;
	int i;
	sfc_key_t redkey, j;
	MPI_Status status;

	/* Get memory for recv count and senddispl array */
	comm = (comm_struct_t *)malloc(sizeof(comm_struct_t)*lb->ncpu);
	if (comm == NULL) {
		io_logging_memfatal(log, "communication array");
		common_terminate(EXIT_FAILURE);
	}

	/* Figure out how many particles are in the boundaries and
	 * hence need to be send */
	for (i=0;i<lb->ncpu; i++) {
		lb->no_parts_loc[i] = UINT64_C(0);
		/* Loop over all cells in the inner boundary to CPU i
		 * and count the particles that need to be send there
		 */
		for (j=UINT64_C(0); j<lb->bound->inner[i].num; j++) {
			lb->no_parts_loc[i] += 
			  lb->bf[(uint64_t)(lb->bound->inner[i].bound[j])];
		}
	}

	/* Telling everyone how many particles are going to be send */
	totrecv = UINT64_C(0);
	totsend = UINT64_C(0);
	for (i=0; i<lb->ncpu; i++) {
		MPI_Scatter((void *)(lb->no_parts_loc), 1, global_mpi.dt_uint64,
		            (void *)&(comm[i].recv), 1, global_mpi.dt_uint64,
		            i, MPI_COMM_WORLD);
		totrecv += comm[i].recv;
		if (i>0)
			comm[i].recv_displ = comm[i-1].recv_displ + comm[i-1].recv;
		else
			comm[i].recv_displ = UINT64_C(0);
		comm[i].send = lb->no_parts_loc[i];
		totsend += comm[i].send;
		if (i>0)
			comm[i].send_displ = comm[i-1].send_displ + comm[i-1].send;
		else
			comm[i].send_displ = UINT64_C(0);
	}

#	ifdef MPI_DEBUG
	local_log_commarr(global_io.log, comm, lb->ncpu);
#	endif

	/* Enlarge the local particle storage */
	new_no_part = *no_part + totrecv;
	newpart = (partptr)realloc(*fst_part,
	                           sizeof(struct particle)*new_no_part);
	if (newpart == NULL) {
		io_logging_memfatal(log, "enlarged local particle storage");
		common_terminate(EXIT_FAILURE);
	}
	*fst_part = newpart;
	
	/* Generate storage for the duplicated particles */
	newpart = (partptr)malloc(sizeof(struct particle)*totsend);
	if (newpart == NULL) {
		io_logging_memfatal(log, "temporary particle send buffer");
		common_terminate(EXIT_FAILURE);
	}

	/* Fill the buffer */
	parts = (partptr *)malloc(sizeof(partptr)*lb->ncpu);
	if (parts == NULL) {
		io_logging_memfatal(log, "particle pointer array");
	}
	parts[0] = newpart;
	for (i=1;i<lb->ncpu;i++) {
		parts[i] = parts[i-1] + comm[i-1].send;
	}
	part = *fst_part;
	for (j=0; j<*no_part; j++, part++) {
    //		redkey = sfc_curve_contract(LOADBALANCE_DOMAIN_LEVEL,
    redkey = sfc_curve_contract(global_io.params->lb_level,
		                            BITS_PER_DIMENSION,
		                            lb->ctype,
		                            part->sfckey);
		for (i=0; i<lb->ncpu; i++) {
			if (sfc_boundary_findKeyPos(log, &(lb->bound->inner[i]),
			                            redkey, NULL)) {
				(void *)memcpy((void *)(parts[i]),
				               (void *)part,
				               sizeof(struct particle));
				parts[i]++;
			}
		}
	}
	free(parts);
	
	/* Fire up the recieving and sending */
	for (i=0; i<lb->ncpu; i++) {
		if (comm[i].recv > UINT64_C(0)) {
			if (comm[i].send*sizeof(struct particle) > LOCAL_MAX_MSGSIZE) {
				io_logging_fatal(log,
				                 "Trying to receive too many particles "
				                 "from process %i",
				                 i);
				common_terminate(EXIT_FAILURE);
			}
			MPI_Irecv((void *)(*fst_part+*no_part+comm[i].recv_displ),
			          comm[i].recv*sizeof(struct particle),
			          MPI_BYTE,
			          i, 0, MPI_COMM_WORLD, &(comm[i].recv_req[0]));
		}
		if (comm[i].send > UINT64_C(0)) {
			MPI_Isend((void *)(newpart+comm[i].send_displ),
			          comm[i].send*sizeof(struct particle),
			          MPI_BYTE,
			          i, 0, MPI_COMM_WORLD, &(comm[i].send_req[0]));
		}
	}

	/* Now sit there and wait for everything to finish */
	for (i=0; i<lb->ncpu; i++) {
		if (comm[i].recv > UINT64_C(0))
			MPI_Wait(&(comm[i].recv_req[0]), &status);
		if (comm[i].send > UINT64_C(0))
			MPI_Wait(&(comm[i].send_req[0]), &status);
	}

	/* Clean up */
	free(newpart);
	free(comm);
	*no_part = new_no_part;

	return totrecv;
}

extern void
comm_update_bound(io_logging_t log,
                  loadbalance_t lb,
                  gridls *curgrid,
                  int offset)
{
	nptr node;
	comm_t comm;
	/* TODO This works more efficient as the information where to put
	 *      the data is stored in the bound-structure attached to the
	 *      loadbalacing module.
	 */
	comm_bound_t send_buff;
	comm_bound_t recv_buff;
	int j;
	uint64_t i, totrecv;
	MPI_Status status;
	double total_masstransfer;

	/* Get memory for the communication array */
	comm = (comm_t)malloc(sizeof(comm_struct_t)*lb->ncpu);
	if (comm == NULL) {
		io_logging_memfatal(log, "communication array");
		common_terminate(EXIT_FAILURE);
	}

	/* Initialise the communication array with all known or safe
	 * initial values
	 */
	comm[0].recv = (uint64_t)lb->bound->inner[0].num;
	comm[0].send = UINT64_C(0);
	totrecv = comm[0].recv;
	comm[0].recv_displ = UINT64_C(0);
	comm[0].send_displ = UINT64_C(0);
	for (j=1; j<lb->ncpu; j++) {
		comm[j].recv = (uint64_t)lb->bound->inner[j].num;
		totrecv += comm[j].recv;
		comm[j].recv_displ =   comm[j-1].recv_displ
		                     + comm[j-1].recv;
		comm[j].send = UINT64_C(0);
		comm[j].send_displ = UINT64_C(0);
	}

	/* Generate the buffer for the sending */
	send_buff = (comm_bound_t)malloc( sizeof(comm_bound_struct_t)
	                                 *lb->bound->outer->num);
	if (send_buff == NULL) {
		io_logging_memfatal(log,
		                    "buffer for sending boundary (%fkB)",
		                     sizeof(comm_bound_struct_t)
		                    *lb->bound->outer->num / 1024.);
		common_terminate(EXIT_FAILURE);
	}

	/* Generate the buffer for the recieving */
	recv_buff = (comm_bound_t)malloc( sizeof(comm_bound_struct_t)
	                                 *totrecv);
	if (recv_buff == NULL) {
		io_logging_memfatal(log,
		                    "buffer for recieving boundary (%fkB)",
		                     sizeof(comm_bound_struct_t)
		                    *totrecv / 1024.);
		common_terminate(EXIT_FAILURE);
	}

	/* Quickly fire up the recieving */
	for (j=0; j<lb->ncpu; j++) {
		if (comm[j].recv > UINT64_C(0)) {
			MPI_Irecv((void *)(recv_buff+comm[j].recv_displ),
			          comm[j].recv*sizeof(comm_bound_struct_t),
			          MPI_BYTE,
			          j, 0, MPI_COMM_WORLD, &(comm[j].recv_req[0]));
# ifdef MPI_DEBUG
			io_logging_msg(global_io.log, 6,
			               "Started to listen for %" PRIu64
			               " elementes from %i",
			               comm[j].recv, j);
# endif
		}
	}

	/* Fill the sending buffer */
	j = 0;
	total_masstransfer = 0.0;
	for (i=0; i<lb->bound->outer->num; i++) {
		/* Figure out to which CPU we have to send this cell */
		while (lb->lstkey[j] < lb->bound->outer->bound[i])
			j++;

		/* Now get the value from the node boundary */
		node = curgrid->bound+i;
		(send_buff+i)->key = lb->bound->outer->bound[i];
		(send_buff+i)->val = *((fpv_t *)(((char *)node)+offset));
		total_masstransfer += (send_buff+i)->val;
		
		/* Count this cell */
		comm[j].send++;
	}
#ifdef MPI_DEBUG
	io_logging_msg(global_io.log, 6,
	               "Total mass transfer (sending): %g",
	               total_masstransfer);
#endif

	/* Fill in the rest of the communication structure and start
	 * sending.
	 */
	comm[0].recv = (uint64_t)lb->bound->inner[0].num;
	comm[0].recv_displ = UINT64_C(0);
	comm[0].send_displ = UINT64_C(0);
	if (comm[0].send > UINT64_C(0))
		MPI_Isend((void *)(send_buff+comm[0].send_displ),
		          comm[0].send*sizeof(comm_bound_struct_t),
		          MPI_BYTE,
		          0, 0, MPI_COMM_WORLD, &(comm[0].send_req[0]));
	for (j=1; j<lb->ncpu; j++) {
		comm[j].send_displ =   comm[j-1].send_displ
		                     + comm[j-1].send;
		if (comm[j].send > UINT64_C(0))
			MPI_Isend((void *)(send_buff+comm[j].send_displ),
			          comm[j].send*sizeof(comm_bound_struct_t),
			          MPI_BYTE,
			          j, 0, MPI_COMM_WORLD, &(comm[j].send_req[0]));
	}
#ifdef MPI_DEBUG
	local_log_commarr(global_io.log, comm, lb->ncpu);
	for (j=0; j<lb->ncpu; j++) {
		io_logging_msg(global_io.log, INT32_C(20),
		               "Sending to %i:", (int)j);
		for (i=0; i<comm[j].send; i++) {
			io_logging_msg(global_io.log, INT32_C(20),
			               "  key: %" SFC_PRIkey "  val: %e",
			               (send_buff+i+comm[j].send_displ)->key,
			               (send_buff+i+comm[j].send_displ)->val);
		}
	}
#endif

	/* Sit there and wait for everything to finish */
	for (j=0; j<lb->ncpu; j++) {
		if (comm[j].recv > UINT64_C(0)) {
			MPI_Wait(&(comm[j].recv_req[0]), &status);
		}
		if (comm[j].send > UINT64_C(0)) {
			MPI_Wait(&(comm[j].send_req[0]), &status);
		}
	}

	total_masstransfer = 0.0;
	/* Copy the boundary information back to the grid */
	for (i=0; i<totrecv; i++) {
		/* Now locate the node in the grid */
		node = get_node_from_key(curgrid,
		                         recv_buff[i].key,
		                         lb->ctype,
		                         lb->level);
		*((fpv_t *)(((char *)node)+offset)) += recv_buff[i].val;
		total_masstransfer += recv_buff[i].val;
	}
#ifdef MPI_DEBUG
	io_logging_msg(global_io.log, 6,
	               "Total mass transfer (recieving): %g",
	               total_masstransfer);
#endif

	/* Clean up */
	free(send_buff);
	free(recv_buff);
	free(comm);

	return;
}

#ifdef WITH_FFTW
extern void
comm_fill_fft_slab(io_logging_t log,
                   gridls *curgrid,
                   fft_t fft,
                   loadbalance_t lb)
{
	int i;
	sfc_key_t key;
	nptr node;
	uint32_t pos[3];
	int ijk;

	for (i=0; i<fft->local_nx*fft->dom_grid*fft->dom_grid; i++)
		fft->data[i] = fft->work[i] = 0.0;

# ifdef DOIT
	/* Loop over the grid and copy the value to the FFT slab */
	for (key=lb->fstkey[global_mpi.rank];
	     key<=lb->lstkey[global_mpi.rank];
	     key++) {
		node = get_node_from_key(curgrid, key, lb->ctype, lb->level);
		sfc_curve_calcPos(lb->ctype, key, lb->level, pos);
		if (    (pos[0] >= fft->local_x_start)
		     && (pos[0] < fft->local_x_start + fft->local_nx)) {
			ijk = pos[2] + fft->dom_grid*(
			      fft->dom_grid*(pos[0]-fft->local_x_start) + pos[1]);
			fft->data[ijk] = node->dens;
		}
	}
# endif

	return;
}

extern void 
comm_read_fft_slab(io_logging_t log,
                   gridls *curgrid,
                   fft_t fft,
                   loadbalance_t lb)
{
	

	return;
}
#endif /* WITH_FFTW */


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
#ifdef MPI_DEBUG
void
local_log_commarr(io_logging_t log, comm_t comm, int ncpu)
{
	int i;

	for (i=0; i<ncpu; i++) {
		io_logging_msg(global_io.log, INT32_C(6),
		               "comm[%i].recv = %" PRIu64,
		               i, comm[i].recv);
		io_logging_msg(global_io.log, INT32_C(6),
		               "comm[%i].recv_displ = %" PRIu64,
		               i, comm[i].recv_displ);
		io_logging_msg(global_io.log, INT32_C(6),
		               "comm[%i].send = %" PRIu64,
		               i, comm[i].send);
		io_logging_msg(global_io.log, INT32_C(6),
		               "comm[%i].send_displ = %" PRIu64,
		               i, comm[i].send_displ);
	}

	return;
}
#endif

 #endif /* WITH_MPI */

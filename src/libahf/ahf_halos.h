#ifndef AHF_HALOS_H
#define AHF_HALOS_H

#  include <stddef.h>
#  include <stdlib.h>
#  include <stdio.h>
#  include <math.h>
#  include <string.h>
#  include "../tdef.h"

/* Those are the necessary conversion factors */
extern double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, Hubble;

int
HaloProfiles(HALO *);

void
rem_nothing(HALO *);

void
rem_unbound(HALO *);

void
rem_outsideRvir(HALO *, int);

void
merge_ll_and_ipart(HALO *);

void
sort_halo_particles(HALO *);


#  if (defined AHFrestart || defined WITH_MPI)
void
rem_boundary_haloes(void);
void
flag_boundary_haloes(void);

#  endif

#  ifdef AHFphspdens
/** Does the phase-space stuff */
int
HaloProfilesPhaseSpace(HALO *);
#  endif

#  ifdef AHFdisks
int
HaloProfilesDisk(HALO *);
#  endif
#endif /* AHF_HALOS_H */

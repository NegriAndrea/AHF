#ifndef AHF_HALOS_SFC_H
#define AHF_HALOS_SFC_H

#include "../tdef.h"

/**
 * \brief  This is the function that will build the haloes from the
 *         given centres and gather radii.
 *
 * \param  halo  The halo to build.
 *
 * \return  Returns nothing, however the external halo structure will be
 *          updated to reflect the changes done in this routine.
 */
void
ahf_halos_sfc_constructHalo(HALO *halo);

/**
 * \brief  This will gather all particles for a given halo which are
 *         inside its  gathering radius.
 *
 * \param  halo  The halo to work on.
 *
 * \return  Returns nothing, however the external halo structure will be
 *          updated to reflect the changes done in this routine.
 */
void
ahf_halos_sfc_gatherParts(HALO *halo);


#endif  /* AHF_HALOS_SFC_H */

#ifndef AHF_IO_SQL_H
#define AHF_IO_SQL_H

/*--- Includes ----------------------------------------------------------*/
#include <stdint.h>
#include "../tdef.h"


/*--- ADT handle --------------------------------------------------------*/
typedef struct ahf_io_sql_struct *ahf_io_sql_t;


/*--- Exported Types ----------------------------------------------------*/
typedef enum {
	AHF_IO_SQL_MODE_SINGLE,
	AHF_IO_SQL_MODE_MULTIPLE
} ahf_io_sql_mode_t;


/*--- Prototypes of exported functions ----------------------------------*/
extern ahf_io_sql_t
ahf_io_sql_new(const char        *dbFileNamePrefix,
               ahf_io_sql_mode_t mode,
               uint64_t          haloIDOffset);

extern void
ahf_io_sql_del(ahf_io_sql_t *sql);

extern void
ahf_io_sql_createTables(ahf_io_sql_t sql);

extern void
ahf_io_sql_writeHalos(ahf_io_sql_t  sql,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos);

extern void
ahf_io_sql_writeProfiles(ahf_io_sql_t  sql,
                         HALO          *halos,
                         unsigned long *idx,
                         int           numHalos);

extern void
ahf_io_sql_writeParticles(ahf_io_sql_t  sql,
                          HALO          *halos,
                          unsigned long *idx,
                          int           numHalos);

#endif

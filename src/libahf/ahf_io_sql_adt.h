#ifndef AHF_IO_SQL_ADT_H
#define AHF_IO_SQL_ADT_H


/*--- Includes ----------------------------------------------------------*/
#include <stdint.h>
#include <sqlite3.h>


/*--- Useful typedesfs --------------------------------------------------*/
typedef sqlite3      *sqlite3_t;
typedef sqlite3_stmt *sqlite3_stmt_t;


/*--- ADT implementation ------------------------------------------------*/
struct ahf_io_sql_struct {
	uint64_t          haloIDOffset;
	ahf_io_sql_mode_t mode;
	char              *dbFileNameHalos;
	char              *dbFileNameProfiles;
	char              *dbFileNameParticles;
	sqlite3_t         dbHalos;
	sqlite3_t         dbProfiles;
	sqlite3_t         dbParticles;
};


#endif

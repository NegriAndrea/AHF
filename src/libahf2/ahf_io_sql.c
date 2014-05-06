#ifdef AHF2

#ifdef AHF_SQL
/*--- Includes ----------------------------------------------------------*/
#define _POSIX_C_SOURCE 199309L

/*************************************************************************\
 * WARNING  The ORDER of the headers is IMPORTANT here:
 * WARNING    sqlite3.h MUST be included before ../tdef.h
 * WARNING    As ahf_io_sql.h does include ../tdef.h, sqlite3.h has been
 * WARNING    moved in front of it.
 * WARNING
 * WARNING  The reason for this messup is unclear at the moment, most
 * WARNING  likely tdef.h defines something that confuses sqlite3.h.
 * WARNING  I will investigate eventually what is going wrong there,
 * WARNING  but for the moment lets be happy with fixing the order of
 * WARNING  the headers.
 * WARNING
 *************************************************************************/
#include <sqlite3.h>
#include "ahf_io_sql.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>
#include "../tdef.h"
#include "../define.h"
#include "../param.h"
#include "../common.h"
#include "../libutility/cosmology.h"
#include "../libutility/utility.h"


/*--- ADT implementation ------------------------------------------------*/
#include "ahf_io_sql_adt.h"


/*--- Dealing with global variables -------------------------------------*/
extern double r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, u_fac, Hubble;


/*--- Local defines -----------------------------------------------------*/
#define LOCAL_SNOOZETIME_IN_NANOSECONDS 25000

#define diediedie(errCode)                                \
	{                                                     \
		fprintf(stderr,                                   \
		        "FATAL:  Deathtrap in %s() at %s:%i\n"    \
		        "Terminating with error code: %i\n\n",    \
		        __func__, __FILE__, __LINE__, (errCode)); \
		exit(errCode);                                    \
	}

#define local_dealWithSQLiteOpenV2Error(returnCode, db)      \
	{                                                        \
		if (returnCode != SQLITE_OK)                         \
		{                                                    \
			fprintf(stderr, "sqlite3_open_v2 failed:  %s\n", \
			        sqlite3_errmsg(db));                     \
			diediedie(EXIT_FAILURE);                         \
		}                                                    \
	}

#define local_dealWithSQLiteExecError(returnCode, message)                \
	{                                                                     \
		if (returnCode != SQLITE_OK)                                      \
		{                                                                 \
			if (message != NULL)                                          \
				fprintf(stderr, "sqlite3_exec failed:  %s\n",             \
				        message);                                         \
			else                                                          \
				fprintf(stderr, "sqlite3_exec failed:  Unknown Error\n"); \
			diediedie(EXIT_FAILURE);                                      \
		}                                                                 \
	}

#define LOCAL_WRITEHALOS_STATEMENT_SQL \
    "INSERT INTO halos VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
#define LOCAL_WRITEPROFILES_STATEMENT_SQL \
    "INSERT INTO profiles VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
#define LOCAL_WRITEPARTICLES_STATEMENT_SQL \
    "INSERT INTO particles VALUES(?,?)"
#define LOCAL_TABLE_HALOS_NAME     "halos"
#define LOCAL_TABLE_PROFILES_NAME  "profiles"
#define LOCAL_TABLE_PARTICLES_NAME "particles"


/*--- Prototypes of local functions -------------------------------------*/
static char *
local_getDBFileName(const char *prefix, const char *qualifier);

static sqlite3_t
local_openConnection(const char *dbFileName);

static void
local_tweakConnection(sqlite3_t db);

static int
local_waitForDBAccess(void *unused, int numCalled);

static void
local_beginTransaction(sqlite3_t db);

static void
local_endTransaction(sqlite3_t db);

static sqlite3_stmt_t
local_getPreparedStatement(sqlite3_t db, const char *stmtSQL);

inline static int64_t
local_getHaloIDFromIdx(int idx, uint64_t offset);

static void
local_createTableHalos(sqlite3_t db);

static void
local_createTableProfiles(sqlite3_t db);

static void
local_createTableParticles(sqlite3_t db);


/*--- Implementations of exported functions -----------------------------*/
extern ahf_io_sql_t
ahf_io_sql_new(const char        *dbFileNamePrefix,
               ahf_io_sql_mode_t mode,
               uint64_t          haloIDOffset)
{
	ahf_io_sql_t sql;

	assert(dbFileNamePrefix != NULL);

	sql = malloc(sizeof(struct ahf_io_sql_struct));
	if (sql == NULL) {
		fprintf(stderr, "Could not allocate memory for sql structure\n");
		exit(EXIT_FAILURE);
	}
	sql->haloIDOffset = haloIDOffset;
	sql->mode         = mode;
	if (sql->mode == AHF_IO_SQL_MODE_MULTIPLE) {
		sql->dbFileNameHalos     = local_getDBFileName(dbFileNamePrefix,
		                                               ".halos");
		sql->dbFileNameProfiles  = local_getDBFileName(dbFileNamePrefix,
		                                               ".profiles");
		sql->dbFileNameParticles = local_getDBFileName(dbFileNamePrefix,
		                                               ".particles");
	} else {
		sql->dbFileNameHalos     = local_getDBFileName(dbFileNamePrefix,
		                                               NULL);
		sql->dbFileNameProfiles  = sql->dbFileNameHalos;
		sql->dbFileNameParticles = sql->dbFileNameHalos;
	}
	sql->dbHalos = local_openConnection(sql->dbFileNameHalos);
	local_tweakConnection(sql->dbHalos);
	if (sql->mode == AHF_IO_SQL_MODE_MULTIPLE) {
		sql->dbProfiles  = local_openConnection(sql->dbFileNameProfiles);
		sql->dbParticles = local_openConnection(sql->dbFileNameParticles);
		local_tweakConnection(sql->dbParticles);
		local_tweakConnection(sql->dbProfiles);
	} else {
		sql->dbProfiles  = sql->dbHalos;
		sql->dbParticles = sql->dbHalos;
	}

	return sql;
} /* ahf_io_sql_new */

extern void
ahf_io_sql_del(ahf_io_sql_t *sql)
{
	assert(sql != NULL && *sql != NULL);

	free((*sql)->dbFileNameHalos);
	sqlite3_close((*sql)->dbHalos);
	if ((*sql)->mode == AHF_IO_SQL_MODE_MULTIPLE) {
		free((*sql)->dbFileNameProfiles);
		sqlite3_close((*sql)->dbProfiles);
		free((*sql)->dbFileNameParticles);
		sqlite3_close((*sql)->dbParticles);
	}
	free(*sql);

	*sql = NULL;
}

extern void
ahf_io_sql_createTables(ahf_io_sql_t sql)
{
	assert(sql != NULL);

	local_beginTransaction(sql->dbHalos);
	{
		local_createTableHalos(sql->dbHalos);
	}
	local_endTransaction(sql->dbHalos);
	local_beginTransaction(sql->dbProfiles);
	{
		local_createTableProfiles(sql->dbProfiles);
	}
	local_endTransaction(sql->dbProfiles);
	local_beginTransaction(sql->dbParticles);
	{
		local_createTableParticles(sql->dbParticles);
	}
	local_endTransaction(sql->dbParticles);
}

extern void
ahf_io_sql_writeHalos(ahf_io_sql_t  sql,
                      HALO          *halos,
                      unsigned long *idx,
                      int           numHalos)
{
	assert(sql != NULL);
	if (numHalos <= 0)
		return;

	assert(halos != NULL);
	assert(idx != NULL);

	local_beginTransaction(sql->dbHalos);
	{
		sqlite3_stmt_t stmt;
		int            rc;
		int            snapID = -1;
		stmt = local_getPreparedStatement(sql->dbHalos,
		                                  LOCAL_WRITEHALOS_STATEMENT_SQL);
		for (int i = 0; i < numHalos; i++) {
			int     j      = idx[i];
			int64_t haloID = local_getHaloIDFromIdx(j, sql->haloIDOffset);

			if ((halos[j].npart < simu.AHF_MINPART))
				continue;

			sqlite3_bind_int64(stmt, 1, haloID);

			sqlite3_bind_int64(stmt, 2, halos[j].npart);
			sqlite3_bind_double(stmt, 3, halos[j].pos.x * x_fac);
			sqlite3_bind_double(stmt, 4, halos[j].pos.y * x_fac);
			sqlite3_bind_double(stmt, 5, halos[j].pos.z * x_fac);
			sqlite3_bind_double(stmt, 6, halos[j].vel.x * v_fac);
			sqlite3_bind_double(stmt, 7, halos[j].vel.y * v_fac);
			sqlite3_bind_double(stmt, 8, halos[j].vel.z * v_fac);


			rc = sqlite3_step(stmt);
			if (rc != SQLITE_DONE) {
				fprintf(stderr, "sqlite3_step failed:  %s\n",
				        sqlite3_errmsg(sql->dbHalos));
				exit(EXIT_FAILURE);
			}
			sqlite3_reset(stmt);
		}
		rc = sqlite3_finalize(stmt);
		if (rc != SQLITE_OK) {
			fprintf(stderr, "sqlite3_finalize failed: %s\n",
			        sqlite3_errmsg(sql->dbHalos));
			exit(EXIT_FAILURE);
		}
	}
	local_endTransaction(sql->dbHalos);
} /* ahf_io_sql_writeHalos */

extern void
ahf_io_sql_writeProfiles(ahf_io_sql_t  sql,
                         HALO          *halos,
                         unsigned long *idx,
                         int           numHalos)
{
	assert(sql != NULL);
	if (numHalos <= 0)
		return;

	assert(halos != NULL);
	assert(idx != NULL);

	local_beginTransaction(sql->dbProfiles);
	{
		sqlite3_stmt_t stmt;
		int            rc;
		int            snapID = -1;
		stmt = local_getPreparedStatement(sql->dbProfiles,
		                                  LOCAL_WRITEPROFILES_STATEMENT_SQL);
		for (int i = 0; i < numHalos; i++) {
			int     j = idx[i];
			int64_t haloID;
			double  rad;
			int     radiiConverged = 0;
			haloID = local_getHaloIDFromIdx(j, sql->haloIDOffset);

			for (int k = 0; k < halos[j].prof.nbins; k++) {
				double t_relax, age;
				/* check for converged radius (Power et al. 2003) */
				rad     = halos[j].prof.r[k];
				t_relax = halos[j].prof.npart[k]
				          / log(rad / halos[j].spaRes)
				          * rad / sqrt(halos[j].prof.v2_circ[k]);

				/* convert to (Mpc/h) / (km/sec) ~ (Mpc/h) / (kpc/Gyr) */
				t_relax *= r_fac / sqrt(phi_fac);

				/* convert to Gyr/h */
				t_relax *= 1E3;

				/* age of the Universe in Gyr/h */
				age = calc_t(global.a) * simu.t_unit * Mpc / Gyr;

				/* if not converged, write negative radius into
				 *.AHF_profiles */
				if (t_relax < 0.9 * age)
					/* The +1 is merely to keep the name consistent:
					 * r_conv_i should give the smallest converged radius,
					 * not the bin before it. Hence to check for converged
					 * bin or not, i<r_conv_i is the thing to do. */
					radiiConverged = k + 1;
			}
			for (int k = 0; k < halos[j].prof.nbins; k++) {
				double flip = k < radiiConverged ? -1.0 : 1.0;
				double rMin, rMax;

				if ((halos[j].npart < simu.AHF_MINPART))
					continue;

				sqlite3_bind_int64(stmt, 1, snapID);
				sqlite3_bind_int64(stmt, 2, haloID);
				rMin = (k > 0) ? flip * halos[j].prof.r[k - 1] : flip * 0.0;
				rMax = flip * halos[j].prof.r[k];
				sqlite3_bind_double(stmt, 3, rMin * x_fac);
				sqlite3_bind_double(stmt, 4, rMax * x_fac);


				//TODO: decide what to dump

				/*------------------------------------------------------------
				 * THIS IS THE PART WHERE YOU DECIDE WHAT INFORMATION TO DUMP
				 *
				 * VERY IMPORTANT THOUGH:
				 * make sure to synchronize this with:
				 *
				 *         local_createTableProfiles(sqlite3_t db)
				 *
				 *------------------------------------------------------------*/


				sqlite3_bind_double(stmt, 5, halos[j].prof.dens[k]   * rho_fac / global.rho_vir);
				sqlite3_bind_double(stmt, 6, halos[j].prof.ovdens[k] * rho_fac / global.rho_vir);


				rc = sqlite3_step(stmt);
				if (rc != SQLITE_DONE) {
					fprintf(stderr, "sqlite3_step failed:  %s\n",
					        sqlite3_errmsg(sql->dbProfiles));
					exit(EXIT_FAILURE);
				}
				sqlite3_reset(stmt);
			}
		}
		rc = sqlite3_finalize(stmt);
		if (rc != SQLITE_OK) {
			fprintf(stderr, "sqlite3_finalize failed: %s\n",
			        sqlite3_errmsg(sql->dbProfiles));
			exit(EXIT_FAILURE);
		}
	}
	local_endTransaction(sql->dbProfiles);
} /* ahf_io_sql_writeProfiles */

extern void
ahf_io_sql_writeParticles(ahf_io_sql_t  sql,
                          HALO          *halos,
                          unsigned long *idx,
                          int           numHalos)
{
	assert(sql != NULL);
	if (numHalos <= 0)
		return;

	assert(halos != NULL);
	assert(idx != NULL);

	local_beginTransaction(sql->dbParticles);
	{
		sqlite3_stmt_t stmt;
		int            rc;
		int            snapID = -1;
		stmt = local_getPreparedStatement(
		    sql->dbParticles,
		    LOCAL_WRITEPARTICLES_STATEMENT_SQL);
		for (int i = 0; i < numHalos; i++) {
			int     j      = idx[i];
			int64_t haloID = local_getHaloIDFromIdx(j, sql->haloIDOffset);

			if ((halos[j].npart < simu.AHF_MINPART))
				continue;

			for (int k = 0; k < halos[j].npart; k++) {
				partptr current = global.fst_part + halos[j].ipart[k];
				sqlite3_bind_int(stmt, 1, snapID);
				sqlite3_bind_int(stmt, 2, haloID);
				sqlite3_bind_int64(stmt, 3, current->id);
				rc = sqlite3_step(stmt);
				if (rc != SQLITE_DONE) {
					fprintf(stderr, "sqlite3_step failed:  %s\n",
					        sqlite3_errmsg(sql->dbParticles));
					exit(EXIT_FAILURE);
				}
				sqlite3_reset(stmt);
			}
		}

		rc = sqlite3_finalize(stmt);
		if (rc != SQLITE_OK) {
			fprintf(stderr, "sqlite3_finalize failed: %s\n",
			        sqlite3_errmsg(sql->dbParticles));
			exit(EXIT_FAILURE);
		}
	}
	local_endTransaction(sql->dbParticles);
} /* ahf_io_sql_writeParticles */

/*--- Implementations of local functions --------------------------------*/
static char *
local_getDBFileName(const char *prefix, const char *qualifier)
{
	char *dbFileName;
	int  lengthOfDBFileName;

	lengthOfDBFileName = 3 + strlen(prefix); // ".db" are the 3
	if (qualifier != NULL)
		lengthOfDBFileName += strlen(qualifier);

	dbFileName = malloc(sizeof(char) * (lengthOfDBFileName + 1));
	sprintf(dbFileName, "%s%s.db",
	        prefix,
	        qualifier != NULL ? qualifier : "");

	return dbFileName;
}

static sqlite3_t
local_openConnection(const char *dbFileName)
{
	int       rc;
	sqlite3_t db;

	rc = sqlite3_open_v2(dbFileName, &db,
	                     SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
	                     NULL);
	local_dealWithSQLiteOpenV2Error(rc, db);

	return db;
}

static void
local_tweakConnection(sqlite3_t db)
{
	int  rc;
	char *errorMessage = NULL;

	sqlite3_busy_handler(db, &local_waitForDBAccess, NULL);
	rc = sqlite3_exec(db, "PRAGMA synchronous=OFF",
	                  NULL, NULL, &errorMessage);
	local_dealWithSQLiteExecError(rc, errorMessage);
	rc = sqlite3_exec(db, "PRAGMA default_temp_store=2",
	                  NULL, NULL, &errorMessage);
	local_dealWithSQLiteExecError(rc, errorMessage);
}

static int
local_waitForDBAccess(void *unused, int numCalled)
{
	static struct timespec sleepTime;

	sleepTime.tv_sec  = 0;
	sleepTime.tv_nsec = LOCAL_SNOOZETIME_IN_NANOSECONDS;

	nanosleep(&sleepTime, NULL);

	return 1;
}

static void
local_beginTransaction(sqlite3_t db)
{
	int  rc;
	char *errmsg = NULL;

	rc = sqlite3_exec(db, "BEGIN EXCLUSIVE TRANSACTION",
	                  NULL, NULL, &errmsg);
	local_dealWithSQLiteExecError(rc, errmsg);
}

static void
local_endTransaction(sqlite3_t db)
{
	int  rc;
	char *errmsg = NULL;

	rc = sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &errmsg);
	local_dealWithSQLiteExecError(rc, errmsg);
}

static sqlite3_stmt_t
local_getPreparedStatement(sqlite3_t db, const char *stmtSQL)
{
	sqlite3_stmt_t stmt;
	int            rc;

	rc = sqlite3_prepare_v2(db, stmtSQL, -1, &stmt, NULL);
	if (rc != SQLITE_OK) {
		fprintf(stderr, "sqlite3 failure:  %s\n", sqlite3_errmsg(db));
		exit(EXIT_FAILURE);
	}

	return stmt;
}

inline static int64_t
local_getHaloIDFromIdx(int idx, uint64_t offset)
{
	return (int64_t)(idx + offset);
}

static void
local_createTableHalos(sqlite3_t db)
{
	int  rc;
	char *errmsg = NULL;

	rc = sqlite3_exec(db,
	                  "CREATE TABLE IF NOT EXISTS 'halos' (\n"
	                  "    haloID     INTEGER NOT NULL,\n"
	                  "    N          UNSIGNED INTEGER NOT NULL,\n"
	                  "    x          DOUBLE NOT NULL,\n"
	                  "    y          DOUBLE NOT NULL,\n"
	                  "    z          DOUBLE NOT NULL,\n"
	                  "    vx         DOUBLE NOT NULL,\n"
	                  "    vy         DOUBLE NOT NULL,\n"
	                  "    vz         DOUBLE NOT NULL,\n"
	                  "    Mvir       DOUBLE NOT NULL,\n"
	                  "    rvir       DOUBLE NOT NULL,\n"
	                  "    vmax       DOUBLE NOT NULL,\n"
	                  "    rmax       DOUBLE NOT NULL,\n"
	                  "    sigv       DOUBLE NOT NULL,\n"
	                  "    lambda     DOUBLE NOT NULL,\n"
	                  "    Lx         DOUBLE NOT NULL,\n"
	                  "    Ly         DOUBLE NOT NULL,\n"
	                  "    Lz         DOUBLE NOT NULL,\n"
	                  "    a          DOUBLE NOT NULL,\n"
	                  "    Eax        DOUBLE NOT NULL,\n"
	                  "    Eay        DOUBLE NOT NULL,\n"
	                  "    Eaz        DOUBLE NOT NULL,\n"
	                  "    b          DOUBLE NOT NULL,\n"
	                  "    Ebx        DOUBLE NOT NULL,\n"
	                  "    Eby        DOUBLE NOT NULL,\n"
	                  "    Ebz        DOUBLE NOT NULL,\n"
	                  "    c          DOUBLE NOT NULL,\n"
	                  "    Ecx        DOUBLE NOT NULL,\n"
	                  "    Ecy        DOUBLE NOT NULL,\n"
	                  "    Ecz        DOUBLE NOT NULL,\n"
	                  "    ovdens     DOUBLE NOT NULL,\n"
	                  "    redge      DOUBLE NOT NULL,\n"
	                  "    Ekin       DOUBLE NOT NULL,\n"
	                  "    Epot       DOUBLE NOT NULL,\n"
	                  "    mbpOffset  DOUBLE NOT NULL,\n"
	                  "    comOffset  DOUBLE NOT NULL,\n"
	                  "    r2         DOUBLE NOT NULL,\n"
	                  "    lambdaE    DOUBLE NOT NULL,\n"
	                  "    vesc       DOUBLE NOT NULL,\n"
	                  "    phi0       DOUBLE NOT NULL,\n"
	                  "    isSubHalo  BOOLEAN\n"
	                  ")",
	                  NULL, NULL, &errmsg);
	local_dealWithSQLiteExecError(rc, errmsg);
} /* local_createTableHalos */

static void
local_createTableProfiles(sqlite3_t db)
{
	int  rc;
	char *errmsg = NULL;

	rc = sqlite3_exec(db,
	                  "CREATE TABLE IF NOT EXISTS 'profiles' (\n"
	                  "    halo_id     INTEGER NOT NULL,\n"
	                  "    rmin        DOUBLE NOT NULL,\n"
	                  "    rmax        DOUBLE NOT NULL,\n"
	                  "    N_sphere    DOUBLE NOT NULL,\n"
	                  "    rho_shell   DOUBLE NOT NULL,\n"
	                  "    rho_sphere  DOUBLE NOT NULL,\n"
	                  "    vcirc       DOUBLE NOT NULL,\n"
	                  "    sigv_sphere DOUBLE NOT NULL,\n"
	                  "    Lx_sphere   DOUBLE NOT NULL,\n"
	                  "    Ly_sphere   DOUBLE NOT NULL,\n"
	                  "    Lz_sphere   DOUBLE NOT NULL,\n"
	                  "    a_sphere    DOUBLE NOT NULL,\n"
	                  "    Eax_sphere  DOUBLE NOT NULL,\n"
	                  "    Eay_sphere  DOUBLE NOT NULL,\n"
	                  "    Eaz_sphere  DOUBLE NOT NULL,\n"
	                  "    b_sphere    DOUBLE NOT NULL,\n"
	                  "    Ebx_sphere  DOUBLE NOT NULL,\n"
	                  "    Eby_sphere  DOUBLE NOT NULL,\n"
	                  "    Ebz_sphere  DOUBLE NOT NULL,\n"
	                  "    c_sphere    DOUBLE NOT NULL,\n"
	                  "    Ecx_sphere  DOUBLE NOT NULL,\n"
	                  "    Ecy_sphere  DOUBLE NOT NULL,\n"
	                  "    Ecz_sphere  DOUBLE NOT NULL,\n"
	                  "    Ekin_sphere DOUBLE NOT NULL,\n"
	                  "    Epot_sphere DOUBLE NOT NULL,\n"
	                  "    vesc_sphere DOUBLE NOT NULL\n"
	                  ")",
	                  NULL, NULL, &errmsg);
	local_dealWithSQLiteExecError(rc, errmsg);
} /* local_createTableProfiles */

static void
local_createTableHalosParticlesRel(sqlite3_t db)
{
	int  rc;
	char *errmsg = NULL;

	rc = sqlite3_exec(db,
	                  "CREATE TABLE IF NOT EXISTS 'halosParticlesRel' (\n"
	                  "    halo_id       INT NOT NULL,\n"
	                  "    particle_id   UNSIGNED BIGINT NOT NULL,\n"
	                  ")",
	                  NULL, NULL, &errmsg);
	local_dealWithSQLiteExecError(rc, errmsg);
}

#endif /* AHF_SQL */

#endif // AHF2

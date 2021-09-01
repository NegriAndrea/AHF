#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hgn.h"
#include <ctype.h>
#include <inttypes.h>

/** \file hgn.c
 *
 * Parse a HGN text file - incrementally.
 */



static const char g_aqarius[] = "Aquarius subhalo comparison project";
static const char g_aqarius2[] = "Aquarius comparison project";
static const char g_ghalo[] = "GHalo subhalo comparison project";
static const char g_level[] = "level:";
static const char g_author[] = "author:";
static const char g_halo[] = "halo finder:";
static const char g_date[] = "date:";


static char *trimcopy(char *str)
{
	char *ep = str + strlen(str);
	ep --;
	while(ep > str && isspace(*ep))
		*ep -- = 0;
	while(isspace(*str) && str < ep)
		str++;
	return strdup(str);
}

static HGN *readprologue(HGN *hgn) 
{
	int c;
	char linebuf[1024];

	/* skip over the opening comments picking up interesting ones */
	while ((c = getc(hgn->fp)) != EOF && c == '#') {
		char *cp;

		if (fgets(linebuf, sizeof linebuf, hgn->fp) == NULL)
		    break;
		for(cp = linebuf; *cp && isspace(*cp); cp++) 
			/*spin*/;
		if(strncasecmp(cp, g_aqarius, sizeof(g_aqarius)-1) == 0 ||
		   strncasecmp(cp, g_aqarius2, sizeof(g_aqarius2)-1) == 0) {
			hgn->idline = strdup(g_aqarius);
			/* good! */
		}
		else if(strncasecmp(cp, g_ghalo, sizeof(g_ghalo)-1) == 0 ||
			strncasecmp(cp, "GHalo", 5) == 0) {
			hgn->idline = strdup(g_ghalo);
			/* good! */
		}
		else if(strncasecmp(cp, g_level, sizeof(g_level)-1) == 0) {
			hgn->level = atoi(cp + sizeof(g_level)-1);
		}
		else if(strncasecmp(cp, g_author, sizeof(g_author)-1) == 0) {
			cp += sizeof(g_author)-1;
			hgn->author = trimcopy(cp);
		}
		else if(strncasecmp(cp, g_date, sizeof(g_date)-1) == 0) {
			cp += sizeof(g_date)-1;
			hgn->date = trimcopy(cp);
		}
		else if(strncasecmp(cp, g_halo, sizeof(g_halo)-1) == 0) {
			cp += sizeof(g_halo)-1;
			hgn->hf = trimcopy(cp);
		}
		else {
			/* other comments */
		}
	}
	ungetc(c, hgn->fp);
	if (fgets(linebuf, sizeof linebuf, hgn->fp) == NULL ||
	    sscanf(linebuf, "%" PRIu64, &hgn->count) != 1) {
	    closehgn(hgn);
	    return NULL;
	}

	return hgn;
}

HGN *openhaloidfile(const char *filename)
{
	HGN *hgn;
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) return NULL;

	hgn = calloc(1, sizeof *hgn);
	if (hgn == NULL) {
		fclose(fp);
		return NULL;
	}
	hgn->fp = fp;
	return readprologue(hgn);
}


HGN *openhaloidstream(FILE *fp)
{
	HGN *hgn;

	if (fp == NULL) return NULL;

	hgn = calloc(1, sizeof *hgn);
	if (hgn == NULL) return NULL;
	hgn->fp = fp;
	return readprologue(hgn);
}



static int fetchids(HGN *hgn, HGNHalo *halo) 
{
	char linebuf[1024];
	int n = halo->partcount;
	PartId partid;
	int ind = 0;

	while (n > 0 && fgets(linebuf, sizeof linebuf, hgn->fp) != NULL) {
		if(linebuf[0] == '#' || linebuf[0] == '\n') 
			continue;
		PartDesc *p = &(halo->partlist[ind]);
		if (sscanf(linebuf, "%" PRIu64 " %lf %lg %lg",
			   &p->id, &p->mass, &p->age, &p->z) == 4) {
		    n--;
		    ind ++;
		}
	}
	return n == 0;
}

/** \brief Read and allocate the next halo dscription
 * Reads the next halo in the file
 * Allocates a new HGNHalo structure which the caller must free
 * \return NULL on failure
 * \return HGNHalo on success
 */
HGNHalo *getnexthalo(HGN *hgn)
{
	HGNHalo *halo = NULL;
	char linebuf[1024];
	PartId haloid = 0, partcount = 0;

	while (fgets(linebuf, sizeof linebuf, hgn->fp) != NULL) {
		if(linebuf[0] == '#' || linebuf[0] == '\n') 
			continue;
		if (sscanf(linebuf, "%" PRIu64 "%" PRIu64, &partcount, &haloid) != 2)
			return NULL;
		break;
	}
	if (feof(hgn->fp) && partcount == 0) return NULL;
	halo = calloc(1, sizeof *halo);
	halo->haloid = haloid;
	halo->partcount = partcount;
	halo->partlist = calloc(partcount, sizeof(*(halo->partlist)));
	if (!fetchids(hgn, halo)) {
		freeHGNHalo(halo);
		return NULL;
	}
	return halo;
}

void freeHGNHalo(HGNHalo *halo)
{
	if (halo == NULL) return;
	if (halo->partlist) free(halo->partlist);
	free(halo);
}
				

void closehgn(HGN *hgn)
{
	if (hgn == NULL) return;
	if (hgn->idline) free(hgn->idline);
	if (hgn->author) free(hgn->author);
	if (hgn->date) free(hgn->date);
	if (hgn->hf) free(hgn->hf);
	if (hgn->fp) fclose(hgn->fp);
	free(hgn);
}

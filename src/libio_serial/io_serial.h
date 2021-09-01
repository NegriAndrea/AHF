#include "../param.h"
#include "../tdef.h"

#ifndef IO_SERIAL_INCLUDED
#define IO_SERIAL_INCLUDED

/* these are the major routines for writing/reading output files in BINARY */
void input            (char *infile_name);
void output           (char *outfile_name, int dumpflag);
void output_grid      (gridls *cur_grid, int dumpflag);

/* wrapper for the actual reading routine */
void read_data        (FILE *icfile);

/* ...the actual reading routines */
void read_amiga       (FILE *icfile);
void read_art         (FILE *icfile);
void read_mlapm       (FILE *icfile);
void read_tipsy       (FILE *icfile);
void read_mare_nostrum(FILE *icfile);
void read_deva        (FILE *icfile);
void read_gadget      (FILE *icfile);
void skim_gadget      (FILE *icfile);

/* these routines serve mainly debugging purposes and hence write information in ASCII */
void write_boundary         (gridls *cur_grid);
void write_density          (gridls *grid_ptr, char *prefix);
void write_nodepart         (gridls *grid);
void write_nodes            (gridls *grid, char *prefix);
void write_positions        (gridls *grid);
void write_residual         (gridls *grid_ptr, char *prefix);

#endif

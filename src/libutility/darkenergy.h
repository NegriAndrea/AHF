#ifndef DARKENERGY_H
#define DARKENERGY_H
/* 
* Dynamical dark energy-related features
*/

struct dark_energy {
	int n_entries;
	double *a;
	double *Hubble;
	double *OmegaM;
} DarkEnergy;

void read_dark_energy_table(char*);

double Hubble_DE(double);

double OmegaM_DE(double);

#endif

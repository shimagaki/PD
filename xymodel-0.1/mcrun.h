#ifndef MCRUN_H
#define MCRUN_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <argp.h>
#include <error.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <sys/time.h>
#include <time.h>
#include <assert.h>

#include <dlfcn.h>

#include "lattice.h"

struct statistics
{
	/* Energy */
	double H;
	double sumH;
	double sumH2;

	/* Magnetisation */
	double mag[2];
	double sumMag;
	double sumMag2;
};

struct arguments
{
	char *shm_file;
	int verbose;
	int average;
	double T;
	double iterations_per_spin;
	
	char *plugin;
	int plugin_samples;
	char *plugin_out_file;
	/* Apparently plugin_out can not be set in parse_opt, since the file descriptor
	 * is closed at exit. */
	FILE *plugin_out;
};

inline double getH(double (*lattice)[SIZE], int x, int y, double spin);
void update_stat(double (*lattice)[SIZE], double T);

#endif

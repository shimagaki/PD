#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "lattice.h"

#include "plugin/plugin.h"

/* Calculate helicity modulus (aka spin stiffness)
 */
size_t get_length()
{
	return 1;
}
static double avg_cos(double (*lattice)[SIZE], size_t x, size_t y);
static double avg_sin(double (*lattice)[SIZE], size_t x, size_t y);

/* secret magic formula :) */
void get_data(double *data, void *l, double T)
{
	double (*lattice)[SIZE] = l;

	data[0] = 0;
	double sin_term = 0;

	size_t x;
	for (x=0; x<SIZE; x++)
	{
		size_t y;
		for (y=0; y<SIZE; y++)
		{
			data[0] -= 0.5*avg_cos(lattice, x,y); /* -1/2 <U> */
			sin_term += avg_sin(lattice, x,y);
		}
	}

	data[0] -= 1/T*pow(sin_term,2);
	data[0] /= SIZE*SIZE;
}

static double avg_cos(double (*lattice)[SIZE], size_t x, size_t y)
{
	const double spin = lattice[x][y];
	double H;
	H = cos(spin - lattice[(x-1)&PBC][y]);
	H += cos(spin - lattice[(x+1)&PBC][y]);
	H += cos(spin - lattice[x][(y+1)&PBC]);
	H += cos(spin - lattice[x][(y-1)&PBC]);

	return -0.5*H;
}
static double avg_sin(double (*lattice)[SIZE], size_t x, size_t y)
{
	const double spin = lattice[x][y];
	double H;
	H = -sin(spin - lattice[(x-1)&PBC][y]);
	H += sin(spin - lattice[(x+1)&PBC][y]);

	return -0.5*H;
}

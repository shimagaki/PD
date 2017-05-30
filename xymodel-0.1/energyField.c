#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "lattice.h"

static const double J = -1; //TODO use argument

//TODO duplicated code
double getH(double (*lattice)[SIZE], int x, int y, double spin)
{
	double H;
	H = cos(spin - lattice[(x-1)&PBC][y]);
	H += cos(spin - lattice[(x+1)&PBC][y]);
	H += cos(spin - lattice[x][(y+1)&PBC]);
	H += cos(spin - lattice[x][(y-1)&PBC]);

	H *= -J;

	return H;
}

int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		printf("Usage: eneryField FILE\n");
		return 1;
	}

	void * shmaddr = get_shm_addr(argv[argc-1]);

	readXY(shmaddr);
	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);



	int i;
	for (i=0; i < SIZE; i++) 
	{
		size_t j;
		for (j=0; j< SIZE; j++) 
			printf("%d\t%d\t%f\n", i, j, getH(lattice, i, j, lattice[i][j]));
	}

	return 0;
}

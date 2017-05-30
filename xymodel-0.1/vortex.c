#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "lattice.h"

#include "vortex.h"
#include "plugin/plugin.h"

size_t get_length()
{
	return 3;
}

/* Count vortices */
void get_data(double *cor, void *l, double T)
{
	double (*lattice)[SIZE] = l;

	cor[0] = 0;
	cor[1] = 0;

	struct point p = next_vortex(lattice);
	while ( p.type != -1)
	{
		if (CLOCKWISE == p.type)
			cor[0]++;
		else
		{
			assert(COUNTERCLOCKWISE == p.type);
			cor[1]++;
		}
		p = next_vortex(lattice);
	}

	cor[2] = abs(cor[0]-cor[1]);
}

	__attribute__ ((visibility("hidden")))
int main(int argc, char *argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: ./vortex file\n");
		exit(1);
	}

	void *shmaddr = get_shm_addr(argv[argc-1]);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	struct point p = next_vortex(lattice);
	while (p.type != -1)
	{
		if (CLOCKWISE == p.type)
			printf("%d\t%d\t\t\n", p.x+1, p.y+1); // clockwise
		else
		{
			assert(COUNTERCLOCKWISE == p.type);
			printf("\t\t%d\t%d\n", p.x+1, p.y+1); // counterclockwise
		}
		p = next_vortex(lattice);
	}

	exit(0);
}

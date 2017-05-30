#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <error.h>

#include "lattice.h"

int main(int argc, char *argv[])
{
	if (argc != 2 && argc != 3)
	{
		fprintf(stderr, "Usage: writelattice SHM_FILE [OUT_FILE]\n"
		"\nWrite lattice to stdout.\n");
		exit(1);
	}

	void * shmaddr = get_shm_addr_rw(argv[1]);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	FILE* out = stdout;
	if (3 == argc)
	{
		out = fopen(argv[2], "w");
		if (NULL == out)
			error(1, errno, "open %s", argv[2]);
	}

	fprintf(out, "%d\n", SIZE);
	size_t x,y;
	for (x=0; x<SIZE; x++) 
	{
		for (y=0; y<SIZE; y++) 
			fprintf(out, "%a\t", lattice[x][y]);
		fprintf(out, "\n");
	}

	if (fclose(out))
	{
		if (out != stdout)
			error(1, errno, "close %s", argv[2]);
		else
			error(1, errno, "close");
	}



	exit(0);
}

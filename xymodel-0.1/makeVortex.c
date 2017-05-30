#include <stdlib.h>
#include <argp.h>
#include <sys/time.h>
#include <time.h>

#include "lattice.h"

const char *argp_program_version = "makeVortex 0.1";

static char doc[] = "makes vortices ;)";

static char args_doc[] = "FILE";

struct arguments
{
	char* shm_file;
	int verbose;
	int pos;
	int neg;
};

static struct argp_option options[] = {
	{"verbose",    'v', 0,       0, "Produce verbose output" },
	{"pos",        '+', 0,       0, "makes a pos vortex" },
	{"neg",        '-', 0,       0, "makes a neg vortex" },
	{ 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
		case '+':
			arguments->pos = 1;
			break;
		case '-':
			arguments->neg = 1;
			break;
		case 'v':
			arguments->verbose = 1;
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num > 0)
				argp_usage(state);
			arguments->shm_file = arg;
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
				argp_usage(state);
			break;
	}

	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

/* Returns H for one spin with its neighbours, but be careful
 * if you use this with a neighbour; you will count double. */
double getH(double (*lattice)[SIZE], int x, int y, double spin)
{
    double H;
	double J = 1;

    H = cos(spin - lattice[(x-1)&PBC][y]);
    H += cos(spin - lattice[(x+1)&PBC][y]);
    H += cos(spin - lattice[x][(y+1)&PBC]);
    H += cos(spin - lattice[x][(y-1)&PBC]);

    H *= -J;

    return H;
}

void makeVortex(double (*lattice)[SIZE], int (*fixMap)[SIZE], size_t x, size_t y, int type)
{
	lattice[x][y+1] = 1 * M_PI/4;
	lattice[x+1][y] = 5 * M_PI/4;
	
	if (1 == type)
	{
		lattice[x][y] = 3 * M_PI/4; 
		lattice[x+1][y+1] = 7 * M_PI/4;
	}
	else
	{
		lattice[x][y] = 7 * M_PI/4; 
		lattice[x+1][y+1] = 3 * M_PI/4;
	}

	fixMap[x  ][y  ] = 1;
	fixMap[x  ][y+1] = 1;
	fixMap[x+1][y  ] = 1;
	fixMap[x+1][y+1] = 1;
}

int main(int argc, char *argv[])
{
	struct arguments arguments;

	arguments.verbose = 0;
	arguments.pos = 0;
	arguments.neg = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	void * shmaddr = get_shm_addr_rw(arguments.shm_file);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	int fixMap[SIZE][SIZE];
	
	if (arguments.verbose)
		printf("Writing to lattice. Size: %d\n", SIZE);

	// init random first
	size_t x,y;
	for (x=0; x< SIZE; x++)
		for (y=0; y< SIZE; y++)
		{	
				int dx = x-SIZE/2;
				int dy = y-SIZE/2;
				double sign = 0;

				if(dx >= 0) dx++;
				if(dy >= 0) dy++;
				
				sign -= M_PI/2;
				if (dx < 0)
					sign += M_PI;
				
				if (arguments.neg)
				{
					if(dx*dy> 0)
					sign += M_PI;
				}

				
				lattice[x][y] =  atan(((float)dy) / dx) + sign;
				if (lattice[x][y] < 0)
					lattice[x][y] += 2*M_PI;
				if (lattice[x][y] >= 2*M_PI)
					lattice[x][y] -= 2*M_PI;
				fixMap[x][y] = 0;
		}
//lazy ;)
	for (x=0; x< SIZE; x++)
	        for (y=0; y< SIZE/5; y++)
				lattice[x][y] = 0;
	for (x=0; x< SIZE; x++)
	        for (y=SIZE; y> SIZE-SIZE/5; y--)
				lattice[x][y] = 0;
	for (x=0; x< SIZE/5; x++)
	        for (y=0; y< SIZE; y++)
				lattice[x][y] = 0;
	for (x=SIZE; x> SIZE-SIZE/5; x--)
	        for (y=0; y< SIZE; y++)
				lattice[x][y] = 0;
	
	// make vortices
	/*
	//if (arguments.pos)
		makeVortex(lattice, fixMap, 2 * SIZE/5, SIZE/2, 1);
	//else
		makeVortex(lattice, fixMap, 3 * SIZE/5, SIZE/2, -1);

	// clean up
	size_t iterations = 1000 * SIZE * SIZE;
	double T;
	for (T=.1; T <1.5;  T+=.1)
	{
    	size_t j;
	    for (j=0; j<iterations; j++)
    	{
        	const int rnd = rand();
    	    const int x =  rnd & PBC;
	        const int y =  (rnd>>LATTICE_POW) & PBC;

    	    if (!fixMap[x][y])
	        {
        	    const double try_spin = rand_angle();
    	        const double dH = getH(lattice, x, y, try_spin) -
	                getH(lattice, x, y, lattice[x][y]);
            	if (dH < 0 || (T != 0 && drand48() < exp(-dH/T)) )
        	        lattice[x][y] = try_spin;
    	    }
		}
	}
*/
	return 0;
}

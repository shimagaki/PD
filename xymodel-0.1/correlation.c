#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>
#include <argp.h>
#include <error.h>

#include "lattice.h"
#include "plugin/plugin.h"


static const int RADIUS = 20;

const char *argp_program_version = "correlation 0.1";
static char doc[] = 
"Read square lattice from shared memory segment of FILE "
"and write correlation in to stdout.\n"
"Format: Distance Correlation Variance";


static char args_doc[] = "FILE";

struct arguments
{
	char* shm_file;
	int verbose;
};

static struct argp_option options[] = {
	{"verbose",    'v', 0,       0, "Produce verbose output" },
	{ 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
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


static int get_nr_of_nearest(int radius)
{
	int count = 1; // don't forget yourself
	size_t i;
	for(i=1; i<=radius; i++)
	{
		size_t j;
		for(j=0; j*j <= radius*radius - i*i; j++)
			count += 4;
	}

	return count;
}

/* not fast, but only called once per run.
 * mallocs - a memory leak waiting to happen! */
static int (*get_nearest(int radius))[2]
{
	int count = get_nr_of_nearest(radius);

	int (*points)[2] = malloc(count*2*sizeof(int));

	points[0][0] = 0;
	points[0][1] = 0;
	count = 1;
	int i, j;
	for(i=1; i <= radius; i++)
		for(j=0; j*j <= radius*radius - i*i; j++)
		{
			// 1st quadrant
			points[count][0]   = i;
			points[count++][1] = j;
			// 2nd quadrant
			points[count][0]   = -j;
			points[count++][1] = i;
			// 3rd quadrant
			points[count][0]   = -i;
			points[count++][1] = -j;
			// 4th quadrant
			points[count][0]   = j;
			points[count++][1] = -i;
		}

	return points;
}

static int get_radius()
{
	const int max_rad = (SIZE < SIZE)? (SIZE-1)/2 : (SIZE-1)/2;
	return (RADIUS < max_rad)? RADIUS : max_rad;
}

size_t get_length()
{
	return get_nr_of_nearest(get_radius());
}

static int (*nearest)[2] = NULL;
//TODO OMFG copy & paste!
void get_data(double *cor, void *l, double T)
{
	double (*lattice)[SIZE] = l;
	int radius = get_radius(); 

	int nr_nearest = get_nr_of_nearest(radius);
	if (NULL == nearest)
		nearest = get_nearest(radius);

	size_t i;
	// initialize
	for (i=0; i < nr_nearest; i++)
		cor[i] = 0;

	size_t x,y;
	for (x=0; x < SIZE; x++)
		for (y=0; y < SIZE; y++)
			for(i=0; i < nr_nearest; i++) {
				double c = cos(lattice[x][y] - lattice[(x+nearest[i][0])&PBC][(y+nearest[i][1])&PBC]);
				cor[i] += c;
			}

	for (i=0; i<nr_nearest; i++)
		cor[i] /= SIZE*SIZE;

}

void get_x(double *x)
{
	if (NULL == nearest)
		nearest = get_nearest(get_radius());
	size_t i;
	for (i=0; i < get_length(); i++)
		x[i] = sqrt(nearest[i][0]*nearest[i][0] + nearest[i][1]*nearest[i][1]);
}

int main(int argc, char *argv[])
{
	struct arguments arguments;
	arguments.verbose = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	void * shmaddr = get_shm_addr(arguments.shm_file);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);


	struct timeval start;
	gettimeofday(&start, NULL);

	int radius = get_radius();

	int nr_nearest = get_nr_of_nearest(radius);
	int (*nearest)[2] = get_nearest(radius);
	double cor[nr_nearest];
	double cor2[nr_nearest];
	double dist[nr_nearest];
	int n = SIZE * SIZE;

	size_t i;
	// initialize
	for(i=0; i < nr_nearest; i++)
	{
		cor[i] = 0;
		cor2[i] = 0;
	}

	// distances
	for(i = 0; i < nr_nearest; i++)
		dist[i] = sqrt(nearest[i][0]*nearest[i][0] + nearest[i][1]*nearest[i][1]);

	size_t x,y;
	for (x=0; x < SIZE; x++)
		for (y=0; y < SIZE; y++)
			for(i=0; i < nr_nearest; i++)
			{
				double c = cos(lattice[x][y] - lattice[(x+nearest[i][0])&PBC][(y+nearest[i][1])&PBC]);
				cor[i] += c;
				cor2[i] += c*c;
			}

	// cleanup
	free(nearest);

	if (arguments.verbose)
	{
		struct timeval end;
		gettimeofday(&end, NULL);
		double time = end.tv_sec+ end.tv_usec * 1E-6
			-start.tv_sec -start.tv_usec * 1E-6;

		fprintf(stderr, "Correlation calculation took %g seconds. %g spins/sec\n", time, (n * radius)/time);
	}

	for (i=0; i<nr_nearest; i++)
		printf("%e\t%e\t%e\n", dist[i], cor[i]/n, 1.0/(n-1) * (cor2[i] - 1.0/n * cor[i]*cor[i]) );

	return 0;
}

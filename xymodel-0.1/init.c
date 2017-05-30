#include <stdlib.h>
#include <argp.h>
#include <sys/time.h>
#include <time.h>

#include "lattice.h"

const char *argp_program_version = "init 0.2";

static char doc[] = "Write spins of a square lattice to the shared memory segment of FILE";

static char args_doc[] = "FILE";

struct arguments
{
	char* shm_file;
	int verbose;
	int random;
	double theta;
};

static struct argp_option options[] = {
	{"verbose",    'v', 0,       0, "Produce verbose output" },
	{"random",     'r', 0,       0, "makes random spins" },
	{"theta",      't', "VALUE", 0, "makes pins turned to the given angle (in degrees)" },
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
		case 'r':
			arguments->random = 1;
			break;
		case 't':
			errno = 0;
			arguments->theta = strtod(arg, NULL);
			if (errno == 0)
			{
				arguments->theta *= 2.f*M_PI/360.f;
			}
			else
			{
				perror("theta");
				argp_usage(state);
			}
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

int main(int argc, char *argv[])
{
	struct arguments arguments;

	arguments.verbose = 0;
	arguments.random = 0;
	arguments.theta = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	void *shmaddr = get_shm_addr_rw(arguments.shm_file);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	if (arguments.verbose)
		printf("Writing to lattice. Size: %d\n", SIZE);

	size_t x,y;
	for (x=0; x< SIZE; x++)
		for (y=0; y< SIZE; y++)
			if (arguments.random)
				lattice[x][y] = rand_angle();
			else 
			{
				double theta = (double)x/(SIZE-1) * arguments.theta;
				lattice[x][y] = (theta < 0) ? theta + 2*M_PI : theta;
			}

	return 0;
}

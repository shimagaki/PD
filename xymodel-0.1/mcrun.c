#include "mcrun.h"
#include "mc_plugs.h"

#include <string.h>
//TODO better variable names.
static struct statistics fstat;

const char *argp_program_version = "mcrun 0.1";
static char doc[] = 
"Run Monte Carlo simulation on lattice associated with FILE "
"with ITERATIONS per spin.";

static char args_doc[] = "ITERATIONS FILE";

static struct argp_option options[] = {
	{"verbose",     'v', 0,          0, "Produce verbose output" },
	{"temperature", 't', "VALUE",    0, "Set temperature" },
	{"average",     'a', 0,          0, "Average" },
	{"plug-in",     'p', "PLUG-IN",  0, "Load a plug-in (.so file)" },
	{"plug-in-out", 'o', "FILE",     0, "Output statistics to FILE" },
	{"plug-in-samples",'s',"SAMPLES",0, "Number of samples to take with plug-in (default: 50)"},
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
		case 'a':
			arguments->average = 1;
			break;
		case 'p':
			arguments->plugin = arg;
			break;
		case 't':
			errno = 0;
			arguments->T = strtod(arg, NULL);
			if (errno != 0)
			{
				perror("temperature");
				argp_usage(state);
			}
			break;
		case 'o':
			arguments->plugin_out_file = arg;
			break;
		case 's':
			errno = 0;
			arguments->plugin_samples = strtol(arg, NULL, 10);
			if (errno != 0)
			{
				perror("plug-in-samples");
				argp_usage(state);
			}
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num > 1)
				argp_usage(state);
			if (0 == state->arg_num)
			{
				errno = 0;
				arguments->iterations_per_spin = strtod(arg, NULL);
				if (errno != 0)
				{
					perror("iterations");
					argp_usage(state);
				}
			}
			else if (1 == state->arg_num)
				arguments->shm_file = arg;
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 2)
				argp_usage(state);
			break;
	}

	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

static const double J = 1;


inline void run(double (*lattice)[SIZE], const size_t iterations,
		const double T, const int average)
{
	/* Make sure random() is big enough for both x and y
	 * max SIZE is something like 2^15 = 32768 */
	assert(RAND_MAX>>LATTICE_POW > SIZE);

	size_t i;
	for (i=0; i<iterations; i++)
	{
		const int rnd = random();
		const int x =  rnd & PBC;
		const int y =  (rnd>>LATTICE_POW) & PBC;

		const double try_spin = rand_angle();
		const double dH = getH(lattice, x, y, try_spin) -
			getH(lattice, x, y, lattice[x][y]);

		if (dH <= 0 || (T != 0 && drand48() < exp(-dH/T)) )
		{
			if (average)
			{
				fstat.mag[0] -= sin(lattice[x][y]);
				fstat.mag[1] -= cos(lattice[x][y]);
			}

			lattice[x][y] = try_spin;

			if (average)
			{
				fstat.mag[0] += sin(lattice[x][y]);
				fstat.mag[1] += cos(lattice[x][y]);

				fstat.H += dH;
			}
		}

		if (average)
		{
			fstat.sumH += fstat.H;
			fstat.sumH2 += pow(fstat.H,2);

			double mag = sqrt(pow(fstat.mag[0],2) + pow(fstat.mag[1],2));
			fstat.sumMag += mag;
			fstat.sumMag2 += pow(mag,2);

		}
	}

}
double T=1;

int main(int argc, char *argv[])
{
	struct arguments arguments;
	arguments.shm_file = NULL;
	arguments.verbose = 0;
	arguments.T = 0;
	arguments.average = 0;

	arguments.plugin = NULL;
	arguments.plugin_out_file = NULL;
	arguments.plugin_out = stdout;
	arguments.plugin_samples = 50;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	if (arguments.plugin_out_file != NULL) {
		arguments.plugin_out = fopen(arguments.plugin_out_file, "w");
		if (arguments.plugin_out == NULL)
			error(1, errno, "open plugin output file %s failed", arguments.plugin_out_file);
	}

	void * shmaddr = get_shm_addr_rw(arguments.shm_file);

	readXY(shmaddr);
	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	const size_t iterations = arguments.iterations_per_spin * SIZE * SIZE;

	struct timeval start;
	gettimeofday(&start, NULL);

	if (arguments.average)
	{
		size_t i,j;
		memset((void *)&fstat, 0, sizeof(struct statistics));
		for (i = 0; i < SIZE; i++)
		{
			for (j = 0; j < SIZE; j++)
			{
				fstat.H += 0.5*getH(lattice, i, j, lattice[i][j]);	// don't count twice
				fstat.mag[0] += sin(lattice[i][j]);
				fstat.mag[1] += cos(lattice[i][j]);
			}
		}
	}

	if (arguments.average || arguments.plugin)
	{
		if (NULL == arguments.plugin)
			run(lattice, iterations, arguments.T, 1);
		else
		{
			if (arguments.verbose)
				fprintf(stderr, "Loading libs\n");
			init_stat(arguments.plugin);

			int i;
			for (i=0; i<arguments.plugin_samples; i++)
			{
				run(lattice, iterations/arguments.plugin_samples, arguments.T, arguments.average);
				update_stat(lattice, arguments.T);
			}
		}
	}
	else
		run(lattice, iterations, arguments.T, 0);

	struct timeval end;
	gettimeofday(&end, NULL);

	if (arguments.verbose)
	{
		double time = end.tv_sec+ end.tv_usec * 1E-6
			-start.tv_sec -start.tv_usec * 1E-6;

		fprintf(stderr, "Simulation took %G seconds. %G spins/sec\n", time, iterations/time);
	}

	if (arguments.average)
	{
		fprintf(stderr, "n = %d\n", iterations);

		fstat.sumH /= SIZE*SIZE;
		fstat.sumH2 /= SIZE*SIZE * SIZE*SIZE;
		printf("%G\t%G", fstat.sumH/iterations,
				(fstat.sumH2 - pow(fstat.sumH,2)/iterations)/(iterations-1));

		fstat.sumMag /= SIZE*SIZE;
		fstat.sumMag2 /= SIZE*SIZE * SIZE*SIZE;
		printf("\t%G\t%G\n", fstat.sumMag/iterations,
				(fstat.sumMag2 - pow(fstat.sumMag,2)/iterations)/(iterations-1));
	}

	if (NULL != arguments.plugin)
	{
		output_stat(arguments.plugin_out, &pstat, arguments.plugin_samples);

		if (arguments.plugin_out_file != NULL)
			if (fclose(arguments.plugin_out)) 
				error(1, errno, "failed to close output file %s", arguments.plugin_out_file);
	}

	return 0;

}

/* Returns H for one spin with its neighbours, but be careful
 * if you use this with a neighbour; you will count double. */
inline double getH(double (*lattice)[SIZE], int x, int y, double spin)
{
	double H;
	H = cos(spin - lattice[(x-1)&PBC][y]);
	H += cos(spin - lattice[(x+1)&PBC][y]);
	H += cos(spin - lattice[x][(y+1)&PBC]);
	H += cos(spin - lattice[x][(y-1)&PBC]);

	H *= -J;

	return H;
}

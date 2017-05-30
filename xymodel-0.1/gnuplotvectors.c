#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <argp.h>
#include <signal.h>

#include "lattice.h"

const char *argp_program_version = "gnuplotvectors 0.1";

static char doc[] = 
"Read square lattice from shm of FILE and write vector data, "
"to be used by gnuplot (plot \"datafile\" with vectors\"), "
"to stdout\n"
"Deprecated, use glvec";

static char args_doc[] = "FILE";

struct arguments
{
	char* shm_file;
	int x11;
	double update_intervall; 
};

static struct argp_option options[] = {
	{"x11", 'x', "UPDATE",  OPTION_ARG_OPTIONAL,
		"Plot in a X11 window. Optionally update every UPDATE seconds." },
	{ 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
		case 'x':
			arguments->x11 = 1;
			if (arg != NULL)
			{
				errno = 0;
				arguments->update_intervall = strtod(arg, NULL);
				if (errno || arguments->update_intervall < 0)
				{
					perror ("Update interval");
					argp_usage(state);
				}
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
	arguments.x11 = 0;
	arguments.update_intervall = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	FILE *stream = stdout;
	if (arguments.x11)
	{
		stream = popen("gnuplot", "w");
		if (NULL == stream)
		{
			perror("gnuplot");
			exit(1);
		}
	}

	void * shmaddr = get_shm_addr(arguments.shm_file);

	readXY(shmaddr);
	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);


	if(arguments.x11)
	{
		fprintf(stream, "set nokey\n");
		if (arguments.update_intervall != 0)
			fprintf(stream, "set term x11 title \"XY model (%s)\" noraise\n",
					arguments.shm_file);
		else
			fprintf(stream, "set term x11 title \"XY model (%s)\" noraise persist\n",
					arguments.shm_file);
	}


	int i;
	do
	{
		if(arguments.x11)
			fprintf(stream, "\n\nplot \'-\' with vectors\n");

		for (i=0; i < SIZE; i++) 
		{
			size_t j;
			for (j=0; j< SIZE; j++) 
			{
				double angle = lattice[i][j];

				fprintf(stream, "%e\t%e\t%e\t%e\n",
						i + 0.5 - 0.4 * cos(angle),
						j + 0.5 - 0.4 * sin(angle),
						0.8*cos(angle),
						0.8*sin(angle));
			}
		}

		if(arguments.x11)
		{
			fprintf(stream, "e\n");
			usleep(arguments.update_intervall);
		}
		else
			break;
	}
	while(arguments.update_intervall != 0);

	return 0;
}

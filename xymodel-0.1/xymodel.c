#include "xymodel.h"

const char *argp_program_version = "xymodel 0.1";
static char doc[] = 
"Create a lattice of a XY model in a shared memory segment "
"(shm) associated with FILE.\n"

"Write lattice to FILE at exit (SIGTERM || SIGINT) when option -w is given.\n"

"If no size is specified, the lattice is read from FILE, "
"or from stdin if '-c' is given.\n"
"\nKeep this program running to prevent the shm from being freed.\n";

static char args_doc[] = "FILE";


static struct argp_option options[] = {
	{"verbose",    'v', 0,      0, "Produce verbose output" },
	{"stdin",      'c', 0,      0, "Read lattice from stdin" },
	{"size",       's', "SIZE", 0, "Set lattice size" },
	{"write",      'w', 0,      0, "Write lattice to FILE" },
	{"background", 'b', 0,      0, "Run in background TODO" },
	{ 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
		case 'c':
			arguments->input = stdin;
			break;
		case 'v':
			arguments->verbose = 1;
			break;
		case 'w':
			arguments->write = 1;
			break;
		case 'b':
			arguments->bg = 1;
			break;
		case 's':
			SIZE = atoi(arg);
			LATTICE_POW = 0;
			check_dimensions();
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
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
	arguments.input = NULL;
	arguments.shm_file = NULL;
	arguments.verbose = 0;
	arguments.write = 0;
	arguments.bg = 0;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);


	if (NULL == arguments.input && 0 == SIZE)
	{
		arguments.input = fopen(arguments.shm_file, "r");
		if (NULL == arguments.input)
			error(1, errno, "error open %s", arguments.shm_file);
	}

	if (arguments.input != NULL)
		get_dimensions(arguments.input);
	if (arguments.verbose)
		fprintf(stderr, "lattice size %d\n", SIZE);

	FILE* out_file=NULL;
	if (arguments.write)
	{
		out_file = fopen(arguments.shm_file, "w");
		if (NULL == out_file)
			error(1, errno, "error open %s for writing", arguments.shm_file);
	}

	void* shmaddr = get_shm_addr_rw(arguments.shm_file);

	double (*lattice)[SIZE];
	lattice = get_shm_lattice(shmaddr);

	if (arguments.input)
	{
		read_lattice_data(lattice, arguments.input);
		if (arguments.input != stdin)
			fclose(arguments.input);
	}


	struct sigaction sa;
	sigemptyset(&sa.sa_mask);
	sa.sa_handler = quit_signal;
	sa.sa_flags = SA_ONESHOT;

	sigaction(SIGTERM, &sa, NULL);
	sigaction(SIGINT, &sa, NULL);

	if (arguments.bg)
	{
		printf("Detaching\n");
		if (setsid() < 0)
			error(1, errno, "setsid");
	}
	else
		printf("Press Ctrl-c to quit\n");

	pause();

	if (arguments.write)
	{
		if (arguments.verbose)
			printf("writing lattice\n");

		fprintf(out_file, "%d\n", SIZE);
		size_t x,y;
		for (x=0; x< SIZE; x++) 
		{
			for (y=0; y< SIZE; y++) 
				fprintf(out_file, "%a\t", lattice[x][y]);
			fprintf(out_file, "\n");
		}

		if(fclose(out_file))
			error(0, errno, "error close file %s", arguments.shm_file);
	}

	shmctl(shmid, IPC_RMID, NULL);//delete shm when last process detaches.
	if (shmdt(shmaddr))
		error(1, errno, "error delete shm");

	exit(0);
}

void quit_signal(int arg)
{
	//do nothing (i.e. just catch signals)
}

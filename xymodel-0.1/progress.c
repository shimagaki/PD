#define _GNU_SOURCE //getline() support

#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_LIBCURSES
#include <curses.h>
#include <term.h>
#endif

#include <errno.h>
#include <error.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <argp.h>

/*
 * TODO:
 * - predict ETA better, by taking slowdown for
 *   high temperatures in account
 * - handle resize
 */

#define TERMINAL_LINE_LENTH 80


const char *argp_program_version = "progress 0.1";

static char doc[] = "Print progress bar and optionally ETA";

static char args_doc[] = "FIFO TOTAL [OFFSET]";

struct progress_bar
{
	double pos;
	double total;
	double offset;
	double start_time;
	int length;
	int ETA;
};

struct arguments
{
	char* fifo;
	struct progress_bar* bar;
};

static struct argp_option options[] = {
	{"eta", 'e', 0, 0, "Print ETA" },
	{ 0 }
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
		case 'e':
			arguments->bar->ETA = 1;
			break;
		case ARGP_KEY_ARG:
			switch(state->arg_num)
			{
				case 0:
					arguments->fifo = arg;
					break;
				case 1:
					errno = 0;
					arguments->bar->total = strtod(arg, NULL);
					if (errno != 0)
						error(1, errno, "total");
					break;
				case 2:
					arguments->bar->offset = strtod(arg, NULL);
					if (errno != 0)
						error(1, errno, "offset");
					arguments->bar->total -= arguments->bar->offset;
					break;

				default:
					argp_usage(state);
			}
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 2)
				argp_usage(state);
			break;
	}

	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};



static void print_initial(const struct progress_bar* bar)
{
	printf("  0%% [");
	int i;
	for (i=0; i<bar->length; i++)
		printf(".");
	printf("]");
}

static void print_time(double sec)
{
	time_t time = sec;
	char *time_str = asctime(localtime(&time));
	time_str[strlen(time_str)-1] = '\0';
	/* asctime() & localtime() return a pointer to statically allocated
	 * memory. i.e: no memleak :) */
	printf(time_str);
}

static void print_eta(const struct progress_bar* bar)
{
	struct timeval cur_time;
	if (gettimeofday(&cur_time, NULL))
		error(1, errno, "gettimeofday");

	double dT = cur_time.tv_sec+cur_time.tv_usec * 1E-6 - bar->start_time;


	printf(" ETA: ");
	print_time(cur_time.tv_sec + dT/bar->pos*(bar->total-bar->pos));
}

void print_bar(const struct progress_bar* bar)
{
	const int percent = 100*bar->pos / bar->total;

	if(percent >= 100)
		printf("100%% [");
	else
	{
		if (percent < 10)
			printf(" ");
		printf(" %d%% [", percent);
	}

	int length = bar->pos * bar->length / bar->total;
	if (length > bar->length)
		length = bar->length;
	int i;
	for (i=0; i<length; i++)
		printf("=");
	for (;i<bar->length; i++)
		printf(".");
	printf("]");

}

void start_bar(struct progress_bar* bar)
{
#ifdef HAVE_LIBCURSES
	if (setupterm(NULL, fileno(stdout), NULL) != OK)
		error(1, errno, "setupterm");

	bar->length = tigetnum("cols");
	if (ERR == bar->length)
		bar->length = TERMINAL_LINE_LENTH;
#else
	bar->length = TERMINAL_LINE_LENTH;
#endif
	bar->length -= 7;
	if (bar->ETA)
		bar->length -= 30;

	print_initial(bar);

	if (bar->ETA)
	{
		struct timeval start;
		if (gettimeofday(&start, NULL))
			error(1, errno, "gettimeofday(start)");

		bar->start_time = start.tv_sec+start.tv_usec * 1E-6;
	}
}

int main(int argc, char* argv[])
{
	struct progress_bar bar;
	bar.ETA = 0;
	bar.total = 0;
	bar.offset = 0;

	struct arguments arguments;
	arguments.fifo = NULL;
	arguments.bar = &bar;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	FILE *in = fopen(arguments.fifo, "r");
	if (NULL == in)
		error(1, errno, "open %s", argv[1]);

	start_bar(&bar);

	size_t len = 0;
	char* line = NULL;
	while (getline(&line, &len, in) != -1)
	{
		errno = 0;
		bar.pos = strtod(line, NULL) - bar.offset;
		if (errno != 0)
			error(1, errno, "readline");

		if (bar.pos < 0)
		{
			printf("\n");
			break;
		}

		printf("\r");
		print_bar(&bar);

		if (bar.pos > 0 && bar.ETA)
			print_eta(&bar);

		if(bar.pos >= bar.total)
		{
			printf("\n");
			break;
		}

		fflush(stdout);
	}
	exit(0);
}

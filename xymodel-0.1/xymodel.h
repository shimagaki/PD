#ifndef XYMODEL_H
#define XYMODEL


#define _GNU_SOURCE //getline() support
#include <errno.h>
#include <error.h>
#include <stdlib.h>
#include <wchar.h>
#include <signal.h>
#include <unistd.h>
#include <argp.h>
#include <string.h>

#include <sys/shm.h>

#include <fcntl.h>

#include "lattice.h"

struct arguments
{
	char* shm_file;
	FILE* input;
	int verbose;
	int write;
	int bg;
};

void get_dimensions(FILE* in);
void read_lattice(double (*lattice)[SIZE], FILE* in);
void quit_signal(int arg);
void check_dimensions();

#endif

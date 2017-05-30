#define _GNU_SOURCE //getline() support
#include <math.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <error.h>

#include "lattice.h"


size_t SIZE = 0;
size_t LATTICE_POW = 0;
size_t PBC;

int shmid;

static void *get_shm_addr2(const char* path, int rw);

struct lattice_property
{
	size_t size_pow;
};

void *get_shm_addr_rw(const char *path)
{
	return get_shm_addr2(path, 1);
}
void *get_shm_addr(const char *path)
{
	return get_shm_addr2(path, 0);
}

void *get_shm_lattice(void *shmaddr)
{
	return shmaddr + sizeof(struct lattice_property);
}

void readXY(void *shmaddr)
{
	const struct lattice_property *args = shmaddr;
	LATTICE_POW = args->size_pow;
	SIZE = 1<<LATTICE_POW;
	PBC = SIZE-1;
}

void writeXY(void *shmaddr)
{
	struct lattice_property *args = shmaddr;
	if (SIZE != 0 && LATTICE_POW != 0)
		assert(1<<LATTICE_POW == SIZE);
	else if (SIZE != 0)
	{
		assert(0 == (SIZE &(SIZE-1)));
		assert(0 == LATTICE_POW);
		size_t tmp_size = SIZE>>1;
		while (tmp_size)
		{
			tmp_size = tmp_size >> 1;
			LATTICE_POW++;
		}
		assert(1<<LATTICE_POW == SIZE);
	}
	args->size_pow = LATTICE_POW;
}

void *get_shmaddr(key_t key, int flags, int rw)
{
	size_t size = sizeof(struct lattice_property)+SIZE*SIZE*sizeof(double);

	if (-1 == (shmid = shmget(key, size, 0644 | flags))) {
		switch (errno)
		{
			case EEXIST:
				error(1, 0,"It's likely that the shm is already in use. "
						"Try to kill the offending programs. (mainly ./xymodel)");
			case ENOENT:
				error(0, 0, "No shm available - start ./xymodel first to use it.");
				return NULL;
			default:
				error(1, errno, "shmget");
		}
	}

	int shmflg = SHM_RDONLY;
	if (rw)
		shmflg = 0;

	void *shmaddr = shmat(shmid, (void *) 0, shmflg);
	if ((void *)-1 == shmaddr)
		error(1, errno, "shmat");

	return shmaddr;
}

static void *read_lattice(const char* path)
{
	FILE *in = fopen(path, "r");
	if (NULL == in)
		error(1, errno, "open %s", path);

	get_dimensions(in);

	void *fake_shm_addr = malloc(sizeof(struct lattice_property)
			+ SIZE*SIZE*sizeof(double));

	writeXY(fake_shm_addr);

	read_lattice_data(get_shm_lattice(fake_shm_addr), in);

	if (fclose(in))
		error(0, errno, "read lattice");

	return fake_shm_addr;
}

static void *get_shm_addr2(const char* path, int rw)
{
	key_t key;

	if ( -1 == (key = ftok(path, 'X')) )
		error(1, errno, "error accessing %s", path);

	int flg = 0;
	if (SIZE != 0)
	{
		assert(1 == rw);
		flg = IPC_CREAT | IPC_EXCL;
	}
	void *shmaddr = get_shmaddr(key, 0644 | flg, rw);

	if (NULL == shmaddr)
	{
		error(0, 0, "Accessing file directly.");
		shmaddr = read_lattice(path);
	}

	if (0 == SIZE) 
	{
		readXY(shmaddr);
		shmdt(shmaddr);

		shmaddr = get_shmaddr(key, 0644, rw);
	}
	else
		writeXY(shmaddr);

	return shmaddr;
}

void get_dimensions(FILE* in)
{
	size_t len = 0;
	char *line = NULL;

	if (-1 == getline(&line, &len, in))
		error(1, 0, "Error: no lines read\n");


	errno = 0;
	char* next = NULL;
	SIZE = strtol(line, &next, 10);
	if (errno != 0)
		error(1, errno, "strtol SIZE");

	if (*next != '\n')
		error(0, 0, "Warning: old data format\n");

	check_dimensions();

	free(line);
}

void check_dimensions() 
{
	if (SIZE & (SIZE-1))
		error(1, 0, "size (%d) has to be power of two", SIZE);
}

void read_lattice_data(double (*lattice)[SIZE], FILE* in)
{
	size_t len = 0;
	char* line = NULL;

	size_t x=0;
	while (getline(&line, &len, in) != -1)
	{
		size_t y=0;

		char* cur = line;
		char* next = NULL;

		while(1)
		{
			errno = 0;
			double val = strtod(cur, &next);
			if (next == cur)
				break;

			if (errno != 0)
				error(1, errno, "strtol read lattice");

			if (y >= SIZE)
				error(1, 0, "line %d to long.", x+1);

			lattice[x][y] = val;

			cur = next;
			next = NULL;
			y++;
		}
		if (y != SIZE)
			error(1, 0, "line %d to short.\n", x+1);

		if (x >= SIZE)
			error(1, 0, "to may lines");

		x++;
	}
	if (x != SIZE)
		error(1, 0, "to few lines (read %d).", x);

	free(line);
}

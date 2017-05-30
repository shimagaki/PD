#ifndef LATTICE_H
#define LATTICE_H

#include <math.h>

extern size_t SIZE;
extern size_t LATTICE_POW;
extern size_t PBC;

extern int shmid;

static inline double rand_angle()
{
	return 2*M_PI*drand48();
}

void *get_shm_addr_rw(const char *path);
void *get_shm_addr(const char *path);

void *get_shm_lattice(void *shmaddr);

void readXY(void *shmaddr);
void writeXY(void *shmaddr);

void *get_shmaddr(key_t key, int flags, int rw);

void check_dimensions();

void get_dimensions(FILE* in);
void read_lattice_data(double (*lattice)[SIZE], FILE* in);

#endif

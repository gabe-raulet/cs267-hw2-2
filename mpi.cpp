#include "common.h"
#include "ParticleStore.h"
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <math.h>
#include <assert.h>
#include <mpi.h>


ParticleStore store;

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    store = ParticleStore(parts, n, size);
    MPI_Barrier(MPI_COMM_WORLD);
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{
    store.compute_forces();
    store.move_particles();
    store.communicate_particles();
}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{
    store.gather_particles(parts);
}



#include "common.h"
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include <math.h>
#include <assert.h>
#include <mpi.h>

int getmyrank(MPI_Comm comm)
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    return myrank;
}

int getnprocs(MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    return nprocs;
}

#define MPI_ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "[rank %d] assertion (%s) failed at %s:%d\n", getmyrank(MPI_COMM_WORLD), #cond, __FILE__, __LINE__); \
        MPI_Barrier(MPI_COMM_WORLD); \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    } \
} while (0);

void apply_force(particle_t& target, const particle_t& ref)
{
    double dx = ref.x - target.x;
    double dy = ref.y - target.y;
    double r2 = dx * dx + dy * dy;

    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    double coef = (1 - cutoff / r) / r2 / mass;
    target.ax += coef * dx;
    target.ay += coef * dy;
}

void move_particle(particle_t& p, double size)
{
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    while (p.x < 0 || p.x > size)
    {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size)
    {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);
    MPI_ASSERT(myrank == rank && nprocs == procs);
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{

}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{

}



#ifndef PARTICLE_STORE_H_
#define PARTICLE_STORE_H_

#include "common.h"
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

typedef struct { particle_t p; int id; } named_particle_t;

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

template <class T>
std::string vector_string(const std::vector<T>& v)
{
    std::ostringstream os;
    os << "{";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, ","));
    os << "}";
    return os.str();
}

#define MPI_ASSERT(cond) do { \
    if (!(cond)) { \
        int _myrank_; \
        MPI_Comm_rank(MPI_COMM_WORLD, &_myrank_); \
        fprintf(stderr, "[rank %d] assertion (%s) failed at %s:%d\n", _myrank_, #cond, __FILE__, __LINE__); \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    } \
} while (0);

class ParticleStore
{
public:
    ParticleStore(const particle_t *particles, int numparts, double size);

    void compute_forces();
    void move_particles();
    void communicate_particles();

    int mynumparts() const;
    int get_particle_rank(const particle_t& p) const;
    void gather_particles(particle_t *parts) const;

    std::vector<named_particle_t> get_my_particles() const;
    std::vector<named_particle_t> gather_neighbor_particles() const;

    void print_info() const;

private:
    MPI_Comm gridcomm;
    int procdim;
    int nneighbors;
    int numparts;
    int myrank, nprocs;
    double procwidth;
    double gridsize;
    std::vector<named_particle_t> parts;
    MPI_Datatype NAMED_PARTICLE;

    void gather_named_particles(std::vector<named_particle_t>& particles, int root) const;
    void disjoint_partition_check() const;
};

ParticleStore::ParticleStore(const particle_t *particles, int numparts, double size)
    : numparts(numparts), gridsize(size)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    procdim = static_cast<int>(std::sqrt(nprocs+0.0));
    procwidth = size / procdim;

    MPI_ASSERT(procdim*procdim == nprocs);
    MPI_ASSERT(procwidth >= 2*cutoff + 1e-16);

    for (int i = 0; i < numparts; ++i)
        if (get_particle_rank(particles[i]) == myrank)
            parts.push_back({particles[i], i});

    disjoint_partition_check();

    std::vector<int> dests, weights;
    int reorder = 0;
    int procrow = myrank / procdim;
    int proccol = myrank % procdim;

    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        {
            int dest = myrank + (dx*procdim + dy);
            int row = procrow + dx;
            int col = proccol + dy;

            if (row >= 0 && row < procdim && col >= 0 && col < procdim)
            {
                dests.push_back(dest);
                weights.push_back(1);
            }
        }

    nneighbors = static_cast<int>(dests.size());
    MPI_Dist_graph_create(MPI_COMM_WORLD, 1, &myrank, &nneighbors, dests.data(), weights.data(), MPI_INFO_NULL, reorder, &gridcomm);

    int indegree, outdegree, weighted;
    MPI_Dist_graph_neighbors_count(gridcomm, &indegree, &outdegree, &weighted);
    MPI_ASSERT(indegree == outdegree && indegree == nneighbors);

    int nitems = 2;
    int blocklens[2] = {1,1};
    MPI_Datatype types[2] = {PARTICLE, MPI_INT};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(named_particle_t, p);
    offsets[1] = offsetof(named_particle_t, id);
    MPI_Type_create_struct(nitems, blocklens, offsets, types, &NAMED_PARTICLE);
    MPI_Type_commit(&NAMED_PARTICLE);
}

int ParticleStore::mynumparts() const
{
    return parts.size();
}

int ParticleStore::get_particle_rank(const particle_t& p) const
{
    int rowid = static_cast<int>(p.x / procwidth);
    int colid = static_cast<int>(p.y / procwidth);

    return rowid*procdim + colid;
}

void ParticleStore::disjoint_partition_check() const
{
    #ifndef NDEBUG
    std::vector<named_particle_t> allparts;
    gather_named_particles(allparts, 0);

    if (myrank == 0)
    {
        int k = 0;
        std::vector<int> allids(allparts.size());
        std::transform(allparts.begin(), allparts.end(), allids.begin(), [](named_particle_t item) { return item.id; });
        std::sort(allids.begin(), allids.end());
        std::vector<int> cmp(numparts);
        std::generate(cmp.begin(), cmp.end(), [&]() { return k++; });
        MPI_ASSERT(cmp == allids);
    }
    #endif
}

void ParticleStore::print_info() const
{
    if (myrank == 0)
    {
        fprintf(stderr, "procdim=%d\n", procdim);
        fprintf(stderr, "numparts=%d\n", numparts);
        fprintf(stderr, "procwidth=%.4f\n", procwidth);
    }

    for (int i = 0; i < nprocs; ++i)
    {
        if (myrank == i)
        {
            int procrow = myrank / procdim;
            int proccol = myrank % procdim;
            fprintf(stderr, "P(%d, %d) currently stores %d particles and has %d neighbor processors\n", procrow, proccol, mynumparts(), nneighbors);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void ParticleStore::gather_named_particles(std::vector<named_particle_t>& particles, int root) const
{
    particles.clear();

    int totcount;
    int mycount = mynumparts();

    std::vector<int> recvcounts, displs;

    if (myrank == root)
    {
        recvcounts.resize(nprocs);
        displs.resize(nprocs);
    }

    MPI_Gather(&mycount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, root, MPI_COMM_WORLD);

    if (myrank == root)
    {
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        particles.resize(numparts);
    }

    MPI_Gatherv(parts.data(), mycount, NAMED_PARTICLE, particles.data(), recvcounts.data(), displs.data(), NAMED_PARTICLE, root, MPI_COMM_WORLD);
}

void ParticleStore::gather_particles(particle_t *parts) const
{
    std::vector<named_particle_t> allparts;
    gather_named_particles(allparts, 0);

    if (myrank == 0)
        for (int i = 0; i < numparts; ++i)
            parts[allparts[i].id] = allparts[i].p;

    MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<named_particle_t> ParticleStore::get_my_particles() const
{
    return parts;
}

std::vector<named_particle_t> ParticleStore::gather_neighbor_particles() const
{
    int mycount = mynumparts();
    std::vector<int> neighbor_counts(nneighbors);

    MPI_Neighbor_allgather(&mycount, 1, MPI_INT, neighbor_counts.data(), 1, MPI_INT, gridcomm);

    std::vector<int> displs(nneighbors);
    displs[0] = 0;
    std::partial_sum(neighbor_counts.begin(), neighbor_counts.end()-1, displs.begin()+1);

    int totrecv = neighbor_counts.back() + displs.back();
    std::vector<named_particle_t> recvbuf(totrecv);

    MPI_Neighbor_allgatherv(parts.data(), mycount, NAMED_PARTICLE, recvbuf.data(), neighbor_counts.data(), displs.data(), NAMED_PARTICLE, gridcomm);

    return recvbuf;
}

void ParticleStore::communicate_particles()
{
    disjoint_partition_check();

    std::vector<named_particle_t> named_neighbors = gather_neighbor_particles();

    parts.clear();

    for (auto it = named_neighbors.begin(); it != named_neighbors.end(); ++it)
        if (get_particle_rank(it->p) == myrank)
            parts.push_back(*it);
}

void ParticleStore::compute_forces()
{
    disjoint_partition_check();

    for (auto it = parts.begin(); it != parts.end(); ++it)
        it->p.ax = it->p.ay = 0;

    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<named_particle_t> refparts = gather_neighbor_particles();

    for (auto it1 = parts.begin(); it1 != parts.end(); ++it1)
        for (auto it2 = refparts.begin(); it2 != refparts.end(); ++it2)
            apply_force(it1->p, it2->p);

    MPI_Barrier(MPI_COMM_WORLD);
}

void ParticleStore::move_particles()
{
    disjoint_partition_check();

    for (auto it = parts.begin(); it != parts.end(); ++it)
        move_particle(it->p, gridsize);

    MPI_Barrier(MPI_COMM_WORLD);
}

#endif

#ifndef PARTICLE_GRID_H_
#define PARTICLE_GRID_H_

#include "common.h"
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <memory>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#define MPI_ASSERT(cond) do {                                                                               \
    if (!(cond)) {                                                                                          \
        int _myrank_;                                                                                       \
        MPI_Comm_rank(MPI_COMM_WORLD, &_myrank_);                                                           \
        fprintf(stderr, "[rank %d] assertion (%s) failed at %s:%d\n", _myrank_, #cond, __FILE__, __LINE__); \
        MPI_Abort(MPI_COMM_WORLD, 1);                                                                       \
    }                                                                                                       \
} while (0);

typedef struct { particle_t p; int id; } named_particle_t;

class ParticleGrid
{
public:
    ParticleGrid(const particle_t *parts, int n, double size);

    void compute_forces();
    void move_particles();
    void update_grid();
    void gather_particles(particle_t *parts);

    int binsperrow() const { return (bindim+numrowprocs-1) / numrowprocs; }
    int binspercol() const { return (bindim+numcolprocs-1) / numcolprocs; }

    int mynumbinrows() const { return myrowproc != numrowprocs-1? binsperrow() : bindim - (numrowprocs-1)*binsperrow(); }
    int mynumbincols() const { return mycolproc != numcolprocs-1? binspercol() : bindim - (numcolprocs-1)*binspercol(); }

    double rowprocwidth() const { return binwidth * binsperrow(); }
    double colprocwidth() const { return binwidth * binspercol(); }

private:
    int procdim, bindim;
    int numrowprocs, numcolprocs;
    int myrowproc, mycolproc;
    double binwidth;
    int numparts;
    double gridsize;
    int nneighbors;
    std::vector<named_particle_t> myparts;
    MPI_Comm gridcomm;
    MPI_Datatype NAMED_PARTICLE;

    void apply_force(particle_t& target, const particle_t& ref);
    void move_particle(particle_t& p);
    int get_particle_rank(const particle_t& p);
    std::vector<named_particle_t> gather_neighbor_particles();
};

ParticleGrid::ParticleGrid(const particle_t *parts, int n, double size)
    : numparts(n), gridsize(size)
{
    int myrank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    procdim = static_cast<int>(std::sqrt(nprocs+0.0));
    numrowprocs = numcolprocs = procdim;

    while (numrowprocs*numcolprocs <= nprocs)
        numcolprocs++;
    numcolprocs--;

    myrowproc = myrank / numcolprocs;
    mycolproc = myrank % numcolprocs;

    int binsperrow, binspercol;
    bindim = static_cast<int>(std::min(gridsize / (cutoff + 1e-16), static_cast<double>(1<<12)));
    binwidth = gridsize / bindim;

    if (myrank == 0) fprintf(stderr, "Running with %d processor rows and %d processor columns (%d processors unused)\n", numrowprocs, numcolprocs, nprocs - numrowprocs*numcolprocs);

    MPI_ASSERT(std::min(rowprocwidth(), colprocwidth()) >= 2*cutoff + 1e-16);

    for (int i = 0; i < numparts; ++i)
        if (get_particle_rank(parts[i]) == myrank)
            myparts.push_back({parts[i], i});

    std::vector<int> dests, weights;

    if (myrank < numrowprocs*numcolprocs)
    {
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
            {
                int dest = myrank + (dx*numcolprocs + dy);
                int row = myrowproc + dx;
                int col = mycolproc + dy;

                if (row >= 0 && row < numrowprocs && col >= 0 && col < numcolprocs)
                {
                    dests.push_back(dest);
                    weights.push_back(1);
                }
            }
    }

    int reorder = 0; /* TODO: Will want to make reorder=1 for improved performance, but will need to be careful! */
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
    MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<named_particle_t> ParticleGrid::gather_neighbor_particles()
{
    if (nneighbors == 0)
    {
        return std::vector<named_particle_t>();
    }

    int mycount = static_cast<int>(myparts.size());
    std::vector<int> recvcounts(nneighbors);
    std::vector<int> displs(nneighbors);

    MPI_Neighbor_allgather(&mycount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, gridcomm);

    displs[0] = 0;
    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    int totrecv = recvcounts.back() + displs.back();
    std::vector<named_particle_t> recvbuf(totrecv);

    MPI_Neighbor_allgatherv(myparts.data(), mycount, NAMED_PARTICLE, recvbuf.data(), recvcounts.data(), displs.data(), NAMED_PARTICLE, gridcomm);

    return recvbuf;
}

int ParticleGrid::get_particle_rank(const particle_t& p)
{
    int rowid = static_cast<int>(p.x / rowprocwidth());
    int colid = static_cast<int>(p.y / colprocwidth());
    return rowid*numcolprocs + colid;
}

void ParticleGrid::compute_forces()
{
    for (auto it = myparts.begin(); it != myparts.end(); ++it)
        it->p.ax = it->p.ay = 0;

    auto neighparts = gather_neighbor_particles();

    for (auto it1 = myparts.begin(); it1 != myparts.end(); ++it1)
        for (auto it2 = neighparts.begin(); it2 != neighparts.end(); ++it2)
            apply_force(it1->p, it2->p);
}

void ParticleGrid::update_grid()
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    auto neighparts = gather_neighbor_particles();
    myparts.clear();

    for (auto it = neighparts.begin(); it != neighparts.end(); ++it)
        if (get_particle_rank(it->p) == myrank)
            myparts.push_back(*it);
}

void ParticleGrid::gather_particles(particle_t *parts)
{
    int myrank, nprocs;
    int totcount, mycount;
    std::vector<named_particle_t> allparts;
    std::vector<int> recvcounts, displs;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myrank == 0)
    {
        recvcounts.resize(nprocs);
        displs.resize(nprocs);
    }

    mycount = static_cast<int>(myparts.size());
    MPI_Gather(&mycount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        allparts.resize(numparts);
    }

    MPI_Gatherv(myparts.data(), mycount, NAMED_PARTICLE, allparts.data(), recvcounts.data(), displs.data(), NAMED_PARTICLE, 0, MPI_COMM_WORLD);

    if (myrank == 0)
        for (int i = 0; i < numparts; ++i)
            parts[allparts[i].id] = allparts[i].p;
}

void ParticleGrid::move_particles()
{
    for (auto it = myparts.begin(); it != myparts.end(); ++it)
        move_particle(it->p);
}

void ParticleGrid::apply_force(particle_t& target, const particle_t& ref)
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

void ParticleGrid::move_particle(particle_t& p)
{
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    while (p.x < 0 || p.x > gridsize)
    {
        p.x = p.x < 0 ? -p.x : 2 * gridsize - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > gridsize)
    {
        p.y = p.y < 0 ? -p.y : 2 * gridsize - p.y;
        p.vy = -p.vy;
    }
}

#endif

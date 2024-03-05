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

#define MPI_ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "[rank %d] assertion (%s) failed at %s:%d\n", getmyrank(MPI_COMM_WORLD), #cond, __FILE__, __LINE__); \
        MPI_Abort(MPI_COMM_WORLD, 1); \
    } \
} while (0);

class ParticleStore
{
public:
    ParticleStore();
    ParticleStore(const ParticleStore& rhs);
    ParticleStore(const particle_t *parts, int numparts, double size);
    ParticleStore& operator=(ParticleStore other);

    int myprocrow() const;
    int myproccol() const;

    int get_particle_rank(const particle_t& p) const;

    friend void swap(ParticleStore& lhs, ParticleStore& rhs);

    static int getmyrank(MPI_Comm comm);
    static int getnprocs(MPI_Comm comm);

    void print_info() const;

private:
    MPI_Comm gridcomm;
    int procdim;
    int nneighbors;
    int numparts;
    double procwidth;
    std::vector<int> myids;
    std::vector<particle_t> myparts;
    bool initialized;

    void disjoint_partition_check();
};

ParticleStore::ParticleStore()
    : initialized(false) {}

ParticleStore::ParticleStore(const ParticleStore& rhs)
    : gridcomm(rhs.gridcomm),
      procdim(rhs.procdim),
      nneighbors(rhs.nneighbors),
      numparts(rhs.numparts),
      procwidth(rhs.procwidth),
      myids(rhs.myids),
      myparts(rhs.myparts),
      initialized(rhs.initialized) {}

ParticleStore::ParticleStore(const particle_t *parts, int numparts, double size)
    : numparts(numparts), initialized(true)
{
    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    procdim = static_cast<int>(std::sqrt(nprocs+0.0));
    procwidth = size / procdim;

    MPI_ASSERT(procdim*procdim == nprocs);
    MPI_ASSERT(procwidth >= 2*cutoff + 1e-16);

    for (int i = 0; i < numparts; ++i)
        if (get_particle_rank(parts[i]) == myrank)
        {
            myids.push_back(i);
            myparts.push_back(parts[i]);
        }

    disjoint_partition_check();

    std::vector<int> dests, weights;
    int reorder = 0;
    int procrow = myprocrow();
    int proccol = myproccol();

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
}

ParticleStore& ParticleStore::operator=(ParticleStore other)
{
    swap(*this, other);
    return *this;
}

void swap(ParticleStore& lhs, ParticleStore& rhs)
{
    std::swap(lhs.gridcomm, rhs.gridcomm);
    std::swap(lhs.procdim, rhs.procdim);
    std::swap(lhs.numparts, rhs.numparts);
    std::swap(lhs.nneighbors, rhs.nneighbors);
    std::swap(lhs.procwidth, rhs.procwidth);
    std::swap(lhs.myids, rhs.myids);
    std::swap(lhs.myparts, rhs.myparts);
    std::swap(lhs.initialized, rhs.initialized);
}

int ParticleStore::getmyrank(MPI_Comm comm)
{
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    return myrank;
}

int ParticleStore::getnprocs(MPI_Comm comm)
{
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    return nprocs;
}

int ParticleStore::myprocrow() const
{
    return getmyrank(MPI_COMM_WORLD) / procdim;
}

int ParticleStore::myproccol() const
{
    return getmyrank(MPI_COMM_WORLD) % procdim;
}

int ParticleStore::get_particle_rank(const particle_t& p) const
{
    MPI_ASSERT(initialized);

    int rowid = static_cast<int>(p.x / procwidth);
    int colid = static_cast<int>(p.y / procwidth);

    return rowid*procdim + colid;
}

void ParticleStore::disjoint_partition_check()
{
    MPI_ASSERT(initialized);

    #ifdef NDEBUG
    return;
    #endif

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    int totcount;
    int mycount = static_cast<int>(myids.size());
    MPI_Reduce(&mycount, &totcount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0) MPI_ASSERT(totcount == numparts);

    std::vector<int> recvcounts, displs;
    std::vector<int> allids;

    if (myrank == 0)
    {
        recvcounts.resize(nprocs);
        displs.resize(nprocs);
    }

    MPI_Gather(&mycount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        allids.resize(totcount);
    }

    MPI_Gatherv(myids.data(), mycount, MPI_INT, allids.data(), recvcounts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        int k = 0;
        std::sort(allids.begin(), allids.end());
        std::vector<int> cmp(totcount);
        std::generate(cmp.begin(), cmp.end(), [&]() { return k++; });
        MPI_ASSERT(cmp == allids);

    }
}

void ParticleStore::print_info() const
{
    MPI_ASSERT(initialized);

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

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
            fprintf(stderr, "P(%d, %d) currently stores %d particles and has %d neighbor processors\n", myprocrow(), myproccol(), static_cast<int>(myids.size()), nneighbors);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}


#endif

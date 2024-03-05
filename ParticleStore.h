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
    void gather_particles(particle_t *parts) const;

    std::vector<named_particle_t> get_my_named_particles() const;
    std::vector<named_particle_t> gather_neighbor_particles() const;

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

    void gather_items(std::vector<int>& allids, std::vector<particle_t>& allparts, int root) const;
    void disjoint_partition_check() const;
    void sanity_check() const;
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

void ParticleStore::sanity_check() const
{
    MPI_ASSERT(initialized);
    MPI_ASSERT(myids.size() == myparts.size());
}

int ParticleStore::get_particle_rank(const particle_t& p) const
{
    sanity_check();

    int rowid = static_cast<int>(p.x / procwidth);
    int colid = static_cast<int>(p.y / procwidth);

    return rowid*procdim + colid;
}

void ParticleStore::disjoint_partition_check() const
{
    sanity_check();

    #ifdef NDEBUG
    return;
    #endif

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    std::vector<int> allids;
    std::vector<particle_t> allparts;

    gather_items(allids, allparts, 0);

    if (myrank == 0)
    {
        int k = 0;
        std::sort(allids.begin(), allids.end());
        std::vector<int> cmp(numparts);
        std::generate(cmp.begin(), cmp.end(), [&]() { return k++; });
        MPI_ASSERT(cmp == allids);
    }
}

void ParticleStore::print_info() const
{
    sanity_check();

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

void ParticleStore::gather_items(std::vector<int>& allids, std::vector<particle_t>& allparts, int root) const
{
    sanity_check();

    allids.clear();
    allparts.clear();

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    int totcount;
    int mycount = static_cast<int>(myids.size());
    MPI_Reduce(&mycount, &totcount, 1, MPI_INT, MPI_SUM, root, MPI_COMM_WORLD);
    if (myrank == root) MPI_ASSERT(totcount == numparts);

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
        allids.resize(numparts);
        allparts.resize(numparts);
    }

    MPI_Gatherv(myids.data(), mycount, MPI_INT, allids.data(), recvcounts.data(), displs.data(), MPI_INT, root, MPI_COMM_WORLD);
    MPI_Gatherv(myparts.data(), mycount, PARTICLE, allparts.data(), recvcounts.data(), displs.data(), PARTICLE, root, MPI_COMM_WORLD);
}

void ParticleStore::gather_particles(particle_t *parts) const
{
    sanity_check();

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    std::vector<int> allids;
    std::vector<particle_t> allparts;

    gather_items(allids, allparts, 0);

    MPI_ASSERT(allids.size() == numparts && allparts.size() == numparts);

    for (int i = 0; i < numparts; ++i)
    {
        parts[allids[i]] = allparts[i];
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<named_particle_t> ParticleStore::get_my_named_particles() const
{
    sanity_check();

    std::vector<named_particle_t> parts;
    parts.reserve(myids.size());

    for (int i = 0; i < myids.size(); ++i)
    {
        named_particle_t p = {myparts[i], myids[i]};
        parts.push_back(p);
    }

    return parts;
}

std::vector<named_particle_t> ParticleStore::gather_neighbor_particles() const
{
    sanity_check();

    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);

    int mycount = static_cast<int>(myparts.size());
    std::vector<int> neighbor_counts(nneighbors);

    MPI_Neighbor_allgather(&mycount, 1, MPI_INT, neighbor_counts.data(), 1, MPI_INT, gridcomm);

    std::vector<int> displs(nneighbors);
    displs[0] = 0;
    std::partial_sum(neighbor_counts.begin(), neighbor_counts.end()-1, displs.begin()+1);

    int totrecv = neighbor_counts.back() + displs.back();


    std::vector<named_particle_t> sendbuf = get_my_named_particles();
    std::vector<named_particle_t> recvbuf(totrecv);

    MPI_Datatype NAMED_PARTICLE;
    int nitems = 2;
    int blocklens[2] = {1,1};
    MPI_Datatype types[2] = {PARTICLE, MPI_INT};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(named_particle_t, p);
    offsets[1] = offsetof(named_particle_t, id);
    MPI_Type_create_struct(nitems, blocklens, offsets, types, &NAMED_PARTICLE);
    MPI_Type_commit(&NAMED_PARTICLE);

    MPI_Neighbor_allgatherv(myparts.data(), mycount, NAMED_PARTICLE, recvbuf.data(), neighbor_counts.data(), displs.data(), NAMED_PARTICLE, gridcomm);

    MPI_Type_free(&NAMED_PARTICLE);

    return recvbuf;
}


#endif

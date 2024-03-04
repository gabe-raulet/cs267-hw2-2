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

static int procdim, procrow, proccol;
static int numparts;
static double procwidth;
static std::vector<int> myids;
static std::vector<particle_t> myparts;

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
        MPI_Abort(MPI_COMM_WORLD, 1); \
    } \
} while (0);

template <class T>
std::string vector_string(const std::vector<T>& v)
{
    std::ostringstream os;
    os << "{";
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, ","));
    os << "}";
    return os.str();
}

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

int get_particle_rank(const particle_t& p)
{
    int rowid = static_cast<int>(p.x / procwidth);
    int colid = static_cast<int>(p.y / procwidth);

    return rowid*procdim + colid;
}

void disjoint_partition_check()
{
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

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    int myrank = getmyrank(MPI_COMM_WORLD);
    int nprocs = getnprocs(MPI_COMM_WORLD);
    MPI_ASSERT(myrank == rank && nprocs == procs);

    procdim = static_cast<int>(std::sqrt(nprocs+0.0));
    procwidth = size / procdim;
    procrow = myrank / procdim;
    proccol = myrank % procdim;
    numparts = n;

    MPI_ASSERT(procwidth >= 2*cutoff + 1e-16);
    MPI_ASSERT(procdim*procdim == nprocs);

    for (int i = 0; i < n; ++i)
        if (get_particle_rank(parts[i]) == myrank)
        {
            myids.push_back(i);
            myparts.push_back(parts[i]);
        }

    disjoint_partition_check();
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{

}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{

}



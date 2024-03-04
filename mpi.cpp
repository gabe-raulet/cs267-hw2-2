#include "common.h"
#include <list>
#include <vector>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include <mpi.h>

int numparts;
int procdim;
double gridsize;
double procwidth;
std::vector<int> myids;

MPI_Comm grid_comm;
int myprocrow;
int myproccol;
int nneighbors;

void get_comm_info(MPI_Comm comm, int *myrank, int *nprocs)
{
    MPI_Comm_rank(comm, myrank);
    MPI_Comm_size(comm, nprocs);
}

int get_particle_rank(const particle_t& p)
{
    int rowid = static_cast<int>(std::floor(p.x / procwidth));
    int colid = static_cast<int>(std::floor(p.y / procwidth));

    return rowid*procdim + colid;
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

void move_particle(particle_t& p)
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

void clear_my_acceleration(particle_t *parts)
{
    for (auto it = myids.begin(); it != myids.end(); ++it)
    {
        particle_t& p = parts[*it];
        p.ax = p.ay = 0;
    }
}

void apply_forces(particle_t *parts)
{
    int myrank, nprocs;
    get_comm_info(grid_comm, &myrank, &nprocs);

    std::vector<int> recvcounts(nneighbors);
    std::vector<int> displs(nneighbors);

    int sendcount = recvcounts[myrank];
    MPI_Neighbor_allgather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, grid_comm);

    displs[0] = 0;
    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    std::vector<particle_t> neighbors(recvcounts.back() + displs.back());
    MPI_Neighbor_allgatherv(myids.data(), sendcount, PARTICLE, neighbors.data(), recvcounts.data(), displs.data(), PARTICLE, grid_comm);

    for (auto tgtit = myids.begin(); tgtit != myids.end(); ++tgtit)
        for (auto refit = neighbors.begin(); refit != neighbors.end(); ++refit)
            apply_force(parts[*tgtit], *refit);
}

void move_particles(particle_t *parts)
{
    for (auto it = myids.begin(); it != myids.end(); ++it)
        move_particle(parts[*it]);
}

void reassign_particles(particle_t *parts)
{
    int myrank, nprocs;
    get_comm_info(MPI_COMM_WORLD, &myrank, &nprocs);
    std::vector<particle_t> psendbuf, precvbuf;
    std::vector<int> isendbuf, irecvbuf;
    std::vector<int> sendcounts(nprocs, 0);
    std::vector<int> sdispls(nprocs);

    for (auto it = myids.begin(); it != myids.end(); ++it)
        sendcounts[get_particle_rank(parts[*it])]++;

    psendbuf.reserve(sendcounts[myrank]);
    isendbuf.reserve(sendcounts[myrank]);

    for (auto it = myids.begin(); it != myids.end(); ++it)
    {
        int prank = get_particle_rank(parts[*it]);
        isendbuf.push_back(*it);
        psendbuf.push_back(parts[*it]);
    }

    std::vector<int> recvcounts(nprocs);
    std::vector<int> rdispls(nprocs);

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    sdispls[0] = rdispls[0] = 0;
    std::partial_sum(sendcounts.begin(), sendcounts.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, rdispls.begin()+1);

    irecvbuf.resize(recvcounts.back() + rdispls.back());
    precvbuf.resize(recvcounts.back() + rdispls.back());

    MPI_Alltoallv(isendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT,  irecvbuf.data(), recvcounts.data(), rdispls.data(), MPI_INT,  MPI_COMM_WORLD);
    MPI_Alltoallv(psendbuf.data(), sendcounts.data(), sdispls.data(), PARTICLE, precvbuf.data(), recvcounts.data(), rdispls.data(), PARTICLE, MPI_COMM_WORLD);

    myids.resize(recvcounts[myrank]);
    std::copy(irecvbuf.begin(), irecvbuf.end(), myids.begin());

    for (auto it = myids.begin(); it != myids.end(); ++it)
        parts[*it] = precvbuf[it - myids.begin()];
}

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    int myrank, nprocs;
    get_comm_info(MPI_COMM_WORLD, &myrank, &nprocs);
    assert(myrank == rank && nprocs == procs); /* sanity check */

    procdim = static_cast<int>(std::floor(std::sqrt(static_cast<double>(nprocs))));
    numparts = n;
    gridsize = size;
    procwidth = gridsize / procdim;
    assert(procdim*procdim == nprocs);
    assert(procwidth >= 2*cutoff + 1e-16);

    myprocrow = myrank / procdim;
    myproccol = myrank % procdim;

    std::vector<int> dests;
    int degree;

    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        {
            int dest = myrank + (dx*procdim + dy);
            int procrow = myprocrow + dx;
            int proccol = myproccol + dy;

            if (procrow >= 0 && procrow < procdim && proccol >= 0 && proccol < procdim)
                dests.push_back(dest);
        }

    degree = static_cast<int>(dests.size());
    std::vector<int> weights(degree, 1);
    MPI_Dist_graph_create(MPI_COMM_WORLD, 1, &myrank, &degree, dests.data(), weights.data(), MPI_INFO_NULL, 0, &grid_comm);

    int indegree, outdegree, weighted;
    MPI_Dist_graph_neighbors_count(grid_comm, &indegree, &outdegree, &weighted);
    assert(indegree == outdegree);
    nneighbors = indegree;

    for (int i = 0; i < numparts; ++i)
        if (get_particle_rank(parts[i]) == myrank)
            myids.push_back(i);
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{
    /*
     * Clear the acceleration vectors of all particles on my
     * processor's partition.
     */
    clear_my_acceleration(parts);

    /*
     * Apply forces to all particles.
     */
    apply_forces(parts);

    /*
     * Move particles based on their current velocity vector (this step
     * may cause some particles to move out of my processor's assigned
     * square of space.)
     */
    move_particles(parts);

    /*
     * Reassign lost particles (particles whose current position has become
     * incongruent with their currently assigned processor).
     */
    reassign_particles(parts);
}

template <class T>
void reorder_array(T *v, const std::vector<int>& indices)
{
    std::vector<T> reordered;
    reordered.reserve(indices.size());

    for (auto it = indices.begin(); it != indices.end(); ++it)
        reordered.push_back(v[*it]);

    std::copy(reordered.begin(), reordered.end(), v);
}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{
    int myrank, nprocs;
    get_comm_info(MPI_COMM_WORLD, &myrank, &nprocs);

    std::vector<particle_t> myparts;
    std::vector<int> myidvec(myids.begin(), myids.end());
    int sendcount = static_cast<int>(myidvec.size());

    myparts.reserve(sendcount);
    for (auto it = myidvec.begin(); it != myidvec.end(); ++it)
        myparts.push_back(parts[*it]);

    std::vector<int> idvec;
    std::vector<int> recvcounts;
    std::vector<int> displs;

    if (myrank == 0)
    {
        idvec.resize(numparts);
        recvcounts.resize(nprocs);
        displs.resize(nprocs);
    }

    MPI_Gather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        displs[0] = 0;
        std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);
        assert(displs.back() + recvcounts.back() == numparts);
    }

    MPI_Gatherv(myidvec.data(), sendcount, MPI_INT, idvec.data(), recvcounts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(myparts.data(), sendcount, PARTICLE, parts, recvcounts.data(), displs.data(), PARTICLE, 0, MPI_COMM_WORLD);

    if (myrank == 0)
    {
        reorder_array(parts, idvec);
    }
}



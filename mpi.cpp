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

int procdim;
int numrowprocs, numcolprocs;
int myrowproc, mycolproc;
double rowprocwidth, colprocwidth;
int numparts;
double gridsize;
int nneighbors;
std::vector<named_particle_t> myparts;
MPI_Comm gridcomm;
MPI_Datatype NAMED_PARTICLE;

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

int get_particle_rank(const particle_t& p)
{
    /*
     * The processor responsible for owning a particle is determined
     * by the x position coordinate (determining the processor row) and
     * the y position coordinate (determining the processor column).
     */
    int rowid = static_cast<int>(p.x / rowprocwidth);
    int colid = static_cast<int>(p.y / colprocwidth);

    return rowid*numcolprocs + colid;
}

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    int myrank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_ASSERT(rank == myrank && procs == nprocs);

    /*
     * Every processor shares the following global variables:
     *
     *  int numparts - total number of particles
     *  int procdim - number of processors in each processor row/column
     *  double gridsize - side length of square simulation space
     *  double procwidth - side length of each processor's assigned square space
     */
    numparts = n;
    gridsize = size;
    procdim = static_cast<int>(std::sqrt(nprocs+0.0));
    numrowprocs = numcolprocs = procdim;

    while (numrowprocs*numcolprocs <= nprocs)
        numcolprocs++;
    numcolprocs--;

    rowprocwidth = size / numrowprocs;
    colprocwidth = size / numcolprocs;
    myrowproc = myrank / numcolprocs;
    mycolproc = myrank % numcolprocs;

    if (myrank == 0) fprintf(stderr, "Running with %d processor rows and %d processor columns (%d processors unused)\n", numrowprocs, numcolprocs, nprocs - numrowprocs*numcolprocs);

    /*
     * Make sure that each processor square is not too small.
     */
    MPI_ASSERT(std::min(rowprocwidth, colprocwidth) >= 2*cutoff + 1e-16);

    /* Every processor has the following rank-local variables:
     *
     *  std::vector<named_particle_t> myparts - local set of particles
     *  int nneighbors - number of neighboring processors in the processor grid
     *
     * Since each processor starts off with a copy of all the particles,
     * each processor simply finds which ones it is responsible for storing
     * and adds it to its local storage @myparts. In order to keep track
     * of which particles are which, we store particles in a structure
     * @named_particle_t which couples the particle with its global index.
     */
    for (int i = 0; i < numparts; ++i)
        if (get_particle_rank(parts[i]) == myrank)
            myparts.push_back({parts[i], i});

    /*
     * Communication between neighboring processors is done using a communicator
     * with a virtual distributed graph topology. This topology is similar to a
     * a Cartesian topology, except that we include the diagonal relation. For example,
     * a 3-by-3 processor grid would have the following structure:
     *
     * P(0,0) has neighbors P(0,0), P(0,1), P(1,1), P(1,0)
     * P(0,1) has neighbors P(0,1), P(0,0), P(0,2), P(1,2), P(1,1), P(1,0)
     * P(0,2) has neighbors P(0,2), P(0,1), P(1,2), P(1,1)
     * P(1,0) has neighbors P(1,0), P(0,0), P(0,1), P(1,1)
     * P(1,1) has neighbors P(1,1), P(0,0), P(0,1), P(0,2), P(1,0), P(1,2), P(2,0), P(2,1), P(2,2)
     *
     * And so on ...
     *
     * Notice that each processor has itself as a neighbor to simplify computation.
     *
     * We build the distributed graph topology by constructing a graph whose vertices are
     * processors and whose edges represent neighbor processors.
     */
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

    /*
     * Make sure that the grid communicator returns the correct number of neighbors
     * for each processor. Note that our graph is perfectly symmetric.
     */
    int indegree, outdegree, weighted;
    MPI_Dist_graph_neighbors_count(gridcomm, &indegree, &outdegree, &weighted);
    MPI_ASSERT(indegree == outdegree && indegree == nneighbors);

    /*
     * In order to commmunicate particles coupled with their global ids, we
     * create an MPI_Datatype for the @named_particle_t structure.
     */
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

std::vector<named_particle_t> gather_neighbor_particles()
{
    /*
     * This function gathers all particles stored in neighboring
     * processors to this processor.
     */

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

void compute_forces()
{
    /*
     * Clear acceleration vectors.
     */
    for (auto it = myparts.begin(); it != myparts.end(); ++it)
        it->p.ax = it->p.ay = 0;

    /*
     * Gather neighboring particles and use them to compute the
     * new acceleration vectors for this processor's particles.
     */
    auto neighparts = gather_neighbor_particles();

    for (auto it1 = myparts.begin(); it1 != myparts.end(); ++it1)
        for (auto it2 = neighparts.begin(); it2 != neighparts.end(); ++it2)
            apply_force(it1->p, it2->p);
}

void move_particles()
{
    for (auto it = myparts.begin(); it != myparts.end(); ++it)
        move_particle(it->p);
}

void communicate_particles()
{
    /*
     * After particles have been moved, it is necessary to reassign particles
     * that may have crossed their storage boundary. This can quickly be done
     * by gathering all neighboring particles (which includes local particles since
     * the graph topology has self-loops for all vertices) and storing those
     * that should be local.
     *
     * Note that the correctness of this function relies on particles not
     * moving far enough to cross over a processor's square.
     */

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    auto neighparts = gather_neighbor_particles();
    myparts.clear();

    for (auto it = neighparts.begin(); it != neighparts.end(); ++it)
        if (get_particle_rank(it->p) == myrank)
            myparts.push_back(*it);
}

void gather_particles(particle_t *parts)
{
    /*
     * In order to get a global view of particles at the root rank, we first
     * need to gather all particles to the root. Because the order of particles
     * is likely to change significantly over the course of the simulation, we
     * will need to reorder them according to their global ids once they have
     * been gathered. We therefore use the @named_particle_t structure to
     * help reorder the particles in correct sorted order.
     */

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

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{
    compute_forces();
    move_particles();
    communicate_particles();
}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{
    gather_particles(parts);
}



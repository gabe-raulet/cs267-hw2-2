#include "common.h"
#include "ParticleGrid.h"

std::unique_ptr<ParticleGrid> grid;

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    grid.reset(new ParticleGrid(parts, n, size));
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{
    grid->compute_forces();
    grid->move_particles();
    grid->update_grid();
}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{
    grid->gather_particles(parts);
}



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

static int bindim;
static int numbins;
static int numparts;
static double binwidth;
static double gridsize;

static constexpr int maxbindim = 1<<12;
static constexpr int bincap = 16;
static int bins[maxbindim*maxbindim][bincap];
static int binsizes[maxbindim*maxbindim];
static std::vector<int> bintable;

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

void compute_forces(particle_t *parts)
{
    for (int i = 0; i < numparts; ++i)
        parts[i].ax = parts[i].ay = 0;

    for (int i = 0; i < numparts; ++i)
    {
        int target_binid = bintable[i];
        particle_t& target = parts[i];

        for (int dx = -bindim; dx <= bindim; dx += bindim)
            for (int dy = -1; dy <= 1; ++dy)
            {
                int ref_binid = target_binid + dx + dy;

                if (ref_binid >= 0 && ref_binid < numbins)
                {
                    int *refbin = bins[ref_binid];
                    for (int j = 0; j < binsizes[ref_binid]; ++j)
                        apply_force(target, parts[refbin[j]]);
                }
            }
    }
}

void move_particles(particle_t *parts)
{
    for (int i = 0; i < numparts; ++i)
        move_particle(parts[i]);
}

void update_bins(particle_t *parts)
{
    //std::fill(binsizes, binsizes + numbins, 0);
    for (int i = 0; i < numbins; ++i)
        binsizes[i] = 0;

    for (int i = 0; i < numparts; ++i)
    {
        int rowid = static_cast<int>(parts[i].x / binwidth);
        int colid = static_cast<int>(parts[i].y / binwidth);
        int binid = rowid*bindim + colid;
        int pos = binsizes[binid]++;
        if (pos >= bincap) fprintf(stderr, "ERROR!!\n");
        bins[binid][pos] = i;
        bintable[i] = binid;
    }
}

void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
{
    gridsize = size;
    numparts = n;
    bindim = std::min(static_cast<int>(gridsize / (cutoff + 1e-16)), maxbindim);
    numbins = bindim * bindim;
    binwidth = gridsize / bindim;
    bintable.resize(numparts);
    update_bins(parts);
}

void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
{
    compute_forces(parts);
    move_particles(parts);
    update_bins(parts);
}

void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
{

}

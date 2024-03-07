// #include "common.h"
// #include <cmath>
// #include <cstddef>
// #include <numeric>
// #include <vector>
// #include <memory>
// #include <iostream>

// typedef struct { particle_t p; int id; } named_particle_t;

// class ParticleGrid
// {
// public:
//     ParticleGrid(const particle_t *parts, int n, double size);

//     void compute_forces();
//     void move_particles();
//     void update_grid();
//     void gather_particles(particle_t *parts);
//     void print_time() {
//         if (myrank == 0) {
//             std::cout << "COMPTIME: " << computation_time << " seconds" << std::endl;
//             std::cout << "ATATIME: " << communication_ATA_time << " seconds" << std::endl;
//             std::cout << "SRTIME: " << communication_SR_time << " seconds" << std::endl;
//         }
//     }

// private:
//     int myrank, nprocs;
//     int numrowprocs, numcolprocs;
//     int myrowproc, mycolproc;
//     double rowsize, colsize;
//     int binrowdim, bincoldim;
//     double binrowwidth, bincolwidth;
    
//     int numparts;
//     double gridsize, computation_time, communication_SR_time, communication_ATA_time;
//     double computation_st, communication_SR_st, communication_ATA_st;
//     std::vector<std::vector<named_particle_t>> mybins, tmpbins;
//     std::vector<std::vector<named_particle_t>> tomove;
//     std::vector<int> sendcounts, recvcounts, senddispls, recvdispls;
//     std::vector<named_particle_t> sendtoothers, recvfromothers;
//     MPI_Datatype NAMED_PARTICLE;

//     void apply_force(particle_t& target, const particle_t& ref);
//     void move_particle(particle_t& p);
//     int get_particle_rank(const particle_t& p);
//     int get_particle_bin(const particle_t& p);
//     void commit_named_particle_type();
//     void add_particle(const particle_t& p, int id);
//     inline int get_bin_id(const int& row, const int& col) const {
//         return row * (bincoldim + 2) + col;
//     }
//     inline void communication_SR_start() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         communication_SR_st = MPI_Wtime();
//     }
//     inline void communication_SR_end() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         communication_SR_time += MPI_Wtime() - communication_SR_st;
//     }
//     inline void communication_ATA_start() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         communication_ATA_st = MPI_Wtime();
//     }
//     inline void communication_ATA_end() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         communication_ATA_time += MPI_Wtime() - communication_ATA_st;
//     }
//     inline void computation_start() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         computation_st = MPI_Wtime();
//     }
//     inline void computation_end() {
//         MPI_Barrier(MPI_COMM_WORLD);
//         computation_time += MPI_Wtime() - computation_st;
//     }
// };

// ParticleGrid::ParticleGrid(const particle_t *parts, int n, double size)
//     : numparts(n), gridsize(size)
// {
//     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//     computation_time = communication_SR_time = communication_ATA_time = 0;

//     numrowprocs = (int)floor(std::sqrt(nprocs+0.5));
//     numcolprocs = nprocs / numrowprocs;

//     myrowproc = myrank / numcolprocs;
//     mycolproc = myrank % numcolprocs;

//     rowsize = size / numrowprocs;
//     colsize = size / numcolprocs;

//     binrowdim = (int)floor(rowsize / (cutoff + 1e-16));
//     bincoldim = (int)floor(colsize / (cutoff + 1e-16));

//     binrowwidth = rowsize / binrowdim;
//     bincolwidth = colsize / bincoldim;

//     mybins.resize((binrowdim + 2) * (bincoldim + 2));
//     tmpbins.resize((binrowdim + 2) * (bincoldim + 2));
//     tomove.resize(nprocs);
//     sendcounts.resize(nprocs);
//     recvcounts.resize(nprocs);
//     senddispls.resize(nprocs);
//     recvdispls.resize(nprocs);

//     for (int i = 0; i < numparts; ++i) {
//         if(get_particle_rank(parts[i]) == myrank)
//             add_particle(parts[i], i);
//     }   

//     commit_named_particle_type();
// }

// void ParticleGrid::commit_named_particle_type()
// {
//     int nitems = 2;
//     int blocklens[2] = {1,1};
//     MPI_Datatype types[2] = {PARTICLE, MPI_INT};
//     MPI_Aint offsets[2];
//     offsets[0] = offsetof(named_particle_t, p);
//     offsets[1] = offsetof(named_particle_t, id);
//     MPI_Type_create_struct(nitems, blocklens, offsets, types, &NAMED_PARTICLE);
//     MPI_Type_commit(&NAMED_PARTICLE);
// }


// int ParticleGrid::get_particle_bin(const particle_t& p)
// {
//     double x = p.x - rowsize*myrowproc;
//     double y = p.y - colsize*mycolproc;

//     int rowid = static_cast<int>(x / binrowwidth);
//     int colid = static_cast<int>(y / bincolwidth);

//     return get_bin_id(rowid + 1, colid + 1);
// }

// int ParticleGrid::get_particle_rank(const particle_t& p)
// {
//     int rowid = static_cast<int>(p.x / rowsize);
//     int colid = static_cast<int>(p.y / colsize);
//     return rowid*numcolprocs + colid;
// }

// void ParticleGrid::compute_forces()
// {
//     computation_start();
//     std::vector<named_particle_t> toneighbor[9], fromneighbor[9];
//     std::vector<MPI_Request> requests;
//     int toneighborsize[9], fromneighborsize[9];
//     toneighbor[0] = mybins[get_bin_id(1, 1)];
//     {
//         int cnt = 0;
//         for (int j = 1;j <= bincoldim; ++j)
//             cnt += mybins[get_bin_id(1, j)].size();
//         toneighbor[1].reserve(cnt);
//         for (int j = 1;j <= bincoldim; ++j) 
//             for (auto &it : mybins[get_bin_id(1, j)])
//                 toneighbor[1].emplace_back(it);
//     }
//     toneighbor[2] = mybins[get_bin_id(1, bincoldim)];
//     {
//         int cnt = 0;
//         for (int i = 1;i <= binrowdim; ++i)
//             cnt += mybins[get_bin_id(i, 1)].size();
//         toneighbor[3].reserve(cnt);
//         for (int i = 1;i <= binrowdim; ++i)
//             for (auto &it : mybins[get_bin_id(i, 1)])
//                 toneighbor[3].emplace_back(it);
//     }
//     {
//         int cnt = 0;
//         for (int i = 1;i <= binrowdim; ++i)
//             cnt += mybins[get_bin_id(i, bincoldim)].size();
//         toneighbor[5].reserve(cnt);
//         for (int i = 1;i <= binrowdim; ++i)
//             for (auto &it : mybins[get_bin_id(i, bincoldim)])
//                 toneighbor[5].emplace_back(it);
//     }
//     toneighbor[6] = mybins[get_bin_id(binrowdim, 1)];
//     {
//         int cnt = 0;
//         for (int j = 1;j <= bincoldim; ++j)
//             cnt += mybins[get_bin_id(binrowdim, j)].size();
//         toneighbor[7].reserve(cnt);
//         for (int j = 1;j <= bincoldim; ++j) 
//             for (auto &it : mybins[get_bin_id(binrowdim, j)])
//                 toneighbor[7].emplace_back(it);
//     }
//     toneighbor[8] = mybins[get_bin_id(binrowdim, bincoldim)];
//     for (int i=0;i<9;++i)
//         toneighborsize[i] = toneighbor[i].size();
//     computation_end();
//     communication_SR_start();

//     for (int di = -1; di <= 1; ++di) {
//         for (int dj = -1; dj <= 1; ++dj) {
//             if (di == 0 && dj == 0)   continue;
//             if (myrowproc >= numrowprocs || mycolproc >= numcolprocs)   continue;
//             int nr = myrowproc + di, nc = mycolproc + dj;
//             if (nr >= 0 && nr < numrowprocs && nc >= 0 && nc < numcolprocs) {
//                 int direction = (di + 1) * 3 + (dj + 1);
//                 requests.push_back(MPI_Request());
//                 MPI_Isend(toneighborsize + direction, 1, MPI_INT, nr * numcolprocs + nc, 0, MPI_COMM_WORLD, &requests.back());
//                 requests.push_back(MPI_Request());
//                 MPI_Irecv(fromneighborsize + direction, 1, MPI_INT, nr * numcolprocs + nc, 0, MPI_COMM_WORLD, &requests.back());
//             }
//         }
//     }
//     MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
//     requests.clear();

//     for (int di = -1; di <= 1; ++di) {
//         for (int dj = -1; dj <= 1; ++dj) {
//             if (di == 0 && dj == 0)   continue;
//             if (myrowproc >= numrowprocs || mycolproc >= numcolprocs)   continue;
//             int nr = myrowproc + di, nc = mycolproc + dj;
//             if (nr >= 0 && nr < numrowprocs && nc >= 0 && nc < numcolprocs) {
//                 int direction = (di + 1) * 3 + (dj + 1);
//                 requests.push_back(MPI_Request());
//                 MPI_Isend(toneighbor[direction].data(), toneighborsize[direction], NAMED_PARTICLE, nr * numcolprocs + nc, 0, MPI_COMM_WORLD, &requests.back());
//                 fromneighbor[direction].resize(fromneighborsize[direction]);
//                 requests.push_back(MPI_Request());
//                 MPI_Irecv(fromneighbor[direction].data(), fromneighborsize[direction], NAMED_PARTICLE, nr * numcolprocs + nc, 0, MPI_COMM_WORLD, &requests.back());
//             }
//         }
//     }
//     MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
//     communication_SR_end();
//     computation_start();
    
//     for (int di = -1; di <= 1; ++di) {
//         for (int dj = -1; dj <= 1; ++dj) {
//             int direction = (di + 1) * 3 + (dj + 1);
//             for (auto &it : fromneighbor[direction]) {
//                 int row = (di == -1 ? 0 : (di == 1 ? binrowdim + 1 : static_cast<int>((it.p.x - rowsize*myrowproc) / binrowwidth) + 1));
//                 int col = (dj == -1 ? 0 : (dj == 1 ? bincoldim + 1 : static_cast<int>((it.p.y - colsize*mycolproc) / bincolwidth) + 1));
//                 int binid = get_bin_id(row, col);
//                 mybins[binid].emplace_back(it);
//             }
//         }
//     }

//     for (int i = 1; i <= binrowdim; ++i)
//         for (int j = 1; j <= bincoldim; ++j) {
//             int binid = get_bin_id(i, j);
//             for (auto &it1 : mybins[binid]) {
//                 it1.p.ax = it1.p.ay = 0;
//                 for (int di = -1; di <= 1; ++di) {
//                     for (int dj = -1; dj <= 1; ++dj) {
//                         int otherbinid = get_bin_id(i + di, j + dj);
//                         for (auto &it2 : mybins[otherbinid]) {
//                             apply_force(it1.p, it2.p);
//                         }
//                     }
//                 }
//             }
//         }
//     computation_end();
// }

// void ParticleGrid::add_particle(const particle_t& p, int id)
// {
//     int binid = get_particle_bin(p);
//     mybins[binid].push_back({p, id});
// }

// void ParticleGrid::update_grid()
// {
//     computation_start();
//     mybins.swap(tmpbins);  
//     for (int i=0;i<mybins.size();++i)
//         mybins[i].clear();
//     for (int i=0;i<nprocs;++i)
//         tomove[i].clear();

//     for (int i=1;i<=binrowdim;++i)
//         for (int j=1;j<=bincoldim;++j) {
//             int binid = get_bin_id(i, j);
//             for (auto &it : tmpbins[binid]) {
//                 int procid = get_particle_rank(it.p);
//                 if (procid == myrank)
//                     add_particle(it.p, it.id);
//                 else {
//                     tomove[procid].emplace_back(it);
//                 }
//             }
//         }

//     for(int i=0;i<nprocs;++i)
//         sendcounts[i] = tomove[i].size();
//     senddispls[0] = 0;
//     std::partial_sum(sendcounts.begin(), sendcounts.end()-1, senddispls.begin()+1);

//     computation_end();
//     communication_ATA_start();

//     MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

//     communication_ATA_end();
//     computation_start();
    
//     recvdispls[0] = 0;
//     std::partial_sum(recvcounts.begin(), recvcounts.end()-1, recvdispls.begin()+1);
    
//     recvfromothers.resize(recvdispls.back() + recvcounts.back());
//     sendtoothers.resize(senddispls.back() + sendcounts.back());

//     for(int i=0;i<nprocs;++i)
//         std::copy(tomove[i].begin(), tomove[i].end(), sendtoothers.begin() + senddispls[i]);

//     computation_end();
//     communication_ATA_start();
    
//     MPI_Alltoallv(sendtoothers.data(), sendcounts.data(), senddispls.data(), NAMED_PARTICLE, 
//                   recvfromothers.data(), recvcounts.data(), recvdispls.data(), NAMED_PARTICLE, MPI_COMM_WORLD);
    
//     communication_ATA_end();

//     computation_start();
//     for (auto &it : recvfromothers)
//         add_particle(it.p, it.id);
//     computation_end();
// }

// void ParticleGrid::gather_particles(particle_t *parts)
// {
//     std::vector<named_particle_t> allparts, myparts;

//     for (int i = 1; i <= binrowdim; ++i)
//         for (int j = 1; j <= bincoldim; ++j) {
//             int binid = get_bin_id(i, j);
//             for (auto &it : mybins[binid])
//                 myparts.emplace_back(it);
//         }

//     int mycount = myparts.size();
//     MPI_Gather(&mycount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

//     if (myrank == 0)
//     {
//         recvdispls[0] = 0;
//         std::partial_sum(recvcounts.begin(), recvcounts.end()-1, recvdispls.begin()+1);
//         allparts.resize(numparts);
//     }

//     MPI_Gatherv(myparts.data(), mycount, NAMED_PARTICLE, allparts.data(), recvcounts.data(), recvdispls.data(), NAMED_PARTICLE, 0, MPI_COMM_WORLD);

//     if (myrank == 0)
//         for (int i = 0; i < numparts; ++i)
//             parts[allparts[i].id] = allparts[i].p;
// }

// void ParticleGrid::move_particles()
// {
//     computation_start();
//     for (int i = 1; i <= binrowdim; ++i)
//         for (int j = 1; j <= bincoldim; ++j) {
//             int binid = get_bin_id(i, j);
//             for (auto &it : mybins[binid])
//                 move_particle(it.p);
//         }
//     computation_end();
// }

// void ParticleGrid::apply_force(particle_t& target, const particle_t& ref)
// {
//     double dx = ref.x - target.x;
//     double dy = ref.y - target.y;
//     double r2 = dx * dx + dy * dy;

//     if (r2 > cutoff * cutoff)
//         return;

//     r2 = fmax(r2, min_r * min_r);
//     double r = sqrt(r2);

//     double coef = (1 - cutoff / r) / r2 / mass;
//     target.ax += coef * dx;
//     target.ay += coef * dy;
// }

// void ParticleGrid::move_particle(particle_t& p)
// {
//     p.vx += p.ax * dt;
//     p.vy += p.ay * dt;
//     p.x += p.vx * dt;
//     p.y += p.vy * dt;

//     while (p.x < 0 || p.x > gridsize)
//     {
//         p.x = p.x < 0 ? -p.x : 2 * gridsize - p.x;
//         p.vx = -p.vx;
//     }

//     while (p.y < 0 || p.y > gridsize)
//     {
//         p.y = p.y < 0 ? -p.y : 2 * gridsize - p.y;
//         p.vy = -p.vy;
//     }
// }



// std::unique_ptr<ParticleGrid> grid;

// void init_simulation(particle_t *parts, int n, double size, int rank, int procs)
// {
//     grid.reset(new ParticleGrid(parts, n, size));
// }
// int cnt = 0;
// void simulate_one_step(particle_t *parts, int n, double size, int rank, int procs)
// {
//     grid->compute_forces();
//     grid->move_particles();
//     grid->update_grid();
//     if (++cnt == nsteps) {
//         grid -> print_time();
//     }
// }

// void gather_for_save(particle_t *parts, int n, double size, int rank, int procs)
// {
//     grid->gather_particles(parts);
// }

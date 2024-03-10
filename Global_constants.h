#ifndef NBODY_CONSTANTS_H
#define NBODY_CONSTANTS_H

const int NforDisc = 16384;
const int NforBulge = 16384;
const int NforHalo = 16384;
const int NforSatellite = 0;
const int NParticles = (NforDisc + NforHalo + NforBulge + NforSatellite);
static const int BLOCK_SIZE=512;
const float eps_for_disk = 0.08;
const float eps_for_halo = 0.4;
const float eps_for_bulge = 0.06;
const float eps_for_sat = 0.08;
#endif
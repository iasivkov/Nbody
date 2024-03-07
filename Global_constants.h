#ifndef NBODY_CONSTANTS_H
#define NBODY_CONSTANTS_H

const int NforDisc = 16384;
const int NforBulge = 16384;
const int NforHalo = 16384;
const int NforSatellite =0;
const int NParticles = (NforDisc + NforHalo + NforBulge + NforSatellite);
static const int BLOCK_SIZE=512;
#endif
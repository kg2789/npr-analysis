#ifndef __Zq_H__

#define __Zq_H__

#include "data.h"
#include "metadata.h"

// Wavefunction renormalization schemes
enum class WavefunctionScheme {
    GAMMA_MU,
    QSLASH,
};

ComplexD Zq_from_ZV(TrajectoryData &data, Metadata &meta, double ZV,
        WavefunctionScheme scheme);
ComplexD Zq_from_ZA(TrajectoryData &data, Metadata &meta, double ZA,
        WavefunctionScheme scheme);


ComplexD Zq_naive_pin(TrajectoryData &data, Metadata &meta);
ComplexD Zq_naive_pout(TrajectoryData &data, Metadata &meta);

#endif

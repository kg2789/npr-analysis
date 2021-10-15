#include "Zq.h"

static ComplexD full_trace(SpinColourMatrix &mat) {
    return TensorRemove(trace(mat));
}

static ComplexD Zq_over_ZV_gamma_mu(TrajectoryData &data) {
    ComplexD ret = 0.0;
    for (int mu = 0; mu < Nd; mu++) {
        Gamma gmu = Gamma::gmu[mu];
        SpinColourMatrix &Lambda_mu = data.current.vector[mu];
        SpinColourMatrix tmp = Lambda_mu * gmu;
        ret += full_trace(tmp);
    }
    return 1.0 / 48.0 * ret;
}

static ComplexD Zq_over_ZV_qslash(TrajectoryData &data, std::vector<double> &q) {
    double qsq = vector_magnitude(q);

    SpinColourMatrix Lambda = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        SpinColourMatrix &Lambda_mu = data.current.vector[mu];
        Lambda = Lambda + Lambda_mu * q[mu];
    }
    Lambda = Lambda * vector_slash(q);
    return 1.0 / (12.0 * qsq) * full_trace(Lambda);
}

ComplexD Zq_from_ZV(TrajectoryData &data, Metadata &meta, double ZV,
        WavefunctionScheme scheme) {
    switch (scheme) {
        case WavefunctionScheme::GAMMA_MU:
            return ZV * Zq_over_ZV_gamma_mu(data);
        case WavefunctionScheme::QSLASH:
            std::vector<double> q = get_lattice_q(meta);
            return ZV * Zq_over_ZV_qslash(data, q);
    }
    return 0;
}

static ComplexD Zq_over_ZA_gamma_mu(TrajectoryData &data) {
    ComplexD ret = 0.0;
    for (int mu = 0; mu < Nd; mu++) {
        Gamma gmu = Gamma::gmu[mu];
        Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
        SpinColourMatrix &Lambda_mu = data.current.axial[mu];
        SpinColourMatrix tmp = Lambda_mu * g5 * gmu;
        ret += full_trace(tmp);
    }
    return 1.0 / 48.0 * ret;
}

static ComplexD Zq_over_ZA_qslash(TrajectoryData &data, std::vector<double> &q) {
    double qsq = vector_magnitude(q);

    SpinColourMatrix Lambda = Zero();
    for (int mu = 0; mu < Nd; mu++) {
        SpinColourMatrix &Lambda_mu = data.current.axial[mu];
        Lambda = Lambda + Lambda_mu * q[mu];
    }
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    Lambda = Lambda * g5 * vector_slash(q);
    return 1.0 / (12.0 * qsq) * full_trace(Lambda);
}

ComplexD Zq_from_ZA(TrajectoryData &data, Metadata &meta, double ZA,
        WavefunctionScheme scheme) {
    switch (scheme) {
        case WavefunctionScheme::GAMMA_MU:
            return ZA * Zq_over_ZA_gamma_mu(data);
        case WavefunctionScheme::QSLASH:
            std::vector<double> q = get_lattice_q(meta);
            return ZA * Zq_over_ZA_qslash(data, q);
    }
    return 0;
}

ComplexD Zq_naive_pin(TrajectoryData &data, Metadata &meta) {
    SpinColourMatrix &Sin = data.external_legs.Sin_average;
    std::vector<RealD> pin = momentum_to_lattice_units(meta.p2, meta);
    SpinMatrix pslash = vector_slash(pin);
    RealD psq = vector_magnitude(pin);

    SpinColourMatrix tmp = invert_sc(Sin) * pslash;
    ComplexD amplitude = 1.0 / (12 * psq) * full_trace(tmp);
    ComplexD imag(0, 1);
    return imag * amplitude;
}

ComplexD Zq_naive_pout(TrajectoryData &data, Metadata &meta) {
    SpinColourMatrix &Sout = data.external_legs.Sout_average;
    std::vector<RealD> pout = momentum_to_lattice_units(meta.p1, meta);
    SpinMatrix pslash = vector_slash(pout);
    RealD psq = vector_magnitude(pout);

    SpinColourMatrix tmp = invert_sc(Sout) * pslash;
    ComplexD amplitude = 1.0 / (12 * psq) * full_trace(tmp);
    ComplexD imag(0, 1);
    return imag * amplitude;
}

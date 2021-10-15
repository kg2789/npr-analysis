#include "msbar.h"

/* This file mainly contains the constants needed to convert RI/SMOM into
 * MSbar. The results are taken from the 1-loop calculation given by Lehner and
 * Sturm here: https://arxiv.org/pdf/1104.4948v1.pdf */

// Taken from https://arxiv.org/pdf/1104.4948v1.pdf, Table IX
const double gamma_mu_gamma_mu_r_ij[7][7] =
{
    { 0.21184, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    { 0.0, -0.10592, 0.31777, 0.0, 0.0, 0.0, 0.0},
    { 0.0, 0.09554, -0.62444, 0.07407, -0.22222, 0.0, 0.0},
    { 0.0, 0.0, 0.0, 0.04319, -1.66667, 0.0, 0.0},
    { 0.0, 0.0, -3.88889, -1.05815, 2.82894, 0.0, 0.0},
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.04319, -0.12957},
    { 0.0, 0.0, 0.0, 0.0, 0.0, -1.61371, 4.49561},
};

// Taken from https://arxiv.org/pdf/1104.4948v1.pdf, Table X
const double qslash_gamma_mu_r_ij[7][7] =
{
    { 2.21184, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    { 0.0, 5.29408, 4.91777, 0.0, 0.0, 0.0, 0.0},
    { 0.0, -4.50446, -6.02444, 0.07407, -0.22222, 0.0, 0.0},
    { 0.0, 0.0, 0.0, 2.70986, -0.12957, 0.0, 0.0},
    { 0.0, -1.66667, -3.88889, -0.05815, 2.49561, 0.0, 0.0},
    { 0.0, 0.0, 0.0, 0.0, 0.0, 2.70986, -0.12957},
    { 0.0, 0.0, 0.0, 0.0, 0.0, -0.61371, 4.16227},
};

template <int N>
Eigen::MatrixXcd array_to_eigen(const double values[N][N]) {
    Eigen::MatrixXcd ret(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ret(i, j) = values[i][j];
        }
    }
    return ret;
}

// Computes alpha_s at a given energy scale from alpha_s at the mass of the Z
// boson and the 1-loop running of alpha_s.
// All units of energy here are in GeV
static double get_alpha_s(double scale) {
    const double M_z = 91.1876;
    const double alpha_s_at_Mz = 0.1184;
    const int n_f = 3;
    const double beta_0 = 11 - 2.0 / 3.0 * n_f;
    const double Lambda = M_z * exp(- 2 * M_PI / (beta_0 * alpha_s_at_Mz));

    return 2 * M_PI / (beta_0 * log(scale / Lambda));
}

Eigen::MatrixXcd RI_SMOM_to_MSbar_Zfactor(ProjectorBasis basis,
        Metadata &metadata) {
    Eigen::MatrixXcd ret;
    switch (basis) {
        case ProjectorBasis::GREG:
            ret = array_to_eigen(gamma_mu_gamma_mu_r_ij);
            break;
        case ProjectorBasis::GREG_QSLASH:
            ret = array_to_eigen(qslash_gamma_mu_r_ij);
            break;
        case ProjectorBasis::MASAAKI:
            std::cerr << "Projector basis MASAAKI not supported "
                << "for conversion to MSbar" << std::endl;
            exit(1);
        case ProjectorBasis::MASAAKI_QSLASH:
        	std::cerr << "Projector basis MASAAKI_QSLASH not supported"
        		<< "for conversion to MSbar" << std::endl;
        	exit(1);
    }
    assert(ret.rows() == ret.cols());

    double scale = get_scale(metadata);
    std::cout << "scale = " << scale << std::endl;
    double alpha_s = get_alpha_s(scale);
    std::cout << "alpha_s = " << alpha_s << std::endl;
    ret = ret * alpha_s / (4 * M_PI);

    // Z_{ij} = \delta_{ij} + r_{ij}
    for (int i = 0; i < ret.rows(); i++) {
        ret(i,i) += 1.0;
    }

    return ret;
}

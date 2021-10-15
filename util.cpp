#include "util.h"

using Grid::Nc;
using Grid::Ns;
using Grid::Nd;
using Grid::vTComplex;
using Grid::ComplexD;

template <typename T>
static void swap(T &a, T &b) {
    T tmp = a;
    a = b;
    b = tmp;
}

// Iterate over all entries of a SpinColourSpinColourMatrix
#define FOR_SCSC \
    for (int is = 0; is < Ns; is++) \
    for (int js = 0; js < Ns; js++) \
    for (int ks = 0; ks < Ns; ks++) \
    for (int ls = 0; ls < Ns; ls++) \
    for (int ic = 0; ic < Nc; ic++) \
    for (int jc = 0; jc < Nc; jc++) \
    for (int kc = 0; kc < Nc; kc++) \
    for (int lc = 0; lc < Nc; lc++) \

static void swap02(SpinColourSpinColourMatrix& mat) {
    SpinColourSpinColourMatrix ret;
    FOR_SCSC {
        // Swap i <-> k
        ret()(is,js)(ic,jc)(ks,ls)(kc,lc) = mat()(ks,js)(kc,jc)(is,ls)(ic,lc);
    }
    mat = ret;
}

static void swap13(SpinColourSpinColourMatrix& mat) {
    SpinColourSpinColourMatrix ret;
    FOR_SCSC {
        // Swap j <-> l
        ret()(is,js)(ic,jc)(ks,ls)(kc,lc) = mat()(is,ls)(ic,lc)(ks,js)(kc,jc);
    }
    mat = ret;
}


void scsc_index_swap(SpinColourSpinColourMatrix& mat, int a, int b) {
    // Ensure a < b
    if (b < a) swap(a,b);
    if (a == 0 && b == 2) {
        swap02(mat);
    }
    else if (a == 1 && b == 3) {
        swap13(mat);
    }
    else {
        std::cerr << "Error: unimplemented swap values "
            << a << ", " << b << std::endl;
        assert(0);
    }
}
SpinColourMatrix spin_colour_identity() {
    SpinColourMatrix ret = ComplexD(1.0);
    return ret;
}
SpinMatrix spin_identity() {
    SpinMatrix ret = ComplexD(1.0);
    return ret;
}


SpinMatrix gamma_matrix(int mu) {
    using Grid::Gamma;
    Gamma gamma_mu = Gamma::gmu[mu];
    return gamma_mu * spin_identity();
}

SpinMatrix gamma5() {
    using Grid::Gamma;
    Gamma g5(Gamma::Algebra::Gamma5);
    return g5 * spin_identity();
}

SpinMatrix vector_slash(std::vector<RealD> &v) {
    SpinMatrix ret;
    assert(v.size() == 4);
    for (int mu = 0; mu < 4; mu++) {
        ret += v[mu] * gamma_matrix(mu);
    }
    return ret;
}

SpinColourMatrix vector_slash_sc(std::vector<RealD> &v) {
    return spin_colour_identity() * vector_slash(v);
}

const int Nsc = Ns * Nc;
typedef Eigen::Matrix<Complex, Nsc, Nsc> EigenSpinColourMatrix;

static SpinColourMatrix eigen_to_grid(const EigenSpinColourMatrix &eigen_mat) {
    SpinColourMatrix mat;
    for (int is = 0; is < Ns; is++) {
    for (int js = 0; js < Ns; js++) {
        for (int ic = 0; ic < Nc; ic++) {
        for (int jc = 0; jc < Nc; jc++) {
            mat()(is,js)(ic,jc)
                = eigen_mat(Ns * ic + is, Ns * jc + js);
        }
        }
    }
    }
    return mat;
}

static EigenSpinColourMatrix grid_to_eigen(const SpinColourMatrix &mat) {
    EigenSpinColourMatrix eigen_mat;
    for (int is = 0; is < Ns; is++) {
    for (int js = 0; js < Ns; js++) {
        for (int ic = 0; ic < Nc; ic++) {
        for (int jc = 0; jc < Nc; jc++) {
            eigen_mat(Ns * ic + is, Ns * jc + js)
                = mat()(is,js)(ic,jc);
        }
        }
    }
    }
    return eigen_mat;
}

SpinColourMatrix invert_sc(const SpinColourMatrix &mat) {
    EigenSpinColourMatrix eigen_mat = grid_to_eigen(mat);
    return eigen_to_grid(eigen_mat.inverse());
}

std::vector<double> momentum_to_lattice_units(const std::vector<double> &n_mu,
        const Metadata &meta) {
    assert(n_mu.size() == meta.lattice_dimensions.size());
    std::vector<double> ret;
    ret.resize(n_mu.size());
    for (int mu = 0; mu < n_mu.size(); mu++) {
        const int Lmu = meta.lattice_dimensions[mu];
        ret[mu] = 2 * M_PI / Lmu * n_mu[mu];
    }
    return ret;
}

double vector_magnitude(const std::vector<double> &v) {
    double ret = 0.0;
    for (auto &val: v)
        ret += val * val;
    return ret;
}

std::vector<double> get_lattice_q(const Metadata &meta) {
    std::vector<double> p1 = momentum_to_lattice_units(meta.p1, meta);
    std::vector<double> p2 = momentum_to_lattice_units(meta.p2, meta);
    assert(p1.size() == p2.size());

    std::vector<double> q = p1;
    for (int mu = 0; mu < p2.size(); mu++) {
        q[mu] -= p2[mu];
    }
    return q;
}

static SpinColourMatrix spin_structure_matrix(const SpinStructure &spin, int mu) {
    using Grid::Gamma;
    Gamma g5(Gamma::Algebra::Gamma5);

    SpinColourMatrix identity = spin_colour_identity();
    Gamma gamma_mu = Gamma::gmu[mu];
    switch (spin) {
        case SpinStructure::LEFT:
            return gamma_mu * (identity - g5 * identity);
            break;
        case SpinStructure::RIGHT:
            return gamma_mu * (identity + g5 * identity);
            break;
        case SpinStructure::VECTOR:
            return gamma_mu * identity;
            break;
        case SpinStructure::AXIAL:
            return gamma_mu * g5 * identity;
        default:
            return identity;
    }
}

void tensor_prod(SpinColourSpinColourMatrix &lret,
        const SpinColourMatrix &a,
        const SpinColourMatrix &b) {
    // Copied from Julia Kettle's code on bilinears with some modification:
    // https://github.com/Julia-Kettle/Grid/blob/c410bdd439c8445d50a231e49dc7135c09e7de9e/Hadrons/Modules/MNPR/Bilinear.hpp#L142
    ComplexD left;
    for(int si=0; si < Ns; ++si){
        for(int sj=0; sj < Ns; ++sj){
            for (int ci=0; ci < Nc; ++ci){
                for (int cj=0; cj < Nc; ++cj){
                    left = a()(si,sj)(ci,cj);
                    lret()(si,sj)(ci,cj)=left*b();
               }
            }
        }
    }
}

// As SpinColourMatrices we consider the gamma-matrices to be color diagonal;
// this function lets us change the color structure to be mixed
// (\delta^{ad} delta^{bc})
static void impose_color_structure(SpinColourSpinColourMatrix &lret,
        ColorStructure structure) {
    if (structure == ColorStructure::COLOR_DIAG) {
        return;
    }
    if (structure == ColorStructure::UNCONTRACTED) {
        std::cerr << "Error: unsupported structure UNCONTRACTED "
            << "in impose_color_structure" << std::endl;
        assert(0);
        return;
    }
    assert(structure == ColorStructure::COLOR_MIXED);
    // This could probably be more efficient, but this should work
    for (int is = 0; is < Ns; is++) {
    for (int js = 0; js < Ns; js++) {
    for (int ks = 0; ks < Ns; ks++) {
    for (int ls = 0; ls < Ns; ls++) {
        ComplexD value = lret()(is,js)(0,0)(ks,ls)(0,0);
        for (int ic = 0; ic < Nc; ic++) {
        for (int jc = 0; jc < Nc; jc++) {
            lret()(is,js)(ic,ic)(ks,ls)(jc,jc) = 0;
            lret()(is,js)(ic,jc)(ks,ls)(jc,ic) = value;
        }}
    }}}}
}

// Computes \Gamma_a \otimes \Gamma_b with the specified color structure, where
// \Gamma_a and \Gamma_b are the gamma-matrix structures specified by a and b
SpinColourSpinColourMatrix
vertex_structure(const SpinStructure &a, const SpinStructure &b,
        const ColorStructure& structure) {
    SpinColourSpinColourMatrix ret;
    for (int mu = 0; mu < Nd; mu++) {
        SpinColourMatrix mat_a = spin_structure_matrix(a, mu);
        SpinColourMatrix mat_b = spin_structure_matrix(b, mu);
        SpinColourSpinColourMatrix tmp;

        tensor_prod(tmp, mat_a, mat_b);
        impose_color_structure(tmp, structure);
        ret += tmp;
    }
    return ret;
}

SpinColourSpinColourMatrix
tensor_prod_with_color_structure(const SpinMatrix &a, const SpinMatrix &b,
        ColorStructure color) {
    SpinColourSpinColourMatrix ret = Grid::Zero();
    for (int is = 0; is < Ns; is++) {
    for (int js = 0; js < Ns; js++) {
    for (int ks = 0; ks < Ns; ks++) {
    for (int ls = 0; ls < Ns; ls++) {
        ComplexD value = a()(is,js)() * b()(ks,ls)();
        for (int ic = 0; ic < Nc; ic++) {
        for (int jc = 0; jc < Nc; jc++) {
            switch (color) {
                case ColorStructure::COLOR_DIAG:
                    ret()(is,js)(ic,ic)(ks,ls)(jc,jc) = value;
                    break;
                case ColorStructure::COLOR_MIXED:
                    ret()(is,js)(ic,jc)(ks,ls)(jc,ic) = value;
                    break;
                default:
                    std::cerr << "Error: color structure " << color
                        << " not supported for tensor_prod_with_color_structure"
                        << std::endl;
                    exit(1);
            }
        }}
    }}}}

    return ret;
}

SpinColourSpinColourMatrix
spin_mixed_tensor_prod_with_color_structure(const SpinMatrix &a, const SpinMatrix &b,
        ColorStructure color) {
    SpinColourSpinColourMatrix ret = Grid::Zero();
    for (int is = 0; is < Ns; is++) {
    for (int js = 0; js < Ns; js++) {
    for (int ks = 0; ks < Ns; ks++) {
    for (int ls = 0; ls < Ns; ls++) {
        ComplexD value = a()(is,js)() * b()(ks,ls)();
        for (int ic = 0; ic < Nc; ic++) {
        for (int jc = 0; jc < Nc; jc++) {
            switch (color) {
                case ColorStructure::COLOR_DIAG:
                    ret()(is,ls)(ic,ic)(ks,js)(jc,jc) = value;
                    break;
                case ColorStructure::COLOR_MIXED:
                    ret()(is,ls)(ic,jc)(ks,js)(jc,ic) = value;
                    break;
                default:
                    std::cerr << "Error: color structure " << color
                        << " not supported for tensor_prod_with_color_structure"
                        << std::endl;
                    exit(1);
            }
        }}
    }}}}

    return ret;
}

void print_scsc(const SpinColourSpinColourMatrix &mat) {
    int i = 0;
    const int width = Ns * Nd * Ns * Nd;
    FOR_SCSC {
        std::cout << mat()(is,js)(ic,jc)(ks,ls)(kc,lc) << " ";
        if (i++ % width == 0 && i != 0) {
            std::cout << std::endl;
        }
    }
}

void print_scsc_real_part(const SpinColourSpinColourMatrix &mat) {
    int i = 0;
    const int width = Ns * Nd * Ns * Nd;
    FOR_SCSC {
        std::cout << mat()(is,js)(ic,jc)(ks,ls)(kc,lc).real() << " ";
        if (i % width == 0 && i != 0) {
            std::cout << std::endl;
        }
        i++;
    }
}

void print_scsc_nonzero(const SpinColourSpinColourMatrix &mat) {
    std::cout << "s, s, c, c, s, s, c, c" << std::endl;
    FOR_SCSC {
        RealD value = mat()(is,js)(ic,jc)(ks,ls)(kc,lc).real();
        if (value != 0) {
            std::cout << is
                << ", " << js
                << ", " << ic
                << ", " << jc
                << ", " << ks
                << ", " << ls
                << ", " << kc
                << ", " << lc
                << ": " << value << std::endl;
        }
    }
}

void print_sc(const SpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int c1 = 0; c1 < Nc; c1++) {
        for (int s2 = 0; s2 < Ns; s2++) {
        for (int c2 = 0; c2 < Nc; c2++) {
            std::cout << mat()(s1,s2)(c1,c2) << " ";
        }
        }
        std::cout << std::endl;
    }
    }
}

void print_sc_real(const SpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int c1 = 0; c1 < Nc; c1++) {
        for (int s2 = 0; s2 < Ns; s2++) {
        for (int c2 = 0; c2 < Nc; c2++) {
            std::cout << mat()(s1,s2)(c1,c2).real() << " ";
        }
        }
        std::cout << std::endl;
    }
    }
}

// Assumes the matrix is color diagonal, so we can ignore the color indices
void print_sc_spin(const SpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
        for (int s2 = 0; s2 < Ns; s2++) {
            std::cout << mat()(s1,s2)(0,0) << " ";
        }
        std::cout << std::endl;
    }
}

// Assumes the matrix is color diagonal, so we can ignore the color indices
void print_sc_spin_real(const SpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
        for (int s2 = 0; s2 < Ns; s2++) {
            std::cout << mat()(s1,s2)(0,0).real() << " ";
        }
        std::cout << std::endl;
    }
}

// Assumes the matrix is color diagonal, so we can ignore the color indices
void print_scsc_spin(const SpinColourSpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
        for (int s3 = 0; s3 < Ns; s3++) {
    for (int s2 = 0; s2 < Ns; s2++) {
        for (int s4 = 0; s4 < Ns; s4++) {
            std::cout << mat()(s1,s2)(0,0)(s3,s4)(0,0) << " ";
        }
    }
    std::cout << std::endl;
        }
    }
}

// Assumes the matrix is color diagonal, so we can ignore the color indices
void print_scsc_spin_real(const SpinColourSpinColourMatrix &mat) {
    for (int s1 = 0; s1 < Ns; s1++) {
        for (int s3 = 0; s3 < Ns; s3++) {
    for (int s2 = 0; s2 < Ns; s2++) {
        for (int s4 = 0; s4 < Ns; s4++) {
            std::cout << mat()(s1,s2)(0,0)(s3,s4)(0,0).real() << " ";
        }
        }
        std::cout << std::endl;
    }
    }
}


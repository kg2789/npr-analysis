#include "parity.h"
#include "util.h"
#include "data.h"

using Grid::Gamma;

template <typename T>
static void parity_transform_matrix(T &mat) {
    Gamma gt = Gamma(Gamma::Algebra::GammaT);
    mat = gt * mat * gt;
}

void parity_transform(SpinColourMatrix &mat) {
    parity_transform_matrix(mat);
}

static SpinColourSpinColourMatrix 
make_tensor_prod(const SpinColourMatrix &a,
        const SpinColourMatrix &b) {
    SpinColourSpinColourMatrix ret;
    tensor_prod(ret, a, b);
    return ret;
}

void parity_transform(SpinColourSpinColourMatrix &mat) {
    using Grid::Ns;
    using Grid::Nc;

    SpinColourMatrix gt = spin_colour_identity() * gamma_matrix(Grid::Tp);
    static const SpinColourSpinColourMatrix gtgt = make_tensor_prod(gt,gt);

    mat = gtgt * mat * gtgt;
    /*for (int s1 = 0; s1 < Ns; s1++) {
    for (int s2 = 0; s2 < Ns; s2++) {
        for (int c1 = 0; c1 < Nc; c1++) {
        for (int c2 = 0; c2 < Nc; c2++) {
            parity_tranform_matrix(mat()(s1,s2)(c1,c2));
        }
        }
    }
    }
    parity_tranform_matrix(mat);*/
}

SpinColourMatrix parity_project(const SpinColourMatrix &mat, Parity parity) {
    SpinColourMatrix tmp = mat;
    parity_transform(tmp);
    switch (parity) {
        case Parity::EVEN:
            return mat + tmp;
        case Parity::ODD:
            return mat - tmp;
        case Parity::BOTH:
            return mat;
    }
}

SpinColourSpinColourMatrix
parity_project(const SpinColourSpinColourMatrix &mat, Parity parity) {
    SpinColourSpinColourMatrix tmp = mat;
    parity_transform(tmp);
    switch (parity) {
        case Parity::EVEN:
            return (1.0 / 2) * (mat + tmp);
        case Parity::ODD:
            return (1.0 / 2) * (mat - tmp);
        case Parity::BOTH:
            return mat;
    }
}

template <typename T>
void fourq_op_parity_project(T &measurements, Parity parity) {
    for (LoopStructure loop: all_loop_structures) {
        for (ColorStructure color: all_color_structures) {
            for (SpinStructure gammaA: va_spin_structures) {
                for (SpinStructure gammaB: va_spin_structures) {
                    switch (parity) {
                        case Parity::EVEN:
                            if (gammaA == gammaB) continue;
                            break;
                        case Parity::ODD:
                            if (gammaA != gammaB) continue;
                            break;
                        case Parity::BOTH:
                            break;
                    }
                    measurements.get_measurement(loop, color, gammaA, gammaB) = Zero();
                }
            }
        }
    }
}

template <typename T>
void twoq_op_parity_project(T &measurements, Parity parity) {
    for (TwoQuarkOp op: all_subtraction_ops) {
        int idx = (int) op;
        switch (parity) {
            case Parity::BOTH:
                break;
            case Parity::EVEN:
                measurements.psuedoscalar_values[idx] = Zero();
                break;
            case Parity::ODD:
                measurements.values[idx] = Zero();
        }
    }
}

void project_parity(TrajectoryData &data, Parity parity) {
    fourq_op_parity_project(data.fourq_op.fourq_ext, parity);
    fourq_op_parity_project(data.fourq_op.twoq_ext, parity);
    twoq_op_parity_project(data.twoq_op.fourq_ext, parity);
    twoq_op_parity_project(data.twoq_op.twoq_ext, parity);
}

RunData read_data_with_parity(const std::string &dir, Parity parity) {
    RunData ret = read_data(dir, Format::RABBOTT, false);
    for (TrajectoryData &data: ret.trajectory_data) {
        project_parity(data, parity);
        data.switch_modes(FourQDataMode::LEFT_RIGHT);
    }
    return ret;
}

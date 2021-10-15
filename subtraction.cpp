#include "subtraction.h"
#include "projector.h"

using Eigen::MatrixXcd;
using Eigen::VectorXcd;

static VectorXcd compute_coefficients(SpinColourMatrix &fourq_diagram,
        std::vector<SpinColourMatrix> &subtraction_diagrams,
        std::vector<TwoQuarkProjector> &projectors) {
    // We want to compute the coefficients x_j such that
    // projector[i].project(fourq_diagram - \sum_j x_j subtraction_diagrams[j]) = 0
    // rearranging we can put this in the form Ax = b where
    //
    // A_ij = projector[i].project(subtraction_diagrams[j])
    // b_i = projector[i].project(fourq_diagram)
    //
    // Note that in order for the solution to be unique A must be an invertible
    // (and hence square matrix). This reflects the requirement that we must
    // have the same number of projectors as sutraction ops.
    assert(projectors.size() == subtraction_diagrams.size());
    MatrixXcd A(projectors.size(), subtraction_diagrams.size());
    VectorXcd b(projectors.size());
    for (int i = 0; i < projectors.size(); i++) {
        for (int j = 0; j < subtraction_diagrams.size(); j++) {
            A(i,j) = projectors[i].project(subtraction_diagrams[j]);
        }
        b(i) = projectors[i].project(fourq_diagram);
    }
    // Solving the linear equation Ax = b; taken from here:
    // https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    VectorXcd x = A.colPivHouseholderQr().solve(b);
    return x;
}

MatrixXcd subtraction_coefficients(TrajectoryData *data,
        std::vector<FourQuarkOperator> &operators,
        std::vector<TwoQuarkOp> &subtraction_ops,
        std::vector<TwoQuarkProjector> &projectors) {
    MatrixXcd ret = MatrixXcd(operators.size(), subtraction_ops.size());

    BilinearOperator external_state
        = BilinearTerm(Quark::S, Quark::DBAR, SpinStructure::UNCONTRACTED);
    BilinearOperator subtraction_state
        = BilinearTerm(Quark::SBAR, Quark::D, SpinStructure::UNCONTRACTED);

    for (int i = 0; i < operators.size(); i++) {
        FourQuarkOperator &op = operators[i];
        // Compute 4-quark diagram
        MatrixElement<BilinearTerm, FourQuarkTerm> matrix_element
            = make_matrix_element(external_state, op);
        SpinColourMatrix fourq_diagram = evaluate_matrix_element(*data, matrix_element);

        // Compute 2-quark diagrams
        std::vector<SpinColourMatrix> twoq_diagrams;
        twoq_diagrams.reserve(subtraction_ops.size());
        // As evaluate_matrix_element() is currently implemented for diagrams
        // with bilinear opeartors in 2-quark external states this will mostly be ignored
        MatrixElement<BilinearTerm, BilinearTerm> twoq_matrix_elem
            = make_matrix_element(external_state, subtraction_state);
        for (TwoQuarkOp &sub_op: subtraction_ops) {
            twoq_diagrams.push_back(evaluate_matrix_element(*data, twoq_matrix_elem, sub_op));
        }

        // Computing the subtraction coefficients for the given operator
        VectorXcd coefficients = compute_coefficients(fourq_diagram, twoq_diagrams,
                                                      projectors);
        for (int j = 0; j < subtraction_ops.size(); j++) {
            ret(i, j) = coefficients(j);
        }
    }
    return ret;
}

// By definition Q_j^sub = Q_j - \sum_k coeff[j][k] subtaction_op[k], where
// coeff[j][k] are the subtraction coefficients determined above. We want to
// compute:
//
// projector[i].project( <Q_j^sub> )
// = projector[i].project( <Q_j - \sum_k coeff[j][k] subtraction_op[k] > )
// = projector[i].project(<Q_j>)
//          - \sum_k coeff[j][k] projector[i].project(<subtraction_op[k]>)
//
// (the projectors here are 4-quark projectors)
//
// This function computes the second term, which gives the contribution of the
// subtraction operators to the mixing matrix.
MatrixXcd subtraction_correction(TrajectoryData *data,
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &fourq_projectors,
        std::vector<OperatorRepresentation> &fourq_op_reps,
        std::vector<TwoQuarkOp> &subtraction_ops,
        std::vector<TwoQuarkProjector> &twoq_projectors,
        std::vector<OperatorRepresentation> &subtraction_op_reps) {
    MatrixXcd coeff = subtraction_coefficients(data,
            operators, subtraction_ops, twoq_projectors);

    // Enforcing non-mixing between different representations
    assert(fourq_op_reps.size() == coeff.rows());
    assert(subtraction_op_reps.size() == coeff.cols());
    for (int i = 0; i < fourq_op_reps.size(); i++) {
        for (int j = 0; j < subtraction_op_reps.size(); j++) {
            if (fourq_op_reps[i] != subtraction_op_reps[j]) {
                coeff(i, j) = 0.0;
            }
        }
    }

    // Matrix for computing projector[i].project(<sub_op[k]>)
    MatrixXcd projected_subtractions(fourq_projectors.size(),
                                     subtraction_ops.size());

    BilinearOperator twoq_structure = BilinearTerm(Quark::SBAR, Quark::D,
                                                   SpinStructure::UNCONTRACTED);
    for (int i = 0; i < fourq_projectors.size(); i++) {
        FourQuarkOperator &ext = external_states[i];
        FourQuarkProjector &proj = fourq_projectors[i];
        for (int k = 0; k < subtraction_ops.size(); k++) {
            // Enforcing non-mixing between representations
            if (fourq_op_reps[i] != subtraction_op_reps[k]) {
                projected_subtractions(i,k) = 0.0;
                continue;
            }

            MatrixElement<FourQuarkTerm, BilinearTerm> matrix_elem
                = make_matrix_element(ext, twoq_structure);
            SpinColourSpinColourMatrix amplitude
                = evaluate_matrix_element(*data, matrix_elem, subtraction_ops[k]);

            projected_subtractions(i,k) = proj.project(amplitude);
        }
    }

    // ret[j][i] = \sum_k coeff[j][k] projector[i].project(<subtraction_op[k]>)
    return coeff * projected_subtractions.transpose();
}


std::vector<OperatorRepresentation> get_subtraction_reps(ProjectorBasis basis) {
    switch (basis) {
        case ProjectorBasis::GREG:
            return std::vector<OperatorRepresentation>(3,
                    OperatorRepresentation::NONE);
        case ProjectorBasis::GREG_QSLASH:
            return std::vector<OperatorRepresentation>(3,
                    OperatorRepresentation::REP8_1);
        case ProjectorBasis::MASAAKI:
            return std::vector<OperatorRepresentation>(3,
                    OperatorRepresentation::NONE);
        case ProjectorBasis::MASAAKI_QSLASH:
            return std::vector<OperatorRepresentation>(3,
                    OperatorRepresentation::NONE);           
    }
}

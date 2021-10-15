#include "../evaluate.h"
#include "../data.h"

SpinColourMatrix evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<BilinearTerm, BilinearTerm> &matrix_element,
        TwoQuarkOp op) {
    SpinColourMatrix ret;
    TwoQOpTwoQExtMeasurements& measurements = data.twoq_op.twoq_ext;
    for (auto &p: matrix_element.term_coefficients) {
        RealD& coeff = p.second;

        // There's only 1 possible contraction with a 2-quark op in a 2-quark
        // external state, so we don't actually need to use the diagram
        //
        // The minus sign is a consequence of Greg's sign choices for the
        // 4-quark operators. In particular his code counts loops including
        // external lines, which means we have 1 loop here and hence have an
        // overall minus sign. Greg's code does not have this minus sign; this
        // problem was discovered when Greg's and Daiaqian's data disagreed.
        ret += -coeff * measurements.get_measurement(op);
    }
    return ret;
}

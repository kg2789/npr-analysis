#include "evaluate.h"
#include "projector.h"
#include "data.h"

Eigen::MatrixXcd compute_mixings(TrajectoryData &data,
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &projection_ops,
        std::vector<OperatorRepresentation> &representations) {
    using Eigen::MatrixXcd;

    assert(operators.size() == external_states.size()
            && external_states.size() == projection_ops.size());
    assert(operators.size() == representations.size());

    const int num_ops = operators.size();
    MatrixXcd ret(num_ops, num_ops);
    for (int i = 0; i < num_ops; i++) {
        FourQuarkOperator &op = operators[i];
        for (int j = 0; j < num_ops; j++) {
            // See note in projector.cpp under get_representations()
            if (representations[i] != representations[j]) {
                ret(i,j) = 0.0;
                continue;
            }

            FourQuarkOperator &ext = external_states[j];
            FourQuarkProjector &proj = projection_ops[j];

            MatrixElement<FourQuarkTerm, FourQuarkTerm> matrix_element
                = make_matrix_element(ext, op);
            SpinColourSpinColourMatrix amplitude
                = evaluate_matrix_element(data, matrix_element);
            ret(i,j) = proj.project(amplitude);
        }
    }
    return ret;
}

#include "tree.h"
#include "evaluate.h"
#include "projector.h"
#include "parity.h"
#include "data.h"


static SpinColourSpinColourMatrix
get_tree_level_measurment(LoopStructure loop, ColorStructure color,
        SpinStructure gammaA, SpinStructure gammaB) {
    SpinColourSpinColourMatrix ret;

    // Only fully-connected diagrams are tree-level
    if (loop != LoopStructure::FULLY_CONNECTED) {
        return ret;
    }

    return vertex_structure(gammaA, gammaB, color);
}

TrajectoryData *tree_level_measurements(Parity parity) {
    FourQOpFourQExtMeasurements *fourq = new FourQOpFourQExtMeasurements();

    std::vector<SpinStructure> spin_structures = {
        SpinStructure::VECTOR, SpinStructure::AXIAL,
        //SpinStructure::LEFT, SpinStructure::RIGHT
    };

    for (int a = 0; a < NUM_LOOP_STRUCTURES; a++) {
        for (int b = 0; b < NUM_COLOR_STRUCTURES; b++) {
            for (SpinStructure &gammaA: spin_structures) {
                for (SpinStructure &gammaB: spin_structures) {
                    LoopStructure loop = (LoopStructure) a;
                    ColorStructure color = (ColorStructure) b;

                    if (parity == Parity::EVEN && gammaA != gammaB)
                        continue;
                    if (parity == Parity::ODD && gammaA == gammaB)
                        continue;

                    SpinColourSpinColourMatrix tmp
                        = get_tree_level_measurment(loop, color, gammaA, gammaB);
                    SpinColourSpinColourMatrix &out
                        = fourq->get_measurement(loop, color, gammaA, gammaB);
                    out  = tmp;
                }
            }
        }
    }

    TrajectoryData *ret = new TrajectoryData();
    ret->fourq_op.fourq_ext = *fourq;
    ret->switch_modes(FourQDataMode::LEFT_RIGHT);

    return ret;
}

Eigen::MatrixXcd tree_level_mixings(
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &projection_ops,
        std::vector<OperatorRepresentation> &reps,
        Parity parity) {
    using Eigen::MatrixXcd;

    assert(operators.size() == external_states.size()
            && external_states.size() == projection_ops.size());

    TrajectoryData *data = tree_level_measurements(parity);
    return compute_mixings(*data, operators, external_states, projection_ops, reps);
}


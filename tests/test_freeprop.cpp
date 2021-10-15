#include <string>
#include "../data.h"
#include "../amputate.h"
#include "../tree.h"
#include "../data/debug.h"

// If we set every gauge link to the identity matrix, then the plane wave
// source exp(i p \cdot x) we use is an eigenfunction of the Dirac operator.
// This implies that the propagator we compute should have the form A exp(i p
// \cdot x), where A is some (constant) matrix. Averaging over the lattice we
// should get this matrix A, so if we amputate then the external propagators
// become the identity. By definition a tree-level result has all external
// propagators equal to the identity, so the result for a unit gauge should be
// the same as tree-level for the nonzero tree-level diagrams. This provides a
// check of both the measurement code and the implementation of the Fierz
// identity to produce the mixed-color-structure diagrams


const RealD tolerance = 0.001;


// Returns success
bool test_trajectory(TrajectoryData &data) {
    // Zero out any non tree-level values
    for (int a = 0; a < NUM_LOOP_STRUCTURES; a++) {
        if ((LoopStructure) a == LoopStructure::FULLY_CONNECTED)
            continue;
        for (int b = 0; b < NUM_COLOR_STRUCTURES; b++) {
            for (int c = 0; c < NUM_SPIN_STRUCTURES; c++) {
                for (int d = 0; d < NUM_SPIN_STRUCTURES; d++) {
                    LoopStructure loop = (LoopStructure) a;
                    ColorStructure color = (ColorStructure) b;
                    SpinStructure gammaA = (SpinStructure) c;
                    SpinStructure gammaB = (SpinStructure) d;

                    if (!data.fourq_op.fourq_ext.is_measurement_stored(gammaA, gammaB)) {
                        continue;
                    }

                    auto &val = data.fourq_op.fourq_ext
                        .get_measurement(loop, color, gammaA, gammaB);

                    val = Zero();
                }
            }
        }
    }

    // Get external propagators
    SpinColourMatrix Sin = data.external_legs.Sin_average;
    SpinColourMatrix Sout = data.external_legs.Sout_average;

    // Ampuate and compare with tree level
    amputate(data, Sin, Sout);
    TrajectoryData *tree_level = tree_level_measurements(Parity::BOTH);

    RealD diff = relative_difference(tree_level->fourq_op.fourq_ext,
                                     data.fourq_op.fourq_ext);

    std::cout << "Realtive difference for trajectory "
        << data.trajectory << ": " << diff << std::endl;

    if (diff >= tolerance) {
        print_diff_profile(*tree_level, data);
    }

    return diff < tolerance;
}

int main() {
    std::string dir = "../../../test_data/unit_gauge_4x4x4x4";

    RunData data = read_data(dir);

    for (TrajectoryData &traj : data.trajectory_data) {
        if (!test_trajectory(traj))
            return 1;
    }

    return 0;
}

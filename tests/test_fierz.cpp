#include <string>
#include "../data.h"

void print_diff_profile(TrajectoryData *data_a, TrajectoryData *data_b) {
    RealD twoq_op_diff = relative_difference(data_a->twoq_op, data_b->twoq_op);
    RealD fourq_op_diff = relative_difference(data_a->fourq_op, data_b->fourq_op);

    std::cout << "twoq_op: " << twoq_op_diff << std::endl;
    std::cout << "fourq_op: " << fourq_op_diff << std::endl;
    std::cout << "fourq_op.fourq_ext: " << relative_difference(data_a->fourq_op.fourq_ext,
                                                               data_b->fourq_op.fourq_ext) << std::endl;
    std::cout << "fourq_op.twoq_ext: " << relative_difference(data_a->fourq_op.twoq_ext,
                                                              data_b->fourq_op.twoq_ext) << std::endl;
}

const double EPSILON = 0.001;

bool test_trajectory(TrajectoryData *raw, TrajectoryData *raw_with_mixed) {
    RealD diff = relative_difference(*raw, *raw_with_mixed);

    std::cout << "Diff = " << diff << std::endl;

    if (diff > EPSILON) {
        std::cout << "Error: diff = " << diff << " > "
            << EPSILON << std::endl;

        print_diff_profile(raw, raw_with_mixed);
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Error: no data directory given" << std::endl;
        return 1;
    }
    std::string dir(argv[1]);

    RunData data = read_data(dir, Format::RAW);
    RunData data_with_mixed = read_data(dir, Format::RAW_WITH_MIXED);

    for (int i = 0; i < data.trajectory_data.size(); i++) {
        TrajectoryData *raw = &data.trajectory_data[i];
        TrajectoryData *raw_with_mixed = &data_with_mixed.trajectory_data[i];

        std::cout << "Testing trajectory " << raw->trajectory << std::endl;
        if (!test_trajectory(raw, raw_with_mixed))
            return 1;
    }

    return 0;
}

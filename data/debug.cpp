#include "debug.h"

#include <iostream>
#include "../data.h"

void print_diff_profile(const TrajectoryData &traj_A, const TrajectoryData &traj_B) {
    std::cout << "norm A: " << l2_norm(traj_A) << std::endl;
    std::cout << "norm B " << l2_norm(traj_B) << std::endl;

    std::cout << "fourq_op.fourq_ext: " << std::endl;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.fourq_op.fourq_ext) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.fourq_op.fourq_ext) << std::endl;
    for (LoopStructure loop: all_loop_structures) {
        std::cout << "\t" << loop << ": " <<
            l2_diagram_difference(traj_A.fourq_op.fourq_ext,
                    traj_B.fourq_op.fourq_ext, loop) << std::endl;
    }
    std::cout << "fourq_op.twoq_ext: " << std::endl;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.fourq_op.twoq_ext) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.fourq_op.twoq_ext) << std::endl;
    for (LoopStructure loop: all_loop_structures) {
        std::cout << "\t" << loop << ": " <<
            l2_diagram_difference(traj_A.fourq_op.twoq_ext,
                    traj_B.fourq_op.twoq_ext, loop) << std::endl;
    }
    std::cout << "twoq_op.fourq_ext: "
        << l2_difference(traj_A.twoq_op.fourq_ext,
                traj_B.twoq_op.fourq_ext) << std::endl;;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.twoq_op.fourq_ext) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.twoq_op.fourq_ext) << std::endl;
    std::cout << "twoq_op.twoq_ext: "
        << l2_difference(traj_A.twoq_op.twoq_ext,
                traj_B.twoq_op.twoq_ext) << std::endl;;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.twoq_op.twoq_ext) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.twoq_op.twoq_ext) << std::endl;

    std::cout << "Incoming external leg: "
        << l2_difference(traj_A.external_legs.Sin_average,
                traj_B.external_legs.Sin_average) << std::endl;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.external_legs.Sin_average) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.external_legs.Sin_average) << std::endl;

    std::cout << "Outgoing external leg: "
        << l2_difference(traj_A.external_legs.Sout_average,
                traj_B.external_legs.Sout_average) << std::endl;
    std::cout << "\tnorm A: "
        << l2_norm(traj_A.external_legs.Sout_average) << std::endl;
    std::cout << "\tnorm B: "
        << l2_norm(traj_B.external_legs.Sout_average) << std::endl;

    std::cout << "G1: " << l2_difference(traj_A.g1, traj_B.g1) << std::endl;
    std::cout << "\tnorm A: " << l2_norm(traj_A.g1) << std::endl;
    std::cout << "\tnorm B: " << l2_norm(traj_B.g1) << std::endl;
    std::cout << "\ttwoq_scalar: "
        << l2_difference(traj_A.g1.twoq_scalar, traj_B.g1.twoq_scalar)
        << std::endl;
    std::cout << "\t\tnorm A: " << l2_norm(traj_A.g1.twoq_scalar) << std::endl;
    std::cout << "\t\tnorm B: " << l2_norm(traj_B.g1.twoq_scalar) << std::endl;
    std::cout << "\ttwoq_gamma5: "
        << l2_difference(traj_A.g1.twoq_gamma5, traj_B.g1.twoq_gamma5)
        << std::endl;
    std::cout << "\t\tnorm A: " << l2_norm(traj_A.g1.twoq_gamma5) << std::endl;
    std::cout << "\t\tnorm B: " << l2_norm(traj_B.g1.twoq_gamma5) << std::endl;
    std::cout << "\tfourq_scalar: "
        << l2_difference(traj_A.g1.fourq_scalar, traj_B.g1.fourq_scalar)
        << std::endl;
    std::cout << "\t\tnorm A: " << l2_norm(traj_A.g1.fourq_scalar) << std::endl;
    std::cout << "\t\tnorm B: " << l2_norm(traj_B.g1.fourq_scalar) << std::endl;
    std::cout << "\tfourq_gamma5: "
        << l2_difference(traj_A.g1.fourq_gamma5, traj_B.g1.fourq_gamma5)
        << std::endl;
    std::cout << "\t\tnorm A: " << l2_norm(traj_A.g1.fourq_gamma5) << std::endl;
    std::cout << "\t\tnorm B: " << l2_norm(traj_B.g1.fourq_gamma5) << std::endl;
}

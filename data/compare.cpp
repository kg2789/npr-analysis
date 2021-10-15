#include <algorithm>
#include <string>

#include <tclap/CmdLine.h>

#include "format.h"
#include "io.h"
#include "../data.h"
#include "debug.h"

int main(int argc, char **argv) {
    Format input_a_format, input_b_format;
    std::string input_a_dir, input_b_dir;

    // Using tclap to parse command-line arguments; see here for more detail:
    // https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
    try {
        TCLAP::CmdLine cmd("Data conversion program for NPR data");

        TCLAP::ValueArg<std::string>
            inputAArg("a", "input-a", "Input 'A' directory",
                    true, "", "filename", cmd);
        TCLAP::ValueArg<std::string>
            inputBArg("b", "input-b", "Output directory",
                    true, "", "filename", cmd);
        TCLAP::ValueArg<std::string>
            inputAFormatArg("f", "format-a", "Format for input 'A'",
                    false, "raw_v2", "format", cmd);
        TCLAP::ValueArg<std::string>
            inputBFormatArg("r", "format-b", "Format for input 'B'",
                    false, "raw_v2", "format", cmd);

        cmd.parse(argc, argv);

        input_a_dir = inputAArg.getValue();
        input_b_dir = inputBArg.getValue();
        input_a_format = parse_format(inputAFormatArg.getValue());
        input_b_format = parse_format(inputBFormatArg.getValue());
    }
    catch (TCLAP::ArgException &e) { // catch any exceptions
        std::cerr << "error: " << e.error()
            << " for arg " << e.argId() << std::endl;
    }

    std::cout << "Reading data from '" << input_a_dir << "'"
        << " in format " << input_a_format << std::endl;
    RunData data_A = read_data(input_a_dir, input_a_format);
    std::cout << "Reading data from '" << input_b_dir << "'"
        << " in format " << input_b_format << std::endl;
    RunData data_B = read_data(input_b_dir, input_b_format);

    if (!(data_A.metadata == data_B.metadata)) {
        std::cerr << "Warning: Metadata mismatch" << std::endl;
    }
    for (TrajectoryData &traj_A: data_A.trajectory_data) {
        int idx = -1;
        int trajectory = traj_A.trajectory;
        for (int i = 0; i < data_B.trajectory_data.size(); i++) {
            if (data_B.trajectory_data[i].trajectory == trajectory) {
                idx = i;
                break;
            }
        }
        if (idx == -1) {
            std::cerr << "Warning: trajectory " << trajectory
                << " not found in input B" << std::endl;
            continue;
        }
        TrajectoryData &traj_B = data_B.trajectory_data[idx];

        RealD diff = l2_difference(traj_A, traj_B);
        std::cout << "Trajectory " << trajectory << ": "
            << diff << std::endl;
        if (diff != 0) {
            print_diff_profile(traj_A, traj_B);
        }
    }

    return 0;
}

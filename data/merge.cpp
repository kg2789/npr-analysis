#include <algorithm>
#include <string>

#include <tclap/CmdLine.h>

#include "format.h"
#include "io.h"
#include "../data.h"

static void print_usage() {
    std::ostream &os = std::cerr;
    os << "Usage: ./merge out out_format in1 in1_format in2 in2_format ..." << std::endl;
}

static void die_usage(int exit_code = 1) {
    print_usage();
    exit(exit_code);
}

Metadata merge_metadata(const Metadata &m1, const Metadata &m2) {
    for (int trajA: m1.trajectories) {
        for (int trajB: m2.trajectories) {
            if (trajA == trajB) {
                std::cerr << "Error: both input and output contain trajectory "
                    << trajA << std::endl;
                exit(-1);
            }
        }
    }

    // Merge all trajectories together into one vector
    std::vector<int> all_trajectories = m1.trajectories;
    all_trajectories.insert(all_trajectories.end(),
            m2.trajectories.begin(),
            m2.trajectories.end());
    std::sort(all_trajectories.begin(), all_trajectories.end());

    // Check that the metadata agree aside from the trajectories
    Metadata merged1 = m1;
    Metadata merged2 = m2;
    merged1.trajectories = all_trajectories;
    merged2.trajectories = all_trajectories;

    if (!(merged1 == merged2)) {
        std::cerr << "Warning: metadata mismatch" << std::endl;
    }
    return merged1;
}

// The point of this program is to merge together datasets which were both
// measured on the same lattice with the same parameters but on different
// trajectories. There are two inputs, called the 'input' and 'output' dataset.
// If the metadata of these datasets checks out, then the program will convert
// the input dataset into the format of the output dataset and write the input
// dataset into the output.
int main(int argc, char **argv) {
    std::cout << "argc = " << argc << std::endl;

    if (argc < 5) {
        die_usage(1);
    }

    if (argc % 2 != 1) {
        std::cerr << "Error: odd number of arguments" << std::endl;
        die_usage(1);
    }

    std::string output_dir = argv[1];
    Format output_format = parse_format(argv[2]);

    std::vector<std::string> in_dirs;
    std::vector<Format> in_format;
    for (int i = 3; i < argc; i += 2) {
        in_dirs.push_back(argv[i]);
        in_format.push_back(parse_format(argv[i + 1]));
    }
    assert(in_dirs.size() == in_format.size());
    const int num_input = in_dirs.size();

    // Read and merge metadata
    std::cout << "Reading metadata from '" << in_dirs[0] << "'" << std::endl;
    Metadata metadata = read_metadata(in_dirs[0] + "/metadata.xml");
    for (int i = 1; i < num_input; i++) {
        std::cout << "Reading metadata from '" << in_dirs[i] << "'" << std::endl;
        Metadata new_metadata = read_metadata(in_dirs[i] + "/metadata.xml");
        metadata = merge_metadata(metadata, new_metadata);
    }

    std::cout << "Done merging metadata; now reading data" << std::endl;

    std::vector<RunData> data;
    data.reserve(in_dirs.size());
    for (int i = 0; i < in_dirs.size(); i++) {
        std::cout << "Reading data from '" << in_dirs[i] << "' in format "
            << in_format[i] << std::endl;
        data.push_back(read_data(in_dirs[i], in_format[i], false));
    }

    RunData output;
    output.metadata = metadata;
    output.trajectory_data.reserve(metadata.trajectories.size());
    for (RunData &run: data) {
        output.trajectory_data.insert(output.trajectory_data.end(),
                run.trajectory_data.begin(),
                run.trajectory_data.end());
    }
    std::sort(output.trajectory_data.begin(),
            output.trajectory_data.end(),
            [] (TrajectoryData& a, TrajectoryData& b) {
            return a.trajectory < b.trajectory;
            });

    std::cout << "Writing data to '" << output_dir << "'"
        << "in format " << output_format << std::endl;
    write_data(output_dir, output, output_format);
}

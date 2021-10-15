#include <cerrno>

#include "io.h"

#ifndef RABBOTT_READER
#define RABBOTT_READER Hdf5Reader
#endif
#ifndef RABBOTT_WRITER
#define RABBOTT_WRITER Hdf5Writer
#endif


static std::string short_name(const SpinStructure &s) {
    switch (s) {
        case SpinStructure::LEFT:   return "L";
        case SpinStructure::RIGHT:  return "R";
        case SpinStructure::VECTOR: return "V";
        case SpinStructure::AXIAL:  return "A";
        default:
            std::cerr << "Error: incorrect SpinStructure '" << s
                << "' supplied to short_name in rabbott.cpp" << std::endl;
            exit(1);
    }
}

static std::string color_structure_name(const ColorStructure& col) {
    switch (col) {
        case ColorStructure::COLOR_DIAG: return "_DIAG";
        case ColorStructure::COLOR_MIXED: return "_MIXED";
        default:
            std::cerr << "Error: incorrect ColorStructure '" << col
                << "' supplied to color_structure_name in rabbott.cpp"
                << std::endl;
            exit(1);
    }
}

template <typename Tval, typename T>
static void rw_fourq_with_structure(bool write_, const std::string& dir,
        T &data, int trajectory, LoopStructure loop, ColorStructure color,
        SpinStructure s1, SpinStructure s2) {
    using namespace Grid;
    std::string subdir = dir + "/" + short_name(s1) + short_name(s2)
        + color_structure_name(color);
    std::string filename = subdir + "/" + std::to_string(trajectory) + ".bin";

    if (write_) {
        make_dir(dir);
        make_dir(subdir);
    }

    Tval &value = data.get_measurement(loop, color, s1, s2);

    if (write_) {
        RABBOTT_WRITER writer(filename);
        write(writer, "data", value);
    }
    else {
        RABBOTT_READER reader(filename);
        read(reader, "data", value);
    }
}

static const std::vector<SpinStructure> spin_structures = {
    SpinStructure::VECTOR, SpinStructure::AXIAL
};
static const std::vector<ColorStructure> color_structures = {
    ColorStructure::COLOR_DIAG, ColorStructure::COLOR_MIXED
};

static const std::vector<LoopStructure> fourq_loop_structures = {
    LoopStructure::FULLY_CONNECTED, LoopStructure::CONNECTED_LOOP,
    LoopStructure::DISCONNECTED_LOOP
};

static const std::vector<std::string> fourq_subdirs = {
    "fully_connected", "connected_loop", "disconnected_loop"
};


// General structure of the 4-quark directories:
// data (top-level directory)
//   |-- fourq_op ('dir' in the function below)
//       |-- fourq_ext
//          |-- fully_connected
//          |-- connected_loop
//          |-- disconnected_loop
//       |-- twoq_ext
//          |-- connected_loop
//          |-- disconnected_loop
//
// In each instance of fully_connnected, connected_loop, disconnected_loop, we
// have the following structure:
// fully_connected
//     |-- VV_MIXED
//     |-- VA_MIXED
//     |-- AV_MIXED
//     |-- AA_MIXED
//     |-- VV_DIAG
//     |-- VA_DIAG
//     |-- AV_DIAG
//     |-- AA_DIAG
//
// Each of the lowest-level subdirs (e.g. VV_MIXED) contains files of the form
// '1500.bin' which contain the for the given trajectory in binary format.

static void rw_fourq_op(bool write, const std::string& dir,
        TrajectoryData *data, int trajectory) {
    if (write) {
        make_dir(dir);
        make_dir(dir + "/fourq_ext");
        make_dir(dir + "/twoq_ext");
    }

    assert(fourq_loop_structures.size() == fourq_subdirs.size());
    for (int i = 0; i < fourq_loop_structures.size(); i++) {
        LoopStructure loop = fourq_loop_structures[i];
        for (SpinStructure s1: spin_structures) {
            for (SpinStructure s2: spin_structures) {
                for (ColorStructure col: color_structures) {
                    // Four-quark external states
                    rw_fourq_with_structure<SpinColourSpinColourMatrix>(write,
                            dir + "/fourq_ext/" + fourq_subdirs[i],
                            data->fourq_op.fourq_ext, trajectory,
                            loop, col, s1, s2);
                    if (loop == LoopStructure::FULLY_CONNECTED) {
                        // There are no fully connected diagrams with two-quark
                        // external states
                        continue;
                    }
                    // Two-quark external states
                    rw_fourq_with_structure<SpinColourMatrix>(write,
                            dir + "/twoq_ext/" + fourq_subdirs[i],
                            data->fourq_op.twoq_ext, trajectory,
                            loop, col, s1, s2);
                }
            }
        }
    }
}

template <typename T>
static void rw_individual_twoq_op(bool write_, const std::string& dir,
        T& value, int trajectory) {
    using namespace Grid;
    std::string filename = dir + "/" + std::to_string(trajectory) + ".bin";

    if (write_) {
        make_dir(dir);
        RABBOTT_WRITER writer(filename);
        write(writer, "data", value);
    }
    else {
        RABBOTT_READER reader(filename);
        read(reader, "data", value);
    }
}

std::vector<TwoQuarkOp> twoq_operators= {
    TwoQuarkOp::SCALAR, TwoQuarkOp::DSLASH_LEFT, TwoQuarkOp::DSLASH_RIGHT
};
std::vector<std::string> twoq_dir_names = {
    "scalar", "dslash_left", "dslash_right"
};

// General structure of the 2-quark directories:
// data (top-level directory)
//   |-- twoq_op (dir in the function below)
//       |-- twoq_ext
//           |-- scalar
//           |-- dslash_left
//           |-- dslash_right
//           |-- scalar_gamma5
//           |-- dslash_left_gamma5
//           |-- dslash_right_gamma5
//       |-- fourq_ext
//           |-- scalar
//           |-- dslash_left
//           |-- dslash_right
//           |-- scalar_gamma5
//           |-- dslash_left_gamma5
//           |-- dslash_right_gamma5
//
// Each of the lowest-level directories has the same format as with the 4-quark
// operators, containing files of the form '1500.bin'

static void rw_twoq_op(bool write, const std::string& dir,
        TrajectoryData *data, int trajectory) {
    assert(twoq_operators.size() == twoq_dir_names.size());
    const int num_ops = twoq_operators.size();
    for (int i = 0; i < num_ops; i++) {
        TwoQuarkOp op = twoq_operators[i];
        std::string name = twoq_dir_names[i];
        std::string name_gamma5 = name + "_gamma5";

        if (write) {
            make_dir(dir);
            make_dir(dir + "/twoq_ext");
            make_dir(dir + "/fourq_ext");
        }

        rw_individual_twoq_op(write, dir + "/twoq_ext/" + name,
                data->twoq_op.twoq_ext.values[(int) op], trajectory);
        rw_individual_twoq_op(write, dir + "/twoq_ext/" + name_gamma5,
                data->twoq_op.twoq_ext.psuedoscalar_values[(int) op], trajectory);

        rw_individual_twoq_op(write, dir + "/fourq_ext/" + name,
                data->twoq_op.fourq_ext.values[(int) op], trajectory);
        rw_individual_twoq_op(write, dir + "/fourq_ext/" + name_gamma5,
                data->twoq_op.fourq_ext.psuedoscalar_values[(int) op], trajectory);
    }
}

static void rw_external_leg(bool write_, const std::string& dir,
        SpinColourMatrix& leg, int trajectory) {
    std::string filename = dir + "/" + std::to_string(trajectory) + ".bin";

    if (write_) {
        make_dir(dir);
        RABBOTT_WRITER writer(filename);
        write(writer, "data", leg);
    }
    else {
        RABBOTT_READER reader(filename);
        read(reader, "data", leg);
    }
}

static void rw_external_legs(bool write, const std::string& dir,
        TrajectoryData* data, int trajectory) {
    if (write) {
        make_dir(dir);
        make_dir(dir + "/in");
        make_dir(dir + "/out");
    }
    rw_external_leg(write, dir + "/in", data->external_legs.Sin_average, trajectory);
    rw_external_leg(write, dir + "/out", data->external_legs.Sout_average, trajectory);
}

static void rw_g1(bool write_, const std::string &dir,
        TrajectoryData *data, int trajectory) {
    if (write_)
        make_dir(dir);

    std::string filename = dir + "/G1." + std::to_string(trajectory) + ".h5";
    if (write_) {
        RABBOTT_WRITER writer(filename);
        write(writer, "G1", data->g1);
    }
    else {
        RABBOTT_READER reader(filename);
        read(reader, "G1", data->g1);
    }
}

static void rw_current(bool write_, const std::string &dir,
        TrajectoryData *data, int trajectory) {
    if (write_)
        make_dir(dir);

    std::string filename = dir + "/" + std::to_string(trajectory) + ".h5";
    if (write_) {
        RABBOTT_WRITER writer(filename);
        write(writer, "current", data->current);
    }
    else {
        RABBOTT_READER reader(filename);
        read(reader, "current", data->current);
    }
}

static void rw_trajectory(bool write, const std::string& dir,
        TrajectoryData* data, int trajectory) {
    if (write) make_dir(dir);
    //// Read/write data for four-quark opeartors
    std::string fourq_dir = dir + "/fourq_op";
    rw_fourq_op(write, fourq_dir, data, trajectory);

    //// Read/write data for two-quark opeartors
    std::string twoq_dir = dir + "/twoq_op";
    rw_twoq_op(write, twoq_dir, data, trajectory);

    //// Read/write external legs
    std::string external_leg_dir = dir + "/external_legs";
    rw_external_legs(write, external_leg_dir, data, trajectory);

    //// Read/write G1 measurments
    std::string g1_dir = dir + "/G1";
    rw_g1(write, g1_dir, data, trajectory);

    //// Read/write vector/axial current measurements
    std::string current_dir = dir + "/current";
    rw_current(write, current_dir, data, trajectory);
}

void write_trajectory_rabbott(const std::string& dir, TrajectoryData* data,
        int trajectory) {
    rw_trajectory(true, dir, data, trajectory);
}

void read_trajectory_rabbott(const std::string& dir, TrajectoryData* data,
        int trajectory) {
    rw_trajectory(false, dir, data, trajectory);
}

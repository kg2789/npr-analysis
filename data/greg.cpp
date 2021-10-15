#include <sstream>
#include "io.h"
#include "greg/DoubleWilsonMatrix.h"
#include "greg/convert.hpp"

static bool file_exists(const std::string &filename) {
    std::ifstream f(filename.c_str());
    return f.good();
}

static void rw_wilson_matrix(bool write, std::string filestem,
        SpinColourMatrix &mat) {
    std::string filename = filestem + ".dat";
    if (write) {
        WilsonMatrix wm = sc_to_wilson(mat);
        wm.WriteToFile(filename.c_str());
    }
    else {
        WilsonMatrix wm(filename.c_str());
        mat = wilson_to_sc(wm);
    }
}

bool EndsWithBin(const char* str);
// Greg's code had a problem with endianness: namely, BG/Q has different
// endianness from qcdserver, so simply reading the binary file wouldn't work.
// To solve this his code reverses endianness whenever it detects when it's
// reading a binary file by checking if the extension is '.bin'. Consequently we must also reverse endianness here.
//
// Only DoubleWilsonMatrix's are ever written to binary files, so we don't need
// to implement the same check above
static void write_double_wilson_matrix(std::string &filename,
        SpinColourSpinColourMatrix &mat) {
    DoubleWilsonMatrix dwm = scsc_to_double_wilson(mat);
    if (EndsWithBin(filename.c_str())) {
        dwm.WriteToBinaryFileReversingEndianness(filename.c_str());
    }
    else {
        dwm.WriteToFile(filename.c_str());
    }
}

static void read_double_wilson_matrix(std::string &filename,
        SpinColourSpinColourMatrix &mat) {
    DoubleWilsonMatrix dwm;
    if (EndsWithBin(filename.c_str())) {
        dwm.LoadFromBinaryFileReversingEndianness(filename.c_str());
    }
    else {
        dwm.LoadFromFile(filename.c_str());
    }
    mat = double_wilson_to_scsc(dwm);
}

// This is somewhat of a hack to get around the fact that Greg's code seems to // have inconsistent conventions on whether files are stored in '.dat' or
// '.bin' files. Rather than try to figure out these conventions, instead we
// just try everything. For reading this means trying to read with both file
// extensions, while for writing this means creating a separate file for every
// extension. In the latter case we can cheat a bit and use hardlinks to avoid
// taking up the extra disk space.
static void rw_double_wilson_matrix(bool write,
        std::string &filestem, SpinColourSpinColourMatrix &mat) {
    std::vector<std::string> file_extensions = { ".bin", ".dat" };
    if (write) {
        for (std::string &ext: file_extensions) {
            std::string filename = filestem + ext;
            write_double_wilson_matrix(filename, mat);
        }
    }
    else {
        for (std::string &ext: file_extensions) {
            std::string filename = filestem + ext;
            if (file_exists(filename)) {
                read_double_wilson_matrix(filename, mat);
            }
        }
    }
}


// Greg's code expects the momenta to be integers, but our mometa can be
// general doubles. When testing with Greg's code we should therefore enforce
// that the momenta are actually integers, but floating point errors can make
// this difficult; to compensate we just round to the nearest integer
static std::vector<int> lround_vector(const std::vector<RealD> &v) {
    std::vector<int> ret;
    ret.reserve(v.size());
    for (const RealD &x: v) {
        ret.push_back(std::lround(x));
    }
    return ret;
}

static std::vector<int> get_momentum_in(Metadata &metadata) {
    return lround_vector(metadata.p2);
}

static std::vector<int> get_momentum_out(Metadata &metadata) {
    return lround_vector(metadata.p1);
}

static std::string momentum_to_string(std::vector<int> p) {
    std::ostringstream os;
    if (p.size() != 4) {
        std::cerr << "Error: momentum has size " << p.size()
            << " (expected 4)" << std::endl;
        exit(1);
    }
    os << p[0] << "_" << p[1] << "_" << p[2] << "_" << p[3];
    return os.str();
}

static std::string str_momentum_in(Metadata &metadata) {
    return momentum_to_string(get_momentum_in(metadata));
}

static std::string str_momentum_out(Metadata &metadata) {
    return momentum_to_string(get_momentum_out(metadata));
}

static std::string subtraction_op_filename(const std::string &dir,
        std::string external_state_tag,
        Metadata &metadata, int trajectory,
        TwoQuarkOp subtraction_op, Parity parity) {
    std::ostringstream out;
    out << dir << "/" << external_state_tag << "_";

    // We only read/write the data one parity at a time, so it makes no sense
    // to have both parities here
    assert(parity == Parity::EVEN || parity == Parity::ODD);

    // We need to match NPRSettings::subtraction_op_names. The names are given
    // in NPRSettings.h lines 82-91 in Greg's code
    switch(subtraction_op) {
        case TwoQuarkOp::SCALAR:
            if (parity == Parity::ODD) {
                out << "pseudo";
            }
            out << "scalar";
            break;
        case TwoQuarkOp::DSLASH_RIGHT:
            // I'm not entirely certain on the conventions here (what does
            // "forward" mean here?), but in the worst case this would just
            // mean accidentally swapping two of the subtraction operators,
            // which shouldn't change the final result
            out << "forward_covariant_dslash";
            if (parity == Parity::ODD) {
                out << "_g5";
            }
            break;
        case TwoQuarkOp::DSLASH_LEFT:
            out << "backward_covariant_dslash";
            if (parity == Parity::ODD) {
                out << "_g5";
            }
            break;
    }

    out << "_pa" << str_momentum_out(metadata) << "_pb"
        << str_momentum_in(metadata) << "_traj"
        << trajectory;

    return out.str();
}

/* General rule for file extensions: if the output is a WilsonMatrix, then the
 * file is a '.dat', while if the output is a DoubleWilsonMatrix, then the
 * files is a '.bin'. Greg's code seems to follow this convention, so we have
 * to follow it here to be able to interface with Greg's code.  */
static void rw_subtraction_ops(bool write, const std::string &dir,
        TrajectoryData *data, Metadata metadata, int trajectory) {
    for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
        TwoQuarkOp op = (TwoQuarkOp) i;
        // Four-quark external state
        std::string parity_even_fourq_filename
            = subtraction_op_filename(dir, "fourq", metadata,
                    trajectory, op, Parity::EVEN);
        rw_double_wilson_matrix(write, parity_even_fourq_filename,
                data->twoq_op.fourq_ext.values[i]);

        std::string parity_odd_fourq_filename
            = subtraction_op_filename(dir, "fourq", metadata,
                    trajectory, op, Parity::ODD);
        rw_double_wilson_matrix(write, parity_odd_fourq_filename,
                data->twoq_op.fourq_ext.psuedoscalar_values[i]);

        // Two-quark external state
        std::string parity_even_twoq_filename
            = subtraction_op_filename(dir, "twoq", metadata,
                    trajectory, op, Parity::EVEN);
        rw_wilson_matrix(write, parity_even_twoq_filename,
                data->twoq_op.twoq_ext.values[i]);

        std::string parity_odd_twoq_filename
            = subtraction_op_filename(dir, "twoq", metadata,
                    trajectory, op, Parity::ODD);
        rw_wilson_matrix(write, parity_odd_twoq_filename,
                data->twoq_op.twoq_ext.psuedoscalar_values[i]);
    }
}


static std::string fourq_op_filename(const std::string &dir,
        std::string external_state_tag,
        Metadata &metadata, int trajectory,
        LoopStructure loop, ColorStructure color,
        SpinStructure gammaA, SpinStructure gammaB) {
    std::ostringstream out;
    out << dir << "/" << external_state_tag << "_";

    // See Greg's code, util.cpp line 14 (definition of variable
    // diagram_names[]). We need to replicate the names given there to match
    // the given loop & color structure here
    switch (loop) {
        case LoopStructure::FULLY_CONNECTED:
            out << "fully_connected";
            break;
        case LoopStructure::DISCONNECTED_LOOP:
            out << "disconnected_loop";
            break;
        case LoopStructure::CONNECTED_LOOP:
            out << "connected_loop";
            break;
        default:
            std::cerr << "Error: invalid loop structure "
                << loop << " while reading/writing in format GREG"
                << std::endl;
            exit(1);
    }

    out << "_";

    switch (color) {
        case ColorStructure::COLOR_DIAG:
            out << "color_diag";
            break;
        case ColorStructure::COLOR_MIXED:
            out << "color_mixed";
            break;
        default:
            std::cerr << "Error: invalid color structure "
                << color << " while reading/writing in format GREG"
                << std::endl;
            exit(1);
    }

    // Writing gammaA, gammaB to the stream should output V for vector
    // structure and A for axial (see variable VA_names[] in Greg's code). As
    // of the writing of this comment, that is what operator << does for
    // SpinStructures
    out << "_" << gammaA << "_" << gammaB << "_pa"
        << str_momentum_out(metadata) << "_pb" << str_momentum_in(metadata)
        << "_traj" << trajectory;

    return out.str();
}

static void rw_fourq_op(bool write, const std::string &dir,
        TrajectoryData *data, Metadata metadata, int trajectory) {
    for (int a = 0; a < NUM_LOOP_STRUCTURES; a++) {
        for (int b = 0; b < NUM_COLOR_STRUCTURES; b++) {
            for (int c = 0; c < NUM_SPIN_STRUCTURES; c++) {
                for (int d = 0; d < NUM_SPIN_STRUCTURES; d++) {
                    LoopStructure loop = (LoopStructure) a;
                    ColorStructure color = (ColorStructure) b;
                    SpinStructure gammaA = (SpinStructure) c;
                    SpinStructure gammaB = (SpinStructure) d;

                    const SpinStructure L = SpinStructure::LEFT;
                    const SpinStructure R = SpinStructure::RIGHT;

                    if (gammaA == L || gammaB == L
                            || gammaA == R || gammaB == R) {
                        // We don't want to read/write diagrams with left/right
                        // structures; we compute the left/right structures
                        // when we need them
                        continue;
                    }

                    // Four-quark external states
                    std::string fourq_ext_filename
                       = fourq_op_filename(dir, "fourq",
                               metadata, trajectory,
                               loop, color, gammaA, gammaB);
                    rw_double_wilson_matrix(write, fourq_ext_filename,
                            data->fourq_op.fourq_ext.get_measurement(loop, color,
                                gammaA, gammaB));

                    // Two-quark external states
                    if (loop == LoopStructure::FULLY_CONNECTED) {
                        // There are no fully-connected diagrams with a 2-quark
                        // external state
                        continue;
                    }
                    std::string twoq_ext_filename
                        = fourq_op_filename(dir, "twoq",
                                metadata, trajectory,
                                loop, color, gammaA, gammaB);
                    rw_wilson_matrix(write, twoq_ext_filename,
                            data->fourq_op.twoq_ext.get_measurement(loop, color,
                                gammaA, gammaB));
                }
            }
        }
    }
}

static std::string external_leg_filename(const std::string& dir,
       std::string momentum, int trajectory) {
    std::ostringstream ret;
    ret << dir << "/external_leg_p" << momentum
        << "_traj" << trajectory;
    return ret.str();
}

static void rw_external_legs(bool write, const std::string& dir,
        TrajectoryData *data, Metadata &metadata, int trajectory) {
    std::string Sin_filename
        = external_leg_filename(dir, str_momentum_in(metadata), trajectory);
    rw_wilson_matrix(write, Sin_filename, data->external_legs.Sin_average);
    std::string Sout_filename
        = external_leg_filename(dir, str_momentum_out(metadata), trajectory);
    rw_wilson_matrix(write, Sout_filename, data->external_legs.Sout_average);
}

static std::string g1_filename(const std::string &dir,
        std::string num_external_quarks, std::string c1_str, bool gamma5,
        std::string pa, std::string pb, int trajectory) {
    std::ostringstream ret;
    ret << dir << "/" << num_external_quarks << "_old_G1";
    if (gamma5)
        ret << "g5";
    ret << c1_str;
    ret << "_pa" << pa << "_pb" << pb
        << "_traj" << trajectory;
    return ret.str();
}

static void rw_g1(bool write, const std::string &dir,
        TrajectoryData *data, Metadata &metadata, int trajectory) {
    std::string pa = str_momentum_out(metadata);
    std::string pb = str_momentum_in(metadata);
    std::string filename;

    //std::string c1_str = "_c1_-0.331";
    std::string c1_str = "";
    filename = g1_filename(dir, "twoq", c1_str, false, pa, pb, trajectory);
    rw_wilson_matrix(write, filename, data->g1.twoq_scalar);
    filename = g1_filename(dir, "twoq", c1_str, true, pa, pb, trajectory);
    rw_wilson_matrix(write, filename, data->g1.twoq_gamma5);

    filename = g1_filename(dir, "fourq", c1_str, false, pa, pb, trajectory);
    rw_double_wilson_matrix(write, filename, data->g1.fourq_scalar);
    filename = g1_filename(dir, "fourq", c1_str, true, pa, pb, trajectory);
    rw_double_wilson_matrix(write, filename, data->g1.fourq_gamma5);
}

static std::string current_filename(const std::string &dir,
        int mu, bool gamma5, std::string pa, std::string pb, int trajectory) {
    std::ostringstream ret;
    ret << dir << "/twoq_gamma" << mu;
    if (gamma5)
        ret << "_gamma5";
    ret << "_pa" << pa << "_pb" << pb
        << "_traj" << trajectory;
    return ret.str();
}

static void rw_current(bool write, const std::string &dir,
        TrajectoryData *data, Metadata &metadata, int trajectory) {
    std::string pa = str_momentum_out(metadata);
    std::string pb = str_momentum_in(metadata);
    if (!write) {
        data->current.vector.resize(4);
        data->current.axial.resize(4);
    }
    for (int mu = 0; mu < 4; mu++) {
        std::string vector_filename
            = current_filename(dir, mu, false, pa, pb, trajectory);
        rw_wilson_matrix(write, vector_filename, data->current.vector[mu]);
        std::string axial_filename
            = current_filename(dir, mu, true, pa, pb, trajectory);
        rw_wilson_matrix(write, axial_filename, data->current.axial[mu]);
    }
}

static void rw_trajectory(bool write, const std::string& dir,
        TrajectoryData *data, Metadata &metadata, int trajectory) {
    rw_fourq_op(write, dir, data, metadata, trajectory);
    rw_subtraction_ops(write, dir, data, metadata, trajectory);
    rw_external_legs(write, dir, data, metadata, trajectory);
    rw_g1(write, dir, data, metadata, trajectory);
    rw_current(write, dir, data, metadata, trajectory);
}

void write_trajectory_greg(const std::string &dir, TrajectoryData *data,
        Metadata &metadata, int trajectory) {
    rw_trajectory(true, dir, data, metadata, trajectory);
}

void read_trajectory_greg(const std::string &dir, TrajectoryData *data,
        Metadata &metadata, int trajectory) {
    rw_trajectory(false, dir, data, metadata, trajectory);
}

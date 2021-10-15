#include "io.h"
#include "fierz.h"
#include <fstream>

using namespace Grid;

// Taken directly from measurement code with a modified name
class SubtractionResult: Serializable
{
    public:
        // Contains information on both the 2 and 4-quark external state
        // diagrams with the given subtraction operator
        class OperatorResult : Serializable {
            public:
                GRID_SERIALIZABLE_CLASS_MEMBERS(OperatorResult,
                        SpinColourSpinColourMatrix, fourq,
                        SpinColourMatrix, twoq);
        };
        GRID_SERIALIZABLE_CLASS_MEMBERS(SubtractionResult,
                OperatorResult, dslash_left,
                OperatorResult, dslash_gamma5_left,
                OperatorResult, dslash_right,
                OperatorResult, dslash_gamma5_right,
                OperatorResult, scalar,
                OperatorResult, psuedoscalar);
};

class FourQuarkResult: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(FourQuarkResult,
                std::vector<SpinColourSpinColourMatrix>, fully_connected_fourq,
                std::vector<SpinColourMatrix>, connected_loop_fourq,
                std::vector<SpinColourMatrix>, disconnected_loop_fourq,
                std::vector<SpinColourMatrix>, connected_loop_bilinear,
                std::vector<SpinColourMatrix>, disconnected_loop_bilinear,
                std::vector<SpinColourMatrix>, fully_connected_bilinear,
                SpinColourMatrix, spectator,
                std::vector<Gamma::Algebra>, gammaA,
                std::vector<Gamma::Algebra>, gammaB);
};

class G1Result: Serializable
{
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(G1Result,
                SpinColourMatrix, fourq_scalar,
                SpinColourMatrix, fourq_gamma5,
                SpinColourMatrix, twoq_scalar,
                SpinColourMatrix, twoq_gamma5);
};


static bool file_exists(const std::string &filename) {
    std::ifstream f(filename.c_str());
    return f.good();
}

template <typename T>
T read_result(const std::string &filestem, const std::string &tag,
        int trajectory) {

    std::string full_filestem = filestem + "." + std::to_string(trajectory);
    std::string hdf5_filename = full_filestem + ".h5";
    std::string xml_filename = full_filestem + ".xml";

    T result;
    if (file_exists(hdf5_filename)) {
        Hdf5Reader reader(hdf5_filename);
        read(reader, tag, result);
    }
    else if (file_exists(xml_filename)) {
        XmlReader reader(xml_filename);
        read(reader, tag, result);
    }
    else {
        std::cerr << "Error: file '" << filestem << "' for trajectory "
            << trajectory << " does not exist" << std::endl;
        exit(1);
    }

    return result;
}

void read_g1(const std::string &filename, TrajectoryData *data, int trajectory,
        const SpinColourMatrix &spectator) {
    // Read result from file
    G1Result result = read_result<G1Result>(filename, "G1", trajectory);

    data->g1.twoq_scalar = result.twoq_scalar;
    data->g1.twoq_gamma5 = result.twoq_gamma5;
    tensor_prod(data->g1.fourq_scalar, result.fourq_scalar, spectator);
    tensor_prod(data->g1.fourq_gamma5, result.fourq_gamma5, spectator);
}

void read_subtraction_diagrams(const std::string &filename,
        TrajectoryData *data, int trajectory) {
    //// Read file into memory
    std::string full_filename = filename + "."
        + std::to_string(trajectory) + ".h5";

    SubtractionResult *result = new SubtractionResult();
    *result = read_result<SubtractionResult>(filename, "subtraction_ops", trajectory);

    //// Move data into output data structure
    auto convert_result = [&](SubtractionResult::OperatorResult &res,
            TwoQuarkOp subtraction_op, bool psuedoscalar) {
        const int idx = (int) subtraction_op;
        if (psuedoscalar) {
            data->twoq_op.twoq_ext.psuedoscalar_values[idx] = res.twoq;
            data->twoq_op.fourq_ext.psuedoscalar_values[idx] = res.fourq;
        }
        else {
            data->twoq_op.twoq_ext.values[idx] = res.twoq;
            data->twoq_op.fourq_ext.values[idx] = res.fourq;
        }
    };

    convert_result(result->dslash_left, TwoQuarkOp::DSLASH_LEFT, false);
    convert_result(result->dslash_gamma5_left, TwoQuarkOp::DSLASH_LEFT, true);

    convert_result(result->dslash_right, TwoQuarkOp::DSLASH_RIGHT, false);
    convert_result(result->dslash_gamma5_right, TwoQuarkOp::DSLASH_RIGHT, true);

    convert_result(result->scalar, TwoQuarkOp::SCALAR, false);
    convert_result(result->psuedoscalar, TwoQuarkOp::SCALAR, true);

    delete result;
}

const std::vector<SpinStructure> spin_structures = {
    SpinStructure::VECTOR, SpinStructure::AXIAL
};

template <typename T>
static void read_loop_info(const std::string& dir,
        std::string subdir, std::string tag,
        int trajectory, std::vector<T>& out) {
    std::string filestem = dir + "/" + subdir + "/" + tag;
    out = read_result<std::vector<T>>(filestem, tag, trajectory);
}

// Helper function for process_fourq_results above; find_gammas is a function
// which takes in 2 'Gammas' and gives us the corresponding index in vals
template <typename T, typename F>
void compute_va_structure(T &lret, std::vector<T> &vals,
        SpinStructure gammaA, SpinStructure gammaB, F find_gammas) {
    lret = Zero();
    assert(gammaA == SpinStructure::VECTOR || gammaA == SpinStructure::AXIAL);
    assert(gammaB == SpinStructure::VECTOR || gammaB == SpinStructure::AXIAL);
    for (int mu = 0; mu < Nd; mu++) {
        Gamma gmu = Gamma::gmu[mu];
        Gamma gmug5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
        Gamma gammaA_mu = (gammaA == SpinStructure::VECTOR) ? gmu : gmug5;
        Gamma gammaB_mu = (gammaB == SpinStructure::VECTOR) ? gmu : gmug5;

        lret += vals[find_gammas(gammaA_mu, gammaB_mu)];
    }
}

// FourQuarkResult has contains all of the relevant info to compute the 4-quark
// diagrams, but not in the format we need. We need structures like \sum_\mu
// (qbar \gamma^\mu q) (qbar \gamma_\mu q), but the measurment code does not
// perform the sum over mu; we perform that sum here, along with putting the
// result in the correct spot.
//
// The parameter 'color' tells us where to put the output; usually this should
// be COLOR_DIAG, except when we're testing against explicitly computed
// color-mixed diagrams
//
// The parameter 'apply_fierz' tells us whether to apply the Fierz identity to
// obtain the color-mixed strucutres. Setting this to true only makes sense if
// color == COLOR_DIAG, as otherwise we won't even be filling in the data
// needed to apply the Fierz identity.
void process_fourq_results(FourQuarkResult *result,
        TrajectoryData *data, ColorStructure color, bool apply_fierz) {
    // Sanity check (see above)
    assert(!apply_fierz || color == ColorStructure::COLOR_DIAG);

    std::vector<Gamma> gammaA;
    std::vector<Gamma> gammaB;
    assert(result->gammaA.size() == result->gammaB.size());
    for (int i = 0; i < result->gammaA.size(); i++) {
        gammaA.push_back(Gamma(result->gammaA[i]));
        gammaB.push_back(Gamma(result->gammaB[i]));
    }

    // We could make this a hashmap or somethething else more efficient, but gammaA
    // and gammaB should be so small (at most 256 elements) that this might even be
    // faster
    auto find_gammas = [&] (const Gamma &ga, const Gamma &gb) {
        for (int i = 0; i < gammaA.size(); i++) {
            if (gammaA[i].g == ga.g && gammaB[i].g == gb.g) {
                return i;
            }
        }
        std::cerr << "Error reading raw data: could not find measurment"
            << "with gammaA = " << ga.g << " and gammaB = " << gb.g << std::endl;
        exit(1);
    };

    //// Reading matrix elements with 4-quark external states
    auto read_fourq_loop_structures
        = [&] (std::vector<SpinColourMatrix>& vals, LoopStructure loop) {
        for (SpinStructure gammaA: va_spin_structures) {
            for (SpinStructure gammaB: va_spin_structures) {
                SpinColourMatrix tmp;
                compute_va_structure(tmp, vals, gammaA, gammaB, find_gammas);
                SpinColourSpinColourMatrix &out =
                    data->fourq_op.fourq_ext.get_measurement(loop, color,
                            gammaA, gammaB);
                tensor_prod(out, tmp, result->spectator);
            }
        }
    };

    read_fourq_loop_structures(result->connected_loop_fourq, LoopStructure::CONNECTED_LOOP);
    read_fourq_loop_structures(result->disconnected_loop_fourq, LoopStructure::DISCONNECTED_LOOP);

    // Reading fully-connected diagrams
    for (SpinStructure gammaA: va_spin_structures) {
        for (SpinStructure gammaB: va_spin_structures) {
            const LoopStructure loop = LoopStructure::FULLY_CONNECTED;
            SpinColourSpinColourMatrix &out =
                data->fourq_op.fourq_ext.get_measurement(loop, color,
                        gammaA, gammaB);
            compute_va_structure(out, result->fully_connected_fourq,
                    gammaA, gammaB, find_gammas);
        }
    }

    if (apply_fierz) {
        // Fills in the color mixed diagrams using the Fierz identity
        apply_fierz_fourq_ext_loops(data->fourq_op.fourq_ext,
                result->connected_loop_fourq, result->disconnected_loop_fourq,
                result->spectator, find_gammas);
        apply_fierz_fully_connected(data->fourq_op.fourq_ext,
                result->fully_connected_fourq, find_gammas);
    }

    //// Reading matrix elements with 2-quark external states

    auto read_twoq_structures = [&] (std::vector<SpinColourMatrix>& vals,
            LoopStructure loop) {
        for (SpinStructure gammaA: va_spin_structures) {
            for (SpinStructure gammaB: va_spin_structures) {
                SpinColourMatrix tmp;
                compute_va_structure(tmp, vals, gammaA, gammaB, find_gammas);
                data->fourq_op.twoq_ext
                    .get_measurement(loop, color, gammaA, gammaB) = tmp;
            }
        }
    };

    read_twoq_structures(result->connected_loop_bilinear, LoopStructure::CONNECTED_LOOP);
    read_twoq_structures(result->disconnected_loop_bilinear, LoopStructure::DISCONNECTED_LOOP);

    if (apply_fierz) {
        apply_fierz_twoq_ext_loops(data->fourq_op.twoq_ext,
                result->connected_loop_bilinear,
                result->disconnected_loop_bilinear,
                find_gammas);
    }

    //// Extract vector and axial currents from computed bilinears
    auto find_bilinear = [&] (Gamma gA) {
        for (int i = 0; i < gammaA.size(); i++) {
            if (gammaA[i].g == gA.g) {
                return result->fully_connected_bilinear[i];
            }
        }
        std::cerr << "Error reading raw data: could not find bilinear "
            "measurment with gammaA = " << gA.g << std::endl;
        exit(1);
    };

    // Vector current J_V^\mu = <qbar \gamma^\mu q>
    std::vector<Gamma> vector_structure = {
        Gamma(Gamma::Algebra::GammaX),
        Gamma(Gamma::Algebra::GammaY),
        Gamma(Gamma::Algebra::GammaZ),
        Gamma(Gamma::Algebra::GammaT)
    };

    data->current.vector.clear();
    std::transform(vector_structure.begin(), vector_structure.end(),
            std::back_inserter(data->current.vector),
            find_bilinear);

    // Axial current J_A^\mu = <qbar \gamma^\mu \gamma^5 q>
    std::vector<Gamma> axial_structure = {
        Gamma(Gamma::Algebra::GammaXGamma5),
        Gamma(Gamma::Algebra::GammaYGamma5),
        Gamma(Gamma::Algebra::GammaZGamma5),
        Gamma(Gamma::Algebra::GammaTGamma5)
    };

    data->current.axial.clear();
    std::transform(axial_structure.begin(), axial_structure.end(),
            std::back_inserter(data->current.axial),
            find_bilinear);
}

void read_fourq_mixed_diagrams(const std::string &filename, TrajectoryData *data,
        int trajectory) {
    FourQuarkResult *result = new FourQuarkResult;
    *result = read_result<FourQuarkResult>(filename, "fourq_mixed_diagrams", trajectory);

    process_fourq_results(result, data, ColorStructure::COLOR_MIXED, false);
    delete result;
}

// Reads in any diagram containing a four-quark operator to 'data'. Includes
// color-mixed diagrams using the Fierz identity
void read_fourq_diagrams(const std::string& filename, TrajectoryData* data,
        int trajectory) {
    FourQuarkResult *result = new FourQuarkResult();
    *result = read_result<FourQuarkResult>(filename, "fourq_diagrams", trajectory);

    process_fourq_results(result, data, ColorStructure::COLOR_DIAG, true);
    delete result;
}

static void read_external_leg(const std::string &filename,
        SpinColourMatrix &leg, int trajectory) {
    leg = read_result<SpinColourMatrix>(filename, "external_leg", trajectory);
}

static void read_external_legs(const std::string &dir, TrajectoryData *data,
        int trajectory) {
    read_external_leg((dir + "/in"), data->external_legs.Sin_average, trajectory);
    read_external_leg((dir + "/out"), data->external_legs.Sout_average, trajectory);
}

void read_trajectory_raw(const std::string& dir, TrajectoryData* data,
        int trajectory) {
    //// Read four-quark data
    std::string fourq_filename = dir + "/fourq_diagrams";
    read_fourq_diagrams(fourq_filename, data, trajectory);


    //// Read two-quark data
    std::string subtraction_filename = dir + "/subtraction_ops";
    read_subtraction_diagrams(subtraction_filename, data, trajectory);

    //// Read external legs
    std::string external_leg_filename = dir + "/external_legs";
    read_external_legs(external_leg_filename, data, trajectory);

    // Note: we don't read G1 data to avoid difficulties with the spectator
    // quark (it would be possible to read the data, but most of the data
    // involving G1 should be in raw_v2 anyways)
}

void read_trajectory_raw_with_mixed(const std::string& dir, TrajectoryData* data,
        int trajectory) {
    read_trajectory_raw(dir, data, trajectory);

    // Extra check: if there's an error in the function below, it's possible
    // the color-mixed diagrams will not be overwritten entirely, in which case
    // we could have a false-positive when we compare to the color-mixed case.
    // By zeroing out preemptively, we avoid this possiblity.
    for (int a = 0; a < NUM_LOOP_STRUCTURES; a++) {
        for (int c = 0; c < NUM_SPIN_STRUCTURES; c++) {
            for (int d = 0; d < NUM_SPIN_STRUCTURES; d++) {
                LoopStructure loop = (LoopStructure) a;
                const ColorStructure color = ColorStructure::COLOR_MIXED;
                SpinStructure gammaA = (SpinStructure) c;
                SpinStructure gammaB = (SpinStructure) d;

                if (!data->fourq_op.fourq_ext.is_measurement_stored(gammaA, gammaB)) {
                    continue;
                }

                SpinColourSpinColourMatrix &fourq_ext_val
                    = data->fourq_op.fourq_ext.get_measurement(loop, color,
                            gammaA, gammaB);
                SpinColourMatrix &twoq_ext_val
                    = data->fourq_op.twoq_ext.get_measurement(loop, color,
                            gammaA, gammaB);
                fourq_ext_val = Zero();
                twoq_ext_val = Zero();
            }
        }
    }

    read_fourq_mixed_diagrams((dir + "/fourq_mixed_diagrams"), data, trajectory);
}

/* RAW_V2 -- the original output format is fairly monolithic, and consequently
 * not suited for integration into Grid/Hadrons. This new format is intended to
 * be the output of a more modular program which separates the
 * four-quark-operator diagrams into separate files */

static void read_fourq_diagrams_v2(const std::string &dir, TrajectoryData *data,
        int trajectory, RealD volume);
static SpinColourMatrix read_spectator(const std::string &dir, int trajectory);

void read_trajectory_raw_v2(const std::string &dir, TrajectoryData *data,
        Metadata &metadata, int trajectory) {
    RealD volume = 1.0;
    for (int Lmu: metadata.lattice_dimensions) {
        volume *= Lmu;
    }
    //// Read four-quark data
    std::string fourq_dir = dir + "/fourq_diagrams";
    read_fourq_diagrams_v2(fourq_dir, data, trajectory, volume);

    //// Read two-quark data
    std::string subtraction_filename = dir + "/subtraction_ops";
    read_subtraction_diagrams(subtraction_filename, data, trajectory);

    //// Read external legs
    std::string external_leg_filename = dir + "/external_legs";
    read_external_legs(external_leg_filename, data, trajectory);

    //// Read G1 data
    std::string g1_filename = dir + "/g1";
    SpinColourMatrix spectator = read_spectator(dir + "/fourq_diagrams", trajectory);
    spectator = spectator * volume;
    read_g1(g1_filename, data, trajectory, spectator);
}

// FIXME: This function and make_filename() below could both be merged into
// read_result<T>() above by refactoring the code further below
template<typename T>
std::vector<T> read_vector(const std::string &filename, const std::string &tag) {
    std::vector<T> result;
    std::string hdf5_filename = filename + ".h5";
    std::string xml_filename = filename + ".xml";
    if (file_exists(hdf5_filename)) {
        Hdf5Reader reader(hdf5_filename);
        read(reader, tag, result);
    }
    else if (file_exists(xml_filename)) {
        XmlReader reader(xml_filename);
        read(reader, tag, result);
    }
    else {
        std::cerr << "Error: file '" << filename
            << "' does not exist" << std::endl;
        exit(1);
    }
    return result;
}

static std::string make_filename(const std::string &dir,
        const std::string &filename, int trajectory) {
    return dir + "/" + filename + "." + std::to_string(trajectory);
}

static SpinColourMatrix read_spectator(const std::string &dir, int trajectory) {
    std::string filestem = dir + "/spectator";
    return read_result<SpinColourMatrix>(filestem, "external_leg", trajectory);
}

static void read_gammas(const std::string &dir, int trajectory,
        std::vector<Gamma::Algebra> &gammaA,
        std::vector<Gamma::Algebra> &gammaB) {
    std::string filename;

    filename = make_filename(dir, "fully_connected/gammaA", trajectory);
    gammaA = read_vector<Gamma::Algebra>(filename, "gammaA");
    filename = make_filename(dir, "fully_connected/gammaB", trajectory);
    gammaB = read_vector<Gamma::Algebra>(filename, "gammaB");

    assert (gammaA.size() == gammaB.size());

    // The logic of the code above assumes that every diagram shares a common
    // gammaA/gammaB vector, so we should check that this is the case
    auto check_gammas = [&] (const std::string &subdir) {
        std::vector<Gamma::Algebra> gammaA_tmp, gammaB_tmp;
        filename = make_filename(dir, subdir + "/gammaA", trajectory);
        gammaA_tmp = read_vector<Gamma::Algebra>(filename, "gammaA");
        filename = make_filename(dir, subdir + "/gammaB", trajectory);
        gammaB_tmp = read_vector<Gamma::Algebra>(filename, "gammaB");

        assert(gammaA_tmp == gammaA && gammaB_tmp == gammaB);
    };

    check_gammas("connected_loop");
    check_gammas("disconnected_loop");
}

static void read_fourq_diagrams_v2(const std::string &dir, TrajectoryData *data,
        int trajectory, RealD volume) {
    // Rather than replicate all of the logic above, we avoid some work here by
    // simply rearranging the data into the format used above
    FourQuarkResult *result = new FourQuarkResult();
    read_gammas(dir, trajectory, result->gammaA, result->gammaB);

    std::string filename = make_filename(dir, "fully_connected/fourq", trajectory);
    result->fully_connected_fourq
        = read_vector<SpinColourSpinColourMatrix>(filename, "fully_connected_fourq");

    filename = make_filename(dir, "fully_connected/twoq", trajectory);
    result->fully_connected_bilinear
        = read_vector<SpinColourMatrix>(filename, "fully_connected_twoq");

    filename = make_filename(dir, "connected_loop/fourq", trajectory);
    result->connected_loop_fourq
        = read_vector<SpinColourMatrix>(filename, "connected_loop_fourq");

    filename = make_filename(dir, "connected_loop/twoq", trajectory);
    result->connected_loop_bilinear
        = read_vector<SpinColourMatrix>(filename, "connected_loop_twoq");

    filename = make_filename(dir, "disconnected_loop/fourq", trajectory);
    result->disconnected_loop_fourq
        = read_vector<SpinColourMatrix>(filename, "disconnected_loop_fourq");

    filename = make_filename(dir, "disconnected_loop/twoq", trajectory);
    result->disconnected_loop_bilinear
        = read_vector<SpinColourMatrix>(filename, "disconnected_loop_twoq");

    result->spectator = read_spectator(dir, trajectory);
    // Since we compute the spectator quarks using the external leg module what
    // we have here is an average over the volume. However, we already have a
    // 1.0 / volume factor from the bilinear, so we need to remote this factor
    // from the spectator
    result->spectator = volume * result->spectator;

    process_fourq_results(result, data, ColorStructure::COLOR_DIAG, true);
    delete result;
}

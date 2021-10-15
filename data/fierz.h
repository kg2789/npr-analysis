#include "../data.h"
#include "../util.h"

template <typename T>
class BilinearStructures {
    public:
    std::vector<T> values;
    T& VV() { return values[0]; };
    T& VA() { return values[1]; };
    T& AV() { return values[2]; };
    T& AA() { return values[3]; };
    T& SS() { return values[4]; };
    T& SP() { return values[5]; };
    T& PS() { return values[6]; };
    T& PP() { return values[7]; };

    BilinearStructures() {
        values.resize(8);
        for (T& val: values) {
            val = Zero();
        }
    }
};

inline BilinearStructures<SpinColourSpinColourMatrix>
tensor_prod(const BilinearStructures<SpinColourMatrix> &left,
        const SpinColourMatrix &right) {
    BilinearStructures<SpinColourSpinColourMatrix> ret
        = BilinearStructures<SpinColourSpinColourMatrix>();
    for (int i = 0; i <  left.values.size(); i++) {
        tensor_prod(ret.values[i], left.values[i], right);
    }
    return ret;
}

template <typename T, typename F>
static inline BilinearStructures<T>
compute_bilinear_structures(std::vector<T> &vals, F find_gammas) {
    BilinearStructures<T> ret = BilinearStructures<T>();
    ret.VV() = Zero();
    ret.AA() = Zero();
    ret.AV() = Zero();
    ret.VA() = Zero();
    ret.SS() = Zero();
    ret.PP() = Zero();
    ret.SP() = Zero();
    ret.PS() = Zero();

    for (int mu = 0; mu < Nd; mu++) {
        Gamma gmu = Gamma::gmu[mu];
        Gamma gmu_g5 = Gamma::mul[gmu.g][Gamma::Algebra::Gamma5];
        ret.VV() += vals[find_gammas(gmu, gmu)];
        ret.VA() += vals[find_gammas(gmu, gmu_g5)];
        ret.AV() += vals[find_gammas(gmu_g5, gmu)];
        ret.AA() += vals[find_gammas(gmu_g5, gmu_g5)];
    }
    Gamma identity = Gamma(Gamma::Algebra::Identity);
    Gamma g5 = Gamma(Gamma::Algebra::Gamma5);
    ret.SS() = vals[find_gammas(identity, identity)];
    ret.SP() = vals[find_gammas(identity, g5)];
    ret.PS() = vals[find_gammas(g5, identity)];
    ret.PP() = vals[find_gammas(g5, g5)];

    return ret;
}

// Uses the Fierz identity to compute the off-diagonal color diagrams. For the
// loop diagrams this can change the toplogical structure of the diagram, so we
// allow the results to be written to an arbitrary loop structure.
//
// Reference for the parity-odd version of the Fierz identity:
// https://arxiv.org/abs/hep-ph/0306087
//
// T should be one of FourQOpFourQExtMeasurements or FourQOpTwoQExtMeasurements
// and S should be SpinColourMatrix or SpinColourSpinColourMatrix, respectively
template <typename T, typename Tval>
void __apply_fierz(T &measurements,
        BilinearStructures<Tval> &input_values,
        LoopStructure output_structure) {
    using namespace Grid;

    Tval& VV = input_values.VV();
    Tval& AA = input_values.AA();
    Tval& AV = input_values.AV();
    Tval& VA = input_values.VA();
    Tval& SS = input_values.SS();
    Tval& PP = input_values.PP();
    Tval& SP = input_values.SP();
    Tval& PS = input_values.PS();

    SpinStructure V = SpinStructure::VECTOR;
    SpinStructure A = SpinStructure::AXIAL;
    LoopStructure loop = output_structure;
    ColorStructure mixed = ColorStructure::COLOR_MIXED;

    // There are 2 minus signs at play here:
    //
    // 1. Typically the Fierz identities are given for the case of commuting
    // spinors, but we have anti-commuting spinors
    //
    // 2. By applying the Fierz identity and swapping around the quarks we
    // change the topology of the graph, potentially introducing a minus sign
    // (in fact, for the cases we consider this always occurs)
    //
    // To account for (1) the values in parenthesis are the opposite of the
    // values cited in the paper, and to account for (2) we give an overall
    // minus sign
    measurements.get_measurement(loop, mixed, V, V) = -(-SS + 0.5 * VV + 0.5 * AA + PP);
    measurements.get_measurement(loop, mixed, A, A) = -(SS + 0.5 * VV + 0.5 * AA - PP);

    measurements.get_measurement(loop, mixed, V, A) = -(SP + 0.5 * VA + 0.5 * AV - PS);
    measurements.get_measurement(loop, mixed, A, V) = -(-SP + 0.5 * VA + 0.5 * AV + PS);

}

template <typename T, typename F>
void apply_fierz_fourq_ext_loops(T &measurements,
        std::vector<SpinColourMatrix> &connected_loop,
        std::vector<SpinColourMatrix> &disconnected_loop,
        SpinColourMatrix &spectator,
        F find_gammas) {
    BilinearStructures<SpinColourSpinColourMatrix> connected_loop_structures
        = tensor_prod(compute_bilinear_structures(connected_loop, find_gammas),
                spectator);
    __apply_fierz(measurements, connected_loop_structures,
            LoopStructure::DISCONNECTED_LOOP);
    BilinearStructures<SpinColourSpinColourMatrix> disconnected_loop_structures
        = tensor_prod(compute_bilinear_structures(disconnected_loop, find_gammas),
                spectator);
    __apply_fierz(measurements, disconnected_loop_structures,
            LoopStructure::CONNECTED_LOOP);
}

template <typename T, typename F>
void apply_fierz_twoq_ext_loops(T &measurements,
        std::vector<SpinColourMatrix> &connected_loop,
        std::vector<SpinColourMatrix> &disconnected_loop,
        F find_gammas) {
    BilinearStructures<SpinColourMatrix> connected_loop_structures
        = compute_bilinear_structures(connected_loop, find_gammas);
    __apply_fierz(measurements, connected_loop_structures,
            LoopStructure::DISCONNECTED_LOOP);
    BilinearStructures<SpinColourMatrix> disconnected_loop_structures
        = compute_bilinear_structures(disconnected_loop, find_gammas);
    __apply_fierz(measurements, disconnected_loop_structures,
            LoopStructure::CONNECTED_LOOP);
}

template <typename T, typename Tval, typename F>
void apply_fierz_fully_connected(T &measurements,
        std::vector<Tval> &values,
        F find_gammas) {
    BilinearStructures<SpinColourSpinColourMatrix> structures
        = compute_bilinear_structures(values, find_gammas);
    __apply_fierz(measurements, structures, LoopStructure::FULLY_CONNECTED);

    for (int c = 0; c < NUM_SPIN_STRUCTURES; c++) {
        for (int d = 0; d < NUM_SPIN_STRUCTURES; d++) {
            LoopStructure loop = LoopStructure::FULLY_CONNECTED;
            ColorStructure color = ColorStructure::COLOR_MIXED;
            SpinStructure gammaA = (SpinStructure) c;
            SpinStructure gammaB = (SpinStructure) d;

            if (!measurements.is_measurement_stored(gammaA, gammaB)) {
                continue;
            }

            SpinColourSpinColourMatrix &val
                = measurements.get_measurement(loop, color, gammaA, gammaB);
            scsc_index_swap(val, 1, 3);
        }
    }
}


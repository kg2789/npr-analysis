#ifndef __EVALUATE_H__

#define __EVALUATE_H__
#include "operator.h"
#include "contraction.h"
#include "projector.h"

using Grid::SpinColourSpinColourMatrix;
using Grid::SpinColourMatrix;

// Quarks come in pairs (qbar q') in every diagram; if position points to one
// of qbar, q', then this function returns the index of the other
static inline int adjacent_position(int position) {
    return position ^ 0x1;
}


enum class LoopStructure {
    FULLY_CONNECTED,
    DISCONNECTED_LOOP,
    CONNECTED_LOOP,
    // Used to keep track of the number of values in this enum; do not add
    // values after this
    NUM_LOOP_STRUCTURES,
};
const int NUM_LOOP_STRUCTURES = (int) LoopStructure::NUM_LOOP_STRUCTURES;

static const LoopStructure all_loop_structures[] = {
    LoopStructure::FULLY_CONNECTED,
    LoopStructure::DISCONNECTED_LOOP,
    LoopStructure::CONNECTED_LOOP
};

inline std::ostream& operator<<(std::ostream& os, LoopStructure &loop) {
    switch (loop) {
        case LoopStructure::FULLY_CONNECTED:
            os << "FULLY_CONNECTED"; break;
        case LoopStructure::DISCONNECTED_LOOP:
            os << "DISCONNECTED_LOOP"; break;
        case LoopStructure::CONNECTED_LOOP:
            os << "CONNECTED_LOOP"; break;
        case LoopStructure::NUM_LOOP_STRUCTURES:
            os << "NUM_LOOP_STRUCTURES"; break;
    }
    return os;
}

// In order to reduce the amount of data we have at any given point, we don't
// want to store both diagrams with left/right spin structures and those with
// vector/axial spin structures at the same time. To get around this, we only
// include one type of data in FourQOpMeasurements below and use a mode to
// identify which one is stored
enum class FourQDataMode {
    VECTOR_AXIAL,
    LEFT_RIGHT,
};

// Contains the relevant measurment information calcluated on the lattice for
// every 4-quark operator we need; the templating is so that we can use this
// for both 2-quark and 4-quark external states, which produce different
// measurement types
template <typename T>
class FourQOpMeasurements {
    public:
    FourQDataMode mode;
    std::vector<T> values;

    // Since we don't actually store data on every possible spin structures
    static const int NUM_STORED_SPIN
        = sizeof(va_spin_structures) / sizeof(SpinStructure);
    static_assert(sizeof(va_spin_structures) == sizeof(lr_spin_structures),
            "Both these should be 2 * sizeof(SpinStructure)");
    static const int DATA_SIZE = NUM_LOOP_STRUCTURES * NUM_COLOR_STRUCTURES
        * NUM_STORED_SPIN * NUM_STORED_SPIN;

    FourQOpMeasurements() {
        values.resize(DATA_SIZE);
        for (T &val: values) {
            val = Zero();
        }
        mode = FourQDataMode::VECTOR_AXIAL;
    }

    bool is_measurement_stored(SpinStructure gammaA, SpinStructure gammaB) const {
        SpinStructure x, y;
        switch (mode) {
            case FourQDataMode::VECTOR_AXIAL:
                x = SpinStructure::VECTOR;
                y = SpinStructure::AXIAL;
                break;
            case FourQDataMode::LEFT_RIGHT:
                x = SpinStructure::LEFT;
                y = SpinStructure::RIGHT;
                break;
        }
        return (gammaA == x || gammaA == y)
            && (gammaB == x || gammaB == y);
    }

    // We store the data somewhat opqauely through an indexing scheme which
    // depends on 'mode'; all attempts to access a particular measurement
    // should use this function
    int get_index(LoopStructure loop, ColorStructure color,
            SpinStructure gammaA, SpinStructure gammaB) const {

        //const int num_loops = NUM_LOOP_STRUCTURES;
        const int num_color = NUM_COLOR_STRUCTURES;
        const int num_spin = NUM_STORED_SPIN;
        // a,b,c,d are effectively the indices to a multi-dimensional array
        int a = (int) loop;
        int b = (int) color;
        int c,d;

        auto va_index = [this] (SpinStructure gamma) {
            switch (gamma) {
                case SpinStructure::VECTOR: return 0;
                case SpinStructure::AXIAL: return 1;
                default: std::cerr << "Error: invalid gamma "
                         << gamma << " in mode "
                             << (mode == FourQDataMode::VECTOR_AXIAL ? "VA" : "LR")
                             << std::endl;
                         exit(1);
            }
        };

        auto lr_index = [this] (SpinStructure gamma) {
            switch (gamma) {
                case SpinStructure::LEFT: return 0;
                case SpinStructure::RIGHT: return 1;
                default: std::cerr << "Error: invalid gamma "
                         << gamma << " in mode "
                             << (mode == FourQDataMode::VECTOR_AXIAL ? "VA" : "LR")
                             << std::endl;
                         exit(1);
            }
        };

        switch (mode) {
            case FourQDataMode::VECTOR_AXIAL:
                c = va_index(gammaA);
                d = va_index(gammaB);
                break;
            case FourQDataMode::LEFT_RIGHT:
                c = lr_index(gammaA);
                d = lr_index(gammaB);
                break;
        }

        return a * num_spin * num_spin * num_color
            + b * num_spin * num_spin
            + c * num_spin + d;
    }

    T& get_measurement(LoopStructure loop,
            ColorStructure color,
            SpinStructure gammaA,
            SpinStructure gammaB) {
        assert(is_measurement_stored(gammaA, gammaB));
        int idx = get_index(loop, color, gammaA, gammaB);
        return values[idx];
    }

    const T& get_measurement_const(LoopStructure loop,
            ColorStructure color,
            SpinStructure gammaA,
            SpinStructure gammaB) const {
        assert(is_measurement_stored(gammaA, gammaB));
        int idx = get_index(loop, color, gammaA, gammaB);
        return values[idx];
    }

    void switch_modes(FourQDataMode new_mode) {
        if (new_mode == mode)
            return;

        FourQDataMode old_mode = mode;

        const SpinStructure V = SpinStructure::VECTOR;
        const SpinStructure A = SpinStructure::AXIAL;
        const SpinStructure L = SpinStructure::LEFT;
        const SpinStructure R = SpinStructure::RIGHT;
        std::vector<T> new_values;
        new_values.resize(values.size());

        std::vector<T>& va_values = (old_mode == FourQDataMode::VECTOR_AXIAL) ? values : new_values;
        std::vector<T>& lr_values = (old_mode == FourQDataMode::LEFT_RIGHT) ? values : new_values;

        for (LoopStructure loop: all_loop_structures) {
            for (ColorStructure color: all_color_structures) {
                // get_index relies on the mode for determining the index, so
                // we have to be careful with 'mode' here
                mode = FourQDataMode::VECTOR_AXIAL;
                T &VV = va_values[get_index(loop, color, V, V)];
                T &VA = va_values[get_index(loop, color, V, A)];
                T &AV = va_values[get_index(loop, color, A, V)];
                T &AA = va_values[get_index(loop, color, A, A)];

                mode = FourQDataMode::LEFT_RIGHT;
                T &LL = lr_values[get_index(loop, color, L, L)];
                T &LR = lr_values[get_index(loop, color, L, R)];
                T &RL = lr_values[get_index(loop, color, R, L)];
                T &RR = lr_values[get_index(loop, color, R, R)];

                switch (new_mode) {
                    case FourQDataMode::LEFT_RIGHT:
                        LL = VV + AA - VA - AV;
                        LR = VV - AA + VA - AV;
                        RL = VV - AA - VA + AV;
                        RR = VV + AA + VA + AV;
                        break;
                    case FourQDataMode::VECTOR_AXIAL:
                        VV = 0.25 * (LL + LR + RL + RR);
                        AA = 0.25 * (LL - LR - RL + RR);
                        VA = 0.25 * (-LL - LR + RL + RR);
                        AV = 0.25 * (-LL + LR - RL + RR);
                        break;
                }
            }
        }

        values = new_values;
        mode = new_mode;
    }

    void operator +=(const FourQOpMeasurements<T> &other) {
        assert(values.size() == other.values.size());
        for (int i = 0; i < values.size(); i++) {
            values[i] += other.values[i];
        }
    }

    void operator *=(const Complex &c) {
        for (int i = 0; i < values.size(); i++) {
            values[i] *= c;
        }
    }

    void operator -=(const FourQOpMeasurements<T> &other) {
        assert(values.size() == other.values.size());
        for (int i = 0; i < values.size(); i++) {
            values[i] -= other.values[i];
        }
    }

    FourQOpMeasurements<T> operator +(const FourQOpMeasurements<T> &other) const {
        FourQOpMeasurements ret = *this;
        ret += other;
        return ret;
    }
};

typedef FourQOpMeasurements<SpinColourSpinColourMatrix> FourQOpFourQExtMeasurements;
typedef FourQOpMeasurements<SpinColourMatrix> FourQOpTwoQExtMeasurements;


enum class TwoQuarkOp {
    SCALAR,
    DSLASH_LEFT,
    DSLASH_RIGHT,
};
const int NUM_TWO_QUARK_OPS = 3;

const TwoQuarkOp all_subtraction_ops[] = {
    TwoQuarkOp::SCALAR,
    TwoQuarkOp::DSLASH_LEFT,
    TwoQuarkOp::DSLASH_RIGHT
};

class TwoQOpTwoQExtMeasurements {
    public:
        SpinColourMatrix values[NUM_TWO_QUARK_OPS];
        // Contains the values for < sbar op \gamma^5 d >
        // TODO: possibly get a better name for this
        SpinColourMatrix psuedoscalar_values[NUM_TWO_QUARK_OPS];

        TwoQOpTwoQExtMeasurements() = default;

        SpinColourMatrix get_measurement(TwoQuarkOp &op) {
            return values[(int) op] + psuedoscalar_values[(int) op];
        }

        void operator +=(const TwoQOpTwoQExtMeasurements &other) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] += other.values[i];
                psuedoscalar_values[i] += other.psuedoscalar_values[i];
            }
        }

        void operator -=(const TwoQOpTwoQExtMeasurements &other) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] -= other.values[i];
                psuedoscalar_values[i] -= other.psuedoscalar_values[i];
            }
        }

        void operator *=(const Complex &c) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] = values[i] * c;
                psuedoscalar_values[i] = psuedoscalar_values[i] * c;
            }
        }

        TwoQOpTwoQExtMeasurements operator +(const TwoQOpTwoQExtMeasurements &other) const {
            TwoQOpTwoQExtMeasurements ret = *this;
            ret += other;
            return ret;
        }
};

class TwoQOpFourQExtMeasurements {
    public:
        SpinColourSpinColourMatrix values[NUM_TWO_QUARK_OPS];
        // Contains the values for < sbar op \gamma^5 d >
        SpinColourSpinColourMatrix psuedoscalar_values[NUM_TWO_QUARK_OPS];

        TwoQOpFourQExtMeasurements() = default;

        SpinColourSpinColourMatrix get_measurement(TwoQuarkOp &op) {
            return values[(int) op] + psuedoscalar_values[(int) op];
        }

        void operator +=(const TwoQOpFourQExtMeasurements &other) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] += other.values[i];
                psuedoscalar_values[i] += other.psuedoscalar_values[i];
            }
        }

        void operator -=(const TwoQOpFourQExtMeasurements &other) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] -= other.values[i];
                psuedoscalar_values[i] -= other.psuedoscalar_values[i];
            }
        }

        void operator *=(const Complex &c) {
            for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
                values[i] = values[i] * c;
                psuedoscalar_values[i] = psuedoscalar_values[i] * c;
            }
        }

        TwoQOpFourQExtMeasurements operator +(const TwoQOpFourQExtMeasurements &other) const {
            TwoQOpFourQExtMeasurements ret = *this;
            ret += other;
            return ret;
        }
};

class TrajectoryData;

SpinColourSpinColourMatrix evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<FourQuarkTerm, FourQuarkTerm> &matrix_element);
SpinColourSpinColourMatrix evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<FourQuarkTerm, BilinearTerm> &matrix_element,
        TwoQuarkOp op);
SpinColourMatrix evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<BilinearTerm, FourQuarkTerm> &matrix_element);
SpinColourMatrix evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<BilinearTerm, BilinearTerm> &matrix_element,
        TwoQuarkOp op);

Eigen::MatrixXcd compute_mixings(TrajectoryData &data,
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &projection_ops,
        std::vector<OperatorRepresentation> &reps);
#endif

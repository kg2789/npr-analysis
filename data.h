#ifndef __DATA_H__

#define __DATA_H__

#include <Grid/Grid.h>
#include <Grid/GridCore.h>

#include "evaluate.h"
#include "metadata.h"
#include "data/io.h"

class FourQuarkOpData {
    public:
        FourQOpFourQExtMeasurements fourq_ext;
        FourQOpTwoQExtMeasurements twoq_ext;

        void operator +=(const FourQuarkOpData &other) {
            twoq_ext += other.twoq_ext;
            fourq_ext += other.fourq_ext;
        }

        void operator -=(const FourQuarkOpData &other) {
            twoq_ext -= other.twoq_ext;
            fourq_ext -= other.fourq_ext;
        }

        void operator *=(const Complex &c) {
            twoq_ext *= c;
            fourq_ext *= c;
        }
};

class TwoQuarkOpData {
    public:
        TwoQOpTwoQExtMeasurements twoq_ext;
        TwoQOpFourQExtMeasurements fourq_ext;

        void operator +=(const TwoQuarkOpData &other) {
            twoq_ext += other.twoq_ext;
            fourq_ext += other.fourq_ext;
        }

        void operator -=(const TwoQuarkOpData &other) {
            twoq_ext -= other.twoq_ext;
            fourq_ext -= other.fourq_ext;
        }

        void operator *=(const Complex &c) {
            twoq_ext *= c;
            fourq_ext *= c;
        }
};

class ExternalLegData {
    public:
        SpinColourMatrix Sin_average, Sout_average;

        void operator +=(const ExternalLegData &other) {
            Sin_average += other.Sin_average;
            Sout_average += other.Sout_average;
        }

        void operator -=(const ExternalLegData &other) {
            Sin_average -= other.Sin_average;
            Sout_average -= other.Sout_average;
        }

        void operator *=(const Complex &c) {
            Sin_average *= c;
            Sout_average *= c;
        }
};


class G1Data : Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(G1Data,
                SpinColourMatrix, twoq_scalar,
                SpinColourMatrix, twoq_gamma5,
                SpinColourSpinColourMatrix, fourq_scalar,
                SpinColourSpinColourMatrix, fourq_gamma5);

        void operator +=(const G1Data &other) {
            twoq_scalar  += other.twoq_scalar;
            twoq_gamma5  += other.twoq_gamma5;
            fourq_scalar += other.fourq_scalar;
            fourq_gamma5 += other.fourq_gamma5;
        }

        void operator -=(const G1Data &other) {
            twoq_scalar  -= other.twoq_scalar;
            twoq_gamma5  -= other.twoq_gamma5;
            fourq_scalar -= other.fourq_scalar;
            fourq_gamma5 -= other.fourq_gamma5;
        }

        void operator *=(const Complex &c) {
            twoq_scalar  *= c;
            twoq_gamma5  *= c;
            fourq_scalar *= c;
            fourq_gamma5 *= c;
        }
};

class CurrentData : Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CurrentData,
                std::vector<SpinColourMatrix>, vector,
                std::vector<SpinColourMatrix>, axial);

        void operator +=(const CurrentData &other) {
            assert(vector.size() == other.vector.size());
            assert(axial.size() == other.axial.size());
            assert(vector.size() == axial.size());

            int dims = vector.size();
            for (int mu = 0; mu < dims; mu++) {
                vector[mu] += other.vector[mu];
                axial[mu] += other.axial[mu];
            }
        }

        void operator -=(const CurrentData &other) {
            assert(vector.size() == other.vector.size());
            assert(axial.size() == other.axial.size());
            assert(vector.size() == axial.size());

            int dims = vector.size();
            for (int mu = 0; mu < dims; mu++) {
                vector[mu] -= other.vector[mu];
                axial[mu] -= other.axial[mu];
            }
        }

        void operator *=(const Complex &c) {
            assert(vector.size() == axial.size());

            int dims = vector.size();
            for (int mu = 0; mu < dims; mu++) {
                vector[mu] *= c;
                axial[mu] *= c;
            }
        }
};

// Stores the results from exactly one trajectory
class TrajectoryData {
    public:
        int trajectory;
        FourQuarkOpData fourq_op;
        TwoQuarkOpData twoq_op;
        ExternalLegData external_legs;
        G1Data g1;
        CurrentData current;

        TrajectoryData() = default;

        void operator +=(const TrajectoryData& other) {
            fourq_op += other.fourq_op;
            twoq_op += other.twoq_op;
            external_legs += other.external_legs;
            g1 += other.g1;
            current += other.current;
        }

        void operator -=(const TrajectoryData& other) {
            fourq_op -= other.fourq_op;
            twoq_op -= other.twoq_op;
            external_legs -= other.external_legs;
            g1 -= other.g1;
            current -= other.current;
        }

        void operator *=(const Complex &c) {
            fourq_op *= c;
            twoq_op *= c;
            external_legs *= c;
            g1 *= c;
            current *= c;
        }


        // Most of the time when we convert the data we don't need the
        // Left/Right spin structures (they can be derived from the VV, VA, AV,
        // AA structures), so we usually don't include these structures. For
        // acutual calculations, however, we require these structures, so this
        // method will calculate and fill in the left/right structures
        void switch_modes(FourQDataMode mode) {
            fourq_op.fourq_ext.switch_modes(mode);
            fourq_op.twoq_ext.switch_modes(mode);
        }
};

class RunData {
    public:
        Metadata metadata;
        std::vector<TrajectoryData> trajectory_data;
};

RunData read_data(const std::string &dir, Format format = Format::RABBOTT, bool full = true);
void write_data(const std::string &dir, RunData &data, Format format = Format::RABBOTT);

template <typename T> RealD l2_norm(const T&);

template <typename T> RealD l2_difference(const T& a, const T& b) {
    // Some of these types are too large to allocate on the stack, so we avoid
    // stack allocation when possible
    T *tmp = new T();

    *tmp = a;
    *tmp -= b;
    RealD ret = l2_norm(*tmp);

    delete tmp;

    return ret;
}

template <typename T, typename F>
RealD filtered_l2_difference(const T& a, const T& b, F filter) {
    RealD diff = 0.0;
    assert(a.mode == b.mode);
    for (LoopStructure loop : all_loop_structures) {
        for (ColorStructure color : all_color_structures) {
            for (SpinStructure gammaA : all_spin_structures) {
                for (SpinStructure gammaB : all_spin_structures) {
                    if (!a.is_measurement_stored(gammaA, gammaB)
                     || !b.is_measurement_stored(gammaA, gammaB)) {
                        continue;
                    }
                    const auto &matA
                        = a.get_measurement_const(loop, color, gammaA, gammaB);
                    const auto &matB
                        = b.get_measurement_const(loop, color, gammaA, gammaB);
                    if (filter(loop, color, gammaA, gammaB))
                        diff += l2_difference(matA, matB);
                }
            }
        }
    }
    return diff;
}

// Computes the difference present in diagrams with the given loop structure.
// Useful when debugging to check where differences come from
template <typename T>
RealD l2_diagram_difference(const T& a, const T& b, LoopStructure loop) {
    auto filter
        = [loop](LoopStructure l, ColorStructure, SpinStructure, SpinStructure) {
        return l == loop;
    };
    return filtered_l2_difference(a, b, filter);
}

template <typename T> RealD relative_difference(const T& a, const T& b) {
    RealD average_norm = (l2_norm(a) + l2_norm(b)) / 2.0;
    return l2_difference(a, b) / average_norm;
}

#endif


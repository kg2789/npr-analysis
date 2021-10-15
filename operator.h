#ifndef __OPERATOR__DEFS__

#define __OPERATOR__DEFS__

#include <Grid/Grid.h>

#include <unordered_map>
#include <iostream>
#include <iterator>

using Grid::RealD;

enum class Quark: int {
    U = 1,
    D = 2,
    S = 3,
    C = 4,
    UBAR = -1,
    DBAR = -2,
    SBAR = -3,
    CBAR = -4,
};

inline Quark antiparticle(const Quark q) {
    return (Quark)(-(int)q);
}

inline std::ostream& operator<< (std::ostream &os, const Quark &q) {
    Quark tmp;
    if ((int)q < 0) {
        tmp = antiparticle(q);
    }
    else {
        tmp = q;
    }
    switch (tmp) {
        case Quark::U: os << "u"; break;
        case Quark::D: os << "d"; break;
        case Quark::S: os << "s"; break;
        case Quark::C: os << "c"; break;
        default: os << "Unkown Quark " << (int)q;
    }
    if ((int)q < 0) {
        os << "bar";
    }
    return os;
}

inline bool can_contract(const Quark &q1, const Quark &q2) {
    return q1 == antiparticle(q2);
}

template <typename T>
class LinearCombination {
    public:
        std::unordered_map<T, RealD> term_coefficients;

        LinearCombination(const T &term) {
            term_coefficients.insert({term, 1});
        }

        LinearCombination() {}

        void operator+= (const LinearCombination<T> &other) {
            LinearCombination<T> ret = other;

            for (auto &it: other.term_coefficients) {
                const T &term = it.first;
                const RealD &coeff = it.second;

                auto entry = term_coefficients.find(term);
                if (entry != ret.term_coefficients.end()) {
                    entry->second += coeff;
                }
                else {
                    term_coefficients.insert({ term, coeff });
                }
            }
        }

        LinearCombination<T> operator- () const {
            LinearCombination<T> ret = *this;
            for (auto &it: ret.term_coefficients) {
                it.second = -it.second;
            }
            return ret;
        }


        // TODO: This could avoid a copy (see += above)
        void operator -= (const LinearCombination<T> &other) {
            *this += -other;
        }

        LinearCombination<T> operator-(const LinearCombination<T> &other) {
            LinearCombination<T> ret = -other;
            ret += *this;
            return ret;
        }

        LinearCombination<T> operator+ (const LinearCombination<T> &other) const {
            LinearCombination<T> ret = *this;
            ret += other;
            return ret;
        }

        void operator*= (const RealD &coeff) {
            for (auto &it: term_coefficients) {
                it.second *= coeff;
            }

        }

        LinearCombination<T> operator* (const RealD &coeff) const {
            LinearCombination<T> ret = *this;
            ret *= coeff;
            return ret;
        }

        T sum() {
            T ret;
            for (auto &it: term_coefficients) {
                ret += it.first * it.second;
            }
            return ret;
        }

        friend std::ostream& operator<< (std::ostream& os, const LinearCombination<T> &lc) {
            bool first = true;
            for (auto &it: lc.term_coefficients) {
                if (!first) {
                    os << " + ";
                }
                else {
                    first = false;
                }
                os << it.second << " * " << it.first;
            }
            return os;
        }
};

// Allows left multiplication
template <typename T>
inline LinearCombination<T> operator*(const RealD coeff, const LinearCombination<T> &lc) {
    return lc * coeff;
}

enum class SpinStructure {
    LEFT,
    RIGHT,
    VECTOR, // \psibar \gamma^\mu \psi
    AXIAL, // \psibar \gamma^\mu \gamma^5 \psi
    UNCONTRACTED, // Used to mark external states
};
// Note: we don't include UNCONTRACTED since it's a special case we usually
// don't want to include
const int NUM_SPIN_STRUCTURES = 4;

const SpinStructure all_spin_structures[] = {
    SpinStructure::LEFT,
    SpinStructure::RIGHT,
    SpinStructure::VECTOR,
    SpinStructure::AXIAL
};

// Contains all spin structures with only vector or axial structures
const SpinStructure va_spin_structures[] = {
    SpinStructure::VECTOR,
    SpinStructure::AXIAL
};

const SpinStructure lr_spin_structures[] = {
    SpinStructure::LEFT,
    SpinStructure::RIGHT
};

inline std::ostream& operator<< (std::ostream &os, const SpinStructure &s) {
    switch (s) {
        case SpinStructure::LEFT: os << "L"; break;
        case SpinStructure::RIGHT: os << "R"; break;
        case SpinStructure::VECTOR: os << "V"; break;
        case SpinStructure::AXIAL: os << "A"; break;
        case SpinStructure::UNCONTRACTED: break;
        default: os << "Unknown spin structure";
    }
    return os;
}

class BilinearTerm {
    public:
        Quark q1, q2;
        SpinStructure spin_structure;

        BilinearTerm(Quark q1_, Quark q2_, SpinStructure structure)
            : q1(q1_), q2(q2_), spin_structure(structure) {}

        friend std::ostream& operator << (std::ostream& os, const BilinearTerm &term) {
            os << term.q1 << " " << term.spin_structure << " " << term.q2;
            return os;
        }

        bool operator== (const BilinearTerm &other) const {
            return q1 == other.q1
                && q2 == other.q2
                && spin_structure == other.spin_structure;
        }

        const int num_quarks() const { return 2; }

        Quark get_quark(int i) const {
            switch(i) {
                case 0: return q1;
                case 1: return q2;
            }
            assert(0);
        }
};

using BilinearOperator = LinearCombination<BilinearTerm>;

enum class ColorStructure {
    COLOR_MIXED,
    COLOR_DIAG,
    UNCONTRACTED, // used to mark external states
};
// As with NUM_SPIN_STRUCTURES above, we don't include UNCONTRACTED since it is
// a special case
const int NUM_COLOR_STRUCTURES = 2;

const ColorStructure all_color_structures[] = {
    ColorStructure::COLOR_MIXED,
    ColorStructure::COLOR_DIAG
};

inline std::ostream& operator<< (std::ostream &os, const ColorStructure &s) {
    switch (s) {
        case ColorStructure::COLOR_DIAG: os << "COLOR_DIAG"; break;
        case ColorStructure::COLOR_MIXED: os << "COLOR_MIXED"; break;
        case ColorStructure::UNCONTRACTED: os << "UNCONTRACTED"; break;
        default: os << "UNKOWN_COLOR_STRUCTURE"; break;
    }
    return os;
}

class FourQuarkTerm {
    public:
    BilinearTerm bilinear1, bilinear2;
    ColorStructure color_structure;

    FourQuarkTerm(const BilinearTerm &b1, const BilinearTerm &b2, const ColorStructure &s)
        : bilinear1(b1), bilinear2(b2), color_structure(s) {}

    friend std::ostream& operator << (std::ostream& os, const FourQuarkTerm &term) {
        os << "(" << term.bilinear1 << ") (" << term.bilinear2 << ")";
        return os;
    }

    bool operator== (const FourQuarkTerm &other) const {
        return bilinear1 == other.bilinear1
            && bilinear2 == other.bilinear2
            && color_structure == other.color_structure;
    }

    const int num_quarks() const { return 4; }

    Quark get_quark(int i) const {
        switch(i) {
            case 0: return bilinear1.q1;
            case 1: return bilinear1.q2;
            case 2: return bilinear2.q1;
            case 3: return bilinear2.q2;
        }
        assert(0);
    }
};

using FourQuarkOperator = LinearCombination<FourQuarkTerm>;

// Hashing functions needed for hashing for unordered_map/LinearCombination
// Primarily taken from here: 
// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
namespace std {
    template <>
        struct hash<BilinearTerm> {
            std::size_t operator()(const BilinearTerm& b) const {
                int x = (int)b.q1 + 100 * (int)b.q2 + 10000 * (int)b.spin_structure;
                return std::hash<int>()(x);
            }
        };

    template <>
        struct hash<FourQuarkTerm> {
            std::size_t operator()(const FourQuarkTerm& t) const {
                using std::size_t;

                size_t res = 17;
                res = res * 31 + hash<BilinearTerm>()(t.bilinear1);
                res = res * 31 + hash<BilinearTerm>()(t.bilinear2);
                res = res * 31 + hash<int>()((int)t.color_structure);
                return res;
            }
        };
}


inline FourQuarkTerm tensor_prod(const BilinearTerm &b1, 
        const BilinearTerm &b2, const ColorStructure &structure) {
    return FourQuarkTerm(b1, b2, structure);
}

FourQuarkOperator tensor_prod(const BilinearOperator &b1, 
        const BilinearOperator &b2, const ColorStructure structure);

enum class OperatorBasis {
    GREG,
    MASAAKI,
    MASAAKI_QSLASH
};

std::vector<FourQuarkOperator> get_operators(OperatorBasis basis);

enum class ProjectorBasis;
std::vector<FourQuarkOperator> get_external_states(ProjectorBasis basis);

enum class TwoQuarkOp;
std::vector<TwoQuarkOp> get_subtraction_operators();

FourQuarkTerm make_external_state(Quark q1, Quark q2, Quark q3, Quark q4);


#endif

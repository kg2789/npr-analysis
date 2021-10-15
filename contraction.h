#ifndef __CONTRACTION_H__

#define __CONTRACTION_H__

#include "operator.h"

class Contraction {
    public:
        std::vector<int> contractions;

        Contraction(int num_quarks) {
            contractions.resize(num_quarks);
            for (int &x: contractions) {
                x = -1;
            }
        }

        void add_contraction(int i, int j) {
            assert(i != j);
            contractions[i] = j;
            contractions[j] = i;
        }

        void remove_contraction(int i) {
            int j = contractions[i];
            assert (j >= 0);
            assert (contractions[j] == i);
            contractions[i] = -1;
            contractions[j] = -1;

        }

        // Given an index i, finds the index that i has been contracted to
        int get_contracted_pair(int i) const {
            assert(0 <= i && i < contractions.size());
            return contractions[i];
        }

        friend std::ostream& operator<< (std::ostream& os, const Contraction &c) {
            std::vector<bool> printed;
            printed.reserve(c.contractions.size());
            for (int i = 0; i < c.contractions.size(); i++) {
                printed.push_back(false);
            }

            for (int i = 0; i < c.contractions.size(); i++) {
                int j = c.contractions[i];
                if (printed[i] && printed[j])
                    continue;
                // We never should have printed half a contraction, which is
                // the only way printed[i] and printed[j] can differ
                assert(!printed[i] && !printed[j]);

                printed[i] = true;
                printed[j] = true;

                os << i << "-" << j;

                bool done = true;
                for (bool p: printed) {
                    if (!p) {
                        // We still have more to print
                        os << " ";
                        done = false;
                        break;
                    }
                }
                if (done) break;
            }
            return os;
        }
};

// A diagram consisting of an external state (Text) along some internal
// operator (Tint) with an 4-quark operator. For instance, a
// Diagram<FourQuarkTerm, FourQuarkTerm> this might represent a term like
// <(s dbar) (u ubar) x (sbar L d) (ubar L u)>
//
// The quarks are numbered from left to right in order to facilitate contractions
template <typename Text, typename Tint>
class Diagram {
    public:
        Text external_state;
        Tint internal_op;

        Diagram(Text _external_state, Tint _internal_op)
            : external_state(_external_state), internal_op(_internal_op) {}

        bool is_external(int i) const {
            return i < external_state.num_quarks();
        }

        bool is_internal(int i) const {
            return !is_external(i);
        }

        const int num_quarks() const {
            return external_state.num_quarks() + internal_op.num_quarks();
        }

        Quark get_quark(int i) const {
            assert(0 <= i && i < num_quarks());
            if (is_external(i)) {
                return external_state.get_quark(i);
            }
            else {
                int idx = i - external_state.num_quarks();
                return internal_op.get_quark(idx);
            }
        }

        bool operator== (const Diagram<Text, Tint>& other) const {
            return external_state == other.external_state
                && internal_op == other.internal_op;
        }

        friend std::ostream& operator<< (std::ostream& os, const Diagram<Text, Tint>& t) {
            os << "< " << t.external_state << " x " << t.internal_op << " >";
            return os;
        }
};


// Hashing functions needed for hashing for unordered_map/LinearCombination
// Primarily taken from here: 
// https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
namespace std {
    template <typename Text, typename Tint>
        struct hash< Diagram<Text, Tint> > {
            std::size_t operator()(const Diagram<Text, Tint> &t) const {
                std::size_t res = 17;
                res = res * 31 + hash<Text>()(t.external_state);
                res = res * 31 + hash<Tint>()(t.internal_op);
                return res;
            }
        };
}

template <typename Text, typename Tint>
class MatrixElement: public LinearCombination< Diagram<Text, Tint> > {
    public:
        MatrixElement(const Diagram<Text, Tint> &d)
            : LinearCombination< Diagram<Text, Tint> >(d) {}
        MatrixElement() {}
};

template <typename Text, typename Tint>
inline MatrixElement<Text, Tint> make_matrix_element(const LinearCombination<Text>& external_state,
        const LinearCombination<Tint> &internal_op) {
    MatrixElement<Text, Tint> ret;
    for (auto &i: external_state.term_coefficients) {
        const Text &external_term = i.first;
        const RealD &c1 = i.second;
        for (auto &j: internal_op.term_coefficients) {
            const Tint &internal_term = j.first;
            const RealD &c2 = j.second;

            MatrixElement<Text, Tint> new_term
                = Diagram<Text, Tint>(external_term, internal_term);
            RealD coeff = c1 * c2;
            ret += coeff * new_term;
        }
    }
    return ret;
}

// In some cases we may have external states where the quarks don't allow any
// contractions. By (ab)using the way we store quarks as integers we note that
// the sum of two contracted quarks is 0, so if the sum of all quarks is
// nonzero then there cannot possibly be a contraction. The converse is not
// necessarily true though since there can be some miraculous cancellation, so
// this serves as a quick check which removes most pathological cases
template <typename Text, typename Tint>
static bool can_possibly_contract(const Diagram<Text, Tint> &diagram) {
    int sum = 0;
    for (int i = 0; i < diagram.num_quarks(); i++) {
        sum += (int) diagram.get_quark(i);
    }
    return sum == 0;
}

template <typename Text, typename Tint>
static void get_contractions_recursive(const Diagram<Text, Tint> &diagram,
        Contraction &current_contraction,
        std::vector<Contraction> &contractions,
        bool *contracted) {
    /*std::cout << "contracted: ";
    for (int i = 0; i < DIAGRAM_NUM_QUARKS; i++) {
        if (i != 0) std::cout << " ";
        std::cout << contracted[i];
    }
    std::cout << std::endl;*/
    int first_uncontracted = -1;
    for (int i = 0; i < diagram.num_quarks(); i++) {
        if (!contracted[i]) {
            first_uncontracted = i;
            break;
        }
    }
    if (first_uncontracted < 0) {
        contractions.push_back(current_contraction);
        return;
    }

    Quark uncontracted_quark = diagram.get_quark(first_uncontracted);
    contracted[first_uncontracted] = true;
    for (int i = first_uncontracted + 1; i < diagram.num_quarks(); i++) {
        if (contracted[i]) continue;

        Quark q = diagram.get_quark(i);
        if (!can_contract(uncontracted_quark, q)) continue;

        contracted[i] = true;
        current_contraction.add_contraction(first_uncontracted, i);

        get_contractions_recursive(diagram,
                current_contraction, contractions, contracted);

        contracted[i] = false;
        current_contraction.remove_contraction(first_uncontracted);
    }
    contracted[first_uncontracted] = false;
}


template <typename Text, typename Tint>
std::vector<Contraction> get_contractions(const Diagram<Text, Tint> &diagram) {
    std::vector<Contraction> ret;

    if (!can_possibly_contract(diagram))
        return ret;

    bool *contracted = new bool[diagram.num_quarks()];
    for (int i = 0; i < diagram.num_quarks(); i++)
        contracted[i] = false;
    Contraction current_contraction(diagram.num_quarks());
    get_contractions_recursive(diagram, current_contraction, ret, contracted);

    delete[] contracted;

    return ret;
}

#endif


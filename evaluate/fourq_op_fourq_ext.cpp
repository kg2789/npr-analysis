#include "../operator.h"
#include "../contraction.h"
#include "../evaluate.h"
#include "../util.h"
#include "../data.h"


// All of the relevant structure we need out of a diagram for evaluation
// purposes
struct FourQInfo {
    ColorStructure color_structure;
    // Indicates the toplogical structure of the diagram
    LoopStructure loop_structure;
    // These both indicate insertions of left/right-handed spin structures,
    // with gammaA denoting the spin structure in the internal (sbar d) pair
    // and gammaB denoting the other insertion
    SpinStructure gammaA, gammaB;
    // Every internal 4-quark op has the structure (sbar d); these tell us
    // which external operators the sbar and d are contracted to
    int external_s, external_dbar;
    // Fermionic sign; always +/- 1
    int sign;
};

static std::vector<std::vector<int>> get_loops(const Contraction &contraction,
        int num_quarks) {
    std::vector<std::vector<int>> loops;
    bool *visited = new bool[num_quarks];
    for (int i = 0; i < num_quarks; i++)
        visited[i] = false;

    while(true) {
        int unvisited = -1;
        for (int i = 0; i < num_quarks; i++) {
            if (!visited[i]) {
                unvisited = i;
                break;
            }
        }
        if (unvisited < 0)
            break;

        std::vector<int> loop;
        int curr = unvisited;
        do {
            visited[curr] = true;
            loop.push_back(curr);
            curr = adjacent_position(curr);

            visited[curr] = true;
            loop.push_back(curr);
            curr = contraction.get_contracted_pair(curr);
        } while (curr != unvisited);
        loops.push_back(loop);
    }

    delete[] visited;
    return loops;
}

// Finds the loop which contains the quark at index idx
static std::vector<int>& find_loop(std::vector<std::vector<int>> &loops,
        int idx) {
    for (int i = 0; i < loops.size(); i++) {
        for (const int &index: loops[i]) {
            if (index == idx) {
                return loops[i];
            }
        }
    }
    assert(0);
}

static LoopStructure
classify_loop_structure_four_q_op(std::vector<std::vector<int>> &loops,
        int internal_sbar, int internal_d, const Contraction &contraction,
        const Diagram<FourQuarkTerm, FourQuarkTerm> &diagram) {
    // internal_q is part of the (q qbar) pair that we always have in addition
    // to the (sbar d) pair
    int internal_q;
    switch (internal_d) {
        case 5: internal_q = 7; break;
        case 7: internal_q = 5; break;
        default:
                // internal_d should always be 5 or 7
                std::cerr << "Error: internal_d = "
                    << internal_d << std::endl;
                assert(0);
    }

    const std::vector<int> &sd_path = find_loop(loops, internal_sbar);
    const std::vector<int> &internal_loop = find_loop(loops, internal_q);

    if (internal_loop.size() == 2) {
        return LoopStructure::DISCONNECTED_LOOP;
    }
    if (&internal_loop == &sd_path) {
        int qbar_contracted = contraction.get_contracted_pair(internal_q);
        if (diagram.is_internal(qbar_contracted)) {
            return LoopStructure::CONNECTED_LOOP;
        }
    }

    return LoopStructure::FULLY_CONNECTED;
}

static FourQInfo classify(const Diagram<FourQuarkTerm, FourQuarkTerm> &diagram,
        const Contraction &contraction) {
    const FourQuarkTerm &op = diagram.internal_op;
    const int num_quarks = diagram.num_quarks();

    // internal_sbar and internal_d are the quarks to which the external s and
    // dbar quarks are contracted. Usually these should be the internal (sbar d)
    // pair, but if the other pair is (sbar s) or (dbar d) then we might have
    // to rearrange
    int internal_sbar = 4;
    if (diagram.is_internal(contraction.get_contracted_pair(internal_sbar))) {
        internal_sbar = 6;
    }

    int internal_d = 5;
    if (diagram.is_internal(contraction.get_contracted_pair(internal_d))) {
        internal_d = 7;
    }

    FourQInfo ret;

    std::vector<std::vector<int>> loops = get_loops(contraction, num_quarks);
    ret.color_structure = op.color_structure;
    ret.loop_structure = classify_loop_structure_four_q_op(loops,
            internal_sbar, internal_d, contraction, diagram);
    ret.external_s = contraction.get_contracted_pair(internal_sbar);
    ret.external_dbar = contraction.get_contracted_pair(internal_d);
    if (internal_sbar == 4) {
        ret.gammaA = op.bilinear1.spin_structure;
        ret.gammaB = op.bilinear2.spin_structure;
    }
    else { // internal_sbar == 6
        ret.gammaA = op.bilinear2.spin_structure;
        ret.gammaB = op.bilinear1.spin_structure;
    }
    ret.sign = (loops.size() % 2 == 0) ? 1 : -1;

    return ret;
}

static SpinColourSpinColourMatrix evaluate(FourQOpFourQExtMeasurements &measurments,
        const FourQInfo &info) {
    SpinColourSpinColourMatrix ret
        = measurments.get_measurement(info.loop_structure,
                                      info.color_structure,
                                      info.gammaA,
                                      info.gammaB);

    if (info.external_s != 0) {
        scsc_index_swap(ret, 0, info.external_s);
    }
    if (info.external_dbar != 1) {
        scsc_index_swap(ret, 1, info.external_dbar);
    }
    if (info.sign != 1) {
        ret = ret * (double)info.sign;
    }

    return ret;
}

static SpinColourSpinColourMatrix
evaluate_diagram_with_contraction(FourQOpFourQExtMeasurements &measurments,
        const Diagram<FourQuarkTerm, FourQuarkTerm> &diagram,
        const Contraction &contraction) {
    FourQInfo info = classify(diagram, contraction);
    return evaluate(measurments, info);
}

static SpinColourSpinColourMatrix
evaluate_diagram(FourQOpFourQExtMeasurements &measurments,
        const Diagram<FourQuarkTerm, FourQuarkTerm> &diagram) {
    SpinColourSpinColourMatrix ret;

    std::vector<Contraction> contractions = get_contractions(diagram);
    for (const Contraction &c: contractions) {
        ret += evaluate_diagram_with_contraction(measurments, diagram, c);
    }
    return ret;
}


static SpinColourSpinColourMatrix
__evaluate_matrix_element(FourQOpFourQExtMeasurements &measurments,
        MatrixElement<FourQuarkTerm, FourQuarkTerm> &matrix_element) {
    SpinColourSpinColourMatrix ret;
    for (auto &p: matrix_element.term_coefficients) {
        RealD coeff = p.second;
        const Diagram<FourQuarkTerm, FourQuarkTerm>& diagram = p.first;

        ret += coeff * evaluate_diagram(measurments, diagram);
    }
    return ret;
}

SpinColourSpinColourMatrix
evaluate_matrix_element(TrajectoryData& data,
        MatrixElement<FourQuarkTerm, FourQuarkTerm> &matrix_element) {
    return __evaluate_matrix_element(data.fourq_op.fourq_ext, matrix_element);
}

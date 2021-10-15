#include "../evaluate.h"
#include "../data.h"
#include "../util.h"

struct DiagramInfo {
    // No loop structure: every diagram has the same structure
    // No color structure: every diagram is color diagonal
    int external_s, external_dbar;
    TwoQuarkOp op;
    int sign;
};

static DiagramInfo classify(const Diagram<FourQuarkTerm, BilinearTerm> &diagram,
        const Contraction& contraction, TwoQuarkOp op) {
    const int internal_sbar = 4;
    const int internal_d = 5;

    DiagramInfo ret;
    ret.external_s = contraction.get_contracted_pair(internal_sbar);
    ret.external_dbar = contraction.get_contracted_pair(internal_d);
    ret.op = op;
    ret.sign = ret.external_s == adjacent_position(ret.external_dbar) ? 1 : -1;

    return ret;
}

static SpinColourSpinColourMatrix
evaluate_info(TwoQOpFourQExtMeasurements &measurements, const DiagramInfo& info) {
    TwoQuarkOp op = info.op;
    SpinColourSpinColourMatrix ret
        = measurements.get_measurement(op);

    if (info.external_s != 0) {
        scsc_index_swap(ret, 0, info.external_s);
    }
    if (info.external_dbar != 1) {
        scsc_index_swap(ret, 1, info.external_dbar);
    }
    if (info.sign != 1) {
        ret = (double) info.sign * ret;
    }

    return ret;
}

//Everything below here is almost identical to the other cases; perhaps there
//is a way to combine these
static SpinColourSpinColourMatrix
evaluate_diagram_with_contraction(TwoQOpFourQExtMeasurements& measurements,
        const Diagram<FourQuarkTerm, BilinearTerm> &diagram,
        const Contraction& contraction, TwoQuarkOp op) {
   DiagramInfo info = classify(diagram, contraction, op);
   return evaluate_info(measurements, info);
}

static SpinColourSpinColourMatrix
evaluate_diagram(TwoQOpFourQExtMeasurements& measurements,
        const Diagram<FourQuarkTerm, BilinearTerm> &diagram, TwoQuarkOp op) {
    SpinColourSpinColourMatrix ret;

    std::vector<Contraction> contractions = get_contractions(diagram);
    for (const Contraction& c: contractions) {
        ret += evaluate_diagram_with_contraction(measurements,
                diagram, c, op);
    }
    return ret;
}

static SpinColourSpinColourMatrix
__evaluate_matrix_element(TwoQOpFourQExtMeasurements &measurements,
        MatrixElement<FourQuarkTerm, BilinearTerm> &matrix_element,
        TwoQuarkOp op) {
    SpinColourSpinColourMatrix ret;
    for (auto &p: matrix_element.term_coefficients) {
        RealD coeff = p.second;
        const Diagram<FourQuarkTerm, BilinearTerm>& diagram = p.first;

        ret += coeff * evaluate_diagram(measurements, diagram, op);
    }
    return ret;
}

SpinColourSpinColourMatrix
evaluate_matrix_element(TrajectoryData& data,
        MatrixElement<FourQuarkTerm, BilinearTerm> &matrix_element,
        TwoQuarkOp op) {
    return __evaluate_matrix_element(data.twoq_op.fourq_ext, matrix_element, op);
}

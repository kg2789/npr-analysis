#include "../evaluate.h"
#include "../data.h"

struct DiagramInfo {
    LoopStructure loop;
    ColorStructure color;
    SpinStructure gammaA, gammaB;
    int sign;
};

static DiagramInfo classify(const Diagram<BilinearTerm, FourQuarkTerm> &diagram,
        const Contraction &contraction) {
    const FourQuarkTerm &op = diagram.internal_op;

    const int external_s = 0;
    const int external_dbar = 1;

    const int internal_sbar = contraction.get_contracted_pair(external_s);
    const int internal_d = contraction.get_contracted_pair(external_dbar);

    DiagramInfo ret;
    ret.color = op.color_structure;
    if (internal_sbar == adjacent_position(internal_d)) {
        ret.loop = LoopStructure::DISCONNECTED_LOOP;
        ret.sign = 1;
    }
    else {
        ret.loop = LoopStructure::CONNECTED_LOOP;
        ret.sign = -1;
    }
    if (internal_sbar == 2) {
        ret.gammaA = op.bilinear1.spin_structure;
        ret.gammaB = op.bilinear2.spin_structure;
    }
    else {
        assert(internal_sbar == 4);
        ret.gammaA = op.bilinear2.spin_structure;
        ret.gammaB = op.bilinear1.spin_structure;
    }
    return ret;
}

static SpinColourMatrix evaluate(FourQOpTwoQExtMeasurements &measurements,
        const DiagramInfo &info) {
    SpinColourMatrix ret
        = measurements.get_measurement(info.loop, info.color,
                                       info.gammaA, info.gammaB);
    if (info.sign != 1) {
        ret = (double)info.sign * ret;
    }
    return ret;
}

static SpinColourMatrix
evaluate_diagram_with_contraction(FourQOpTwoQExtMeasurements &measurements,
        const Diagram<BilinearTerm, FourQuarkTerm> &diagram,
        const Contraction &contraction) {
    DiagramInfo info = classify(diagram, contraction);
    return evaluate(measurements, info);
}

static SpinColourMatrix
evaluate_diagram(FourQOpTwoQExtMeasurements &measurements,
        const Diagram<BilinearTerm, FourQuarkTerm> &diagram) {
    SpinColourMatrix ret;

    for (const Contraction &c: get_contractions(diagram)) {
        ret += evaluate_diagram_with_contraction(measurements, diagram, c);
    }
    return ret;
}

static SpinColourMatrix
__evaluate_matrix_element(FourQOpTwoQExtMeasurements &measurements,
        MatrixElement<BilinearTerm, FourQuarkTerm> &matrix_element) {
    SpinColourMatrix ret;
    for (auto &p: matrix_element.term_coefficients) {
        RealD &coeff = p.second;
        const Diagram<BilinearTerm, FourQuarkTerm> &diagram = p.first;

        ret += coeff * evaluate_diagram(measurements, diagram);
    }
    return ret;
}

SpinColourMatrix
evaluate_matrix_element(TrajectoryData &data,
        MatrixElement<BilinearTerm, FourQuarkTerm> &matrix_element) {
    return __evaluate_matrix_element(data.fourq_op.twoq_ext, matrix_element);
}

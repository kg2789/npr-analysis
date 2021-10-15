#include "amputate.h"
#include "util.h"


AmputationData::AmputationData(const SpinColourMatrix &prop1, const SpinColourMatrix &prop2) {
    prop1_inv = invert_sc(prop1);
    prop2_inv = invert_sc(prop2);

    // left = S(p_1)^{-1} \otimes S(p_1)^{-1}
    tensor_prod(left, prop1_inv, prop1_inv);
    // right = S(p_2)^{-1} \otimes S(p_2)^{-1}
    tensor_prod(right, prop2_inv, prop2_inv);
}

template <>
void amputate<SpinColourMatrix>(SpinColourMatrix &mat,
        const AmputationData &amp) {
    mat = amp.prop1_inv * mat * amp.prop2_inv;
}

// TODO: check that this matches conventions (seems to match Greg's thesis eq. 7.26)
template <>
void amputate<SpinColourSpinColourMatrix>(SpinColourSpinColourMatrix &mat,
        const AmputationData &amp) {
    mat = amp.left * mat * amp.right;
}

template <typename T> void fourq_op_amputate(FourQOpMeasurements<T>& measurements,
        const AmputationData &amp) {
    for (auto &val: measurements.values) {
        amputate(val, amp);
    }
}

template <>
void amputate<FourQOpFourQExtMeasurements>(FourQOpFourQExtMeasurements &m,
        const AmputationData &amp) {
    fourq_op_amputate<SpinColourSpinColourMatrix>(m, amp);
}

template <>
void amputate<FourQOpTwoQExtMeasurements>(FourQOpTwoQExtMeasurements &m,
        const AmputationData &amp) {
    fourq_op_amputate<SpinColourMatrix>(m, amp);
}

template <typename T> void twoq_op_amputate(T& measurements,
        const AmputationData &amp) {
    for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
        amputate(measurements.values[i], amp);
        amputate(measurements.psuedoscalar_values[i], amp);
    }
}

template <>
void amputate<TwoQOpFourQExtMeasurements>(TwoQOpFourQExtMeasurements &m,
        const AmputationData &amp) {
    twoq_op_amputate(m, amp);
}

template <>
void amputate<TwoQOpTwoQExtMeasurements>(TwoQOpTwoQExtMeasurements &m,
        const AmputationData &amp) {
    twoq_op_amputate(m, amp);
}

template <>
void amputate<FourQuarkOpData>(FourQuarkOpData &data,
        const AmputationData &amp) {
    amputate(data.fourq_ext, amp);
    amputate(data.twoq_ext, amp);
}

template <>
void amputate<TwoQuarkOpData>(TwoQuarkOpData &data,
        const AmputationData &amp) {
    amputate(data.fourq_ext, amp);
    amputate(data.twoq_ext, amp);
}

template <>
void amputate<G1Data>(G1Data &data,
        const AmputationData &amp) {
    amputate(data.twoq_scalar, amp);
    amputate(data.twoq_gamma5, amp);
    amputate(data.fourq_scalar, amp);
    amputate(data.fourq_gamma5, amp);
}

template <>
void amputate<CurrentData>(CurrentData &data,
        const AmputationData &amp) {
    for (SpinColourMatrix& val : data.vector) {
        amputate(val, amp);
    }
    for (SpinColourMatrix& val : data.axial) {
        amputate(val, amp);
    }
}

template <>
void amputate<TrajectoryData>(TrajectoryData &data,
        const AmputationData &amp) {
    amputate(data.fourq_op, amp);
    amputate(data.twoq_op, amp);
    amputate(data.g1, amp);
    amputate(data.current, amp);
}

void amputate_trajectory_data(TrajectoryData &data) {
    SpinColourMatrix &Sin = data.external_legs.Sin_average;
    SpinColourMatrix &Sout = data.external_legs.Sout_average;
    amputate(data, Sin, Sout);
}

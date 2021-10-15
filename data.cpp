#include "data.h"
#include "metadata.h"
#include "data/io.h"

template <>
RealD l2_norm<SpinColourMatrix>(const SpinColourMatrix &mat) {
    return norm2(mat);
}

template <>
RealD l2_norm<SpinColourSpinColourMatrix>(const SpinColourSpinColourMatrix &mat) {
    return norm2(mat);
}

template <typename T> RealD fourq_op_norm(const FourQOpMeasurements<T>& measurements) {
    RealD ret = 0.0;
    for (const T& val: measurements.values) {
        ret += l2_norm(val);
    }
    return ret;
}

template <>
RealD l2_norm<FourQOpFourQExtMeasurements>(const FourQOpFourQExtMeasurements &m) {
    return fourq_op_norm<SpinColourSpinColourMatrix>(m);
}

template <>
RealD l2_norm<FourQOpTwoQExtMeasurements>(const FourQOpTwoQExtMeasurements &m) {
    return fourq_op_norm<SpinColourMatrix>(m);
}

template <typename T> RealD twoq_op_norm(const T& measurements) {
    RealD ret = 0.0;
    for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
        ret += l2_norm(measurements.values[i]);
        ret += l2_norm(measurements.psuedoscalar_values[i]);
    }
    return ret;
}

template <>
RealD l2_norm<TwoQOpTwoQExtMeasurements>(const TwoQOpTwoQExtMeasurements &m) {
    return twoq_op_norm(m);
}

template <>
RealD l2_norm<TwoQOpFourQExtMeasurements>(const TwoQOpFourQExtMeasurements &m) {
    return twoq_op_norm(m);
}

template <>
RealD l2_norm<FourQuarkOpData>(const FourQuarkOpData &data) {
    return l2_norm(data.fourq_ext) + l2_norm(data.twoq_ext);
}

template <>
RealD l2_norm<TwoQuarkOpData>(const TwoQuarkOpData &data) {
    return l2_norm(data.fourq_ext) + l2_norm(data.twoq_ext);
}

template <>
RealD l2_norm<ExternalLegData>(const ExternalLegData &data) {
    return l2_norm(data.Sin_average) + l2_norm(data.Sout_average);
}

template <>
RealD l2_norm<G1Data>(const G1Data &data) {
    return l2_norm(data.twoq_scalar) + l2_norm(data.twoq_gamma5)
        + l2_norm(data.fourq_scalar) + l2_norm(data.fourq_gamma5);
}

template <>
RealD l2_norm<TrajectoryData>(const TrajectoryData &data) {
    return l2_norm(data.fourq_op) + l2_norm(data.twoq_op)
        + l2_norm(data.external_legs) + l2_norm(data.g1);
}

RunData read_data(const std::string &dir, Format format, bool full) {
    RunData ret;
    ret.metadata = read_metadata(dir + "/metadata.xml");

    ret.trajectory_data.reserve(ret.metadata.trajectories.size());
    TrajectoryData *tmp = new TrajectoryData();
    for (int traj : ret.metadata.trajectories) {
        // All read and write operations are done with vector/axial structures
        // in mind; if we don't set this here then read_trajectory might throw
        // an error on an invalid mode
        tmp->fourq_op.fourq_ext.mode = FourQDataMode::VECTOR_AXIAL;
        tmp->fourq_op.twoq_ext.mode = FourQDataMode::VECTOR_AXIAL;

        read_trajectory(dir, tmp, ret.metadata, traj, format);
        if (full) {
            // For moving around data we don't need to compute the left/right
            // spin-structure diagrams, but we do need this for analysis. This
            // fills in the missing data when we need it, but by setting 'full'
            // to false we can avoid the extra work.
            tmp->switch_modes(FourQDataMode::LEFT_RIGHT);
        }
        ret.trajectory_data.push_back(*tmp);
    }

    return ret;
}

void write_data(const std::string &dir, RunData &data, Format format) {
    make_dir(dir);
    write_metadata(dir + "/metadata.xml", data.metadata);

    int i = 0;
    for (int traj : data.metadata.trajectories) {
        TrajectoryData *tmp = &data.trajectory_data[i++];
        if (tmp->trajectory != traj) {
            std::cerr << "Error: trajectory data does not match metadata"
                << "(metadata has trajectory " << traj << "while trajectory data "
                << "has " << tmp->trajectory << ")" << std::endl;
        }
        write_trajectory(dir, tmp, data.metadata, traj, format);
    }
}

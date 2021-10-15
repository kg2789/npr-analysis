#include "greg/WilsonMatrix.h"
#include "greg/DoubleWilsonMatrix.h"
#include "greg/convert.hpp"
#include "../data.h"
//#include "util.h"

#include <fstream>
#include <array>

/* This file was taken from Greg's code (ConvertZiyuanData.cpp) and modified to
 * match up with the rest of my code */

static WilsonMatrix ReadZiyuanExternalLeg(const char* filename)
{
    WilsonMatrix ret;
    std::ifstream f(filename);
    if (!f) {
	printf("Couldn't open %s\n", filename);
	exit(-1);
    }
    for (int s1 = 0; s1 < 4; s1++) {
	for (int c1 = 0; c1 < 3; c1++) {
	    for (int s2 = 0; s2 < 4; s2++) {
		for (int c2 = 0; c2 < 3; c2++) {
		    double real;
		    double imag;
		    f >> real;
		    f >> imag;
		    if (f.fail()) {
			printf("Failure reading %s\n", filename);
			exit(-1);
		    }
		    ret(s1, c1, s2, c2) = complex<double>(real, imag);
		}
	    }
	}
    }
    f.close();
    return ret;
}

static DoubleWilsonMatrix ReadZiyuanDoubleWilsonMatrixFromStream(std::ifstream &f)
{
    DoubleWilsonMatrix ret;
    for (int i = 0; i < 144 * 144; ++i) {
	double real, imag;
	f >> real;
	f >> imag;
	// Qi/Ziyuan/Daiqian have the tensor indices in a different order than I do
	int qi_indices[4];
	int j = i;
	qi_indices[0] = j % 12; j /= 12;
	qi_indices[1] = j % 12; j /= 12;
	qi_indices[2] = j % 12; j /= 12;
	qi_indices[3] = j;
	int my_index =
	    qi_indices[3] +
	    12 * (qi_indices[2] +
	    12 * (qi_indices[1] +
	    12 * qi_indices[0]));

	ret.data[my_index] = complex<double>(real, imag);
    }

    return ret;
}

static WilsonMatrix ReadZiyuanOpLineFromStream(std::ifstream &f)
{
    WilsonMatrix ret;
    for (int i = 0; i < 12 * 12; ++i) {
	double real, imag;
	f >> real;
	f >> imag;
	// Qi/Ziyuan/Daiqian have the tensor indices in a different order than I do
	int my_index = 12 * (i % 12) + (i / 12);
	ret.data[my_index] = complex<double>(real, imag);
    }

    return ret;
}

enum VAStructure {
    GammaVV,
    GammaAA,
    GammaVA,
    GammaAV,
};

static const std::array<ColorStructure, 4> zbai_color_structure_order{{
    ColorStructure::COLOR_DIAG,
    ColorStructure::COLOR_MIXED,
    ColorStructure::COLOR_DIAG,
    ColorStructure::COLOR_MIXED
}};

// order in which VA structures appear in EVEN PARITY Ziyuan format bulk files
static const std::array<VAStructure, 4> zbai_even_VA_order{{
    GammaVV, GammaVV, GammaAA, GammaAA
}};
// order in which VA structures appear in ODD PARITY Ziyuan format bulk files
static const std::array<VAStructure, 4> zbai_odd_VA_order{{
    GammaAV, GammaAV, GammaVA, GammaVA
}};

static SpinStructure get_gammaA(VAStructure va) {
    const SpinStructure V = SpinStructure::VECTOR;
    const SpinStructure A = SpinStructure::AXIAL;
    switch (va) {
        case GammaVV: return V;
        case GammaAA: return A;
        case GammaVA: return V;
        case GammaAV: return A;
    }
}

static SpinStructure get_gammaB(VAStructure va) {
    const SpinStructure V = SpinStructure::VECTOR;
    const SpinStructure A = SpinStructure::AXIAL;
    switch (va) {
        case GammaVV: return V;
        case GammaAA: return A;
        case GammaVA: return A;
        case GammaAV: return V;
    }
}

static void ReadZiyuanFourqBulkFile(TrajectoryData &data,
        const char* filename, Parity parity)
{
    std::ifstream f(filename);
    if (!f) {
	printf("couldn't open bulk file %s\n", filename);
	exit(-1);
    }

    const std::array<VAStructure, 4> &VA_order
        = (parity == Parity::EVEN ? zbai_even_VA_order : zbai_odd_VA_order);

    // read fully connected diagrams
    LoopStructure loop = LoopStructure::FULLY_CONNECTED;
    for (int i = 0; i < 4; ++i) {
        DoubleWilsonMatrix measurement
            = ReadZiyuanDoubleWilsonMatrixFromStream(f);

        SpinStructure gammaA = get_gammaA(VA_order[i]);
        SpinStructure gammaB = get_gammaB(VA_order[i]);
        ColorStructure color = zbai_color_structure_order[i];

        SpinColourSpinColourMatrix &out
            = data.fourq_op.fourq_ext.get_measurement(loop, color,
                    gammaA, gammaB);
        out = double_wilson_to_scsc(measurement);
    }

    // read disconnected diagrams
    WilsonMatrix spectator = ReadZiyuanOpLineFromStream(f);
    std::vector<LoopStructure> disconnected_loop_structures = {
        LoopStructure::DISCONNECTED_LOOP,
        LoopStructure::CONNECTED_LOOP
    };
    for (LoopStructure loop: disconnected_loop_structures) {
        for (int i = 0; i < 4; ++i) {
            WilsonMatrix op_line = ReadZiyuanOpLineFromStream(f);
            DoubleWilsonMatrix measurement
                = DoubleWilsonMatrix::OuterProduct(op_line, spectator);

            SpinStructure gammaA = get_gammaA(VA_order[i]);
            SpinStructure gammaB = get_gammaB(VA_order[i]);
            ColorStructure color = zbai_color_structure_order[i];

            SpinColourSpinColourMatrix &out
                = data.fourq_op.fourq_ext.get_measurement(loop, color,
                        gammaA, gammaB);
            out = double_wilson_to_scsc(measurement);
        }
    }

    // Read subtraction operators

    SpinColourSpinColourMatrix (&values)[NUM_TWO_QUARK_OPS]
        = (parity == Parity::EVEN ? data.twoq_op.fourq_ext.values
                : data.twoq_op.fourq_ext.psuedoscalar_values);
    std::vector<SpinColourSpinColourMatrix> amplitudes;
    amplitudes.reserve(NUM_TWO_QUARK_OPS);
    for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
        WilsonMatrix op_line = ReadZiyuanOpLineFromStream(f);
        DoubleWilsonMatrix measurement
            = DoubleWilsonMatrix::OuterProduct(op_line, spectator);
        amplitudes.push_back(double_wilson_to_scsc(measurement));
    }
    values[(int) TwoQuarkOp::SCALAR] = amplitudes[0];
    values[(int) TwoQuarkOp::DSLASH_RIGHT] = amplitudes[1];
    values[(int) TwoQuarkOp::DSLASH_LEFT] = amplitudes[2];

    f.close();
}

static void ReadZiyuanTwoqBulkFile(TrajectoryData &data,
        const char* filename, Parity parity)
{
    std::ifstream f(filename);
    if (!f) {
	printf("couldn't open bulk file %s\n", filename);
	exit(-1);
    }

    const std::array<VAStructure, 4> &VA_order
        = (parity == Parity::EVEN ? zbai_even_VA_order : zbai_odd_VA_order);

    // read disconnected diagrams
    std::vector<LoopStructure> disconnected_loop_structures = {
        LoopStructure::DISCONNECTED_LOOP,
        LoopStructure::CONNECTED_LOOP
    };
    for (LoopStructure loop: disconnected_loop_structures) {
        for (int i = 0; i < 4; ++i) {
            WilsonMatrix measurement = ReadZiyuanOpLineFromStream(f);

            SpinStructure gammaA = get_gammaA(VA_order[i]);
            SpinStructure gammaB = get_gammaB(VA_order[i]);
            ColorStructure color = zbai_color_structure_order[i];

            SpinColourMatrix &out
                = data.fourq_op.twoq_ext.get_measurement(loop, color,
                        gammaA, gammaB);
            out = wilson_to_sc(measurement);
        }
    }

    // Read subtraction operators

    SpinColourMatrix (&values)[NUM_TWO_QUARK_OPS]
        = (parity == Parity::EVEN ? data.twoq_op.twoq_ext.values
                : data.twoq_op.twoq_ext.psuedoscalar_values);
    std::vector<SpinColourMatrix> amplitudes;
    amplitudes.reserve(NUM_TWO_QUARK_OPS);
    for (int i = 0; i < NUM_TWO_QUARK_OPS; i++) {
        WilsonMatrix measurement = ReadZiyuanOpLineFromStream(f);
        amplitudes.push_back(wilson_to_sc(measurement));
    }

    values[(int) TwoQuarkOp::SCALAR] = amplitudes[0];
    values[(int) TwoQuarkOp::DSLASH_RIGHT] = amplitudes[1];
    values[(int) TwoQuarkOp::DSLASH_LEFT] = amplitudes[2];

    f.close();
}

void ConvertZiyuanNPRData(TrajectoryData &data, const char* zbai_dir,
        int conf, const int mom1[4], const int mom2[4])
{
    // convert external legs

    char zbai_prop1_filename[512];
    char zbai_prop2_filename[512];
    sprintf(zbai_prop1_filename, "%s/traj_%d_legs_p0.txt", zbai_dir, conf);
    sprintf(zbai_prop2_filename, "%s/traj_%d_legs_p1.txt", zbai_dir, conf);
    WilsonMatrix zbai_prop1 = ReadZiyuanExternalLeg(zbai_prop1_filename);
    WilsonMatrix zbai_prop2 = ReadZiyuanExternalLeg(zbai_prop2_filename);

    data.external_legs.Sout_average = wilson_to_sc(zbai_prop1);
    data.external_legs.Sin_average = wilson_to_sc(zbai_prop2);

    // convert fourq diagrams

    char even_fourq_bulk_file[512];
    char odd_fourq_bulk_file[512];
    sprintf(even_fourq_bulk_file, "%s/traj_%d_k2pipiNPR_even_p0_p1.txt", zbai_dir, conf);
    sprintf(odd_fourq_bulk_file, "%s/traj_%d_k2pipiNPR_odd_p0_p1.txt", zbai_dir, conf);

    ReadZiyuanFourqBulkFile(data, even_fourq_bulk_file, Parity::EVEN);
    ReadZiyuanFourqBulkFile(data, odd_fourq_bulk_file, Parity::ODD);

    // convert twoq diagrams

    char even_twoq_bulk_file[512];
    char odd_twoq_bulk_file[512];
    sprintf(even_twoq_bulk_file, "%s/traj_%d_subtractionCoeff_even_p0_p1.txt", zbai_dir, conf);
    sprintf(odd_twoq_bulk_file, "%s/traj_%d_subtractionCoeff_odd_p0_p1.txt", zbai_dir, conf);

    ReadZiyuanTwoqBulkFile(data, even_twoq_bulk_file, Parity::EVEN);
    ReadZiyuanTwoqBulkFile(data, odd_twoq_bulk_file, Parity::ODD);
}

static void read_Zq(const std::string &dir, TrajectoryData &data, int trajectory) {
    std::string filename = dir + "/traj_" + std::to_string(trajectory) + "_bilinear_p0_p1.txt";
    std::ifstream f(filename);

    WilsonMatrix scalar = ReadZiyuanOpLineFromStream(f);
    WilsonMatrix psuedoscalar = ReadZiyuanOpLineFromStream(f);

    data.current.vector.resize(4);
    data.current.axial.resize(4);
    for (int mu = 0; mu < 4; mu++) {
        WilsonMatrix mat = ReadZiyuanOpLineFromStream(f);
        data.current.vector[mu] = wilson_to_sc(mat);
    }
    for (int mu = 0; mu < 4; mu++) {
        WilsonMatrix mat = ReadZiyuanOpLineFromStream(f);
        // Looking at the results for Zq, it seems we need a minus sign here;
        // this sounds reasonable, since the axial current can be either
        // \gamma^\mu \gammma^5 or \gamma^\mu \gamma^5, and the difference
        // between the two is a minus sign.
        data.current.axial[mu] = -wilson_to_sc(mat);
    }

    f.close();
}

void read_trajectory_ziyuan(const std::string &dir, TrajectoryData *data,
        Metadata &metadata, int trajectory) {
    int mom1[4];
    int mom2[4];
    for (int i = 0; i < 4; i++) {
        mom1[i] = metadata.p1[i];
        mom2[i] = metadata.p2[i];
    }
    ConvertZiyuanNPRData(*data, dir.c_str(), trajectory, mom1, mom2);

    read_Zq(dir, *data, trajectory);
}

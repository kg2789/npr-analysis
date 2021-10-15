#ifndef __PROJECTOR_H__

#define __PROJECTOR_H__

#include <Grid/Grid.h>

#include "util.h"
#include "parity.h"
#include "metadata.h"

using Grid::SpinColourSpinColourMatrix;
using Grid::SpinColourMatrix;
using Grid::ComplexD;

class FourQuarkProjector {
    public:
        SpinColourSpinColourMatrix mat;

        FourQuarkProjector() {}
        FourQuarkProjector(const SpinColourSpinColourMatrix &mat_)
            : mat(mat_) {}
        FourQuarkProjector(const SpinStructure &a, const SpinStructure &b,
                ColorStructure structure);

        ComplexD project(SpinColourSpinColourMatrix &op);
};

class TwoQuarkProjector {
    public:
        SpinColourMatrix mat;

        TwoQuarkProjector() = default;
        TwoQuarkProjector(const SpinColourMatrix &mat_)
            : mat(mat_) {}

        ComplexD project(SpinColourMatrix &op);
};

enum class ProjectorBasis {
    GREG,
    MASAAKI,
    GREG_QSLASH,
    MASAAKI_QSLASH
};

enum class OperatorRepresentation {
	REP84_1,
    REP27_1,
    REP20_1,
    REP15_15,
    REP15_1,
    REP8_8,
    REP8_1,
    NONE
};

std::vector<OperatorRepresentation> get_representations(ProjectorBasis basis);

std::vector<FourQuarkProjector> get_projectors(ProjectorBasis basis, Parity parity,
        Metadata metadata);
// Returns the projectors which we use to determine the subtraction
// coefficients
std::vector<TwoQuarkProjector> get_subtraction_projectors(Parity parity,
        Metadata metadata);

#endif

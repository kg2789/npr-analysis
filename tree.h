#ifndef __TREE_H__

#define __TREE_H__

#include "contraction.h"
#include "util.h"
#include "projector.h"
#include "parity.h"
#include "data.h"

SpinColourSpinColourMatrix
evaluate_matrix_elem_tree_level(MatrixElement<FourQuarkTerm, FourQuarkTerm> &matrix_element);

TrajectoryData *tree_level_measurements(Parity parity);

Eigen::MatrixXcd tree_level_mixings(
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &projection_ops,
        std::vector<OperatorRepresentation> &reps,
        Parity parity);

#endif

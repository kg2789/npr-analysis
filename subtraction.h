#ifndef __SUBTRACTION_H__

#define __SUBTRACTION_H__

#include <Grid/Grid.h>
#include "evaluate.h"
#include "projector.h"

// Assumes data has already been amputated
Eigen::MatrixXcd subtraction_coefficients(TrajectoryData *data,
        std::vector<FourQuarkOperator> &opeartors,
        std::vector<TwoQuarkOp> &subtraction_operators,
        std::vector<TwoQuarkProjector> &projectors);

// Assumes data has already been amputated
Eigen::MatrixXcd subtraction_correction(TrajectoryData *data,
        std::vector<FourQuarkOperator> &operators,
        std::vector<FourQuarkOperator> &external_states,
        std::vector<FourQuarkProjector> &fourq_projectors,
        std::vector<OperatorRepresentation> &fourq_op_reps,
        std::vector<TwoQuarkOp> &subtraction_ops,
        std::vector<TwoQuarkProjector> &twoq_projectors,
        std::vector<OperatorRepresentation> &subtraction_op_reps);

std::vector<OperatorRepresentation> get_subtraction_reps(ProjectorBasis basis);

#endif

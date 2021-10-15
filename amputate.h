#ifndef __AMPUTATE_H__

#define __AMPUTATE_H__

#include "data.h"

// Contains the relevant inverses of matrices to amputate any diagram; serves
// as a cache so we don't recompute the inverse for every diagram
class AmputationData {
    public:
        // Inverses needed to amputate bilinears
        SpinColourMatrix prop1_inv, prop2_inv;
        // Inverses needed to amputate 4-quark diagrsm
        SpinColourSpinColourMatrix left, right;

        AmputationData(const SpinColourMatrix &prop1, const SpinColourMatrix &prop2);
};

template <typename T>
void amputate(T&, const AmputationData&);

template <typename T>
void amputate(T& vals, SpinColourMatrix &Sin, SpinColourMatrix &Sout) {
    AmputationData amp = AmputationData(Sout, Sin);
    return amputate(vals, amp);
}

void amputate_trajectory_data(TrajectoryData &data);

#endif

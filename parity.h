#ifndef __PARITY_H__

#define __PARITY_H__

#include <Grid/Grid.h>
#include <iostream>

using Grid::SpinColourMatrix;
using Grid::SpinColourSpinColourMatrix;

enum class Parity {
    EVEN,
    ODD,
    BOTH,
};

inline std::ostream& operator<< (std::ostream& os, const Parity& p) {
    switch(p) {
        case Parity::EVEN: os << "EVEN"; break;
        case Parity::ODD:  os << "ODD";  break;
        case Parity::BOTH: os << "BOTH"; break;
    }
    return os;
}

void parity_tranform(SpinColourSpinColourMatrix &mat);
void parity_tranform(SpinColourMatrix &mat);
SpinColourMatrix parity_project(const SpinColourMatrix &mat, Parity parity);
SpinColourSpinColourMatrix
parity_project(const SpinColourSpinColourMatrix &mat, Parity parity);

class RunData;
RunData read_data_with_parity(const std::string &dir, Parity parity);

#endif

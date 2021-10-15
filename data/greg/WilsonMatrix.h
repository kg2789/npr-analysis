#ifndef WILSON_MATRIX_H
#define WILSON_MATRIX_H

#include <vector>
#include <complex>
//#include "Matrix.h"
//#include "SpinMatrix.h"

class WilsonMatrix
{
public:
    static const int DATA_SIZE = 12 * 12;
    std::vector<std::complex<double> > data;

    int index(int s1, int c1, int s2, int c2) const;

    WilsonMatrix() : data(DATA_SIZE, std::complex<double>(0, 0)) {}
    WilsonMatrix(const char* filename) : data(DATA_SIZE) { LoadFromFile(filename); }

    void Zero();

    std::complex<double>& operator()(int s1, int c1, int s2, int c2);
    const std::complex<double>& operator()(int s1, int c1, int s2, int c2) const;

    WilsonMatrix& operator+=(const WilsonMatrix &rhs);
    WilsonMatrix& operator-=(const WilsonMatrix &rhs);
    WilsonMatrix operator*(const WilsonMatrix &rhs) const;
    /*WilsonMatrix operator*(const SpinMatrix &rhs) const;*/
    WilsonMatrix& operator*=(std::complex<double> coeff);
    WilsonMatrix operator+(const WilsonMatrix &rhs) const;
    WilsonMatrix operator-(const WilsonMatrix &rhs) const;

    WilsonMatrix Inverse() const;

    double Norm2() const;

    std::complex<double> Trace() const;
    /*SpinMatrix ColorTrace() const;
    SpinMatrix ColorTraceOverThree() const;*/

    void WriteToFile(const char *filename) const;
    void LoadFromFile(const char *filename);
};

WilsonMatrix operator*(std::complex<double> coeff, const WilsonMatrix &mat);
/*WilsonMatrix operator*(const SpinMatrix &sm, const WilsonMatrix &wm);*/

WilsonMatrix operator-(const WilsonMatrix m);

WilsonMatrix Amputate(const WilsonMatrix &unamputated_diagram,
    const WilsonMatrix &prop1, const WilsonMatrix &prop2);


#endif

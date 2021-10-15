#ifndef DOUBLE_WILSON_MATRIX_H
#define DOUBLE_WILSON_MATRIX_H

#include "WilsonMatrix.h"
//#include "SpinMatrix.h"

#include <vector>
#include <complex>

/*enum ColorStructure
{
    COLOR_STRUCTURE_DIAGONAL,
    COLOR_STRUCTURE_MIXED
};*/

// Object with the tensor structure of the outer product
// of two Wilson matrices. It has 4 Dirac indices and 4 color
// indices for a total of (4*3)^4 = 20736 complex entries.
//
// A DoubleWilsonMatrix can represent the tensor structure of a four-quark operator
// or the result of a four-quark contraction.
class DoubleWilsonMatrix
{
public:
    static const int DATA_SIZE = 12 * 12 * 12 * 12;
    std::vector<std::complex<double> > data;

    int index(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4) const;


    DoubleWilsonMatrix() : data(DATA_SIZE, std::complex<double>(0, 0)) {}

    DoubleWilsonMatrix(const char* filename) : data(DATA_SIZE) { LoadFromFile(filename); }
    
    void Zero();

    std::complex<double>& operator()(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4);
    const std::complex<double>& operator()(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4) const;

    static DoubleWilsonMatrix OuterProduct(const WilsonMatrix &a, const WilsonMatrix &b);

    // Contract a Wilson matrix with the one set of spin/color indices
    // index = 0: result_abcd = sum_x (m_ax this_xbcd)
    // index = 1: result_abcd = sum_x (this_axcd m_xb)
    // index = 2: result_abcd = sum_x (m_cx this_abxd)
    // index = 3: result_abcd = sum_x (this_abcx m_xd)
    /*DoubleWilsonMatrix Contract(const WilsonMatrix &m, int index) const;

    static DoubleWilsonMatrix Construct(const SpinMatrix &spin_structure_1,
	const SpinMatrix &spin_structure_2,
	ColorStructure color_structure);*/

    void Swap02();
    void Swap13();

    void Swap01();
    void Swap23();

    std::complex<double> Project(DoubleWilsonMatrix dwm) const;

    DoubleWilsonMatrix operator+(const DoubleWilsonMatrix &rhs) const;
    DoubleWilsonMatrix operator-(const DoubleWilsonMatrix &rhs) const;
    DoubleWilsonMatrix operator-() const; // unary minus
    DoubleWilsonMatrix& operator+=(const DoubleWilsonMatrix &rhs);
    DoubleWilsonMatrix& operator-=(const DoubleWilsonMatrix &rhs);
    DoubleWilsonMatrix& operator*=(double rhs);
    DoubleWilsonMatrix& operator*=(std::complex<double> rhs);

    double norm() const;

    void WriteToFile(const char* filename) const;
    void LoadFromFile(const char *filename);
    void WriteToBinaryFileReversingEndianness(const char* filename) const;
    void LoadFromBinaryFileReversingEndianness(const char* filename);
};

DoubleWilsonMatrix operator*(std::complex<double> coeff, const DoubleWilsonMatrix &dwm);

/*DoubleWilsonMatrix Amputate(const DoubleWilsonMatrix &unamputated, const WilsonMatrix &prop1, const WilsonMatrix &prop2,
    const WilsonMatrix &prop3, const WilsonMatrix &prop4);*/

#endif

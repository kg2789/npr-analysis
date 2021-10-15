#include "WilsonMatrix.h"
//#include "Matrix.h"

#include <fstream>
#include <cassert>

int WilsonMatrix::index(int s1, int c1, int s2, int c2) const
{
    assert(s1 >= 0 && c1 >= 0 && s2 >= 0 && c2 >= 0);
    assert(s1 < 4 && c1 < 3 && s2 < 4 && c2 < 3);
    int ind1 = 3 * s1 + c1;
    int ind2 = 3 * s2 + c2;
    return ind1 + 12 * ind2;
}

void WilsonMatrix::Zero()
{
    for (int i = 0; i < DATA_SIZE; ++i) data[i] = 0;
}

std::complex<double>& WilsonMatrix::operator()(int s1, int c1, int s2, int c2)
{
    return this->data[index(s1, c1, s2, c2)];
}

const std::complex<double>& WilsonMatrix::operator()(int s1, int c1, int s2, int c2) const
{
    return this->data[index(s1, c1, s2, c2)];
}


WilsonMatrix& WilsonMatrix::operator+=(const WilsonMatrix &rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] += rhs.data[i];
    }
    return *this;
}

WilsonMatrix& WilsonMatrix::operator-=(const WilsonMatrix &rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] -= rhs.data[i];
    }
    return *this;
}

WilsonMatrix WilsonMatrix::operator*(const WilsonMatrix &rhs) const
{
    WilsonMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int x = 0; x < 4; ++x) {
			for (int y = 0; y < 3; ++y) {
			    ret(s1, c1, s2, c2) += (*this)(s1, c1, x, y) * rhs(x, y, s2, c2);
			}
		    }
		}
	    }
	}
    }
    return ret;
}

/*WilsonMatrix WilsonMatrix::operator*(const SpinMatrix &rhs) const
{
    WilsonMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int x = 0; x < 4; ++x) {
			ret(s1, c1, s2, c2) += (*this)(s1, c1, x, c2) * rhs(x, s2);
		    }
		}
	    }
	}
    }
    return ret;
}*/

WilsonMatrix& WilsonMatrix::operator*=(std::complex<double> coeff)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] *= coeff;
    }
    return *this;
}

WilsonMatrix WilsonMatrix::operator+(const WilsonMatrix &rhs) const
{
    WilsonMatrix ret;
    for (int i = 0; i < DATA_SIZE; ++i) ret.data[i] = this->data[i] + rhs.data[i];
    return ret;
}

WilsonMatrix WilsonMatrix::operator-(const WilsonMatrix &rhs) const
{
    WilsonMatrix ret;
    for (int i = 0; i < DATA_SIZE; ++i) ret.data[i] = this->data[i] - rhs.data[i];
    return ret;
}

/*WilsonMatrix WilsonMatrix::Inverse() const
{
    Matrix<std::complex<double>, 12> mat;
    for (int s1 = 0; s1 < 4; s1++) {
	for (int c1 = 0; c1 < 3; c1++) {
	    for (int s2 = 0; s2 < 4; s2++) {
		for (int c2 = 0; c2 < 3; c2++) {
		    mat[3 * s1 + c1][3 * s2 + c2] = (*this)(s1, c1, s2, c2);
		}
	    }
	}
    }
    Matrix<std::complex<double>, 12> inv_mat = mat.Inverse();
    WilsonMatrix ret;
    for (int s1 = 0; s1 < 4; s1++) {
	for (int c1 = 0; c1 < 3; c1++) {
	    for (int s2 = 0; s2 < 4; s2++) {
		for (int c2 = 0; c2 < 3; c2++) {
		    ret(s1, c1, s2, c2) = inv_mat[3 * s1 + c1][3 * s2 + c2];
		}
	    }
	}
    }
    return ret;
}*/

double WilsonMatrix::Norm2() const
{
    double ret = 0.0;
    for (int i = 0; i < DATA_SIZE; ++i) ret += data[i].real() * data[i].real() + data[i].imag() * data[i].imag();
    return ret;
}

std::complex<double> WilsonMatrix::Trace() const
{
    std::complex<double> ret(0, 0);
    for (int s = 0; s < 4; ++s) {
	for (int c = 0; c < 3; ++c) {
	    ret += (*this)(s, c, s, c);
	}
    }
    return ret;
}

/*SpinMatrix WilsonMatrix::ColorTrace() const
{
    SpinMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int s2 = 0; s2 < 4; ++s2) {
	    for (int c = 0; c < 3; ++c) {
		ret(s1, s2) += (*this)(s1, c, s2, c);
	    }
	}
    }
    return ret;
}

SpinMatrix WilsonMatrix::ColorTraceOverThree() const
{
    return (1.0 / 3.0) * ColorTrace();
}*/

void WilsonMatrix::WriteToFile(const char *filename) const
{
    FILE * f = fopen(filename, "w");
    for (int s1 = 0; s1 < 4; s1++) {
	for (int c1 = 0; c1 < 3; c1++) {
	    for (int s2 = 0; s2 < 4; s2++) {
		for (int c2 = 0; c2 < 3; c2++) {
		    double re = (*this)(s1, c1, s2, c2).real();
		    double im = (*this)(s1, c1, s2, c2).imag();
		    fprintf(f, "%0.16e %0.16e ", re, im);
		}
	    }
	}
    }
    fprintf(f, "\n");
    fclose(f);
}

void WilsonMatrix::LoadFromFile(const char *filename)
{
    //printf("Going to load WilsonMatrix from %s\n", filename);
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
		    (*this)(s1, c1, s2, c2) = std::complex<double>(real, imag);
		}
	    }
	}
    }
    f.close();
}


WilsonMatrix operator*(std::complex<double> coeff, const WilsonMatrix &mat)
{
    WilsonMatrix ret;
    for (int i = 0; i < WilsonMatrix::DATA_SIZE; ++i) {
	ret.data[i] = coeff * mat.data[i];
    }
    return ret;
}

/*WilsonMatrix operator*(const SpinMatrix &sm, const WilsonMatrix &wm)
{
    WilsonMatrix ret;
    ret.Zero();
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int x = 0; x < 4; ++x) {
			ret(s1, c1, s2, c2) += sm(s1, x) * wm(x, c1, s2, c2);
		    }
		}
	    }
	}
    }
    return ret;
}*/

WilsonMatrix operator-(const WilsonMatrix m)
{
    WilsonMatrix ret;
    for (int i = 0; i < WilsonMatrix::DATA_SIZE; ++i) {
	ret.data[i] = -m.data[i];
    }
    return ret;
}


/*WilsonMatrix Amputate(const WilsonMatrix &unamputated_diagram,
    const WilsonMatrix &prop1, const WilsonMatrix &prop2)
{
    return prop1.Inverse() * unamputated_diagram * prop2.Inverse();
}*/

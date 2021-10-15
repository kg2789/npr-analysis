#include "DoubleWilsonMatrix.h"
#include "WilsonMatrix.h"
#include "util.h"

#include <fstream>
#include <cstdio>
#include <cassert>

int DoubleWilsonMatrix::index(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4) const
{
    assert(s1 >= 0 && c1 >= 0 && s2 >= 0 && c2 >= 0 && s3 >= 0 && c3 >= 0 && s4 >= 0 && c4 >= 0);
    assert(s1 < 4 && c1 < 3 && s2 < 4 && c2 < 3 && s3 < 4 && c3 < 3 && s4 < 4 && c4 < 3);
    int ind1 = 3 * s1 + c1;
    int ind2 = 3 * s2 + c2;
    int ind3 = 3 * s3 + c3;
    int ind4 = 3 * s4 + c4;
    return ind1 + 12 * (ind2 + 12 * (ind3 + 12 * ind4));
}

void DoubleWilsonMatrix::Zero()
{
    for (int i = 0; i < DATA_SIZE; ++i) data[i] = 0.0;
}

std::complex<double>& DoubleWilsonMatrix::operator()(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4)
{
    return data[index(s1, c1, s2, c2, s3, c3, s4, c4)];
}

const std::complex<double>& DoubleWilsonMatrix::operator()(int s1, int c1, int s2, int c2, int s3, int c3, int s4, int c4) const
{
    return data[index(s1, c1, s2, c2, s3, c3, s4, c4)];
}

DoubleWilsonMatrix DoubleWilsonMatrix::OuterProduct(const WilsonMatrix &a, const WilsonMatrix &b)
{
    DoubleWilsonMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    ret(s1, c1, s2, c2, s3, c3, s4, c4) = a(s1, c1, s2, c2) * b(s3, c3, s4, c4);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return ret;
}

// Contract a Wilson matrix with the one set of spin/color indices
// index = 0: result_abcd = sum_x (m_ax this_xbcd)
// index = 1: result_abcd = sum_x (this_axcd m_xb)
// index = 2: result_abcd = sum_x (m_cx this_abxd)
// index = 3: result_abcd = sum_x (this_abcx m_xd)
/*DoubleWilsonMatrix DoubleWilsonMatrix::Contract(const WilsonMatrix &m, int index) const
{
    DoubleWilsonMatrix ret;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {

				    for (int x = 0; x < 4; ++x) { // x and y are the contracted indices
					for (int y = 0; y < 3; ++y) {

					    switch (index) {
						case 0:
						    ret(s1, c1, s2, c2, s3, c3, s4, c4) += m(s1, c1, x, y) * (*this)(x, y, s2, c2, s3, c3, s4, c4);
						    break;

						case 1:
						    ret(s1, c1, s2, c2, s3, c3, s4, c4) += (*this)(s1, c1, x, y, s3, c3, s4, c4) * m(x, y, s2, c2);
						    break;

						case 2:
						    ret(s1, c1, s2, c2, s3, c3, s4, c4) += m(s3, c3, x, y) * (*this)(s1, c1, s2, c2, x, y, s4, c4);
						    break;

						case 3:
						    ret(s1, c1, s2, c2, s3, c3, s4, c4) += (*this)(s1, c1, s2, c2, s3, c3, x, y) * m(x, y, s4, c4);
						    break;

						case 4:
						    printf("invalid index: %d\n", index);
						    exit(-1);
						    break;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return ret;
}*/


/*DoubleWilsonMatrix DoubleWilsonMatrix::Construct(const SpinMatrix &spin_structure_1,
    const SpinMatrix &spin_structure_2,
    ColorStructure color_structure)
{
    DoubleWilsonMatrix ret; // initially all zeros
    for (int c1 = 0; c1 < 3; ++c1) {
	for (int c2 = 0; c2 < 3; ++c2) {
	    for (int c3 = 0; c3 < 3; ++c3) {
		for (int c4 = 0; c4 < 3; ++c4) {
		    if (color_structure == COLOR_STRUCTURE_DIAGONAL) {
			if (!(c1 == c2 && c3 == c4)) continue;
		    } else { // color-mixed
			if (!(c1 == c4 && c3 == c2)) continue;
		    }
		    for (int s1 = 0; s1 < 4; ++s1) {
			for (int s2 = 0; s2 < 4; ++s2) {
			    for (int s3 = 0; s3 < 4; ++s3) {
				for (int s4 = 0; s4 < 4; ++s4) {
				    ret(s1, c1, s2, c2, s3, c3, s4, c4) = spin_structure_1(s1, s2) * spin_structure_2(s3, s4);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return ret;
}*/

void DoubleWilsonMatrix::Swap02()
{
    DoubleWilsonMatrix swapped;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    swapped(s1, c1, s2, c2, s3, c3, s4, c4) = (*this)(s3, c3, s2, c2, s1, c1, s4, c4);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    *this = swapped;
}

void DoubleWilsonMatrix::Swap13()
{
    DoubleWilsonMatrix swapped;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    swapped(s1, c1, s2, c2, s3, c3, s4, c4) = (*this)(s1, c1, s4, c4, s3, c3, s2, c2);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    *this = swapped;
}

void DoubleWilsonMatrix::Swap01()
{
    DoubleWilsonMatrix swapped;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    swapped(s1, c1, s2, c2, s3, c3, s4, c4) = (*this)(s2, c2, s1, c1, s3, c3, s4, c4);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    *this = swapped;
}

void DoubleWilsonMatrix::Swap23()
{
    DoubleWilsonMatrix swapped;
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    swapped(s1, c1, s2, c2, s3, c3, s4, c4) = (*this)(s1, c1, s2, c2, s4, c4, s3, c3);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    *this = swapped;
}


std::complex<double> DoubleWilsonMatrix::Project(DoubleWilsonMatrix dwm) const
{
    std::complex<double> ret = 0.0;
    // WRONG for (int i = 0; i < DATA_SIZE; ++i) ret += this->data[i] * dwm.data[i];
    
    for (int s1 = 0; s1 < 4; ++s1) {
	for (int c1 = 0; c1 < 3; ++c1) {
	    for (int s2 = 0; s2 < 4; ++s2) {
		for (int c2 = 0; c2 < 3; ++c2) {
		    for (int s3 = 0; s3 < 4; ++s3) {
			for (int c3 = 0; c3 < 3; ++c3) {
			    for (int s4 = 0; s4 < 4; ++s4) {
				for (int c4 = 0; c4 < 3; ++c4) {
				    // Note the swapped indices on the projector
				    ret += (*this)(s2, c2, s1, c1, s4, c4, s3, c3) * dwm(s1, c1, s2, c2, s3, c3, s4, c4);
				}
			    }
			}
		    }
		}
	    }
	}
    }
    return ret;
}

DoubleWilsonMatrix DoubleWilsonMatrix::operator+(const DoubleWilsonMatrix &rhs) const
{
    DoubleWilsonMatrix ret;
    for (int i = 0; i < DATA_SIZE; ++i) {
	ret.data[i] = this->data[i] + rhs.data[i];
    }
    return ret;
}

DoubleWilsonMatrix DoubleWilsonMatrix::operator-(const DoubleWilsonMatrix &rhs) const
{
    DoubleWilsonMatrix ret;
    for (int i = 0; i < DATA_SIZE; ++i) {
	ret.data[i] = this->data[i] - rhs.data[i];
    }
    return ret;
}

// unary minus
DoubleWilsonMatrix DoubleWilsonMatrix::operator-() const
{
    DoubleWilsonMatrix ret;
    for (int i = 0; i < DATA_SIZE; ++i) {
	ret.data[i] = -this->data[i];
    }
    return ret;
}

DoubleWilsonMatrix& DoubleWilsonMatrix::operator+=(const DoubleWilsonMatrix &rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] += rhs.data[i];
    }
    return *this;
}

DoubleWilsonMatrix& DoubleWilsonMatrix::operator-=(const DoubleWilsonMatrix &rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] -= rhs.data[i];
    }
    return *this;
}

DoubleWilsonMatrix& DoubleWilsonMatrix::operator*=(double rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] *= rhs;
    }
    return *this;
}

DoubleWilsonMatrix& DoubleWilsonMatrix::operator*=(std::complex<double> rhs)
{
    for (int i = 0; i < DATA_SIZE; ++i) {
	this->data[i] *= rhs;
    }
    return *this;
}

DoubleWilsonMatrix operator*(std::complex<double> coeff, const DoubleWilsonMatrix &dwm)
{
    DoubleWilsonMatrix ret = dwm;
    ret *= coeff;
    return ret;
}

double DoubleWilsonMatrix::norm() const
{
    double ret = 0;
    for (int i = 0; i < DATA_SIZE; ++i) {
	std::complex<double> val = data[i];
	ret += val.real() * val.real() + val.imag() * val.imag();
    }
    return ret;
}

void DoubleWilsonMatrix::WriteToFile(const char* filename) const
{
    FILE* f = OpenFile(filename, "w");
    for (int i = 0; i < DATA_SIZE; ++i) {
	fprintf(f, "%0.16e %0.16e\n", data[i].real(), data[i].imag());
    }
    fclose(f);
}

void DoubleWilsonMatrix::LoadFromFile(const char* filename)
{
    //printf("Going to load DoubleWilsonMatrix from %s\n", filename);
    std::ifstream f(filename);
    if (!f) {
	printf("Couldn't open file %s\n", filename);
	exit(-1);
    }
    for (int i = 0; i < DATA_SIZE; ++i) {
	double real, imag;
	f >> real;
	f >> imag;
	data[i] = std::complex<double>(real, imag);
    }
    f.close();
}

void DoubleWilsonMatrix::WriteToBinaryFileReversingEndianness(const char* filename) const
{
    // usually we read files written on BG/Q, which writes doubles with reversed endianness
    // compared to qcdserver. So if we write out a file on qcdserver, we had better reverse
    // the endianness manually so that the result looks like a file written on BG/Q.
    std::vector<std::complex<double>> reversed_data(this->data);
    ReverseDoubleEndianness((double *)reversed_data.data(), 2 * DATA_SIZE);

    FILE* f = OpenFile(filename, "wb");
    fwrite(reversed_data.data(), sizeof(std::complex<double>), DATA_SIZE, f);
    fclose(f);
}

void DoubleWilsonMatrix::LoadFromBinaryFileReversingEndianness(const char* filename)
{
    FILE* f = OpenFile(filename, "rb");
    fread(data.data(), sizeof(std::complex<double>), DATA_SIZE, f);
    fclose(f);

    ReverseDoubleEndianness((double *)data.data(), 2 * DATA_SIZE); // account for difference between qcdserver and BG/Q
}

/*DoubleWilsonMatrix Amputate(const DoubleWilsonMatrix &unamputated, const WilsonMatrix &prop1, const WilsonMatrix &prop2,
    const WilsonMatrix &prop3, const WilsonMatrix &prop4)
{
    DoubleWilsonMatrix amp1 = unamputated.Contract(prop1.Inverse(), 0);
    DoubleWilsonMatrix amp2 = amp1.Contract(prop2.Inverse(), 1);
    DoubleWilsonMatrix amp3 = amp2.Contract(prop3.Inverse(), 2);
    DoubleWilsonMatrix ret = amp3.Contract(prop4.Inverse(), 3);
    return ret;
}*/

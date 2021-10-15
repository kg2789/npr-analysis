#ifndef __JACKKNIFE_H__

#define __JACKKNIFE_H__

#include <Grid/Grid.h>
#include <vector>
#include "output.h"

using Grid::RealD;
using Grid::iScalar;
using Grid::iVector;
using Grid::iMatrix;

// Squares/take the square-root of a vector/matrix element by element, which is
// necessary in calculating the jackknife errors
//template <typename T> strong_inline T elementwise_square(const T&);
//template <typename T> T elementwise_sqrt(const T&);

RealD elementwise_square(const RealD& val) {
    return val * val;
}

RealD elementwise_sqrt(const RealD& val) {
    return std::sqrt(val);
}

template <typename T>
std::complex<T> elementwise_square(const std::complex<T> &val) {
    std::complex<T> ret;
    ret.real(elementwise_square(val.real()));
    ret.imag(elementwise_square(val.imag()));
    return ret;
}

template <typename T>
std::complex<T> elementwise_sqrt(const std::complex<T> &val) {
    std::complex<T> ret;
    ret.real(elementwise_sqrt(val.real()));
    ret.imag(elementwise_sqrt(val.imag()));
    return ret;
}

using Eigen::MatrixXcd;

MatrixXcd elementwise_square(const MatrixXcd &val) {
    MatrixXcd ret = Eigen::MatrixXcd(val.rows(), val.cols());
    for (int i = 0; i < val.rows(); i++) {
        for (int j = 0; j < val.cols(); j++) {
            ret(i,j) = elementwise_square(val(i,j));
        }
    }
    return ret;
}

MatrixXcd elementwise_sqrt(const MatrixXcd &val) {
    MatrixXcd ret = Eigen::MatrixXcd(val.rows(), val.cols());
    for (int i = 0; i < val.rows(); i++) {
        for (int j = 0; j < val.cols(); j++) {
            ret(i,j) = elementwise_sqrt(val(i,j));
        }
    }
    return ret;
}


template <typename T, int N>
iMatrix<T, N> elementwise_square(const iMatrix<T, N> &val);
template <typename T, int N>
iVector<T, N> elementwise_square(const iVector<T, N> &val);

template <typename T, int N>
iMatrix<T, N> elementwise_sqrt(const iMatrix<T, N> &val);
template <typename T, int N>
iVector<T, N> elementwise_sqrt(const iVector<T, N> &val);

template <typename T>
iScalar<T> elementwise_square(const iScalar<T> &val) {
    iScalar<T> ret;
    ret._internal = elementwise_square(val._internal);
    return ret;
}

template <typename T>
iScalar<T> elementwise_sqrt(const iScalar<T> &val) {
    iScalar<T> ret;
    ret._internal = elementwise_sqrt(val._internal);
    return ret;
}

template <typename T, int N>
iVector<T, N> elementwise_square(const iVector<T, N> &val) {
    iVector<decltype(elementwise_square(val._internal[0])), N> ret;
    for (int i = 0; i < N; i++) {
        ret._internal[i] = elementwise_square(val._internal[i]);
    }
    return ret;
}

template <typename T, int N>
iVector<T, N> elementwise_sqrt(const iVector<T, N> &val) {
    iVector<decltype(elementwise_sqrt(val._internal[0])), N> ret;
    for (int i = 0; i < N; i++) {
        ret._internal[i] = elementwise_sqrt(val._internal[i]);
    }
    return ret;
}

template <typename T, int N>
iMatrix<T, N> elementwise_square(const iMatrix<T, N> &val) {
    iMatrix<T, N> ret;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ret._internal[i][j] = elementwise_square(val._internal[i][j]);
        }
    }
    return ret;
}

template <typename T, int N>
iMatrix<T, N> elementwise_sqrt(const iMatrix<T, N> &val) {
    iMatrix<T, N> ret;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ret._internal[i][j] = elementwise_sqrt(val._internal[i][j]);
        }
    }
    return ret;
}

template <typename T>
class JackknifeDatabase {
    public:
        std::vector<T> values;
        T *_central_value;

        JackknifeDatabase() {
            _central_value = new T();
        };
        JackknifeDatabase(const std::vector<T> &v): values(v) {
            _central_value = new T();
            *_central_value = compute_central_value();
        }
        JackknifeDatabase(const JackknifeDatabase &other): values(other.values) {
            _central_value = new T();
            *_central_value = *other._central_value;
        }

        void operator =(const JackknifeDatabase &other) {
            values = other.values;
            *_central_value = *other._central_value;
        }

        void operator +=(const JackknifeDatabase &other) {
            assert(values.size() == other.values.size());
            parallel_for(int i = 0; i < values.size(); i++) {
                values[i] += other.values[i];
            }
            *_central_value += *other._central_value;
        }

        void operator -=(const JackknifeDatabase &other) {
            assert(values.size() == other.values.size());
            parallel_for(int i = 0; i < values.size(); i++) {
                values[i] -= other.values[i];
            }
            *_central_value -= *other._central_value;
        }

        template <typename S>
        void operator *=(const JackknifeDatabase<S> &other) {
            assert(values.size() == other.values.size());
            parallel_for(int i = 0; i < values.size(); i++) {
                values[i] *= other.values[i];
            }
            *_central_value *= *other._central_value;
        }

        JackknifeDatabase operator +(const JackknifeDatabase &other) const {
            JackknifeDatabase ret = *this;
            ret += other;
            return ret;
        }

        JackknifeDatabase operator -(const JackknifeDatabase &other) const {
            JackknifeDatabase ret = *this;
            ret -= other;
            return ret;
        }

        template <typename F>
        auto fmap(F f) -> JackknifeDatabase<decltype(f(values[0]))> {
            auto ret = JackknifeDatabase<decltype(f(values[0]))>();
            ret.values.resize(values.size());
            parallel_for(int i = 0; i < values.size(); i++) {
                ret.values[i] = f(values[i]);
            }
            *ret._central_value = f(*_central_value);
            return ret;
        }

        // Used for applying functions which return void; if we try to use the
        // above implementation then the compiler complains about taking a void
        // reference
        template <typename F>
        void fmap_void(F f) {
            parallel_for(int i = 0; i < values.size(); i++) {
                f(values[i]);
            }
            f(*_central_value);
        }

        T compute_central_value() {
            T ret = values[0];
            for (const T& val: values) {
                if (&val == &values[0])
                    continue;
                ret += val;
            }
            ret *= 1.0 / values.size();
            return ret;
        }

        T &central_value(bool reaverage = false) {
            if (reaverage) {
                *_central_value = compute_central_value();
            }
            return *_central_value;
        }

        T error() {
            T &mean = central_value();
            T ret = elementwise_square(values[0] - mean);
            for (const T& val: values) {
                if (&val == &values[0])
                    continue;
                ret += elementwise_square(val - mean);
            }
            const RealD N = values.size();
            ret *= (N - 1) / N;
            return elementwise_sqrt(ret);
        }
};

template <typename T>
JackknifeDatabase<T> make_jackknife(const std::vector<T> &vals, int bin_size = 1) {
    assert(vals.size() % bin_size == 0 && bin_size > 0);
    int Nbins = vals.size() / bin_size;

    RealD N = vals.size();
    T total = vals[0];
    for (int i = 1; i < vals.size(); i++) {
        total += vals[i];
    }

    std::vector<T> averages;
    averages.resize(Nbins);

    if (vals.size() == 1) {
        // When we only have one sample the below formulae make everything NaN.
        // By setting things up this way we end up with a correct result with 0
        // error, which is a more reasonable result.
        averages[0] = total;
        return JackknifeDatabase<T>(averages);
    }

    parallel_for (int i = 0; i < Nbins; i++) {
        averages[i] = total;
        for (int j = 0; j < bin_size; j++) {
            averages[i] -= vals[i * bin_size + j];
        }
        averages[i] *= 1.0 / (N - bin_size);
    }

    return JackknifeDatabase<T>(averages);
}

template <int N>
inline OutputData<N> convert_to_output(JackknifeDatabase<Eigen::MatrixXcd> &jack,
        std::vector<int> &trajectories) {
    OutputData<N> ret;
    ret.central_value = jack.central_value();
    ret.values = jack.values;
    ret.trajectories = trajectories;
    return ret;
}

template <int N>
inline JackknifeDatabase<Eigen::MatrixXcd> convert_from_output(OutputData<N> &output) {
    JackknifeDatabase<MatrixXcd> ret;
    *ret._central_value = output.central_value;
    ret.values = output.values;
    return ret;
}

#endif

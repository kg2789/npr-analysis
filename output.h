#ifndef __OUTPUT_H__

#define __OUTPUT_H__

#include <iomanip>
#include <limits>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include <Grid/Eigen/Dense>


using Eigen::MatrixXcd;

template <int N>
inline void write_matrix(std::ostream &os, const MatrixXcd &mat) {
    assert(mat.rows() == N && mat.cols() == N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            os << mat(i,j);
        }
    }
}

template <int N>
inline MatrixXcd read_matrix(std::istream &is) {
    MatrixXcd ret(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            is >> ret(i,j);
        }
    }
    return ret;
}

template <int N>
class OutputData {
    public:
        MatrixXcd central_value;
        std::vector<MatrixXcd> values;
        std::vector<int> trajectories;

        OutputData() {}

        void write_to_file(const std::string &filename) {
            assert(values.size() == trajectories.size());

            std::ofstream file(filename);
            int precision = std::numeric_limits<double>::max_digits10;
            file << std::setprecision(precision);
            write_matrix<N>(file, central_value);
            for (size_t i = 0; i < values.size(); i++) {
                auto &val = values[i];
                int traj = trajectories[i];

                file << traj;
                write_matrix<N>(file, val);
            }
        }

        void read_from_file(const std::string &filename) {
            std::ifstream file(filename);
            central_value = read_matrix<N>(file);
            MatrixXcd tmp;
            while (true) {
                int traj;
                file >> traj;
                if (file.eof()) break;
                tmp = read_matrix<N>(file);
                // If the file ends here then we found a trajectory number
                // without a corresponding matrix, which means something is
                // wrong with the data
                assert(!file.eof());

                values.push_back(tmp);
                trajectories.push_back(traj);
            }
        }
};

#endif

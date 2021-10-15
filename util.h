#ifndef __UTIL_H__

#define __UTIL_H__

#include <Grid/Grid.h>
#include "operator.h"
#include "metadata.h"

#define parallel_for _Pragma("openmp parallel for") for

using Grid::SpinColourSpinColourMatrix;
using Grid::SpinColourMatrix;
using Grid::SpinMatrix;

void scsc_index_swap(SpinColourSpinColourMatrix& mat, int a, int b);

SpinColourMatrix spin_colour_identity();
SpinMatrix spin_identity();
SpinMatrix gamma_matrix(int mu);
SpinMatrix gamma5();
// Computes v_\mu \gamma^\mu
SpinMatrix vector_slash(std::vector<RealD> &v);
SpinColourMatrix vector_slash_sc(std::vector<RealD> &v);

SpinColourMatrix invert_sc(const SpinColourMatrix &mat);

// Takes a momentum in units of 2 pi / L_mu and returns the momentum in lattice
// units
std::vector<double> momentum_to_lattice_units(const std::vector<double> &n_mu,
        const Metadata &meta);
double vector_magnitude(const std::vector<double> &v);
// Computes q = p_1 - p_2 in physical units
std::vector<double> get_lattice_q(const Metadata &meta);

SpinColourSpinColourMatrix
vertex_structure(const SpinStructure &a, const SpinStructure &b,
        const ColorStructure &structure);
SpinColourSpinColourMatrix
tensor_prod_with_color_structure(const SpinMatrix &a, const SpinMatrix &b,
        ColorStructure color);
SpinColourSpinColourMatrix
spin_mixed_tensor_prod_with_color_structure(const SpinMatrix &a, const SpinMatrix &b,
        ColorStructure color);

void tensor_prod(SpinColourSpinColourMatrix &lret,
        const SpinColourMatrix &a,
        const SpinColourMatrix &b);

void print_scsc(const SpinColourSpinColourMatrix &mat);
void print_scsc_real_part(const SpinColourSpinColourMatrix &mat);
void print_scsc_nonzero(const SpinColourSpinColourMatrix &mat);
void print_sc(const SpinColourMatrix &mat);
void print_sc_real(const SpinColourMatrix &mat);
void print_scsc_spin(const SpinColourSpinColourMatrix &mat);
void print_scsc_spin_real(const SpinColourSpinColourMatrix &mat);

// Primarily used for reading in matrices from Greg's code for debug purposes
inline std::istream& operator>> (std::istream& istr,
        SpinColourSpinColourMatrix& mat) {
    using Grid::Ns;
    using Grid::Nc;

    for (int is = 0; is < Ns; is++) {
        for (int ic = 0; ic < Nc; ic++) {
    for (int js = 0; js < Ns; js++) {
        for (int jc = 0; jc < Nc; jc++) {
    for (int ks = 0; ks < Ns; ks++) {
        for (int kc = 0; kc < Nc; kc++) {
    for (int ls = 0; ls < Ns; ls++) {
        for (int lc = 0; lc < Nc; lc++) {
            istr >> mat()(is,js)(ic,jc)(ks,ls)(kc,lc);
        }}}}
    }}}}

    return istr;
}

#endif

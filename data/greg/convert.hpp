#include "DoubleWilsonMatrix.h"
#include "WilsonMatrix.h"
#include <Grid/Grid.h>

using Grid::SpinColourSpinColourMatrix;
using Grid::SpinColourMatrix;
using Grid::Nc;
using Grid::Ns;

inline SpinColourSpinColourMatrix
double_wilson_to_scsc(const DoubleWilsonMatrix &dwm) {
    using Grid::Nc;
    using Grid::Ns;

    SpinColourSpinColourMatrix ret;
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int s2 = 0; s2 < Ns; s2++) {
    for (int c1 = 0; c1 < Nc; c1++) {
    for (int c2 = 0; c2 < Nc; c2++) {
    for (int s3 = 0; s3 < Ns; s3++) {
    for (int s4 = 0; s4 < Ns; s4++) {
    for (int c3 = 0; c3 < Nc; c3++) {
    for (int c4 = 0; c4 < Nc; c4++) {
        ret()(s1,s2)(c1,c2)(s3,s4)(c3,c4)
            = dwm(s1, c1, s2, c2, s3, c3, s4, c4);
    }}}}}}}}
    return ret;
}

inline DoubleWilsonMatrix
scsc_to_double_wilson(const SpinColourSpinColourMatrix &scsc) {
    using Grid::Nc;
    using Grid::Ns;

    DoubleWilsonMatrix ret;
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int s2 = 0; s2 < Ns; s2++) {
    for (int c1 = 0; c1 < Nc; c1++) {
    for (int c2 = 0; c2 < Nc; c2++) {
    for (int s3 = 0; s3 < Ns; s3++) {
    for (int s4 = 0; s4 < Ns; s4++) {
    for (int c3 = 0; c3 < Nc; c3++) {
    for (int c4 = 0; c4 < Nc; c4++) {
        ret(s1, c1, s2, c2, s3, c3, s4, c4)
            = scsc()(s1,s2)(c1,c2)(s3,s4)(c3,c4);
    }}}}}}}}
    return ret;
}

inline SpinColourMatrix
wilson_to_sc(const WilsonMatrix &mat) {
    using Grid::Nc;
    using Grid::Ns;

    SpinColourMatrix ret;
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int s2 = 0; s2 < Ns; s2++) {
    for (int c1 = 0; c1 < Nc; c1++) {
    for (int c2 = 0; c2 < Nc; c2++) {
        ret()(s1,s2)(c1,c2) = mat(s1, c1, s2, c2);
    }}}}
    return ret;
}

inline WilsonMatrix
sc_to_wilson(const SpinColourMatrix &mat) {
    using Grid::Nc;
    using Grid::Ns;

    WilsonMatrix ret;
    for (int s1 = 0; s1 < Ns; s1++) {
    for (int s2 = 0; s2 < Ns; s2++) {
    for (int c1 = 0; c1 < Nc; c1++) {
    for (int c2 = 0; c2 < Nc; c2++) {
        ret(s1, c1, s2, c2) = mat()(s1,s2)(c1,c2);
    }}}}
    return ret;
}


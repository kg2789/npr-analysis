#ifndef UTIL_H
#define UTIL_H

#include <cstdio>

extern const double PI;

enum Parity { POSITIVE_PARITY, NEGATIVE_PARITY, BOTH_PARITIES };

enum LRStructure { GammaLL = 0, GammaRR = 1, GammaLR = 2, GammaRL = 3 };
extern const char *LR_names[];

enum VAStructure { GammaVV = 0, GammaAA = 1, GammaVA = 2, GammaAV = 3 };
extern const char *VA_names[];

enum Diagram
{
    FULLY_CONNECTED_COLOR_DIAG, FULLY_CONNECTED_COLOR_MIXED,
    DISCONNECTED_LOOP_COLOR_DIAG, DISCONNECTED_LOOP_COLOR_MIXED,
    CONNECTED_LOOP_COLOR_DIAG, CONNECTED_LOOP_COLOR_MIXED
};
extern const char *diagram_names[];

FILE* OpenFile(const char* filename, const char* mode);

void SprintfValueAndError(char *str, double val, double err);

void ReverseDoubleEndianness(double *arr, int N);


bool EndsWithBin(const char* str);


#endif

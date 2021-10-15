#include "util.h"

#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstring>

const double PI = 3.1415926535897932384626433832;

const char *LR_names[] = { "L_L", "R_R", "L_R", "R_L" };

const char *VA_names[] = { "V_V", "A_A", "V_A", "A_V" };

const char *diagram_names[] = { "fully_connected_color_diag", "fully_connected_color_mixed",
"disconnected_loop_color_diag", "disconnected_loop_color_mixed",
"connected_loop_color_diag", "connected_loop_color_mixed" };

FILE* OpenFile(const char* filename, const char* mode)
{
    FILE *f = fopen(filename, mode);
    if (!f) {
	printf("couldn't open %s\n", filename);
	exit(-1);
    }
    return f;
}


void SprintfValueAndError(char *str, double val, double err)
{
    if (err == 0) {
	sprintf(str, "%0.6f(0)", val);
    } else if (err >= 10) {
	sprintf(str, "%d(%d)", (int)round(val), (int)round(err));
    } else if (err >= 1.0) {
	sprintf(str, "%0.1f(%0.1f)", val, err);
    } else {
	int num_leading_zeros_in_err = -(int)std::floor(std::log10(err)) - 1;
	int err_digits = (int)round(err * std::pow(10, num_leading_zeros_in_err + 2));
	int num_digits_past_decimal = num_leading_zeros_in_err + 2;
	char format_string[128];
	sprintf(format_string, "%%0.%df(%%d)", num_digits_past_decimal);
	sprintf(str, format_string, val, err_digits);
    }
}



static void swap(unsigned char *a, unsigned char *b)
{
    char tmp = *a;
    *a = *b;
    *b = tmp;
}

void ReverseDoubleEndianness(double *arr, int N)
{
    unsigned char *v = (unsigned char *)arr;

    for (int i = 0; i < N; ++i) {
	unsigned char* base = v + 8 * i;
	swap(base, base + 7);
	swap(base + 1, base + 6);
	swap(base + 2, base + 5);
	swap(base + 3, base + 4);
    }
}



bool EndsWithBin(const char* str)
{
    int len = strlen(str);
    if (len < 3) return false;
    return str[len - 3] == 'b' && str[len - 2] == 'i' && str[len - 1] == 'n';
}

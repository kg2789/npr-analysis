#include "../jackknife.h"

const double EPSILON = 0.0001;

template <typename T>
bool float_equals(const T&a, const T&b) {
    return std::abs(a - b) < EPSILON;
}

template <typename T>
bool assert_eq(const T&a, const T&b) {
    if (float_equals(a, b))
        return true;

    std::cout << "Error: " << a << " != " << b << std::endl;
    return false;
}

bool test_jackknife() {
    std::vector<double> vals = { 1.0, 6.5, 2.8, 4.2, 6.8 };

    JackknifeDatabase<double> jack = make_jackknife(vals);

    double avg = jack.central_value();
    double err = jack.error();

    // All values are calculated by hand
    if (!assert_eq(avg, 4.26))
        return false;

    if (!assert_eq(jack.values[0], 5.075))
        return false;

    if (!assert_eq(jack.values[1], 3.7))
        return false;

    if (!assert_eq(jack.values[2], 4.625))
        return false;

    if (!assert_eq(jack.values[3], 4.275))
        return false;

    if (!assert_eq(jack.values[4], 3.625))
        return false;

    if (!assert_eq(err, 1.1007270325))
        return false;

    return true;
}

int main(int argc, char **argv) {
    if (!test_jackknife())
        return 1;

    return 0;
}

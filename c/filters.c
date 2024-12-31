// Copyright (c) 2024 silverslither.

#include "filters.h"
#include <math.h>
#include <stdint.h>

double Triangle(double x) {
    return 1.0 - x;
}

double Hermite(double x) {
    return 1.0 - x * x * (3.0 - 2.0 * x);
}

double Lagrange2(double x) {
    if (x <= 0.5)
        return 1.0 - x * x;
    return 1.0 - x * (1.5 - 0.5 * x);
}

// Normalized to 3.0
double Lagrange3(double x) {
    if (x <= 1.0)
        return 3.0 - x * (1.5 + x * (3.0 - 1.5 * x));
    return 3.0 - x * (5.5 - x * (3.0 - 0.5 * x));
}

double BSpline2(double x) {
    if (x <= 0.5)
        return 0.75 - x * x;
    const double x_ = x - 1.5;
    return 0.5 * x_ * x_;
}

// Normalized to 6.0
double BSpline3(double x) {
    if (x <= 1.0)
        return 4.0 - x * x * (6.0 - 3.0 * x);
    const double x_ = 2.0 - x;
    return x_ * x_ * x_;
}

// Normalized to 9.0
double MitNet(double x) {
    if (x <= 1.0)
        return 8.0 - x * x * (18.0 - 10.5 * x);
    return 16.0 - x * (30.0 - x * (18.0 - 3.5 * x));
}

double CatRom(double x) {
    if (x <= 1.0)
        return 1.0 - x * x * (2.5 - 1.5 * x);
    return 2.0 - x * (4.0 - x * (2.5 - 0.5 * x));
}

static inline double lsin3(double x) {
    const double x2 = x * x;
    const double v = 2.829828552115177 - 1.2490259408560183 * x2;
    return x * (1.7320508075688772 - v * x2);
}
static inline double lsin3w(double x) {
    const double x2 = x * x;
    const double v = 0.10480846489315472 - 0.005140024447967154 * x2;
    return x * (0.5773502691896257 - v * x2);
}
double Lanczos3(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int round_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(round_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)round_x);
    return lsin3(poly_x) * lsin3w(x > 1.5 ? 3.0 - x : x) / (x * x);
}

static inline double lsin4(double x) {
    const double x2 = x * x;
    const double v = 3.2676045526483732 - 1.4422509263560956 * x2;
    return x * (2.0 - v * x2);
}
static inline double lsin4w(double x) {
    const double x2 = x * x;
    const double v = 0.05105632113513083 - 0.0014084481702696247 * x2;
    return x * (0.5 - v * x2);
}
double Lanczos4(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int round_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(round_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)round_x);
    return lsin4(poly_x) * lsin4w(x > 2.0 ? 4.0 - x : x) / (x * x);
}

static inline double hsin(double x) {
    const double x2 = x * x;
    const double v = 1.6338022763241866 - 0.7211254631780478 * x2;
    return x * (1.0 - v * x2);
}

static inline double hcos3w(double x) {
    const double x2 = x * x;
    const double v = 0.24920390748853422 - 0.019569144068978167 * x2;
    return 0.46164 - v * x2;
}
double Hamming3(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int round_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(round_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)round_x);
    const double w = copysign(hcos3w(x > 1.5 ? x - 3.0 : x), 1.5 - x);
    return hsin(poly_x) * (0.53836 + w) / x;
}

static inline double hcos4w(double x) {
    const double x2 = x * x;
    const double v = 0.1401771979623005 - 0.006191799490575123 * x2;
    return 0.46164 - v * x2;
}
double Hamming4(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int round_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(round_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)round_x);
    const double w = copysign(hcos4w(x > 2.0 ? x - 4.0 : x), 2.0 - x);
    return hsin(poly_x) * (0.53836 + w) / x;
}

double L_bspline2i[11] = {
    0.14583333333333334,
    0.17142857142857143,
    0.1715686274509804,
    0.17157275021026072,
    0.17157287157287157,
    0.1715728751454532,
    0.17157287525062018,
    0.171572875253716,
    0.17157287525380713,
    0.17157287525380982,
    0.1715728752538099
};

double L_bspline3i[14] = {
    0.20833333333333334,
    0.26666666666666666,
    0.26785714285714285,
    0.2679425837320574,
    0.26794871794871794,
    0.2679491583648231,
    0.26794918998527245,
    0.26794919225551855,
    0.2679491924185149,
    0.26794919243021753,
    0.2679491924310577,
    0.26794919243111803,
    0.26794919243112236,
    0.2679491924311227
};

double L_omoms3[18] = {
    0.2490842490842491,
    0.33986928104575165,
    0.34362717574396406,
    0.34407148031876356,
    0.34412408743959544,
    0.34413031736062233,
    0.3441310551448089,
    0.34413114251779625,
    0.34413115286505125,
    0.3441311540904378,
    0.3441311542355558,
    0.34413115425274154,
    0.3441311542547768,
    0.34413115425501783,
    0.34413115425504637,
    0.34413115425504975,
    0.34413115425505014,
    0.3441311542550502
};

// Normalized to 5.25
double OMOMS3(double x) {
    if (x <= 1.0)
        return 3.25 + x * (0.375 - x * (5.25 - 2.625 * x));
    return 7.25 - x * (10.625 - x * (5.25 - 0.875 * x));
}

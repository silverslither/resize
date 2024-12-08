// Copyright (c) 2024 silverslither.

#include "filters.h"
#include <stdint.h>

double Triangle(double x) {
    return 1.0 - x;
}

double Hermite(double x) {
    return 1.0 - x * x * (3.0 - 2.0 * x);
}

double BSpline2(double x) {
    if (x <= 0.5)
        return 0.75 - x * x;
    const double x_ = x - 1.5;
    return 0.5 * x_ * x_;
}

double BSpline3(double x) {
    if (x <= 1.0)
        return 0.6666666666666666 - x * x * (1.0 - 0.5 * x);
    const double x_ = x - 2;
    return -0.16666666666666666 * x_ * x_ * x_;
}

double MitNet(double x) {
    if (x <= 1.0)
        return 0.8888888888888888 - x * x * (2.0 - 1.1666666666666667 * x);
    return 1.7777777777777777 - x * (3.3333333333333333 - x * (2.0 - 0.3888888888888889 * x));
}

double CatRom(double x) {
    if (x <= 1.0)
        return 1.0 - x * x * (2.5 - 1.5 * x);
    return 2.0 - x * (4.0 - x * (2.5 - 0.5 * x));
}

double MKS2013(double x) {
    if (x <= 0.5)
        return 1.0625 - 1.75 * x * x;
    if (x <= 1.5)
        return 1.75 - x * (2.75 - x);
    const double x_ = x - 2.5;
    return -0.125 * x_ * x_;
}

static inline double sinc3(double x) {
    const double x2 = x * x;
    const double v = 2.829828552115177 - 1.2490259408560183 * x2;
    return x * (1.7320508075688772 - v * x2);
}
static inline double sinc3w(double x) {
    const double x2 = x * x;
    const double v = 0.10480846489315472 - 0.005140024447967154 * x2;
    return x * (0.5773502691896257 - v * x2);
}
double Lanczos3(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int floor_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(~floor_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)floor_x);
    return sinc3(poly_x) * sinc3w(x > 1.5 ? 3.0 - x : x) / (x * x);
}

static inline double sinc4(double x) {
    const double x2 = x * x;
    const double v = 3.2676045526483732 - 1.4422509263560956 * x2;
    return x * (2.0 - v * x2);
}
static inline double sinc4w(double x) {
    const double x2 = x * x;
    const double v = 0.05105632113513083 - 0.0014084481702696247 * x2;
    return x * (0.5 - v * x2);
}
double Lanczos4(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int floor_x = (int)(x + 0.5);
    const uint64_t sign = ((uint64_t)(~floor_x) << 63) | 0x3ff0000000000000;
    const double poly_x = *(const double *)&sign * (x - (double)floor_x);
    return sinc4(poly_x) * sinc4w(x > 2.0 ? 4.0 - x : x) / (x * x);
}

double L_bspline3i[15] = {
    0.2,
    0.2631578947368421,
    0.2676056338028169,
    0.2679245283018868,
    0.2679474216380182,
    0.2679490652939583,
    0.2679491833030853,
    0.2679491917757591,
    0.2679491923840697,
    0.26794919242774445,
    0.26794919243088017,
    0.26794919243110527,
    0.2679491924311215,
    0.26794919243112264,
    0.2679491924311227
};
double pNormBSpline3(double x) {
    if (x <= 1.0)
        return 4.0 - x * x * (6.0 - 3.0 * x);
    const double x_ = 2 - x;
    return x_ * x_ * x_;
}

double L_omoms3[18] = {
    0.23529411764705882,
    0.33170731707317075,
    0.34266610948600085,
    0.34395774192389234,
    0.34411061890821637,
    0.34412872234809083,
    0.34413086625366646,
    0.34413112014812935,
    0.34413115021589546,
    0.3441311537767083,
    0.34413115419840196,
    0.34413115424834156,
    0.3441311542542557,
    0.3441311542549561,
    0.34413115425503904,
    0.34413115425504887,
    0.34413115425505003,
    0.3441311542550502
};
double pNormOMOMS3(double x) {
    if (x <= 1.0)
        return 3.25 + x * (0.375 - x * (5.25 - 2.625 * x));
    return 7.25 - x * (10.625 - x * (5.25 - 0.875 * x));
}

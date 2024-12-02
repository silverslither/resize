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

double KeysHalf(double x) {
    if (x <= 1.0)
        return 0.8333333333333333 - x * x * (1.75 - x);
    return 1.6666666666666667 - x * (3.0 - x * (1.75 - 0.3333333333333333 * x));
}

double MitNet(double x) {
    if (x <= 1.0)
        return 0.8888888888888888 - x * x * (2.0 - 1.1666666666666667 * x);
    return 1.7777777777777777 - x * (3.3333333333333333 - x * (2.0 - 0.3888888888888889 * x));
}

double MitNetSharp(double x) {
    if (x <= 1.0)
        return 0.9166666666666666 - x * x * (2.125 - 1.25 * x);
    return 1.8333333333333333 - x * (3.5 - x * (2.125 - 0.4166666666666667 * x));
}

double CatRom(double x) {
    if (x <= 1.0)
        return 1.0 - x * x * (2.5 - 1.5 * x);
    return 2.0 - x * (4.0 - x * (2.5 - 0.5 * x));
}

double CatRomSharp(double x) {
    if (x <= 1.0)
        return 1.0 - x * x * (2.4 - 1.4 * x);
    return 2.4 - x * (4.8 - x * (3 - 0.6 * x));
}

double MagicKernelSharp2013(double x) {
    if (x <= 0.5)
        return 1.0625 - 1.75 * x * x;
    if (x <= 1.5)
        return 1.75 - x * (2.75 - x);
    const double x_ = x - 2.5;
    return -0.125 * x_ * x_;
}

static double sinc3(double x) {
    const double x2 = x * x;
    const double v = 0.19010152698956836 - 0.00932297307587116 * x2;
    return 1.0471975511965979 - v * x2;
}
static double lnorm3(double x) {
    const double c0 = x;
    const double c0_ = 3.0 - x;
    const double c1 = x >= 0.5 ? 2.0 - x : x + 1.0;
    const double c1_ = 3.0 - c1;
    const double c2 = 1.0 - x;
    const double c2_ = x + 2.0;
    const double o0 = sinc3(c0);
    const double o1 = sinc3(c1);
    const double o2 = sinc3(c2);
    double v = 0.0;
    v -= o2 / c2;
    v += o1 / c1;
    v -= o0 / c0;
    v -= c0 * o0 / (c0_ * c0_);
    v += c1 * o1 / (c1_ * c1_);
    v -= c2 * o2 / (c2_ * c2_);
    return v;
}
double Lanczos3(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int floor_x = x;
    const uint64_t sign = ((uint64_t)(~floor_x) << 63) | 0x3ff0000000000000;
    const double poly_x = x > 1.5 ? 3.0 - x : x;
    const double numerator = sinc3(poly_x) * poly_x / (x * x);
    const double denominator = lnorm3(x - floor_x);
    return *(double *)&sign * numerator / denominator;
}

static double sinc4(double x) {
    const double x2 = x * x;
    const double v = 0.08019908169872415 - 0.002212385212340519 * x2;
    return 0.7853981633974483 - v * x2;
}
static double lnorm4(double x) {
    const double c0 = x;
    const double c0_ = 4.0 - x;
    const double c1 = x + 1.0;
    const double c1_ = 3.0 - x;
    const double c2 = 2.0 - x;
    const double c2_ = x + 2.0;
    const double c3 = 1.0 - x;
    const double c3_ = x + 3.0;
    const double o0 = sinc4(c0);
    const double o1 = sinc4(c1);
    const double o2 = sinc4(c2);
    const double o3 = sinc4(c3);
    double v = 0.0;
    v -= o3 / c3;
    v += o2 / c2;
    v -= o0 / c0;
    v += o1 / c1;
    v -= c1 * o1 / (c1_ * c1_);
    v += c0 * o0 / (c0_ * c0_);
    v -= c2 * o2 / (c2_ * c2_);
    v += c3 * o3 / (c3_ * c3_);
    return v;
}
double Lanczos4(double x) {
    if (x < 1.0e-8)
        return 1.0;
    const int floor_x = x;
    const uint64_t sign = ((uint64_t)(~floor_x) << 63) | 0x3ff0000000000000;
    const double poly_x = x > 2.0 ? 4.0 - x : x;
    const double numerator = sinc4(poly_x) * poly_x / (x * x);
    const double denominator = lnorm4(x - floor_x);
    return *(double *)&sign * numerator / denominator;
}

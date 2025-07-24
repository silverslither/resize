#include "colour.h"
#include <math.h>

void mul_alpha(double *arr, size_t len) {
    const double *end = arr + len;
    while (arr < end) {
        arr[0] *= arr[3];
        arr[1] *= arr[3];
        arr[2] *= arr[3];
        arr += 4;
    }
}

void div_alpha(double *arr, size_t len) {
    const double *end = arr + len;
    while (arr < end) {
        if (arr[3] != 0) {
            arr[0] /= arr[3];
            arr[1] /= arr[3];
            arr[2] /= arr[3];
        }
        arr += 4;
    }
}

static inline double _srgb_encode_polyapprox(double x) {
    x = sqrt(sqrt(x));
    double v = -0.9379296312853258;
    v = v * x + 9.2771292175317921;
    v = v * x + -42.538447026088356;
    v = v * x + 120.03418911775015;
    v = v * x + -233.36351766806072;
    v = v * x + 331.64911228504167;
    v = v * x + -356.81113272500056;
    v = v * x + 297.04753668361309;
    v = v * x + -194.09552975820117;
    v = v * x + 100.55408719624158;
    v = v * x + -41.732365455317492;
    v = v * x + 14.171258154510459;
    v = v * x + -4.2475018941351363;
    v = v * x + 1.9789474424876097;
    v = v * x + 0.070073818045790701;
    return v * x + -0.055909757142145876;
}

static inline double _srgb_decode_polyapprox(double x) {
    double v = -35.214090542702046;
    v = v * x + 381.58899386463918;
    v = v * x + -1928.8060440754307;
    v = v * x + 6042.0071934476673;
    v = v * x + -13144.742650923334;
    v = v * x + 21097.422483855;
    v = v * x + -25907.037146644285;
    v = v * x + 24918.290546697084;
    v = v * x + -19071.343465115136;
    v = v * x + 11742.31511660925;
    v = v * x + -5863.0107704639386;
    v = v * x + 2390.382658684478;
    v = v * x + -802.17473499797052;
    v = v * x + 224.50267306872107;
    v = v * x + -53.813907435619484;
    v = v * x + 11.736268792279825;
    v = v * x + -2.7110681752833026;
    v = v * x + 1.1071900562656445;
    v = v * x + 0.46354401297592546;
    v = v * x + 0.036375405673038159;
    return v * x + 0.00083387965612594117;
}

double srgb_encode(double x) {
    if (x > 0.003130668442500634)
        return _srgb_encode_polyapprox(x);
    return 12.92 * x;
}

double srgb_decode(double x) {
    if (x > 0.04044823627710819)
        return _srgb_decode_polyapprox(x);
    return 0.07739938080495357 * x;
}

SigmoidizationParams get_sigmoidization_params(double beta) {
    beta *= 1.4426950408889634;
    SigmoidizationParams params;
    double exp_beta = exp2(beta);
    double sqrt_exp_beta = sqrt(exp_beta);
    double inv_denominator = 1.0 / (sqrt_exp_beta - 1.0);
    params.beta = beta;
    params.inv_beta = 1.0 / beta;
    params.a = sqrt_exp_beta * inv_denominator;
    params.b = (exp_beta + sqrt_exp_beta) * inv_denominator;
    params.c = sqrt_exp_beta;
    return params;
}

double sigmoidal_contrast_increase(double x, SigmoidizationParams *params) {
    return params->a - params->b / (exp2(params->beta * x) + params->c);
}

double sigmoidal_contrast_decrease(double x, SigmoidizationParams *params) {
    return log2(params->b / (params->a - x) - params->c) * params->inv_beta;
}

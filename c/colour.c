#include "colour.h"
#include <math.h>

double srgb_encode(double x) {
    if (x > 0.003130668442500634)
        return 1.055 * pow(x, 0.4166666666666667) - 0.055;
    return 12.92 * x;
}

double srgb_decode(double x) {
    if (x > 0.04044823627710819)
        return pow(0.9478672985781991 * (x + 0.055), 2.4);
    return 0.07739938080495357 * x;
}

void get_sigmoidization_params(double beta, SigmoidizationParams *params) {
    double exp_beta = exp2(beta);
    double sqrt_exp_beta = sqrt(exp_beta);
    double inv_denominator = 1.0 / (sqrt_exp_beta - 1.0);
    params->beta = beta;
    params->inv_beta = 1.0 / beta;
    params->a = sqrt_exp_beta * inv_denominator;
    params->b = (exp_beta + sqrt_exp_beta) * inv_denominator;
    params->c = sqrt_exp_beta;
}

double sigmoidal_contrast_increase(double x, SigmoidizationParams *params) {
    return params->a - params->b / (exp2(params->beta * x) + params->c);
}

double sigmoidal_contrast_decrease(double x, SigmoidizationParams *params) {
    return log2(params->b / (params->a - x) - params->c) * params->inv_beta;
}

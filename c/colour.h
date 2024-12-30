#ifndef RESIZE_COLOUR
#define RESIZE_COLOUR

#include <stddef.h>

typedef struct SigmoidizationParams {
    double beta;
    double inv_beta;
    double a;
    double b;
    double c;
} SigmoidizationParams;

void mul_alpha(double *arr, size_t len);
void div_alpha(double *arr, size_t len);

double srgb_encode(double x);
double srgb_decode(double x);

SigmoidizationParams get_sigmoidization_params(double beta);
double sigmoidal_contrast_increase(double x, SigmoidizationParams *params);
double sigmoidal_contrast_decrease(double x, SigmoidizationParams *params);

#endif

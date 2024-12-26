#ifndef RESIZE_COLOUR
#define RESIZE_COLOUR

typedef struct SigmoidizationParams {
    double beta;
    double inv_beta;
    double a;
    double b;
    double c;
} SigmoidizationParams;

double srgb_encode(double x);
double srgb_decode(double x);

void get_sigmoidization_params(double beta, SigmoidizationParams *params);
double sigmoidal_contrast_increase(double x, SigmoidizationParams *params);
double sigmoidal_contrast_decrease(double x, SigmoidizationParams *params);

#endif

// Copyright (c) 2024 silverslither.

#include "colour.h"
#include "include/lodepng/lodepng.h"
#include "resize.h"
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static double *u8_to_f64(const unsigned char *arr, size_t len, int linearize, SigmoidizationParams *params) {
    double *out = malloc(len << 3);
    for (size_t i = 0; i < len; i++) {
        double temp = 0.00392156862745098 * (double)arr[i];
        if (i % 4 != 3) {
            if (linearize)
                temp = srgb_decode(temp);
            if (params)
                temp = sigmoidal_contrast_decrease(temp, params);
        }
        out[i] = temp;
    }
    return out;
}

static unsigned char *f64_to_u8(const double *arr, size_t len, int linearize, SigmoidizationParams *params) {
    unsigned char *out = malloc(len);
    for (size_t i = 0; i < len; i++) {
        double temp = arr[i];
        if (i % 4 != 3) {
            if (params)
                temp = sigmoidal_contrast_increase(temp, params);
            if (linearize)
                temp = srgb_encode(temp);
        }
        temp = temp > 0.0 ? temp : 0.0;
        temp = temp < 1.0 ? temp : 1.0;
        out[i] = (unsigned char)__builtin_roundeven(255.0 * temp);
    }
    return out;
}

int parseDimension(char *str, unsigned src_dimension) {
    int len = strlen(str);
    double num;
    int pos;

    if (!sscanf(str, "%lf%n", &num, &pos) || num <= 0 || pos < len - 1)
        return -1;

    if (str[len - 1] == '%')
        num *= 0.01 * (double)src_dimension;
    if (tolower(str[len - 1]) == 'x')
        num *= (double)src_dimension;

    if (num >= INT_MAX + 0.5)
        return -1;

    return (int)__builtin_roundeven(num);
}

Filter parseFilter(char *str) {
    static char filters[][20] = {
        "DEFAULT",
        "NEAREST",
        "AREA",
        "TRIANGLE",
        "HERMITE",
        "LAGRANGE2",
        "LAGRANGE3",
        "BSPLINE2",
        "BSPLINE3",
        "MITNET",
        "CATROM",
        "MKS2013",
        "LANCZOS3",
        "LANCZOS4",
        "HAMMING3",
        "HAMMING4",
        "BSPLINE3I",
        "OMOMS3I"
    };
    static int len = sizeof(filters) / 20;

    char *p = str - 1;
    while (*++p)
        *p = toupper(*p);

    for (int i = 0; i < len; i++)
        if (strcmp(str, filters[i]) == 0)
            return i;

    fprintf(stderr, "warning: invalid filter '%s'\n", str);

    return 0;
}

SigmoidizationParams *parseSigmoidizationBeta(char *str, SigmoidizationParams *params) {
    int len = strlen(str);
    double beta;
    int pos;

    if (!sscanf(str, "%lf%n", &beta, &pos) || beta < 0 || pos != len) {
        fprintf(stderr, "warning: invalid sigmoidization contrast '%s'\n", str);
        return NULL;
    }

    if (beta == 0)
        return NULL;

    get_sigmoidization_params(beta, params);
    return params;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: resize <input> <width> <height> <output> [-f filter] [-l] [-s contrast]\n");
        return 1;
    }

    // options declaration
    Filter filter = DEFAULT;
    int linearize = 0;
    SigmoidizationParams *sigParamsPtr = NULL;
    SigmoidizationParams _sigParams;

    if (argc == 5)
        goto no_options;

    // options parsing
    for (int i = 5; i < argc; i++) {
        const char *str = argv[i];
        if (str[0] != '-' || str[1] == 0 || str[2] != 0) {
invalid:
            fprintf(stderr, "warning: invalid argument '%s'\n", str);
            continue;
        }

        switch (str[1]) {
        case 'f':
            if (i++ == argc - 1) {
                fprintf(stderr, "warning: no filter given\n");
                break;
            }
            filter = parseFilter(argv[i]);
            break;
        case 'l':
            linearize = 1;
            break;
        case 's':
            if (i++ == argc - 1) {
                fprintf(stderr, "warning: no sigmoidization contrast given\n");
                break;
            }
            sigParamsPtr = parseSigmoidizationBeta(argv[i], &_sigParams);
            break;
        default:
            goto invalid;
        }
    }

no_options:;
    // input img
    unsigned char *_img;
    unsigned src_width, src_height;
    unsigned err = lodepng_decode32_file(&_img, &src_width, &src_height, argv[1]);
    if (err) {
        fprintf(stderr, "error: lodepng %u\n", err);
        return 1;
    }
    size_t area = ((size_t)src_width * (size_t)src_height) << 2;

    // dimensions parsing
    int dst_width = parseDimension(argv[2], src_width);
    if (dst_width == -1) {
        fprintf(stderr, "error: invalid dimension '%s'\n", argv[2]);
        return 1;
    }
    int dst_height = parseDimension(argv[3], src_height);
    if (dst_height == -1) {
        fprintf(stderr, "error: invalid dimension '%s'\n", argv[3]);
        return 1;
    }

    if (dst_width > 65536 || dst_height > 65536) {
        printf("dst image is very large, are you sure? [y/N] ");
        if (tolower(getchar()) != 'y')
            return 0;
    }

    // resize
    double *img = u8_to_f64(_img, area, linearize, sigParamsPtr);
    free(_img);
    mul_alpha(img, area);

    double *resized = resize(img, src_width, src_height, dst_width, dst_height, filter);
    free(img);

    area = ((size_t)dst_width * (size_t)dst_height) << 2;
    div_alpha(resized, area);

    _img = f64_to_u8(resized, area, linearize, sigParamsPtr);
    free(resized);

    // output img
    err = lodepng_encode32_file(argv[4], _img, dst_width, dst_height);
    if (err) {
        fprintf(stderr, "error: lodepng %u\n", err);
        return 1;
    }
    free(_img);

    return 0;
}

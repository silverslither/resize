// Copyright (c) 2024 silverslither.

#include "include/lodepng/lodepng.h"
#include "resize.h"
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double *u8_to_f64(const unsigned char *arr, size_t len) {
    double *out = malloc(len << 3);
    for (size_t i = 0; i < len; i++)
        out[i] = (double)arr[i];
    return out;
}

unsigned char *f64_to_u8(const double *arr, size_t len) {
    unsigned char *out = malloc(len);
    for (size_t i = 0; i < len; i++) {
        const double temp = fmin(fmax(arr[i], 0.0), 255.0);
        out[i] = nearbyint(temp);
    }
    return out;
}

Filter parseFilter(char *str) {
    static char filters[][20] = {
        "DEFAULT",
        "NEAREST",
        "AREA",
        "TRIANGLE",
        "HERMITE",
        "BSPLINE2",
        "BSPLINE3",
        "KEYSHALF",
        "MITNET",
        "MITNETSHARP",
        "CATROM",
        "CATROMSHARP",
        "LANCZOS3",
        "LANCZOS4"
    };
    static int len = sizeof(filters) / 20;

    char *p = str - 1;
    while (*++p)
        *p = toupper(*p);

    for (int i = 0; i < len; i++) {
        if (strcmp(str, filters[i]) == 0)
            return i;
    }

    return 0;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        printf("usage: resize <input> <width> <height> <output> [filter name]\n");
        return 0;
    }

    Filter filter = DEFAULT;
    if (argc >= 6)
        filter = parseFilter(argv[5]);

    unsigned dst_width = abs(atoi(argv[2]));
    unsigned dst_height = abs(atoi(argv[3]));

    if (dst_width > 65536 || dst_height > 65536) {
        printf("dst image is very large, are you sure? [y/N] ");
        if (getchar() != 'y')
            return 0;
    }

    unsigned char *_img;
    unsigned src_width, src_height;
    unsigned err = lodepng_decode32_file(&_img, &src_width, &src_height, argv[1]);
    if (err) {
        fprintf(stderr, "lodepng error %u\n", err);
        return err;
    }

    size_t area = ((size_t)src_width * (size_t)src_height) << 2;

    double *img = u8_to_f64(_img, area);
    free(_img);

    double *resized = resize(img, src_width, src_height, dst_width, dst_height, filter);
    area = ((size_t)dst_width * (size_t)dst_height) << 2;
    free(img);

    _img = f64_to_u8(resized, area);
    free(resized);

    err = lodepng_encode32_file(argv[4], _img, dst_width, dst_height);
    if (err) {
        fprintf(stderr, "lodepng error %u\n", err);
        return err;
    }

    return 0;
}

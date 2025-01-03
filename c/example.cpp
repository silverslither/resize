// Copyright (c) 2024-2025 silverslither.

extern "C" {
#include "colour.h"
#include "resize.h"
}

#include <cctype>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define WUFFS_IMPLEMENTATION
#define WUFFS_CONFIG__STATIC_FUNCTIONS

#define WUFFS_CONFIG__MODULES
#define WUFFS_CONFIG__MODULE__AUX__BASE
#define WUFFS_CONFIG__MODULE__AUX__IMAGE
#define WUFFS_CONFIG__MODULE__BASE
#define WUFFS_CONFIG__MODULE__CRC32
#define WUFFS_CONFIG__MODULE__JPEG
#define WUFFS_CONFIG__MODULE__ADLER32
#define WUFFS_CONFIG__MODULE__BMP
#define WUFFS_CONFIG__MODULE__DEFLATE
#define WUFFS_CONFIG__MODULE__ETC2
#define WUFFS_CONFIG__MODULE__GIF
#define WUFFS_CONFIG__MODULE__NETPBM
#define WUFFS_CONFIG__MODULE__NIE
#define WUFFS_CONFIG__MODULE__PNG
#define WUFFS_CONFIG__MODULE__QOI
#define WUFFS_CONFIG__MODULE__TARGA
#define WUFFS_CONFIG__MODULE__THUMBHASH
#define WUFFS_CONFIG__MODULE__VP8
#define WUFFS_CONFIG__MODULE__WBMP
#define WUFFS_CONFIG__MODULE__WEBP
#define WUFFS_CONFIG__MODULE__ZLIB

#define WUFFS_CONFIG__DST_PIXEL_FORMAT__ENABLE_ALLOWLIST
#define WUFFS_CONFIG__DST_PIXEL_FORMAT__ALLOW_RGBA_NONPREMUL

#include "include/fpng.cpp"
#include "include/wuffs-v0.4.c"

class RGBA_NonPremul_DecodeImageCallbacks : public wuffs_aux::DecodeImageCallbacks {
  private:
    wuffs_base__pixel_format SelectPixfmt(const wuffs_base__image_config &image_config) override {
        (void)image_config;
        return wuffs_base__make_pixel_format(WUFFS_BASE__PIXEL_FORMAT__RGBA_NONPREMUL);
    }
};

#undef assert
#define assert(cond, ...)             \
    if (!(cond)) {                    \
        fprintf(stderr, __VA_ARGS__); \
        return 1;                     \
    }                                 \
    static_assert(true, "")

static double *u8_to_f64(const uint8_t *arr, size_t len, bool linearize, SigmoidizationParams *params) {
    double *out = static_cast<double *>(malloc(len << 3));
    for (size_t i = 0; i < len; i++) {
        double temp = 0.00392156862745098 * static_cast<double>(arr[i]);
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

static uint8_t *f64_to_u8(const double *arr, size_t len, bool linearize, SigmoidizationParams *params) {
    uint8_t *out = static_cast<uint8_t *>(malloc(len));
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
        out[i] = static_cast<uint8_t>(__builtin_roundeven(255.0 * temp));
    }
    return out;
}

int32_t parseDimension(char *str, unsigned src_dimension) {
    int len = strlen(str);
    double num;
    int pos;

    if (!sscanf(str, "%lf%n", &num, &pos) || pos < len - 1)
        return -1;

    if (pos == len - 1) {
        if (str[len - 1] == '%')
            num *= 0.01 * static_cast<double>(src_dimension);
        else if (tolower(str[len - 1]) == 'x')
            num *= static_cast<double>(src_dimension);
        else
            return -1;
    }

    if (num <= 0.5 || num >= INT_MAX + 0.5)
        return -1;

    return static_cast<int32_t>(__builtin_roundeven(num));
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
        "MITNET",
        "CATROM",
        "HAMMING3",
        "HAMMING4",
        "BSPLINE2I",
        "BSPLINE3I",
        "OMOMS3I"
    };
    static int len = sizeof(filters) / 20;

    char *p = str - 1;
    while (*++p)
        *p = toupper(*p);

    for (int i = 0; i < len; i++)
        if (strcmp(str, filters[i]) == 0)
            return static_cast<Filter>(i);

    fprintf(stderr, "warning: invalid filter '%s'\n", str);

    return DEFAULT;
}

SigmoidizationParams *parseSigmoidizationBeta(char *str, SigmoidizationParams *params) {
    int len = strlen(str);
    double beta;
    int pos;

    if (!sscanf(str, "%lf%n", &beta, &pos) || beta < 0 || pos != len) {
        fprintf(stderr, "warning: invalid sigmoidization contrast '%s'\n", str);
        return nullptr;
    }

    if (beta == 0)
        return nullptr;

    *params = get_sigmoidization_params(beta);
    return params;
}

int main(int argc, char **argv) {
    assert(argc >= 5, "usage: resize <input> <width> <height> <output> [-f filter] [-l] [-s contrast]\n");

    // default options
    Filter filter = DEFAULT;
    bool linearize = false;
    SigmoidizationParams *sigParamsPtr = nullptr;
    SigmoidizationParams _sigParams;

    if (argc == 5)
        goto no_options;

    // options
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
            linearize = true;
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

no_options:
    // input
    FILE *input_file = fopen(argv[1], "rb");
    assert(input_file, "error: failed to open '%s' for reading\n", argv[2]);

    RGBA_NonPremul_DecodeImageCallbacks callbacks;
    wuffs_aux::sync_io::FileInput input(input_file);
    wuffs_aux::DecodeImageResult result = wuffs_aux::DecodeImage(callbacks, input);
    assert(result.error_message.empty(), "error: wuffs '%s'\n", result.error_message.c_str());
    fclose(input_file);

    uint32_t src_width = result.pixbuf.pixcfg.width();
    uint32_t src_height = result.pixbuf.pixcfg.height();
    uint8_t *input_u8 = result.pixbuf.plane(0).ptr;
    size_t area = (static_cast<size_t>(src_width) * static_cast<size_t>(src_height)) << 2;

    // dimensions
    int32_t dst_width = parseDimension(argv[2], src_width);
    assert(dst_width != -1, "error: invalid width '%s'\n", argv[2]);

    int32_t dst_height = parseDimension(argv[3], src_height);
    assert(dst_height != -1, "error: invalid height '%s'\n", argv[3]);

    if (dst_width > 65536 || dst_height > 65536) {
        printf("dst image is very large, are you sure? [y/N] ");
        if (tolower(getchar()) != 'y')
            return 0;
    }

    // resize
    double *input_f64 = u8_to_f64(input_u8, area, linearize, sigParamsPtr);
    mul_alpha(input_f64, area);

    double *output_f64 = resize(input_f64, src_width, src_height, dst_width, dst_height, filter);
    free(input_f64);

    area = (static_cast<size_t>(dst_width) * static_cast<size_t>(dst_height)) << 2;
    div_alpha(output_f64, area);

    uint8_t *output_u8 = f64_to_u8(output_f64, area, linearize, sigParamsPtr);
    free(output_f64);

    // output
    fpng::fpng_init();
    assert(fpng::fpng_encode_image_to_file(argv[4], output_u8, dst_width, dst_height, 4, fpng::FPNG_ENCODE_SLOWER), "error: failed to write to '%s'\n", argv[4]);
    free(output_u8);

    return 0;
}

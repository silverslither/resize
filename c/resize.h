// Copyright (c) 2024 silverslither.

#ifndef RESIZE_RESIZE
#define RESIZE_RESIZE

#include <stdint.h>

typedef int32_t s32;

typedef enum FILTER {
    DEFAULT,
    NEAREST,
    AREA,
    TRIANGLE,
    HERMITE,
    B_SPLINE_2,
    B_SPLINE_3,
    KEYS_HALF,
    MITNET,
    MITNET_SHARP,
    CATROM,
    CATROM_SHARP,
    LANCZOS_3,
    LANCZOS_4
} FILTER;

double *sample(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

double *scale(double *data, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

double *resize(double *data, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, FILTER filter_name);

#endif

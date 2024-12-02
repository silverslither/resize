// Copyright (c) 2024 silverslither.

#ifndef RESIZE_RESIZE
#define RESIZE_RESIZE

#include <stdint.h>

typedef int32_t s32;

typedef enum Filter {
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
    MKS2013,
    LANCZOS_3,
    LANCZOS_4
} Filter;

/**
 * \brief Resample an image using nearest neighbor interpolation.
 * \param src Source image in RGBA format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \return Destination image in RGBA format, or null pointer if OOM.
 */
double *sample(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

/**
 * \brief Resize an image using area averaging / pixel mixing.
 * \param src Source image in RGBA format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \return Destination image in RGBA format, or null pointer if OOM.
 */
double *scale(double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

/**
 * \brief Resample an image using a reconstruction filter. Also acts as a wrapper for `sample` and `scale`.
 * \param src Source image in RGBA format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Reconstruction filter to be used. `NEAREST` acts as a wrapper for `sample`, and `AREA` acts as a wrapper for `scale`. The default filter used is Mitchell-Netravali.
 * \return Destination image in RGBA format, or null pointer if OOM.
 */
double *resize(double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, Filter filter);

#endif

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
    BSPLINE2,
    BSPLINE3,
    MITNET,
    CATROM,
    HAMMING3,
    HAMMING4,
    BSPLINE2I,
    BSPLINE3I,
    OMOMS3I,
} Filter;

/**
 * \brief Resample an image using nearest neighbor interpolation.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *sample(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

/**
 * \brief Resize an image using area averaging / pixel mixing.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *scale(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height);

/**
 * \brief Resize an image using a reconstruction filter.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Reconstruction filter function.
 * \param window Filter function window.
 * \param norm Normalization constant.
 * \param nop Boolean flag for a no-op case.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *reconstruct(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, double (*filter)(double), double window, double norm, int nop);

/**
 * \brief Resize an image using a reconstruction filter and an inverse discrete convolution.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Reconstruction filter function.
 * \param window Filter function window.
 * \param norm Normalization constant.
 * \param L Lower matrix coefficients.
 * \param m Number of lower matrix coefficients.
 * \param c Edge multiplier constant.
 * \param nop Boolean flag for a no-op case.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *reconstruct_iconvolve(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, double (*filter)(double), double window, double norm, const double *L, int m, double c, int nop);

/**
 * \brief Wrapper for `sample`, `scale`, `reconstruct`, and `reconstruct_iconvolve`.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Resizing method (filter) to be used. Defaults to Mitchell-Netravali.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *resize(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, Filter filter);

#endif

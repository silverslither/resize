// Copyright (c) 2024-2025 silversliher.

#ifndef RESIZE_RESIZE
#define RESIZE_RESIZE

#include <stdbool.h>
#include <stddef.h>

typedef ptrdiff_t pdt;

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
    HAMMING8,
    BSPLINE2I,
    BSPLINE3I,
    OMOMS3I,
    OMOMS7I,
    OMOMS11I,
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
double *sample(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height);

/**
 * \brief Resize an image using area averaging / pixel mixing.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *scale(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height);

/**
 * \brief Convolve an image with a horizontal and vertical kernel.
 * \param src Source image in a 4-channel format.
 * \param width Image width.
 * \param height Image height.
 * \param h_kernel Horizontal kernel, or null pointer if no horizontal convolution is desired.
 * \param v_kernel Vertical kernel, or null pointer if no vertical convolution is desired.
 * \param h_support Support window for the horizontal kernel. Must be an odd number.
 * \param v_support Support window for the vertical kernel. Must be an odd number.
 * \param h_support Normalization constant for the horizontal kernel.
 * \param v_support Normalization constant for the vertical kernel.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *convolve(const double *src, pdt width, pdt height, const double *h_kernel, const double *v_kernel, pdt h_support, pdt v_support, double h_norm, double v_norm);

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
double *reconstruct(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, bool nop);

/**
 * \brief Resize an image using a reconstruction filter and an inverse discrete convolution with bandwidth 2.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Reconstruction filter function.
 * \param window Filter function window.
 * \param norm Normalization constant.
 * \param LU LU-decomposition matrix constants.
 * \param m Number of matrix coefficients before convergence.
 * \param nop Boolean flag for a no-op case.
 * \return Destination image in a 4-channel format, or null pointer if OOM.
 */
double *reconstruct_iconvolve_b2(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, const double *LU, pdt m, bool nop);

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
 * \param LU LU-decomposition matrix constants.
 * \param m Number of matrix coefficients before convergence.
 * \param b Inverse discrete convolution bandwidth. Note that the domain in which the inverse convolution is performed must have each dimension be at least 2 * bandwidth - 1.
 * \param nop Boolean flag for a no-op case.
 * \return Destination image in a 4-channel format, or null pointer if OOM or iconvolve error.
 */
double *reconstruct_iconvolve(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, const double *LU, pdt m, pdt b, bool nop);

/**
 * \brief Wrapper for `sample`, `scale`, `reconstruct`, and `reconstruct_iconvolve`.
 * \param src Source image in a 4-channel format.
 * \param src_width Source image width.
 * \param src_height Source image height.
 * \param dst_width Destination image width.
 * \param dst_height Destination image height.
 * \param filter Resizing method (filter) to be used. Defaults to Mitchell-Netravali.
 * \return Destination image in a 4-channel format, or null pointer if OOM or iconvolve error.
 */
double *resize(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, Filter filter);

#endif

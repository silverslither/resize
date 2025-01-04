// Copyright (c) 2024-2025 silverslither.

#include "resize.h"
#include "filters.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static inline double q_fmax(double x, double y) {
    return x > y ? x : y;
}
static inline s32 min(s32 x, s32 y) {
    return x < y ? x : y;
}
static inline s32 max(s32 x, s32 y) {
    return x > y ? x : y;
}

static double *imgcpy(const double *src, s32 width, s32 height) {
    size_t size = ((size_t)width * (size_t)height * sizeof(double)) << 2;
    double *dst = malloc(size);
    if (!dst)
        return NULL;
    memcpy(dst, src, size);
    return dst;
}

double *sample(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;
    if (src_width == dst_width && src_height == dst_height)
        return imgcpy(src, src_width, src_height);

    const double x_factor = (double)src_width / dst_width;
    const double y_factor = (double)src_height / dst_height;

    double *dst = malloc(((size_t)dst_width * (size_t)dst_height * sizeof(double)) << 2);
    if (!dst)
        return NULL;

    s32 dst_pixel = 0;
    for (s32 y = 0; y < dst_height; y++) {
        s32 mapped_y = (s32)ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (s32 x = 0; x < dst_width; x++, dst_pixel += 4) {
            const s32 mapped_x = (s32)ceil(x_factor * (x + 0.5) - 1.0);
            const s32 src_pixel = (mapped_y + mapped_x) << 2;
            dst[dst_pixel + 0] = src[src_pixel + 0];
            dst[dst_pixel + 1] = src[src_pixel + 1];
            dst[dst_pixel + 2] = src[src_pixel + 2];
            dst[dst_pixel + 3] = src[src_pixel + 3];
        }
    }

    return dst;
}

static double *h_filter(const double *src, s32 src_width, s32 height, s32 dst_width, const s32 *bounds, const double *coeffs, size_t support) {
    const s32 adj_src_width = src_width << 2;
    double *dst = calloc(((size_t)dst_width * (size_t)height) << 2, 8);
    if (!dst)
        return NULL;

    const s32 *bounds_ptr;
    const double *coeffs_ptr;
    s32 src_offset = 0;
    s32 dst_pixel = 0;
    for (s32 y = 0; y < height; y++, src_offset += adj_src_width) {
        coeffs_ptr = coeffs;
        bounds_ptr = bounds;

        for (s32 x = 0; x < dst_width; x++, bounds_ptr += 2, coeffs_ptr += support, dst_pixel += 4) {
            const s32 min_x = bounds_ptr[0];
            const s32 max_x = bounds_ptr[1];

            s32 src_pixel = src_offset + (min_x << 2);
            for (s32 s = min_x, i = 0; s <= max_x; s++, i++, src_pixel += 4) {
                const double weight = coeffs_ptr[i];

                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

static double *v_filter(const double *src, s32 width, s32 dst_height, const s32 *bounds, const double *coeffs, size_t support) {
    const s32 adj_width = width << 2;
    double *dst = calloc((size_t)adj_width * (size_t)dst_height, 8);
    if (!dst)
        return NULL;

    const s32 *bounds_ptr = bounds;
    const double *coeffs_ptr = coeffs;
    s32 dst_offset = 0;
    for (s32 y = 0; y < dst_height; y++, bounds_ptr += 2, coeffs_ptr += support, dst_offset += adj_width) {
        const s32 min_y = bounds_ptr[0];
        const s32 max_y = bounds_ptr[1];

        s32 src_pixel = adj_width * min_y;
        for (s32 s = min_y, i = 0; s <= max_y; s++, i++) {
            const double weight = coeffs_ptr[i];

            s32 dst_pixel = dst_offset;
            for (s32 x = 0; x < width; x++, src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

static bool gen_area_filter(s32 **_bounds, double **_coeffs, size_t *_support, s32 src, s32 dst) {
    const double factor = (double)src / dst;
    const double inv_factor = (double)dst / src;
    const size_t support = (size_t)ceil(factor) + 1;

    s32 *bounds = malloc(2 * (size_t)dst * sizeof(s32));
    if (!bounds)
        return true;

    double *coeffs = malloc(support * (size_t)dst * sizeof(double));
    if (!coeffs) {
        free(bounds);
        return true;
    }

    *_bounds = bounds;
    *_coeffs = coeffs;
    *_support = support;

    for (s32 z = 0; z < dst; z++, bounds += 2, coeffs += support) {
        const double min_mapped_z = factor * z;
        const double max_mapped_z = min_mapped_z + factor;
        s32 min_z = (s32)floor(min_mapped_z + 1.0);
        s32 max_z = (s32)ceil(max_mapped_z - 1.0);

        if (max_z < min_z) {
            bounds[1] = bounds[0] = max_z;
            coeffs[0] = 1.0;
            continue;
        }

        bounds[0] = min_z - 1;
        bounds[1] = max_z;
        coeffs[0] = inv_factor * ((double)min_z - min_mapped_z);
        s32 i = 1;
        for (s32 s = min_z; s < max_z; s++, i++)
            coeffs[i] = inv_factor;
        coeffs[i] = inv_factor * (max_mapped_z - (double)max_z);
    }

    return false;
}

static double *h_scale(const double *src, s32 src_width, s32 height, s32 dst_width) {
    s32 *bounds;
    double *coeffs;
    size_t support;
    if (gen_area_filter(&bounds, &coeffs, &support, src_width, dst_width))
        return NULL;

    double *dst = h_filter(src, src_width, height, dst_width, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

static double *v_scale(const double *src, s32 width, s32 src_height, s32 dst_height) {
    s32 *bounds;
    double *coeffs;
    size_t support;
    if (gen_area_filter(&bounds, &coeffs, &support, src_height, dst_height))
        return NULL;

    double *dst = v_filter(src, width, dst_height, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

double *scale(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (dst_width == src_width) {
        if (dst_height == src_height)
            return imgcpy(src, src_width, src_height);
        return v_scale(src, src_width, src_height, dst_height);
    } else if (dst_height == src_height) {
        return h_scale(src, src_width, src_height, dst_width);
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor > y_factor) {
        double *temp = v_scale(src, src_width, src_height, dst_height);
        if (!temp)
            return NULL;

        double *ret = h_scale(temp, src_width, dst_height, dst_width);
        free(temp);
        return ret;
    } else {
        double *temp = h_scale(src, src_width, src_height, dst_width);
        if (!temp)
            return NULL;

        double *ret = v_scale(temp, dst_width, src_height, dst_height);
        free(temp);
        return ret;
    }
}

static bool gen_discrete_filter(s32 **_bounds, double **_coeffs, size_t *_support, s32 src, s32 dst, double (*filter)(double), double window, double norm) {
    const s32 max_s = src - 1;
    const double factor = (double)src / dst;
    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;
    size_t support_val = (size_t)ceil(2.0 * window);

    s32 *bounds = malloc(2 * (size_t)dst * sizeof(s32));
    if (!bounds)
        return true;

    double *coeffs = malloc(support_val * (size_t)dst * sizeof(double));
    if (!coeffs) {
        free(bounds);
        return true;
    }

    *_bounds = bounds;
    *_coeffs = coeffs;
    *_support = support_val;

    for (s32 z = 0; z < dst; z++, bounds += 2, coeffs += support_val) {
        const double mapped_x = factor * (z + 0.5) - 0.5;
        const s32 min_z = max((s32)floor(mapped_x - window + 1.0), 0);
        const s32 max_z = min((s32)ceil(mapped_x + window - 1.0), max_s);
        bounds[0] = min_z;
        bounds[1] = max_z;

        double weight_total = 0.0;
        for (s32 s = min_z, i = 0; s <= max_z; s++, i++) {
            const double weight = filter(fabs(mapped_x - s) * filter_scale);
            coeffs[i] = weight;
            weight_total += weight;
        }

        weight_total = norm / weight_total;
        for (s32 s = min_z, i = 0; s <= max_z; s++, i++)
            coeffs[i] *= weight_total;
    }

    return false;
}

static double *h_reconstruct(const double *src, s32 src_width, s32 height, s32 dst_width, double (*filter)(double), double window, double norm) {
    s32 *bounds;
    double *coeffs;
    size_t support;
    if (gen_discrete_filter(&bounds, &coeffs, &support, src_width, dst_width, filter, window, norm))
        return NULL;

    double *dst = h_filter(src, src_width, height, dst_width, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

static double *v_reconstruct(const double *src, s32 width, s32 src_height, s32 dst_height, double (*filter)(double), double window, double norm) {
    s32 *bounds;
    double *coeffs;
    size_t support;
    if (gen_discrete_filter(&bounds, &coeffs, &support, src_height, dst_height, filter, window, norm))
        return NULL;

    double *dst = v_filter(src, width, dst_height, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

double *reconstruct(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, double (*filter)(double), double window, double norm, int nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return imgcpy(src, src_width, src_height);
            return v_reconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        } else if (dst_height == src_height) {
            return h_reconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        }
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor > y_factor) {
        double *temp = v_reconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        if (!temp)
            return NULL;

        double *ret = h_reconstruct(temp, src_width, dst_height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        double *temp = h_reconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        if (!temp)
            return NULL;

        double *ret = v_reconstruct(temp, dst_width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    }
}

static void h_iconvolve(double *img, s32 width, s32 height, const double *L, int m, double c) {
    int n = m - 1;
    if (m >= width)
        n = m = width - 1;

    const double L_0 = L[0];
    const double L_inf = L[m - 1];
    const double L_inf_mul = L_inf * c;
    const double L_inf_div = L[n] / c;
    const s32 f_adj = (width << 2) + 8;
    double *f = img + 4;

    for (s32 y = 0; y < height; y++) {
        s32 x;

        for (x = 1; x < m; x++, f += 4) {
            const double L_x = L[x - 1];
            f[0] -= L_x * f[-4];
            f[1] -= L_x * f[-3];
            f[2] -= L_x * f[-2];
            f[3] -= L_x * f[-1];
        }

        for (; x < width - 1; x++, f += 4) {
            f[0] -= L_inf * f[-4];
            f[1] -= L_inf * f[-3];
            f[2] -= L_inf * f[-2];
            f[3] -= L_inf * f[-1];
        }

        f[0] = L_inf_div * (f[0] - L_inf_mul * f[-4]);
        f[1] = L_inf_div * (f[1] - L_inf_mul * f[-3]);
        f[2] = L_inf_div * (f[2] - L_inf_mul * f[-2]);
        f[3] = L_inf_div * (f[3] - L_inf_mul * f[-1]);
        f -= 4;

        for (x = width - 2; x >= m; x--, f -= 4) {
            f[0] = L_inf * (f[0] - f[4]);
            f[1] = L_inf * (f[1] - f[5]);
            f[2] = L_inf * (f[2] - f[6]);
            f[3] = L_inf * (f[3] - f[7]);
        }

        for (; x > 0; x--, f -= 4) {
            const double L_x = L[x];
            f[0] = L_x * (f[0] - f[4]);
            f[1] = L_x * (f[1] - f[5]);
            f[2] = L_x * (f[2] - f[6]);
            f[3] = L_x * (f[3] - f[7]);
        }

        f[0] = L_0 * (f[0] - c * f[4]);
        f[1] = L_0 * (f[1] - c * f[5]);
        f[2] = L_0 * (f[2] - c * f[6]);
        f[3] = L_0 * (f[3] - c * f[7]);

        f += f_adj - 4;
    }
}

static void v_iconvolve(double *img, s32 width, s32 height, const double *L, int m, double c) {
    int n = m - 1;
    if (m >= height)
        n = m = height - 1;

    const double L_0 = L[0];
    const double L_inf = L[m - 1];
    const double L_inf_mul = L_inf * c;
    const double L_inf_div = L[n] / c;
    const s32 adj_width = width << 2;
    double *pf = img;
    double *f = pf + adj_width;

    for (s32 x = 0; x < width; x++) {
        s32 y;

        for (y = 1; y < m; y++, pf = f, f += adj_width) {
            const double L_y = L[y - 1];
            f[0] -= L_y * pf[0];
            f[1] -= L_y * pf[1];
            f[2] -= L_y * pf[2];
            f[3] -= L_y * pf[3];
        }

        for (; y < height - 1; y++, pf = f, f += adj_width) {
            f[0] -= L_inf * pf[0];
            f[1] -= L_inf * pf[1];
            f[2] -= L_inf * pf[2];
            f[3] -= L_inf * pf[3];
        }

        f[0] = L_inf_div * (f[0] - L_inf_mul * pf[0]);
        f[1] = L_inf_div * (f[1] - L_inf_mul * pf[1]);
        f[2] = L_inf_div * (f[2] - L_inf_mul * pf[2]);
        f[3] = L_inf_div * (f[3] - L_inf_mul * pf[3]);
        pf = f;
        f = pf - adj_width;

        for (y = height - 2; y >= m; y--, pf = f, f -= adj_width) {
            f[0] = L_inf * (f[0] - pf[0]);
            f[1] = L_inf * (f[1] - pf[1]);
            f[2] = L_inf * (f[2] - pf[2]);
            f[3] = L_inf * (f[3] - pf[3]);
        }

        for (; y > 0; y--, pf = f, f -= adj_width) {
            const double L_y = L[y];
            f[0] = L_y * (f[0] - pf[0]);
            f[1] = L_y * (f[1] - pf[1]);
            f[2] = L_y * (f[2] - pf[2]);
            f[3] = L_y * (f[3] - pf[3]);
        }

        f[0] = L_0 * (f[0] - c * pf[0]);
        f[1] = L_0 * (f[1] - c * pf[1]);
        f[2] = L_0 * (f[2] - c * pf[2]);
        f[3] = L_0 * (f[3] - c * pf[3]);

        pf = f + 4;
        f = pf + adj_width;
    }
}

static double *h_reconstruct_iconvolve(const double *src, s32 src_width, s32 height, s32 dst_width, double (*filter)(double), double window, double norm, const double *L, int m, double c) {
    if (dst_width > src_width) {
        if (src_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *temp = imgcpy(src, src_width, height);
        if (!temp)
            return NULL;
        h_iconvolve(temp, src_width, height, L, m, c);

        double *ret = h_reconstruct(temp, src_width, height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *ret = h_reconstruct(src, src_width, height, dst_width, filter, window, norm);
        if (!ret)
            return NULL;

        h_iconvolve(ret, dst_width, height, L, m, c);
        return ret;
    }
}

static double *v_reconstruct_iconvolve(const double *src, s32 width, s32 src_height, s32 dst_height, double (*filter)(double), double window, double norm, const double *L, int m, double c) {
    if (dst_height > src_height) {
        if (src_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *temp = imgcpy(src, width, src_height);
        if (!temp)
            return NULL;
        v_iconvolve(temp, width, src_height, L, m, c);

        double *ret = v_reconstruct(temp, width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *ret = v_reconstruct(src, width, src_height, dst_height, filter, window, norm);
        if (!ret)
            return NULL;

        v_iconvolve(ret, width, dst_height, L, m, c);
        return ret;
    }
}

double *reconstruct_iconvolve(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, double (*filter)(double), double window, double norm, const double *L, int m, double c, int nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return imgcpy(src, src_width, src_height);
            return v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        } else if (dst_height == src_height) {
            return h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        }
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor > y_factor) {
        double *temp = v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        if (!temp)
            return NULL;

        double *ret = h_reconstruct_iconvolve(temp, src_width, dst_height, dst_width, filter, window, norm, L, m, c);
        free(temp);
        return ret;
    } else {
        double *temp = h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        if (!temp)
            return NULL;

        double *ret = v_reconstruct_iconvolve(temp, dst_width, src_height, dst_height, filter, window, norm, L, m, c);
        free(temp);
        return ret;
    }
}

double *resize(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, Filter filter) {
    switch (filter) {
    case NEAREST:
        return sample(src, src_width, src_height, dst_width, dst_height);
    case AREA:
        return scale(src, src_width, src_height, dst_width, dst_height);
    case TRIANGLE:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Triangle, 1.0, 1.0, 1);
    case HERMITE:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hermite, 1.0, 1.0, 1);
    case BSPLINE2:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, BSpline2, 1.5, 1.0, 0);
    case BSPLINE3:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, BSpline3, 2.0, 1.0, 0);
    default:
    case MITNET:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, MitNet, 2.0, 1.0, 0);
    case CATROM:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, CatRom, 2.0, 1.0, 1);
    case HAMMING3:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hamming3, 3.0, 1.0, 1);
    case HAMMING4:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hamming4, 4.0, 1.0, 1);
    case BSPLINE2I:
        return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, BSpline2, 1.5, 8.0, L_bspline2i, 11, 1.1428571428571428, 1);
    case BSPLINE3I:
        return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, BSpline3, 2.0, 6.0, L_bspline3i, 14, 1.2, 1);
    case OMOMS3I:
        return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, OMOMS3, 2.0, 5.25, L_omoms3, 18, 1.2352941176470589, 1);
    }
}

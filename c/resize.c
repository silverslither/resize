// Copyright (c) 2024 silverslither.

#include "resize.h"
#include "filters.h"
#include <math.h>
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
    size_t size = ((size_t)width * (size_t)height) << 5;
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

    double *dst = malloc(((size_t)dst_width * (size_t)dst_height) << 5);
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

static double *hscale(const double *src, s32 src_width, s32 height, s32 dst_width) {
    const s32 max_high_x = src_width - 1;
    const s32 adj_src_width = src_width << 2;
    const s32 adj_dst_width = dst_width << 2;
    const s32 adj_src_area = adj_src_width * height;
    const double factor = (double)src_width / dst_width;
    const double inv_factor = (double)dst_width / src_width;

    double *dst = calloc((size_t)height * (size_t)adj_dst_width, 8);
    if (!dst)
        return NULL;

    s32 dst_offset = 0;
    for (s32 x = 0; x < dst_width; x++, dst_offset += 4) {
        const double min_mapped_x = factor * x;
        const double max_mapped_x = min_mapped_x + factor;
        s32 min_x = (s32)ceil(min_mapped_x);
        s32 max_x = (s32)floor(max_mapped_x);

        s32 dst_pixel = dst_offset;

        if (max_x < min_x) {
            for (s32 src_pixel = max_x << 2; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] = src[src_pixel + 0];
                dst[dst_pixel + 1] = src[src_pixel + 1];
                dst[dst_pixel + 2] = src[src_pixel + 2];
                dst[dst_pixel + 3] = src[src_pixel + 3];
            }
            continue;
        }

        const s32 low_x = max(min_x - 1, 0) << 2;
        const s32 high_x = min(max_x, max_high_x) << 2;
        const double low_mult = (double)min_x - min_mapped_x;
        const double high_mult = max_mapped_x - (double)max_x;
        min_x <<= 2;
        max_x <<= 2;

        for (s32 y = 0; y < adj_src_area; y += adj_src_width, dst_pixel += adj_dst_width) {
            s32 src_pixel = y + low_x;
            dst[dst_pixel + 0] += src[src_pixel + 0] * low_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * low_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * low_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * low_mult;

            const s32 src_max = y + max_x;
            for (src_pixel = y + min_x; src_pixel < src_max; src_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }

            src_pixel = y + high_x;
            dst[dst_pixel + 0] += src[src_pixel + 0] * high_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * high_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * high_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * high_mult;
            dst[dst_pixel + 0] *= inv_factor;
            dst[dst_pixel + 1] *= inv_factor;
            dst[dst_pixel + 2] *= inv_factor;
            dst[dst_pixel + 3] *= inv_factor;
        }
    }

    return dst;
}

static double *vscale(const double *src, s32 width, s32 src_height, s32 dst_height) {
    const s32 max_high_y = src_height - 1;
    const s32 adj_width = width << 2;
    const double factor = (double)src_height / dst_height;
    const double inv_factor = (double)dst_height / src_height;

    double *dst = calloc((size_t)adj_width * (size_t)dst_height, 8);
    if (!dst)
        return NULL;

    s32 dst_pixel = 0;
    for (s32 y = 0; y < dst_height; y++) {
        const double min_mapped_y = factor * y;
        const double max_mapped_y = min_mapped_y + factor;
        s32 min_y = (s32)ceil(min_mapped_y);
        s32 max_y = (s32)floor(max_mapped_y);

        if (max_y < min_y) {
            const s32 src_offset = adj_width * max_y;
            const s32 src_max = src_offset + adj_width;
            for (s32 src_pixel = src_offset; src_pixel < src_max; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] = src[src_pixel + 0];
                dst[dst_pixel + 1] = src[src_pixel + 1];
                dst[dst_pixel + 2] = src[src_pixel + 2];
                dst[dst_pixel + 3] = src[src_pixel + 3];
            }
            continue;
        }

        const s32 low_y = adj_width * max(min_y - 1, 0);
        const s32 high_y = adj_width * min(max_y, max_high_y);
        const double low_mult = (double)min_y - min_mapped_y;
        const double high_mult = max_mapped_y - (double)max_y;
        min_y *= adj_width;
        max_y *= adj_width;

        for (s32 x = 0; x < adj_width; x += 4, dst_pixel += 4) {
            s32 src_pixel = x + low_y;
            dst[dst_pixel + 0] += src[src_pixel + 0] * low_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * low_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * low_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * low_mult;

            const s32 src_max = x + max_y;
            for (src_pixel = x + min_y; src_pixel < src_max; src_pixel += adj_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }

            src_pixel = x + high_y;
            dst[dst_pixel + 0] += src[src_pixel + 0] * high_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * high_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * high_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * high_mult;
            dst[dst_pixel + 0] *= inv_factor;
            dst[dst_pixel + 1] *= inv_factor;
            dst[dst_pixel + 2] *= inv_factor;
            dst[dst_pixel + 3] *= inv_factor;
        }
    }

    return dst;
}

double *scale(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (dst_width == src_width) {
        if (dst_height == src_height)
            return imgcpy(src, src_width, src_height);
        return vscale(src, src_width, src_height, dst_height);
    } else if (dst_height == src_height) {
        return hscale(src, src_width, src_height, dst_width);
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor >= y_factor) {
        double *temp = vscale(src, src_width, src_height, dst_height);
        if (!temp)
            return NULL;

        double *ret = hscale(temp, src_width, dst_height, dst_width);
        free(temp);
        return ret;
    } else {
        double *temp = hscale(src, src_width, src_height, dst_width);
        if (!temp)
            return NULL;

        double *ret = vscale(temp, dst_width, src_height, dst_height);
        free(temp);
        return ret;
    }
}

static double *hreconstruct(const double *src, s32 src_width, s32 height, s32 dst_width, double (*filter)(double), double window, double norm) {
    const s32 max_s = src_width - 1;
    const s32 adj_src_width = src_width << 2;
    const s32 adj_dst_width = dst_width << 2;
    const s32 adj_dst_area = adj_dst_width * height;
    const double factor = (double)src_width / dst_width;

    double *dst = calloc((size_t)height * (size_t)adj_dst_width, 8);
    if (!dst)
        return NULL;

    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    s32 dst_offset = 0;
    for (s32 x = 0; x < dst_width; x++, dst_offset += 4) {
        const double mapped_x = factor * (x + 0.5) - 0.5;
        const s32 min_x = max((s32)ceil(mapped_x - window), 0);
        const s32 max_x = min((s32)floor(mapped_x + window), max_s);
        double weight_total = 0;

        s32 src_offset = min_x << 2;
        for (s32 s = min_x; s <= max_x; s++, src_offset += 4) {
            const double weight = filter(fabs(mapped_x - s) * filter_scale);
            weight_total += weight;

            for (s32 src_pixel = src_offset, dst_pixel = dst_offset; dst_pixel < adj_dst_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        weight_total = norm / weight_total;
        for (s32 dst_pixel = dst_offset; dst_pixel < adj_dst_area; dst_pixel += adj_dst_width) {
            dst[dst_pixel + 0] *= weight_total;
            dst[dst_pixel + 1] *= weight_total;
            dst[dst_pixel + 2] *= weight_total;
            dst[dst_pixel + 3] *= weight_total;
        }
    }

    return dst;
}

static double *vreconstruct(const double *src, s32 width, s32 src_height, s32 dst_height, double (*filter)(double), double window, double norm) {
    const s32 max_s = src_height - 1;
    const s32 adj_width = width << 2;
    const double factor = (double)src_height / dst_height;

    double *dst = calloc((size_t)adj_width * (size_t)dst_height, 8);
    if (!dst)
        return NULL;

    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    s32 dst_offset = 0;
    s32 ndst_offset = adj_width;
    for (s32 y = 0; y < dst_height; y++, dst_offset = ndst_offset, ndst_offset += adj_width) {
        const double mapped_y = factor * (y + 0.5) - 0.5;
        const s32 min_y = max((s32)ceil(mapped_y - window), 0);
        const s32 max_y = min((s32)floor(mapped_y + window), max_s);
        double weight_total = 0;

        s32 src_offset = adj_width * min_y;
        for (s32 s = min_y; s <= max_y; s++, src_offset += adj_width) {
            const double weight = filter(fabs(mapped_y - s) * filter_scale);
            weight_total += weight;

            for (s32 src_pixel = src_offset, dst_pixel = dst_offset; dst_pixel < ndst_offset; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        weight_total = norm / weight_total;
        for (s32 dst_pixel = dst_offset; dst_pixel < ndst_offset; dst_pixel += 4) {
            dst[dst_pixel + 0] *= weight_total;
            dst[dst_pixel + 1] *= weight_total;
            dst[dst_pixel + 2] *= weight_total;
            dst[dst_pixel + 3] *= weight_total;
        }
    }

    return dst;
}

double *reconstruct(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, double (*filter)(double), double window, double norm, int nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return imgcpy(src, src_width, src_height);
            return vreconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        } else if (dst_height == src_height) {
            return hreconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        }
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor >= y_factor) {
        double *temp = vreconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        if (!temp)
            return NULL;

        double *ret = hreconstruct(temp, src_width, dst_height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        double *temp = hreconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        if (!temp)
            return NULL;

        double *ret = vreconstruct(temp, dst_width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    }
}

static void hiconvolve(double *img, s32 width, s32 height, const double *L, int m, double c) {
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

static void viconvolve(double *img, s32 width, s32 height, const double *L, int m, double c) {
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

static double *hreconstruct_iconvolve(const double *src, s32 src_width, s32 height, s32 dst_width, double (*filter)(double), double window, double norm, const double *L, int m, double c) {
    if (dst_width > src_width) {
        if (src_width == 1)
            return hreconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *temp = imgcpy(src, src_width, height);
        if (!temp)
            return NULL;
        hiconvolve(temp, src_width, height, L, m, c);

        double *ret = hreconstruct(temp, src_width, height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_width == 1)
            return hreconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *ret = hreconstruct(src, src_width, height, dst_width, filter, window, norm);
        if (!ret)
            return NULL;

        hiconvolve(ret, dst_width, height, L, m, c);
        return ret;
    }
}

static double *vreconstruct_iconvolve(const double *src, s32 width, s32 src_height, s32 dst_height, double (*filter)(double), double window, double norm, const double *L, int m, double c) {
    if (dst_height > src_height) {
        if (src_height == 1)
            return vreconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *temp = imgcpy(src, width, src_height);
        if (!temp)
            return NULL;
        viconvolve(temp, width, src_height, L, m, c);

        double *ret = vreconstruct(temp, width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_height == 1)
            return vreconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *ret = vreconstruct(src, width, src_height, dst_height, filter, window, norm);
        if (!ret)
            return NULL;

        viconvolve(ret, width, dst_height, L, m, c);
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
            return vreconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        } else if (dst_height == src_height) {
            return hreconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        }
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor >= y_factor) {
        double *temp = vreconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        if (!temp)
            return NULL;

        double *ret = hreconstruct_iconvolve(temp, src_width, dst_height, dst_width, filter, window, norm, L, m, c);
        free(temp);
        return ret;
    } else {
        double *temp = hreconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        if (!temp)
            return NULL;

        double *ret = vreconstruct_iconvolve(temp, dst_width, src_height, dst_height, filter, window, norm, L, m, c);
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

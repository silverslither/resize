// Copyright (c) 2024 silverslither.

#include "resize.h"
#include "filters.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static inline s32 q_floor(double x) {
    return (s32)x;
}

static inline s32 q_ceil(double x) {
    const s32 r = (s32)x;
    return r + (x > (double)r);
}

static inline double q_fmax(double x, double y) {
    return x > y ? x : y;
}
static inline s32 min(s32 x, s32 y) {
    return x < y ? x : y;
}
static inline s32 max(s32 x, s32 y) {
    return x > y ? x : y;
}

double *sample(const double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height) {
    if (dst_width <= 0 || dst_height <= 0)
        return NULL;

    const double x_factor = (double)src_width / dst_width;
    const double y_factor = (double)src_height / dst_height;

    double *dst = calloc((dst_width * dst_height) << 2, 8);
    if (!dst)
        return NULL;

    s32 dst_pixel = 0;
    for (s32 y = 0; y < dst_height; y++) {
        s32 mapped_y = q_ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (s32 x = 0; x < dst_width; x++, dst_pixel += 4) {
            const s32 mapped_x = q_ceil(x_factor * (x + 0.5) - 1.0);
            const s32 src_pixel = (mapped_y + mapped_x) << 2;
            dst[dst_pixel + 0] += src[src_pixel + 0];
            dst[dst_pixel + 1] += src[src_pixel + 1];
            dst[dst_pixel + 2] += src[src_pixel + 2];
            dst[dst_pixel + 3] += src[src_pixel + 3];
        }
    }

    return dst;
}

static double *hscale(const double *src, s32 src_width, s32 height, s32 dst_width) {
    const s32 adj_src_width = src_width << 2;
    const s32 adj_dst_width = dst_width << 2;
    const s32 adj_src_area = adj_src_width * height;
    const double factor = (double)src_width / dst_width;
    const double inv_factor = (double)dst_width / src_width;

    double *dst = calloc(height * adj_dst_width, 8);
    if (!dst)
        return NULL;

    s32 dst_offset = 0;
    for (s32 x = 0; x < dst_width; x++, dst_offset += 4) {
        const double min_mapped_x = factor * x;
        const double max_mapped_x = min_mapped_x + factor;
        s32 min_x = q_ceil(min_mapped_x);
        s32 max_x = q_floor(max_mapped_x);

        s32 dst_pixel = dst_offset;

        if (max_x < min_x) {
            for (s32 src_pixel = max_x << 2; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }
            continue;
        }

        const s32 low_x = max(min_x - 1, 0) << 2;
        const s32 high_x = min(max_x, src_width - 1) << 2;
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
            for (s32 src_pixel = y + min_x; src_pixel < src_max; src_pixel += 4) {
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
    const s32 adj_width = width << 2;
    const double factor = (double)src_height / dst_height;
    const double inv_factor = (double)dst_height / src_height;

    double *dst = calloc(adj_width * dst_height, 8);
    if (!dst)
        return NULL;

    s32 dst_pixel = 0;
    for (s32 y = 0; y < dst_height; y++) {
        const double min_mapped_y = factor * y;
        const double max_mapped_y = min_mapped_y + factor;
        s32 min_y = q_ceil(min_mapped_y);
        s32 max_y = q_floor(max_mapped_y);

        if (max_y < min_y) {
            const s32 src_offset = adj_width * max_y;
            const s32 src_max = src_offset + adj_width;
            for (s32 src_pixel = src_offset; src_pixel < src_max; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }
            continue;
        }

        const s32 low_y = adj_width * max(min_y - 1, 0);
        const s32 high_y = adj_width * min(max_y, src_height - 1);
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
            for (s32 src_pixel = x + min_y; src_pixel < src_max; src_pixel += adj_width) {
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

double *scale(double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height) {
    if (dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (dst_width == src_width) {
        if (dst_height == src_height)
            return src;
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

static double *hfilter(const double *src, s32 src_width, s32 height, s32 dst_width, double (*filter)(double), double window) {
    const s32 adj_src_width = src_width << 2;
    const s32 adj_dst_width = dst_width << 2;
    const s32 adj_src_area = adj_src_width * height;
    const double factor = (double)src_width / dst_width;

    double *dst = calloc(height * adj_dst_width, 8);
    if (!dst)
        return NULL;

    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    const s32 width_end = adj_src_width - 4;

    s32 dst_offset = 0;
    for (s32 x = 0; x < dst_width; x++, dst_offset += 4) {
        const double mapped_x = factor * (x + 0.5) - 0.5;
        s32 min_x = q_ceil(mapped_x - window);
        s32 max_x = q_floor(mapped_x + window);

        while (min_x < 0) {
            const double weight = filter((mapped_x - (min_x++)) * filter_scale) * filter_scale;
            for (s32 src_pixel = 0, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        while (max_x >= src_width) {
            const double weight = filter(((max_x--) - mapped_x) * filter_scale) * filter_scale;
            for (s32 src_pixel = width_end, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        s32 src_offset = min_x << 2;
        for (s32 s = min_x; s <= max_x; s++, src_offset += 4) {
            const double weight = filter(fabs(mapped_x - s) * filter_scale) * filter_scale;
            for (s32 src_pixel = src_offset, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

static double *vfilter(const double *src, s32 width, s32 src_height, s32 dst_height, double (*filter)(double), double window) {
    const s32 adj_width = width << 2;
    const s32 adj_src_area = adj_width * src_height;
    const double factor = (double)src_height / dst_height;

    double *dst = calloc(adj_width * dst_height, 8);
    if (!dst)
        return NULL;

    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    const s32 height_end = adj_width * (src_height - 1);

    s32 dst_offset = 0;
    for (s32 y = 0; y < dst_height; y++, dst_offset += adj_width) {
        const double mapped_y = factor * (y + 0.5) - 0.5;
        s32 min_y = q_ceil(mapped_y - window);
        s32 max_y = q_floor(mapped_y + window);

        while (min_y < 0) {
            const double weight = filter((mapped_y - (min_y++)) * filter_scale) * filter_scale;
            for (s32 src_pixel = 0, dst_pixel = dst_offset; src_pixel < adj_width; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        while (max_y >= src_height) {
            const double weight = filter(((max_y--) - mapped_y) * filter_scale) * filter_scale;
            for (s32 src_pixel = height_end, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        s32 src_offset = adj_width * min_y;
        s32 src_max = src_offset + adj_width;
        for (s32 s = min_y; s <= max_y; s++, src_offset = src_max, src_max += adj_width) {
            const double weight = filter(fabs(mapped_y - s) * filter_scale) * filter_scale;
            for (s32 src_pixel = src_offset, dst_pixel = dst_offset; src_pixel < src_max; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

double *resize(double *src, s32 src_width, s32 src_height, s32 dst_width, s32 dst_height, Filter filter) {
    if (dst_width <= 0 || dst_height <= 0)
        return NULL;

    double (*filter_func)(double);
    double window;
    int nop = 0;

    switch (filter) {
    case NEAREST:
        return sample(src, src_width, src_height, dst_width, dst_height);
    case AREA:
        return scale(src, src_width, src_height, dst_width, dst_height);
    case TRIANGLE:
        filter_func = Triangle;
        window = 1.0;
        nop = 1;
        break;
    case HERMITE:
        filter_func = Hermite;
        window = 1.0;
        nop = 1;
        break;
    case B_SPLINE_2:
        filter_func = BSpline2;
        window = 1.5;
        break;
    case B_SPLINE_3:
        filter_func = BSpline3;
        window = 2.0;
        break;
    case KEYS_HALF:
        filter_func = KeysHalf;
        window = 2.0;
        break;
    default:
    case MITNET:
        filter_func = MitNet;
        window = 2.0;
        break;
    case MITNET_SHARP:
        filter_func = MitNetSharp;
        window = 2.0;
        break;
    case CATROM:
        filter_func = CatRom;
        window = 2.0;
        nop = 1;
        break;
    case CATROM_SHARP:
        filter_func = CatRomSharp;
        window = 2.0;
        nop = 1;
        break;
    case LANCZOS_3:
        filter_func = Lanczos3;
        window = 3.0;
        nop = 1;
        break;
    case LANCZOS_4:
        filter_func = Lanczos4;
        window = 4.0;
        nop = 1;
        break;
    }

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return src;
            return vfilter(src, src_width, src_height, dst_height, filter_func, window);
        } else if (dst_height == src_height) {
            return hfilter(src, src_width, src_height, dst_width, filter_func, window);
        }
    }

    const double x_factor = (double)dst_width / src_width;
    const double y_factor = (double)dst_height / src_height;

    if (x_factor >= y_factor) {
        double *temp = vfilter(src, src_width, src_height, dst_height, filter_func, window);
        if (!temp)
            return NULL;
        double *ret = hfilter(temp, src_width, dst_height, dst_width, filter_func, window);
        free(temp);
        return ret;
    } else {
        double *temp = hfilter(src, src_width, src_height, dst_width, filter_func, window);
        if (!temp)
            return NULL;
        double *ret = vfilter(temp, dst_width, src_height, dst_height, filter_func, window);
        free(temp);
        return ret;
    }
}

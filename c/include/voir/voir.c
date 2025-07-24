// Copyright (c) 2024-2025 silverslither.

#include "voir.h"
#include "filters.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static inline double q_fmax(double x, double y) {
    return x > y ? x : y;
}
static inline pdt min(pdt x, pdt y) {
    return x < y ? x : y;
}
static inline pdt max(pdt x, pdt y) {
    return x > y ? x : y;
}

static double *imgcpy(const double *src, pdt width, pdt height) {
    double *dst;
    pdt size = (width * height * sizeof(*dst)) << 2;
    dst = malloc(size);
    if (!dst)
        return NULL;
    memcpy(dst, src, size);
    return dst;
}

double *sample(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;
    if (src_width == dst_width && src_height == dst_height)
        return imgcpy(src, src_width, src_height);

    const double x_factor = (double)src_width / dst_width;
    const double y_factor = (double)src_height / dst_height;

    double *dst = malloc((dst_width * dst_height * sizeof(*dst)) << 2);
    if (!dst)
        return NULL;

    pdt dst_pixel = 0;
    for (pdt y = 0; y < dst_height; y++) {
        pdt mapped_y = (pdt)ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (pdt x = 0; x < dst_width; x++, dst_pixel += 4) {
            const pdt mapped_x = (pdt)ceil(x_factor * (x + 0.5) - 1.0);
            const pdt src_pixel = (mapped_y + mapped_x) << 2;
            dst[dst_pixel + 0] = src[src_pixel + 0];
            dst[dst_pixel + 1] = src[src_pixel + 1];
            dst[dst_pixel + 2] = src[src_pixel + 2];
            dst[dst_pixel + 3] = src[src_pixel + 3];
        }
    }

    return dst;
}

static double *h_filter(const double *src, pdt src_width, pdt height, pdt dst_width, const pdt *bounds, const double *coeffs, pdt support) {
    const pdt adj_src_width = src_width << 2;
    double *dst = calloc((dst_width * height) << 2, sizeof(*dst));
    if (!dst)
        return NULL;

    const pdt *bounds_ptr;
    const double *coeffs_ptr;
    pdt src_offset = 0;
    pdt dst_pixel = 0;
    for (pdt y = 0; y < height; y++, src_offset += adj_src_width) {
        coeffs_ptr = coeffs;
        bounds_ptr = bounds;

        for (pdt x = 0; x < dst_width; x++, bounds_ptr += 2, coeffs_ptr += support, dst_pixel += 4) {
            const pdt min_x = bounds_ptr[0];
            const pdt max_x = bounds_ptr[1];

            pdt src_pixel = src_offset + (min_x << 2);
            for (pdt s = min_x, i = 0; s <= max_x; s++, i++, src_pixel += 4) {
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

static double *v_filter(const double *src, pdt width, pdt dst_height, const pdt *bounds, const double *coeffs, pdt support) {
    const pdt adj_width = width << 2;
    double *dst = calloc(adj_width * dst_height, sizeof(*dst));
    if (!dst)
        return NULL;

    const pdt *bounds_ptr = bounds;
    const double *coeffs_ptr = coeffs;
    pdt dst_offset = 0;
    for (pdt y = 0; y < dst_height; y++, bounds_ptr += 2, coeffs_ptr += support, dst_offset += adj_width) {
        const pdt min_y = bounds_ptr[0];
        const pdt max_y = bounds_ptr[1];

        pdt src_pixel = adj_width * min_y;
        for (pdt s = min_y, i = 0; s <= max_y; s++, i++) {
            const double weight = coeffs_ptr[i];

            pdt dst_pixel = dst_offset;
            for (pdt x = 0; x < width; x++, src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

static bool gen_area_filter(pdt **_bounds, double **_coeffs, pdt *_support, pdt src, pdt dst) {
    const double factor = (double)src / dst;
    const double inv_factor = (double)dst / src;
    const pdt support = (pdt)ceil(factor) + 1;

    pdt *bounds = malloc(2 * dst * sizeof(*bounds));
    if (!bounds)
        return true;

    double *coeffs = malloc(support * dst * sizeof(*coeffs));
    if (!coeffs) {
        free(bounds);
        return true;
    }

    *_bounds = bounds;
    *_coeffs = coeffs;
    *_support = support;

    for (pdt z = 0; z < dst; z++, bounds += 2, coeffs += support) {
        const double min_mapped_z = factor * z;
        const double max_mapped_z = min_mapped_z + factor;
        const pdt min_z = (pdt)floor(min_mapped_z + 1.0);
        const pdt max_z = (pdt)ceil(max_mapped_z - 1.0);

        if (max_z < min_z) {
            bounds[1] = bounds[0] = max_z;
            coeffs[0] = 1.0;
            continue;
        }

        bounds[0] = min_z - 1;
        bounds[1] = max_z;
        coeffs[0] = inv_factor * ((double)min_z - min_mapped_z);
        pdt i = 1;
        for (pdt s = min_z; s < max_z; s++, i++)
            coeffs[i] = inv_factor;
        coeffs[i] = inv_factor * (max_mapped_z - (double)max_z);
    }

    return false;
}

static double *h_scale(const double *src, pdt src_width, pdt height, pdt dst_width) {
    pdt *bounds;
    double *coeffs;
    pdt support;
    if (gen_area_filter(&bounds, &coeffs, &support, src_width, dst_width))
        return NULL;

    double *dst = h_filter(src, src_width, height, dst_width, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

static double *v_scale(const double *src, pdt width, pdt src_height, pdt dst_height) {
    pdt *bounds;
    double *coeffs;
    pdt support;
    if (gen_area_filter(&bounds, &coeffs, &support, src_height, dst_height))
        return NULL;

    double *dst = v_filter(src, width, dst_height, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

double *scale(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (dst_width == src_width) {
        if (dst_height == src_height)
            return imgcpy(src, src_width, src_height);
        return v_scale(src, src_width, src_height, dst_height);
    } else if (dst_height == src_height) {
        return h_scale(src, src_width, src_height, dst_width);
    }

    const pdt h_intermediate = dst_width * src_height;
    const pdt v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
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

static bool gen_kernel_filter(pdt **_bounds, double **_coeffs, pdt len, const double *kernel, pdt support, double norm) {
    const pdt offset = support >> 1;
    pdt *bounds = malloc(2 * len * sizeof(*bounds));
    if (!bounds)
        return true;

    double *coeffs = malloc(support * len * sizeof(*coeffs));
    if (!coeffs) {
        free(bounds);
        return true;
    }

    *_bounds = bounds;
    *_coeffs = coeffs;

    for (pdt z = 0; z < len; z++, bounds += 2, coeffs += support) {
        const pdt min_z = z - offset;
        const pdt max_z = min(z + offset, len - 1);
        bounds[0] = max(min_z, 0);
        bounds[1] = max_z;

        double weight_total = 0.0;
        for (pdt s = min_z, i = 0; s <= max_z; s++, i++) {
            if (s < 0) {
                i--;
                continue;
            }
            const double weight = kernel[i];
            coeffs[i] = weight;
            weight_total += weight;
        }

        if (fabs(weight_total) > 2.3283064365386963e-10) {
            weight_total = norm / weight_total;
            for (pdt s = min_z, i = 0; s <= max_z; s++, i++)
                coeffs[i] *= weight_total;
        }
    }

    return false;
}

static double *h_convolve(const double *src, pdt width, pdt height, const double *kernel, pdt support, double norm) {
    pdt *bounds;
    double *coeffs;
    if (gen_kernel_filter(&bounds, &coeffs, width, kernel, support, norm))
        return NULL;

    double *dst = h_filter(src, width, height, width, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

static double *v_convolve(const double *src, pdt width, pdt height, const double *kernel, pdt support, double norm) {
    pdt *bounds;
    double *coeffs;
    if (gen_kernel_filter(&bounds, &coeffs, height, kernel, support, norm))
        return NULL;

    double *dst = v_filter(src, width, height, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

double *convolve(const double *src, pdt width, pdt height, const double *h_kernel, const double *v_kernel, pdt h_support, pdt v_support, double h_norm, double v_norm) {
    if (width <= 0 || height <= 0)
        return NULL;

    if (!h_kernel) {
        if (!v_kernel)
            return imgcpy(src, width, height);
        return v_convolve(src, width, height, v_kernel, v_support, v_norm);
    } else if (!v_kernel) {
        return h_convolve(src, width, height, h_kernel, h_support, h_norm);
    }

    double *temp = h_convolve(src, width, height, h_kernel, h_support, h_norm);
    if (!temp)
        return NULL;

    double *ret = v_convolve(temp, width, height, v_kernel, v_support, v_norm);
    free(temp);
    return ret;
}

static bool gen_discrete_filter(pdt **_bounds, double **_coeffs, pdt *_support, pdt src, pdt dst, double (*filter)(double), double window, double norm) {
    const pdt max_s = src - 1;
    const double factor = (double)src / dst;
    const double inv_filter_scale = q_fmax(factor, 1.0);
    const double filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;
    const pdt support = (pdt)ceil(2.0 * window);

    pdt *bounds = malloc(2 * dst * sizeof(*bounds));
    if (!bounds)
        return true;

    double *coeffs = malloc(support * dst * sizeof(*coeffs));
    if (!coeffs) {
        free(bounds);
        return true;
    }

    *_bounds = bounds;
    *_coeffs = coeffs;
    *_support = support;

    for (pdt z = 0; z < dst; z++, bounds += 2, coeffs += support) {
        const double mapped_z = factor * (z + 0.5) - 0.5;
        const pdt min_z = max((pdt)floor(mapped_z - window + 1.0), 0);
        const pdt max_z = min((pdt)ceil(mapped_z + window - 1.0), max_s);
        bounds[0] = min_z;
        bounds[1] = max_z;

        double weight_total = 0.0;
        for (pdt s = min_z, i = 0; s <= max_z; s++, i++) {
            const double weight = filter(fabs(mapped_z - s) * filter_scale);
            coeffs[i] = weight;
            weight_total += weight;
        }

        if (fabs(weight_total) > 2.3283064365386963e-10) {
            weight_total = norm / weight_total;
            for (pdt s = min_z, i = 0; s <= max_z; s++, i++)
                coeffs[i] *= weight_total;
        }
    }

    return false;
}

static double *h_reconstruct(const double *src, pdt src_width, pdt height, pdt dst_width, double (*filter)(double), double window, double norm) {
    pdt *bounds;
    double *coeffs;
    pdt support;
    if (gen_discrete_filter(&bounds, &coeffs, &support, src_width, dst_width, filter, window, norm))
        return NULL;

    double *dst = h_filter(src, src_width, height, dst_width, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

static double *v_reconstruct(const double *src, pdt width, pdt src_height, pdt dst_height, double (*filter)(double), double window, double norm) {
    pdt *bounds;
    double *coeffs;
    pdt support;
    if (gen_discrete_filter(&bounds, &coeffs, &support, src_height, dst_height, filter, window, norm))
        return NULL;

    double *dst = v_filter(src, width, dst_height, bounds, coeffs, support);

    free(coeffs);
    free(bounds);
    return dst;
}

double *reconstruct(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, bool nop) {
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

    const pdt h_intermediate = dst_width * src_height;
    const pdt v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
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

static void h_iconvolve_b2_ip(double *img, pdt width, pdt height, const double *LU, pdt m) {
    const double c = LU[m];
    pdt n = m - 1;
    if (m >= width)
        n = m = width - 1;

    const double L_0 = LU[0];
    const double L_inf = LU[m - 1];
    const double L_inf_mul = L_inf * c;
    const double L_inf_div = LU[n] / c;
    const pdt f_adj = (width << 2) + 4;
    double *f = img + 4;

    for (pdt y = 0; y < height; y++) {
        pdt x;

        for (x = 1; x < m; x++, f += 4) {
            const double L_x = LU[x - 1];
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
            const double L_x = LU[x];
            f[0] = L_x * (f[0] - f[4]);
            f[1] = L_x * (f[1] - f[5]);
            f[2] = L_x * (f[2] - f[6]);
            f[3] = L_x * (f[3] - f[7]);
        }

        f[0] = L_0 * (f[0] - c * f[4]);
        f[1] = L_0 * (f[1] - c * f[5]);
        f[2] = L_0 * (f[2] - c * f[6]);
        f[3] = L_0 * (f[3] - c * f[7]);

        f += f_adj;
    }
}

static void v_iconvolve_b2_ip(double *img, pdt width, pdt height, const double *LU, pdt m) {
    const double c = LU[m];
    pdt n = m - 1;
    if (m >= height)
        n = m = height - 1;

    const double L_0 = LU[0];
    const double L_inf = LU[m - 1];
    const double L_inf_mul = L_inf * c;
    const double L_inf_div = LU[n] / c;
    const pdt adj_width = width << 2;
    double *pf = img;
    double *f = pf + adj_width;

    pdt y;

    for (y = 1; y < m; y++) {
        const double L_y = LU[y - 1];
        for (pdt x = 0; x < adj_width; x++, f++, pf++)
            f[0] -= L_y * pf[0];
    }

    for (; y < height - 1; y++)
        for (pdt x = 0; x < adj_width; x++, f++, pf++)
            f[0] -= L_inf * pf[0];

    for (pdt x = 0; x < adj_width; x++, f++, pf++)
        f[0] = L_inf_div * (f[0] - L_inf_mul * pf[0]);

    pf = f - 1;
    f = pf - adj_width;

    for (y = height - 2; y >= m; y--)
        for (pdt x = 0; x < adj_width; x++, f--, pf--)
            f[0] = L_inf * (f[0] - pf[0]);

    for (; y > 0; y--) {
        const double L_y = LU[y];
        for (pdt x = 0; x < adj_width; x++, f--, pf--)
            f[0] = L_y * (f[0] - pf[0]);
    }

    for (pdt x = 0; x < adj_width; x++, f--, pf--)
        f[0] = L_0 * (f[0] - c * pf[0]);
}

static double *h_reconstruct_iconvolve_b2(const double *src, pdt src_width, pdt height, pdt dst_width, double (*filter)(double), double window, double norm, const double *L, pdt m) {
    if (dst_width > src_width) {
        if (src_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *temp = imgcpy(src, src_width, height);
        if (!temp)
            return NULL;
        h_iconvolve_b2_ip(temp, src_width, height, L, m);

        double *ret = h_reconstruct(temp, src_width, height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *ret = h_reconstruct(src, src_width, height, dst_width, filter, window, norm);
        if (!ret)
            return NULL;

        h_iconvolve_b2_ip(ret, dst_width, height, L, m);
        return ret;
    }
}

static double *v_reconstruct_iconvolve_b2(const double *src, pdt width, pdt src_height, pdt dst_height, double (*filter)(double), double window, double norm, const double *L, pdt m) {
    if (dst_height > src_height) {
        if (src_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *temp = imgcpy(src, width, src_height);
        if (!temp)
            return NULL;
        v_iconvolve_b2_ip(temp, width, src_height, L, m);

        double *ret = v_reconstruct(temp, width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *ret = v_reconstruct(src, width, src_height, dst_height, filter, window, norm);
        if (!ret)
            return NULL;

        v_iconvolve_b2_ip(ret, width, dst_height, L, m);
        return ret;
    }
}

double *reconstruct_iconvolve_b2(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, const double *L, pdt m, bool nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return imgcpy(src, src_width, src_height);
            return v_reconstruct_iconvolve_b2(src, src_width, src_height, dst_height, filter, window, norm, L, m);
        } else if (dst_height == src_height) {
            return h_reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, filter, window, norm, L, m);
        }
    }

    const pdt h_intermediate = dst_width * src_height;
    const pdt v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        double *temp = v_reconstruct_iconvolve_b2(src, src_width, src_height, dst_height, filter, window, norm, L, m);
        if (!temp)
            return NULL;

        double *ret = h_reconstruct_iconvolve_b2(temp, src_width, dst_height, dst_width, filter, window, norm, L, m);
        free(temp);
        return ret;
    } else {
        double *temp = h_reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, filter, window, norm, L, m);
        if (!temp)
            return NULL;

        double *ret = v_reconstruct_iconvolve_b2(temp, dst_width, src_height, dst_height, filter, window, norm, L, m);

        free(temp);
        return ret;
    }
}

static void h_filter_ip(double *img, pdt width, pdt height, const pdt *bounds, const double *coeffs, pdt support, bool reverse) {
    const pdt adj_width = width << 2;
    const pdt *bounds_ptr;
    const double *coeffs_ptr;

    const pdt dst_change = reverse ? -4 : 4;
    pdt src_offset = 0;
    pdt dst_offset = reverse ? adj_width - 4 : 0;
    for (pdt y = 0; y < height; y++, src_offset += adj_width, dst_offset += adj_width) {
        coeffs_ptr = coeffs;
        bounds_ptr = bounds;

        pdt dst_pixel = dst_offset;
        for (pdt x = 0; x < width; x++, bounds_ptr += 2, coeffs_ptr += support, dst_pixel += dst_change) {
            const pdt min_x = bounds_ptr[0];
            const pdt max_x = bounds_ptr[1];

            pdt src_pixel = src_offset + (min_x << 2);
            for (pdt s = min_x, i = 0; s <= max_x; s++, i++, src_pixel += 4) {
                const double weight = coeffs_ptr[i];

                if (dst_pixel == src_pixel) {
                    img[dst_pixel + 0] = img[src_pixel + 0] * weight;
                    img[dst_pixel + 1] = img[src_pixel + 1] * weight;
                    img[dst_pixel + 2] = img[src_pixel + 2] * weight;
                    img[dst_pixel + 3] = img[src_pixel + 3] * weight;
                } else {
                    img[dst_pixel + 0] += img[src_pixel + 0] * weight;
                    img[dst_pixel + 1] += img[src_pixel + 1] * weight;
                    img[dst_pixel + 2] += img[src_pixel + 2] * weight;
                    img[dst_pixel + 3] += img[src_pixel + 3] * weight;
                }
            }
        }
    }
}

static void v_filter_ip(double *img, pdt width, pdt height, const pdt *bounds, const double *coeffs, pdt support, bool reverse) {
    const pdt adj_width = width << 2;
    const pdt *bounds_ptr = bounds;
    const double *coeffs_ptr = coeffs;

    const pdt dst_change = reverse ? -adj_width : adj_width;
    pdt dst_offset = reverse ? adj_width * (height - 1) : 0;
    for (pdt y = 0; y < height; y++, bounds_ptr += 2, coeffs_ptr += support, dst_offset += dst_change) {
        const pdt min_y = bounds_ptr[0];
        const pdt max_y = bounds_ptr[1];

        pdt src_pixel = adj_width * min_y;
        for (pdt s = min_y, i = 0; s <= max_y; s++, i++) {
            const double weight = coeffs_ptr[i];

            pdt dst_pixel = dst_offset;
            for (pdt x = 0; x < width; x++, src_pixel += 4, dst_pixel += 4) {
                if (dst_pixel == src_pixel) {
                    img[dst_pixel + 0] = img[src_pixel + 0] * weight;
                    img[dst_pixel + 1] = img[src_pixel + 1] * weight;
                    img[dst_pixel + 2] = img[src_pixel + 2] * weight;
                    img[dst_pixel + 3] = img[src_pixel + 3] * weight;
                } else {
                    img[dst_pixel + 0] += img[src_pixel + 0] * weight;
                    img[dst_pixel + 1] += img[src_pixel + 1] * weight;
                    img[dst_pixel + 2] += img[src_pixel + 2] * weight;
                    img[dst_pixel + 3] += img[src_pixel + 3] * weight;
                }
            }
        }
    }
}

static bool gen_iconvolve_filters(pdt **_forward_bounds, double **_forward_coeffs, pdt **_backward_bounds, double **_backward_coeffs, pdt len, const double *LU, pdt m, pdt b) {
    if (len < 2 * b - 1)
        return true;

#define L(x, y) LU[(x) * m + min(y, m - 1)]
#define C(x, y) _C[(x) * special + (y)]
    const pdt special = b - 1;
    const double *_C = LU + special * m;
    const pdt cutoff = len - special;

    pdt *forward_bounds = malloc(2 * len * sizeof(*forward_bounds));
    if (!forward_bounds)
        return true;
    pdt *backward_bounds = malloc(2 * len * sizeof(*backward_bounds));
    if (!backward_bounds) {
        free(forward_bounds);
        return true;
    }

    pdt support = special;
    double *forward_coeffs = malloc(support * len * sizeof(*forward_coeffs));
    if (!forward_coeffs) {
        free(forward_bounds);
        free(backward_bounds);
        return true;
    }
    double *backward_coeffs = malloc(b * len * sizeof(*backward_coeffs));
    if (!backward_bounds) {
        free(forward_bounds);
        free(backward_bounds);
        free(forward_coeffs);
        return true;
    }

    *_forward_bounds = forward_bounds;
    *_forward_coeffs = forward_coeffs;
    *_backward_bounds = backward_bounds;
    *_backward_coeffs = backward_coeffs;

    pdt z;

    forward_bounds[0] = 0;
    forward_bounds[1] = -1;
    forward_bounds += 2;
    forward_coeffs += support;

    for (z = 1; z < cutoff; z++, forward_bounds += 2, forward_coeffs += support) {
        const pdt min_z = max(z - support, 0);
        const pdt max_z = z - 1;
        forward_bounds[0] = min_z;
        forward_bounds[1] = max_z;

        for (pdt s = max_z, i = 0; s >= min_z; s--, i++)
            forward_coeffs[max_z - min_z - i] = -L(i, (z - 1) - i);
    }

    for (; z < len; z++, forward_bounds += 2, forward_coeffs += support) {
        const pdt min_z = z - support;
        const pdt max_z = z - 1;
        forward_bounds[0] = min_z;
        forward_bounds[1] = max_z;

        for (pdt s = max_z, i = 0; s >= min_z; s--, i++)
            forward_coeffs[support - 1 - i] = -C(i, (len - 1) - z) * L(i, (z - 1) - i);
    }

    support = b;

    for (z = len - 1; z >= cutoff; z--, backward_bounds += 2, backward_coeffs += support) {
        const double M = L(special - 1, z) / C(special - 1, (len - 1) - z);
        const pdt min_z = z;
        const pdt max_z = min(z + support - 1, len - 1);
        backward_bounds[0] = min_z;
        backward_bounds[1] = max_z;

        backward_coeffs[0] = M;
        for (pdt s = min_z + 1, i = 1; s <= max_z; s++, i++)
            backward_coeffs[i] = -L(i - 1, z);
    }

    for (; z >= special; z--, backward_bounds += 2, backward_coeffs += support) {
        const double M = L(special - 1, z);
        const pdt min_z = z;
        const pdt max_z = z + support - 1;
        backward_bounds[0] = min_z;
        backward_bounds[1] = max_z;

        backward_coeffs[0] = M;
        for (pdt s = min_z + 1, i = 1; s < max_z; s++, i++)
            backward_coeffs[i] = -L(i - 1, z);
        backward_coeffs[support - 1] = -M;
    }

    for (; z >= 0; z--, backward_bounds += 2, backward_coeffs += support) {
        const double M = L(special - 1, z);
        const pdt min_z = z;
        const pdt max_z = z + support - 1;
        backward_bounds[0] = min_z;
        backward_bounds[1] = max_z;

        backward_coeffs[0] = M;
        for (pdt s = min_z + 1, i = 1; s < max_z; s++, i++)
            backward_coeffs[i] = -C(i - 1, z) * L(i - 1, z);
        backward_coeffs[support - 1] = -M * C(special - 1, z);
    }

#undef L
#undef C
    return false;
}

static bool h_iconvolve_ip(double *img, pdt width, pdt height, const double *LU, pdt m, pdt b) {
    pdt *forward_bounds;
    pdt *backward_bounds;
    double *forward_coeffs;
    double *backward_coeffs;
    if (gen_iconvolve_filters(&forward_bounds, &forward_coeffs, &backward_bounds, &backward_coeffs, width, LU, m, b))
        return true;

    h_filter_ip(img, width, height, forward_bounds, forward_coeffs, b - 1, false);
    h_filter_ip(img, width, height, backward_bounds, backward_coeffs, b, true);

    free(forward_bounds);
    free(forward_coeffs);
    free(backward_bounds);
    free(backward_coeffs);
    return false;
}

static bool v_iconvolve_ip(double *img, pdt width, pdt height, const double *LU, pdt m, pdt b) {
    pdt *forward_bounds;
    pdt *backward_bounds;
    double *forward_coeffs;
    double *backward_coeffs;
    if (gen_iconvolve_filters(&forward_bounds, &forward_coeffs, &backward_bounds, &backward_coeffs, height, LU, m, b))
        return true;

    v_filter_ip(img, width, height, forward_bounds, forward_coeffs, b - 1, false);
    v_filter_ip(img, width, height, backward_bounds, backward_coeffs, b, true);

    free(forward_bounds);
    free(forward_coeffs);
    free(backward_bounds);
    free(backward_coeffs);
    return false;
}

static double *h_reconstruct_iconvolve(const double *src, pdt src_width, pdt height, pdt dst_width, double (*filter)(double), double window, double norm, const double *LU, pdt m, pdt b) {
    if (dst_width > src_width) {
        if (src_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *temp = imgcpy(src, src_width, height);
        if (!temp)
            return NULL;
        if (h_iconvolve_ip(temp, src_width, height, LU, m, b)) {
            free(temp);
            return NULL;
        }

        double *ret = h_reconstruct(temp, src_width, height, dst_width, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_width == 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        double *ret = h_reconstruct(src, src_width, height, dst_width, filter, window, norm);
        if (!ret)
            return NULL;

        if (h_iconvolve_ip(ret, dst_width, height, LU, m, b)) {
            free(ret);
            return NULL;
        }
        return ret;
    }
}

static double *v_reconstruct_iconvolve(const double *src, pdt width, pdt src_height, pdt dst_height, double (*filter)(double), double window, double norm, const double *LU, pdt m, pdt b) {
    if (dst_height > src_height) {
        if (src_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *temp = imgcpy(src, width, src_height);
        if (!temp)
            return NULL;
        if (v_iconvolve_ip(temp, width, src_height, LU, m, b)) {
            free(temp);
            return NULL;
        }

        double *ret = v_reconstruct(temp, width, src_height, dst_height, filter, window, norm);
        free(temp);
        return ret;
    } else {
        if (dst_height == 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        double *ret = v_reconstruct(src, width, src_height, dst_height, filter, window, norm);
        if (!ret)
            return NULL;

        if (v_iconvolve_ip(ret, width, dst_height, LU, m, b)) {
            free(ret);
            return NULL;
        }
        return ret;
    }
}

double *reconstruct_iconvolve(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, double (*filter)(double), double window, double norm, const double *LU, pdt m, pdt b, bool nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return NULL;

    if (nop) {
        if (dst_width == src_width) {
            if (dst_height == src_height)
                return imgcpy(src, src_width, src_height);
            return v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, LU, m, b);
        } else if (dst_height == src_height) {
            return h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, LU, m, b);
        }
    }

    const pdt h_intermediate = dst_width * src_height;
    const pdt v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        double *temp = v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, LU, m, b);
        if (!temp)
            return NULL;

        double *ret = h_reconstruct_iconvolve(temp, src_width, dst_height, dst_width, filter, window, norm, LU, m, b);
        free(temp);
        return ret;
    } else {
        double *temp = h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, LU, m, b);
        if (!temp)
            return NULL;

        double *ret = v_reconstruct_iconvolve(temp, dst_width, src_height, dst_height, filter, window, norm, LU, m, b);
        free(temp);
        return ret;
    }
}

double *resize(const double *src, pdt src_width, pdt src_height, pdt dst_width, pdt dst_height, Filter filter) {
    switch (filter) {
    case NEAREST:
        return sample(src, src_width, src_height, dst_width, dst_height);
    case AREA:
        return scale(src, src_width, src_height, dst_width, dst_height);
    case TRIANGLE:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Triangle, 1.0, 1.0, true);
    case HERMITE:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hermite, 1.0, 1.0, true);
    case BSPLINE2:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, BSpline2, 1.5, 1.0, false);
    case BSPLINE3:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, BSpline3, 2.0, 1.0, false);
    default:
    case MITNET:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, MitNet, 2.0, 1.0, false);
    case CATROM:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, CatRom, 2.0, 1.0, true);
    case HAMMING3:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hamming3, 3.0, 1.0, true);
    case HAMMING4:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hamming4, 4.0, 1.0, true);
    case HAMMING8:
        return reconstruct(src, src_width, src_height, dst_width, dst_height, Hamming8, 8.0, 1.0, true);
    case BSPLINE2I:
        return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, BSpline2, 1.5, 8.0, LU_bspline2i, 11, true);
    case BSPLINE3I:
        return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, BSpline3, 2.0, 6.0, LU_bspline3i, 14, true);
    case OMOMS3I:
        return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, OMOMS3, 2.0, 5.25, LU_omoms3, 18, true);
    case OMOMS7I:
        return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, OMOMS7, 4.0, 1952.8179190751446, LU_omoms7, 35, 4, true);
    case OMOMS11I:
        return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, OMOMS11, 6.0, 5148314.024847966, LU_omoms11, 54, 6, true);
    }
}

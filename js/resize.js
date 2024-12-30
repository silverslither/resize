// Copyright (c) 2024 silverslither.

import * as Filters from "./filters.js";

/**
 * Resample an image using nearest neighbor interpolation.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @returns {TypedArray} Destination image in a 4-channel format, with the same type as source image.
 */
export function sample(src, src_width, src_height, dst_width, dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new (Object.getPrototypeOf(src).constructor)();

    if (src_width === dst_width && src_height === dst_height)
        return new (Object.getPrototypeOf(src).constructor)(src);

    const x_factor = src_width / dst_width;
    const y_factor = src_height / dst_height;

    const dst = new (Object.getPrototypeOf(src).constructor)((dst_width * dst_height) << 2);

    let dst_pixel = 0;
    for (let y = 0; y < dst_height; y++) {
        const mapped_y = Math.ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (let x = 0; x < dst_width; x++, dst_pixel += 4) {
            const mapped_x = Math.ceil(x_factor * (x + 0.5) - 1.0);
            const src_pixel = (mapped_y + mapped_x) << 2;
            dst[dst_pixel + 0] = src[src_pixel + 0];
            dst[dst_pixel + 1] = src[src_pixel + 1];
            dst[dst_pixel + 2] = src[src_pixel + 2];
            dst[dst_pixel + 3] = src[src_pixel + 3];
        }
    }

    return dst;
}

function hscale(src, src_width, height, dst_width) {
    const max_high_x = src_width - 1;
    const adj_src_width = src_width << 2;
    const adj_dst_width = dst_width << 2;
    const adj_src_area = adj_src_width * height;
    const factor = src_width / dst_width;
    const inv_factor = dst_width / src_width;

    const dst = new Float64Array(height * adj_dst_width);

    let dst_offset = 0;
    for (let x = 0; x < dst_width; x++, dst_offset += 4) {
        const min_mapped_x = factor * x;
        const max_mapped_x = min_mapped_x + factor;
        let min_x = Math.ceil(min_mapped_x);
        let max_x = Math.floor(max_mapped_x);

        let dst_pixel = dst_offset;

        if (max_x < min_x) {
            for (let src_pixel = max_x << 2; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] = src[src_pixel + 0];
                dst[dst_pixel + 1] = src[src_pixel + 1];
                dst[dst_pixel + 2] = src[src_pixel + 2];
                dst[dst_pixel + 3] = src[src_pixel + 3];
            }
            continue;
        }

        const low_x = Math.max(min_x - 1, 0) << 2;
        const high_x = Math.min(max_x, max_high_x) << 2;
        const low_mult = min_x - min_mapped_x;
        const high_mult = max_mapped_x - max_x;
        min_x <<= 2;
        max_x <<= 2;

        for (let y = 0; y < adj_src_area; y += adj_src_width, dst_pixel += adj_dst_width) {
            let src_pixel = y + low_x;
            dst[dst_pixel + 0] += src[src_pixel + 0] * low_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * low_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * low_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * low_mult;

            const src_max = y + max_x;
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

function vscale(src, width, src_height, dst_height) {
    const max_high_y = src_height - 1;
    const adj_width = width << 2;
    const factor = src_height / dst_height;
    const inv_factor = dst_height / src_height;

    const dst = new Float64Array(adj_width * dst_height);

    let dst_pixel = 0;
    for (let y = 0; y < dst_height; y++) {
        const min_mapped_y = factor * y;
        const max_mapped_y = min_mapped_y + factor;
        let min_y = Math.ceil(min_mapped_y);
        let max_y = Math.floor(max_mapped_y);

        if (max_y < min_y) {
            const src_offset = adj_width * max_y;
            const src_max = src_offset + adj_width;
            for (let src_pixel = src_offset; src_pixel < src_max; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] = src[src_pixel + 0];
                dst[dst_pixel + 1] = src[src_pixel + 1];
                dst[dst_pixel + 2] = src[src_pixel + 2];
                dst[dst_pixel + 3] = src[src_pixel + 3];
            }
            continue;
        }

        const low_y = adj_width * Math.max(min_y - 1, 0);
        const high_y = adj_width * Math.min(max_y, max_high_y);
        const low_mult = min_y - min_mapped_y;
        const high_mult = max_mapped_y - max_y;
        min_y *= adj_width;
        max_y *= adj_width;

        for (let x = 0; x < adj_width; x += 4, dst_pixel += 4) {
            let src_pixel = x + low_y;
            dst[dst_pixel + 0] += src[src_pixel + 0] * low_mult;
            dst[dst_pixel + 1] += src[src_pixel + 1] * low_mult;
            dst[dst_pixel + 2] += src[src_pixel + 2] * low_mult;
            dst[dst_pixel + 3] += src[src_pixel + 3] * low_mult;

            const src_max = x + max_y;
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

/**
 * Resize an image using area averaging / pixel mixing.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @returns {Float64Array} Destination image in a 4-channel format.
 */
export function scale(src, src_width, src_height, dst_width, dst_height) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (dst_width === src_width) {
        if (dst_height === src_height)
            return new Float64Array(src);
        return vscale(src, src_width, src_height, dst_height);
    } else if (dst_height === src_height) {
        return hscale(src, src_width, src_height, dst_width);
    }

    const x_factor = dst_width / src_width;
    const y_factor = dst_height / src_height;

    if (x_factor >= y_factor) {
        const temp = vscale(src, src_width, src_height, dst_height);
        return hscale(temp, src_width, dst_height, dst_width);
    } else {
        const temp = hscale(src, src_width, src_height, dst_width);
        return vscale(temp, dst_width, src_height, dst_height);
    }
}

function hreconstruct(src, src_width, height, dst_width, filter, window, norm) {
    const max_s = src_width - 1;
    const adj_src_width = src_width << 2;
    const adj_dst_width = dst_width << 2;
    const adj_dst_area = adj_dst_width * height;
    const factor = src_width / dst_width;

    const dst = new Float64Array(adj_dst_width * height);

    const inv_filter_scale = Math.max(factor, 1.0);
    const filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    let dst_offset = 0;
    for (let x = 0; x < dst_width; x++, dst_offset += 4) {
        const mapped_x = factor * (x + 0.5) - 0.5;
        const min_x = Math.max(Math.ceil(mapped_x - window), 0);
        const max_x = Math.min(Math.floor(mapped_x + window), max_s);
        let weight_total = 0;

        let src_offset = min_x << 2;
        for (let s = min_x; s <= max_x; s++, src_offset += 4) {
            const weight = filter(Math.abs(mapped_x - s) * filter_scale);
            weight_total += weight;

            for (let src_pixel = src_offset, dst_pixel = dst_offset; dst_pixel < adj_dst_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        weight_total = norm / weight_total;
        for (let dst_pixel = dst_offset; dst_pixel < adj_dst_area; dst_pixel += adj_dst_width) {
            dst[dst_pixel + 0] *= weight_total;
            dst[dst_pixel + 1] *= weight_total;
            dst[dst_pixel + 2] *= weight_total;
            dst[dst_pixel + 3] *= weight_total;
        }
    }

    return dst;
}

function vreconstruct(src, width, src_height, dst_height, filter, window, norm) {
    const max_s = src_height - 1;
    const adj_width = width << 2;
    const factor = src_height / dst_height;

    const dst = new Float64Array(adj_width * dst_height);

    const inv_filter_scale = Math.max(factor, 1.0);
    const filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    let dst_offset = 0;
    let ndst_offset = adj_width;
    for (let y = 0; y < dst_height; y++, dst_offset = ndst_offset, ndst_offset += adj_width) {
        const mapped_y = factor * (y + 0.5) - 0.5;
        const min_y = Math.max(Math.ceil(mapped_y - window), 0);
        const max_y = Math.min(Math.floor(mapped_y + window), max_s);
        let weight_total = 0;

        let src_offset = adj_width * min_y;
        for (let s = min_y; s <= max_y; s++, src_offset += adj_width) {
            const weight = filter(Math.abs(mapped_y - s) * filter_scale);
            weight_total += weight;

            for (let src_pixel = src_offset, dst_pixel = dst_offset; dst_pixel < ndst_offset; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        weight_total = norm / weight_total;
        for (let dst_pixel = dst_offset; dst_pixel < ndst_offset; dst_pixel += 4) {
            dst[dst_pixel + 0] *= weight_total;
            dst[dst_pixel + 1] *= weight_total;
            dst[dst_pixel + 2] *= weight_total;
            dst[dst_pixel + 3] *= weight_total;
        }
    }

    return dst;
}

/**
 * Resize an image using a reconstruction filter.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {(x: number) => number} filter Reconstruction filter function.
 * @param {number} window Filter function window.
 * @param {number} norm Normalization constant.
 * @param {boolean} nop Boolean flag for a no-op case.
 * @returns {Float64Array} Destination image in a 4-channel format.
 */
export function reconstruct(src, src_width, src_height, dst_width, dst_height, filter, window, norm, nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (nop) {
        if (dst_width === src_width) {
            if (dst_height === src_height)
                return new Float64Array(src);
            return vreconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        } else if (dst_height === src_height) {
            return hreconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        }
    }

    const x_factor = dst_width / src_width;
    const y_factor = dst_height / src_height;

    if (x_factor >= y_factor) {
        const temp = vreconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        return hreconstruct(temp, src_width, dst_height, dst_width, filter, window, norm);
    } else {
        const temp = hreconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        return vreconstruct(temp, dst_width, src_height, dst_height, filter, window, norm);
    }
}

function hiconvolve(img, width, height, L, m, c) {
    let n = m - 1;
    if (m >= width)
        n = m = width - 1;

    const L_0 = L[0];
    const L_inf = L[m - 1];
    const L_inf_mul = L_inf * c;
    const L_inf_div = L[n] / c;
    const f_adj = (width << 2) + 8;
    let f = 4;

    for (let y = 0; y < height; y++) {
        let x;

        for (x = 1; x < m; x++, f += 4) {
            const L_x = L[x - 1];
            img[f + 0] -= L_x * img[f - 4];
            img[f + 1] -= L_x * img[f - 3];
            img[f + 2] -= L_x * img[f - 2];
            img[f + 3] -= L_x * img[f - 1];
        }

        for (; x < width - 1; x++, f += 4) {
            img[f + 0] -= L_inf * img[f - 4];
            img[f + 1] -= L_inf * img[f - 3];
            img[f + 2] -= L_inf * img[f - 2];
            img[f + 3] -= L_inf * img[f - 1];
        }

        img[f + 0] = L_inf_div * (img[f + 0] - L_inf_mul * img[f - 4]);
        img[f + 1] = L_inf_div * (img[f + 1] - L_inf_mul * img[f - 3]);
        img[f + 2] = L_inf_div * (img[f + 2] - L_inf_mul * img[f - 2]);
        img[f + 3] = L_inf_div * (img[f + 3] - L_inf_mul * img[f - 1]);
        f -= 4;

        for (x = width - 2; x >= m; x--, f -= 4) {
            img[f + 0] = L_inf * (img[f + 0] - img[f + 4]);
            img[f + 1] = L_inf * (img[f + 1] - img[f + 5]);
            img[f + 2] = L_inf * (img[f + 2] - img[f + 6]);
            img[f + 3] = L_inf * (img[f + 3] - img[f + 7]);
        }

        for (; x > 0; x--, f -= 4) {
            const L_x = L[x];
            img[f + 0] = L_x * (img[f + 0] - img[f + 4]);
            img[f + 1] = L_x * (img[f + 1] - img[f + 5]);
            img[f + 2] = L_x * (img[f + 2] - img[f + 6]);
            img[f + 3] = L_x * (img[f + 3] - img[f + 7]);
        }

        img[f + 0] = L_0 * (img[f + 0] - c * img[f + 4]);
        img[f + 1] = L_0 * (img[f + 1] - c * img[f + 5]);
        img[f + 2] = L_0 * (img[f + 2] - c * img[f + 6]);
        img[f + 3] = L_0 * (img[f + 3] - c * img[f + 7]);

        f += f_adj - 4;
    }
}

function viconvolve(img, width, height, L, m, c) {
    let n = m - 1;
    if (m >= height)
        n = m = height - 1;

    const L_0 = L[0];
    const L_inf = L[m - 1];
    const L_inf_mul = L_inf * c;
    const L_inf_div = L[n] / c;
    const adj_width = width << 2;
    let pf = 0;
    let f = adj_width;

    for (let x = 0; x < width; x++) {
        let y;

        for (y = 1; y < m; y++, pf = f, f += adj_width) {
            const L_y = L[y - 1];
            img[f + 0] -= L_y * img[pf + 0];
            img[f + 1] -= L_y * img[pf + 1];
            img[f + 2] -= L_y * img[pf + 2];
            img[f + 3] -= L_y * img[pf + 3];
        }

        for (; y < height - 1; y++, pf = f, f += adj_width) {
            img[f + 0] -= L_inf * img[pf + 0];
            img[f + 1] -= L_inf * img[pf + 1];
            img[f + 2] -= L_inf * img[pf + 2];
            img[f + 3] -= L_inf * img[pf + 3];
        }

        img[f + 0] = L_inf_div * (img[f + 0] - L_inf_mul * img[pf + 0]);
        img[f + 1] = L_inf_div * (img[f + 1] - L_inf_mul * img[pf + 1]);
        img[f + 2] = L_inf_div * (img[f + 2] - L_inf_mul * img[pf + 2]);
        img[f + 3] = L_inf_div * (img[f + 3] - L_inf_mul * img[pf + 3]);
        pf = f;
        f = pf - adj_width;

        for (y = height - 2; y >= m; y--, pf = f, f -= adj_width) {
            img[f + 0] = L_inf * (img[f + 0] - img[pf + 0]);
            img[f + 1] = L_inf * (img[f + 1] - img[pf + 1]);
            img[f + 2] = L_inf * (img[f + 2] - img[pf + 2]);
            img[f + 3] = L_inf * (img[f + 3] - img[pf + 3]);
        }

        for (; y > 0; y--, pf = f, f -= adj_width) {
            const L_y = L[y];
            img[f + 0] = L_y * (img[f + 0] - img[pf + 0]);
            img[f + 1] = L_y * (img[f + 1] - img[pf + 1]);
            img[f + 2] = L_y * (img[f + 2] - img[pf + 2]);
            img[f + 3] = L_y * (img[f + 3] - img[pf + 3]);
        }

        img[f + 0] = L_0 * (img[f + 0] - c * img[pf + 0]);
        img[f + 1] = L_0 * (img[f + 1] - c * img[pf + 1]);
        img[f + 2] = L_0 * (img[f + 2] - c * img[pf + 2]);
        img[f + 3] = L_0 * (img[f + 3] - c * img[pf + 3]);

        pf = f + 4;
        f = pf + adj_width;
    }
}

function hreconstruct_iconvolve(src, src_width, height, dst_width, filter, window, norm, L, m, c) {
    if (dst_width > src_width) {
        if (src_width === 1)
            return hreconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const temp = new Float64Array(src);
        hiconvolve(temp, src_width, height, L, m, c);
        return hreconstruct(temp, src_width, height, dst_width, filter, window, norm);
    } else {
        if (dst_width === 1)
            return hreconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const ret = hreconstruct(src, src_width, height, dst_width, filter, window, norm);
        hiconvolve(ret, dst_width, height, L, m, c);
        return ret;
    }
}

function vreconstruct_iconvolve(src, width, src_height, dst_height, filter, window, norm, L, m, c) {
    if (dst_height > src_height) {
        if (src_height === 1)
            return vreconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const temp = new Float64Array(src);
        viconvolve(temp, width, src_height, L, m, c);
        return vreconstruct(temp, width, src_height, dst_height, filter, window, norm);
    } else {
        if (dst_height === 1)
            return vreconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const ret = vreconstruct(src, width, src_height, dst_height, filter, window, norm);
        viconvolve(ret, width, dst_height, L, m, c);
        return ret;
    }
}

/**
 * Resize an image using a reconstruction filter and an inverse discrete convolution.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {(x: number) => number} filter Reconstruction filter function.
 * @param {number} window Filter function window.
 * @param {number} norm Normalization constant.
 * @param {number[]} L Lower matrix coefficients.
 * @param {number} m Number of lower matrix coefficients.
 * @param {number} c Edge multiplier constant.
 * @param {boolean} nop Boolean flag for a no-op case.
 * @returns {Float64Array} Destination image in a 4-channel format.
 */
export function reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, filter, window, norm, L, m, c, nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (nop) {
        if (dst_width === src_width) {
            if (dst_height === src_height)
                return new Float64Array(src);
            return vreconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        } else if (dst_height === src_height) {
            return hreconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        }
    }

    const x_factor = dst_width / src_width;
    const y_factor = dst_height / src_height;

    if (x_factor >= y_factor) {
        const temp = vreconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, L, m, c);
        return hreconstruct_iconvolve(temp, src_width, dst_height, dst_width, filter, window, norm, L, m, c);
    } else {
        const temp = hreconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, L, m, c);
        return vreconstruct_iconvolve(temp, dst_width, src_height, dst_height, filter, window, norm, L, m, c);
    }
}

/**
 * @enum {number}
 */
export const Filter = {
    DEFAULT: 0,
    NEAREST: 1,
    AREA: 2,
    TRIANGLE: 3,
    HERMITE: 4,
    LAGRANGE_2: 5,
    LAGRANGE_3: 6,
    B_SPLINE_2: 7,
    B_SPLINE_3: 8,
    MITNET: 9,
    CATROM: 10,
    MKS_2013: 11,
    LANCZOS_3: 12,
    LANCZOS_4: 13,
    HAMMING_3: 14,
    HAMMING_4: 15,
    B_SPLINE_3_I: 16,
    O_MOMS_3_I: 17
};

/**
 * Wrapper for `sample`, `scale`, `reconstruct`, and `reconstruct_iconvolve`.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {Filter} filter Resizing method (filter) to be used. Defaults to Mitchell-Netravali.
 * @returns {TypedArray | Float64Array} Destination image in a 4-channel format. Returns the same type as source image if `filter` is `NEAREST`, otherwise returns Float64Array.
 */
export function resize(src, src_width, src_height, dst_width, dst_height, filter) {
    switch (filter) {
        case Filter.NEAREST:
            return sample(src, src_width, src_height, dst_width, dst_height);
        case Filter.AREA:
            return scale(src, src_width, src_height, dst_width, dst_height);
        case Filter.TRIANGLE:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Triangle, 1.0, 1.0, 1);
        case Filter.HERMITE:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hermite, 1.0, 1.0, 1);
        case Filter.LAGRANGE_2:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Lagrange2, 1.5, 1.0, 1);
        case Filter.LAGRANGE_3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Lagrange3, 2.0, 1.0, 1);
        case Filter.B_SPLINE_2:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.BSpline2, 1.5, 1.0, 0);
        case Filter.B_SPLINE_3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.BSpline3, 2.0, 1.0, 0);
        default:
        case Filter.MITNET:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.MitNet, 2.0, 1.0, 0);
        case Filter.CATROM:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.CatRom, 2.0, 1.0, 1);
        case Filter.MKS_2013:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.MKS2013, 2.5, 1.0, 0);
        case Filter.LANCZOS_3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Lanczos3, 3.0, 1.0, 1);
        case Filter.LANCZOS_4:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Lanczos4, 4.0, 1.0, 1);
        case Filter.HAMMING_3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hamming3, 3.0, 1.0, 1);
        case Filter.HAMMING_4:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hamming4, 4.0, 1.0, 1);
        case Filter.B_SPLINE_3_I:
            return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, Filters.BSpline3, 2.0, 6.0, Filters.L_bspline3i, 14, 1.2, 1);
        case Filter.O_MOMS_3_I:
            return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, Filters.OMOMS3, 2.0, 5.25, Filters.L_omoms3, 18, 1.2352941176470589, 1);
    }
}

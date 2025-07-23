// Copyright (c) 2024-2025 silverslither.

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

    const dst = new (Object.getPrototypeOf(src).constructor)((dst_width * dst_height) * 4);

    let dst_pixel = 0;
    for (let y = 0; y < dst_height; y++) {
        const mapped_y = Math.ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (let x = 0; x < dst_width; x++, dst_pixel += 4) {
            const mapped_x = Math.ceil(x_factor * (x + 0.5) - 1.0);
            const src_pixel = (mapped_y + mapped_x) * 4;
            dst[dst_pixel + 0] = src[src_pixel + 0];
            dst[dst_pixel + 1] = src[src_pixel + 1];
            dst[dst_pixel + 2] = src[src_pixel + 2];
            dst[dst_pixel + 3] = src[src_pixel + 3];
        }
    }

    return dst;
}

function h_filter(src, src_width, height, dst_width, bounds, coeffs, support) {
    const adj_src_width = src_width * 4;
    const dst = new Float64Array((dst_width * height) * 4);

    let bounds_ptr;
    let coeffs_ptr;
    let src_offset = 0;
    let dst_pixel = 0;
    for (let y = 0; y < height; y++, src_offset += adj_src_width) {
        coeffs_ptr = 0;
        bounds_ptr = 0;

        for (let x = 0; x < dst_width; x++, bounds_ptr += 2, coeffs_ptr += support, dst_pixel += 4) {
            const min_x = bounds[bounds_ptr + 0];
            const max_x = bounds[bounds_ptr + 1];

            let src_pixel = src_offset + (min_x * 4);
            for (let s = min_x, i = 0; s <= max_x; s++, i++, src_pixel += 4) {
                const weight = coeffs[coeffs_ptr + i];

                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

function v_filter(src, width, dst_height, bounds, coeffs, support) {
    const adj_width = width * 4;
    const dst = new Float64Array(adj_width * dst_height);

    let bounds_ptr = 0;
    let coeffs_ptr = 0;
    let dst_offset = 0;
    for (let y = 0; y < dst_height; y++, bounds_ptr += 2, coeffs_ptr += support, dst_offset += adj_width) {
        const min_y = bounds[bounds_ptr + 0];
        const max_y = bounds[bounds_ptr + 1];

        let src_pixel = adj_width * min_y;
        for (let s = min_y, i = 0; s <= max_y; s++, i++) {
            const weight = coeffs[coeffs_ptr + i];

            let dst_pixel = dst_offset;
            for (let x = 0; x < width; x++, src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

function gen_area_filter(src, dst) {
    const factor = src / dst;
    const inv_factor = dst / src;
    const support = Math.ceil(factor) + 1;

    const bounds = new Int32Array(2 * dst);
    const coeffs = new Float64Array(support * dst);

    let bounds_ptr = 0;
    let coeffs_ptr = 0;
    for (let z = 0; z < dst; z++, bounds_ptr += 2, coeffs_ptr += support) {
        const min_mapped_z = factor * z;
        const max_mapped_z = min_mapped_z + factor;
        const min_z = Math.floor(min_mapped_z + 1.0);
        const max_z = Math.ceil(max_mapped_z - 1.0);

        if (max_z < min_z) {
            bounds[bounds_ptr + 1] = bounds[bounds_ptr + 0] = max_z;
            coeffs[coeffs_ptr + 0] = 1.0;
            continue;
        }

        bounds[bounds_ptr + 0] = min_z - 1;
        bounds[bounds_ptr + 1] = max_z;
        coeffs[coeffs_ptr + 0] = inv_factor * (min_z - min_mapped_z);
        let i = 1;
        for (let s = min_z; s < max_z; s++, i++)
            coeffs[coeffs_ptr + i] = inv_factor;
        coeffs[coeffs_ptr + i] = inv_factor * (max_mapped_z - max_z);
    }

    return { bounds, coeffs, support };
}

function h_scale(src, src_width, height, dst_width) {
    const { bounds, coeffs, support } = gen_area_filter(src_width, dst_width);
    return h_filter(src, src_width, height, dst_width, bounds, coeffs, support);
}

function v_scale(src, width, src_height, dst_height) {
    const { bounds, coeffs, support } = gen_area_filter(src_height, dst_height);
    return v_filter(src, width, dst_height, bounds, coeffs, support);
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
        return v_scale(src, src_width, src_height, dst_height);
    } else if (dst_height === src_height) {
        return h_scale(src, src_width, src_height, dst_width);
    }

    const h_intermediate = dst_width * src_height;
    const v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        const temp = v_scale(src, src_width, src_height, dst_height);
        return h_scale(temp, src_width, dst_height, dst_width);
    } else {
        const temp = h_scale(src, src_width, src_height, dst_width);
        return v_scale(temp, dst_width, src_height, dst_height);
    }
}

function gen_kernel_filter(len, kernel, support, norm) {
    const offset = Math.floor(0.5 * support);

    const bounds = new Int32Array(2 * len);
    const coeffs = new Float64Array(support * len);

    let bounds_ptr = 0;
    let coeffs_ptr = 0;
    for (let z = 0; z < len; z++, bounds_ptr += 2, coeffs_ptr += support) {
        const min_z = z - offset;
        const max_z = Math.min(z + offset, len - 1);
        bounds[bounds_ptr + 0] = Math.max(min_z, 0);
        bounds[bounds_ptr + 1] = max_z;

        let weight_total = 0.0;
        for (let s = min_z, i = 0; s <= max_z; s++, i++) {
            if (s < 0) {
                i--;
                continue;
            }
            const weight = kernel[i];
            coeffs[coeffs_ptr + i] = weight;
            weight_total += weight;
        }

        if (Math.abs(weight_total) > 2.3283064365386963e-10) {
            weight_total = norm / weight_total;
            for (let s = min_z, i = 0; s <= max_z; s++, i++)
                coeffs[coeffs_ptr + i] *= weight_total;
        }
    }

    return { bounds, coeffs };
}

function h_convolve(src, width, height, kernel, support, norm) {
    const { bounds, coeffs } = gen_kernel_filter(width, kernel, support, norm);
    return h_filter(src, width, height, width, bounds, coeffs, support);
}

function v_convolve(src, width, height, kernel, support, norm) {
    const { bounds, coeffs } = gen_kernel_filter(height, kernel, support, norm);
    return v_filter(src, width, height, bounds, coeffs, support);
}

/**
 * Convolve an image with a horizontal and vertical kernel.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} width Image width.
 * @param {number} height Image height.
 * @param {Float64Array | null} h_kernel Horizontal kernel, or null if no horizontal convolution is desired.
 * @param {Float64Array | null} v_kernel Vertical kernel, or null if no vertical convolution is desired.
 * @param {number} h_support Support window for the horizontal kernel. Must be an odd number.
 * @param {number} v_support Support window for the vertical kernel. Must be an odd number.
 * @param {number} h_support Normalization constant for the horizontal kernel.
 * @param {number} v_support Normalization constant for the vertical kernel.
 * @returns {Float64Array} Destination image in a 4-channel format.
 */
export function convolve(src, width, height, h_kernel, v_kernel, h_support, v_support, h_norm, v_norm) {
    if (width <= 0 || height <= 0)
        return new Float64Array();

    if (h_kernel == null) {
        if (v_kernel == null)
            return new Float64Array(src);
        return v_convolve(src, width, height, v_kernel, v_support, v_norm);
    } else if (v_kernel == null) {
        return h_convolve(src, width, height, h_kernel, h_support, h_norm);
    }

    const temp = h_convolve(src, width, height, h_kernel, h_support, h_norm);
    return v_convolve(temp, width, height, v_kernel, v_support, v_norm);
}

function gen_discrete_filter(src, dst, filter, window, norm) {
    const max_s = src - 1;
    const factor = src / dst;
    const inv_filter_scale = Math.max(factor, 1.0);
    const filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;
    const support = Math.ceil(2.0 * window);

    const bounds = new Int32Array(2 * dst);
    const coeffs = new Float64Array(support * dst);

    let bounds_ptr = 0;
    let coeffs_ptr = 0;
    for (let z = 0; z < dst; z++, bounds_ptr += 2, coeffs_ptr += support) {
        const mapped_x = factor * (z + 0.5) - 0.5;
        const min_z = Math.max(Math.floor(mapped_x - window + 1.0), 0);
        const max_z = Math.min(Math.ceil(mapped_x + window - 1.0), max_s);
        bounds[bounds_ptr + 0] = min_z;
        bounds[bounds_ptr + 1] = max_z;

        let weight_total = 0.0;
        for (let s = min_z, i = 0; s <= max_z; s++, i++) {
            const weight = filter(Math.abs(mapped_x - s) * filter_scale);
            coeffs[coeffs_ptr + i] = weight;
            weight_total += weight;
        }

        if (Math.abs(weight_total) > 2.3283064365386963e-10) {
            weight_total = norm / weight_total;
            for (let s = min_z, i = 0; s <= max_z; s++, i++)
                coeffs[coeffs_ptr + i] *= weight_total;
        }
    }

    return { bounds, coeffs, support };
}

function h_reconstruct(src, src_width, height, dst_width, filter, window, norm) {
    const { bounds, coeffs, support } = gen_discrete_filter(src_width, dst_width, filter, window, norm);
    return h_filter(src, src_width, height, dst_width, bounds, coeffs, support);
}

function v_reconstruct(src, width, src_height, dst_height, filter, window, norm) {
    const { bounds, coeffs, support } = gen_discrete_filter(src_height, dst_height, filter, window, norm);
    return v_filter(src, width, dst_height, bounds, coeffs, support);
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
            return v_reconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        } else if (dst_height === src_height) {
            return h_reconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        }
    }

    const h_intermediate = dst_width * src_height;
    const v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        const temp = v_reconstruct(src, src_width, src_height, dst_height, filter, window, norm);
        return h_reconstruct(temp, src_width, dst_height, dst_width, filter, window, norm);
    } else {
        const temp = h_reconstruct(src, src_width, src_height, dst_width, filter, window, norm);
        return v_reconstruct(temp, dst_width, src_height, dst_height, filter, window, norm);
    }
}

function h_iconvolve_b2_ip(img, width, height, LU, m) {
    const c = LU[m];
    let n = m - 1;
    if (m >= width)
        n = m = width - 1;

    const L_0 = LU[0];
    const L_inf = LU[m - 1];
    const L_inf_mul = L_inf * c;
    const L_inf_div = LU[n] / c;
    const f_adj = (width * 4) + 4;
    let f = 4;

    for (let y = 0; y < height; y++) {
        let x;

        for (x = 1; x < m; x++, f += 4) {
            const L_x = LU[x - 1];
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
            const L_x = LU[x];
            img[f + 0] = L_x * (img[f + 0] - img[f + 4]);
            img[f + 1] = L_x * (img[f + 1] - img[f + 5]);
            img[f + 2] = L_x * (img[f + 2] - img[f + 6]);
            img[f + 3] = L_x * (img[f + 3] - img[f + 7]);
        }

        img[f + 0] = L_0 * (img[f + 0] - c * img[f + 4]);
        img[f + 1] = L_0 * (img[f + 1] - c * img[f + 5]);
        img[f + 2] = L_0 * (img[f + 2] - c * img[f + 6]);
        img[f + 3] = L_0 * (img[f + 3] - c * img[f + 7]);

        f += f_adj;
    }
}

function v_iconvolve_b2_ip(img, width, height, LU, m) {
    const c = LU[m];
    let n = m - 1;
    if (m >= height)
        n = m = height - 1;

    const L_0 = LU[0];
    const L_inf = LU[m - 1];
    const L_inf_mul = L_inf * c;
    const L_inf_div = LU[n] / c;
    const adj_width = width * 4;
    let pf = 0;
    let f = adj_width;

    let y;

    for (y = 1; y < m; y++) {
        const L_y = LU[y - 1];
        for (let x = 0; x < adj_width; x++, f++, pf++)
            img[f] -= L_y * img[pf];
    }

    for (; y < height - 1; y++)
        for (let x = 0; x < adj_width; x++, f++, pf++)
            img[f] -= L_inf * img[pf];

    for (let x = 0; x < adj_width; x++, f++, pf++)
        img[f] = L_inf_div * (img[f] - L_inf_mul * img[pf]);

    pf = f - 1;
    f = pf - adj_width;

    for (y = height - 2; y >= m; y--)
        for (let x = 0; x < adj_width; x++, f--, pf--)
            img[f] = L_inf * (img[f] - img[pf]);

    for (; y > 0; y--) {
        const L_y = LU[y];
        for (let x = 0; x < adj_width; x++, f--, pf--)
            img[f] = L_y * (img[f] - img[pf]);
    }

    for (let x = 0; x < adj_width; x++, f--, pf--)
        img[f] = L_0 * (img[f] - c * img[pf]);
}

function h_reconstruct_iconvolve_b2(src, src_width, height, dst_width, filter, window, norm, L, m) {
    if (dst_width > src_width) {
        if (src_width === 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const temp = new Float64Array(src);
        h_iconvolve_b2_ip(temp, src_width, height, L, m);
        return h_reconstruct(temp, src_width, height, dst_width, filter, window, norm);
    } else {
        if (dst_width === 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const ret = h_reconstruct(src, src_width, height, dst_width, filter, window, norm);
        h_iconvolve_b2_ip(ret, dst_width, height, L, m);
        return ret;
    }
}

function v_reconstruct_iconvolve_b2(src, width, src_height, dst_height, filter, window, norm, L, m) {
    if (dst_height > src_height) {
        if (src_height === 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const temp = new Float64Array(src);
        v_iconvolve_b2_ip(temp, width, src_height, L, m);
        return v_reconstruct(temp, width, src_height, dst_height, filter, window, norm);
    } else {
        if (dst_height === 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const ret = v_reconstruct(src, width, src_height, dst_height, filter, window, norm);
        v_iconvolve_b2_ip(ret, width, dst_height, L, m);
        return ret;
    }
}

/**
 * Resize an image using a reconstruction filter and an inverse discrete convolution with bandwidth 2.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {(x: number) => number} filter Reconstruction filter function.
 * @param {number} window Filter function window.
 * @param {number} norm Normalization constant.
 * @param {Float64Array} LU LU-decomposition matrix constants.
 * @param {number} m Number of matrix coefficients before convergence.
 * @param {boolean} nop Boolean flag for a no-op case.
 * @returns {Float64Array} Destination image in a 4-channel format.
 */
export function reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, filter, window, norm, LU, m, nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (nop) {
        if (dst_width === src_width) {
            if (dst_height === src_height)
                return new Float64Array(src);
            return v_reconstruct_iconvolve_b2(src, src_width, src_height, dst_height, filter, window, norm, LU, m);
        } else if (dst_height === src_height) {
            return h_reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, filter, window, norm, LU, m);
        }
    }

    const h_intermediate = dst_width * src_height;
    const v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        const temp = v_reconstruct_iconvolve_b2(src, src_width, src_height, dst_height, filter, window, norm, LU, m);
        return h_reconstruct_iconvolve_b2(temp, src_width, dst_height, dst_width, filter, window, norm, LU, m);
    } else {
        const temp = h_reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, filter, window, norm, LU, m);
        return v_reconstruct_iconvolve_b2(temp, dst_width, src_height, dst_height, filter, window, norm, LU, m);
    }
}

function h_filter_ip(img, width, height, bounds, coeffs, support, reverse) {
    const adj_width = width * 4;
    let bounds_ptr;
    let coeffs_ptr;

    const dst_change = reverse ? -4 : 4;
    let src_offset = 0;
    let dst_offset = reverse ? adj_width - 4 : 0;
    for (let y = 0; y < height; y++, src_offset += adj_width, dst_offset += adj_width) {
        coeffs_ptr = 0;
        bounds_ptr = 0;

        let dst_pixel = dst_offset;
        for (let x = 0; x < width; x++, bounds_ptr += 2, coeffs_ptr += support, dst_pixel += dst_change) {
            const min_x = bounds[bounds_ptr + 0];
            const max_x = bounds[bounds_ptr + 1];

            let src_pixel = src_offset + (min_x * 4);
            for (let s = min_x, i = 0; s <= max_x; s++, i++, src_pixel += 4) {
                const weight = coeffs[coeffs_ptr + i];

                if (dst_pixel === src_pixel) {
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

function v_filter_ip(img, width, height, bounds, coeffs, support, reverse) {
    const adj_width = width * 4;
    let bounds_ptr = 0;
    let coeffs_ptr = 0;

    const dst_change = reverse ? -adj_width : adj_width;
    let dst_offset = reverse ? adj_width * (height - 1) : 0;
    for (let y = 0; y < height; y++, bounds_ptr += 2, coeffs_ptr += support, dst_offset += dst_change) {
        const min_y = bounds[bounds_ptr + 0];
        const max_y = bounds[bounds_ptr + 1];

        let src_pixel = adj_width * min_y;
        for (let s = min_y, i = 0; s <= max_y; s++, i++) {
            const weight = coeffs[coeffs_ptr + i];

            let dst_pixel = dst_offset;
            for (let x = 0; x < width; x++, src_pixel += 4, dst_pixel += 4) {
                if (dst_pixel === src_pixel) {
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

function gen_iconvolve_filters(len, LU, m, b) {
    if (len < 2 * b - 1)
        return true;

    const L = (x, y) => LU[(x) * m + Math.min(y, m - 1)];
    const C = (x, y) => LU[((x) + m) * special + (y)];

    const special = b - 1;
    const cutoff = len - special;

    const forward_bounds = new Int32Array(2 * len);
    const backward_bounds = new Int32Array(2 * len);
    let support = special;
    const forward_coeffs = new Float64Array(support * len);
    const backward_coeffs = new Float64Array(b * len);

    let bounds_ptr = 0;
    let coeffs_ptr = 0;

    let z;

    forward_bounds[0] = 0;
    forward_bounds[1] = -1;
    bounds_ptr += 2;
    coeffs_ptr += support;

    for (z = 1; z < cutoff; z++, bounds_ptr += 2, coeffs_ptr += support) {
        const min_z = Math.max(z - support, 0);
        const max_z = z - 1;
        forward_bounds[bounds_ptr + 0] = min_z;
        forward_bounds[bounds_ptr + 1] = max_z;

        for (let s = max_z, i = 0; s >= min_z; s--, i++)
            forward_coeffs[coeffs_ptr + max_z - min_z - i] = -L(i, (z - 1) - i);
    }

    for (; z < len; z++, bounds_ptr += 2, coeffs_ptr += support) {
        const min_z = z - support;
        const max_z = z - 1;
        forward_bounds[bounds_ptr + 0] = min_z;
        forward_bounds[bounds_ptr + 1] = max_z;

        for (let s = max_z, i = 0; s >= min_z; s--, i++)
            forward_coeffs[coeffs_ptr + support - 1 - i] = -C(i, (len - 1) - z) * L(i, (z - 1) - i);
    }

    support = b;
    bounds_ptr = 0;
    coeffs_ptr = 0;

    for (z = len - 1; z >= cutoff; z--, bounds_ptr += 2, coeffs_ptr += support) {
        const M = L(special - 1, z) / C(special - 1, (len - 1) - z);
        const min_z = z;
        const max_z = Math.min(z + support - 1, len - 1);
        backward_bounds[bounds_ptr + 0] = min_z;
        backward_bounds[bounds_ptr + 1] = max_z;

        backward_coeffs[coeffs_ptr + 0] = M;
        for (let s = min_z + 1, i = 1; s <= max_z; s++, i++)
            backward_coeffs[coeffs_ptr + i] = -L(i - 1, z);
    }

    for (; z >= special; z--, bounds_ptr += 2, coeffs_ptr += support) {
        const M = L(special - 1, z);
        const min_z = z;
        const max_z = z + support - 1;
        backward_bounds[bounds_ptr + 0] = min_z;
        backward_bounds[bounds_ptr + 1] = max_z;

        backward_coeffs[coeffs_ptr + 0] = M;
        for (let s = min_z + 1, i = 1; s < max_z; s++, i++)
            backward_coeffs[coeffs_ptr + i] = -L(i - 1, z);
        backward_coeffs[coeffs_ptr + support - 1] = -M;
    }

    for (; z >= 0; z--, bounds_ptr += 2, coeffs_ptr += support) {
        const M = L(special - 1, z);
        const min_z = z;
        const max_z = z + support - 1;
        backward_bounds[bounds_ptr + 0] = min_z;
        backward_bounds[bounds_ptr + 1] = max_z;

        backward_coeffs[coeffs_ptr + 0] = M;
        for (let s = min_z + 1, i = 1; s < max_z; s++, i++)
            backward_coeffs[coeffs_ptr + i] = -C(i - 1, z) * L(i - 1, z);
        backward_coeffs[coeffs_ptr + support - 1] = -M * C(special - 1, z);
    }

    return { forward_bounds, forward_coeffs, backward_bounds, backward_coeffs };
}

function h_iconvolve_ip(img, width, height, LU, m, b) {
    const { forward_bounds, forward_coeffs, backward_bounds, backward_coeffs } = gen_iconvolve_filters(width, LU, m, b);
    h_filter_ip(img, width, height, forward_bounds, forward_coeffs, b - 1, false);
    h_filter_ip(img, width, height, backward_bounds, backward_coeffs, b, true);
}

function v_iconvolve_ip(img, width, height, LU, m, b) {
    const { forward_bounds, forward_coeffs, backward_bounds, backward_coeffs } = gen_iconvolve_filters(height, LU, m, b);
    v_filter_ip(img, width, height, forward_bounds, forward_coeffs, b - 1, false);
    v_filter_ip(img, width, height, backward_bounds, backward_coeffs, b, true);
}

function h_reconstruct_iconvolve(src, src_width, height, dst_width, filter, window, norm, LU, m, b) {
    if (src == null)
        return null;

    if (dst_width > src_width) {
        if (src_width === 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const temp = new Float64Array(src);
        if (h_iconvolve_ip(temp, src_width, height, LU, m, b))
            return null;
        return h_reconstruct(temp, src_width, height, dst_width, filter, window, norm);
    } else {
        if (dst_width === 1)
            return h_reconstruct(src, src_width, height, dst_width, filter, window, 1.0);

        const ret = h_reconstruct(src, src_width, height, dst_width, filter, window, norm);
        if (h_iconvolve_ip(ret, dst_width, height, LU, m, b))
            return null;
        return ret;
    }
}

function v_reconstruct_iconvolve(src, width, src_height, dst_height, filter, window, norm, LU, m, b) {
    if (src == null)
        return null;

    if (dst_height > src_height) {
        if (src_height === 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const temp = new Float64Array(src);
        if (v_iconvolve_ip(temp, width, src_height, LU, m, b))
            return null;
        return v_reconstruct(temp, width, src_height, dst_height, filter, window, norm);
    } else {
        if (dst_height === 1)
            return v_reconstruct(src, width, src_height, dst_height, filter, window, 1.0);

        const ret = v_reconstruct(src, width, src_height, dst_height, filter, window, norm);
        if (v_iconvolve_ip(ret, width, dst_height, LU, m, b))
            return null;
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
 * @param {Float64Array} LU LU-decomposition matrix constants.
 * @param {number} m Number of matrix coefficients before convergence.
 * @param {number} b Inverse discrete convolution bandwidth. Note that the domain in which the inverse convolution is performed must have each dimension be at least 2 * bandwidth - 1.
 * @param {boolean} nop Boolean flag for a no-op case.
 * @returns {Float64Array | null} Destination image in a 4-channel format, or null if iconvolve error.
 */
export function reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, filter, window, norm, LU, m, b, nop) {
    if (src_width <= 0 || src_height <= 0 || dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (nop) {
        if (dst_width === src_width) {
            if (dst_height === src_height)
                return new Float64Array(src);
            return v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, LU, m, b);
        } else if (dst_height === src_height) {
            return h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, LU, m, b);
        }
    }

    const h_intermediate = dst_width * src_height;
    const v_intermediate = dst_height * src_width;

    if (h_intermediate > v_intermediate) {
        const temp = v_reconstruct_iconvolve(src, src_width, src_height, dst_height, filter, window, norm, LU, m, b);
        return h_reconstruct_iconvolve(temp, src_width, dst_height, dst_width, filter, window, norm, LU, m, b);
    } else {
        const temp = h_reconstruct_iconvolve(src, src_width, src_height, dst_width, filter, window, norm, LU, m, b);
        return v_reconstruct_iconvolve(temp, dst_width, src_height, dst_height, filter, window, norm, LU, m, b);
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
    BSPLINE2: 5,
    BSPLINE3: 6,
    MITNET: 7,
    CATROM: 8,
    HAMMING3: 9,
    HAMMING4: 10,
    HAMMING8: 11,
    BSPLINE2I: 12,
    BSPLINE3I: 13,
    OMOMS3I: 14,
    OMOMS7I: 15,
    OMOMS11I: 16,
};

/**
 * Wrapper for `sample`, `scale`, `reconstruct`, and `reconstruct_iconvolve`.
 * @param {TypedArray} src Source image in a 4-channel format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {Filter} filter Resizing method (filter) to be used. Defaults to Mitchell-Netravali.
 * @returns {TypedArray | null} Destination image in a 4-channel format, or null if iconvolve error. Returns the same type as source image if `filter` is `NEAREST`, otherwise returns Float64Array.
 */
export function resize(src, src_width, src_height, dst_width, dst_height, filter) {
    switch (filter) {
        case Filter.NEAREST:
            return sample(src, src_width, src_height, dst_width, dst_height);
        case Filter.AREA:
            return scale(src, src_width, src_height, dst_width, dst_height);
        case Filter.TRIANGLE:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Triangle, 1.0, 1.0, true);
        case Filter.HERMITE:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hermite, 1.0, 1.0, true);
        case Filter.BSPLINE2:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.BSpline2, 1.5, 1.0, false);
        case Filter.BSPLINE3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.BSpline3, 2.0, 1.0, false);
        default:
        case Filter.MITNET:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.MitNet, 2.0, 1.0, false);
        case Filter.CATROM:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.CatRom, 2.0, 1.0, true);
        case Filter.HAMMING3:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hamming3, 3.0, 1.0, true);
        case Filter.HAMMING4:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hamming4, 4.0, 1.0, true);
        case Filter.HAMMING8:
            return reconstruct(src, src_width, src_height, dst_width, dst_height, Filters.Hamming8, 8.0, 1.0, true);
        case Filter.BSPLINE2I:
            return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, Filters.BSpline2, 1.5, 8.0, Filters.LU_bspline2i, 11, 1.1428571428571428, true);
        case Filter.BSPLINE3I:
            return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, Filters.BSpline3, 2.0, 6.0, Filters.LU_bspline3i, 14, 1.2, true);
        case Filter.OMOMS3I:
            return reconstruct_iconvolve_b2(src, src_width, src_height, dst_width, dst_height, Filters.OMOMS3, 2.0, 5.25, Filters.LU_omoms3, 18, 1.2352941176470589, true);
        case Filter.OMOMS7I:
            return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, Filters.OMOMS7, 4.0, 1952.8179190751446, Filters.LU_omoms7, 35, 4, true);
        case Filter.OMOMS11I:
            return reconstruct_iconvolve(src, src_width, src_height, dst_width, dst_height, Filters.OMOMS11, 6.0, 5148314.024847966, Filters.LU_omoms11, 54, 6, true);
    }
}

/** Simple class abstrating a pool of Web Workers. */
export class Resizer {
    /**
     * Create a resizer with a pool of Web Workers.
     * @param {number} pool Number of Web Workers in pool. Defaults to 4.
     */
    constructor(pool = 4) {
        this.workers = [];
        this.locks = Array(pool).fill(false);
        this.queue = [];
        for (let i = 0; i < pool; i++)
            this.workers.push(new Worker(import.meta.url, { type: "module" }));
    }

    /**
     * Wrapper for `sample`, `scale`, `reconstruct`, and `reconstruct_iconvolve`. Compute is done in a Web Worker.
     * @param {TypedArray} src Source image in a 4-channel format. Since the image is transferred to a Worker, it will become inaccessible. It is returned once the operation completes
     * @param {number} src_width Source image width.
     * @param {number} src_height Source image height.
     * @param {number} dst_width Destination image width.
     * @param {number} dst_height Destination image height.
     * @param {Filter} filter Resizing method (filter) to be used. Defaults to Mitchell-Netravali.
     * @returns {Promise<{src: TypedArray, dst: TypedArray | null}>} Destination image in a 4-channel format, or null if iconvolve error. Returns the same type as source image if `filter` is `NEAREST`, otherwise returns Float64Array.
     */
    resize(src, src_width, src_height, dst_width, dst_height, filter) {
        return new Promise((resolve) => {
            this.queue.push([[src, src_width, src_height, dst_width, dst_height, filter], resolve]);
            this.#check();
        });
    }

    #check() {
        if (this.queue.length === 0)
            return;
        for (let i = 0; i < this.workers.length; i++) {
            if (!this.locks[i]) {
                const [args, resolve] = this.queue.shift();
                this.locks[i] = true;
                this.workers[i].onmessage = (event) => {
                    this.locks[i] = false;
                    this.#check();
                    resolve(event.data);
                };
                this.workers[i].postMessage(args, [args[0].buffer]);
                return;
            }
        }
    }
}

// Code to be executed in a worker
if (typeof WorkerGlobalScope !== "undefined" && self instanceof WorkerGlobalScope) {
    onmessage = (event) => {
        const [src, src_width, src_height, dst_width, dst_height, filter] = event.data;
        const dst = resize(src, src_width, src_height, dst_width, dst_height, filter);
        postMessage({ src, dst }, dst?.buffer ? [src.buffer, dst.buffer] : [src.buffer]);
    };
}

// Copyright (c) 2024 silverslither.

import * as FilterFunctions from "./filters.js";

/**
 * Resample an image using nearest neighbor interpolation.
 * @param {TypedArray} src Source image in RGBA format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @returns {TypedArray} Destination image in RGBA format, with the same type as source image.
 */
export function sample(src, src_width, src_height, dst_width, dst_height) {
    if (dst_width <= 0 || dst_height <= 0)
        return new (Object.getPrototypeOf(src).constructor)();

    const x_factor = src_width / dst_width;
    const y_factor = src_height / dst_height;

    const dst = new (Object.getPrototypeOf(src).constructor)((dst_width * dst_height) << 2);

    let dst_pixel = 0;
    for (let y = 0; y < dst_height; y++) {
        const mapped_y = Math.ceil(y_factor * (y + 0.5) - 1.0) * src_width;
        for (let x = 0; x < dst_width; x++, dst_pixel += 4) {
            const mapped_x = Math.ceil(x_factor * (x + 0.5) - 1.0);
            const src_pixel = (mapped_y + mapped_x) << 2;
            dst[dst_pixel + 0] += src[src_pixel + 0];
            dst[dst_pixel + 1] += src[src_pixel + 1];
            dst[dst_pixel + 2] += src[src_pixel + 2];
            dst[dst_pixel + 3] += src[src_pixel + 3];
        }
    }

    return dst;
}

function hscale(src, src_width, height, dst_width) {
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
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }
            continue;
        }

        const low_x = Math.max(min_x - 1, 0) << 2;
        const high_x = Math.min(max_x, src_width - 1) << 2;
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
            for (let src_pixel = y + min_x; src_pixel < src_max; src_pixel += 4) {
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
                dst[dst_pixel + 0] += src[src_pixel + 0];
                dst[dst_pixel + 1] += src[src_pixel + 1];
                dst[dst_pixel + 2] += src[src_pixel + 2];
                dst[dst_pixel + 3] += src[src_pixel + 3];
            }
            continue;
        }

        const low_y = adj_width * Math.max(min_y - 1, 0);
        const high_y = adj_width * Math.min(max_y, src_height - 1);
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
            for (let src_pixel = x + min_y; src_pixel < src_max; src_pixel += adj_width) {
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
 * @param {TypedArray} src Source image in RGBA format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @returns {Float64Array} Destination image in RGBA format.
 */
export function scale(src, src_width, src_height, dst_width, dst_height) {
    if (dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    if (dst_width === src_width) {
        if (dst_height === src_height)
            return src;
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

function hfilter(src, src_width, height, dst_width, filter, window) {
    const adj_src_width = src_width << 2;
    const adj_dst_width = dst_width << 2;
    const adj_src_area = adj_src_width * height;
    const factor = src_width / dst_width;

    const dst = new Float64Array(height * adj_dst_width);

    const inv_filter_scale = Math.max(factor, 1.0);
    const filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    const width_end = adj_src_width - 4;

    let dst_offset = 0;
    for (let x = 0; x < dst_width; x++, dst_offset += 4) {
        const mapped_x = factor * (x + 0.5) - 0.5;
        let min_x = Math.ceil(mapped_x - window);
        let max_x = Math.floor(mapped_x + window);

        while (min_x < 0) {
            const weight = filter((mapped_x - (min_x++)) * filter_scale) * filter_scale;
            for (let src_pixel = 0, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        while (max_x >= src_width) {
            const weight = filter(((max_x--) - mapped_x) * filter_scale) * filter_scale;
            for (let src_pixel = width_end, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        let src_offset = min_x << 2;
        for (let s = min_x; s <= max_x; s++, src_offset += 4) {
            const weight = filter(Math.abs(mapped_x - s) * filter_scale) * filter_scale;
            for (let src_pixel = src_offset, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += adj_src_width, dst_pixel += adj_dst_width) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

function vfilter(src, width, src_height, dst_height, filter, window) {
    const adj_width = width << 2;
    const adj_src_area = adj_width * src_height;
    const factor = src_height / dst_height;

    const dst = new Float64Array(adj_width * dst_height);

    const inv_filter_scale = Math.max(factor, 1.0);
    const filter_scale = 1.0 / inv_filter_scale;
    window *= inv_filter_scale;

    const height_end = adj_width * (src_height - 1);

    let dst_offset = 0;
    for (let y = 0; y < dst_height; y++, dst_offset += adj_width) {
        const mapped_y = factor * (y + 0.5) - 0.5;
        let min_y = Math.ceil(mapped_y - window);
        let max_y = Math.floor(mapped_y + window);

        while (min_y < 0) {
            const weight = filter((mapped_y - (min_y++)) * filter_scale) * filter_scale;
            for (let src_pixel = 0, dst_pixel = dst_offset; src_pixel < adj_width; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        while (max_y >= src_height) {
            const weight = filter(((max_y--) - mapped_y) * filter_scale) * filter_scale;
            for (let src_pixel = height_end, dst_pixel = dst_offset; src_pixel < adj_src_area; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }

        let src_offset = adj_width * min_y;
        let src_max = src_offset + adj_width;
        for (let s = min_y; s <= max_y; s++, src_offset = src_max, src_max += adj_width) {
            const weight = filter(Math.abs(mapped_y - s) * filter_scale) * filter_scale;
            for (let src_pixel = src_offset, dst_pixel = dst_offset; src_pixel < src_max; src_pixel += 4, dst_pixel += 4) {
                dst[dst_pixel + 0] += src[src_pixel + 0] * weight;
                dst[dst_pixel + 1] += src[src_pixel + 1] * weight;
                dst[dst_pixel + 2] += src[src_pixel + 2] * weight;
                dst[dst_pixel + 3] += src[src_pixel + 3] * weight;
            }
        }
    }

    return dst;
}

/**
 * @enum {number}
 */
export const Filters = {
    DEFAULT: 0,
    NEAREST: 1,
    AREA: 2,
    TRIANGLE: 3,
    HERMITE: 4,
    B_SPLINE_2: 5,
    B_SPLINE_3: 6,
    KEYS_HALF: 7,
    MITNET: 8,
    MITNET_SHARP: 9,
    CATROM: 10,
    CATROM_SHARP: 11,
    LANCZOS_3: 12,
    LANCZOS_4: 13
};

/**
 * Resample an image using a reconstruction filter. Also acts as a wrapper for `sample` and `scale`.
 * @param {TypedArray} src Source image in RGBA format.
 * @param {number} src_width Source image width.
 * @param {number} src_height Source image height.
 * @param {number} dst_width Destination image width.
 * @param {number} dst_height Destination image height.
 * @param {Filters} filter Reconstruction filter to be used. `NEAREST` acts as a wrapper for `sample`, and `AREA` acts as a wrapper for `scale`. The default filter used is Mitchell-Netravali.
 * @returns {Float64Array} Destination image in RGBA format.
 */
export function resize(src, src_width, src_height, dst_width, dst_height, filter = 0) {
    if (dst_width <= 0 || dst_height <= 0)
        return new Float64Array();

    let filter_func, window, nop = false;

    switch (filter) {
        case Filters.NEAREST:
            return sample(src, src_width, src_height, dst_width, dst_height);
        case Filters.AREA:
            return scale(src, src_width, src_height, dst_width, dst_height);
        case Filters.TRIANGLE:
            filter_func = FilterFunctions.Triangle;
            window = 1.0;
            nop = true;
            break;
        case Filters.HERMITE:
            filter_func = FilterFunctions.Hermite;
            window = 1.0;
            nop = true;
            break;
        case Filters.B_SPLINE_2:
            filter_func = FilterFunctions.BSpline2;
            window = 1.5;
            break;
        case Filters.B_SPLINE_3:
            filter_func = FilterFunctions.BSpline3;
            window = 2.0;
            break;
        case Filters.KEYS_HALF:
            filter_func = FilterFunctions.KeysHalf;
            window = 2.0;
            break;
        default:
        case Filters.MITNET:
            filter_func = FilterFunctions.MitNet;
            window = 2.0;
            break;
        case Filters.MITNET_SHARP:
            filter_func = FilterFunctions.MitNetSharp;
            window = 2.0;
            break;
        case Filters.CATROM:
            filter_func = FilterFunctions.CatRom;
            window = 2.0;
            nop = true;
            break;
        case Filters.CATROM_SHARP:
            filter_func = FilterFunctions.CatRomSharp;
            window = 2.0;
            nop = true;
            break;
        case Filters.LANCZOS_3:
            filter_func = FilterFunctions.Lanczos3;
            window = 3.0;
            nop = true;
            break;
        case Filters.LANCZOS_4:
            filter_func = FilterFunctions.Lanczos4;
            window = 4.0;
            nop = true;
            break;
    }

    if (nop) {
        if (dst_width === src_width) {
            if (dst_height === src_height)
                return src;
            return vfilter(src, src_width, src_height, dst_height, filter_func, window);
        } else if (dst_height === src_height) {
            return hfilter(src, src_width, src_height, dst_width, filter_func, window);
        }
    }

    const x_factor = dst_width / src_width;
    const y_factor = dst_height / src_height;

    if (x_factor >= y_factor) {
        const temp = vfilter(src, src_width, src_height, dst_height, filter_func, window);
        return hfilter(temp, src_width, dst_height, dst_width, filter_func, window);
    } else {
        const temp = hfilter(src, src_width, src_height, dst_width, filter_func, window);
        return vfilter(temp, dst_width, src_height, dst_height, filter_func, window);
    }
}

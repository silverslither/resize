// Copyright (c) 2024-2025 silverslither.

import { mul_alpha, div_alpha, srgb_encode, srgb_decode, get_sigmoidization_params, sigmoidal_contrast_increase, sigmoidal_contrast_decrease } from "./voir/colour.js";
import { Filter, Voir } from "./voir/voir.js";
const voir = new Voir();

let input, width, height, filterName, smoothFilterName, gradientFilterName, gmultiplier, rawGradientBox, linearizeBox, beta, submit, err;

document.addEventListener("DOMContentLoaded", () => {
    input = document.getElementById("input");
    width = document.getElementById("width");
    height = document.getElementById("height");
    filterName = document.getElementById("filter");
    smoothFilterName = document.getElementById("smooth-filter");
    gradientFilterName = document.getElementById("gradient-filter");
    gmultiplier = document.getElementById("gradient-multiplier");
    rawGradientBox = document.getElementById("raw-gradient");
    linearizeBox = document.getElementById("linearize");
    beta = document.getElementById("beta");
    submit = document.querySelector("button");
    err = document.getElementById("err");
    submit.addEventListener("click", listener);
    input.value = "";
});

async function listener() {
    try {
        err.innerText = "";
        const file = input.files[0];
        if (file == null) {
            err.innerText += "error: no file selected";
            return;
        }
        await main(file);
    } catch (e) {
        err.innerText += e;
    }
}

async function main(file) {
    const src = await decode(file);
    src.data = new Float64Array(src.data);

    const dst_width = parseDimension(width.value, src.width);
    if (dst_width === -1) {
        err.innerText += `error: invalid width '${width.value}'\n`;
        return;
    }
    const dst_height = parseDimension(height.value, src.height);
    if (dst_height === -1) {
        err.innerText += `error: invalid height '${height.value}'\n`;
        return;
    }
    const filter = parseFilter(filterName.value);
    let smoothFilter = parseFilter(smoothFilterName.value);
    let gradientFilter = parseFilter(gradientFilterName.value);
    const haloMinimize = smoothFilter !== Filter.DEFAULT;
    const gradientMultiplier = parseGradientMultiplier(gmultiplier.value);
    const rawGradient = rawGradientBox.checked;
    const linearize = linearizeBox.checked;
    const sigParams = parseSigmoidizationBeta(beta.value);

    preprocess(src.data, linearize, sigParams);
    mul_alpha(src.data);

    if (smoothFilter == Filter.DEFAULT)
        smoothFilter = filter;
    if (gradientFilter == Filter.DEFAULT)
        gradientFilter = smoothFilter;

    let dst;
    if (rawGradient) {
        const { dst: temp } = await voir.resize(src.data, src.width, src.height, dst_width, dst_height, gradientFilter);
        const xmult = gradientMultiplier * Math.max(dst_width / src.width, 1.0);
        const ymult = gradientMultiplier * Math.max(dst_height / src.height, 1.0);
        ({ dst } = await gradientMagnitude(temp, dst_width, dst_height, xmult, ymult, false));
    } else if (haloMinimize) {
        ({ dst } = await haloMinimizedResize(src.data, src.width, src.height, dst_width, dst_height, filter, smoothFilter, gradientFilter, gradientMultiplier));
    } else {
        ({ dst } = await voir.resize(src.data, src.width, src.height, dst_width, dst_height, filter));
    }

    div_alpha(dst);
    postprocess(dst, linearize, sigParams);

    writeFileSync(`RESIZED_${file.name}`, await encode(dst, dst_width, dst_height));
}

function preprocess(arr, linearize, params) {
    for (let i = 0; i < arr.length; i++) {
        let temp = 0.00392156862745098 * arr[i];
        if (i % 4 != 3) {
            if (linearize)
                temp = srgb_decode(temp);
            if (params != null)
                temp = sigmoidal_contrast_decrease(temp, params);
        }
        arr[i] = temp;
    }
}

function postprocess(arr, linearize, params) {
    for (let i = 0; i < arr.length; i++) {
        let temp = arr[i];
        if (i % 4 != 3) {
            if (params != null)
                temp = sigmoidal_contrast_increase(temp, params);
            if (linearize)
                temp = srgb_encode(temp);
        }
        arr[i] = 255.0 * temp;
    }
}

function parseDimension(str, src_dimension) {
    let num = Number(str);
    if (num !== num) {
        num = Number(str.slice(0, -1));
        if (num !== num)
            return -1;
        const last = str.at(-1);
        if (last === "%")
            num *= 0.01 * src_dimension;
        else if (last.toLowerCase() === "x")
            num *= src_dimension;
        else
            return -1;
    }

    if (num <= 0.5 || num >= 2147483647.5)
        return -1;

    return Number.isInteger(num - 0.5) ? 2 * Math.round(0.5 * num) : Math.round(num);
}

function parseFilter(str) {
    if (str === "")
        return 0;

    const filters = Object.keys(Filter);
    str = str.toUpperCase();

    for (let i = 0; i < filters.length; i++)
        if (str === filters[i])
            return Filter[filters[i]];

    err.innerText += `warning: invalid filter '${str}'\n`;

    return 0;
}

function parseSigmoidizationBeta(str) {
    if (str === "")
        return null;

    const beta = Number(str);
    if (beta !== beta || beta <= 0) {
        err.innerText += `warning: invalid sigmoidization contrast '${str}'\n`;
        return null;
    }

    return get_sigmoidization_params(beta);
}

function parseGradientMultiplier(str) {
    if (str === "")
        return 2.0;

    const mult = Number(str);
    if (mult !== mult || mult <= 0) {
        err.innerText += `warning: invalid gradient magnitude multiplier '${str}'\n`;
        return 2.0;
    }

    return mult;
}

async function haloMinimizedResize(src, src_width, src_height, dst_width, dst_height, sharpFilter, smoothFilter, gradientFilter, multiplier) {
    let sharp, smooth;
    ({ src, dst: sharp } = await voir.resize(src, src_width, src_height, dst_width, dst_height, sharpFilter));
    ({ src, dst: smooth } = await voir.resize(src, src_width, src_height, dst_width, dst_height, smoothFilter));

    const xmult = Math.max(dst_width / src_width, 1.0);
    const ymult = Math.max(dst_height / src_height, 1.0);
    let gradient;
    if (gradientFilter === sharpFilter) {
        ({ src: sharp, dst: gradient } = await gradientMagnitude(sharp, dst_width, dst_height, xmult, ymult));
    } else if (gradientFilter === smoothFilter) {
        ({ src: smooth, dst: gradient } = await gradientMagnitude(smooth, dst_width, dst_height, xmult, ymult));
    } else {
        let temp;
        ({ src, dst: temp } = await voir.resize(src, src_width, src_height, dst_width, dst_height, gradientFilter));
        ({ dst: gradient } = await gradientMagnitude(temp, dst_width, dst_height, xmult, ymult));
    }

    const length = dst_width * dst_height << 2;
    for (let i = 0; i < length; i++) {
        const c = multiplier * gradient[i];
        gradient[i] = c * sharp[i] + (1.0 - c) * smooth[i];
    }

    return { src, dst: gradient };
}

async function gradientMagnitude(src, width, height, xmult, ymult, alpha = true) {
    const CDiffKernel = new Float64Array([0.5, 0, -0.5]);

    let Gx, Gy;
    ({ src, dst: Gx } = await voir.convolve(src, width, height, CDiffKernel, null, 3, 3, 0.0, 0.0));
    ({ src, dst: Gy } = await voir.convolve(src, width, height, null, CDiffKernel, 3, 3, 0.0, 0.0));

    xmult *= xmult;
    ymult *= ymult;
    const length = width * height << 2;
    if (alpha) {
        for (let i = 0; i < length; i++)
            Gx[i] = Math.sqrt(xmult * Gx[i] * Gx[i] + ymult * Gy[i] * Gy[i]);
    } else {
        for (let i = 0; i < length; i++) {
            if (i % 4 == 3) {
                Gx[i] = 1.0;
                continue;
            }
            Gx[i] = Math.sqrt(xmult * Gx[i] * Gx[i] + ymult * Gy[i] * Gy[i]);
        }
    }

    return { src, dst: Gx };
}

function decode(file) {
    return new Promise((resolve) => {
        const reader = new FileReader();
        reader.addEventListener("load", async () => {
            const image = new Image();
            image.src = reader.result;
            await image.decode();
            const canvas = new OffscreenCanvas(image.naturalWidth, image.naturalHeight);
            const context = canvas.getContext("2d");
            context.drawImage(image, 0, 0);
            const data = context.getImageData(0, 0, image.naturalWidth, image.naturalHeight);
            resolve({ data: data.data, width: data.width, height: data.height });
        });
        reader.readAsDataURL(file);
    });
}

function encode(dataArray, width, height) {
    if (!(dataArray instanceof Uint8ClampedArray))
        dataArray = new Uint8ClampedArray(dataArray);
    const data = new ImageData(dataArray, width, height);
    const canvas = new OffscreenCanvas(width, height);
    const context = canvas.getContext("2d");
    context.putImageData(data, 0, 0);
    return canvas.convertToBlob();
}

function writeFileSync(name, blob) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
}

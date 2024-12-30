// Copyright (c) 2024 silverslither.

import { resize } from "./resize.js";
import { mul_alpha, div_alpha, srgb_encode, srgb_decode, get_sigmoidization_params, sigmoidal_contrast_increase, sigmoidal_contrast_decrease } from "./colour.js";
import { encode, decode } from "https://cdn.jsdelivr.net/npm/fast-png@6.2.0/+esm"

let input, output, width, height, filterName, linearizeBox, beta, submit, err;

document.addEventListener("DOMContentLoaded", () => {
    input = document.getElementById("input");
    output = document.querySelector("a");
    width = document.getElementById("width");
    height = document.getElementById("height");
    filterName = document.getElementById("filter");
    linearizeBox = document.getElementById("linearize");
    beta = document.getElementById("beta");
    submit = document.querySelector("button");
    err = document.getElementById("err");
    submit.addEventListener("click", main);
    input.value = "";
});

async function main() {
    try {
        err.innerText = "";
        const _file = input.files[0];
        if (_file == null) {
            err.innerText += "error: no file selected";
            return;
        }

        const src = expandChannels(decode(await _file.arrayBuffer()));

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
        const sigParams = parseSigmoidizationBeta(beta.value);
        const linearize = linearizeBox.checked;

        preprocess(src.data, linearize, sigParams);
        mul_alpha(src.data);

        const dst = resize(src.data, src.width, src.height, dst_width, dst_height, filter);

        div_alpha(dst);
        postprocess(dst, linearize, sigParams);

        writeFileSync(`RESIZED_${_file.name}`, encodeArray(dst, dst_width, dst_height));
    } catch (e) {
        err.innerText += e;
    }
}

function expandChannels(img) {
    if (img.channels === 4)
        return new Float64Array(img);
    const area = img.width * img.height;
    const newData = new Float64Array(area << 2);
    for (let i = 0; i < area; i++) {
        const offset_4 = i << 2;
        const offset_3 = i * 3;
        newData[offset_4 + 0] = img.data[offset_3 + 0];
        newData[offset_4 + 1] = img.data[offset_3 + 1];
        newData[offset_4 + 2] = img.data[offset_3 + 2];
        newData[offset_4 + 3] = 255.0;
    }
    return { data: newData, width: img.width, height: img.height };
}

function preprocess(arr, linearize, params) {
    for (let i = 0; i < arr.length; i++) {
        let temp = 0.00392156862745098 * arr[i];
        if (i % 4 != 3) {
            if (linearize)
                temp = srgb_decode(temp);
            if (params != null)
                temp = sigmoidal_contrast_decrease(temp, params)
        }
        arr[i] = temp;
    }
}

function postprocess(arr, linearize, params) {
    for (let i = 0; i < arr.length; i++) {
        let temp = arr[i];
        if (i % 4 != 3) {
            if (params != null)
                temp = sigmoidal_contrast_increase(temp, params)
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

    const filters = [
        "DEFAULT",
        "NEAREST",
        "AREA",
        "TRIANGLE",
        "HERMITE",
        "LAGRANGE2",
        "LAGRANGE3",
        "BSPLINE2",
        "BSPLINE3",
        "MITNET",
        "CATROM",
        "MKS2013",
        "LANCZOS3",
        "LANCZOS4",
        "HAMMING3",
        "HAMMING4",
        "BSPLINE3I",
        "OMOMS3I"
    ];
    str = str.toUpperCase();

    for (let i = 0; i < filters.length; i++)
        if (str === filters[i])
            return i;

    err.innerText += `warning: invalid filter '${str}'\n`;

    return 0;
}

function parseSigmoidizationBeta(str) {
    const beta = Number(str);
    if (beta !== beta || beta < 0) {
        err.innerText += `warning: invalid sigmoidization contrast '${str}'\n`;
        return null;
    }

    if (beta === 0)
        return null;

    return get_sigmoidization_params(beta);
}

function encodeArray(typedArray, width, height) {
    if (!(typedArray instanceof Uint8ClampedArray) && !(typedArray instanceof Uint8Array))
        typedArray = new Uint8ClampedArray(typedArray);
    return new Blob([encode({ data: typedArray, width, height })]);
}

function writeFileSync(name, blob) {
    output.href = URL.createObjectURL(blob);
    output.download = name;
    output.click();
}

// Copyright (c) 2024-2025 silverslither.

import { Filter, resize } from "./resize.js";
import { mul_alpha, div_alpha, srgb_encode, srgb_decode, get_sigmoidization_params, sigmoidal_contrast_increase, sigmoidal_contrast_decrease } from "./colour.js";

let input, width, height, filterName, linearizeBox, beta, submit, err;

document.addEventListener("DOMContentLoaded", () => {
    input = document.getElementById("input");
    width = document.getElementById("width");
    height = document.getElementById("height");
    filterName = document.getElementById("filter");
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
    const sigParams = parseSigmoidizationBeta(beta.value);
    const linearize = linearizeBox.checked;

    preprocess(src.data, linearize, sigParams);
    mul_alpha(src.data);

    const dst = resize(src.data, src.width, src.height, dst_width, dst_height, filter);

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

    const filters = Object.keys(Filter);
    str = str.toUpperCase();

    for (let i = 0; i < filters.length; i++)
        if (str === filters[i])
            return Filter[filters[i]];

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

async function decode(file) {
    return new Promise((resolve) => {
        const reader = new FileReader();
        reader.addEventListener("load", async () => {
            const image = new Image();
            image.src = reader.result;
            await image.decode();
            const canvas = document.createElement("canvas");
            const context = canvas.getContext("2d");
            canvas.width = image.naturalWidth;
            canvas.height = image.naturalHeight;
            context.drawImage(image, 0, 0);
            const data = context.getImageData(0, 0, image.naturalWidth, image.naturalHeight);
            resolve({ data: data.data, width: data.width, height: data.height });
        });
        reader.readAsDataURL(file);
    });
}

async function encode(dataArray, width, height) {
    return new Promise((resolve, reject) => {
        if (!(dataArray instanceof Uint8ClampedArray))
            dataArray = new Uint8ClampedArray(dataArray);
        const data = new ImageData(dataArray, width, height);
        const canvas = document.createElement("canvas");
        const context = canvas.getContext("2d");
        canvas.width = width;
        canvas.height = height;
        context.putImageData(data, 0, 0);
        canvas.toBlob((blob) => {
            if (blob == null)
                reject();
            else
                resolve(blob);
        });
    });
}

function writeFileSync(name, blob) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
}

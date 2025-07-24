// Copyright (c) 2024-2025 silverslither.

import { Filter, Voir } from "./voir/voir.js";
const pica = window.pica();
const voir = new Voir();

let input, width, height, picaFilter, voirFilter, submit, err;
let lock = false, start = NaN;

document.addEventListener("DOMContentLoaded", () => {
    input = document.getElementById("input");
    width = document.getElementById("width");
    height = document.getElementById("height");
    picaFilter = document.getElementById("pica-filter");
    voirFilter = document.getElementById("voir-filter");
    submit = document.querySelector("button");
    err = document.getElementById("err");
    submit.addEventListener("click", listener);
    input.value = "";
});

async function listener() {
    if (lock)
        return;
    lock = true;

    try {
        err.innerText = "";
        const file = input.files[0];
        if (file == null) {
            err.innerText += "error: no file selected";
            lock = false;
            return;
        }
        await main(file);
    } catch (e) {
        err.innerText += e;
    }

    lock = false;
}

async function main(file) {
    const src = await decode(file);

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

    time();
    await pica.resize(src.canvas, new OffscreenCanvas(dst_width, dst_height), { filter: picaFilter.value || "mks2013" });
    err.innerText += `pica: ${time()} ms\n`;

    time();
    await voir.resize(src.data, src.width, src.height, dst_width, dst_height, parseFilter(voirFilter.value));
    err.innerText += `voir: ${time()} ms\n`;
}

function time() {
    if (start !== start) {
        start = performance.now();
    } else {
        const result = Math.round(performance.now() - start);
        start = NaN;
        return result;
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
            resolve({ canvas, data: data.data, width: data.width, height: data.height });
        });
        reader.readAsDataURL(file);
    });
}

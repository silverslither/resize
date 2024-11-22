// Copyright (c) 2024 silverslither.

import { resize, FILTERS } from "./resize.js";
import { encode, decode } from "https://cdn.jsdelivr.net/npm/fast-png@6.2.0/+esm"

main();

async function main() {
    const dst_width = 4500;
    const dst_height = 1500;

    const src = expandChannels(await getImage("in.png"));
    console.log(src);
    const dst = resize(src.data, src.width, src.height, dst_width, dst_height, FILTERS.MITNET_SHARP);

    writeFileSync("mine.png", encodeArray(dst, dst_width, dst_height));
}

async function getImage(file) {
    const res = await fetch(file);
    return decode(await res.arrayBuffer());
}

function expandChannels(img) {
    if (img.channels === 4)
        return img;
    const area = img.width * img.height;
    const newData = new Uint8Array(area << 2);
    for (let i = 0; i < area; i++) {
        newData[(i << 2) + 0] = img.data[3 * i + 0];
        newData[(i << 2) + 1] = img.data[3 * i + 1];
        newData[(i << 2) + 2] = img.data[3 * i + 2];
        newData[(i << 2) + 3] = 255;
    }
    return { data: newData, width: img.width, height: img.height };
}

function encodeArray(typedArray, width, height) {
    if (!(typedArray instanceof Uint8ClampedArray) && !(typedArray instanceof Uint8Array))
        typedArray = new Uint8ClampedArray(typedArray);
    return new Blob([encode({ data: typedArray, width, height })]);
}

function writeFileSync(name, blob) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
}

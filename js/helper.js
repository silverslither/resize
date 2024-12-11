function multiplyAlpha(arr, len) {
    for (let i = 0; i < len; i += 4) {
        arr[i + 0] *= arr[i + 3];
        arr[i + 1] *= arr[i + 3];
        arr[i + 2] *= arr[i + 3];
    }
}

function divideAlpha(arr, len) {
    for (let i = 0; i < len; i += 4) {
        if (arr[i + 3] === 0) {
            arr += 4;
            continue;
        }
        arr[i + 0] /= arr[i + 3];
        arr[i + 1] /= arr[i + 3];
        arr[i + 2] /= arr[i + 3];
    }
}

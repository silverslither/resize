export function multiplyAlpha(arr) {
    for (let i = 0; i < arr.length; i += 4) {
        arr[i + 0] *= arr[i + 3];
        arr[i + 1] *= arr[i + 3];
        arr[i + 2] *= arr[i + 3];
    }
}

export function divideAlpha(arr) {
    for (let i = 0; i < arr.length; i += 4) {
        if (arr[i + 3] !== 0) {
            arr[i + 0] /= arr[i + 3];
            arr[i + 1] /= arr[i + 3];
            arr[i + 2] /= arr[i + 3];
        }
    }
}

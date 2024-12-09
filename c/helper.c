#include "helper.h"

void multiplyAlpha(double *arr, size_t len) {
    const double *end = arr + len;
    while (arr < end) {
        arr[0] *= arr[3];
        arr[1] *= arr[3];
        arr[2] *= arr[3];
        arr += 4;
    }
}

void divideAlpha(double *arr, size_t len) {
    const double *end = arr + len;
    while (arr < end) {
        if (arr[3] == 0)
            continue;
        arr[0] /= arr[3];
        arr[1] /= arr[3];
        arr[2] /= arr[3];
        arr += 4;
    }
}

export function mul_alpha(arr) {
    for (let i = 0; i < arr.length; i += 4) {
        arr[i + 0] *= arr[i + 3];
        arr[i + 1] *= arr[i + 3];
        arr[i + 2] *= arr[i + 3];
    }
}

export function div_alpha(arr) {
    for (let i = 0; i < arr.length; i += 4) {
        if (arr[i + 3] !== 0) {
            arr[i + 0] /= arr[i + 3];
            arr[i + 1] /= arr[i + 3];
            arr[i + 2] /= arr[i + 3];
        }
    }
}

export function srgb_encode(x) {
    if (x > 0.003130668442500634)
        return 1.055 * x ** 0.4166666666666667 - 0.055;
    return 12.92 * x;
}

export function srgb_decode(x) {
    if (x > 0.04044823627710819)
        return (0.9478672985781991 * (x + 0.055)) ** 2.4;
    return 0.07739938080495357 * x;
}

export function get_sigmoidization_params(beta) {
    beta *= 1.4426950408889634;
    const params = {};
    const exp_beta = 2.0 ** beta;
    const sqrt_exp_beta = Math.sqrt(exp_beta);
    const inv_denominator = 1.0 / (sqrt_exp_beta - 1.0);
    params.beta = beta;
    params.inv_beta = 1.0 / beta;
    params.a = sqrt_exp_beta * inv_denominator;
    params.b = (exp_beta + sqrt_exp_beta) * inv_denominator;
    params.c = sqrt_exp_beta;
    return params;
}

export function sigmoidal_contrast_increase(x, params) {
    return params.a - params.b / (2.0 ** (params.beta * x) + params.c);
}

export function sigmoidal_contrast_decrease(x, params) {
    return Math.log2(params.b / (params.a - x) - params.c) * params.inv_beta;
}

// Copyright (c) 2024 silverslither.

export const Triangle = (x) => 1.0 - x;

export const Hermite = (x) => {
    return 1.0 - x * x * (3.0 - 2.0 * x);
};

export const BSpline2 = (x) => {
    if (x <= 0.5)
        return 0.75 - x * x;
    const x_ = x - 1.5;
    return 0.5 * x_ * x_;
};

export const BSpline3 = (x) => {
    if (x <= 1.0)
        return 0.6666666666666666 - x * x * (1.0 - 0.5 * x);
    const x_ = x - 2;
    return -0.16666666666666666 * x_ * x_ * x_;
};

export const KeysHalf = (x) => {
    if (x <= 1.0)
        return 0.8333333333333333 - x * x * (1.75 - x);
    return 1.6666666666666667 - x * (3.0 - x * (1.75 - 0.3333333333333333 * x));
};

export const MitNet = (x) => {
    if (x <= 1.0)
        return 0.8888888888888888 - x * x * (2.0 - 1.1666666666666667 * x);
    return 1.7777777777777777 - x * (3.3333333333333333 - x * (2.0 - 0.3888888888888889 * x));
};

export const MitNetSharp = (x) => {
    if (x <= 1.0)
        return 0.9166666666666666 - x * x * (2.125 - 1.25 * x);
    return 1.8333333333333333 - x * (3.5 - x * (2.125 - 0.4166666666666667 * x));
}

export const CatRom = (x) => {
    if (x <= 1.0)
        return 1.0 - x * x * (2.5 - 1.5 * x);
    return 2.0 - x * (4.0 - x * (2.5 - 0.5 * x));
};

export const CatRomSharp = (x) => {
    if (x <= 1.0)
        return 1.0 - x * x * (2.4 - 1.4 * x);
    return 2.4 - x * (4.8 - x * (3 - 0.6 * x));
};

export const MagicKernelSharp2013 = (x) => {
    if (x <= 0.5)
        return 1.0625 - 1.75 * x * x;
    if (x <= 1.5)
        return 1.75 - x * (2.75 - x);
    const x_ = x - 2.5;
    return -0.125 * x_ * x_;
};

const sinc3 = (x) => {
    const x2 = x * x;
    const v = 0.19010152698956836 - 0.00932297307587116 * x2;
    return 1.0471975511965979 - v * x2;
};
const lnorm3 = (x) => {
    const c0 = x;
    const c0_ = 3.0 - x;
    const c1 = x >= 0.5 ? 2.0 - x : x + 1.0;
    const c1_ = 3.0 - c1;
    const c2 = 1.0 - x;
    const c2_ = x + 2.0;
    const o0 = sinc3(c0);
    const o1 = sinc3(c1);
    const o2 = sinc3(c2);
    let v = 0.0;
    v -= o2 / c2;
    v += o1 / c1;
    v -= o0 / c0;
    v -= c0 * o0 / (c0_ * c0_);
    v += c1 * o1 / (c1_ * c1_);
    v -= c2 * o2 / (c2_ * c2_);
    return v;
};
export const Lanczos3 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const sign = ((x & 1) - 1) | 1;
    const poly_x = x > 1.5 ? 3.0 - x : x;
    const numerator = sinc3(poly_x) * poly_x / (x * x);
    const denominator = lnorm3(x % 1.0);
    return sign * numerator / denominator;
};

const sinc4 = (x) => {
    const x2 = x * x;
    const v = 0.08019908169872415 - 0.002212385212340519 * x2;
    return 0.7853981633974483 - v * x2;
};
const lnorm4 = (x) => {
    const c0 = x;
    const c0_ = 4.0 - x;
    const c1 = x + 1.0;
    const c1_ = 3.0 - x;
    const c2 = 2.0 - x;
    const c2_ = x + 2.0;
    const c3 = 1.0 - x;
    const c3_ = x + 3.0;
    const o0 = sinc4(c0);
    const o1 = sinc4(c1);
    const o2 = sinc4(c2);
    const o3 = sinc4(c3);
    let v = 0.0;
    v -= o3 / c3;
    v += o2 / c2;
    v -= o0 / c0;
    v += o1 / c1;
    v -= c1 * o1 / (c1_ * c1_);
    v += c0 * o0 / (c0_ * c0_);
    v -= c2 * o2 / (c2_ * c2_);
    v += c3 * o3 / (c3_ * c3_);
    return v;
};
export const Lanczos4 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const sign = ((x & 1) - 1) | 1;
    const poly_x = x > 2.0 ? 4.0 - x : x;
    const numerator = sinc4(poly_x) * poly_x / (x * x);
    const denominator = lnorm4(x % 1.0);
    return sign * numerator / denominator;
};

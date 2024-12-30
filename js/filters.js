// Copyright (c) 2024 silverslither.

export const Triangle = (x) => {
    1.0 - x;
};

export const Hermite = (x) => {
    return 1.0 - x * x * (3.0 - 2.0 * x);
};

export const Lagrange2 = (x) => {
    if (x <= 0.5)
        return 1.0 - x * x;
    return 1.0 - x * (1.5 - 0.5 * x);
}

// Normalized to 3.0
export const Lagrange3 = (x) => {
    if (x <= 1.0)
        return 3.0 - x * (1.5 + x * (3.0 - 1.5 * x));
    return 3.0 - x * (5.5 - x * (3.0 - 0.5 * x));
}

export const BSpline2 = (x) => {
    if (x <= 0.5)
        return 0.75 - x * x;
    const x_ = x - 1.5;
    return 0.5 * x_ * x_;
};

// Normalized to 6.0
export const BSpline3 = (x) => {
    if (x <= 1.0)
        return 4.0 - x * x * (6.0 - 3.0 * x);
    const x_ = 2.0 - x;
    return x_ * x_ * x_;
}

// Normalized to 9.0
export const MitNet = (x) => {
    if (x <= 1.0)
        return 8.0 - x * x * (18.0 - 10.5 * x);
    return 16.0 - x * (30.0 - x * (18.0 - 3.5 * x));
}

export const CatRom = (x) => {
    if (x <= 1.0)
        return 1.0 - x * x * (2.5 - 1.5 * x);
    return 2.0 - x * (4.0 - x * (2.5 - 0.5 * x));
};

export const MKS2013 = (x) => {
    if (x <= 0.5)
        return 1.0625 - 1.75 * x * x;
    if (x <= 1.5)
        return 1.75 - x * (2.75 - x);
    const x_ = x - 2.5;
    return -0.125 * x_ * x_;
};

const lsin3 = (x) => {
    const x2 = x * x;
    const v = 2.829828552115177 - 1.2490259408560183 * x2;
    return x * (1.7320508075688772 - v * x2);
};
const lsin3w = (x) => {
    const x2 = x * x;
    const v = 0.10480846489315472 - 0.005140024447967154 * x2;
    return x * (0.5773502691896257 - v * x2);
};
export const Lanczos3 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const round_x = Math.round(x);
    const sign = (round_x << 31 >> 31) | 1;
    const poly_x = sign * (x - round_x);
    return lsin3(poly_x) * lsin3w(x > 1.5 ? 3.0 - x : x) / (x * x);
};

const lsin4 = (x) => {
    const x2 = x * x;
    const v = 3.2676045526483732 - 1.4422509263560956 * x2;
    return x * (2.0 - v * x2);
};
const lsin4w = (x) => {
    const x2 = x * x;
    const v = 0.05105632113513083 - 0.0014084481702696247 * x2;
    return x * (0.5 - v * x2);
};
export const Lanczos4 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const round_x = Math.round(x);
    const sign = (round_x << 31 >> 31) | 1;
    const poly_x = sign * (x - round_x);
    return lsin4(poly_x) * lsin4w(x > 2.0 ? 4.0 - x : x) / (x * x);
};

const hsin = (x) => {
    const x2 = x * x;
    const v = 1.6338022763241866 - 0.7211254631780478 * x2;
    return x * (1.0 - v * x2);
};

const hcos3w = (x) => {
    const x2 = x * x;
    const v = 0.24920390748853422 - 0.019569144068978167 * x2;
    return 0.46164 - v * x2;
};
export const Hamming3 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const round_x = Math.round(x);
    const sign = (round_x << 31 >> 31) | 1;
    const poly_x = sign * (x - round_x);
    const flip = x > 1.5;
    const w = hcos3w(flip ? x - 3.0 : x);
    return hsin(poly_x) * (0.53836 + (flip ? -w : w)) / x;
};

const hcos4w = (x) => {
    const x2 = x * x;
    const v = 0.1401771979623005 - 0.006191799490575123 * x2;
    return 0.46164 - v * x2;
};
export const Hamming4 = (x) => {
    if (x < 1.0e-8)
        return 1.0;
    const round_x = Math.round(x);
    const sign = (round_x << 31 >> 31) | 1;
    const poly_x = sign * (x - round_x);
    const flip = x > 2.0;
    const w = hcos4w(flip ? x - 4.0 : x);
    return hsin(poly_x) * (0.53836 + (flip ? -w : w)) / x;
};

export const L_bspline3i = [
    0.20833333333333334,
    0.26666666666666666,
    0.26785714285714285,
    0.2679425837320574,
    0.26794871794871794,
    0.2679491583648231,
    0.26794918998527245,
    0.26794919225551855,
    0.2679491924185149,
    0.26794919243021753,
    0.2679491924310577,
    0.26794919243111803,
    0.26794919243112236,
    0.2679491924311227
];

export const L_omoms3 = [
    0.2490842490842491,
    0.33986928104575165,
    0.34362717574396406,
    0.34407148031876356,
    0.34412408743959544,
    0.34413031736062233,
    0.3441310551448089,
    0.34413114251779625,
    0.34413115286505125,
    0.3441311540904378,
    0.3441311542355558,
    0.34413115425274154,
    0.3441311542547768,
    0.34413115425501783,
    0.34413115425504637,
    0.34413115425504975,
    0.34413115425505014,
    0.3441311542550502
];

// Normalized to 5.25
export const OMOMS3 = (x) => {
    if (x <= 1.0)
        return 3.25 + x * (0.375 - x * (5.25 - 2.625 * x));
    return 7.25 - x * (10.625 - x * (5.25 - 0.875 * x));
};

// Copyright (c) 2024-2025 silverslither.

export const Triangle = (x) => {
    1.0 - x;
};

export const Hermite = (x) => {
    return 1.0 - x * x * (3.0 - 2.0 * x);
};

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

const hsin = (x) => {
    const x2 = x * x;
    let v = 1.6338022763241866 - 0.7211254631780478 * x2;
    v = 1.0 - v * x2;
    return v * x;
};

const hcos3w = (x) => {
    const x2 = x * x;
    let v = 0.023077941650754795 - 0.0007869230084260478 * x2;
    v = 0.25311490431737477 - v * x2;
    return 0.46164 - v * x2;
};
export const Hamming3 = (x) => {
    if (x < 7.450580596923828e-9)
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
    let v = 0.007302004975434134 - 0.00014005538895082734 * x2;
    v = 0.1423771336785233 - v * x2;
    return 0.46164 - v * x2;
};
export const Hamming4 = (x) => {
    if (x < 7.450580596923828e-9)
        return 1.0;
    const round_x = Math.round(x);
    const sign = (round_x << 31 >> 31) | 1;
    const poly_x = sign * (x - round_x);
    const flip = x > 2.0;
    const w = hcos4w(flip ? x - 4.0 : x);
    return hsin(poly_x) * (0.53836 + (flip ? -w : w)) / x;
};

export const L_bspline2i = new Float64Array([
    0.14583333333333334,
    0.17142857142857143,
    0.1715686274509804,
    0.17157275021026072,
    0.17157287157287157,
    0.1715728751454532,
    0.17157287525062018,
    0.171572875253716,
    0.17157287525380713,
    0.17157287525380982,
    0.1715728752538099
]);

export const L_bspline3i = new Float64Array([
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
]);

export const L_omoms3 = new Float64Array([
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
]);

// Normalized to 5.25
export const OMOMS3 = (x) => {
    if (x <= 1.0)
        return 3.25 + x * (0.375 - x * (5.25 - 2.625 * x));
    return 7.25 - x * (10.625 - x * (5.25 - 0.875 * x));
};

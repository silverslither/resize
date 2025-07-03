// Copyright (c) 2024-2025 silverslither.

#ifndef RESIZE_FILTERS
#define RESIZE_FILTERS

double Triangle(double x);
double Hermite(double x);

double BSpline2(double x);
double BSpline3(double x);

double MitNet(double x);
double CatRom(double x);

double Hamming3(double x);
double Hamming4(double x);
double Hamming8(double x);

double OMOMS3(double x);
double OMOMS7(double x);
double OMOMS11(double x);

extern double LU_bspline2i[1 * 11 + 1];
extern double LU_bspline3i[1 * 14 + 1];

extern double LU_omoms3[1 * 18 + 1];
extern double LU_omoms7[3 * 35 + 9];
extern double LU_omoms11[5 * 54 + 25];

#endif

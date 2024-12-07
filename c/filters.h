// Copyright (c) 2024 silverslither.

#ifndef RESIZE_FILTERS
#define RESIZE_FILTERS

double Triangle(double x);

double Hermite(double x);

double BSpline2(double x);

double BSpline3(double x);

double MitNet(double x);

double CatRom(double x);

double MKS2013(double x);

double Lanczos3(double x);

double Lanczos4(double x);

extern double L_bspline3i[15];
double pNormBSpline3(double x);

extern double L_omoms3[18];
double pNormOMOMS3(double x);

#endif

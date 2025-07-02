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
double OMOMS5(double x);
double OMOMS7(double x);
double OMOMS9(double x);

extern double LU_bspline2i[12];
extern double LU_bspline3i[15];

extern double LU_omoms3[19];
extern double LU_omoms5[62];
extern double LU_omoms7[126];
extern double LU_omoms9[212];

#endif

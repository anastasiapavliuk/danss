#ifndef MUPROCESSING_FIT_H_INCLUDED
#define MUPROCESSING_FIT_H_INCLUDED
#pragma once

#include "data.h"
#include "muprocessing.h"

int FitErr(int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *p0, double *p1, double *chi2, double *ndf);
int MyFit(int eventID, int nsipm, int side, std::vector<double>& xsipmarr, std::vector<double>& zxsipmarr, double* p0, double* p1, double* chi2, double* ndf);
int RotateAndFitErr (int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *old_p0, double *old_p1, double *chi2, double *ndf, int *vertical);
int RotateAndFit (int eventID, int nsize,  double *old_x, double *old_z, double *old_p0, double *old_p1);
double mnk_p0 (int n, double* x, double* y);
double mnk_p1(int n, double* x, double * y);

int FFitErr(int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *p0, double *p1, double *chi2, double *ndf, int *vertical);
int FFit(int eventID, int n, double * x, double * z, double *p0, double *p1);

#endif
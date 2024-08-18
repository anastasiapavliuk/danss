#pragma once

#include "data.h"

Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);
void langaus(TH1D*hSNR, int MCflag, double low, double up, double number, double *energy, double *energy_err, double *sigma, double *sigma_err, double *width, double *width_err, int k, double *chisqr, int *ndf, int units = 1);
Double_t langaufunfix(Double_t *x, Double_t *par);
TF1 *langaufitfix(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Double_t fix);
void langaus_energy(TH1D*hSNR, int MCflag, double low, double up, double number, double *energy, double *energy_err, double *sigma, double *sigma_err, double *width, double *width_err, int k, double *fwhm, double fix_value);

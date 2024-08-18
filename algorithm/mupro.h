#ifndef MUPRO_H
#define MUPRO_H
#pragma once

#include "data.h"
#include "TChain.h"
#include "TLeafD.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"

double PrimX(double x);
double PrimY(double y);
double PrimZ(double z);
int AlexSiPMNumbering(int SiPMID);
double SiPMZCoordinate(int SiPMID);
double SiPMXYCoordinate(int SiPMID);
double meanXY (int n, double* x, double * y);
double meanX (int n, double* x);
double mnk_p0 (int n, double* x, double* y);
double mnk_p1(int n, double* x, double * y);
int RotateAndFit (int eventID, int nsize,  double *old_x, double *old_z, double *old_p0, double *old_p1);
int FFit(int eventID, int n, double * x, double * z, double *p0, double *p1);
int DistanceX (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout);
int DistanceY (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout);
int Cross(double x1, double x2, double z1, double z2, double p0, double p1, double *xin, double *xout, double *zin, double *zout);
int LengthX (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout);
int LengthY (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout);
double Length(int SiPMID, double px0, double px1, double py0, double py1, double *xin, double *xout, double *yin, double *yout, double *zin, double *zout);

#endif
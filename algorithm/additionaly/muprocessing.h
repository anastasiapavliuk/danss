#ifndef MUPROCESSING_H
#define MUPROCESSING_H
#pragma once

#include "data.h"
#include "TChain.h"
#include "TLeafD.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"

void FillLeavesEvent(TChain* parttreechain, TLeafD** sigleafevid, TLeafD** sigleafnsipm, TLeafD** sigleafrealevid);
void FillLeavesSiPM(TChain* parttreechain, TLeafD** sipmleafevid, TLeafD** sipmleafsipmid, TLeafD** sipmleafsig);
int IrinaSiPMNumbering(int SiPMID);
int ReverseIrinaSiPMNumbering(int SiPMID);
int ReverseAlexSiPMNumbering(int SiPMID);
int AlexSiPMNumbering(int SiPMID);
double PrimX(double x);
double PrimY(double y);
double PrimZ(double z);
int SiPMIDbyCoordinates(double layer, double column);
double SiPMZCoordinate(int SiPMID);
double SiPMXYCoordinate(int SiPMID);
int CalculateEventIdByGlobalTime(long global_time);
int Cross(double x1, double x2, double z1, double z2, double p0, double p1, double *xin, double *xout, double *zin, double *zout);
int LengthX (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout);
int LengthY (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout);

double Length(int SiPMID, double px0, double px1, double py0, double py1, double *xin, double *xout, double *yin, double *yout, double *zin, double *zout);

int DistanceX (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout);
int DistanceY (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout);

bool AreNeighbors(TStrip a, TStrip b);
int StripHasNeighbors(TStrip curstrip, const std::vector<TStrip>& strip);
int RecNeighborFunc(int status, TStrip strip, std::vector<TStrip>& stripvec, int n);
int FindMaxConnectedSpace (std::vector<TStrip>& stripvec, double p0, double p1, TH1D *hist = 0, bool TestMode = false);
void HoughTransform(TH2D* hist, const std::vector<TStrip>& strip, double* p0, double* p1);
void RotateAndHoughTransform (TH2D* hist, const std::vector<TStrip>& strip, double* old_p0, double* old_p1);

double meanXY (int n, double* x, double * y);
double meanX (int n, double* x);
double mean_distance(double* x, double* z, double p0, double p1);

double PointLineDistance (double x, double z, double p0, double p1);
double CutDelta(int eventID, const std::vector<TStrip>& strip, int flag);
double maxX (double *x);
double minX (double *x);
double LengthToMev(double l);
double LengthErrorToMev(double l);

double PheToMeV(double phe);
bool StriIsBad (int SiPMID);
//void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//Double_t LineFunc(double x, Double_t *par);
void FillEnergyLengthHists(std::vector <TH1D>& energy_hist_vec, std::vector <TH1D*> length_vec, std::vector <TH1D*> left_border_l, std::vector <TH1D*> right_border_l, double energy_sipm, double length_sipm, int SiPMID);
int FindTracks(std::vector<TStrip>& x_strip_vector, std::vector<TStrip>& y_strip_vector, bool TestMode, std::vector<int> interesting_events, int current_event_ID, int n_sipm, double phe_threshold, bool MCflag, bool DrawMode, TrackStruct &track, TFitResultPtr &xfit, TFitResultPtr &yfit);

// bool StripCondition (double up, double low, int side, double left, double right, int SiPMID);
bool StripCondition (int SiPMID);
void FillEnergyLengthHists(std::vector <TH1D>& energy_hist_vec, std::vector <TH1D*> length_vec, std::vector <double> left_border_l, std::vector <double> right_border_l, double energy_sipm, double length_sipm, int SiPMID);

double AttenuationCorrectionSiPM(double distance, bool mode);

#endif
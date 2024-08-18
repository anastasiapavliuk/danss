#ifndef MUPROCESSING_NAMES_H
#define MUPROCESSING_NAMES_H

#pragma once

#include "data.h"
// #include "muprocessing.h"

std::string DoubleToEventName (double number, int MCflag, int TestMode);
std::string DoubleToSiPMName (double number, int MCflag, int TestMode);
std::string InputFileName (int flag, std::string directory = std::string());
std::string AllSiPMFileName(int MCflag, int TestMode, std::string directory = std::string());
std::string ProcessingResultsFileName(int MCflag, int TestMode, std::string directory = std::string());
std::string DoubleToPNGName (double number, int side, int flag, const std::string& description = std::string());
std::string DoubleToEPSName (double number, int side, int flag, const std::string& description = std::string());
std::string DoubleToTitleName (double number, int side, double chi2=-10, double ndf=-10);
std::string DoubleToTitleNameEnergyHist (double low_e, double up_e, double low_d, double up_d);
std::string TwoDoubleToString(double num1, double num2, const std::string& description, const std::string& suffix = std::string());

#endif
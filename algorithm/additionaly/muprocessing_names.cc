#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

#include "TTree.h"
#include "TBranch.h"
#include "TLeafD.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "Fit/Fitter.h"
#include "TFitter.h"
#include "TStyle.h"
#include "TMarker.h"

#include "data.h"
#include "muprocessing_names.h"
// #include "muprocessing.h"
// #include "muprocessing_fit.h"

using std::cout;
using std::endl;
using std::cin;

std::string DoubleToEventName (double number, int MCflag, int TestMode)
{
    
    //out << std::fixed << std::precision(9) << s << "-" << m << ".root";
    std::ostringstream out;
    if(TestMode){
        out << "../output_test/eventID" <<number << ".root"; 
    }
    else{
        if (MCflag == 0)out << "../output_data/eventID" <<number << ".root"; 
        else out << "../output_mc/eventID" <<number << ".root"; 
    }
    return out.str();
}

std::string DoubleToSiPMName (double number, int MCflag, int TestMode)
{
    std::ostringstream out;
    if(TestMode){
        out << "../output_test/SiPMID" <<number << ".root";
    }
    else{
        if (MCflag == 0)out << "../output_data/SiPMID" <<number << ".root"; 
        else out << "../output_mc/SiPMID" <<number << ".root"; 
    }
    return out.str();
}

std::string InputFileName (int flag, std::string directory)
{
    std::ostringstream out;
    
    if (flag == 0)out << "../input_data/"<<directory<<"*.root"; 
    else out << "../input_mc/"<<directory<<"*.root"; 
    return out.str();
}

std::string AllSiPMFileName(int MCflag, int TestMode, std::string directory)
{
    std::ostringstream out;
    if(TestMode){
        out << "../output_test/AllSiPM.root";
    }
    else{
        if (MCflag == 0)out << "../output_data/"<<directory<<"AllSiPM_0_19.root";
        else out << "../output_mc/"<<directory<<"AllSiPM.root";
    }
    return out.str();
}

std::string ProcessingResultsFileName(int MCflag, int TestMode, std::string directory)
{
    std::ostringstream out;
    if(TestMode){
        out << "../output_test/ProcessingResults.root";
    }
    else{
        if (MCflag == 0)out << "../output_data/"<<directory<<"ProcessingResults_2.root"; 
        else out << "../output_mc/"<<directory<<"ProcessingResults.root"; 
    }
    return out.str();
}

std::string DoubleToPNGName(double number, int side, int flag, const std::string& description ) // 0 -- X, 1 -- Y in my coordinates
{
    std::ostringstream out;
    if (flag == 0)
    {
        if (side == 0)out << "../output_data/" <<description<<number << "x.png";
        else out << "../output_data/" <<description<<number << "y.png";
    }
    else
    {
        if (side == 0)out << "../output_mc/" <<description<<number << "x.png";
        else out << "../output_mc/" <<description<<number << "y.png";
    }
    return out.str();
}

std::string DoubleToEPSName(double number, int side, int flag, const std::string& description ) // 0 -- X, 1 -- Y in my coordinates
{
    std::ostringstream out;
    if (flag == 0)
    {
        if (side == 0)out << "../output_data/" <<description<<number << "x.pdf";
        else out << "../output_data/" <<description<<number << "y.pdf";
    }
    else
    {
        if (side == 0)out << "../output_mc/" <<description<<number << "x.pdf";
        else out << "../output_mc/" <<description<<number << "y.pdf";
    }
    return out.str();
}

std::string DoubleToTitleName (double number, int side, double chi2, double ndf)
{
    std::ostringstream out;
    if (ndf<0)
    {
        if (side == 0)out << "X side, eventID: " <<number; 
        else out << "Y side, eventID: " <<number;
        
    }
    else
    {
        if (side == 0) out << "X side, eventID: "<<number<<", chi: "<<chi2<<", Ndf: "<<ndf; 
        else out << "Y side, eventID: " <<number<<", chi: "<<chi2<<", Ndf: "<<ndf;
    }
    return out.str();
}

std::string TwoDoubleToString(double num1, double num2, const std::string& description, const std::string& suffix)
{
    std::ostringstream out;
    out << description <<num1<<"_"<<num2<<suffix; 
        
    return out.str();
}

std::string DoubleToTitleNameEnergyHist (double low_e, double up_e, double low_d, double up_d)
{
    std::ostringstream out;
    
    if (low_d>=0)out<<"Length: ("<<low_e<<" -- "<<up_e<<")cm, Distance: ("<<low_d<<"-"<<up_d<<") cm";
    else out<<"Length: ("<<low_e<<" -- "<<up_e<<") cm";
    return out.str();
}


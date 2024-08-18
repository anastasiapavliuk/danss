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

#include "data.h"
// #include "muprocessing.h"
// #include "muprocessing_names.h"
#include "danss_constants.h"
#include "muprocessing_fit.h"

using std::cout;
using std::endl;
using std::cin;

int MyFit(int eventID, int nsipm, int side, std::vector<double>& xsipmarr, std::vector<double>& zxsipmarr, double* p0, double* p1, double* chi2, double* ndf)
{
    
    std::vector<double> xerr;
    std::vector<double> zerr;
    double pp0 = *p0;
    double pp1 = *p1;
    for (int i = 0; i <nsipm; i++)
    {
        xerr.push_back(0.5);
        zerr.push_back(2.);
         //xerr[i] = 0.5;
         //zerr[i] = 2.;
    }
    TFitResultPtr ffit;
    double p0limit = 2.;//TODO
    
    if (nsipm < xsipmarr.size())std::cerr << "Array size is less the SiPM number! Event: " << eventID << "\n";
    TGraphErrors* gx = new TGraphErrors(nsipm, xsipmarr.data(), zxsipmarr.data(), xerr.data(), zerr.data()); 
    TF1* fx = new TF1("fx","[0]*x+[1]", 0., 104.); 
    fx->SetParameter(0, pp0);      
    fx->SetParameter(1, pp1);  
    fx->SetParLimits(0, pp0-p0limit,pp0+p0limit);
    TVirtualFitter::Fitter(gx)->SetPrecision(FitPrecision);
    //TVirtualFitter::Fitter(gx)->SetMaxIterations(5000);
    
    
  
    ffit = gx->Fit("fx","SRQ");
    //if (eventID == 35) ffit = gx->Fit("fx","SRV");
    //else 
    *p0 = fx->GetParameter(0);
    *p1 = fx->GetParameter(1);
    
    /*TCanvas* c1 = new TCanvas("c1", "Signal data", 30, 20, 1000, 1000);
    c1->cd();
    if (eventID == 121 && side == 0){
        gx->SetMarkerStyle(20);
        gx->SetMarkerSize(1);
        gx->SetFillColor(kCyan-3);
        gx->SetFillStyle(3014);
        gx->SetLineColor(kCyan+3);
        gx->SetTitle(" ");
        gx->Draw("A5Z");
        gx->GetXaxis()->SetLimits(0, 100);
        gx->SetMinimum(-100);
        gx->SetMaximum(0);
        fx->Draw("Same");
        c1->SaveAs(DoubleToPNGName(eventID, side, 1, "fromfit").c_str());
    }
    delete c1;*/
    //*p0 = mnk_p0(nsipm, xsipmarr, zxsipmarr);
   // *p1 = mnk_p1(nsipm, xsipmarr, zxsipmarr);
    if(int(ffit) == 0){
        *ndf = ffit->Ndf();
        *chi2 = ffit->Chi2();
    }
    delete gx;
    delete fx;
    return int(ffit);
}

int FitErr(int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *p0, double *p1, double *chi2, double *ndf)
{
    std::vector<double> xsipmarr;
    std::vector<double> zxsipmarr;
    int k = 0;
    for (int i = 0; i < strip.size(); i++)
    {
        if (strip[i].status!=0)   
        {
            xsipmarr.push_back(strip[i].XY);
            zxsipmarr.push_back(strip[i].Z);
            k++;
        }
        /*else {xsipmarr[i] = -4000; zxsipmarr[i] = -4000;}*/
        
    }
    if (k!=nsipm)std::cerr<<"k-nsipm: "<<k-nsipm<<"\t"<<eventID<<"\n";
    int ret = MyFit(eventID, nsipm, side, xsipmarr, zxsipmarr, p0, p1, chi2, ndf);
    
    return ret;
}

int RotateAndFitErr (int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *old_p0, double *old_p1, double *chi2, double *ndf, int *vertical)
{
    
    double epsilon = LengthXY/DANSSConstants::DetectorSizeVertical;
    std::vector<double> x;
    std::vector<double> z;
   
    double p0, p1;
   
    //Transform p0 and p1
    p0 = -1./(*old_p0);
    p1 = (*old_p1)/(*old_p0);
   
    int k = 0;
    for (int i = 0; i < strip.size(); i++)
    {
        
        if (strip[i].status!=0)
        {
            x.push_back(strip[i].Z);
            z.push_back(-(strip[i].XY));
            k++;
        }
        /*else {x[i] = -4000; z[i] = -4000;}*/
    }
    
    if (k!=nsipm)std::cerr<<"k-nsipm: "<<k-nsipm<<"\t"<<eventID<<"\n";
    int ret = MyFit(eventID, nsipm, side, x, z, &p0, &p1, chi2, ndf);
    
    //Transform p0 and p1 back
    if (fabs(p0)<epsilon){
        *vertical = 1;
        //cout<<"Vertical track in "<<eventID<< " event\n";
    }
    else {
        *vertical = 0;
    }
    *old_p0 = -1/p0;
    *old_p1 = -p1/p0;
    //cout<<"Converted p0:\t"<<*old_p0<<"\tp1\t"<<*old_p1<<"\n";
   
    return ret;
}


int FFitErr(int eventID, int nsipm, int side, const std::vector<TStrip>& strip, double *p0, double *p1, double *chi2, double *ndf, int *vertical)
{
    if (fabs(*p0) < 1) {
        *vertical = 0;
        return  FitErr(eventID, nsipm, side, strip, p0, p1, chi2, ndf);
    }
    else if (fabs(*p0)>=1)return RotateAndFitErr (eventID, nsipm, side, strip, p0, p1, chi2, ndf, vertical);
    else return -10;
}


double mnk_p0 (int n, double * x, double * y)
{
     return (meanXY(n, x, y)- meanX(n, x)*meanX(n, y))/(meanXY (n, x,x)-meanX(n,x)*meanX(n,x));
}

double mnk_p1(int n, double * x, double * y)
{
    return meanX(n,y) - mnk_p0(n,x,y)*meanX(n,x);
}

int RotateAndFit (int eventID, int nsize,  double *old_x, double *old_z, double *old_p0, double *old_p1)
{
    
    double epsilon = 1E-40;
    double *x = new double[nsize];
    double *z = new double[nsize];
    
    double p0limit = 1.;//TODO
    
    int ret = -500;
    double p0, p1;
    int i = 0;
    int k = 0;
    //Transform p0 and p1
    p0 = -1/(*old_p0);
    p1 = (*old_p1)/(*old_p0);
    
    
    //TFitResultPtr ffit;
    for (i = 0; i < nsize; i++)
    {
        if ((old_x[i] != -10) && (old_z[i] != -10))
        {
            x[k] = (old_z[i]);
            z[k] = (-old_x[i]);
           
            k++;
        } 
    }
    
    
    p0 = mnk_p0(nsize, x, z);
    p1 = mnk_p1(nsize, x, z);

    
    ret = 0;   
    //ret = int(ffit);
    
    //Transform p0 and p1 back
    //*ndf = ffit->Ndf();
    //*chi2 = ffit->Chi2();
    
    *old_p0 = -1/p0;
    *old_p1 = -p1/p0;
    //cout<<"Converted p0:\t"<<*old_p0<<"\tp1\t"<<*old_p1<<"\n";
    delete []x;
    delete []z;
    return ret;
}



int FFit(int eventID, int n, double * x, double * z, double *p0, double *p1)
{
    if (fabs(*p0) < 1) {
        
        *p0 = mnk_p0(n, x, z);
        *p1 = mnk_p1(n, x, z);
        return 0;
    }
    else if (fabs(*p0)>=1)return RotateAndFit (eventID, n, x, z, p0, p1);
    else return -10;
}

#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>


#include "TTree.h"
#include "TChain.h"
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
#include "muprocessing_fit.h"
#include "muprocessing_names.h"
#include "muprocessing.h"

using std::cout;
using std::endl;
using std::cin;

void FillLeavesEvent(TChain* parttreechain, TLeafD** sigleafevid, TLeafD** sigleafnsipm, TLeafD** sigleafrealevid)
{
    *sigleafevid = (TLeafD*) parttreechain->GetLeaf("DANSSSigEventBranch", "EventID");
    *sigleafnsipm = (TLeafD*) parttreechain->GetLeaf("DANSSSigEventBranch", "NSiPM");
    *sigleafrealevid = (TLeafD*) parttreechain->GetLeaf("DANSSSigEventBranch", "RealEventID");
}

void FillLeavesSiPM(TChain* parttreechain, TLeafD** sipmleafevid, TLeafD** sipmleafsipmid, TLeafD** sipmleafsig)
{
    *sipmleafevid = (TLeafD*) parttreechain->GetLeaf("DANSSSigSiPMBranch", "EventID");
    *sipmleafsipmid = (TLeafD*) parttreechain->GetLeaf("DANSSSigSiPMBranch", "SiPMID");
    *sipmleafsig = (TLeafD*) parttreechain->GetLeaf("DANSSSigSiPMBranch", "Signal");
}
double PrimX(double y)
{
    return y/1e4+50;
}
double PrimY(double x)
{
    return 100-(x/1e4+50);
}
double PrimZ(double z)
{
    return z/1e4+52;
}
int CalculateEventIdByGlobalTime(long global_time) {
    int cycles_per_second = 125000000;
    return 1 + std::round(static_cast<long double>(global_time) / cycles_per_second);
}
int IrinaSiPMNumbering(int SiPMID)//from Irina's to mine
{
    int row = -10;
    int column = -10;
    //cout<<"Current SiPMID: "<<SiPMID<<endl;
    if (SiPMID < 1250)
    {
        row = (SiPMID/25)*2;
        column = 25 - SiPMID%25;
        SiPMID = 25*row+column-1;
    }
    else
    {
        SiPMID-=1250;
        row = (SiPMID/25)*2+1;
        column = SiPMID%25;
        SiPMID = 25*row+column; 
    }//cout<< "row: "<<row<<" column: "<<column<<" new ID: "<<SiPMID<<endl;
    return SiPMID;
}

int ReverseIrinaSiPMNumbering(int SiPMID)//to Irina's
{
    int row = -10;
    int column = -10;
    //cout<<"Current SiPMID: "<<SiPMID<<endl;
    if (((SiPMID)/25)%2==0)// ||| xside
    {
        row = SiPMID/25;
        column = SiPMID%25;
        SiPMID = row*25/2+24-column;   
    }
    else
    {
        row = SiPMID/25;
        column = SiPMID%25;
        SiPMID = (row-1)*25/2+column+1250;
    }//cout<< "row: "<<row<<" column: "<<column<<" new ID: "<<SiPMID<<endl;
    return SiPMID;
}
int AlexSiPMNumbering(int SiPMID){//TODO
    int row = -10;
    int column = -10;
    if (((SiPMID)/25)%2==0){ //my x, Sasha y
            
            column = SiPMID%25; //as is
            row = SiPMID/25;
            
        }
        else{
            column = 24-SiPMID%25; // flip
            row = SiPMID/25;
            
        }
    SiPMID = 25*row+column; 
    return SiPMID;    
    
    
}
int ReverseAlexSiPMNumbering(int SiPMID)//to Alex
{
    int row = -10;
    int column = -10;
    if (((SiPMID)/25)%2==0){ //my x, Sasha y
            
            column = SiPMID%25; //as is
            row = SiPMID/25;
            
        }
        else{
            column = 24-SiPMID%25; // flip
            row = SiPMID/25;
            
        }
    SiPMID = 25*row+column; 
    return SiPMID;
}

double SiPMZCoordinate(int SiPMID){
    return double(SiPMID/25)*(LengthZ+GapLength)+HalfLengthZ;   
}

double SiPMXYCoordinate(int SiPMID){
    return double((SiPMID%25)*LengthXY+HalfLengthXY);   
}

int SiPMIDbyCoordinates(double layer, double column)//good only for (0.5, 2), (1.5, 98) and like this
{
    
    return 25*int(layer)+int(column/4);
}

double meanXY (int n, double * x, double * y)
{
    //int n = 2500;
    double RealLength = 0;
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        if ((x[i] != -10) && (y[i] != -10))
        {
            sum+=(x[i]*y[i]);
            RealLength++;
        }

    }
    return sum/RealLength;
}

double meanX(int n, double *x)
{
    //int n = 2500;
    double RealLength = 0;
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        if (x[i] != -10) 
        {
            sum+=(x[i]);
            RealLength++;
        }

    }
    return sum/RealLength;
}

double mean_distance(double *x, double *z, double p0, double p1)
{
    int n = 2500;
    double RealLength = 0;
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        if (z[i] != -10) 
        {
            sum+=fabs(double(z[i] - (p0*x[i]+p1)));
            RealLength++;
        }

    }
    return sum/RealLength;
}


int Cross(double x1, double x2, double z1, double z2, double p0, double p1, double *xin, double *xout, double *zin, double *zout)//fucking conditions!
{
    
    
    //cout <<"****************\n";
    //cout << "p0: "<<p0<<"p1: "<<p1<<endl;
    //cout << "x1: " << x1<<"\nx2: "<<x2<<"\nz1: "<<z1<<"\nz2: "<<z2<<endl;
    double xcross1 = (z1-p1)/p0;
    double xcross2 = (z2-p1)/p0;
    double zcross1 = p0*x1+p1;
    double zcross2 = p0*x2+p1;
    
    //if (xcross1 > 100 || xcross1 <0 || xcross2 > 100 || xcross2 <0) return -1;
    //cout << "xc1: " << xcross1<<"\nxc2: "<<xcross2<<"\nzc1: "<<zcross1<<"\nzc2: "<<zcross2<<endl;


//z2      xcross2
//    _____________
//   |             |
//   |             | zcrsoss2
//   |_____________|
//z1
//
//   x1   xcross1  x2

    if ((zcross1>z1) && (zcross1 < z2))//left hit
    {
        if ((xcross2>x1) && (xcross2<x2))
        {
            *xin = x1;
            *xout = xcross2;
            *zin = zcross1;
            *zout = z2;
            return 0;
            //  ____
            // |
        }
        else if ((zcross2<z2)&&(zcross2>z1))
        {
            *xin = x1;
            *xout = x2;
            *zin = zcross1;
            *zout = zcross2;
            return 0;
            // |    |
        }
        else if ((xcross1>x1)&&(xcross1<x2))
        {
            *xin = x1;
            *xout = xcross1;
            *zin = zcross1;
            *zout = z1;
            return 0;
            // |____
        }
    }
    else if ((xcross2>x1) && (xcross2<x2))//up hit
    {
        if ((zcross2<z2)&&(zcross2>z1))
        {
            *xin = xcross2;
            *xout = x2;
            *zin = z2;
            *zout = zcross2;
            return 0;
            //  ____
            //      |
        }
        else if ((xcross1>x1)&&(xcross1<x2))
        {
            *xin = xcross1;
            *xout = xcross2;
            *zin = z1;
            *zout = z2;
            return 0;
            //  ____
            //  ____
        }
    }
    else if ((zcross2<z2)&&(zcross2>z1))//right hit
    {
        if ((xcross1>x1)&&(xcross1<x2))
        {
            *xin = xcross1;
            *xout = x2;
            *zin = z1;
            *zout = zcross2;
            return 0;
            //  ____|
        }
    }
    if ((zcross1 == z2) && (xcross2 == x1))//up left corner
    {
        cout << "Up left\n";
        if (p0 > 0) return -10;
        if ((zcross2 == z1) && (xcross1 == x2))//down right corner
        {
            cout<<"And down right too\n";
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross1;
            *zout = zcross2;
            return 0;
        }
        else if ((zcross2<z2)&&(zcross2>z1))//right hit
        {
            *xin = xcross2;
            *xout = x2;
            *zin = zcross1;
            *zout = zcross2;
        }
        else if ((xcross1>x1)&&(xcross1<x2))//down hit
        {
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross1;
            *zout = z1;
            return 0;
        }   
        else return -10;
    }
    else if ((zcross2 == z2) && (xcross2 == x2))//up right corner
    {
        cout << "Up right\n";
        if (p0 < 0) return -10;
        if ((zcross1 == z1) && (xcross1 == x1))//down left corner
        {
            cout << "And down left too\n";
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((zcross1>z1) && (zcross1 < z2))//left hit
        {
            *xin = xcross2;
            *xout = x1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((xcross1>x1)&&(xcross1<x2))//down hit
        {
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross2;
            *zout = z1;
            return 0;
        }
        else return -10;
    }
    else if ((zcross2 == z1) && (xcross1 == x2))//down right
    {
        cout << "Down right\n";
        if (p0 > 0) return -10;
        else if ((zcross1>z1) && (zcross1 < z2))//left hit
        {
            *xin = xcross1;
            *xout = x1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((xcross2>x1) && (xcross2<x2))//up hit
        {
            *xin = xcross1;
            *xout = xcross2;
            *zin = zcross2;
            *zout = z2;
            return 0;
        }
        else return -10;
        
    }
    
    else if ((zcross1 == z1) && (xcross1 == x1))//down left
    {
        cout << "Down left\n";
        if (p0 < 0) return -10;
        else if ((xcross2>x1) && (xcross2<x2))//up hit
        {
            *xin = xcross1;
            *xout = xcross2;
            *zin = zcross1;
            *zout = z2;
            return 0;
        }
        else if ((zcross2<z2)&&(zcross2>z1))//right hit
        {
            *xin = xcross1;
            *xout = x2;
            *zin = zcross1;
            *zout = zcross2;
            return 0;
        }
        else return -10;
    }
    
    return -10;
}

int LengthX (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout)
{
    if (((SiPMID)/25)%2!=0)return -100;
    double length = -10;
    double x1 = SiPMXYCoordinate(SiPMID)-HalfLengthXY+NoSenseLength;    
    double x2 = SiPMXYCoordinate(SiPMID)+HalfLengthXY-NoSenseLength;
    double z1 = SiPMZCoordinate(SiPMID)-HalfLengthZ+NoSenseLength;
    double z2 = SiPMZCoordinate(SiPMID)+HalfLengthZ-NoSenseLength;
    
    //z2      xcross2
//    _____________
//   |             |
//   |             | zcrsoss2
//   |_____________|
//z1
//
//   x1   xcross1  x2
    
    return Cross (x1, x2, z1, z2, px0, px1, xin, xout, zin, zout);
}
int LengthY (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout)
{
    if (((SiPMID)/25)%2==0)return -100;//Wrong side!
    double length = -10;
    double y1 = SiPMXYCoordinate(SiPMID)-HalfLengthXY+NoSenseLength;    
    double y2 = SiPMXYCoordinate(SiPMID)+HalfLengthXY-NoSenseLength; 
    double z1 = SiPMZCoordinate(SiPMID)-HalfLengthZ+NoSenseLength;
    double z2 = SiPMZCoordinate(SiPMID)+HalfLengthZ-NoSenseLength;
    
//z2      ycross2
//    _____________
//   |             |
//   |             | zcrsoss2
//   |_____________|
//z1
//
//   y1   ycross1  y2
    return Cross (y1, y2, z1, z2, py0, py1, yin, yout, zin, zout);
}
int DistanceX (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout)
{
    //double distance = -10;
 
    //double z1 = double((SiPMID-1)/25);
    //p0 = p0;
    *yin =  (*zin-py1)/py0;
    *yout = (*zout-py1)/py0;
    if ((*yin < 0) || (*yin >100) || (*yout < 0) || (*yout >100)) return -10;
    else return 0;//Status OK
}

int DistanceY (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout)
{
    *xin =  (*zin-px1)/px0;
    *xout = (*zout-px1)/px0;
    if ((*xin < 0) || (*xin >100) || (*xout < 0) || (*xout >100)) return -10;
    else return 0;//Status OK
}


double Length(int SiPMID, double px0, double px1, double py0, double py1, double *xin, double *xout, double *yin, double *yout, double *zin, double *zout)
{
    double cross = -10;
    double d_status = -10;
    if (((SiPMID)/25)%2==0){ //x side
        cross = LengthX(SiPMID, px0, px1, xin, xout, zin, zout);
        d_status = DistanceX(SiPMID, py0, py1, yin, yout, zin, zout);
        
    }
    else{
        cross = LengthY(SiPMID, py0, py1, yin, yout, zin, zout);
        d_status = DistanceY(SiPMID, px0, px1, xin, xout, zin, zout);
    }
    if (cross == 0 && d_status == 0){
            return sqrt ((*xout-*xin)*(*xout-*xin) + (*yout-*yin)*(*yout-*yin) + (*zout-*zin)*(*zout-*zin));
    }
    return 0;
    
}
bool AreNeighbors(TStrip a, TStrip b)
{
    if (fabs(a.XY-b.XY)<=LengthXY && a.status!=0 && b.status!=0 && fabs(a.Z-b.Z) <= 2*(LengthZ+GapLength)+GapEpsilon) return true;
    
    else return false;
}
int StripHasNeighbors(TStrip curstrip, const std::vector<TStrip>& strip)
{
    int n = 0;
    for (int i = 0; i < strip.size(); i++)
    {
        if ((curstrip.SiPMID!=strip[i].SiPMID) && (AreNeighbors(curstrip, strip[i]))) n++;  
        
    }
    return n;
}

int RecNeighborFunc(int status, TStrip strip, std::vector<TStrip>& stripvec, int n)
{
    
    for (int i = 0; i < stripvec.size(); i++)
    {
       
        if (AreNeighbors(strip, stripvec[i]) && stripvec[i].status != status)
        {
            stripvec[i].status = status;
            n++;
            n = RecNeighborFunc(status, stripvec[i], stripvec, n);
        }
    }
    return n;
}

int FindMaxConnectedSpace (std::vector<TStrip>& stripvec, double p0, double p1, TH1D *hist, bool TestMode)
{
    int nsipm = 0;
    int k = 1;
    int kmax = 1;
    int nmax = 0;    
    for (int i = 0; i<stripvec.size(); i++){
        if (TestMode){
                cout<<"SiPMID: "<<stripvec[i].SiPMID<<", XY: "<< stripvec[i].XY<<", Z: " <<  stripvec[i].Z<<", distance: "<< PointLineDistance(stripvec[i].XY,stripvec[i].Z, p0, p1)<<endl;
        }
        if(PointLineDistance(stripvec[i].XY,stripvec[i].Z, p0, p1)<2.2 && stripvec[i].status == 1){ //sqrt(2^2+0.5^2)==2.0616
            k++;
            if (TestMode){
                cout<<"SiPMID: "<<stripvec[i].SiPMID<<", XY: "<< stripvec[i].XY<<", Z: " <<  stripvec[i].Z<<endl;
                cout<<"k: "<<k<<endl;
            }
            nsipm = RecNeighborFunc(k, stripvec[i], stripvec, 0);
            if (TestMode) cout<<"nsipm: "<<nsipm<<endl;
            if (nsipm>nmax){
                nmax = nsipm;
                kmax = k;
            }
            
               
        }
    }
    if (hist!=0)hist->Fill(nmax);
    nsipm = stripvec.size();
        
        
    for (int i = 0; i < stripvec.size(); i++){
        if (stripvec[i].status!=kmax)
        {
           stripvec[i].status = 0;
           nsipm--;
        }
    }
    
    return nsipm;
}
void HoughTransform (TH2D* hist, const std::vector<TStrip>& strip, double* p0, double* p1)
{
    int n = strip.size();
    double max = 0;
    int xbin, zbin;
    int houghflag = 0;
    double angle;
    double R;
    for (int i = 0; i < n; i++)
    {
        for (int t = 0; t < M_PI*NHoughBinsA; t++)hist->Fill (double(t)/NHoughBinsA, strip[i].XY * cos(double(t)/NHoughBinsA)+strip[i].Z * sin(double(t)/NHoughBinsA));
                
    }
    
   
    for (int i = NHoughBinsA; i>0; i--)
    {
        for (int k = NHoughBinsR; k > 0; k--)
        {
            if (hist->GetBinContent(i,k)>=max)
            {
                max  = hist->GetBinContent(i,k);
                xbin = i;
                zbin = k;
                //houghflag++;
            }
        }
    }
        
        //if (houghflag>1)continue;
        
    angle = M_PI*double(xbin)/NHoughBinsA;
    R = 2.*HoughR*double(zbin)/(NHoughBinsR)-HoughR;
    *p0 = - cos(angle)/sin(angle);
    *p1 = R/sin(angle);   
}

void RotateAndHoughTransform (TH2D* hist, const std::vector<TStrip>& strip, double* old_p0, double* old_p1)
{
    int n = strip.size();
    double max = 0;
    int xbin, zbin;
    int houghflag = 0;
    double angle, R;

    double p0, p1;
    
    for (int i = 0; i < n; i++)
    {
        for (int t = 0; t < M_PI*NHoughBinsA; t++)hist->Fill (double(t)/NHoughBinsA, strip[i].Z * cos(double(t)/NHoughBinsA)-strip[i].XY * sin(double(t)/NHoughBinsA));
                
    }
    
   
    for (int i = NHoughBinsA; i>0; i--)
    {
        for (int k = NHoughBinsR; k > 0; k--)
        {
            if (hist->GetBinContent(i,k)>=max)
            {
                max  = hist->GetBinContent(i,k);
                xbin = i;
                zbin = k;
                //houghflag++;
            }
        }
    }
        
        //if (houghflag>1)continue;
        
    angle = M_PI*double(xbin)/NHoughBinsA;
    R = 2.*HoughR*double(zbin)/(NHoughBinsR)-HoughR;
    p0 = - cos(angle)/sin(angle);
    p1 = R/sin(angle); 
    *old_p0 = -1/p0;
    *old_p1 = -p1/p0;
}




double maxX (double *x)
{
    int n = 2500;
    double max = -20;
    for (int i = 0; i < n; i++)
    {
        if (x[i] > max) max = x[i] ;
        
    }
    return max;
    
}

double minX (double *x)
{
    int n = 2500;
    double min = 110;
    for (int i = 0; i < n; i++)
    {
        if ((x[i] < min) && (x[i]>=0)) min = x[i] ;
        
    }
    return min;
    
}

double CutDelta(int eventID, const std::vector<TStrip>& strip, int flag)
{
    double max = 0;
    double min = 100;
    if (flag == 0){
        for (int i = 0; i<strip.size(); i++){
           if ((strip[i].status !=0)&&(strip[i].XY > max)) max = strip[i].XY;
           if ((strip[i].status !=0)&&(strip[i].XY < min)) min = strip[i].XY;
        }
    }
    else {
        for (int i = 0; i<strip.size(); i++){
           if ((strip[i].status !=0)&&(strip[i].Z > max)) max = strip[i].Z;
           if ((strip[i].status !=0)&&(strip[i].Z < min)) min = strip[i].Z;
        }
    }
    //if (eventID == 4) cout<<"max:\t"<<max<<"\tminn:\t"<<min<<"\t"<<flag<<"\n";
    return max-min;
}
double PointLineDistance (double x, double z, double p0, double p1)
{
    return fabs(p0*x-z+p1)/sqrt(p0*p0+1);    
}

double LengthToMev(double l)
{
    //return l*1.777-0.1529;
//     return l*1.786-0.09708; //new MC //2019
    return l*1.89 - 0.1959;

}

double LengthErrorToMev(double l)
{
    //return (0.0028/1.777+0.001)*l+0.0035;
//     return (0.0054/1.786+0.001)*l+0.006; //new MC //2019
    return (0.005/1.9 + 0.001)*l+0.001;//FIXME

}

double PheToMeV(double phe)
{
    double a = 17.3571;
    double b = -2.96342;
    return (phe-b)/a;
}
bool StriIsBad (int SiPMID)
{   
    int BadStrips[] = {36, 150, 184, 190, 191, 192, 193, 194, 211, 218, 240, 241, 242, 243, 244, 272, 290, 291, 292, 293, 294, 310, 323, 361, 379, 389, 400, 464, 674, 695, 745, 769, 961, 1128, 1204, 1217, 1340, 1365, 1687, 1754, 1775, 1902, 1914, 1937, 1991, 2007, 2040, 2131, 2240, 2497};
    
    int *p;
    p = std::find(BadStrips, BadStrips+50, SiPMID);
    if (p != BadStrips+50) return true;   
    else return false;
}

int FindTracks(std::vector<TStrip>& x_strip_vector, std::vector<TStrip>& y_strip_vector, bool TestMode, std::vector<int> interesting_events, int current_event_ID, int n_sipm, double phe_threshold, bool MCflag, bool DrawMode, TrackStruct &track, TFitResultPtr &xfit, TFitResultPtr &yfit)
{
  
    int n_x_SiPM = x_strip_vector.size();
    int n_y_SiPM = y_strip_vector.size();
    
    Double_t deltaX, deltaXZ, deltaY, deltaYZ; 
    Double_t px0, px1, py0, py1;
    Double_t distance_cut = 2.5;
     
    Deleter<TF1> fx = new TF1("fx","[0]*x+[1]", 0., 100); 
    Deleter<TF1> fy = new TF1("fy","[0]*x+[1]", 0., 100); 
    
    Deleter<TH1D> horizontal_distribution_coord_x = new TH1D("h_d_c_x", "Strip coordinates X", 25, 0, 100);
    Deleter<TH1D> horizontal_distribution_coord_y = new TH1D("h_d_c_y", "Strip coordinates X", 25, 0, 100);
    
    for (int i = 0; i < x_strip_vector.size(); i++){
        horizontal_distribution_coord_x->Fill(x_strip_vector[i].XY);
    }
    
    for (int i = 0; i < y_strip_vector.size(); i++){
        horizontal_distribution_coord_y->Fill(y_strip_vector[i].XY);
    }
    
    if (horizontal_distribution_coord_x->GetMaximum() > 30) cout << "Likely vertical x, event: " << current_event_ID << "\n";
    if (horizontal_distribution_coord_y->GetMaximum() > 30) cout << "Likely vertical y, event: " << current_event_ID << "\n";
    
    TH2D* HoughHistX = new TH2D ("houghx", "Histogram for Hough transform", NHoughBinsA, 0, M_PI, NHoughBinsR, -HoughR, HoughR);
    TH2D* HoughHistY = new TH2D ("houghy", "Histogram for Hough transform", NHoughBinsA, 0, M_PI, NHoughBinsR, -HoughR, HoughR);

        if(TestMode)cout<<"Just loaded event, event:\t"<<current_event_ID<<", nx: "<<n_x_SiPM<<", ny: "<<n_y_SiPM<<endl;
        
       // if (n_sipm != n_x_SiPM+n_y_SiPM && phe_threshold == 0)std::cerr<<"\nError with SiPM number\t"<<current_event_ID<<"\n";
        //TODO
        
        
        if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end())
        {
            Deleter<TH2D> sipmx = new TH2D("sipmx", "", 25, 0, 100, 100, 0, 100);
            Deleter<TH2D> sipmy = new TH2D("sipmy","" , 25, 0, 100, 100, 0, 100);
            Deleter<TCanvas> c11 = new TCanvas("c11", "Some name", 30, 40, 900, 900);
            c11->cd();
            sipmx->Reset();
            sipmy->Reset();
            /*TGraph* gz = new TGraph(arrlen, ysipmarr, zysipmarr);
            gz->Draw("AP*");
            gz->GetXaxis()->SetLimits(0, 100);
            gz->GetYaxis()->SetLimits(0, 100);
            gz->SetTitle(DoubleToTitleName(curevid, 0).c_str());
            c1->SaveAs("./output_data/ygr.png");*/
            
            // sipmx->SetTitle(DoubleToTitleName(current_event_ID, 0, -10, -10).c_str());
            for (int iter = 0; iter < x_strip_vector.size(); iter++)
            {   
                sipmx->Fill(x_strip_vector[iter].XY,x_strip_vector[iter].Z, x_strip_vector[iter].energy);
                
            }
            sipmx->Draw("COLZ");
            sipmx->SetStats(0);
            // c11->SaveAs(DoubleToPNGName(current_event_ID, 0, MCflag).c_str());    
            
            
           
           // cout<<current_event_ID<<"\n";
            //delete c11;
            
        }
        if (n_x_SiPM < 5 || n_y_SiPM < 5){
            delete HoughHistX;
            delete HoughHistY;
            return 5;
            
        }
     
        // if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
        //                 SavePicturePrim(n_x_SiPM, current_event_ID, 0, MCflag, x_strip_vector, 0, 0, 0, 0, -10, -10, 0, "0start");
        //                 SavePicturePrim(n_y_SiPM, current_event_ID, 1, MCflag, y_strip_vector, 0, 0, 0, 0, -10, -10, 0, "0start");
        //             }    
        
        //HoughTransform (HoughHistX, xstripvec, &px0, &px1);
        RotateAndHoughTransform(HoughHistX, x_strip_vector, &track.px0, &track.px1);
        /*c2->cd();
        HoughHistX->Draw("COLZ");
        HoughHistX->SetStats(0);
        c2->SaveAs(DoubleToPNGName(curevid, 0, MCflag, "hough").c_str());*/
       
        
        
        fx->SetParameters(track.px0, track.px1);        
        
        
        
        RotateAndHoughTransform (HoughHistY, y_strip_vector, &track.py0, &track.py1);
        //HoughTransform (HoughHistY, ystripvec, &py0, &py1);
    
        
        fy->SetParameters(track.py0, track.py1);  
        
        // if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
        //                 SavePicturePrim(n_x_SiPM, current_event_ID, 0, MCflag, x_strip_vector, track.px0, track.px1, 0, 0, -10, -10, 0, "1houghline");
        //                 SavePicturePrim(n_y_SiPM, current_event_ID, 1, MCflag, y_strip_vector, track.py0, track.py1, 0, 0, -10, -10, 0, "1houghline");
        //             } 
        
        int n_x_SiPM_old = n_x_SiPM;
        int n_y_SiPM_old = n_y_SiPM;
        
//         n_X_SiPM_all_hist->Fill(n_x_SiPM);
//         n_Y_SiPM_all_hist->Fill(n_y_SiPM);
        
        xfit = 0;
        yfit = 0;
        
//         if(TestMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
//             cout<<"Event ID: "<<current_event_ID<<endl;
//             n_x_SiPM =  FindMaxConnectedSpace (x_strip_vector, *px0, *px1, n_X_SiPM_neighbor_hist, true);
//         }
//         else
            n_x_SiPM =  FindMaxConnectedSpace (x_strip_vector, track.px0, track.px1);
        
        n_y_SiPM =  FindMaxConnectedSpace (y_strip_vector, track.py0, track.py1);
        //if (options.TestMode && curevid<100) SavePicturePrim(nxsipmold, curevid, 0, MCflag, xstripvec, px0, px1, 0, 0, track.chiX, track.ndfX, 0, "NeighborsTest");
        if(TestMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end())
            cout<<"Before first cut, event:\t"<<current_event_ID<<", nx: "<<n_x_SiPM<<", ny: "<<n_y_SiPM<<endl;
        //First cut in fact!
        
        if (n_x_SiPM<5){ 
            delete HoughHistX;
            delete HoughHistY;
            return 5;
        }
        //if(options.TestMode)cout<<"First cut nxsipm\t"<<curevid<<endl;
        
       
        
        if (n_y_SiPM<5){ 
            delete HoughHistX;
            delete HoughHistY;
            return 5;
        }
        //if(options.TestMode)cout<<"First cut nysipm\t"<<curevid<<endl;
       
        // if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
        //     SavePicturePrim(n_x_SiPM_old, current_event_ID, 0, MCflag, x_strip_vector, 0, 0, 0, 0, track.chi_x, track.ndf_x, 0, "2Neighbors");
        //     SavePicturePrim(n_y_SiPM_old, current_event_ID, 1, MCflag, y_strip_vector, 0, 0, 0, 0, track.chi_y, track.ndf_y, 0, "2Neighbors");
        // }
        
       
        deltaX = CutDelta(current_event_ID, x_strip_vector,0);
        deltaXZ = CutDelta(current_event_ID, x_strip_vector,1);
        deltaY = CutDelta(current_event_ID, y_strip_vector,0);
        deltaYZ = CutDelta(current_event_ID, y_strip_vector,1);
        
//         delta_x_hist->Fill(deltaX);
//         delta_xz_hist->Fill(deltaXZ);
//         
//         delta_y_hist->Fill(deltaY);
//         delta_yz_hist->Fill(deltaYZ);
//         
        
        //cout<<"Before delta cut\n";
        if (deltaX<30 && deltaXZ<20){ 
            delete HoughHistX;
            delete HoughHistY;
            return 5;
            
        }
        if (deltaY<30 && deltaYZ<20){ 
            delete HoughHistX;
            delete HoughHistY;
            return 5;
        }
        
       
        //cout<<"Third cut (deltaX<30 && deltaXZ<20)\t"<<curevid<<endl;
        
        //if (curevid  < 20)SavePicturePrim(nxsipmold, curevid, 0, MCflag, xstripvec, px0, px1, 0, 0, track.chiX, track.ndfX, 0, "1_" );
        
        //Fit
        xfit = FFitErr(current_event_ID, n_x_SiPM, 0, x_strip_vector, &track.px0, &track.px1, &track.chi_x, &track.ndf_x, &track.vertical_x);
            
        yfit = FFitErr(current_event_ID, n_y_SiPM, 1, y_strip_vector, &track.py0, &track.py1, &track.chi_y, &track.ndf_y, &track.vertical_y);
        
//        if (curevid == 60 || curevid == 68 || curevid == 72 || curevid == 40){
//                         SavePicturePrim(nxsipmold, curevid, 0, MCflag, xstripvec, px0, px1, 0, 0, chi2x, ndfx, 0, "after_first_fit");
//                         SavePicturePrim(nysipmold, curevid, 1, MCflag, ystripvec, py0, py1, 0, 0, chi2y, ndfy, 0, "after_first_fit");
//                     } 
//        
       //Fit done, looking for good and bad strips
        for (int i = 0; i<x_strip_vector.size(); i++)
        {
            //point_line_hist_x->Fill(PointLineDistance(x_strip_vector[i].XY,x_strip_vector[i].Z, px0, px1));
            if ((PointLineDistance(x_strip_vector[i].XY,x_strip_vector[i].Z, track.px0, track.px1) < (distance_cut-track.vertical_x*0.5)) && x_strip_vector[i].status==0)
            {
                x_strip_vector[i].status = 1;
                n_x_SiPM++;
            }
            if((PointLineDistance(x_strip_vector[i].XY,x_strip_vector[i].Z, track.px0, track.px1) > (distance_cut-track.vertical_x*0.5)) && x_strip_vector[i].status!=0)//TODO 
            {
                x_strip_vector[i].status = 0;
                n_x_SiPM--;
            }
           
        }
       // n_X_SiPM_track_hist->Fill(n_x_SiPM);
        
        for (int i = 0; i<y_strip_vector.size(); i++)
        {
           // point_line_hist_y->Fill(PointLineDistance(y_strip_vector[i].XY,y_strip_vector[i].Z, py0, py1));
            if ((PointLineDistance(y_strip_vector[i].XY,y_strip_vector[i].Z, track.py0, track.py1) < (distance_cut-track.vertical_y*0.5)) && y_strip_vector[i].status==0)
            {
                y_strip_vector[i].status = 1;
                n_y_SiPM++;
            }
            if((PointLineDistance(y_strip_vector[i].XY,y_strip_vector[i].Z, track.py0, track.py1) > (distance_cut-track.vertical_y*0.5)) && y_strip_vector[i].status!=0)//TODO 
            {
                y_strip_vector[i].status = 0;
                n_y_SiPM--;
            }
           
        }
       
      // n_Y_SiPM_track_hist->Fill(n_y_SiPM);
       
    //    if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
    //                     SavePicturePrim(n_x_SiPM_old, current_event_ID, 0, MCflag, x_strip_vector, track.px0, track.px1, 0, 0, track.chi_x, track.ndf_x, 0, "3befor_second_fit");
    //                     SavePicturePrim(n_y_SiPM_old, current_event_ID, 1, MCflag, y_strip_vector, track.py0, track.py1, 0, 0, track.chi_y, track.ndf_y, 0, "3befor_second_fit");
    //                 } 
                    
       //And fit again    
       //cout<<"n SiPM: "<<n_x_SiPM<<", vector size: "<<x_strip_vector.size()<<"\n";
       xfit = FFitErr(current_event_ID, n_x_SiPM, 0, x_strip_vector, &track.px0, &track.px1, &track.chi_x, &track.ndf_x, &track.vertical_x);
            
       yfit = FFitErr(current_event_ID, n_y_SiPM, 1, y_strip_vector, &track.py0, &track.py1, &track.chi_y, &track.ndf_y, &track.vertical_y);
        
       
    //    if (DrawMode && std::find (interesting_events.begin(), interesting_events.end(), current_event_ID)!= interesting_events.end()){
    //                     SavePicturePrim(n_x_SiPM_old, current_event_ID, 0, MCflag, x_strip_vector, track.px0, track.px1, 0, 0, track.chi_x, track.ndf_x, 0, "4after_second_fit");
    //                     SavePicturePrim(n_y_SiPM_old, current_event_ID, 1, MCflag, y_strip_vector, track.py0, track.py1, 0, 0, track.chi_y, track.ndf_y, 0, "4after_second_fit");
                        
    //                     cout<<"vertical status: "<<track.vertical_x<<", " <<track.vertical_y<<"\n";
    //                 }
                    
   

    delete HoughHistX;
    delete HoughHistY;
    if(Int_t(xfit)==0 && Int_t(yfit)==0) return 0;
    else return 5;

}

bool StripCondition (int SiPMID)
{
    //return (((int(SiPMID)/25)+0.5>up || (int(SiPMID)/25)+0.5<low) && side==0);
    //if (SiPMZCoordinate(SiPMID)>=52) return true;
    //if (((SiPMID)/25)%2==1) return true;
    //else return false;
    return true;
    
}
double AttenuationCorrectionSiPM(double distance, bool mode){
    //<func(x)=1>
    if (mode) return 1;
    const double FuncAverage = 1.02208;
    double rez;
    double x = distance;
    rez = (0.00577381*exp(-0.1823*(x-50))+0.999583*exp(-0.0024*(x-50))-8.095E-13*exp(0.5205*(x-50))-0.00535714*exp(-0.1838*(x-50)))
        / FuncAverage;
    return rez;

}
void FillEnergyLengthHists(std::vector <TH1D>& energy_hist_vec, std::vector <TH1D*> length_vec, std::vector <double> left_border_l, std::vector <double> right_border_l, double energy_sipm, double length_sipm, int SiPMID){
    
    for (int i = 0; i < energy_hist_vec.size(); i++){
        if ((length_sipm>=left_border_l[i]) && (length_sipm<right_border_l[i]) && StripCondition(SiPMID)){
            length_vec[i]->Fill(length_sipm);
            energy_hist_vec[i].Fill(energy_sipm);
        }
        
    }
//     if ((length_sipm>=1.0) && (length_sipm<1.2) && StripCondition(SiPMID)){
//         length_vec[0]->Fill(length_sipm);
//         energy_hist_vec[0].Fill(energy_sipm);
//     }
//     if ((length_sipm>=1.2) && (length_sipm<1.4) && StripCondition(SiPMID)){
//         length_vec[1]->Fill(length_sipm);
//         energy_hist_vec[1].Fill(energy_sipm);
//     }
//     if ((length_sipm>=1.4) && (length_sipm<1.6) && StripCondition(SiPMID)){
//         length_vec[2]->Fill(length_sipm);
//         energy_hist_vec[2].Fill(energy_sipm);
//     }
//     if ((length_sipm>=1.6) && (length_sipm<1.8) && StripCondition(SiPMID)){
//         length_vec[3]->Fill(length_sipm);
//         energy_hist_vec[3].Fill(energy_sipm);
//     }
//     if ((length_sipm>=1.8) && (length_sipm<2.) && StripCondition(SiPMID)){
//         length_vec[4]->Fill(length_sipm);
//         energy_hist_vec[4].Fill(energy_sipm);
//     }
//     if ((length_sipm>=2) && (length_sipm<2.5) && StripCondition(SiPMID)){
//         length_vec[5]->Fill(length_sipm);
//         energy_hist_vec[5].Fill(energy_sipm);
//     }
//     if (energy_hist_vec.size() > 6){
//     if ((length_sipm>=2.5) && (length_sipm<3. && StripCondition(SiPMID))){
//         length_vec[6]->Fill(length_sipm);
//         energy_hist_vec[6].Fill(energy_sipm);
//     }
//     if ((length_sipm>=3.) && (length_sipm<4. && StripCondition(SiPMID))){
//         length_vec[7]->Fill(length_sipm);
//         energy_hist_vec[7].Fill(energy_sipm);
//     }
//     if ((length_sipm>=4.) && (length_sipm<5.5 && StripCondition(SiPMID))){
//         length_vec[8]->Fill(length_sipm);
//         energy_hist_vec[8].Fill(energy_sipm);
//     }
//     if ((length_sipm>=5.5) && (length_sipm<7. && StripCondition(SiPMID))){
//         length_vec[9]->Fill(length_sipm);
//         energy_hist_vec[9].Fill(energy_sipm);
//     }
//         
//     }
}

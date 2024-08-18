#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

#include "data.h"
#include "mupro.h"


using std::cout;
using std::endl;
using std::cin;


double PrimX(double y) {
    return y/1e4+50;
}

double PrimY(double x) {
    return 100-(x/1e4+50);
}

double PrimZ(double z) {
    return z/1e4+52;
}

int AlexSiPMNumbering(int SiPMID) {//TODO
    int row = -10;
    int column = -10;
    if (((SiPMID)/25)%2==0) { //my x, Sasha y
            column = SiPMID%25; //as is
            row = SiPMID/25; 
        }
        else {
            column = 24-SiPMID%25; // flip
            row = SiPMID/25;   
        }
    SiPMID = 25*row+column; 
    return SiPMID;       
}

double SiPMZCoordinate(int SiPMID) {
    return double(SiPMID/25)*(LengthZ+GapLength)+HalfLengthZ;   
}

double SiPMXYCoordinate(int SiPMID) {
    return double((SiPMID%25)*LengthXY+HalfLengthXY);   
}

double meanXY (int n, double * x, double * y) {
    //int n = 2500;
    double RealLength = 0;
    double sum = 0;
    for (int i = 0; i < n; i++) {
        if ((x[i] != -10) && (y[i] != -10)) {
            sum+=(x[i]*y[i]);
            RealLength++;
        }

    }
    return sum/RealLength;
}

double meanX(int n, double *x) {
    //int n = 2500;
    double RealLength = 0;
    double sum = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] != -10) {
            sum+=(x[i]);
            RealLength++;
        }

    }
    return sum/RealLength;
}

double mnk_p0 (int n, double * x, double * y) {
     return (meanXY(n, x, y)- meanX(n, x)*meanX(n, y))/(meanXY (n, x,x)-meanX(n,x)*meanX(n,x));
}

double mnk_p1(int n, double * x, double * y) {
    return meanX(n,y) - mnk_p0(n,x,y)*meanX(n,x);
}

int RotateAndFit (int eventID, int nsize,  double *old_x, double *old_z, double *old_p0, double *old_p1) {
    
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
    for (i = 0; i < nsize; i++) {
        if ((old_x[i] != -10) && (old_z[i] != -10)) {
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
    delete []x;
    delete []z;
    return ret;  //???????
}

int FFit(int eventID, int n, double * x, double * z, double *p0, double *p1) {
    if (fabs(*p0) < 1) {
        *p0 = mnk_p0(n, x, z);
        *p1 = mnk_p1(n, x, z);
        return 0;
    }
    else if (fabs(*p0)>=1)return RotateAndFit (eventID, n, x, z, p0, p1);
    else return -10;
}

int DistanceX (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout) {
    //double distance = -10;
    //double z1 = double((SiPMID-1)/25);
    //p0 = p0;
    *yin =  (*zin-py1)/py0;
    *yout = (*zout-py1)/py0;
    if ((*yin < 0) || (*yin >100) || (*yout < 0) || (*yout >100)) return -10;
    else return 0;//Status OK
}

int DistanceY (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout) {
    *xin =  (*zin-px1)/px0;
    *xout = (*zout-px1)/px0;
    if ((*xin < 0) || (*xin >100) || (*xout < 0) || (*xout >100)) return -10;
    else return 0;//Status OK
}

int Cross(double x1, double x2, double z1, double z2, double p0, double p1, double *xin, double *xout, double *zin, double *zout) {
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

    if ((zcross1>z1)&&(zcross1 < z2)) { //left hit
        if ((xcross2>x1)&&(xcross2<x2)) {
            *xin = x1;
            *xout = xcross2;
            *zin = zcross1;
            *zout = z2;
            return 0;
            //  ____
            // |
        }
        else if ((zcross2<z2)&&(zcross2>z1)) {
            *xin = x1;
            *xout = x2;
            *zin = zcross1;
            *zout = zcross2;
            return 0;
            // |    |
        }
        else if ((xcross1>x1)&&(xcross1<x2)) {
            *xin = x1;
            *xout = xcross1;
            *zin = zcross1;
            *zout = z1;
            return 0;
            // |____
        }
    }
    else if ((xcross2>x1)&&(xcross2<x2)) {  //up hit
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
        else if ((xcross1>x1)&&(xcross1<x2)) {
            *xin = xcross1;
            *xout = xcross2;
            *zin = z1;
            *zout = z2;
            return 0;
            //  ____
            //  ____
        }
    }
    else if ((zcross2<z2)&&(zcross2>z1)) {  //right hit
        if ((xcross1>x1)&&(xcross1<x2)) {
            *xin = xcross1;
            *xout = x2;
            *zin = z1;
            *zout = zcross2;
            return 0;
            //  ____|
        }
    }
    if ((zcross1 == z2) && (xcross2 == x1)) {  //up left corner
        cout << "Up left\n";
        if (p0 > 0) return -10;
        if ((zcross2 == z1) && (xcross1 == x2)) {  //down right corner
            cout<<"And down right too\n";
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross1;
            *zout = zcross2;
            return 0;
        }
        else if ((zcross2<z2)&&(zcross2>z1)) {  //right hit
            *xin = xcross2;
            *xout = x2;
            *zin = zcross1;
            *zout = zcross2;
        }
        else if ((xcross1>x1)&&(xcross1<x2)) {  //down hit
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross1;
            *zout = z1;
            return 0;
        }   
        else return -10;
    }
    else if ((zcross2 == z2) && (xcross2 == x2)) {  //up right corner
        cout << "Up right\n";
        if (p0 < 0) return -10;
        if ((zcross1 == z1) && (xcross1 == x1)) {  //down left corner
            cout << "And down left too\n";
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((zcross1>z1) && (zcross1 < z2)) {  //left hit
            *xin = xcross2;
            *xout = x1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((xcross1>x1)&&(xcross1<x2)) {  //down hit
            *xin = xcross2;
            *xout = xcross1;
            *zin = zcross2;
            *zout = z1;
            return 0;
        }
        else return -10;
    }
    else if ((zcross2 == z1) && (xcross1 == x2)) {  //down right
        cout << "Down right\n";
        if (p0 > 0) return -10;
        else if ((zcross1>z1) && (zcross1 < z2)) {  //left hit
            *xin = xcross1;
            *xout = x1;
            *zin = zcross2;
            *zout = zcross1;
            return 0;
        }
        else if ((xcross2>x1) && (xcross2<x2)) {  //up hit
            *xin = xcross1;
            *xout = xcross2;
            *zin = zcross2;
            *zout = z2;
            return 0;
        }
        else return -10;
        
    }
    
    else if ((zcross1 == z1) && (xcross1 == x1)) {  //down left
        cout << "Down left\n";
        if (p0 < 0) return -10;
        else if ((xcross2>x1) && (xcross2<x2)) {  //up hit
            *xin = xcross1;
            *xout = xcross2;
            *zin = zcross1;
            *zout = z2;
            return 0;
        }
        else if ((zcross2<z2)&&(zcross2>z1)) {  //right hit
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

int LengthX (int SiPMID, double px0, double px1, double *xin, double *xout, double *zin, double *zout) {
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

int LengthY (int SiPMID, double py0, double py1, double *yin, double *yout, double *zin, double *zout) {
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

double Length(int SiPMID, double px0, double px1, double py0, double py1, double *xin, double *xout, double *yin, double *yout, double *zin, double *zout) {
    double cross = -10;
    double d_status = -10;
    if (((SiPMID)/25)%2==0) { //x side
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

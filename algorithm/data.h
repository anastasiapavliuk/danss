#pragma once


const int NHoughBinsR = 576; //144*4
const int NHoughBinsA = 400;
const int InterestingSiPMID = 413; //FIXME
#define HoughR 144 //double or int?
//144 == sqrt (100^2+104^2)
const double FitPrecision = 1e-9;
//#define StepSize 1e-3



#include "TBranch.h"
#include "TFile.h"
#include "TLeafF.h"
#include "TLeafI.h"
#include "TLeafL.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <sys/types.h>
#include <dirent.h>
#include "danss_constants.h"

const double LengthXY = 4;
const double HalfLengthXY = 2;
const double LengthZ = 1;
const double HalfLengthZ = 0.5;
const int NEnergyBins = 100;
const double NoSenseLength = 0.03;

const double GapLength  = 4.04040404040404075e-02; // == 4/99
const double DANSSxSize = 100;
const double DANSSySize = 100;
const double DANSSzSize = 104;
const double GapEpsilon = 1e-10;

const double theory_slope_width = 0.307075*0.53768*2*1.06;

struct TStrip {
  Int_t SiPMID;
  Double_t XY;
  Double_t Z;
  Double_t energy;
  Int_t status;

  static Int_t GetZLevelBySiPMID(Int_t SiPMID);
  static Double_t GetZBySiPMID(Int_t SiPMID);
  static Double_t GetXYBySiPMID(Int_t SiPMID);
  static bool IsSideXZ(Int_t SiPMID);
};

struct PMTID {
    Int_t XY;
    Int_t Z;
    
    bool IsSideXZ() const{
    return Z % 2 == 0;
}
    bool ContainsSiPM(Int_t sipm_id) const{
    static constexpr Double_t strips_in_pmt_vertical = static_cast<Double_t>(DANSSConstants::StripsTotalVertical) / 
                                                        DANSSConstants::PMTsTotalVertical;
    static constexpr Double_t pmt_height = (strips_in_pmt_vertical - 1) * DANSSConstants::StripSizeVertical +
                                            (strips_in_pmt_vertical - 2) * DANSSConstants::GapVertical;
    static constexpr Double_t pmt_double_height = strips_in_pmt_vertical * 
                                                    (DANSSConstants::StripSizeVertical + DANSSConstants::GapVertical);
    static constexpr Double_t pmt_width = DANSSConstants::DetectorSizeHorizontal / DANSSConstants::PMTsTotalHorizontal;
    
    if (TStrip::IsSideXZ(sipm_id) != IsSideXZ()) {
        return false;
    }
    Double_t xy_left = XY * pmt_width;
    Double_t z_bottom = Z / 2 * pmt_double_height;
    if (!IsSideXZ()) {
        z_bottom += DANSSConstants::StripSizeVertical + DANSSConstants::GapVertical;
    }
    
    Double_t sipm_xy = TStrip::GetXYBySiPMID(sipm_id);
    if (sipm_xy < xy_left || sipm_xy > xy_left + pmt_width) {
        return false;
    }

    Double_t sipm_z = TStrip::GetZBySiPMID(sipm_id);
    if (sipm_z < z_bottom || sipm_z > z_bottom + pmt_height) {
        return false;
    }
    return true;
}
};



struct PMTInfo {
    PMTID id;
    Double_t energy;
};

struct EventInfo {
    Long64_t event_id;
    Long64_t global_time;
    Long64_t veto_hits;
    Double_t veto_energy;
    std::vector<TStrip> strips_xz;
    std::vector<TStrip> strips_yz;
    std::vector<PMTInfo> pmts_xz;
    std::vector<PMTInfo> pmts_yz;
};

struct MCEventInfo {
    Long64_t event_id;
    Long64_t global_time;
    Long64_t veto_hits;
    Double_t veto_energy;
    std::vector<TStrip> strips_xz;
    std::vector<TStrip> strips_yz;
    std::vector<PMTInfo> pmts_xz;
    std::vector<PMTInfo> pmts_yz;
    Double_t birth_x;
    Double_t birth_y;
    Double_t birth_z;
};
struct OutStrip {
    int64_t SiPMID;
    double row;
    double column;
    double energy;
};


struct TrackStruct {
    double chi_x;
    double chi_y;
    double ndf_x;
    double ndf_y;
    double px0;
    double py0;
    double px1;
    double py1;
    double theta;
    int64_t eventID;
    int64_t RealEventID;
    int vertical_x;
    int vertical_y;
};

struct OutputSiPMStruct{
    int64_t eventID;
    int64_t RealEventID;
    int SiPMID;
    int side; // 0 -- X, 1 -- Y
    double energy;
    double distance;
    double length;
    double width;
    double px0;
    double py0;
    double px1;
    double py1;
    double chiX;
    double chiY;
    double ndfX;
    double ndfY;
    int vertical_x;
    int vertical_y;
};

struct OutputPMTStruct{
    int64_t eventID;
    int64_t RealEventID;
    int side; // 0 -- X, 1 -- Y
    int XY;
    int Z;    
    int nSiPMraw;
    int vertical_x;
    int vertical_y;
    double energy;
    double energySiPMraw;
    double energySiPMtrack;
    double lengthStripSum;

};

template<class T>
class Deleter {
public:
    Deleter(T* ptr = 0)
        : Ptr(ptr)
    {
        //
    }
    T* operator->() {
        return Ptr;
    }
    ~Deleter() {
    
        delete Ptr;
    }
private:
    T* Ptr;
    Deleter<T>& operator=(const Deleter<T>& other);
};

#ifndef DANSS_CONSTANTS_H_INCLUDED
#define DANSS_CONSTANTS_H_INCLUDED

//#include "RtypesCore.h"

namespace DANSSConstants {
    static constexpr Double_t DetectorSizeHorizontal = 100;
    static constexpr Double_t DetectorSizeVertical = 104;
    static constexpr Int_t StripsTotalHorizontal = 25;
    static constexpr Int_t StripsTotalVertical = 100;
    static constexpr Double_t StripSizeHorizontal = DetectorSizeHorizontal / StripsTotalHorizontal;
    static constexpr Double_t StripSizeVertical = 1;
    static constexpr Double_t GapVertical =
        (DetectorSizeVertical - StripSizeVertical * StripsTotalVertical) /
        (StripsTotalVertical - 1);
    static constexpr Int_t PMTsTotalHorizontal = 5;
    static constexpr Int_t PMTsTotalVertical = 5;
};

#endif // DANSS_CONSTANTS_H_INCLUDED

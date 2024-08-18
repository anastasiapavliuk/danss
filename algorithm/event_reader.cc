//#include "data.h"
//#include "danss_constants.h"
//#include "muprocessing.h"
#include "event_reader.h"

const std::string AlekseevEventReader::kNameTree = "DanssEvent";
const std::string AlekseevEventReader::kNameBranchData = "Data";
const std::string AlekseevEventReader::kNameBranchHitType = "HitType";
const std::string AlekseevEventReader::kNameBranchHitEnergy = "HitE";
const std::string AlekseevEventReader::kNameLeafHitsTotal = "NHits";
const std::string AlekseevEventReader::kNameLeafGlobalTime = "globalTime";

Double_t TStrip::GetZBySiPMID(Int_t SiPMID) {
    return (SiPMID / DANSSConstants::StripsTotalHorizontal)
           * (DANSSConstants::StripSizeVertical + DANSSConstants::GapVertical)
           + DANSSConstants::StripSizeVertical / 2;
}
Int_t TStrip::GetZLevelBySiPMID(Int_t SiPMID) {
    return (SiPMID / DANSSConstants::StripsTotalHorizontal)
           * (DANSSConstants::StripSizeVertical)
           + DANSSConstants::StripSizeVertical / 2;
}
Double_t TStrip::GetXYBySiPMID(Int_t SiPMID) {
    return (SiPMID % DANSSConstants::StripsTotalHorizontal)
           * DANSSConstants::StripSizeHorizontal
           + DANSSConstants::StripSizeHorizontal / 2;
}

bool TStrip::IsSideXZ(Int_t SiPMID) {
  return (SiPMID / DANSSConstants::StripsTotalHorizontal) % 2 == 0;
}

AlekseevEventReader::AlekseevEventReader(const std::vector<std::string>& file_names)
    : chain_(new TChain(kNameTree.c_str())) {

    for (const auto& file_name : file_names) {
        chain_->Add(file_name.c_str());
    }
}

Long64_t AlekseevEventReader::GetEntries() const {
    return chain_->GetEntries();
}

EventInfo AlekseevEventReader::ReadEvent(Long64_t event_id) {
    chain_->GetEntry(event_id);
    EventInfo event_info;
    event_info.event_id = event_id;
    TLeafL* leaf_global_time = static_cast<TLeafL*>
        (chain_->GetLeaf(kNameBranchData.c_str(), kNameLeafGlobalTime.c_str()));
    TLeafI* leaf_hits_total = static_cast<TLeafI*>
        (chain_->GetLeaf(kNameBranchData.c_str(), kNameLeafHitsTotal.c_str()));
    TLeafI* leaf_hit_type = static_cast<TLeafI*>
        (chain_->GetLeaf(kNameBranchHitType.c_str(), kNameBranchHitType.c_str()));
    TLeafF* leaf_hit_energy = static_cast<TLeafF*>
        (chain_->GetLeaf(kNameBranchHitEnergy.c_str(), kNameBranchHitEnergy.c_str()));
    event_info.global_time = leaf_global_time->GetValue();
    size_t hits_total = leaf_hits_total->GetValue();
    for (size_t hit_id = 0; hit_id < hits_total; ++hit_id) {
        struct HitTypeStruct {
            char type;
            char z;
            char xy;
            char flag;
        } hit_type_struct;
        Int_t hit_type_value = leaf_hit_type->GetValue(hit_id);
        memcpy(&hit_type_struct, &hit_type_value, sizeof(hit_type_value));
        if (hit_type_struct.z < 0 || hit_type_struct.z >= 100) {
            throw std::runtime_error("Bad z value");
        }
        if (hit_type_struct.xy < 0 || hit_type_struct.xy >= 25) {
            throw std::runtime_error("Bad xy value");
        }
        if (hit_type_struct.type == kTypeSiPM) {
            TStrip strip;
            if(static_cast<Int_t>(hit_type_struct.z)%2==0){ //X 
                strip.SiPMID = static_cast<Int_t>(hit_type_struct.z) *
                DANSSConstants::StripsTotalHorizontal + ( DANSSConstants::StripsTotalHorizontal - 1 - 
                static_cast<Int_t>(hit_type_struct.xy));
                
            }
            else{
                strip.SiPMID = static_cast<Int_t>(hit_type_struct.z) *
                DANSSConstants::StripsTotalHorizontal +
                static_cast<Int_t>(hit_type_struct.xy);
            }
            strip.XY = TStrip::GetXYBySiPMID(strip.SiPMID);
            strip.Z = TStrip::GetZBySiPMID(strip.SiPMID);
            strip.energy = leaf_hit_energy->GetValue(hit_id);
            strip.status = 1;
            if (strip.energy < sipm_energy_threshold_) {
                continue;
            }
            if (TStrip::IsSideXZ(strip.SiPMID)) {
                event_info.strips_xz.push_back(strip);
            } else {
                event_info.strips_yz.push_back(strip);
            }
        } else if (hit_type_struct.type == kTypePMT) {
            PMTInfo pmt;
            pmt.id.XY = hit_type_struct.xy;
            pmt.id.Z = hit_type_struct.z;
            pmt.energy = leaf_hit_energy->GetValue(hit_id);
            if (pmt.id.Z % 2 == 0) { // TO FIX
                pmt.id.XY = DANSSConstants::PMTsTotalHorizontal - 1 - pmt.id.XY;
                event_info.pmts_xz.push_back(pmt);
            } else {
                event_info.pmts_yz.push_back(pmt);
            }
        }
    }
    return event_info;
}

void AlekseevEventReader::SetSiPMEnergyThreshold(Double_t threshold) {
    sipm_energy_threshold_ = threshold;
}

template<class FileNames>
std::unique_ptr<EventReader> InitReader(const FileNames& file_names,
    EventReader::FileType file_type) {

    switch(file_type) {
        case EventReader::FileType::Alekseev:
            return std::unique_ptr<EventReader>(new AlekseevEventReader(file_names));
        default:
            return nullptr;
        // case EventReader::FileType::Irina:
        //     return std::unique_ptr<EventReader>(new IrinaEventReader(file_names));
    }
    return nullptr;
}

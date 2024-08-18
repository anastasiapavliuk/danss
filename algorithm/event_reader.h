#ifndef EVENT_READER_H_INCLUDED
#define EVENT_READER_H_INCLUDED

#include "TChain.h"
#include "TFile.h"
#include "TLeafD.h"
#include "TLeafF.h"
#include "TLeafI.h"
#include "TLeafL.h"

#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "data.h"
#include "danss_constants.h"


class AlekseevEventReader;

class EventReader {
public:
    enum class FileType { Alekseev, Irina };

    virtual ~EventReader() = default;
    virtual Long64_t GetEntries() const = 0;
    virtual EventInfo ReadEvent(Long64_t event_id) = 0;
    virtual void SetSiPMEnergyThreshold(Double_t threshold) = 0;
};

class AlekseevEventReader : public EventReader {
public:
    explicit AlekseevEventReader(const std::vector<std::string>& file_names);
    explicit AlekseevEventReader(const std::string& file_name)
        : AlekseevEventReader(std::vector<std::string>{file_name}) {
    }
    Long64_t GetEntries() const override;
    EventInfo ReadEvent(Long64_t event_id) override;
    void SetSiPMEnergyThreshold(Double_t threshold) override;

private:
    static const std::string kNameTree;  // = "DanssEvent";
    static const std::string kNameBranchData;  // = "Data";
    static const std::string kNameBranchHitType;  // = "HitType";
    static const std::string kNameBranchHitEnergy;  // = "HitE";
    static const std::string kNameLeafHitsTotal;  // = "NHits";
    static const std::string kNameLeafGlobalTime;  // = "globalTime";
    static const std::string kNameLeafVetoHits;  // = "VetoCleanHits";
    static const std::string kNameLeafVetoEnergy;  // = "VetoCleanEnergy";
    static const std::string kNameBranchMC;  // = "MC";
    static const std::string kNameLeafMCX;  // = "McX";
    static constexpr Int_t kTypeSiPM = 0;
    static constexpr Int_t kTypePMT = 1;
    
    Double_t sipm_energy_threshold_ = 0;

    std::unique_ptr<TChain> chain_;
};

template<class FileNames>
std::unique_ptr<EventReader> InitReader(const FileNames& file_names,
    EventReader::FileType file_type);

#endif // EVENT_READER_H_INCLUDED

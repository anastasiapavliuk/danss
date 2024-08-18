#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include "./json.hpp"
#include <algorithm>

#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeafD.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "data.h"
#include "danss_constants.h"
#include "event_reader.h"
#include "event_reader.cc"
#include "mupro.h"
#include "langaus.h"

#include "TStyle.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TMath.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

Double_t n_pix(Double_t energy) {
    Double_t n_pix;
    n_pix = 512 * (1 - exp(-(energy*18.67/512)));
    return n_pix;
}

void MNK (vector<TStrip> v, Double_t *k, Double_t *c) {

    //x = k*z + c

    Int_t n = v.size();

    Double_t s_x = 0;
    Double_t s2_z = 0;
    Double_t s_z = 0; 
    Double_t s_xz = 0;

    for (Int_t i = 0; i < n; i++) {
        s_x += v[i].XY;
        s2_z +=v[i].Z*v[i].Z;
        s_z += v[i].Z;
        s_xz += v[i].XY*v[i].Z;
    }

    *k = (n*s_xz - s_x*s_z) / (n*s2_z - s_z*s_z);
    *c = (s_x - *k*s_z) / n;
}

Double_t everage_div (vector<TStrip> v, Double_t k, Double_t c) {

    Double_t everage_div = 0;
    Int_t n = v.size();
    for (Int_t i = 0; i < n; i++) {
        everage_div += abs(v[i].Z - v[i].XY*k - c);
    }
    everage_div /=n;
    return everage_div;
}

Double_t dist (Double_t k, Double_t c, Double_t z0, Double_t x0) {
    Double_t dist;
    dist = sqrt(pow( ((z0+k*x0-k*c) / (pow(k, 2)+1) - z0), 2 ) + pow( (k*(z0+k*x0-k*c) / (pow(k , 2)+1) + c - x0), 2));
    return dist;
}

void write (std::string filename, vector<TStrip> xz, vector<TStrip> yz, Double_t k_x, Double_t c_x, Double_t k_y, Double_t c_y) {
    int size_x = xz.size();
    int size_y = yz.size();
    
    vector <double_t> x;
    vector <double_t> z_x;
    vector <double_t> e_x;
    vector <double_t> y;
    vector <double_t> z_y;
    vector <double_t> e_y;


    for (Int_t i=0; i<size_x; i++) {
        x.push_back(xz[i].XY);
        z_x.push_back(xz[i].Z);
        e_x.push_back(xz[i].energy);
    }

    for (Int_t j=0; j<size_y; j++) {
        y.push_back(yz[j].XY);
        z_y.push_back(yz[j].Z);
        e_y.push_back(yz[j].energy);
    }

    nlohmann::json j;
    j["x"]   = x;
    j["z_x"] = z_x;
    j["energy_x"] = e_x;
    j["y"]   = y;
    j["z_y"] = z_y;
    j["energy_y"] = e_y;
    j["k_x"] = k_x;
    j["c_x"] = c_x;
    j["k_y"] = k_y;
    j["c_y"] = c_y;

    std::ofstream o(filename);
    o << j;
}




int main(int argc, char* argv[]) 
{

    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // обрабатываем обработанный мк -- строим траектории по стрипам:

    std::stringstream filename;
    filename << "/data/mc_Muons_indLY_transcode_rawProc_pedSim_01_01.root";

    TFile   *file0 = new TFile("/data/mc_Muons_indLY_transcode_rawProc_pedSim_01_01.root", "READ");
    TTree   *evtree0 = (TTree*) file0->Get("DanssEvent");
    TBranch *evbranch0 = evtree0->GetBranch("MCEvent");
    TLeafD  *eventIDleaf0 = (TLeafD*) evbranch0->GetLeaf("EventID");
    Int_t nevent0 = evtree0->GetEntries();

    AlekseevEventReader reader = AlekseevEventReader(filename.str());
    Long64_t Nevent = reader.GetEntries(); 
    Long64_t event_id = 547;                                                            // from 0 to Nevent
    Int_t EventID;
    vector<int> goodeventid;

    for (Long64_t id = 1; id <= nevent0; id++) {

        evbranch0->GetEntry(id);
        EventID = eventIDleaf0->GetValue();

        Long64_t global_time = reader.ReadEvent(EventID).global_time;
        Long64_t veto_hits = reader.ReadEvent(EventID).veto_hits;
        Double_t veto_energy = reader.ReadEvent(EventID).veto_energy;

        std::vector<TStrip> strips_xz = reader.ReadEvent(EventID).strips_xz;
        std::vector<TStrip> strips_yz = reader.ReadEvent(EventID).strips_yz;

        std::vector<PMTInfo> pmts_xz = reader.ReadEvent(EventID).pmts_xz;
        std::vector<PMTInfo> pmts_yz = reader.ReadEvent(EventID).pmts_yz;

        Int_t event_size_x = strips_xz.size();
        Int_t event_size_y = strips_yz.size();
        Int_t event_size = event_size_x + event_size_y;

        if (event_size > 0) {

            //cout << "event number: " << id << "; event size: " << event_size << endl;
            Int_t i=0; 
            while (i<event_size_x) {                                      //выкидываем шум (#1)
                if(n_pix(strips_xz[i].energy) < 3) {
                    //cout << "removing noise" << endl;
                    strips_xz.erase(strips_xz.begin()+i);
                    event_size_x -=1;
                }
                else {i++;}
            }
            i=0; 
            while (i<event_size_y) {
                if(n_pix(strips_yz[i].energy) < 3) {
                    //cout << "removing noise" << endl;
                    strips_yz.erase(strips_yz.begin()+i);
                    event_size_y -=1;
                }
                else {i++;}
            }

            //Time?

            event_size_x = strips_xz.size();
            event_size_y = strips_yz.size();
            event_size = event_size_x + event_size_y;

            if (event_size_x >= 5 && event_size_y >= 5 && event_size >= 20) {             //выкидываем мелкие события(#3)
                
                Double_t k_x0;
                Double_t c_x0;
                MNK(strips_xz, &k_x0, &c_x0);

                Double_t k_y0;
                Double_t c_y0;
                MNK(strips_yz, &k_y0, &c_y0);

                std::vector<TStrip> strips_eject_xz;
                std::vector<TStrip> strips_eject_yz;

                Int_t j = 0; 
                while (j < event_size_x) {                                                  // строим МНК и выкидываем выбивающиеся (#4)
                    if(dist(k_x0, c_x0, strips_xz[j].Z, strips_xz[j].XY) > 6) {             //(надеюсь, что strips_xz[j].XY в сантиметрах)
                        strips_eject_xz.emplace_back(strips_xz[j]);
                        strips_xz.erase(strips_xz.begin()+j);
                        event_size_x -=1;
                    }
                    else {j++;}
                }
                j = 0; 
                while (j < event_size_y) {                              
                    if(dist(k_y0, c_y0, strips_yz[j].Z, strips_yz[j].XY) > 6) {           
                        strips_eject_yz.emplace_back(strips_yz[j]);
                        strips_yz.erase(strips_yz.begin()+j);
                        event_size_y -=1;
                    }
                    else{j++;}
                }

                Double_t k_x1;
                Double_t c_x1;
                MNK(strips_xz, &k_x1, &c_x1);

                Double_t k_y1;
                Double_t c_y1;
                MNK(strips_yz, &k_y1, &c_y1);

                for (Int_t j = 0; j < strips_eject_xz.size(); j++) {                               // перестраиваем МНК и возвращаем  (#5)
                    if(dist(k_x1, c_x1, strips_xz[j].Z, strips_xz[j].XY) <= 6) {           
                        strips_xz.emplace_back(strips_eject_xz[j]);
                    }
                }

                for (Int_t j = 0; j < strips_eject_yz.size(); j++) {                              
                    if(dist(k_y1, c_y1, strips_yz[j].Z, strips_yz[j].XY) <= 6) {           
                        strips_yz.emplace_back(strips_eject_yz[j]);
                    }
                }

                event_size_x = strips_xz.size();
                event_size_y = strips_yz.size();
                event_size = event_size_x + event_size_y;



                // потом поудаляю все лишнее

                Double_t k_x;                                                                       // финальная итерация МНК (#6)
                Double_t c_x;
                MNK(strips_xz, &k_x, &c_x);

                Double_t k_y;
                Double_t c_y;
                MNK(strips_yz, &k_y, &c_y);

                // проверяем на отклонения (#7)

                if ((everage_div(strips_xz, k_x, c_x) < 1.4) && (everage_div(strips_yz, k_y, c_y) < 1.4)) {
                    cout << "draw hist in " << EventID << " event with size of " << event_size << endl;
                    goodeventid.push_back(EventID);
                    std::stringstream dataname;
                    dataname << "mce" << EventID << ".json";                           //file event << id
                    write(dataname.str(), strips_xz, strips_yz, k_x, c_x, k_y, c_y);
                }
                //else {cout << "bad event (adron)" << endl;}
            }
            //else {cout << "too few hits" << endl << endl;}
        }
        //else {cout << id <<" is null event" << endl << endl;}
    }

    
    //-----------------------------------------------------------------------------------------------------------------------------------------------------------------
    // обрабатываем сырой мс -- строим истинные треки

    
    TFile   *file = new TFile("/data/DANSSmod1_1.root", "READ");

    TTree   *primtree = (TTree*) file->Get("DANSSPrimary");
    TBranch *primbranch = primtree->GetBranch("PrimaryData");
    TLeafD  *primXleaf = (TLeafD*)primbranch->GetLeaf("X");
    TLeafD  *primYleaf = (TLeafD*)primbranch->GetLeaf("Y");
    TLeafD  *primZleaf = (TLeafD*)primbranch->GetLeaf("Z");
    TLeafD  *primeventIDleaf = (TLeafD*) primbranch->GetLeaf("EventID");
    Int_t nprim = primtree->GetEntries();
    // cout << "Primary entries: " << nprim << "\n";

    TTree   *evtree = (TTree*) file->Get("DANSSSignal");
    TBranch *evbranch = evtree->GetBranch("RealSignalData");
    TLeafD  *signalleaf = (TLeafD*) evbranch->GetLeaf("RealSignal");
    TBranch *sipmbranch = evtree->GetBranch("SiPMSignalData");
    TLeafD  *sipm_signalleaf = (TLeafD*) sipmbranch->GetLeaf("SiPMSignal");
    TLeafD  *eventIDleaf = (TLeafD*) evbranch->GetLeaf("EventID");
    Int_t nevent = evtree->GetEntries();


    Int_t primID, eventID, currentEventID, foundID, SiPMID, xfit, yfit, rootSiPMID;
    //TCanvas *c11 = new TCanvas("с11", "Canvas_1",258,136,1658,787);
    
    Double_t px0prim, px1prim, py0prim, py1prim, theta;
    Double_t xin, xout, yin, yout, zin, zout;
    Double_t primX, primY, primZ;
    Double_t energy, length, energy_phe;

    std::vector<double> primXvec;
    std::vector<double> primZvec;
    std::vector<double> primYvec;
     
    TStrip curstrip;
    std::vector<TStrip> ystripvec;
    std::vector<TStrip> xstripvec;
    
    px0prim = 0;
    px1prim = 0;
    py0prim = 0;
    py1prim = 0;
    theta = -10;
    foundID = -40;
    
    Int_t primentry = 0;

    currentEventID = -10;
    //eventID = 29; 

    for (int entry = 0; entry < nevent; entry++){
        evbranch->GetEntry(entry);
        currentEventID = eventIDleaf->GetValue();
        if (std::find(goodeventid.begin(), goodeventid.end(), currentEventID) != goodeventid.end()){ 
            eventID = currentEventID;
            for (int i = primentry; i<nprim; i++) {    //Looking for this new event in raw MC
                file->cd();
                primbranch->GetEntry(i);
                primID = primeventIDleaf->GetValue();
                // cout << "index, ID: " << i << ", " << primID << "\n"; 

                if (primID == eventID) {
                    foundID = primID;
                    primX = PrimX(primYleaf->GetValue());
                    primY = PrimY(primXleaf->GetValue());
                    primZ = PrimZ(primZleaf->GetValue());

                    if (primZ>=0 && primZ<=DANSSzSize && primX>=0 && primX<=DANSSxSize && primY>=0 && primY<=DANSSySize) {
                        primXvec.push_back(primX);
                        primZvec.push_back(primZ);
                        primYvec.push_back(primY);
                    }
                }

                if (primID > eventID)
                {
                    primentry = i;
                    break;

                }
            }
            //cout<<" cycle passed! Primary: "<<foundID<<", data: "<<eventID << "\n";
            //cout << primXvec.size() <<"\tvecsize\n";

            if (primXvec.size()>0)
            {
                // cout<<"Hit!\t"<<eventID<<"\n";
                xfit = FFit(eventID, primXvec.size(), primXvec.data(), primZvec.data(), &px0prim, &px1prim);
                if (int(xfit)!=0) cerr<<"Problems in primary fit x, eventID:\t"<<eventID<<"\n";

                yfit = FFit(eventID, primYvec.size(), primYvec.data(), primZvec.data(), &py0prim, &py1prim);
                if (int(yfit)!=0) cerr<<"Problems in primary fit y, eventID:\t"<<eventID<<"\n";

                // cout << "primary x: (p0, p1) " << px0prim << ", " << px1prim << "\n";
                // cout << "primary y: (p0, p1) " << py0prim << ", " << py1prim << "\n";
            }

            else {
                cout<<"Empty!\n";
                px0prim = 0;
                px1prim = 0;
                py0prim = 0;
                py1prim = 0;
            }

            evbranch->GetEntry(entry);
            sipmbranch->GetEntry(entry);

            for (int column = 0; column <25; column++){
                for (int row = 0; row< 100; row++){

                    SiPMID = 25*row+column;
                    rootSiPMID = 100*column+row;
                    energy = signalleaf->GetValue(rootSiPMID);
                    //there should be correction for attenuation
                    energy_phe = sipm_signalleaf->GetValue(rootSiPMID);

                    if (energy>0){
                        SiPMID = AlexSiPMNumbering(SiPMID);

                        if ((int(SiPMID)/25)%2==0){ //X   
                            length = Length(SiPMID, px0prim, px1prim, py0prim, py1prim, &xin, &xout, &yin, &yout, &zin, &zout);
                            curstrip.SiPMID = SiPMID;
                            curstrip.energy = energy;
                            curstrip.XY = SiPMXYCoordinate(SiPMID);
                            curstrip.Z = SiPMZCoordinate(SiPMID); 
                            curstrip.status = 1;
                            xstripvec.push_back(curstrip);
                        }
                        else {  //Y
                            length = Length(SiPMID, px0prim, px1prim, py0prim, py1prim, &xin, &xout, &yin, &yout, &zin, &zout);
                            curstrip.SiPMID = SiPMID;
                            curstrip.energy = energy;
                            curstrip.XY = SiPMXYCoordinate(SiPMID);
                            curstrip.Z = SiPMZCoordinate(SiPMID);
                            curstrip.status = 1;
                            ystripvec.push_back(curstrip);  
                        }
                    }
                }
            }  
            cout << "draw hist in " << eventID << " true event " << endl;
            std::stringstream dataname0;
            dataname0 << "mc0e" << eventID << ".json";                           //file event << id
            write(dataname0.str(), xstripvec, ystripvec, 0, 0, 0, 0);          
        }
    }


    return 0;
}

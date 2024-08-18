#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>

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

using std::cout;
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


int main(int argc, char* argv[]) 
{

AlekseevEventReader reader = AlekseevEventReader("/work/danss/danss_134194.root");

Long64_t Nevent = reader.GetEntries(); 
Long64_t id = 55010;                                                            // from 0 to Nevent

    Long64_t global_time = reader.ReadEvent(id).global_time;
    Long64_t veto_hits = reader.ReadEvent(id).veto_hits;
    Double_t veto_energy = reader.ReadEvent(id).veto_energy;

    std::vector<TStrip> strips_xz = reader.ReadEvent(id).strips_xz;
    std::vector<TStrip> strips_yz = reader.ReadEvent(id).strips_yz;

    std::vector<PMTInfo> pmts_xz = reader.ReadEvent(id).pmts_xz;
    std::vector<PMTInfo> pmts_yz = reader.ReadEvent(id).pmts_yz;

    Int_t event_size_x = strips_xz.size();
    Int_t event_size_y = strips_yz.size();
    Int_t event_size = event_size_x + event_size_y;

    if (event_size > 0) {

        cout << "event number: " << id << "; event size: " << event_size << endl;
        Int_t i=0; 
        while (i<event_size_x) {                                      //выкидываем шум (#1)
            if(n_pix(strips_xz[i].energy) < 3) {
                cout << "removing noise" << endl;
                strips_xz.erase(strips_xz.begin()+i);
                event_size_x -=1;
            }
            else {i++;}
        }
        i=0; 
        while (i<event_size_y) {
            if(n_pix(strips_yz[i].energy) < 3) {
                cout << "removing noise" << endl;
                strips_yz.erase(strips_yz.begin()+i);
                event_size_y -=1;
            }
            else {i++;}
        }

        //Time?

        event_size_x = strips_xz.size();
        event_size_y = strips_yz.size();
        event_size = event_size_x + event_size_y;

        cout << "event size after removing noise " << event_size << endl;

        if (event_size_x >= 5 && event_size_y >= 5 && event_size >= 20) {           //выкидываем мелкие события(#3)
            
            for (Int_t i=0; i<event_size_x; i++) {
                cout << i <<" f(x)=" << strips_xz[i].Z << "; x=" <<strips_xz[i].XY << endl;
            }
            for (Int_t i=0; i<event_size_y; i++) {
                cout << i <<" f(y)=" << strips_yz[i].Z << "; y=" <<strips_yz[i].XY << endl;
            }
            
            Double_t k_x0;
            Double_t c_x0;
            MNK(strips_xz, &k_x0, &c_x0);

            Double_t k_y0;
            Double_t c_y0;
            MNK(strips_yz, &k_y0, &c_y0);


            cout << "k_x=" << k_x0 << "; c_x=" << c_x0 << endl;
            cout << "k_y=" << k_y0 << "; c_y=" << c_y0 << endl;


            std::vector<TStrip> strips_eject_xz;
            std::vector<TStrip> strips_eject_yz;

            Int_t j = 0; 
            while (j < event_size_x) {                                               // строим МНК и выкидываем выбивающиеся (#4)
                if(abs(k_x0*strips_xz[j].XY + c_x0 - strips_xz[j].Z) > 6) {          //(надеюсь, что strips_xz[j].XY в сантиметрах)
                    strips_eject_xz.emplace_back(strips_xz[j]);
                    strips_xz.erase(strips_xz.begin()+j);
                    event_size_x -=1;
                }
                else {j++;}
            }
            j = 0; 
            while (j < event_size_y) {                              
                if(abs(k_y0*strips_yz[j].XY + c_y0 - strips_yz[j].Z) > 6) {           
                    strips_eject_yz.emplace_back(strips_yz[j]);
                    strips_yz.erase(strips_yz.begin()+j);
                    event_size_y -=1;
                }
                else{j++;}
            }

            cout << "size of error is " << strips_eject_xz.size()+strips_eject_yz.size() << endl;

            Double_t k_x1;
            Double_t c_x1;
            MNK(strips_xz, &k_x1, &c_x1);

            Double_t k_y1;
            Double_t c_y1;
            MNK(strips_yz, &k_y1, &c_y1);

            for (Int_t j = 0; j < strips_eject_xz.size(); j++) {                               // перестраиваем МНК и возвращаем  (#5)
                if(abs(k_x1*strips_eject_xz[j].XY + c_x1 - strips_eject_xz[j].Z) <= 6) {           
                    strips_xz.emplace_back(strips_eject_xz[j]);
                }
            }

            for (Int_t j = 0; j < strips_eject_yz.size(); j++) {                              
                if(abs(k_x1*strips_eject_yz[j].XY + c_x1 - strips_eject_yz[j].Z) <= 6) {           
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
                cout << "draw hist in " << id << " event with size of " << event_size << endl;
                //-------//
            }

            else {cout << "bad event (adron)" << endl;}
        }

        else {cout << "too few hits" << endl << endl;}

    }

    return 0;
}

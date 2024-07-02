void energyspectrum(){

gROOT->Reset();


TFile *file = new TFile("99.root", "READ");

TTree *evtree = (TTree*) file->Get("PositronTree");
TBranch *evbranch = evtree->GetBranch("SummaryBranch");

TLeafD  *mix_energy_leaf = (TLeafD*) evbranch->GetLeaf("mix_energy");

Int_t nevent;

TH1D *henergy = new TH1D("henergy", "Energy spectrum", 100, 0, 50);


nevent = evbranch->GetEntries();
cout << "Total events: " << nevent << endl;

for(int i = 0; i < nevent; i++)
{
    evbranch->GetEntry(i);
    henergy->Fill(mix_energy_leaf->GetValue());
       
}

TCanvas* c1 = new TCanvas("c1", "Signal data", 30, 20, 750, 550);
c1->cd();

henergy->Draw();

}

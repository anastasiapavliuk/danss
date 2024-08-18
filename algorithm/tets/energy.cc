void energy(){
    
    gROOT->Reset();
    
    TFile *file = new TFile("danss_134328.root", "READ");

    TTree *evtree = (TTree*) file->Get("DanssEvent");

    TBranch *evbranch1 = (TBranch*) evtree->GetBranch("HitE");
    TLeafF *energy_leaf = (TLeafF*) evbranch1->GetLeaf("HitE");

    TBranch *evbranch2 = (TBranch*) evtree->GetBranch("HitType");
    TLeafI *type_leaf = (TLeafI*) evbranch1->GetLeaf("HitType");

    Int_t nevent;

    nevent = evbranch1->GetEntries();

    TH1D *energy = new TH1D("energy", "Energy spectrum", 100, 0, 50);

    cout << "Total events: " << nevent << endl;

    Int_t j = 0;

    for(int i = 0; i < nevent; i++)
    {
        evbranch1->GetEntry(i);
        energy->Fill(energy_leaf->GetValue());
        
    }
    TCanvas* c1 = new TCanvas("c1", "Signal data", 30, 20, 750, 550);
    c1->cd();

    energy->Draw();

}

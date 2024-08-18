void distribution(){

    gROOT->Reset();


    TFile *file = new TFile("99.root", "READ");

    TTree *evtree = (TTree*) file->Get("PositronTree");
    TBranch *evbranch = evtree->GetBranch("SummaryBranch");

    TLeafD  *time_leaf = (TLeafD*) evbranch->GetLeaf("time");

    Int_t nevent;

    TH1D *hits = new TH1D("hits", "Time Distribution", 100, 0, 1000);


    nevent = evbranch->GetEntries();
    cout << "Total events: " << nevent << endl;

    for(int i = 0; i < nevent; i++)
    {
        evbranch->GetEntry(i);
        hits->Fill(time_leaf->GetValue());
        
    }

    TCanvas* c1 = new TCanvas("c1", "Signal data", 30, 20, 750, 550);
    c1->cd();

    hits->Draw();

}

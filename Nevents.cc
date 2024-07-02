void Nevents(TString filename)
{
    TFile *file = new TFile(filename, "READ");

    TTree *evtree = (TTree*) file->Get("PositronTree");
    TBranch *evbranch = evtree->GetBranch("SummaryBranch");

    Int_t nevent;

    nevent = evbranch->GetEntries();
    cout << "." << nevent;
}

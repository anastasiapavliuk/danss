#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include "TF1.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"

using std::cout;
using std::endl;



int main(int argc, char* argv[])
{
    std::string filename = argv[1];
    TFile* file = new TFile(filename.c_str());
    TTree* evtree = (TTree*) file->Get("SummaryTree");
    TBranch* evbranch = evtree->GetBranch("SummaryBranch");
    int nevent;
    nevent = evbranch->GetEntries();
    cout << nevent << std::endl;
    delete file, evtree, evbranch;
    return 0;
}

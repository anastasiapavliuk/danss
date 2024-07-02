#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>

#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TCanvas.h"

using std::cout;
using std::endl;



int main(int argc, char* argv[])
{

    TH1D* htest = new TH1D("htest", "some hist; x;y", 100, 0, 10);     
    delete htest;
    return 0;
}



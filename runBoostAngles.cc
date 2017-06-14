#include <iostream>
#include <string>
#include <sstream>
using namespace std;


#include "rootIncludes.inc"

#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TFrame.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH1.h"

#include "TH2F.h"
#include "TTree.h"
#include "TDirectory.h"
#include "Riostream.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"

#include "BoostAngles.C"


//====================================
int main(int argc, char** argv){

	Double_t nSigma = 3.;
	
  gROOT->ProcessLine(".L BoostAngles.C+");


  BoostAngles();


}
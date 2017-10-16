#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#include "genDefs_data.h"

#include <string>
#include <iostream>

double g_mass_B(double mass)
{
  // here assumed to be a flat distribution.
  // The (unnormalized) background best-fit function should go here.
  // Used only to calculate proportion of sideband-to-signal-region events.
  // If the BG shape changes with pT, this proportion should be
  // given directly as input to the algorithm as a pT function
  return 1.;
}


void polSub(const std::string& filename)
{
  using namespace std;
  // calculate proportion of sideband-to-signal-region background events
  // This is the ONLY necessary information from the invariant-mass fit to
  // the BG subtraction algorithm; if this fraction depends on pT (due to
  // a significant change in *shape* of the BG function) this pT dependence
  // should be determined and given as input to the algorithm

  double normBG_SR = 0.;
  double normBG_tot  = 0.;

  const int n_intgrsteps = 100000;

  double intgrstep = ( mass_max - mass_min ) / double(n_intgrsteps);

  for ( int istep = 0; istep < n_intgrsteps; istep++ ) {

    double xmass = mass_min + istep * intgrstep;

    normBG_tot = normBG_tot + g_mass_B(xmass)*intgrstep;
    if ( xmass > lftSideEnd && xmass < rgtSideStt ) normBG_SR = normBG_SR + g_mass_B(xmass)*intgrstep;
  }

  double normBG_Side = normBG_tot - normBG_SR;
  double NBGtot_div_NBGside = normBG_tot / normBG_Side;

  //////////////////////////////////////////////////////////////////////////////////////////////////////

  gROOT->Reset();

  // input ntuple

  TFile* genFile = TFile::Open(filename.c_str(), "update");

  TTree* genData = (TTree*)genFile->Get("genData");

  // structure of existing ntuple

  double mass;      genData->SetBranchAddress( "mass",         &mass     );

  double wS;
  TBranch* wSbranch = genData->Branch("wS",  &wS,  "wS/D"  );


  int numEvts = int( genData->GetEntries() );
  int n_step = numEvts/50;


  ///////////////////////////////////////////////////////////////
  // loop over events in the input ntuple, to calculate weight
  ///////////////////////////////////////////////////////////////


  cout << endl;
  cout << "Loop to calculate signal weight"<< endl;
  cout << "Reading " << numEvts << " dilepton events"<< endl;
  cout << "-------------------------------------------------------------" << endl;
  cout << "Progress: ";


  for ( int i = 0; i < numEvts; i++ ) {

    if ( i%n_step == 0 ) cout << "X";

    genData->GetEvent( i );


    ////////////////////////////////////////////////////////////////
    // BG SUBTRACTION ALGORITHM ////////////////////////////////////


    if ( mass > lftSideEnd && mass < rgtSideStt ) { wS = 1.; }
    else { wS = 1. - NBGtot_div_NBGside; }


    ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////
    // use weight wS to plot signal distributions
    //          1-wS to plot BG distributions
    ////////////////////////////////////////////////////////////////


    wSbranch->Fill();

  } // end of weight calculations


  cout << endl << endl;


  /////// end

  genFile->Write(nullptr, TObject::kWriteDelete);
  genFile->Close();

}

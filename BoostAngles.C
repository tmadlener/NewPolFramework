#include "rootIncludes.inc"
//#include "calcPol.C"
#include "TH3D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

enum{L,R};
const Int_t kNbSpecies = 3;
enum{UPS1S, UPS2S, UPS3S};

double contamination2Sin1S;
double contamination1Sin2S;
double contamination3Sin2S;
double contamination2Sin3S;

void BoostAngles(Int_t nSigma=3
                      ){

  Char_t fileNameIn[100];
  sprintf(fileNameIn, "tmpFiles/selEvents_data_Ups.root");
  Char_t fileNameInMass[100];
  sprintf(fileNameInMass, "tmpFiles/massFitParameters_Ups.root");
  //==============================
  //read inputs from input file:
  TFile *fIn = new TFile(fileNameIn);
  TFile *fInmass = new TFile(fileNameInMass);
  TLorentzVector *lepP;
  TLorentzVector *lepN;
  TTree *treeIn = (TTree *) fIn->Get("selectedData");
  TTree *treeIn2 = (TTree *) fInmass->Get("massFitParameters");
  
  if(fIn->Get("selectedData")==NULL){
    printf("\n\n\nMissing data.\n\n\n");
  }
  
  if(fInmass->Get("massFitParameters")==NULL){
    printf("\n\n\nMissing mass params.\n\n\n");
  }

  //==============================

  //definition of output variables 
  Char_t fileNameOut[100];
  sprintf(fileNameOut, "tmpFiles/data.root");
  TFile *fOut = new TFile(fileNameOut, "RECREATE");
  gStyle->SetPadRightMargin(0.2);
  TTree *treeOut =  new TTree ("selectedData", "Boosted events");
  TTree *treeOut2 = treeIn2->CloneTree(0);

  //==========================================================
  //reading fit parameters to establish signal mass window
  //as well as the L and R sideband window for the 3D BG histo
  //==========================================================
  fInmass->cd();
  TTree *treeFitPar = (TTree *) gDirectory->Get("massFitParameters");
  TF1 *fUps[kNbSpecies], *fBG = 0;
  fUps[0] = 0, fUps[1] = 0, fUps[2] = 0;
  treeFitPar->SetBranchAddress("fUps1S", &fUps[0]);
  treeFitPar->SetBranchAddress("fUps2S", &fUps[1]);
  treeFitPar->SetBranchAddress("fUps3S", &fUps[2]);
  treeFitPar->SetBranchAddress("fBG", &fBG);
  treeFitPar->LoadTree(0);
  treeFitPar->GetEntry(0);

  Double_t mass[kNbSpecies], sigma[kNbSpecies];
  for(int iState = 0; iState < kNbSpecies; iState++){
    mass[iState] = fUps[iState]->GetParameter(1);
    sigma[iState] = fUps[iState]->GetParameter(2);
  }
  printf("1S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS1S], sigma[UPS1S]);
  printf("2S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS2S], sigma[UPS2S]);
  printf("3S: mass = %1.3f GeV, sigma = %1.3f GeV\n", mass[UPS3S], sigma[UPS3S]);
  Double_t massMin, massMax, massMin1S, massMax1S, massMin2S, massMax2S, massMin3S, massMax3S;
  massMin = mass[UPS1S] - nSigma*sigma[UPS1S];
  massMax = mass[UPS3S] + nSigma*sigma[UPS3S];
  massMin1S = mass[UPS1S] - nSigma*sigma[UPS1S];
  massMax1S = mass[UPS1S] + nSigma*sigma[UPS1S];
  massMin2S = mass[UPS2S] - nSigma*sigma[UPS2S];
  massMax2S = mass[UPS2S] + nSigma*sigma[UPS2S];
  massMin3S = mass[UPS3S] - nSigma*sigma[UPS3S];
  massMax3S = mass[UPS3S] + nSigma*sigma[UPS3S];
  

  cout<<"massMin = "<<massMin<<endl;
  cout<<"massMax = "<<massMax<<endl;
  cout<<"massMin1S = "<<massMin1S<<endl;
  cout<<"massMax1S = "<<massMax1S<<endl;
  cout<<"massMin2S = "<<massMin2S<<endl;
  cout<<"massMax2S = "<<massMax2S<<endl;
  cout<<"massMin3S = "<<massMin3S<<endl;
  cout<<"massMax3S = "<<massMax3S<<endl;

  printf("--> signal mass window: %1.3f < M < %1.3f GeV\n", massMin, massMax);

  //calculate the L and R mass windows:
  Double_t massMinBG[2], massMaxBG[2];
  massMinBG[L] = 8.6;
  massMaxBG[L] = mass[UPS1S] - 3*sigma[UPS1S];
  massMinBG[R] = mass[UPS3S] + 3*sigma[UPS3S];
  massMaxBG[R] = 11.4;
  printf("--> L mass window: %1.3f < M < %1.3f GeV\n", massMinBG[L], massMaxBG[L]);
  printf("--> R mass window: %1.3f < M < %1.3f GeV\n", massMinBG[R], massMaxBG[R]);
  
  Double_t nBGSB = fBG->Integral(massMinBG[L],massMaxBG[L])+fBG->Integral(massMinBG[R],massMaxBG[R]);
  Double_t nBGSR1S = fBG->Integral(massMin1S,massMax1S);
  Double_t nBGSR2S = fBG->Integral(massMin2S,massMax2S);
  Double_t nBGSR3S = fBG->Integral(massMin3S,massMax3S);

  fInmass->Close();
  fIn->cd();
  
  lepP = 0; lepN = 0;
  Double_t Nch, pT, costh_HX, phi_HX, costh_CS, phi_CS, onia_mass, onia_rap, w_Y1S, w_Y2S, w_Y3S;
  treeIn->SetBranchAddress("Nch", &Nch);
  treeOut->Branch("pT", &pT, "pT/D");
  treeOut->Branch("mass", &onia_mass, "mass/D");
  treeOut->Branch("rap", &onia_rap, "rap/D");
  treeOut->Branch("costh_HX", &costh_HX, "costh_HX/D");
  treeOut->Branch("phi_HX", &phi_HX, "phi_HX/D");
  treeOut->Branch("costh_CS", &costh_CS, "costh_CS/D");
  treeOut->Branch("phi_CS", &phi_CS, "phi_CS/D");
  treeOut->Branch("Nch", &Nch, "Nch/D");
  treeOut->Branch("w_Y1S", &w_Y1S, "w_Y1S/D");
  treeOut->Branch("w_Y2S", &w_Y2S, "w_Y2S/D");
  treeOut->Branch("w_Y3S", &w_Y3S, "w_Y3S/D");
  treeIn->SetBranchAddress("lepP", &lepP);
  treeIn->SetBranchAddress("lepN", &lepN);
  TLorentzVector *onia = new TLorentzVector();

//TLorentzVector lepton_DILEP = *lepP;
const double pbeam_ = 7000.; // exact number irrelevant as long as pbeam >> Mprot
const double Mprot_ = 0.9382720;
const double gPI_ = TMath::Pi();
const double Ebeam_ = sqrt( pbeam_*pbeam_ + Mprot_*Mprot_ );
TLorentzVector beam1_LAB_( 0., 0., pbeam_, Ebeam_ );
TLorentzVector beam2_LAB_( 0., 0., -pbeam_, Ebeam_ );

  
  
  for(int iEn = 0; iEn < treeIn->GetEntries(); iEn++){ 
  
    Long64_t iEntry = treeIn->LoadTree(iEn);
    treeIn->GetEntry(iEntry);
    if(iEn % 100000 == 0)
      cout << "entry " << iEntry << " out of " << treeIn->GetEntries() << endl;

    *onia = *(lepP) + *(lepN);
    onia_mass = onia->M();
    pT = onia->Pt();
    onia_rap  = onia->Rapidity();
    
//Calculating boosted angles
////////////////////////////
////////////////////////////
////////////////////////////

TVector3 lab_to_dilep = -onia->BoostVector();

 TLorentzVector beam1_DILEP = beam1_LAB_;
	    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
	    TLorentzVector beam2_DILEP = beam2_LAB_;
	    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

	    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
	    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
	    TVector3 dilep_direction     = onia->Vect().Unit();
	    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


// all polarization frames have the same Y axis = the normal to the plane formed by
	    // the directions of the colliding hadrons:

	    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

	    // flip of y axis with rapidity:

	    if ( onia_rap < 0. ) Yaxis = - Yaxis;

	    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();
	    
	    TLorentzVector lepton_DILEP = *lepP;
	    lepton_DILEP.Boost(lab_to_dilep);

	    // CS frame angles:

	    TVector3 newZaxis = beam1_beam2_bisect;
	    TVector3 newYaxis = Yaxis;
	    TVector3 newXaxis = newYaxis.Cross( newZaxis );

	    TRotation rotation;
	    rotation.SetToIdentity();
	    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame
	    TVector3 lepton_DILEP_rotated = lepton_DILEP.Vect();
	    lepton_DILEP_rotated.Transform(rotation);

	     costh_CS = lepton_DILEP_rotated.CosTheta();
	     phi_CS   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
	    double phith_CS;
	    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
	    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;
	    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


	    // HELICITY frame angles:

	    newZaxis = dilep_direction;
	    newYaxis = Yaxis;
	    newXaxis = newYaxis.Cross( newZaxis );

	    rotation.SetToIdentity();
	    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
	    rotation.Invert();
	    lepton_DILEP_rotated = lepton_DILEP.Vect();
	    lepton_DILEP_rotated.Transform(rotation);

	     costh_HX = lepton_DILEP_rotated.CosTheta();
	     phi_HX   = lepton_DILEP_rotated.Phi() * 180. / gPI_;
	    double phith_HX;
	    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
	    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;
	    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;
////////////////////////////
////////////////////////////
////////////////////////////
  


  if(onia_mass > massMin1S && onia_mass < massMax1S) {w_Y1S = 1; w_Y2S = 0; w_Y3S = 0;}
  if(onia_mass > massMin2S && onia_mass < massMax2S) {w_Y1S = 0; w_Y2S = 1; w_Y3S = 0;}
  if(onia_mass > massMin3S && onia_mass < massMax3S) {w_Y1S = 0; w_Y2S = 0; w_Y3S = 1;}
  else {w_Y1S = 1 - (nBGSB+nBGSR1S)/nBGSB; w_Y2S = 1 - (nBGSB+nBGSR2S)/nBGSB; w_Y3S = 1 - (nBGSB+nBGSR3S)/nBGSB;}
  

 //   if(onia_mass > massMin && onia_mass < massMax){
      treeOut->Fill(); //stores TLorenzVectors of the two muons

      // //store now the cosTheta and phi distributions of the signal window:
      // calcPol(*lepP, *lepN);

   // }
  }



  //write the output
  fOut->cd();
  treeOut->Write();
  treeOut2->Write();
  fOut->Close();
  fIn->Close();

}

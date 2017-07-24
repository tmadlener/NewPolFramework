#define data_or_ref 0 // 0 = reference, 1 = data

#include "Riostream.h"
#include "TSystem.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"

#if data_or_ref == 1
#include "genDefs_data.h"
#else
#include "genDefs_ref.h"
#endif

const double pbeam = 7000;  // (actually there is no dependence of polarization-frame definitions on beam energy)
const double Mprot = 0.9382720;
const double Mlepton = 0.10566;  // (muon)
double gPI = TMath::Pi();
double Ebeam = sqrt( pbeam*pbeam + Mprot*Mprot );

TLorentzVector beam1_LAB( 0., 0., pbeam, Ebeam );
TLorentzVector beam2_LAB( 0., 0., -pbeam, Ebeam );


void polGen(){


  gROOT->Reset();

  delete gRandom;
  gRandom = new TRandom3(0);

  TF1* pT_distr = new TF1("pT_distr",func_pT_gen,pTdilepton_min,pTdilepton_max,0);
  TF1* pT_distr_BG = new TF1("pT_distr_BG",func_pT_gen_BG,pTdilepton_min,pTdilepton_max,0);

  TF1* Nch_distr = new TF1("Nch_distr",func_Nch_gen,Nch_min,Nch_max,0);
  TF1* Nch_distr_BG = new TF1("Nch_distr_BG",func_Nch_gen_BG,Nch_min,Nch_max,0);

  TF1* rap_distr = new TF1("rap_distr",func_rap_gen,rapdilepton_min,rapdilepton_max,0);

  TFile* hfileout = new TFile(outfile, "RECREATE", "genData");
  TTree* genData = new TTree("genData","genData");

  // Structure of output ntuple

  TLorentzVector* lepP = new TLorentzVector(0.,0.,0.,0.);
  genData->Branch("lepP","TLorentzVector",&lepP);
  TLorentzVector* lepN = new TLorentzVector(0.,0.,0.,0.);
  genData->Branch("lepN","TLorentzVector",&lepN);

  double costh_CS;  genData->Branch("costh_CS",     &costh_CS,     "costh_CS/D");
  double phi_CS;    genData->Branch("phi_CS",       &phi_CS,       "phi_CS/D"  );
  double phith_CS;  genData->Branch("phith_CS",     &phith_CS,     "phith_CS/D");

  double costh_HX;  genData->Branch("costh_HX",     &costh_HX,     "costh_HX/D");
  double phi_HX;    genData->Branch("phi_HX",       &phi_HX,       "phi_HX/D"  );
  double phith_HX;  genData->Branch("phith_HX",     &phith_HX,     "phith_HX/D");

  double costh_PX;  genData->Branch("costh_PX",     &costh_PX,     "costh_PX/D");
  double phi_PX;    genData->Branch("phi_PX",       &phi_PX,       "phi_PX/D"  );
  double phith_PX;  genData->Branch("phith_PX",     &phith_PX,     "phith_PX/D");

  double cosalpha;  genData->Branch("cosalpha",     &cosalpha,     "cosalpha/D");

  double pT;        genData->Branch("pT",           &pT,           "pT/D"  );
  double rap;       genData->Branch("rap",          &rap,          "rap/D" );
  double mass;      genData->Branch("mass",         &mass,         "mass/D");
  double Nch;       genData->Branch("Nch",          &Nch,          "Nch/D" );

  double deltaHXCS; genData->Branch("deltaHXCS",    &deltaHXCS,    "deltaHXCS/D");

  int isBG;         genData->Branch("isBG",         &isBG,         "isBG/I");


  const int n_step = n_events/50;


  cout << endl;
  cout << "Generating " << n_events << " dilepton events"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";



  /////////////////// CYCLE OF EVENTS ////////////////////////
  for(int i_event = 1; i_event <= n_events; i_event++){

    // generation of dilepton in the pp event in the pp CM

    // mass and pT

    isBG = 0;

    if ( gRandom->Uniform() < f_BG ) { isBG = 1;
      mass = gRandom->Uniform(mass_min, mass_max);
      pT = pT_distr_BG->GetRandom();
      Nch = Nch_distr_BG->GetRandom(); }
    else { do { mass = gRandom->Gaus(mass_signal_peak, mass_signal_sigma); }
      while ( mass < mass_min || mass > mass_max );
      pT = pT_distr->GetRandom();
      Nch = Nch_distr->GetRandom();  }

    // pL:

    double rap_sign = gRandom->Uniform(-1., 1.); rap_sign /= fabs(rap_sign);
    rap = rap_distr->GetRandom() * rap_sign;
    double mT = sqrt( mass*mass + pT*pT );
    double pL1 = 0.5 *mT * exp(rap);
    double pL2 = - 0.5 *mT * exp(-rap);
    double pL = pL1 + pL2;

    // Phi:

    double Phi = 2. * gPI * gRandom->Uniform(1.);

    // 4-vector:

    TLorentzVector dilepton;
    dilepton.SetXYZM( pT * cos(Phi) , pT * sin(Phi), pL, mass );


    // generation of polarization (generic reference frame)

    double lambda_theta    = lambda_theta_sig(pT, Nch);

    double lambda_phi      = lambda_phi_sig(pT, Nch);
    double lambda_thetaphi = lambda_thetaphi_sig(pT, Nch);
    bool HX_is_natural = HX_is_natural_sig;
    bool PX_is_natural = PX_is_natural_sig;

    if ( isBG ) {
      lambda_theta    = lambda_theta_bkg(pT, Nch);
      lambda_phi      = lambda_phi_bkg(pT, Nch);
      lambda_thetaphi = lambda_thetaphi_bkg(pT, Nch);
      HX_is_natural = HX_is_natural_bkg;
      PX_is_natural = PX_is_natural_bkg;
    }


    double costhphidistr_max = 1. + fabs(lambda_phi) + fabs(lambda_thetaphi);
    double costhphidistr_rnd;
    double costhphidistr;
    double costh_gen;
    double sinth_gen;
    double phi_gen;

    if ( lambda_theta > 0. ) costhphidistr_max += lambda_theta;

    do { costh_gen = -1. + 2. * gRandom->Uniform(1.);
      phi_gen   = 2. * gPI * gRandom->Uniform(1.);
      sinth_gen = sqrt( 1. - costh_gen*costh_gen );
      costhphidistr_rnd = costhphidistr_max * gRandom->Uniform(1.);
      costhphidistr = 1. + lambda_theta    * costh_gen*costh_gen
        + lambda_phi      * sinth_gen*sinth_gen * cos(2.*phi_gen)
        + lambda_thetaphi * 2.* sinth_gen*costh_gen * cos(phi_gen);
    } while ( costhphidistr_rnd > costhphidistr );


    // lepton momentum in the dilepton rest frame:

    double p_lepton_DILEP = sqrt( 0.25*mass*mass - Mlepton*Mlepton );

    TLorentzVector lepton_DILEP;

    lepton_DILEP.SetXYZM( p_lepton_DILEP * sinth_gen * cos(phi_gen),
                          p_lepton_DILEP * sinth_gen * sin(phi_gen),
                          p_lepton_DILEP * costh_gen,
                          Mlepton );


    // reference directions to calculate angles:

    TVector3 lab_to_dilep = -dilepton.BoostVector();

    TLorentzVector beam1_DILEP = beam1_LAB;
    beam1_DILEP.Boost(lab_to_dilep);         // beam1 in the dilepton rest frame
    TLorentzVector beam2_DILEP = beam2_LAB;
    beam2_DILEP.Boost(lab_to_dilep);         // beam2 in the dilepton rest frame

    TVector3 beam1_direction     = beam1_DILEP.Vect().Unit();
    TVector3 beam2_direction     = beam2_DILEP.Vect().Unit();
    TVector3 dilep_direction     = dilepton.Vect().Unit();
    TVector3 beam1_beam2_bisect  = ( beam1_direction - beam2_direction ).Unit();


    deltaHXCS = dilep_direction.Angle(beam1_beam2_bisect) * 180./gPI;

    // all polarization frames have the same Y axis = the normal to the plane formed by
    // the directions of the colliding hadrons

    TVector3 Yaxis = ( beam1_direction.Cross( beam2_direction ) ).Unit();

    // flip of y axis with rapidity

    if ( rap < 0 ) Yaxis = - Yaxis;

    TVector3 perpendicular_to_beam = ( beam1_beam2_bisect.Cross( Yaxis ) ).Unit();


    // step 1: transform (rotation) lepton momentum components from generation frame
    // to the frame with x,y,z axes as in the laboratory

    TVector3 oldZaxis = beam1_beam2_bisect;
    if ( HX_is_natural ) oldZaxis = dilep_direction;
    if ( PX_is_natural ) oldZaxis = perpendicular_to_beam;

    TVector3 oldYaxis = Yaxis;
    TVector3 oldXaxis = oldYaxis.Cross(oldZaxis);

    TRotation rotation;
    rotation.RotateAxes(oldXaxis, oldYaxis, oldZaxis);
    // transforms coordinates from the "old" frame to the "xyz" frame

    TLorentzVector lepton_DILEP_xyz = lepton_DILEP;

    lepton_DILEP_xyz.Transform(rotation);
    // lepton_DILEP_xyz is the lepton in the dilepton rest frame
    // wrt to the lab axes

    // lepton 4-vectors in the LAB frame:

    TVector3 dilep_to_lab = dilepton.BoostVector();

    *lepP = lepton_DILEP_xyz;
    lepP->Boost(dilep_to_lab);
    lepN->SetPxPyPzE(-lepton_DILEP_xyz.Px(),-lepton_DILEP_xyz.Py(),-lepton_DILEP_xyz.Pz(),lepton_DILEP_xyz.E());
    lepN->Boost(dilep_to_lab);


    /////////////////////////////////////////////////////////////////////
    // CS frame


    TVector3 newZaxis = beam1_beam2_bisect;
    TVector3 newYaxis = Yaxis;
    TVector3 newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();  // transforms coordinates from the "xyz" frame to the new frame

    TVector3 lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    costh_CS = lepton_DILEP_rotated.CosTheta();

    phi_CS = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_CS < 0. ) phith_CS = phi_CS - 135.;
    if ( costh_CS > 0. ) phith_CS = phi_CS - 45.;

    if ( phith_CS < -180. ) phith_CS = 360. + phith_CS;


    /////////////////////////////////////////////////////////////////////
    // HELICITY frame

    newZaxis = dilep_direction;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    costh_HX = lepton_DILEP_rotated.CosTheta();

    phi_HX = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_HX < 0. ) phith_HX = phi_HX - 135.;
    if ( costh_HX > 0. ) phith_HX = phi_HX - 45.;

    if ( phith_HX < -180. ) phith_HX = 360. + phith_HX;


    /////////////////////////////////////////////////////////////////////
    // PERPENDICULAR HELICITY frame

    newZaxis = perpendicular_to_beam;
    newYaxis = Yaxis;
    newXaxis = newYaxis.Cross( newZaxis );

    rotation.SetToIdentity();
    rotation.RotateAxes( newXaxis, newYaxis, newZaxis );
    rotation.Invert();

    lepton_DILEP_rotated = lepton_DILEP_xyz.Vect();

    lepton_DILEP_rotated.Transform(rotation);

    costh_PX = lepton_DILEP_rotated.CosTheta();

    phi_PX = lepton_DILEP_rotated.Phi() * 180. / gPI;

    if ( costh_PX < 0. ) phith_PX = phi_PX - 135.;
    if ( costh_PX > 0. ) phith_PX = phi_PX - 45.;

    if ( phith_PX < -180. ) phith_PX = 360. + phith_PX;


    /////////////////////////////////////////////////////////////////////
    // invariant polarization angle

    cosalpha = sqrt( 1. - pow(costh_PX, 2.) ) * sin( lepton_DILEP_rotated.Phi() );



    //  filling of the ntuple:

    bool acceptedEvent = lepN->Pt() > leptoncut && lepP->Pt() > leptoncut;

    if( acceptedEvent ) genData->Fill();

    if (i_event%n_step == 0) cout << "X";

  } // end of external loop (generated events)

  cout << endl << endl;

  hfileout->Write();
  hfileout->Close();

} // end of main

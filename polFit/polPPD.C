#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"

#include<iostream>
#include<string>
#include<sstream>

const bool is_analysis_vs_pT = true;

const double pTmin_analysis_not_vs_pT = 10.;   // set min and max pT in analysis done vs Nch
const double pTmax_analysis_not_vs_pT = 100.;  //

string nameOfX(is_analysis_vs_pT?"chicPt":"Nch");

// x dependence: grade of polynomial parametrization
const int npar_th = 1;
const int npar_ph = 1;
const int npar_tp = 1;

const int pTreweightingIter = 3; // number of iterations for reweighting of the reference pT distribution
// for measurements vs pT: set to at least 2; ineffective for analysis of relative polarization vs Nch:
const int n_iterations = is_analysis_vs_pT?pTreweightingIter:1;

// for relative Nch measurement, no pT reweighting will be applied to the "reference"
const bool apply_pT_weight_to_ref = is_analysis_vs_pT?true:false;

const double x_unit = 1.; // a "unit" for x; pT -> GeV
const double n_ev_min = 100;  // in the calculation of the extremes in x of the analysis
                              // a minimum event-density of n_ev_min/x_unit will be required:
                              // margins of the distribution are cut out until they have more than n_ev_min events per x unit
const double xmin_imposed =     0.*x_unit;
const double xmax_imposed =  1000.*x_unit;
// min and max x values in the data sample will be determined automatically,
// but the range considered can be further restricted changing these definitions

const int n_xbins = 4;
// preliminary x binning to determine the starting parameter values and sigmas
// This number is not related to npar_th etc, but it must be >= 3
const double kfactor = 2.0; // defines the increase in size of bin i+1 wrt bin i

const int ncellcosth = 12; // for preliminary fits
const int ncellphi = 15;

///// PARAMETER SCANNING STRATEGY /////////////////////////////////////////////

const int n_scan_points = 10000;

const int n_scan_phases = 4;

double scan_progress_vs_phase[n_scan_phases]
= { .50, .80, .90, 1.0 }; // must be increasing sequence;
// the last phase completes (100%) the scan

// How much broader should the Gaussian scan sigma be
// wrt the last temporary calculation of the sigma of a parameter PPD?
double scan_sigma_margin_factor[n_scan_phases]
= { 1.6, 1.4, 1.3, 1.2 };

// physical boundaries for parameter scanning:
// outside them the PPD will not be calculated (= zero)

const double Ath_min = 0.0;
const double Ath_max = 1.0;

const double Aph_min = -1.0;
const double Aph_max = 1.0;

const double Atp_min = -2./3.*sqrt(2.);
const double Atp_max = -Atp_min;

/////////////////////////////////////////////////////////////////////////////////////
// polarization of the reference, vs x:
// IF = 0, the DIFFERENCE wrt the reference polarization is measured.
// Actually this is a first-order approximation, valid when
// A) the polarization difference is not big
// B) OR lambdath_ref is not big
// A better approximation is
// DeltAthambdath = lambdath_meas(value obained from the fit) *
//                * [ 1 + lambdath_ref * (1 + 3/5 * lambdath_ref) / (3 + lambdath_ref)
// Note: it is not possible to determine the correction without an estimate of lambdath_ref
// (its sign even determines the sign of the correction!)
// (similar lambdaphi and lambdathetaphi expressions can be calculated but given
// the smallness of the measured values the 1st order approximation is certainly sufficient)

double lambdath_ref(double x)  // x could here represent pT or any other variable wrt which lambda of reference is known to change
{                              // Nch analysis: set to zero for relative measurement, or to the known Nch-integrated polarization
  return 0;                    // for absolute measurement
}
double lambdaph_ref(double x)
{
  return 0;
}
double lambdatp_ref(double x)
{
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// function fitting the pT-distribution ratio data/reference
// It is necessary to determine the pT-distribution ratio data/reference WITH EQUAL POLARIZATION
// HYPOTHESIS. The the pT distribution of data and reference are different beacause of
// 1) production mechanism (difference at "generation") -> this difference must be equalized
// 2) different polarization, determining different acceptances and different reconstructed distributions
//    -> this effect MUST NOT be equalized out
// To equalize 1 and not 2, pT distributions must fist re-determined AS IF both signal and reference
// had the same polarization
// switch off the data/ref-ratio polarization weight if it is known that
// reference and data have by definition the same pT distribution at "generation" (as in Nch-dependent study)

double func_pT_weight(double* x, double* par)
{
  return par[0] * pow( x[0], par[1]  );
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// parametrization of polarization dependence on x: polynomial with a maximum of 3 support points
// this is used only in the preliminary fit, where the two "IFs" do not cause a problem of speed

double A_vs_x( double* x, double* par )
{                // this function has 1 + 3 + npar_th|npar_ph|npar_tp parameters
  int n_polParams = par[0]; // first parameter indicates the type of polinomial interpolation

  double x1, x2, x3, A1, A2, A3;

  x1 = par[1]; // support points: always 3
  x2 = par[2];
  x3 = par[3];

  A1 = par[4];
  A2 = A1; A3 = A1;
  // constant dependence obtained fixing A[3] = A[2] = A[1]
  // linear dependence obtained imposing A3 = A1 + (x3-x1) * (A2-A1) / (x2-x1)
  if ( n_polParams == 2 ) { A2 = par[5]; A3 = A1 + (x3-x1) * (A2-A1) / (x2-x1); }
  if ( n_polParams == 3 ) { A2 = par[5]; A3 = par[6]; }

  double xval = x[0];

  return   A1*(xval-x2)*(xval-x3)/((x1-x2)*(x1-x3))
    + A2*(xval-x1)*(xval-x3)/((x2-x1)*(x2-x3))
    + A3*(xval-x1)*(xval-x2)/((x3-x1)*(x3-x2));
}

const double gPI = TMath::Pi();

double func_pol_A(double* x, double* par)  // function to fit costh-phi distribution to obtain A's
{
  double costh = x[0];
  double phi = gPI/180.*x[1];
  double Ath = par[1];
  double Aph = par[2];
  double Atp = par[3];
  double sinth2 = 1. - costh*costh;
  double sin2th = 2.*costh*sqrt(sinth2);
  return par[0] * 3./4. * ( 1. + Ath + (1.-3.*Ath) * costh*costh + Aph * sinth2*cos(2*phi) + Atp * sin2th*cos(phi) );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void polPPD(){

  gROOT->Reset();
  gErrorIgnoreLevel = kWarning;

  TStopwatch* timer = new TStopwatch;
  timer->Start();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// open data and reference ntuples ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // TFile* dataFile = TFile::Open("~/work/PhD/data/Upsilon2011_Nch/tests_newFW/Upsilon2011_pt10_rap1p2_sampSplitRand_Seed77_gt/tmpFiles/data.root");
  // TTree* dataSample = (TTree*)dataFile->Get("selectedData");

  TFile* dataFile = TFile::Open("~/cernbox/PhD/data/chic_steppingstone_tests/PR/chic_tuple_PR2.root");
  TTree* dataSample = (TTree*)dataFile->Get("chic_tuple");

  double costh;  dataSample->SetBranchAddress( "costh_HX",     &costh );  // chose frame here
  double phi;    dataSample->SetBranchAddress( "phi_HX",       &phi   );
  double wS;     dataSample->SetBranchAddress( "wChic2",           &wS    );  // BG subtraction weight
  double pT;     dataSample->SetBranchAddress( "chicPt",           &pT    );
  double Nch;    // dataSample->SetBranchAddress( "Nch",          &Nch   );

  // TFile* refFile = TFile::Open("~/work/PhD/data/Upsilon2011_Nch/tests_newFW/Upsilon2011_pt10_rap1p2_sampSplitRand_Seed77_lt/tmpFiles/data.root");
  // TTree* refSample = (TTree*)refFile->Get("selectedData");

  TFile* refFile = TFile::Open("~/cernbox/PhD/data/chic_steppingstone_tests/PR/chic_tuple_PR.root");
  TTree* refSample = (TTree*)refFile->Get("chic_tuple");

  double costh_R;  refSample->SetBranchAddress( "costh_HX",    &costh_R );
  double phi_R;    refSample->SetBranchAddress( "phi_HX",      &phi_R   );
  double wS_R;     refSample->SetBranchAddress( "wChic1",          &wS_R    );
  double pT_R;     refSample->SetBranchAddress( "chicPt",          &pT_R    );
  double Nch_R;     // refSample->SetBranchAddress( "Nch",         &Nch_R   );

  // determine x range and binning

  const int n_cells_per_unit = 10;
  int nxcells = int((xmax_imposed-xmin_imposed)/x_unit)*n_cells_per_unit;

  stringstream command_option;
  command_option << nameOfX << ">>h_x_data(" << nxcells << "," << xmin_imposed << "," << xmax_imposed << ")";

  dataSample->Draw(command_option.str().c_str(),"wChic2","goff");
  TH1D *h_x_data = (TH1D*)gDirectory->Get("h_x_data");

  command_option.str("");
  command_option << nameOfX << ">>h_x_ref(" << nxcells << "," << xmin_imposed << "," << xmax_imposed << ")";

  refSample->Draw(command_option.str().c_str(),"wChic1","goff");
  TH1D *h_x_ref = (TH1D*)gDirectory->Get("h_x_ref");

  // determine total x range as the interesection of data and ref ranges

  double xmin = 1.e20; double xmax =  -1.e20;

  double n_ev_min_per_cell = n_ev_min/double(n_cells_per_unit);

  for ( int i_cell = 1; i_cell <= nxcells; i_cell++ ) {

    double cellcontent_left  = min( h_x_data->GetBinContent(i_cell), h_x_ref->GetBinContent(i_cell) );
    double cellcontent_right = min( h_x_data->GetBinContent(nxcells+1-i_cell), h_x_ref->GetBinContent(nxcells+1-i_cell) );

    double xval_left  = h_x_data->GetXaxis()->GetBinLowEdge(i_cell);
    double xval_right = h_x_data->GetXaxis()->GetBinUpEdge(nxcells+1-i_cell);

    if ( cellcontent_left  > n_ev_min_per_cell && xval_left  < xmin ) { xmin = xval_left;  } // require > n_ev_min events per unit
    if ( cellcontent_right > n_ev_min_per_cell && xval_right > xmax ) { xmax = xval_right; } // in the etreme cells
  }

  cout << endl;
  cout << "Analysis range:  " << xmin << "  <  x  <  " << xmax << endl;

  // calculate preliminary x binning

  double* xbinsize = new double[n_xbins];

  double* xbin = new double[n_xbins+1]; xbin[0] = xmin;

  double xbinsize_norm = 0.;
  for ( int i_bin = 0; i_bin < n_xbins; i_bin++ ) {
    xbinsize[i_bin] = pow( kfactor, double(i_bin) );
    xbinsize_norm += xbinsize[i_bin];
  }
  for ( int i_bin = 0; i_bin < n_xbins; i_bin++ ) {
    xbin[i_bin+1] = xbin[i_bin] + (xmax-xmin) * xbinsize[i_bin] / xbinsize_norm;
  }

  // average x value for each bin

  double* x_i      = new double[n_xbins];   for ( int i = 0; i < n_xbins; i++ ) { x_i[i] = 0.; }
  double* x_i_norm = new double[n_xbins];   for ( int i = 0; i < n_xbins; i++ ) { x_i_norm[i] = 0.; }

  // loop over data histogram cells to calculate centres of x bins

  for ( int i_cell = 1; i_cell <= nxcells; i_cell++ ) {

    double cellcontent = h_x_data->GetBinContent(i_cell);

    double xval = h_x_data->GetXaxis()->GetBinCenter(i_cell);

    for ( int i_bin = 0; i_bin < n_xbins; i_bin++ ) {
      if ( xval > xbin[i_bin] && xval < xbin[i_bin+1] ) {
        x_i[i_bin] += xval * cellcontent;
        x_i_norm[i_bin] += cellcontent;
      }
    }
  }

  for ( int i = 0; i < n_xbins; i++ ) {
    x_i[i] /= x_i_norm[i];
  }

  cout << endl;
  cout << "Binning for preliminary calculations"<< endl;
  cout << endl;
  for ( int i = 0; i < n_xbins; i++ ) {
    cout << "bin " << i << ":  " << xbin[i] << " < x < " << xbin[i+1]
         << ";     <x> = " << x_i[i] << ";     data signal = " << int(x_i_norm[i]) << endl;
  }

  delete h_x_data, h_x_ref;

  int numDataEvtsTot = dataSample->GetEntries();
  int numRefEvtsTot  = refSample->GetEntries();

  command_option.str("");
  command_option << nameOfX << ">" << xmin << "&&"<< nameOfX << "<" << xmax;

  int numDataEvts = dataSample->GetEntries(command_option.str().c_str());  // only the analysis range
  int numRefEvts  = refSample->GetEntries(command_option.str().c_str());

  // Fill arrays of events for data and reference

  double* x_j                 = new double[numDataEvts];
  double* pT_j                = new double[numDataEvts];
  double* costheta_j          = new double[numDataEvts];
  double* phi_j               = new double[numDataEvts];
  double* costheta2_j         = new double[numDataEvts];
  double* sintheta2cos2phi_j  = new double[numDataEvts];
  double* sin2thetacosphi_j   = new double[numDataEvts];
  double* wS_j                = new double[numDataEvts];

  double* x_jR                 = new double[numRefEvts];
  double* pT_jR                = new double[numRefEvts];
  double* costheta_jR          = new double[numRefEvts];
  double* phi_jR               = new double[numRefEvts];
  double* costheta2_jR         = new double[numRefEvts];
  double* sintheta2cos2phi_jR  = new double[numRefEvts];
  double* sin2thetacosphi_jR   = new double[numRefEvts];
  double* wS_jR                = new double[numRefEvts];


  // pT histogram to afterwards calculate pT-weight to be applied to the reference sample
  //////////////////////////////////////////////////
  double pTmin = xmin;
  double pTmax = xmax;

  if ( !is_analysis_vs_pT ) { pTmin = pTmin_analysis_not_vs_pT; pTmax = pTmax_analysis_not_vs_pT; }

  int ncellspT = pTmax-pTmin;
  //////////////////////////////////////////////////
  const double pTstep = (pTmax-pTmin)/double(ncellspT);
  TH1D* h_pT_data  = new TH1D("", "", ncellspT, pTmin, pTmax );
  TH1D* h_pT_ref   = new TH1D("", "", ncellspT, pTmin, pTmax );
  h_pT_ref->Sumw2();
  h_pT_data->Sumw2();
  // for now only the reference histogram will be filled:
  // the data histogram must be filled reweighted inside the iteration loop


  // costh_phi histograms for calculation of preliminary lambda values
  TH2D** h_costhphi_data  = new TH2D*[n_xbins];
  TH2D** h_costhphi_ref   = new TH2D*[n_xbins];
  TH2D** h_costhphi_ratio = new TH2D*[n_xbins];

  for (int i = 0; i < n_xbins; i++) {
    h_costhphi_data[i]  = new TH2D("", "", ncellcosth, -1., 1., ncellphi, -180., 180. );
    h_costhphi_ref[i]   = new TH2D("", "", ncellcosth, -1., 1., ncellphi, -180., 180. );
    h_costhphi_data[i]->Sumw2();
    h_costhphi_ref[i]->Sumw2();
  }
  // for now only the data histograms will be filled:
  // the reference histograms must be filled reweighted inside the iteration loop

  int n_step = numDataEvtsTot/50;

  cout << endl;
  cout << "Filling arrays of data events"<< endl;
  cout << "Reading " << numDataEvtsTot << " dilepton events"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";

  double n_dataSignal = 0.;
  int k = 0;

  for ( int j = 0; j < numDataEvtsTot; j++ ) { //loop over all ntuple events

    if ( (j+1)%n_step == 0 ) cout << "X";

    dataSample->GetEvent( j );

    double xx = is_analysis_vs_pT?pT:Nch;

    if ( xx > xmin && xx < xmax && pT > pTmin && pT < pTmax ) {

      x_j[k] = xx;
      pT_j[k] = pT;
      costheta_j[k] = costh;
      phi_j[k] = phi;
      double costheta2 = costh*costh;
      costheta2_j[k] = costheta2;
      sintheta2cos2phi_j[k] = (1. - costheta2) * cos(gPI/90.*phi);
      sin2thetacosphi_j[k] = 2.*costh*sqrt(1.-costheta2) * cos(gPI/180.*phi);
      wS_j[k] = wS;

      n_dataSignal += wS;

      for (int i = 0; i < n_xbins; i++) {

        if ( xx > xbin[i] && xx < xbin[i+1] ) { h_costhphi_data[i]->Fill(costh, phi, wS); }
      }

      k++;
    }

  } // end of data loop

  // normalize
  for (int i = 0; i < n_xbins; i++) {
    h_costhphi_data[i]->Scale( 1./h_costhphi_data[i]->Integral() );
  }

  cout << endl;
  cout << "Signal events in data sample: " << int(n_dataSignal) << endl;


  n_step = numRefEvtsTot/50;

  cout << endl;
  cout << "Filling arrays of reference events"<< endl;
  cout << "Reading " << numRefEvtsTot << " dilepton events"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";

  double n_refSignal = 0.;
  k = 0;

  for ( int j = 0; j < numRefEvtsTot; j++ ) { //loop over ALL ntuple events

    if ( (j+1)%n_step == 0 ) cout << "X";

    refSample->GetEvent( j );

    double xx = is_analysis_vs_pT?pT_R:Nch_R;

    if ( xx > xmin && xx < xmax  && pT_R > pTmin && pT_R < pTmax ) {

      x_jR[k] = xx;
      pT_jR[k] = pT_R;
      costheta_jR[k] = costh_R;
      phi_jR[k] = phi_R;
      double costheta2 = costh_R*costh_R;
      costheta2_jR[k] = costheta2;
      sintheta2cos2phi_jR[k] = (1. - costheta2) * cos(gPI/90.*phi_R);
      sin2thetacosphi_jR[k] = 2.*costh_R*sqrt(1.-costheta2) * cos(gPI/180.*phi_R);
      wS_jR[k] = wS_R;

      n_refSignal += wS_R;

      double xpol = pT_R; // variable wrt which the reference polarization is known to change
      double wUnpol = ( 3.+lambdath_ref(xpol) ) /3. /(1 + lambdath_ref(xpol) * costheta2
                                                      + lambdaph_ref(xpol) * sintheta2cos2phi_jR[k]
                                                      + lambdatp_ref(xpol) * sin2thetacosphi_jR[k] );
      // weight to eliminate the known polarization of the reference

      h_pT_ref->Fill(pT_R, wS_R * wUnpol); // this histogram is not changed by iterations

      k++;
    }
  } // end of reference loop

  // normalize
  h_pT_ref->Scale( 1./ h_pT_ref->Integral() );

  cout << endl;
  cout << "Signal events in reference sample: " << int(n_refSignal) << endl;


  // x abscissas (npar_th, etc.) of the final results, also used to parametrize
  // lambda(pT) functions for iterative pT reweighting (analysis with x = pT)

  // the choice is obviously arbitrary: the fitting curve will be the same.
  // Chosen: centres of 3 equidistant x bins
  // will use only the first one if npar = 1, only the first two if npar = 2

  double* x   = new double[3];
  double* Ath = new double[3]; // fits done using A parameters
  double* Aph = new double[3];
  double* Atp = new double[3];

  double* dAth = new double[3];
  double* dAph = new double[3];
  double* dAtp = new double[3];


  x[0] = x_i[0];
  x[2] = x_i[n_xbins-1];
  x[1] = 0.5 * ( x[0] + x[2] );


  for ( int i = 0; i < 3; i++ ) {

    Ath[i] = 1./3.; // longitudinal fraction = 1/3 in the unpolarized case
    Aph[i] = 0.;
    Atp[i] = 0.;
  }


  // parameters of the fine-x-binned fit
  // abscissas were defined above as:   double* x_i  = new double[n_xbins];
  double* Ath_i    = new double[n_xbins];
  double* Aph_i    = new double[n_xbins];
  double* Atp_i    = new double[n_xbins];
  double* dAth_i   = new double[n_xbins];
  double* dAph_i   = new double[n_xbins];
  double* dAtp_i   = new double[n_xbins];


  for ( int i = 0; i < n_xbins; i++ ) {

    // starting polarization parameter values for first fit

    Ath_i[i] = 1./3.; // longitudinal fraction = 1/3 in the unpolarized case
    Aph_i[i] = 0.;
    Atp_i[i] = 0.;
  }

  double* ex = new double[20];
  memset(ex, 0., 20*sizeof(double) );

  // polarization fitting functions

  TF2* polFit_A = new TF2("polFit_A",func_pol_A,-1.,1.,-180.,180.,4); // main fitting function

  TF1* Ath_vs_x = new TF1("Ath_vs_x", A_vs_x, xmin, xmax, 4 + npar_th );
  Ath_vs_x->FixParameter(0, npar_th);
  Ath_vs_x->FixParameter(1, x[0]);
  Ath_vs_x->FixParameter(2, x[1]);
  Ath_vs_x->FixParameter(3, x[2]);
  Ath_vs_x->SetParameter(4, Ath[0]);
  if ( npar_th > 1 ) Ath_vs_x->SetParameter(5, Ath[1]);
  if ( npar_th > 2 ) Ath_vs_x->SetParameter(6, Ath[2]);

  TF1* Aph_vs_x = new TF1("Aph_vs_x", A_vs_x, xmin, xmax, 4 + npar_ph );
  Aph_vs_x->FixParameter(0, npar_ph);
  Aph_vs_x->FixParameter(1, x[0]);
  Aph_vs_x->FixParameter(2, x[1]);
  Aph_vs_x->FixParameter(3, x[2]);
  Aph_vs_x->SetParameter(4, Aph[0]);
  if ( npar_ph > 1 ) Aph_vs_x->SetParameter(5, Aph[1]);
  if ( npar_ph > 2 ) Aph_vs_x->SetParameter(6, Aph[2]);

  TF1* Atp_vs_x = new TF1("Atp_vs_x", A_vs_x, xmin, xmax, 4 + npar_tp );
  Atp_vs_x->FixParameter(0, npar_tp);
  Atp_vs_x->FixParameter(1, x[0]);
  Atp_vs_x->FixParameter(2, x[1]);
  Atp_vs_x->FixParameter(3, x[2]);
  Atp_vs_x->SetParameter(4, Atp[0]);
  if ( npar_tp > 1 ) Atp_vs_x->SetParameter(5, Atp[1]);
  if ( npar_tp > 2 ) Atp_vs_x->SetParameter(6, Atp[2]);

  TF1* pTweight = new TF1("pTweight",func_pT_weight,pTmin,pTmax,2);


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////// Preliminary calculations of the polarization parameters, including optional iterative pT reweighting

  for ( int i_iter = 0; i_iter < n_iterations; i_iter++ ) {

    h_pT_data->Reset(); // reset histogram to zero: it is re-filled at each iteration

    //////// fill DATA pT distribution, reweighted for the current best knowledge of the data polarization (step 0: unpolarized)
    for ( int j_ev = 0; j_ev < numDataEvts; j_ev++ ) { //loop only over array of analyzed events

      // weight to eliminate the previously measured polarization (iterative method)
      // used ONLY to produce the pT-distribution ratio data/reference in the zero-polarization hypothesis
      // At first iteration weight is 1; in Nch-dependent analysis there must not be following iterations,
      // so no weight is ever applied to account for data polarization

      double pT = pT_j[j_ev];

      double Ath_curr = Ath_vs_x->Eval(pT);
      double Aph_curr = Aph_vs_x->Eval(pT);
      double Atp_curr = Atp_vs_x->Eval(pT);

      double wUnpol = 4./3. / ( 1. + Ath_curr + (1.-3.*Ath_curr) * costheta2_j[j_ev]
                                + Aph_curr * sintheta2cos2phi_j[j_ev]
                                + Atp_curr * sin2thetacosphi_j[j_ev] );

      h_pT_data->Fill( pT, wS_j[j_ev] * wUnpol );
    } // end of data loop to fill weighted pT distribution

    h_pT_data->Scale( 1./ h_pT_data->Integral() );


    ///////// calculation of pT reweighing function for the reference sample

    h_pT_data->Divide(h_pT_ref);

    pTweight->SetParameter(0, 1.);
    pTweight->SetParameter(1, 0.);

    h_pT_data->SetAxisRange(0.0, 2.0,"Y");
    TCanvas* canvas1 = new TCanvas();
    h_pT_data->Draw(" ");

    h_pT_data->Fit("pTweight","IQ","",pTmin,pTmax);
    //h_pT_data->Fit("pTweight","IQ","",pTmin+pTstep,pTmax-pTstep); //remove first and last cell

    command_option.str("");
    command_option << "pTweight_fitted_iter" << i_iter << ".pdf";

    canvas1->Print( command_option.str().c_str() );
    canvas1->Close();


    ////////// fill REF costh-phi histograms (with full pT-reweighting)

    for (int i = 0; i < n_xbins; i++) {
      h_costhphi_ref[i]->Reset();  // histograms re-filled at each iteration
      h_costhphi_ratio[i] = (TH2D*) h_costhphi_data[i]->Clone(); // also clone data histograms into the histos for the ratio
    }

    for ( int j_ev = 0; j_ev < numRefEvts; j_ev++ ) { //loop only over array of analyzed events

      double xx = x_jR[j_ev];

      double pT = pT_jR[j_ev];

      double xpol = pT; // variable wrt which the reference polarization is known to change
      double wUnpol = ( 3.+lambdath_ref(xpol) ) /3. /(1 + lambdath_ref(xpol) * costheta2_jR[j_ev]
                                                      + lambdaph_ref(xpol) * sintheta2cos2phi_jR[j_ev]
                                                      + lambdatp_ref(xpol) * sin2thetacosphi_jR[j_ev] );
      // weight to eliminate the known polarization of the reference

      double wpT = pTweight->Eval(pT);
      // pT weight

      if ( ! apply_pT_weight_to_ref ) wpT = 1.;

      double weight = wS_jR[j_ev] * wUnpol * wpT;

      for (int i_bin = 0; i_bin < n_xbins; i_bin++) {

        if ( xx > xbin[i_bin] && xx < xbin[i_bin+1] ) {
          h_costhphi_ref[i_bin]->Fill(costheta_jR[j_ev], phi_jR[j_ev], weight);
        }
      }
    } // end of reference-events loop to fill weighted costh-phi histograms

    // normalize reference (data histograms already normalized) and calculate
    // ratios of data/reference histograms to be fitted
    for (int i = 0; i < n_xbins; i++) {
      h_costhphi_ref[i]->Scale( 1./h_costhphi_ref[i]->Integral() );
      h_costhphi_ratio[i]->Divide(h_costhphi_ref[i]);
    }


    //////////// FITS /////////////////////////////////////////////////


    // polarization fits in the n_xbins

    for (int i = 0; i < n_xbins; i++) {

      polFit_A->SetParameter(0, 1.);
      polFit_A->SetParameter(1, Ath_i[i]);
      polFit_A->SetParameter(2, Aph_i[i]);
      polFit_A->SetParameter(3, Atp_i[i]);

      h_costhphi_ratio[i]->Fit("polFit_A","IN0Q");

      Ath_i[i] = polFit_A->GetParameter(1);
      Aph_i[i] = polFit_A->GetParameter(2);
      Atp_i[i] = polFit_A->GetParameter(3);
      dAth_i[i] = polFit_A->GetParError(1);
      dAph_i[i] = polFit_A->GetParError(2);
      dAtp_i[i] = polFit_A->GetParError(3);
    }


    // x dependence

    TCanvas* canvas2 = new TCanvas("canvas2","Polarization",10,20,300,450);

    TGraphErrors* Ath_vs_x_graph = new TGraphErrors( n_xbins, x_i, Ath_i, ex, dAth_i );
    Ath_vs_x_graph->Draw("AP");
    Ath_vs_x_graph->GetYaxis()->SetRangeUser(-1.,1.);

    Ath_vs_x_graph->SetTitle("");
    Ath_vs_x_graph->SetMarkerColor(2);
    Ath_vs_x_graph->SetLineColor(2);
    Ath_vs_x_graph->SetLineWidth(2);
    Ath_vs_x_graph->SetMarkerStyle(21);
    Ath_vs_x_graph->SetMarkerSize(0.3);
    Ath_vs_x_graph->Draw("P");

    TGraphErrors* Aph_vs_x_graph = new TGraphErrors( n_xbins, x_i, Aph_i, ex, dAph_i );
    Aph_vs_x_graph->SetTitle("");
    Aph_vs_x_graph->SetMarkerColor(4);
    Aph_vs_x_graph->SetLineColor(4);
    Aph_vs_x_graph->SetLineWidth(2);
    Aph_vs_x_graph->SetMarkerStyle(21);
    Aph_vs_x_graph->SetMarkerSize(0.3);
    Aph_vs_x_graph->Draw("P same");

    TGraphErrors* Atp_vs_x_graph = new TGraphErrors( n_xbins, x_i, Atp_i, ex, dAtp_i );
    Atp_vs_x_graph->SetTitle("");
    Atp_vs_x_graph->SetMarkerColor(3);
    Atp_vs_x_graph->SetLineColor(3);
    Atp_vs_x_graph->SetLineWidth(2);
    Atp_vs_x_graph->SetMarkerStyle(21);
    Atp_vs_x_graph->SetMarkerSize(0.3);
    Atp_vs_x_graph->Draw("P same");

    Ath_vs_x->SetLineColor(2);
    Ath_vs_x->SetLineStyle(2);
    Ath_vs_x->SetLineWidth(1);
    Ath_vs_x_graph->Fit("Ath_vs_x","Q");

    Aph_vs_x->SetLineColor(4);
    Aph_vs_x->SetLineStyle(2);
    Aph_vs_x->SetLineWidth(1);
    Aph_vs_x_graph->Fit("Aph_vs_x","Q");

    Atp_vs_x->SetLineColor(3);
    Atp_vs_x->SetLineStyle(2);
    Atp_vs_x->SetLineWidth(1);
    Atp_vs_x_graph->Fit("Atp_vs_x","Q");

    command_option.str("");
    command_option << "Polarization_fitted_iter" << i_iter << ".pdf";

    canvas2->Print( command_option.str().c_str() );
    canvas2->Close();


    // re-set Ath[i] from fit results

    Ath[0]  = Ath_vs_x->GetParameter(4);
    dAth[0] = Ath_vs_x->GetParError(4);
    if ( npar_th > 1 ) { Ath[1]  = Ath_vs_x->GetParameter(5);
      dAth[1]  = Ath_vs_x->GetParError(5);  }
    if ( npar_th > 2 ) { Ath[2]  = Ath_vs_x->GetParameter(6);
      dAth[2]  = Ath_vs_x->GetParError(6);  }

    Aph[0]  = Aph_vs_x->GetParameter(4);
    dAph[0] = Aph_vs_x->GetParError(4);
    if ( npar_th > 1 ) { Aph[1]  = Aph_vs_x->GetParameter(5);
      dAph[1]  = Aph_vs_x->GetParError(5);  }
    if ( npar_th > 2 ) { Aph[2]  = Aph_vs_x->GetParameter(6);
      dAph[2]  = Aph_vs_x->GetParError(6);  }

    Atp[0]  = Atp_vs_x->GetParameter(4);
    dAtp[0] = Atp_vs_x->GetParError(4);
    if ( npar_th > 1 ) { Atp[1]  = Atp_vs_x->GetParameter(5);
      dAtp[1]  = Atp_vs_x->GetParError(5);  }
    if ( npar_th > 2 ) { Atp[2]  = Atp_vs_x->GetParameter(6);
      dAtp[2]  = Atp_vs_x->GetParError(6);  }


  } // end of iterations
  ///////////////////////////


  // plot of starting values represented in terms of lambdas
  // this is a pure representation. It will not be used in the fit
  // error calculations for lph and ltp are approximate (neglect correlation with lth)

  // abscissa were defined as  double* x   = new double[3];
  double* lth = new double[3];
  double* lph = new double[3];
  double* ltp = new double[3];

  double* dlth = new double[3];
  double* dlph = new double[3];
  double* dltp = new double[3];

  for ( int i = 0; i < npar_th; i++ ) {

    lth[i]  = (1. - 3.* Ath[i]) / ( 1. + Ath[i] );
    dlth[i] = 4. * dAth[i] / pow( 1. + Ath[i], 2. );
  }

  for ( int i = 0; i < npar_ph; i++ ) {

    double Ath_x = Ath_vs_x->Eval(x[i]);

    lph[i]  =  Aph[i] / ( 1. + Ath_x );
    dlph[i] = dAph[i] / ( 1. + Ath_x );
  }

  for ( int i = 0; i < npar_tp; i++ ) {

    double Ath_x = Ath_vs_x->Eval(x[i]);

    ltp[i]  =  Atp[i] / ( 1. + Ath_x );
    dltp[i] = dAtp[i] / ( 1. + Ath_x );
  }

  TCanvas* canvas3 = new TCanvas("canvas3","Lambdas",10,20,300,450);

  TGraphErrors* lth_vs_x_graph = new TGraphErrors( npar_th, x, lth, ex, dlth );
  lth_vs_x_graph->Draw("AP");
  lth_vs_x_graph->GetYaxis()->SetRangeUser(-1.,1.);
  lth_vs_x_graph->SetTitle("");
  lth_vs_x_graph->SetMarkerColor(2);
  lth_vs_x_graph->SetLineColor(2);
  lth_vs_x_graph->SetLineWidth(2);
  lth_vs_x_graph->SetMarkerStyle(21);
  lth_vs_x_graph->SetMarkerSize(0.3);
  lth_vs_x_graph->Draw("P");

  TGraphErrors* lph_vs_x_graph = new TGraphErrors( npar_ph, x, lph, ex, dlph );
  lph_vs_x_graph->SetTitle("");
  lph_vs_x_graph->SetMarkerColor(4);
  lph_vs_x_graph->SetLineColor(4);
  lph_vs_x_graph->SetLineWidth(2);
  lph_vs_x_graph->SetMarkerStyle(21);
  lph_vs_x_graph->SetMarkerSize(0.3);
  lph_vs_x_graph->Draw("P same");

  TGraphErrors* ltp_vs_x_graph = new TGraphErrors( npar_tp, x, ltp, ex, dltp );
  ltp_vs_x_graph->SetTitle("");
  ltp_vs_x_graph->SetMarkerColor(3);
  ltp_vs_x_graph->SetLineColor(3);
  ltp_vs_x_graph->SetLineWidth(2);
  ltp_vs_x_graph->SetMarkerStyle(21);
  ltp_vs_x_graph->SetMarkerSize(0.3);
  ltp_vs_x_graph->Draw("P same");

  command_option.str("");
  command_option << "Lambdas_starting.pdf";

  canvas3->Print( command_option.str().c_str() );
  canvas3->Close();

  // dataSample->Delete(); // ROOT v6 crashes here if the TTrees are deleted
  dataFile->Close();
  // refSample->Delete();
  refFile->Close();


  cout << endl;
  cout << "Output of preliminary fits:"<< endl;
  cout << endl;
  for ( int i = 0; i < npar_th; i++ ) {
    cout << "x = " << x[i] << ":   lth = " << lth[i] << " +/- " << dlth[i] << endl;
  }
  for ( int i = 0; i < npar_ph; i++ ) {
    cout << "x = " << x[i] << ":   lph = " << lph[i] << " +/- " << dlph[i] << endl;
  }
  for ( int i = 0; i < npar_tp; i++ ) {
    cout << "x = " << x[i] << ":   ltp = " << ltp[i] << " +/- " << dltp[i] << endl;
  }

  ///////////////////////////////////////////////////////////////////
  /////////// PREPARATION OF FINAL MULTIDIMENSIONAL FIT /////////////
  ///////////////////////////////////////////////////////////////////

  const double deltaAth_min = 0.05;   // prevent starting scan sigmas from being too small
  const double deltaAph_min = 0.03;
  const double deltaAtp_min = 0.03;

  // enlarge sigmas for initial scan

  for ( int i = 0; i < npar_th; i++ ) {

    dAth[i] *= scan_sigma_margin_factor[0];
    if ( dAth[i] < deltaAth_min ) dAth[i] = deltaAth_min;
  }

  for ( int i = 0; i < npar_ph; i++ ) {

    dAph[i] *= scan_sigma_margin_factor[0];
    if ( dAph[i] < deltaAph_min ) dAph[i] = deltaAph_min;
  }

  for ( int i = 0; i < npar_tp; i++ ) {

    dAtp[i] *= scan_sigma_margin_factor[0];
    if ( dAtp[i] < deltaAtp_min ) dAtp[i] = deltaAtp_min;
  }

  ////////////////////////////////////////////////////////////////////////////////////
  // loop over reference events to calculate reference averages

  double avg_norm = 0.; // sum of weights


  double avg_x = 0.;
  double avg_x_square = 0.; // averages of "square" quantities needed to calculate uncertainties
                            // obviously in this case avg_x_square = avg_x2, but notation is kept symmetric
  double avg_x2 = 0.;
  double avg_x2_square = 0.;


  double avg_costheta2 = 0.;
  double avg_costheta2_square = 0.;

  double avg_xcostheta2 = 0.;
  double avg_xcostheta2_square = 0.;

  double avg_x2costheta2 = 0.;
  double avg_x2costheta2_square = 0.;


  double avg_sintheta2cos2phi = 0.;
  double avg_sintheta2cos2phi_square = 0.;

  double avg_xsintheta2cos2phi = 0.;
  double avg_xsintheta2cos2phi_square = 0.;

  double avg_x2sintheta2cos2phi = 0.;
  double avg_x2sintheta2cos2phi_square = 0.;


  double avg_sin2thetacosphi = 0.;
  double avg_sin2thetacosphi_square = 0.;

  double avg_xsin2thetacosphi = 0.;
  double avg_xsin2thetacosphi_square = 0.;

  double avg_x2sin2thetacosphi = 0.;
  double avg_x2sin2thetacosphi_square = 0.;


  for ( int j_ev = 0; j_ev < numRefEvts; j_ev++ ) { //loop only over array of analyzed ref events

    double x1 = x_jR[j_ev];
    double x2 = x1*x1;

    double pT = pT_jR[j_ev];


    double costheta2 = costheta2_jR[j_ev];
    double sintheta2cos2phi = sintheta2cos2phi_jR[j_ev];
    double sin2thetacosphi = sin2thetacosphi_jR[j_ev];

    double costheta2_square = costheta2*costheta2;
    double sintheta2cos2phi_square = sintheta2cos2phi*sintheta2cos2phi;
    double sin2thetacosphi_square = sin2thetacosphi*sin2thetacosphi;

    double xpol = pT; // variable wrt which the reference polarization is known to change
    double wUnpol = ( 3.+lambdath_ref(xpol) ) /3. /(1 + lambdath_ref(xpol) * costheta2
                                                    + lambdaph_ref(xpol) * sintheta2cos2phi
                                                    + lambdatp_ref(xpol) * sin2thetacosphi );
    // weight to eliminate the known polarization of the reference

    double wpT = pTweight->Eval(pT);
    // pT weight

    if ( ! apply_pT_weight_to_ref ) wpT = 1.;

    double wTot = wS_jR[j_ev] * wUnpol * wpT;


    avg_norm += wTot;

    avg_x += wTot * x1;
    avg_x2 += wTot * x2;

    avg_costheta2        += wTot * costheta2;
    avg_xcostheta2       += wTot * x1 * costheta2;
    avg_x2costheta2      += wTot * x2 * costheta2;

    avg_sintheta2cos2phi   += wTot * sintheta2cos2phi;
    avg_xsintheta2cos2phi  += wTot * x1 * sintheta2cos2phi;
    avg_x2sintheta2cos2phi += wTot * x2 * sintheta2cos2phi;

    avg_sin2thetacosphi   += wTot * sin2thetacosphi;
    avg_xsin2thetacosphi  += wTot * x1 * sin2thetacosphi;
    avg_x2sin2thetacosphi += wTot * x2 * sin2thetacosphi;


    avg_x_square += wTot * x2;
    avg_x2_square += wTot * x2*x2;

    avg_costheta2_square        += wTot * costheta2_square;
    avg_xcostheta2_square       += wTot * x2 * costheta2_square;
    avg_x2costheta2_square      += wTot * x2*x2 * costheta2_square;

    avg_sintheta2cos2phi_square    += wTot * sintheta2cos2phi_square;
    avg_xsintheta2cos2phi_square   += wTot * x2 * sintheta2cos2phi_square;
    avg_x2sintheta2cos2phi_square  += wTot * x2*x2 * sintheta2cos2phi_square;

    avg_sin2thetacosphi_square    += wTot * sin2thetacosphi_square;
    avg_xsin2thetacosphi_square   += wTot * x2 * sin2thetacosphi_square;
    avg_x2sin2thetacosphi_square  += wTot * x2*x2 * sin2thetacosphi_square;

  } // end of reference-events loop to calculate event probability normalizations


  // averages of the distributions

  avg_x /= avg_norm;
  avg_x2 /= avg_norm;

  avg_costheta2 /= avg_norm;
  avg_xcostheta2 /= avg_norm;
  avg_x2costheta2 /= avg_norm;

  avg_sintheta2cos2phi /= avg_norm;
  avg_xsintheta2cos2phi /= avg_norm;
  avg_x2sintheta2cos2phi /= avg_norm;

  avg_sin2thetacosphi /= avg_norm;
  avg_xsin2thetacosphi /= avg_norm;
  avg_x2sin2thetacosphi /= avg_norm;


  avg_x_square /= avg_norm;
  avg_x2_square /= avg_norm;

  avg_costheta2_square  /= avg_norm;
  avg_xcostheta2_square  /= avg_norm;
  avg_x2costheta2_square  /= avg_norm;

  avg_sintheta2cos2phi_square   /= avg_norm;
  avg_xsintheta2cos2phi_square   /= avg_norm;
  avg_x2sintheta2cos2phi_square   /= avg_norm;

  avg_sin2thetacosphi_square   /= avg_norm;
  avg_xsin2thetacosphi_square   /= avg_norm;
  avg_x2sin2thetacosphi_square   /= avg_norm;

  // std deviations of the distributions (multidimensional fit only)

  double sig_x = sqrt((avg_x_square - avg_x*avg_x)/ avg_norm);  //  division by avg_norm to obtain std dev OF THE MEAN
  double sig_x2 = sqrt((avg_x2_square - avg_x2*avg_x2)/ avg_norm);

  double sig_costheta2 = sqrt((avg_costheta2_square - avg_costheta2*avg_costheta2)/ avg_norm);
  double sig_xcostheta2 = sqrt((avg_xcostheta2_square - avg_xcostheta2*avg_xcostheta2)/ avg_norm);
  double sig_x2costheta2 = sqrt((avg_x2costheta2_square - avg_x2costheta2*avg_x2costheta2)/ avg_norm);

  double sig_sintheta2cos2phi = sqrt((avg_sintheta2cos2phi_square - avg_sintheta2cos2phi*avg_sintheta2cos2phi)/ avg_norm);
  double sig_xsintheta2cos2phi = sqrt((avg_xsintheta2cos2phi_square - avg_xsintheta2cos2phi*avg_xsintheta2cos2phi)/ avg_norm);
  double sig_x2sintheta2cos2phi = sqrt((avg_x2sintheta2cos2phi_square - avg_x2sintheta2cos2phi*avg_x2sintheta2cos2phi)/ avg_norm);

  double sig_sin2thetacosphi = sqrt((avg_sin2thetacosphi_square - avg_sin2thetacosphi*avg_sin2thetacosphi)/ avg_norm);
  double sig_xsin2thetacosphi = sqrt((avg_xsin2thetacosphi_square - avg_xsin2thetacosphi*avg_xsin2thetacosphi)/ avg_norm);
  double sig_x2sin2thetacosphi = sqrt((avg_x2sin2thetacosphi_square - avg_x2sin2thetacosphi*avg_x2sin2thetacosphi)/ avg_norm);



  //////////////////////////////////////////////////////////////////////////////////////////////////
  // definition of the "systematic" (actually statistic) variations of the reference averages

  const int n_ref_quantities = 11; // n. of reference quantities to be varied
  const int n_variations_of_reference = 35; // n. of systematic variations considered for the reference quantities
  // (+-1 sigma for one parameter at a time, plus no variation)
  // a different PPD will be calculated (and independently normalized) for each case

  double* avg_x_var                 = new double[n_variations_of_reference];
  double* avg_x2_var                = new double[n_variations_of_reference];

  double* avg_costheta2_var         = new double[n_variations_of_reference];
  double* avg_xcostheta2_var        = new double[n_variations_of_reference];
  double* avg_x2costheta2_var       = new double[n_variations_of_reference];

  double* avg_sintheta2cos2phi_var  = new double[n_variations_of_reference];
  double* avg_xsintheta2cos2phi_var = new double[n_variations_of_reference];
  double* avg_x2sintheta2cos2phi_var= new double[n_variations_of_reference];

  double* avg_sin2thetacosphi_var   = new double[n_variations_of_reference];
  double* avg_xsin2thetacosphi_var  = new double[n_variations_of_reference];
  double* avg_x2sin2thetacosphi_var = new double[n_variations_of_reference];


  // each "test" (matrix line) = 1 sigma variation of one quantity at a time; test zero = no variation
  double variation[n_variations_of_reference][n_ref_quantities] = { { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    // some correlated variations:
                                                                    { 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    {-1.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 1., 1., 1., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0.,-1.,-1.,-1., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0.,-1.,-1.,-1., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0.,-1.,-1.,-1. },
                                                                    { 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0. },
                                                                    {-1., 0., 0.,-1., 0., 0.,-1., 0., 0.,-1., 0. },
                                                                    { 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1. },
                                                                    { 0.,-1., 0., 0.,-1., 0., 0.,-1., 0., 0.,-1. },
                                                                    // all fully uncorrelated variations:
                                                                    { 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    {-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0.,-1., 0., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0.,-1., 0., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0.,-1., 0., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0.,-1., 0., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0.,-1., 0., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 0.,-1., 0. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1. },
                                                                    { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,-1. } };

  for ( int i_variation = 0; i_variation < n_variations_of_reference; i_variation++ ) {

    avg_x_var                 [i_variation] = avg_x                  + sig_x                  * variation[i_variation][0];
    avg_x2_var                [i_variation] = avg_x2                 + sig_x2                 * variation[i_variation][1];

    avg_costheta2_var         [i_variation] = avg_costheta2          + sig_costheta2          * variation[i_variation][2];
    avg_xcostheta2_var        [i_variation] = avg_xcostheta2         + sig_xcostheta2         * variation[i_variation][3];
    avg_x2costheta2_var       [i_variation] = avg_x2costheta2        + sig_x2costheta2        * variation[i_variation][4];

    avg_sintheta2cos2phi_var  [i_variation] = avg_sintheta2cos2phi   + sig_sintheta2cos2phi   * variation[i_variation][5];
    avg_xsintheta2cos2phi_var [i_variation] = avg_xsintheta2cos2phi  + sig_xsintheta2cos2phi  * variation[i_variation][6];
    avg_x2sintheta2cos2phi_var[i_variation] = avg_x2sintheta2cos2phi + sig_x2sintheta2cos2phi * variation[i_variation][7];

    avg_sin2thetacosphi_var   [i_variation] = avg_sin2thetacosphi    + sig_sin2thetacosphi    * variation[i_variation][8];
    avg_xsin2thetacosphi_var  [i_variation] = avg_xsin2thetacosphi   + sig_xsin2thetacosphi   * variation[i_variation][9];
    avg_x2sin2thetacosphi_var [i_variation] = avg_x2sin2thetacosphi  + sig_x2sin2thetacosphi  * variation[i_variation][10];
  }

  // some variations will not cause any effect.
  // For example, in a linear fit the terms with "x2" do not exist and do not affect the likelihood normalization
  // Each "varied" PPD will be "added" to the sum of the varied PPDs only if the the variation has changed the
  // normalization (and therefore the PPD shape) wrt the unvaried one, otherwise the unvaried likelihood would be counted
  // more times, reducing the smearing effect.


  ////// ARRAYS storing the scanned parameter values and the PPD

  double* Ath1_arr = new double[n_scan_points]; // arrays of parameter values, one for each event
  double* Ath2_arr = new double[n_scan_points];
  double* Ath3_arr = new double[n_scan_points];
  double* Aph1_arr = new double[n_scan_points];
  double* Aph2_arr = new double[n_scan_points];
  double* Aph3_arr = new double[n_scan_points];
  double* Atp1_arr = new double[n_scan_points];
  double* Atp2_arr = new double[n_scan_points];
  double* Atp3_arr = new double[n_scan_points];

  double** lnPPD_arr = new double*[n_scan_points];
  for (int i = 0; i < n_scan_points; ++i) {
    lnPPD_arr[i] = new double[n_variations_of_reference]; // array of PPD values for each type of systematic variation
  }

  double* lnPPD_max = new double[n_variations_of_reference];
  memset(lnPPD_max, -1.E20, n_variations_of_reference*sizeof(double) );
  // maximum value of PPD for each type of systematic variation

  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];

  // constants appearing in the polynomial parametrization A = a + b * x + c * x^2

  double c1 = 1./((x1-x2)*(x1-x3));
  double c2 = 1./((x2-x1)*(x2-x3));
  double c3 = 1./((x3-x1)*(x3-x2));

  double a1 = x2*x3*c1;
  double a2 = x1*x3*c2;
  double a3 = x1*x2*c3;

  double b1 = -(x2+x3)*c1;
  double b2 = -(x1+x3)*c2;
  double b3 = -(x1+x2)*c3;


  ///////////////////////////////////////////////////////////////////
  /////////// FINAL MULTIDIMENSIONAL FIT ////////////////////////////
  ///////////////////////////////////////////////////////////////////

  n_step = n_scan_points/50;

  cout << endl;
  cout << "Parameter scanning"<< endl;
  cout << "Calculating " << n_scan_points << " PPD values"<< endl;
  cout << "------------------------------------------------------------" << endl;
  cout << "Progress: ";


  int i_scan_phase = 0; // fist phase of parameter scan, using broadest Gaussians

  for(int i_scan = 0; i_scan < n_scan_points; i_scan++){

    if ( i_scan%n_step == 0 ) cout << "X";

    // log(weight) for PPD, compensating for the non-flat scanning of the parameter space
    double ln_scan_weight = 0.;

    double Ath1, Ath2, Ath3;
    Ath1 = gRandom->Gaus(Ath[0], dAth[0]); ln_scan_weight -= log( TMath::Gaus(Ath1, Ath[0], dAth[0]) );
    Ath2 = Ath1; Ath3 = Ath1;
    if ( npar_th == 2 ) { do { Ath2 = gRandom->Gaus(Ath[1], dAth[1]); } while ( Ath2 > Ath_max || Ath2 < Ath_min );
      ln_scan_weight -= log( TMath::Gaus(Ath2, Ath[1], dAth[1]) );
      Ath3 = Ath1 + (x3-x1) * (Ath2-Ath1) / (x2-x1); }
    else if ( npar_th == 3 ) { do { Ath2 = gRandom->Gaus(Ath[1], dAth[1]); } while ( Ath2 > Ath_max || Ath2 < Ath_min );
      ln_scan_weight -= log( TMath::Gaus(Ath2, Ath[1], dAth[1]) );
      do { Ath3 = gRandom->Gaus(Ath[2], dAth[2]); } while ( Ath3 > Ath_max || Ath3 < Ath_min );
      ln_scan_weight -= log( TMath::Gaus(Ath3, Ath[2], dAth[2]) ); }
    Ath1_arr[i_scan] = Ath1;
    Ath2_arr[i_scan] = Ath2;
    Ath3_arr[i_scan] = Ath3;

    double Aph1, Aph2, Aph3;
    Aph1 = gRandom->Gaus(Aph[0], dAph[0]); ln_scan_weight -= log( TMath::Gaus(Aph1, Aph[0], dAph[0]) );
    Aph2 = Aph1; Aph3 = Aph1;
    if ( npar_ph == 2 ) { do { Aph2 = gRandom->Gaus(Aph[1], dAph[1]); } while ( Aph2 > Aph_max || Aph2 < Aph_min );
      ln_scan_weight -= log( TMath::Gaus(Aph2, Aph[1], dAph[1]) );
      Aph3 = Aph1 + (x3-x1) * (Aph2-Aph1) / (x2-x1); }
    else if ( npar_ph == 3 ) { do { Aph2 = gRandom->Gaus(Aph[1], dAph[1]); } while ( Aph2 > Aph_max || Aph2 < Aph_min );
      ln_scan_weight -= log( TMath::Gaus(Aph2, Aph[1], dAph[1]) );
      do { Aph3 = gRandom->Gaus(Aph[2], dAph[2]); } while ( Aph3 > Aph_max || Aph3 < Aph_min );
      ln_scan_weight -= log( TMath::Gaus(Aph3, Aph[2], dAph[2]) ); }
    Aph1_arr[i_scan] = Aph1;
    Aph2_arr[i_scan] = Aph2;
    Aph3_arr[i_scan] = Aph3;

    double Atp1, Atp2, Atp3;
    Atp1 = gRandom->Gaus(Atp[0], dAtp[0]); ln_scan_weight -= log( TMath::Gaus(Atp1, Atp[0], dAtp[0]) );
    Atp2 = Atp1; Atp3 = Atp1;
    if ( npar_tp == 2 ) { do { Atp2 = gRandom->Gaus(Atp[1], dAtp[1]); } while ( Atp2 > Atp_max || Atp2 < Atp_min );
      ln_scan_weight -= log( TMath::Gaus(Atp2, Atp[1], dAtp[1]) );
      Atp3 = Atp1 + (x3-x1) * (Atp2-Atp1) / (x2-x1); }
    else if ( npar_tp == 3 ) { do { Atp2 = gRandom->Gaus(Atp[1], dAtp[1]); } while ( Atp2 > Atp_max || Atp2 < Atp_min );
      ln_scan_weight -= log( TMath::Gaus(Atp2, Atp[1], dAtp[1]) );
      do { Atp3 = gRandom->Gaus(Atp[2], dAtp[2]); } while ( Atp3 > Atp_max || Atp3 < Atp_min );
      ln_scan_weight -= log( TMath::Gaus(Atp3, Atp[2], dAtp[2]) ); }
    Atp1_arr[i_scan] = Atp1;
    Atp2_arr[i_scan] = Atp2;
    Atp3_arr[i_scan] = Atp3;

    double a_th = Ath1*a1 + Ath2*a2 + Ath3*a3;
    double b_th = Ath1*b1 + Ath2*b2 + Ath3*b3;
    double c_th = Ath1*c1 + Ath2*c2 + Ath3*c3;

    double a_ph = Aph1*a1 + Aph2*a2 + Aph3*a3;
    double b_ph = Aph1*b1 + Aph2*b2 + Aph3*b3;
    double c_ph = Aph1*c1 + Aph2*c2 + Aph3*c3;

    double a_tp = Atp1*a1 + Atp2*a2 + Atp3*a3;
    double b_tp = Atp1*b1 + Atp2*b2 + Atp3*b3;
    double c_tp = Atp1*c1 + Atp2*c2 + Atp3*c3;

    double lnPPD = 0.;

    double xval, xval2, Ath_j, Aph_j, Atp_j;

    ///////// Event loop ////////////////////////////////////////////////////////////
    for ( int j = 0; j < numDataEvts; j++ ) {

      xval  = x_j[j];
      xval2 = xval*xval;

      Ath_j = a_th  + b_th * xval  + c_th * xval2;
      Aph_j = a_ph  + b_ph * xval  + c_ph * xval2;
      Atp_j = a_tp  + b_tp * xval  + c_tp * xval2;

      lnPPD +=  wS_j[j] * log(     1. + Ath_j
                                   + ( 1. - 3. * Ath_j ) * costheta2_j[j]
                                   +        Aph_j        * sintheta2cos2phi_j[j]
                                   +        Atp_j        * sin2thetacosphi_j[j]  );
    }////////////////////////////////////////////////////////////////////////////////

    lnPPD += ln_scan_weight;

    // lnPPD normalization and calculation, for each kind of systematic variation

    for ( int i_variation = 0; i_variation < n_variations_of_reference; i_variation++ ) {

      double lnPPD_temp = lnPPD
        - n_dataSignal* log(     1. + a_th + b_th * avg_x_var[i_variation] + c_th * avg_x2_var[i_variation]
                                 + ( 1. - 3. * a_th ) * avg_costheta2_var[i_variation]
                                 - 3. * b_th * avg_xcostheta2_var[i_variation]
                                 - 3. * c_th * avg_x2costheta2_var[i_variation]

                                 + a_ph * avg_sintheta2cos2phi_var[i_variation]
                                 + b_ph * avg_xsintheta2cos2phi_var[i_variation]
                                 + c_ph * avg_x2sintheta2cos2phi_var[i_variation]

                                 + a_tp * avg_sin2thetacosphi_var[i_variation]
                                 + b_tp * avg_xsin2thetacosphi_var[i_variation]
                                 + c_tp * avg_x2sin2thetacosphi_var[i_variation]    );

      lnPPD_arr[i_scan][i_variation] = lnPPD_temp;

      if ( lnPPD_temp > lnPPD_max[i_variation] ) lnPPD_max[i_variation] = lnPPD_temp;
    }

    // check change of phase

    bool end_of_phase =

      ( double(i_scan)/double(n_scan_points) > scan_progress_vs_phase[i_scan_phase] );

    // at the end of each scanning phase, calculate averages and sigmas
    // to redefine an optimized Gaussian scanning width

    if ( end_of_phase ) {

      i_scan_phase = i_scan_phase + 1;

      // re-scan so-far calculated PPD values, normalize them and calculate the
      // temporary averages and sigmas of the parameter distributions

      // sums for calculation of temporary centres and sigmas of the 1D PPDs
      double wPPD_sum = 0.;

      double Ath1_avg = 0.;
      double Ath1_sqavg = 0.;
      double Ath2_avg = 0.;
      double Ath2_sqavg = 0.;
      double Ath3_avg = 0.;
      double Ath3_sqavg = 0.;

      double Aph1_avg = 0.;
      double Aph1_sqavg = 0.;
      double Aph2_avg = 0.;
      double Aph2_sqavg = 0.;
      double Aph3_avg = 0.;
      double Aph3_sqavg = 0.;

      double Atp1_avg = 0.;
      double Atp1_sqavg = 0.;
      double Atp2_avg = 0.;
      double Atp2_sqavg = 0.;
      double Atp3_avg = 0.;
      double Atp3_sqavg = 0.;

      for(int i_rescan = 0; i_rescan < i_scan; i_rescan++){

        double PPD0 = exp(lnPPD_arr[i_rescan][0]-lnPPD_max[0]);

        double wPPD_temp = 0.;

        for ( int i_variation = 1; i_variation < n_variations_of_reference; i_variation++ ) {

          double PPDi = exp(lnPPD_arr[i_rescan][i_variation]-lnPPD_max[i_variation]);

          if ( fabs(PPDi - PPD0) > 0.001 * PPD0 ) wPPD_temp += PPDi; // sum of equally normalized PPDs
          // of all non-null systematic variations
        }

        double Ath1_temp = Ath1_arr[i_rescan];
        double Ath2_temp = Ath2_arr[i_rescan];
        double Ath3_temp = Ath3_arr[i_rescan];

        double Aph1_temp = Aph1_arr[i_rescan];
        double Aph2_temp = Aph2_arr[i_rescan];
        double Aph3_temp = Aph3_arr[i_rescan];

        double Atp1_temp = Atp1_arr[i_rescan];
        double Atp2_temp = Atp2_arr[i_rescan];
        double Atp3_temp = Atp3_arr[i_rescan];

        wPPD_sum += wPPD_temp;

        Ath1_avg   += wPPD_temp * Ath1_temp;
        Ath1_sqavg += wPPD_temp * Ath1_temp*Ath1_temp;
        Ath2_avg   += wPPD_temp * Ath2_temp;
        Ath2_sqavg += wPPD_temp * Ath2_temp*Ath2_temp;
        Ath3_avg   += wPPD_temp * Ath3_temp;
        Ath3_sqavg += wPPD_temp * Ath3_temp*Ath3_temp;

        Aph1_avg   += wPPD_temp * Aph1_temp;
        Aph1_sqavg += wPPD_temp * Aph1_temp*Aph1_temp;
        Aph2_avg   += wPPD_temp * Aph2_temp;
        Aph2_sqavg += wPPD_temp * Aph2_temp*Aph2_temp;
        Aph3_avg   += wPPD_temp * Aph3_temp;
        Aph3_sqavg += wPPD_temp * Aph3_temp*Aph3_temp;

        Atp1_avg   += wPPD_temp * Atp1_temp;
        Atp1_sqavg += wPPD_temp * Atp1_temp*Atp1_temp;
        Atp2_avg   += wPPD_temp * Atp2_temp;
        Atp2_sqavg += wPPD_temp * Atp2_temp*Atp2_temp;
        Atp3_avg   += wPPD_temp * Atp3_temp;
        Atp3_sqavg += wPPD_temp * Atp3_temp*Atp3_temp;
      }

      // re-define scan centres and sigmas

      Ath[0] = Ath1_avg / wPPD_sum;
      dAth[0] = sqrt( fabs( Ath1_sqavg / wPPD_sum - Ath[0]*Ath[0] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Ath[1] = Ath2_avg / wPPD_sum;
      dAth[1] = sqrt( fabs( Ath2_sqavg / wPPD_sum - Ath[1]*Ath[1] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Ath[2] = Ath3_avg / wPPD_sum;
      dAth[2] = sqrt( fabs( Ath3_sqavg / wPPD_sum - Ath[2]*Ath[2] ) ) * scan_sigma_margin_factor[i_scan_phase];

      Aph[0] = Aph1_avg / wPPD_sum;
      dAph[0] = sqrt( fabs( Aph1_sqavg / wPPD_sum - Aph[0]*Aph[0] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Aph[1] = Aph2_avg / wPPD_sum;
      dAph[1] = sqrt( fabs( Aph2_sqavg / wPPD_sum - Aph[1]*Aph[1] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Aph[2] = Aph3_avg / wPPD_sum;
      dAph[2] = sqrt( fabs( Aph3_sqavg / wPPD_sum - Aph[2]*Aph[2] ) ) * scan_sigma_margin_factor[i_scan_phase];

      Atp[0] = Atp1_avg / wPPD_sum;
      dAtp[0] = sqrt( fabs( Atp1_sqavg / wPPD_sum - Atp[0]*Atp[0] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Atp[1] = Atp2_avg / wPPD_sum;
      dAtp[1] = sqrt( fabs( Atp2_sqavg / wPPD_sum - Atp[1]*Atp[1] ) ) * scan_sigma_margin_factor[i_scan_phase];
      Atp[2] = Atp3_avg / wPPD_sum;
      dAtp[2] = sqrt( fabs( Atp3_sqavg / wPPD_sum - Atp[2]*Atp[2] ) ) * scan_sigma_margin_factor[i_scan_phase];
      /*
        cout << endl;
        for ( int i = 0; i < npar_th; i++ ) {
        cout << "x = " << x[i] << ":   Ath = " << Ath[i] << " +/- " << dAth[i] << endl;
        }
        for ( int i = 0; i < npar_ph; i++ ) {
        cout << "x = " << x[i] << ":   Aph = " << Aph[i] << " +/- " << dAph[i] << endl;
        }
        for ( int i = 0; i < npar_tp; i++ ) {
        cout << "x = " << x[i] << ":   Atp = " << Atp[i] << " +/- " << dAtp[i] << endl;
        }
      */
    }


  } // end of parameter scan
  //////////////////////////////////////////////////////////////////////
  cout << endl;


  ///////////////////////////////////////////////////////////////////////////////////
  //// PPD (output) ntuple //////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  TFile* ppdFile = new TFile("polPPD.root", "RECREATE", "polPPD");

  // histogram with valus of the 3 support points and qualification of the x variable
  TH1D* xi_values  = new TH1D("xi_values", "xi_values", 4, 0.5, 4.5 );
  xi_values->SetBinContent(1, x[0] );
  xi_values->SetBinContent(2, x[1] );
  xi_values->SetBinContent(3, x[2] );
  xi_values->SetBinContent(4, is_analysis_vs_pT?1:0 ); // pT: bin 4 has content 1; Nch: bin 4 has content 0

  // output ntuple
  TTree* polPPD = new TTree("polPPD","polPPD");

  double Ath1;      polPPD->Branch("Ath1",   &Ath1,   "Ath1/D");
  double Ath2;      polPPD->Branch("Ath2",   &Ath2,   "Ath2/D");
  double Ath3;      polPPD->Branch("Ath3",   &Ath3,   "Ath3/D");

  double Aph1;      polPPD->Branch("Aph1",   &Aph1,   "Aph1/D");
  double Aph2;      polPPD->Branch("Aph2",   &Aph2,   "Aph2/D");
  double Aph3;      polPPD->Branch("Aph3",   &Aph3,   "Aph3/D");

  double Atp1;      polPPD->Branch("Atp1",   &Atp1,   "Atp1/D");
  double Atp2;      polPPD->Branch("Atp2",   &Atp2,   "Atp2/D");
  double Atp3;      polPPD->Branch("Atp3",   &Atp3,   "Atp3/D");

  double wPPD;      polPPD->Branch("wPPD",   &wPPD,   "wPPD/D");


  double  lth1,  lth2,  lth3,  lph1,  lph2,  lph3,  ltp1,  ltp2,  ltp3;
  double dlth1, dlth2, dlth3, dlph1, dlph2, dlph3, dltp1, dltp2, dltp3;

  /////////////////////////////////////////////////////////////////////
  // now re-scan arrays to calculate the PPD, print results
  // and fill ntuple

  double wPPD_sum = 0.;

  double lth1_avg = 0.;
  double lth1_sqavg = 0.;
  double lth2_avg = 0.;
  double lth2_sqavg = 0.;
  double lth3_avg = 0.;
  double lth3_sqavg = 0.;

  double lph1_avg = 0.;
  double lph1_sqavg = 0.;
  double lph2_avg = 0.;
  double lph2_sqavg = 0.;
  double lph3_avg = 0.;
  double lph3_sqavg = 0.;

  double ltp1_avg = 0.;
  double ltp1_sqavg = 0.;
  double ltp2_avg = 0.;
  double ltp2_sqavg = 0.;
  double ltp3_avg = 0.;
  double ltp3_sqavg = 0.;

  for(int i_scan = 0; i_scan < n_scan_points; i_scan++){

    wPPD = 0.;

    double PPD0 = exp(lnPPD_arr[i_scan][0]-lnPPD_max[0]);

    for ( int i_variation = 1; i_variation < n_variations_of_reference; i_variation++ ) {

      double PPDi = exp(lnPPD_arr[i_scan][i_variation]-lnPPD_max[i_variation]);

      if ( fabs(PPDi - PPD0) > 0.001 * PPD0 ) wPPD += PPDi; // sum of equally normalized PPDs
      // of all non-null systematic variations
    }

    //remove smearing due to reference variation
    //wPPD = PPD0;

    Ath1 = Ath1_arr[i_scan];
    Ath2 = Ath2_arr[i_scan];
    Ath3 = Ath3_arr[i_scan];

    Aph1 = Aph1_arr[i_scan];
    Aph2 = Aph2_arr[i_scan];
    Aph3 = Aph3_arr[i_scan];

    Atp1 = Atp1_arr[i_scan];
    Atp2 = Atp2_arr[i_scan];
    Atp3 = Atp3_arr[i_scan];

    lth1 = (1. - 3.* Ath1) / ( 1. + Ath1 );
    lph1 =      Aph1       / ( 1. + Ath1 );
    ltp1 =      Atp1       / ( 1. + Ath1 );

    lth2 = (1. - 3.* Ath2) / ( 1. + Ath2 );
    lph2 =      Aph2       / ( 1. + Ath2 );
    ltp2 =      Atp2       / ( 1. + Ath2 );

    lth3 = (1. - 3.* Ath3) / ( 1. + Ath3 );
    lph3 =      Aph3       / ( 1. + Ath3 );
    ltp3 =      Atp3       / ( 1. + Ath3 );

    wPPD_sum += wPPD;

    lth1_avg   += wPPD * lth1;
    lth1_sqavg += wPPD * lth1*lth1;
    lth2_avg   += wPPD * lth2;
    lth2_sqavg += wPPD * lth2*lth2;
    lth3_avg   += wPPD * lth3;
    lth3_sqavg += wPPD * lth3*lth3;

    lph1_avg   += wPPD * lph1;
    lph1_sqavg += wPPD * lph1*lph1;
    lph2_avg   += wPPD * lph2;
    lph2_sqavg += wPPD * lph2*lph2;
    lph3_avg   += wPPD * lph3;
    lph3_sqavg += wPPD * lph3*lph3;

    ltp1_avg   += wPPD * ltp1;
    ltp1_sqavg += wPPD * ltp1*ltp1;
    ltp2_avg   += wPPD * ltp2;
    ltp2_sqavg += wPPD * ltp2*ltp2;
    ltp3_avg   += wPPD * ltp3;
    ltp3_sqavg += wPPD * ltp3*ltp3;

    polPPD->Fill();
  }

  ////////////////////////////////////////////////////////////////////////////////////

  ppdFile->Write();
  ppdFile->Close();

  lth1 = lth1_avg / wPPD_sum;
  dlth1 = sqrt( fabs( lth1_sqavg / wPPD_sum - lth1*lth1 ) );
  lph1 = lph1_avg / wPPD_sum;
  dlph1 = sqrt( fabs( lph1_sqavg / wPPD_sum - lph1*lph1 ) );
  ltp1 = ltp1_avg / wPPD_sum;
  dltp1 = sqrt( fabs( ltp1_sqavg / wPPD_sum - ltp1*ltp1 ) );

  lth2 = lth2_avg / wPPD_sum;
  dlth2 = sqrt( fabs( lth2_sqavg / wPPD_sum - lth2*lth2 ) );
  lph2 = lph2_avg / wPPD_sum;
  dlph2 = sqrt( fabs( lph2_sqavg / wPPD_sum - lph2*lph2 ) );
  ltp2 = ltp2_avg / wPPD_sum;
  dltp2 = sqrt( fabs( ltp2_sqavg / wPPD_sum - ltp2*ltp2 ) );

  lth3 = lth3_avg / wPPD_sum;
  dlth3 = sqrt( fabs( lth3_sqavg / wPPD_sum - lth3*lth3 ) );
  lph3 = lph3_avg / wPPD_sum;
  dlph3 = sqrt( fabs( lph3_sqavg / wPPD_sum - lph3*lph3 ) );
  ltp3 = ltp3_avg / wPPD_sum;
  dltp3 = sqrt( fabs( ltp3_sqavg / wPPD_sum - ltp3*ltp3 ) );

  cout << endl;
  cout << "x = " << x[0] << ":" << endl;
  cout << "   lth = " << lth1 << " +/- " << dlth1 << endl;
  cout << "   lph = " << lph1 << " +/- " << dlph1 << endl;
  cout << "   ltp = " << ltp1 << " +/- " << dltp1 << endl;
  cout << endl;
  cout << "x = " << x[1] << ":" << endl;
  cout << "   lth = " << lth2 << " +/- " << dlth2 << endl;
  cout << "   lph = " << lph2 << " +/- " << dlph2 << endl;
  cout << "   ltp = " << ltp2 << " +/- " << dltp2 << endl;
  cout << endl;
  cout << "x = " << x[2] << ":" << endl;
  cout << "   lth = " << lth3 << " +/- " << dlth3 << endl;
  cout << "   lph = " << lph3 << " +/- " << dlph3 << endl;
  cout << "   ltp = " << ltp3 << " +/- " << dltp3 << endl;

  timer->Stop();
  cout << endl << "CPU time:  " << int(timer->CpuTime()) << " s" << endl;
  cout         << "Real time: " << int(timer->RealTime()) << " s" << endl << endl;

} // END

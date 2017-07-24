#include "Riostream.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"
#include "TMath.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TFile.h"

// plotting ranges

const double lth_min = -1.;
const double lth_max =  1.;

const double lph_min = -0.3;
const double lph_max =  0.3;

const double ltp_min = -0.1;
const double ltp_max =  0.1;

const double ltd_min = -1.;
const double ltd_max =  1.;

const int nquantiles = 7;
//                                       ______________ 99.7% _________________
//                                      |       ________  95% __________      |
//                                      |      |                        |     |
const double qprobs[nquantiles] = { 0.0015, 0.025, 0.16, 0.50, 0.84, 0.975, 0.9985 };
//                                                   |          |
//                                                   ---- 68% ---

double A_vs_x( double x, double x1, double x2, double x3, double A1, double A2, double A3 )
{
  return   A1*(x-x2)*(x-x3)/((x1-x2)*(x1-x3))
    + A2*(x-x1)*(x-x3)/((x2-x1)*(x2-x3))
    + A3*(x-x1)*(x-x2)/((x3-x1)*(x3-x2));
}

void polPlot(){

  gROOT->Reset();
  gErrorIgnoreLevel = kWarning;

  TFile* PPDfile = new TFile("polPPD.root");

  TH1D* xi_values = (TH1D*)PPDfile->Get("xi_values");

  double* x   = new double[3];

  for (int i = 0; i < 3; i++) {

    x[i] = xi_values->GetBinContent( i+1 );
  }

  double x_definition = xi_values->GetBinContent( 4 ); // check if the variable is pT or Nch

  string xLabel(x_definition>0?"pT [GeV]":"N_{ch}");

  double xmin = x[0] - 0.3*(x[1]-x[0]); if (xmin < 0) xmin = 0.;
  double xmax = x[2] + 0.3*(x[2]-x[1]);

  TTree* polPPD = (TTree*)PPDfile->Get("polPPD");

  double Ath1;      polPPD->SetBranchAddress( "Ath1",   &Ath1 );
  double Ath2;      polPPD->SetBranchAddress( "Ath2",   &Ath2 );
  double Ath3;      polPPD->SetBranchAddress( "Ath3",   &Ath3 );

  double Aph1;      polPPD->SetBranchAddress( "Aph1",   &Aph1 );
  double Aph2;      polPPD->SetBranchAddress( "Aph2",   &Aph2 );
  double Aph3;      polPPD->SetBranchAddress( "Aph3",   &Aph3 );

  double Atp1;      polPPD->SetBranchAddress( "Atp1",   &Atp1 );
  double Atp2;      polPPD->SetBranchAddress( "Atp2",   &Atp2 );
  double Atp3;      polPPD->SetBranchAddress( "Atp3",   &Atp3 );

  double wPPD;      polPPD->SetBranchAddress( "wPPD",   &wPPD );

  int numEvts = int( polPPD->GetEntries() );

  const int nPPDslices = 100;
  double xstep = (xmax-xmin)/double(nPPDslices);

  const int nbinPPD = 10000;

  TH1D** PPDlth = new TH1D*[nPPDslices];
  TH1D** PPDlph = new TH1D*[nPPDslices];
  TH1D** PPDltp = new TH1D*[nPPDslices];
  TH1D** PPDltd = new TH1D*[nPPDslices];

  for ( int i = 0; i < nPPDslices; i++ ) {
    PPDlth[i] = new TH1D("","",nbinPPD,-1.,1.);
    PPDlph[i] = new TH1D("","",nbinPPD,-1.,1.);
    PPDltp[i] = new TH1D("","",nbinPPD,-0.71,0.71);
    PPDltd[i] = new TH1D("","",nbinPPD,-1.,1.);
  }

  TProfile* lth_vs_x = new TProfile("lth_vs_x","lth_vs_x",100,xmin,xmax,"s");
  TProfile* lph_vs_x = new TProfile("lph_vs_x","lph_vs_x",100,xmin,xmax,"s");
  TProfile* ltp_vs_x = new TProfile("ltp_vs_x","ltp_vs_x",100,xmin,xmax,"s");
  TProfile* ltd_vs_x = new TProfile("ltd_vs_x","ltd_vs_x",100,xmin,xmax,"s");


  ///////////////////////////////////////////////////////////////
  // loop over events in the input ntuple
  ///////////////////////////////////////////////////////////////

  for ( int i = 0; i < numEvts; i++ ) {

    polPPD->GetEvent( i );

    for(int j_x = 0; j_x < nPPDslices; j_x++){

      double xval = xmin + xstep * (j_x+0.5);

      double Ath_x = A_vs_x( xval, x[0], x[1], x[2], Ath1, Ath2, Ath3 );
      double Aph_x = A_vs_x( xval, x[0], x[1], x[2], Aph1, Aph2, Aph3 );
      double Atp_x = A_vs_x( xval, x[0], x[1], x[2], Atp1, Atp2, Atp3 );

      double lth = (1. - 3.*Ath_x)/(1. + Ath_x);
      double lph = Aph_x/(1. + Ath_x);
      double ltp = Atp_x/(1. + Ath_x);
      double ltd = (1. - 3.*(Ath_x-Aph_x))/(1. + (Ath_x-Aph_x));

      lth_vs_x->Fill(xval, lth, wPPD);
      lph_vs_x->Fill(xval, lph, wPPD);
      ltp_vs_x->Fill(xval, ltp, wPPD);
      ltd_vs_x->Fill(xval, ltd, wPPD);

      PPDlth[j_x]->Fill(lth, wPPD);
      PPDlph[j_x]->Fill(lph, wPPD);
      PPDltp[j_x]->Fill(ltp, wPPD);
      PPDltd[j_x]->Fill(ltd, wPPD);
    }
  } // end ntuple event loop
  ///////////////////////////////////////////////////////////////

  double* ex = new double[nPPDslices];
  memset(ex, 0., nPPDslices*sizeof(double) );
  double* xx      = new double[nPPDslices];

  double* lth     = new double[nPPDslices];
  double* dlth1up = new double[nPPDslices];
  double* dlth1dw = new double[nPPDslices];
  double* dlth2up = new double[nPPDslices];
  double* dlth2dw = new double[nPPDslices];
  double* dlth3up = new double[nPPDslices];
  double* dlth3dw = new double[nPPDslices];

  double* lph     = new double[nPPDslices];
  double* dlph1up = new double[nPPDslices];
  double* dlph1dw = new double[nPPDslices];
  double* dlph2up = new double[nPPDslices];
  double* dlph2dw = new double[nPPDslices];
  double* dlph3up = new double[nPPDslices];
  double* dlph3dw = new double[nPPDslices];

  double* ltp     = new double[nPPDslices];
  double* dltp1up = new double[nPPDslices];
  double* dltp1dw = new double[nPPDslices];
  double* dltp2up = new double[nPPDslices];
  double* dltp2dw = new double[nPPDslices];
  double* dltp3up = new double[nPPDslices];
  double* dltp3dw = new double[nPPDslices];

  double* ltd     = new double[nPPDslices];
  double* dltd1up = new double[nPPDslices];
  double* dltd1dw = new double[nPPDslices];
  double* dltd2up = new double[nPPDslices];
  double* dltd2dw = new double[nPPDslices];
  double* dltd3up = new double[nPPDslices];
  double* dltd3dw = new double[nPPDslices];

  double* quantile_position = new double[nquantiles];

  for(int j_x = 0; j_x < nPPDslices; j_x++){

    xx[j_x] = xmin + xstep * (j_x+0.5);

    PPDlth[j_x]->GetQuantiles(nquantiles, quantile_position, qprobs);
    lth[j_x] = quantile_position[3]; // median
    dlth1up[j_x] = quantile_position[4] - quantile_position[3];
    dlth1dw[j_x] = quantile_position[3] - quantile_position[2];
    dlth2up[j_x] = quantile_position[5] - quantile_position[3];
    dlth2dw[j_x] = quantile_position[3] - quantile_position[1];
    dlth3up[j_x] = quantile_position[6] - quantile_position[3];
    dlth3dw[j_x] = quantile_position[3] - quantile_position[0];

    PPDlph[j_x]->GetQuantiles(nquantiles, quantile_position, qprobs);
    lph[j_x] = quantile_position[3]; // median
    dlph1up[j_x] = quantile_position[4] - quantile_position[3];
    dlph1dw[j_x] = quantile_position[3] - quantile_position[2];
    dlph2up[j_x] = quantile_position[5] - quantile_position[3];
    dlph2dw[j_x] = quantile_position[3] - quantile_position[1];
    dlph3up[j_x] = quantile_position[6] - quantile_position[3];
    dlph3dw[j_x] = quantile_position[3] - quantile_position[0];

    PPDltp[j_x]->GetQuantiles(nquantiles, quantile_position, qprobs);
    ltp[j_x] = quantile_position[3]; // median
    dltp1up[j_x] = quantile_position[4] - quantile_position[3];
    dltp1dw[j_x] = quantile_position[3] - quantile_position[2];
    dltp2up[j_x] = quantile_position[5] - quantile_position[3];
    dltp2dw[j_x] = quantile_position[3] - quantile_position[1];
    dltp3up[j_x] = quantile_position[6] - quantile_position[3];
    dltp3dw[j_x] = quantile_position[3] - quantile_position[0];

    PPDltd[j_x]->GetQuantiles(nquantiles, quantile_position, qprobs);
    ltd[j_x] = quantile_position[3]; // median
    dltd1up[j_x] = quantile_position[4] - quantile_position[3];
    dltd1dw[j_x] = quantile_position[3] - quantile_position[2];
    dltd2up[j_x] = quantile_position[5] - quantile_position[3];
    dltd2dw[j_x] = quantile_position[3] - quantile_position[1];
    dltd3up[j_x] = quantile_position[6] - quantile_position[3];
    dltd3dw[j_x] = quantile_position[3] - quantile_position[0];

  }

  TGraph* lthMed = new TGraph(nPPDslices, xx, lth);
  TGraphAsymmErrors* lth68 = new TGraphAsymmErrors(nPPDslices,xx,lth,ex,ex,dlth1dw,dlth1up);
  TGraphAsymmErrors* lth95 = new TGraphAsymmErrors(nPPDslices,xx,lth,ex,ex,dlth2dw,dlth2up);
  TGraphAsymmErrors* lth99 = new TGraphAsymmErrors(nPPDslices,xx,lth,ex,ex,dlth3dw,dlth3up);

  TGraph* lphMed = new TGraph(nPPDslices, xx, lph);
  TGraphAsymmErrors* lph68 = new TGraphAsymmErrors(nPPDslices,xx,lph,ex,ex,dlph1dw,dlph1up);
  TGraphAsymmErrors* lph95 = new TGraphAsymmErrors(nPPDslices,xx,lph,ex,ex,dlph2dw,dlph2up);
  TGraphAsymmErrors* lph99 = new TGraphAsymmErrors(nPPDslices,xx,lph,ex,ex,dlph3dw,dlph3up);

  TGraph* ltpMed = new TGraph(nPPDslices, xx, ltp);
  TGraphAsymmErrors* ltp68 = new TGraphAsymmErrors(nPPDslices,xx,ltp,ex,ex,dltp1dw,dltp1up);
  TGraphAsymmErrors* ltp95 = new TGraphAsymmErrors(nPPDslices,xx,ltp,ex,ex,dltp2dw,dltp2up);
  TGraphAsymmErrors* ltp99 = new TGraphAsymmErrors(nPPDslices,xx,ltp,ex,ex,dltp3dw,dltp3up);

  TGraph* ltdMed = new TGraph(nPPDslices, xx, ltd);
  TGraphAsymmErrors* ltd68 = new TGraphAsymmErrors(nPPDslices,xx,ltd,ex,ex,dltd1dw,dltd1up);
  TGraphAsymmErrors* ltd95 = new TGraphAsymmErrors(nPPDslices,xx,ltd,ex,ex,dltd2dw,dltd2up);
  TGraphAsymmErrors* ltd99 = new TGraphAsymmErrors(nPPDslices,xx,ltd,ex,ex,dltd3dw,dltd3up);


  TCanvas *c1;
  TH1F *hframe;

  // lambdas, avg +- stddev

  c1 = new TCanvas("c1","",200,10,700,500);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.07);
  c1->SetTopMargin(0.075);
  c1->SetBottomMargin(0.125);
  c1->SetFrameBorderMode(0);

  hframe = gPad->DrawFrame( xmin, -1., xmax, 1. );
  hframe->GetXaxis()->SetTitle(xLabel.c_str());
  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetTitleOffset(1.08);
  hframe->GetYaxis()->SetTitle("#lambda_{#vartheta}, #lambda_{#varphi}, #lambda_{#vartheta#varphi}, #tilde{#lambda}");
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.16);
  hframe->Draw(" ");

  lth_vs_x->SetTitle("");
  lth_vs_x->SetFillColor(kRed);
  lth_vs_x->SetLineStyle(1);
  lth_vs_x->SetLineWidth(1.);
  lth_vs_x->Draw( "E3 same" );

  lph_vs_x->SetTitle("");
  lph_vs_x->SetFillColor(kBlue);
  lph_vs_x->SetLineStyle(1);
  lph_vs_x->SetLineWidth(1.);
  lph_vs_x->Draw( "E3 same" );

  ltp_vs_x->SetTitle("");
  ltp_vs_x->SetFillColor(kMagenta);
  ltp_vs_x->SetLineStyle(1);
  ltp_vs_x->SetLineWidth(1.);
  ltp_vs_x->Draw( "E3 same" );

  ltd_vs_x->SetTitle("");
  ltd_vs_x->SetFillColor(kGreen+2);
  ltd_vs_x->SetLineStyle(1);
  ltd_vs_x->SetLineWidth(1.);
  ltd_vs_x->Draw( "E3 same" );

  c1->Print( "lambdas_vs_x_avg-sig.pdf" );



  // lambdatheta, with 3 bands

  hframe = gPad->DrawFrame( xmin, lth_min, xmax, lth_max );
  hframe->GetXaxis()->SetTitle(xLabel.c_str());
  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetTitleOffset(1.08);
  hframe->GetYaxis()->SetTitle("#lambda_{#vartheta}");
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.16);
  hframe->Draw(" ");

  lth99->SetTitle("");
  lth99->SetFillColor(kRed-9);
  lth99->Draw( "E3 same" );

  lth95->SetTitle("");
  lth95->SetFillColor(kRed-7);
  lth95->Draw( "E3 same" );

  lth68->SetTitle("");
  lth68->SetFillColor(kRed);
  lth68->Draw( "E3 same" );

  lthMed->SetTitle("");
  lthMed->SetLineStyle(1);
  lthMed->SetLineWidth(1.);
  lthMed->SetLineColor(kRed+2);
  lthMed->Draw( "L same" );

  c1->Print( "lambdatheta_vs_x.pdf" );



  // lambdaphi, with 3 bands

  hframe = gPad->DrawFrame( xmin, lph_min, xmax, lph_max );
  hframe->GetXaxis()->SetTitle(xLabel.c_str());
  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetTitleOffset(1.08);
  hframe->GetYaxis()->SetTitle("#lambda_{#varphi}");
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.16);
  hframe->Draw(" ");

  lph99->SetTitle("");
  lph99->SetFillColor(kBlue-9);
  lph99->Draw( "E3 same" );

  lph95->SetTitle("");
  lph95->SetFillColor(kBlue-7);
  lph95->Draw( "E3 same" );

  lph68->SetTitle("");
  lph68->SetFillColor(kBlue);
  lph68->Draw( "E3 same" );

  lphMed->SetTitle("");
  lphMed->SetLineStyle(1);
  lphMed->SetLineWidth(1.);
  lphMed->SetLineColor(kBlue+2);
  lphMed->Draw( "L same" );

  c1->Print( "lambdaphi_vs_x.pdf" );


  // lambdathetaphi, with 3 bands

  hframe = gPad->DrawFrame( xmin, ltp_min, xmax, ltp_max );
  hframe->GetXaxis()->SetTitle(xLabel.c_str());
  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetTitleOffset(1.08);
  hframe->GetYaxis()->SetTitle("#lambda_{#vartheta#varphi}");
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.16);
  hframe->Draw(" ");

  ltp99->SetTitle("");
  ltp99->SetFillColor(kMagenta-9);
  ltp99->Draw( "E3 same" );

  ltp95->SetTitle("");
  ltp95->SetFillColor(kMagenta-7);
  ltp95->Draw( "E3 same" );

  ltp68->SetTitle("");
  ltp68->SetFillColor(kMagenta);
  ltp68->Draw( "E3 same" );

  ltpMed->SetTitle("");
  ltpMed->SetLineStyle(1);
  ltpMed->SetLineWidth(1.);
  ltpMed->SetLineColor(kMagenta+2);
  ltpMed->Draw( "L same" );

  c1->Print( "lambdathetaphi_vs_x.pdf" );


  // lambdatilde, with 3 bands

  hframe = gPad->DrawFrame( xmin, ltd_min, xmax, ltd_max );
  hframe->GetXaxis()->SetTitle(xLabel.c_str());
  hframe->GetXaxis()->SetLabelOffset(0.007);
  hframe->GetXaxis()->SetLabelSize(0.04);
  hframe->GetXaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetTitleOffset(1.08);
  hframe->GetYaxis()->SetTitle("#tilde{#lambda}");
  hframe->GetYaxis()->SetLabelSize(0.04);
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetYaxis()->SetTitleOffset(1.16);
  hframe->Draw(" ");

  ltd99->SetTitle("");
  ltd99->SetFillColor(kGreen);
  ltd99->Draw( "E3 same" );

  ltd95->SetTitle("");
  ltd95->SetFillColor(kGreen+1);
  ltd95->Draw( "E3 same" );

  ltd68->SetTitle("");
  ltd68->SetFillColor(kGreen+2);
  ltd68->Draw( "E3 same" );

  ltdMed->SetTitle("");
  ltdMed->SetLineStyle(1);
  ltdMed->SetLineWidth(1.);
  ltdMed->SetLineColor(kGreen+3);
  ltdMed->Draw( "L same" );

  c1->Print( "lambdatilde_vs_x.pdf" );

  /////// end
}

const char outfile[256] = "genDataRef.root";

const double pTdilepton_min =  10.;
const double pTdilepton_max =  70.;
const double rapdilepton_min = 0.0; // put here only interval in the positive-rapidity region!
const double rapdilepton_max = 1.0; // e.g.: [0.5, 1.0] to generate in [-1.0, -0.5] and [0.5, 1.0]
const double Nch_min =  0.;
const double Nch_max =  90.;
const double mass_signal_peak  =  3.10;
const double mass_signal_sigma =  0.05;
const double n_sigmas_signal = 4.;
double mass_min = mass_signal_peak - n_sigmas_signal*mass_signal_sigma;
double mass_max = mass_signal_peak + n_sigmas_signal*mass_signal_sigma;

// sideband definition
double lftSideEnd = mass_signal_peak - 3.*mass_signal_sigma;
double rgtSideStt = mass_signal_peak + 3.*mass_signal_sigma;


const long n_events = 6000000; // number of signal + bkg events generated in the entire lepton momentum space

// background fraction
const double f_BG = 0.40; // do not set smaller than 0.001

// pT distribution of the signal
double func_pT_gen(double* x, double* par)
{
  const double mass = 3.1;
  double pTovM = x[0] / mass;
  const double beta = 3.45;
  const double gamma = 0.8;
  return pTovM * pow( 1. + 1./(beta - 2.) * pTovM*pTovM / gamma, -beta  );
}

// pT distribution of the BG
double func_pT_gen_BG(double* x, double* par)
{
  const double mass = 3.1;
  double pTovM = x[0] / mass;
  const double beta = 6.0;
  const double gamma = 0.8;
  return pTovM * pow( 1. + 1./(beta - 2.) * pTovM*pTovM / gamma, -beta  );
}

// Nch distribution of the signal
double func_Nch_gen(double* x, double* par)
{
  double Nch = x[0];
  const double beta = 3.15;
  const double gamma = 900.;
  return Nch * pow( 1. + 1./(beta - 2.) * Nch*Nch / gamma, -beta  );
}

// Nch distribution of the BG
double func_Nch_gen_BG(double* x, double* par)
{
  double Nch = x[0];
  const double beta = 2.7;
  const double gamma = 900.;
  return Nch * pow( 1. + 1./(beta - 2.) * Nch*Nch / gamma, -beta  );
}

// rapidity distribution (of both signal and background)
double func_rap_gen(double* x, double* par)
{
  return   1.;
}

// generated polarization for signal, as a possible function of pT and/or Nch
inline double lambda_theta_sig(double pT, double Nch)
{
  return 0.;
}
inline double lambda_phi_sig(double pT, double Nch)
{
  return 0.;
}
inline double lambda_thetaphi_sig(double pT, double Nch)
{
  return 0.;
}

// natural frame for signal
const bool HX_is_natural_sig = true;  // put both to false to generate in the CS frame
const bool PX_is_natural_sig = false;


// generated polarization for background, as a possible function of pT and/or Nch
inline double lambda_theta_bkg(double pT, double Nch)
{
  return 1.0;
}
inline double lambda_phi_bkg(double pT, double Nch)
{
  return 1.0;
}
inline double lambda_thetaphi_bkg(double pT, double Nch)
{
  return 0.0;
}

// natural frame for background
const bool HX_is_natural_bkg = true;
const bool PX_is_natural_bkg = false;


// acceptance

double leptoncut = 4.0;
//double leptoncut = 0.0;

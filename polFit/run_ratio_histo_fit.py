#!/usr/bin/env python

import argparse
import json


import ROOT as r

# since root associates all histograms to the file currently pointed to by gDirectory,
# TH1s that are created in a function, where the file is only in scope for the function,
# get destroyed once the TFile is closed, resulting in a None, when the histo is
# accessed later.
# To avoid having TFiles going out of scope to early, collect them in this global list.
_open_files = []


def create_histo(filename, weight, treename='genData', name=''):
    """
    Create the costh_PX histogram from the datafile using the passed treename
    """
    from utils.jupyter_helpers import drawVarToHist, setHistOpts
    from utils.miscHelpers import createRandomString

    dataf = r.TFile.Open(filename)
    global _open_files
    _open_files.append(dataf)
    tree = dataf.Get(treename)

    if not name:
        name = createRandomString(16)

    hist = r.TH1D(name, ';cos#theta', 32, -1, 1)
    drawVarToHist(tree, hist, 'costh_PX', '', weight)

    setHistOpts(hist)

    return hist


def get_histos(refn, datan, treen):
    """
    Get the costh_PX histogram from the reference and the data file
    """
    datahist = create_histo(datan, weight='wS', name='h_costh_PX_chic2', treename=treen)
    refhist = create_histo(refn, weight='wS', name='h_costh_PX_chic1', treename=treen)

    return (datahist, refhist)


def set_neg_bins_to_zero(h, verbose=False):
    """Set negative bins of passed histo to zero"""
    neg_bins = [(i, b) for i, b in enumerate(h) if b < 0]
    for nb, cont in neg_bins:
        h.SetBinContent(nb, 0)
        h.SetBinError(nb, 0)

    if verbose:
        print('checked {} for negative bins'.format(h.GetName()))
        for nb, cont in neg_bins:
            print('Set bin {} to 0, content was {}'.format(nb, cont))


def def_fit_func(name='', chic1_limits=True, fix_ref=None, fit_range=[-1,1]):
    """
    Define the fitting function
    """
    if not name:
        name = 'W_costh_ratio'
    f = r.TF1(name,
              '[0] * (1.0 + ([1] + [2]) * x[0]*x[0]) / (1.0 + [1] * x[0]*x[0])',
              fit_range[0], fit_range[1])

    if chic1_limits:
        f.SetParLimits(1, -0.3, 1)

    # don't have to check if chic1_limits is true, this overrides it in any cas
    if fix_ref is not None:
        f.FixParameter(1, fix_ref)

    return f


def divide(datah, refh):
    """
    Clone the datah and divide it by refh
    """
    ratioh = datah.Clone(datah.GetName().replace('chic2', 'ratio')) # tailored to chic2 - chic1 currently
    ratioh.Divide(refh)

    return ratioh


def do_fit(datah, refh, chic1_limits=True, fix_ref=None, fit_range=False):
    """
    Do the fit and return the fit results
    """
    ratio = divide(datah, refh)

    range_to_fit = [-1, 1]
    if fit_range:
        filled_bins = [i for i, b in enumerate(ratio) if b != 0]
        range_to_fit = [ratio.GetBinLowEdge(filled_bins[0]),
                        ratio.GetBinLowEdge(filled_bins[-1]) + ratio.GetBinWidth(filled_bins[-1])]

    fit_func = def_fit_func(chic1_limits=chic1_limits, fix_ref=fix_ref, fit_range=range_to_fit)

    fit_rlt = ratio.Fit(fit_func, 'SR')
    ratio.GetListOfFunctions().Clear()
    if int(fit_rlt) == 0:
        return (fit_rlt, ratio)

    print('Fit returned status {} when fitting the ratio of '
          '{} to {}'.format(int(fit_rlt), datah.GetName(), refh.GetName()))
    return None


def extract_par_from_result(fit_rlt):
    """
    Extract the parameters from the fit results
    """
    return {
        'N': [fit_rlt.Parameter(0), fit_rlt.Error(0)],
        'lth_ref': [fit_rlt.Parameter(1), fit_rlt.Error(1)],
        'delta_lth': [fit_rlt.Parameter(2), fit_rlt.Error(2)],
        'chi2': fit_rlt.Chi2(),
        'ndf': fit_rlt.Ndf()
    }



def run(datafn, reffn, outfn, treen, chic1_limits, fix_ref, fit_range):
    """
    Run the fit for a given data and reference file, store the histograms
    and fit results into a root file and the results in a json file that can be
    read by the analyzer script
    """
    # get the histograms and set negtaive bins to zero
    datah, refh = get_histos(datafn, reffn, treen)

    outfile = r.TFile(outfn, 'recreate')
    outfile.cd()

    datah.Write('_'.join([datah.GetName(), 'raw']))
    refh.Write('_'.join([refh.GetName(), 'raw']))

    set_neg_bins_to_zero(datah, True)
    set_neg_bins_to_zero(refh, True)
    datah.Write()
    refh.Write()

    fit_rlt, ratio = do_fit(datah, refh, chic1_limits, fix_ref, fit_range)
    if fit_rlt is not None:
        fit_rlt.Write('fit_result_costh_ratio')
        ratio.Write()
        with open(outfn.replace('.root', '.json'), 'w') as f:
            json.dump(extract_par_from_result(fit_rlt), f, indent=2)

    outfile.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script that runs a costh ratio fit using histograms')
    parser.add_argument('dataFileName', help='filename of the data tuple')
    parser.add_argument('refFileName', help='filename of the ref tuple')
    parser.add_argument('outFileName', help='output file name')
    parser.add_argument('-t', '--treename', help='treename', dest='tree', default='genData')
    parser.add_argument('-l', '--chic1_limits', action='store_true', default=False, dest='limit',
                        help='Impose a limit on the chic1 polarization (i.e. reference)')
    parser.add_argument('-f', '--fix_ref', type=float, default=None, dest='fix_ref',
                        help='fix ref_lth to this value in fits')
    parser.add_argument('-r', '--fix_range', action='store_true', default=False, dest='fit_range',
                        help='fix the range to the bins that are filled in the histogram.')


    args = parser.parse_args()

    r.gROOT.SetBatch()

    run(args.dataFileName, args.refFileName, args.outFileName,
        args.tree, args.limit, args.fix_ref, args.fit_range)

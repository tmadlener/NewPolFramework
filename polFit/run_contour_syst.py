#!/usr/bin/env python

import pandas as pd
import numpy as np

from utils.ratio_fitting import CosthRatioFit

from run_ratio_histo_fit import get_histos, divide
from run_ratio_fit_new import get_free_norm, error_levels
from runToyMCFits import read_gen_config, get_combinations


def run_one_fit(datafn, reffn, treen, n_bins):
    """Run fit for one data and reference file"""
    import ROOT as r
    norm = get_free_norm(datafn, reffn, treen, n_bins)

    if norm is None:
        return None

    fit_func = r.TF1('fit_func',
                     '[0] * (3 + [1]) / (3 + [1] + [2]) * (1 + ([1] + [2]) * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                     -1, 1)

    datah, refh = get_histos(datafn, reffn, treen, n_bins)
    ratioh = divide(datah, refh)
    ratio_fit = CosthRatioFit(ratioh, fit_func, fix_params=[(0, norm)])

    return ratio_fit.get_results(error_levels)



def run_all_gens(gendir, treename, n_bins):
    """Run all combinations of data and ref files"""
    data_ref_combis = get_combinations(gendir, 'foo', True, ra_lth_ref=[0.5,0.5],
                                       ra_lth_data=[-0.5,-0.5])

    results = []

    for dataf, reff, _ in data_ref_combis:
        fit_res = run_one_fit(dataf, reff, treename, n_bins)
        results.append(fit_res)

    return results

def analyze_one_result(results):
    """Analyze one fit result"""
    def fill_nan(res_d, keys):
        for k in keys:
            res_d[k] = np.nan


    # print(results)
    res_d = {}
    if results is None:
        fill_nan(res_d, ['lth', 'lth_err1D', 'lth_err2D', 'lth_errup1D',
                         'lth_errlow1D', 'lth_errup2D', 'lth_errlow2D',
                         'dlam', 'dlam_err1D', 'dlam_err2D', 'dlam_errup1D',
                         'dlam_errlow1D', 'dlam_errup2D', 'dlam_errlow2D',
                         'N', 'chi2', 'ndf'])
        res_d['valid1D'] = False
        res_d['valid2D'] = False

    else:
        # central values are identical for both results
        params = results[0]['params']
        res_d['N'] = params[0]
        res_d['lth'] = params[1]
        res_d['dlam'] = params[2]

        # errors differ for the two error levels
        errors = results[0]['errors']
        res_d['lth_err1D'] = errors[1]
        res_d['dlam_err1D'] = errors[2]

        errors = results[1]['errors']
        res_d['lth_err2D'] = errors[1]
        res_d['dlam_err2D'] = errors[2]

        # asym errors
        errors = results[0]['low_errors']
        res_d['lth_errlow1D'] = errors[1]
        res_d['dlam_errlow1D'] = errors[2]

        errors = results[0]['up_errors']
        res_d['lth_errup1D'] = errors[1]
        res_d['dlam_errup1D'] = errors[2]

        errors = results[1]['low_errors']
        res_d['lth_errlow2D'] = errors[1]
        res_d['dlam_errlow2D'] = errors[2]

        errors = results[1]['up_errors']
        res_d['lth_errup2D'] = errors[1]
        res_d['dlam_errup2D'] = errors[2]

        res_d['chi2'] = results[0]['chi2']
        res_d['ndf'] = results[0]['ndf']

        res_d['valid1D'] = (results[0]['contour'] is not None)
        res_d['valid2D'] = (results[1]['contour'] is not None)

    return res_d

def analyze_fit_results(fit_results):
    """Cleanup fit_results and put them into a DataFrame"""
    results = []

    for fit_res in fit_results:
        results.append(analyze_one_result(fit_res))

    return pd.DataFrame(results)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='script for running systematic tests on toy mc data')
    parser.add_argument('gendir', help='directory containing the generated events')
    parser.add_argument('-t', '--treename', default='genData',
                        help='treename')
    parser.add_argument('-n', '--nbins', default=32, type=int,
                        help='number of bins to use in costh')

    args = parser.parse_args()

    results = run_all_gens(args.gendir, args.treename, args.nbins)
    df = analyze_fit_results(results)

    df.to_pickle('toymc_contour_results_' + str(args.nbins) + '.pkl')

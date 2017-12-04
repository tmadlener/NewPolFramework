#!/usr/bin/env python

import argparse

import pandas as pd
import ROOT as r
import numpy as np

from scipy.stats import chi2

from utils.ratio_fitting import CosthRatioFit
from utils.plotHelpers import setColor

from run_ratio_histo_fit import get_histos, divide
from runLamRefScan import run_scan

# define error levels for contours
# 1 - corresponds to 68 % confidence intervals on single variable
# ~2.3 - corresponds to 68 % confidence region on bivariate
error_levels = [1, chi2.ppf(chi2.cdf(2,2), 2)]
# error_levels = [1]
line_labels = ['1 #sigma (1D)', '1 #sigma (2D)']

# settings for line drawing
line_settings = [(0, 2), (0, 3), (0, 4), (0, 6)]

# clones of objects to be kept alive for legend
_legend_clones = []

def flat_list_tuple(tup_list):
    """
    From a list of tuples containing two elements, return a list of tuples,
    where the first elements, where the first is the tuple containing all first elemnts
    of the input list and the second contains one of all possible combinations
    of the items in the second elements of the tuple.
    """
    from itertools import product
    fix_funcs, fix_vals = zip(*tup_list) # split apart the fix_funcs and the values

    val_comb_gen = product(*[vals.tolist() for vals in fix_vals])
    return ((fix_funcs, val_comb) for val_comb in val_comb_gen)


def calc_mean_norm(datafn, reffn):
    """calculate the mean normalization"""
    scan_res = run_scan(datafn, reffn, ref_range=[-1,1], n_points=10)
    scan_res = pd.DataFrame(scan_res)
    # get mean normalization and min and max value to check stability
    norm = scan_res.N.apply(lambda x: x[0])
    mean_norm = norm.mean()
    max_norm = norm.max()
    min_norm = norm.min()

    percent = lambda x: (x - mean_norm) / mean_norm * 100

    print('mean,max,min = {}, {}, {}'.format(mean_norm, max_norm, min_norm))
    print('deviations from mean = {} %, {} %'.format(percent(max_norm), percent(min_norm)))

    return mean_norm


def scan_chi2_params(hist, fit_func, scan_params):
    results = []
    par_names = {0: 'N', 1: 'lth_ref', 2: 'delta_lth'}

    for fix_funcs, val_comb in flat_list_tuple(scan_params):
        pars = {}
        for (i, func) in enumerate(fix_funcs):
            func(fit_func, val_comb[i])
            pars[par_names[i]] = val_comb[i]

        pars['chi2'] = hist.Chisquare(fit_func)
        results.append(pars)

    return pd.DataFrame(results)


def get_bin_centers(x_min, x_max, n_bins):
    """Get the values for the evaluation of the parameters in the bin centers"""
    bin_bounds = np.linspace(x_min, x_max, n_bins + 1)
    return 0.5 * (bin_bounds[1:] + bin_bounds[:-1])


def get_best_fit(fit_rlt):
    """Get the best fit and make a TGraph with one point of it"""
    p1_best = np.array([fit_rlt['params'][1]])
    p2_best = np.array([fit_rlt['params'][2]])

    return r.TGraph(1, p1_best, p2_best)


def add_clone_to_leg(leg, elem, label, opt):
    """Add a clone to the legend, with black line and markers for visiblity"""
    global _legend_clones
    e_clone = elem.Clone()
    _legend_clones.append(e_clone)

    setColor(e_clone,1)
    leg.AddEntry(e_clone, label, opt)


def make_paed_plot(ratioh, func, err_lvls, fit_rlts, plotname,
                   n_grid=100, x_ran=[-1, 1], y_ran=[-1, 1]):
    """Make a paedagoical plot"""
    plotHist = r.TH2D('h', ';#lambda_{ref};#Delta_{#lambda};#chi^{2} - #chi^{2}_{min}',
                      n_grid, x_ran[0], x_ran[1], n_grid, y_ran[0], y_ran[1])
    plotHist.SetStats(0)
    plotHist.SetTitleSize(0.04, 'XYZ')
    plotHist.SetTitleOffset(0.9, 'X')
    plotHist.SetTitleOffset(0.95, 'YZ')


    par_scan_settings = (
        (lambda f, x: f.FixParameter(0, x), np.array([fit_rlts[0]['params'][0]])),
        (lambda f, x: f.FixParameter(1, x), get_bin_centers(x_ran[0], x_ran[1], n_grid)),
        (lambda f, x: f.FixParameter(2, x), get_bin_centers(y_ran[0], y_ran[1], n_grid))
    )

    min_chi2 = fit_rlts[0]['chi2']
    best_fit = get_best_fit(fit_rlts[0])

    scan_results = scan_chi2_params(ratioh, func, par_scan_settings)
    for _,vals in scan_results.iterrows():
        if vals.chi2 - min_chi2 > 25:
            continue
        i_bin = plotHist.FindBin(vals.lth_ref, vals.delta_lth)
        plotHist.SetBinContent(i_bin, vals.chi2 - min_chi2)

    # it is possible that the scanning actually finds a (slightly) better minimum
    # than the minimization from the fit, so to have plots with white space, where
    # the histogram has not been filled, set all bins with content zero to -1
    # and adjust the color range to include slightly negative values
    for i,b in enumerate(plotHist):
        if b == 0:
            plotHist.SetBinContent(i, -1)

    plotHist.GetZaxis().SetRangeUser(-0.01, 25.0)

    c = r.TCanvas('c', 'c', 800, 800)
    pad = c.cd()
    c.SetGrid()
    pad.SetRightMargin(0.15)
    pad.SetLeftMargin(0.12)

    leg = r.TLegend(0.12, 0.1, 0.25, 0.25)

    plotHist.Draw('colz')
    for i, rlts in enumerate(fit_rlts):
        if rlts['contour'] is None: # catch cases of failing contour
            continue
        rlts['contour'].SetLineColor(line_settings[i][0])
        rlts['contour'].SetLineStyle(line_settings[i][1])
        rlts['contour'].Draw('sameL')
        add_clone_to_leg(leg, rlts['contour'], line_labels[i], 'L')

        best_fit.SetMarkerStyle(22)
    best_fit.SetMarkerColor(0)
    best_fit.Draw('sameP')
    add_clone_to_leg(leg, best_fit, 'best fit', 'P')

    leg.Draw()
    c.Draw()
    c.SaveAs(plotname)


def get_plotting_range(results):
    """Determine the range for the paedagoical plot"""
    x_cent, y_cent = results['params'][1:] # hardcoded to three params
    x_sym, y_sym = results['errors'][1:]
    x_low, y_low = results['low_errors'][1:]
    x_up, y_up = results['up_errors'][1:]

    x_upper_bounds = np.array([x_cent + x_sym, x_cent + x_up])
    y_upper_bounds = np.array([y_cent + y_sym, y_cent + y_up])
    x_lower_bounds = np.array([x_cent - x_sym, x_cent + x_low])
    y_lower_bounds = np.array([y_cent - y_sym, y_cent + y_low])

    return ([np.min(x_lower_bounds), np.max(x_upper_bounds)],
            [np.min(y_lower_bounds), np.max(y_upper_bounds)])


def main(datafn, reffn, outbase):
    # first run a scan to check how stable the normalization is
    mean_norm = calc_mean_norm(datafn, reffn)
    fit_func = r.TF1('fit_func',
                     '[0] * (3 + [1]) / (3 + [1] + [2]) * (1 + ([1] + [2]) * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                     -1, 1)

    datah, refh = get_histos(datafn, reffn, 'chic_tuple')
    ratioh = divide(datah, refh)
    ratio_fit = CosthRatioFit(ratioh, fit_func, fix_params=[(0, mean_norm)])


    fit_results = ratio_fit.get_results(error_levels)

    x_ran, y_ran = get_plotting_range(fit_results[-1])

    make_paed_plot(ratioh, fit_func, error_levels, fit_results,
                   '.'.join([outbase, 'pdf']),
                   n_grid=100, x_ran=x_ran, y_ran=y_ran)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for running a costh ratio fit '
                                     'returning contours')
    parser.add_argument('data_file_name', help='data file name')
    parser.add_argument('ref_file_name', help='ref file name')
    parser.add_argument('-o', '--output_base', help='output base dir',
                        default='fit_output')

    args = parser.parse_args()

    main(args.data_file_name, args.ref_file_name, args.output_base)

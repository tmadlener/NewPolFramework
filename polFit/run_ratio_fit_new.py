#!/usr/bin/env python

import argparse

import pandas as pd
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
import numpy as np

from scipy.stats import chi2
from scipy.spatial import ConvexHull

from utils.ratio_fitting import CosthRatioFit
from utils.plotHelpers import setColor

from run_ratio_histo_fit import get_histos, divide, run, set_bins_to_zero
from runLamRefScan import run_scan

# define error levels for contours
# 1 - corresponds to 68 % confidence intervals on single variable
# ~2.3 - corresponds to 68 % confidence region on bivariate
error_levels = [1, chi2.ppf(chi2.cdf(1,1), 2)]
# error_levels = [1]
line_labels = ['1 #sigma (1D)', '1 #sigma (2D)']

# settings for line drawing
line_settings = [(0, 2), (0, 3), (0, 4), (0, 6)]

# clones of objects to be kept alive for legend
_legend_clones = []

# clones of scan graphs
_scan_graphs = []

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


def calc_mean_norm(datafn, reffn, treen, n_bins):
    """calculate the mean normalization"""
    scan_res = run_scan(datafn, reffn, treen, ref_range=[-1,1], n_points=10,
                        n_bins=n_bins)
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


def get_free_norm(datafn, reffn, treen, n_bins):
    """calculate the normalization from a fit with all free parameters"""
    fit_results = run(datafn, reffn, '', treen, chic1_limits=False,
                      fix_ref=None, fit_range=False, save=False, fix_norm=False,
                      n_bins=n_bins)

    if fit_results is not None:
        print('Normalization from free fit = {} +/- {}'.format(fit_results['N'][0], fit_results['N'][1]))
        return fit_results['N'][0]

    print('FITERROR: Cannot determine normalization from free fit for inputs:\n'
          'datafile: {}\nreffile: {}'.format(datafn, reffn))
    return None


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

    p1_err_low = -np.array([fit_rlt['low_errors'][1]])
    p2_err_low = -np.array([fit_rlt['low_errors'][2]])

    p1_err_up = np.array([fit_rlt['up_errors'][1]])
    p2_err_up = np.array([fit_rlt['up_errors'][2]])

    # return r.TGraph(1, p1_best, p2_best)
    return r.TGraphAsymmErrors(1, p1_best, p2_best, p1_err_low, p1_err_up,
                               p2_err_low, p2_err_up)


def add_clone_to_leg(leg, elem, label, opt):
    """Add a clone to the legend, with black line and markers for visiblity"""
    global _legend_clones
    e_clone = elem.Clone()
    _legend_clones.append(e_clone)

    setColor(e_clone,1)
    leg.AddEntry(e_clone, label, opt)


def get_contour_scan(scan_results, min_chi2, err_level):
    """
    Get the contour from the scan by simply calculating the convex hull of
    all points with chi2 values smaller than min_chi2 + err_level

    Return a dict with the same structure as the CosthRatioFit.get_results
    """
    # min_chi2 = scan_results.chi2.min()
    cont_sel = np.abs(scan_results.chi2 - min_chi2) < err_level

    sel_points = np.array([scan_results[cont_sel].lth_ref,
                           scan_results[cont_sel].delta_lth]).T

    hull = ConvexHull(sel_points)

    n_points = len(hull.vertices)
    n_x = sel_points[hull.vertices, 0]
    n_y = sel_points[hull.vertices, 1]

    # calculate the rest of the information as well
    # i.e. min and asymm uncertainties
    chi2_min_idx = scan_results.chi2.argmin()
    lth_ref, delta_lth = scan_results.iloc[chi2_min_idx][['lth_ref', 'delta_lth']]

    lth_min, delta_min = hull.min_bound
    lth_max, delta_max = hull.max_bound

    fit_res = {}
    fit_res['chi2'] = min_chi2
    fit_res['params'] = [scan_results.iloc[chi2_min_idx].N, lth_ref, delta_lth]
    fit_res['up_errors'] = [0, lth_max - lth_ref, delta_max - delta_lth]
    fit_res['low_errors'] = [0, lth_min - lth_ref, delta_min - delta_lth]
    fit_res['ndf'] = -1
    fit_res['errors'] = [0, 0, 0]
    fit_res['err_level'] = err_level
    fit_res['contour'] = r.TGraph(n_points, n_x, n_y)

    return fit_res


def make_paed_plot(ratioh, func, err_lvls, fit_rlts, plotname,
                   n_grid=100, x_ran=[-1, 1], y_ran=[-1, 1]):
    """Make a paedagoical plot"""
    plotHist = r.TH2D('h', ';#lambda_{ref};#Delta_{#lambda};#chi^{2} - #chi^{2}_{min}',
                      n_grid, x_ran[0], x_ran[1], n_grid, y_ran[0], y_ran[1])
    plotHist.SetStats(0)
    plotHist.SetTitleSize(0.04, 'XYZ')
    plotHist.SetTitleOffset(0.9, 'X')
    plotHist.SetTitleOffset(0.95, 'Z')
    plotHist.SetTitleOffset(1.4, 'Y')

    par_scan_settings = (
        (lambda f, x: f.FixParameter(0, x), np.array([fit_rlts[0]['params'][0]])),
        (lambda f, x: f.FixParameter(1, x), get_bin_centers(x_ran[0], x_ran[1], n_grid)),
        (lambda f, x: f.FixParameter(2, x), get_bin_centers(y_ran[0], y_ran[1], n_grid))
    )

    min_chi2 = fit_rlts[0]['chi2']
    best_fit = get_best_fit(fit_rlts[0])

    scan_results = scan_chi2_params(ratioh, func, par_scan_settings)
    # it is possible that the scanning actually finds a (slightly) better minimum
    # than the minimization from the fit. If that is the case use the new minimum chi2
    # value for the plot and also indicate the better result in the plot
    scan_smaller_minimization = np.any(scan_results.chi2 < min_chi2)
    if scan_smaller_minimization:
        min_chi2 = scan_results.chi2.min()

    for _,vals in scan_results.iterrows():
        if vals.chi2 - min_chi2 > 25:
            continue
        i_bin = plotHist.FindBin(vals.lth_ref, vals.delta_lth)
        plotHist.SetBinContent(i_bin, vals.chi2 - min_chi2)

    # not all bins get filled above. To have them appear white, change their
    # value to -1
    for i, b in enumerate(plotHist):
        if b == 0:
            plotHist.SetBinContent(i, -1)

    plotHist.GetZaxis().SetRangeUser(-0.01, 25.0)

    c = r.TCanvas('c', 'c', 800, 800)
    pad = c.cd()
    c.SetGrid()
    pad.SetRightMargin(0.15)
    pad.SetLeftMargin(0.12)

    leg = r.TLegend(0.12, 0.1, 0.25, 0.25)

    scan_fit_results = []

    plotHist.Draw('colz')
    for i, rlts in enumerate(fit_rlts):
        scan_cont = False
        # catch cases of failing contour or new minimum in scanning
        if rlts['contour'] is None or scan_smaller_minimization:
            scan_fit_results.append(get_contour_scan(scan_results, min_chi2, rlts['err_level']))
            contour = scan_fit_results[-1]['contour']
            scan_cont = True
        else:
            contour = rlts['contour']
        contour.SetLineColor(line_settings[i][0])
        contour.SetLineStyle(line_settings[i][1])
        contour.Draw('sameL')

        label_add = ' scan' if scan_cont else ''
        add_clone_to_leg(leg, contour, line_labels[i] + label_add, 'L')

    # if scan_fit_results and scan_smaller_minimization:
    if scan_smaller_minimization:
        scan_best_fit = get_best_fit(scan_fit_results[0])
        scan_best_fit.SetMarkerStyle(23)
        scan_best_fit.SetMarkerColor(0)
        scan_best_fit.Draw('sameP')
        add_clone_to_leg(leg, scan_best_fit, 'best fit (scan)', 'P')

    best_fit.SetMarkerStyle(22)
    best_fit.SetMarkerColor(0)
    best_fit.Draw('samePX') # 'X' for not drawing error bars
    add_clone_to_leg(leg, best_fit, 'best fit', 'P')

    leg.Draw()
    c.Draw()
    c.SaveAs(plotname)

    return scan_fit_results


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


def print_result(fit_result):
    """Print one fit result in a (relatively) easyly readable way"""
    print('============================== RESULT ==============================')
    print('error level = {:.2f}'.format(fit_result['err_level']))
    print('chi2 / ndf = {:.2f}/{}'.format(fit_result['chi2'], fit_result['ndf']))
    print('N = {:.4f} +/- {:.4f}'.format(fit_result['params'][0], fit_result['errors'][0]))
    print('lambda_ref = {:.4f} +/- {:.4f} (+{:.4f}, -{:.4f})'.format(fit_result['params'][1],
                                                                     fit_result['errors'][1],
                                                                     fit_result['up_errors'][1],
                                                                     -fit_result['low_errors'][1]))
    print('delta_lam = {:.4f} +/- {:.4f} (+{:.4f}, -{:.4f})'.format(fit_result['params'][2],
                                                                    fit_result['errors'][2],
                                                                    fit_result['up_errors'][2],
                                                                    -fit_result['low_errors'][2]))



def make_fit_plot(ratioh, fit_res, fit_func, outbase):
    """Make a plot showing the ratio and the fit result"""
    from utils.TH_utils import calc_pulls
    from utils.plotHelpers import _defaultColors, mkplot, plotOnCanvas

    colors = _defaultColors()

    fit_func.SetParameters(*fit_res['params'])

    can = r.TCanvas('fit_plot_canvas', '', 800, 800)
    can.cd()
    hpad = r.TPad('hpad', 'hpad', 0, 0.25, 1, 1)
    r.SetOwnership(hpad, False)
    hpad.Draw()
    hpad.SetGrid()

    ratioh.GetYaxis().SetRangeUser(0, 1)
    ratioh.SetYTitle('N_{#chi_{c2}} / N_{#chi_{c1}}')
    plotOnCanvas(hpad, [ratioh, fit_func], drawOpt='E1')

    can.cd()
    ppad = r.TPad('ppad', 'ppad', 0, 0, 1, 0.25)
    r.SetOwnership(ppad, False)
    ppad.Draw()

    ppad.cd()
    ppad.SetGrid()
    plothist = ppad.DrawFrame(-1, -5, 1, 5)
    # plothist.SetXTitle('cos#theta')
    plothist.SetYTitle('(fit - data)/#sigma_{data}')
    plothist.GetYaxis().SetLabelSize(ratioh.GetYaxis().GetLabelSize() * 0.75 / 0.25)
    plothist.GetYaxis().SetNdivisions(205)
    plothist.GetYaxis().SetTitleSize(ratioh.GetYaxis().GetTitleSize() * 0.75 / 0.25)
    plothist.GetYaxis().SetTitleOffset(0.4)

    plothist.GetXaxis().SetLabelSize(ratioh.GetXaxis().GetLabelSize() * 0.75 / 0.25)

    pulls = calc_pulls(ratioh, fit_func)
    pulls.SetMarkerStyle(2)
    pulls.SetMarkerSize(1)

    plotOnCanvas(ppad, [pulls], drawOpt='PE')

    can.Draw()
    plot_name = '_'.join([outbase, 'fit'])
    can.SaveAs(plot_name + '.pdf')


def write_result_file(filename, fit_results):
    """Write the contour and the central results to a root file"""
    contour_file = r.TFile(filename, 'recreate')
    for result in fit_results:
        graph = result['contour']
        fit_graph = get_best_fit(result)
        if graph is not None:
            name = '_'.join(['contour', 'graph', 'errlevel', '{:.2f}'.format(result['err_level'])])
            fit_name = name.replace('contour', 'fit')
            if result['ndf'] == -1: # scan fit results get the ndf set to -1
                name += '_scan'
                fit_name += '_scan'

            graph.SetName(name)
            graph.Write()

            fit_graph.SetName(fit_name)
            fit_graph.Write()

    contour_file.Close()


def main(datafn, reffn, outbase, treen, n_bins, bin_thresh=0, scan_plot=True,
         use_center=False):
    # first run a scan to check how stable the normalization is
    # mean_norm = calc_mean_norm(datafn, reffn, treen, n_bins)
    mean_norm = get_free_norm(datafn, reffn, treen, n_bins)

    fit_func = r.TF1('fit_func',
                     '[0] * (3 + [1]) / (3 + [1] + [2]) * (1 + ([1] + [2]) * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                     -1, 1)

    datah, refh = get_histos(datafn, reffn, treen, n_bins)
    # for (i,b) in enumerate(datah):
    #     print('{}, data = {}, ref = {}'.format(i, b, refh.GetBinContent(i))

    set_bins_to_zero(datah, bin_thresh, True)
    set_bins_to_zero(refh, bin_thresh, True)

    ratioh = divide(datah, refh)
    ratio_fit = CosthRatioFit(ratioh, fit_func, fix_params=[(0, mean_norm)],
                              use_center=use_center)
    fit_results = ratio_fit.get_results(error_levels)

    plotbase = '_'.join([outbase, str(n_bins)])
    make_fit_plot(ratioh, fit_results[0], fit_func, plotbase)

    if scan_plot:
        x_ran, y_ran = get_plotting_range(fit_results[-1])
        scan_fit_result = make_paed_plot(ratioh, fit_func, error_levels, fit_results,
                                         '.'.join([plotbase, 'pdf']),
                                         n_grid=250, x_ran=x_ran, y_ran=y_ran)
        fit_results += scan_fit_result

    for result in fit_results:
        print_result(result)

    cfile_name = '_'.join([plotbase, 'contours']) + '.root'
    write_result_file(cfile_name, fit_results)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for running a costh ratio fit '
                                     'returning contours')
    parser.add_argument('data_file_name', help='data file name')
    parser.add_argument('ref_file_name', help='ref file name')
    parser.add_argument('-o', '--output_base', help='output base dir',
                        default='fit_output')
    parser.add_argument('-t', '--treename', default='chic_tuple',
                        help='name of tree in input files')
    parser.add_argument('-n', '--nbins', default=32, type=int,
                        help='number of bins to use in costh')
    parser.add_argument('-ns', '--noscan', default=False, action='store_true',
                        help='do not produce scan plot')
    parser.add_argument('-bt', '--binthresh', default=0, type=float,
                        help='minimum number of events that have to be in each bin of the'
                        ' data and reference histograms. Bins below this threshold will be'
                        ' set to zero')
    parser.add_argument('-c', '--usecenter', default=False, action='store_true',
                        help='Evaluate the fit function at the bin center instead of '
                        'integrating it over the whole bin in the fit.')

    args = parser.parse_args()

    r.gROOT.SetBatch()
    main(args.data_file_name, args.ref_file_name, args.output_base,
         args.treename, args.nbins, bin_thresh=args.binthresh,
         scan_plot=not args.noscan, use_center=args.usecenter)

#!/usr/bin/env python
"""
Script to plot chi2 values
"""
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import sys
import pandas as pd
import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from scipy.spatial import ConvexHull
from root_pandas import read_root

from run_ratio_fit_new import fix_non_free_params, create_scan_plot, add_clone_to_leg

# dataframe names to axis label conversion
_axis_labels = {
    'N': 'N',
    'lth_ref': '#lambda_{ref}',
    'delta_lth': '#Delta_{#lambda}'
}

def set_plot_range(col_map, xran, yran):
    """
    Set the plotting range on the color map in the x & y direction.
    """
    logging.debug('Setting plotting range in x & y direction')
    x_min = col_map.GetXaxis().GetBinLowEdge(1)
    x_max = col_map.GetXaxis().GetBinUpEdge(col_map.GetXaxis().GetNbins())
    y_min = col_map.GetYaxis().GetBinLowEdge(1)
    y_max = col_map.GetYaxis().GetBinUpEdge(col_map.GetYaxis().GetNbins())
    logging.debug('Original range was: x = [{},{}], y = [{},{}]'
                  .format(x_min, x_max, y_min, y_max))

    repl_none = lambda x, v: v if x is None else x
    x_ran = [repl_none(xran[0], x_min), repl_none(xran[1], x_max)]
    y_ran = [repl_none(yran[0], y_min), repl_none(yran[1], y_max)]

    logging.debug('New range is x = [{},{}], y = [{},{}]'
                  .format(x_ran[0], x_ran[1], y_ran[0], y_ran[1]))

    col_map.GetXaxis().SetRangeUser(x_ran[0], x_ran[1])
    col_map.GetYaxis().SetRangeUser(y_ran[0], y_ran[1])


def get_contour_best_fit(scan_res, x_param, y_param, err_level=1.0):
    """
    Get the contour and the best fit (includeing asymmetric
    uncertainties) from the scanned chi2 values by simply calculating
    the convex hull of all points with chi2 values smaller than
    min_chi2 + err_level.

    Args:
        scan_res (pandas.DataFrame): chi2 values and corresponding
        parameters, where all but the x_param and y_param are fixed to
        their value at the min chic2
        [x,y]_param (str): Names of the x & y parameters for which the
        uncertainty contour should be determined in scan_res
        min_chi2 (float): The minimum observed chi2 value
        err_level (float): The delta chi2 to be used to calculate the
        uncertainty contour. Default = 1.0, corresponding to the 39.9
        % 2D CI, or the 68.8 % 1D CI, when the projecting onto 1
        parameter.

    Returns:
        (contour (ROOT.TGraph), best_fit (ROOT.TGraphAsymmErrors)):
        tuple containing the contour with the x-y coordinates of the
        vertices of the convex hull and the best fit results with
        asymmetric uncertainties.
    """
    logging.debug('Determining uncertainty contour for parameters {} and {} at '
                  'err_level = {}'.format(x_param, y_param, err_level))
    # use np.argmin() to avoid problems with pandas argmin/idxmin
    min_chi2_idx = np.argmin(scan_res.chi2.values)
    min_chi2 = scan_res.iloc[min_chi2_idx]['chi2']
    logging.debug('Found min chi2 = {} at position {}'.format(min_chi2, min_chi2_idx))

    cont_sel = np.abs(scan_res.chi2 - min_chi2) < err_level
    sel_points = np.array([scan_res[cont_sel][x_param],
                           scan_res[cont_sel][y_param]]).T

    logging.debug('Determining convex hull with {} points with delta chi2 < {}'
                  .format(np.sum(cont_sel), err_level))
    hull = ConvexHull(sel_points)

    n_points = len(hull.vertices)
    n_x = sel_points[hull.vertices, 0]
    n_y = sel_points[hull.vertices, 1]
    logging.debug('Found convex hull with {} vertices'.format(n_points))

    logging.debug('Getting best fit results and asymmetric uncertainties')
    x_val, y_val = scan_res.iloc[min_chi2_idx][[x_param, y_param]]
    x_low = np.array([x_val - scan_res[cont_sel][x_param].min()])
    x_up = np.array([scan_res[cont_sel][x_param].max() - x_val])
    y_low = np.array([y_val - scan_res[cont_sel][y_param].min()])
    y_up = np.array([scan_res[cont_sel][y_param].max() - y_val])

    logging.info('Best fit: {} = {:.4f} (+{:.4f} -{:.4f})'.format(x_param, x_val, x_low[0], x_up[0]))
    logging.info('Best fit: {} = {:.4f} (+{:.4f} -{:.4f})'.format(y_param, y_val, y_low[0], y_up[0]))

    return (r.TGraph(n_points, n_x, n_y),
            r.TGraphAsymmErrors(1, np.array(x_val), np.array(y_val),
                                x_low, x_up, y_low, y_up))


def get_dataframe(infile):
    """Get the dataframe from the input file"""
    logging.debug('Getting DataFrame from {}'.format(infile))
    if infile.endswith('.pkl'):
        return pd.read_pickle(infile)
    if infile.endswith('.root'):
        return read_root(infile)
    logging.error('Infile does not have a parseable format: {}'
                  ' Valid formats are .root and .pkl'.format(inputfile))
    sys.exit(1)


def make_two_var_plot(chi2, var_comb, outbase,
                      max_dchi2=25, xran=[None,None], yran=[None,None]):
    """
    Make a 2D plot of the chi2 square difference.
    All other variables are fixed to their values at the minimum
    """
    logging.debug('Making combi plot of {1} vs {0}'.format(*var_comb))
    def get_plot_name(var_comb, outbase):
        """Get the filename for the output plot"""

        return '_'.join([outbase, var_comb[1], 'v', var_comb[0]]) +'.pdf'

    scan_plot = create_scan_plot(chi2, *var_comb, err_lvls=[1, 2.3],
                                 min_chi2=1e9, plot_thresh=25)
    if scan_plot is None:
        return

    logging.debug('Creating canvas for plotting')
    c = r.TCanvas('c', 'c', 800, 800)
    pad = c.cd()
    c.SetGrid()
    pad.SetRightMargin(0.15)
    pad.SetLeftMargin(0.12)

    set_plot_range(scan_plot, xran, yran)
    scan_plot.SetXTitle(_axis_labels[var_comb[0]])
    scan_plot.SetYTitle(_axis_labels[var_comb[1]])
    scan_plot.Draw('colz')

    plot_chi2 = fix_non_free_params(chi2, *var_comb)

    leg = r.TLegend(0.12, 0.1, 0.25, 0.25)

    contour1D, best_fit = get_contour_best_fit(plot_chi2, *var_comb, err_level=1.0)
    contour1D.SetLineColor(0)
    contour1D.SetLineStyle(2)
    contour1D.Draw('sameL')
    add_clone_to_leg(leg, contour1D, '39.9 % CI', 'L')

    best_fit.SetMarkerStyle(23)
    best_fit.SetMarkerColor(0)
    best_fit.Draw('samePX')
    add_clone_to_leg(leg, best_fit, 'best fit', 'P')

    contour2D, _ = get_contour_best_fit(plot_chi2, *var_comb, err_level=2.3)
    contour2D.SetLineColor(0)
    contour2D.SetLineStyle(3)
    contour2D.Draw('sameL')
    add_clone_to_leg(leg, contour2D, '68.8 % CI', 'L')

    leg.Draw()
    plot_name = get_plot_name(var_comb, outbase)
    logging.debug('Saving plot {}'.format(plot_name))
    c.SaveAs(plot_name)


def parse_range(range_str):
    """Parse the range from a passed string"""
    return [float(s) if s != 'None' else None for s in range_str.split(',')]


def main(args):
    """Main"""
    chi2_vals = get_dataframe(args.chi2file)
    if args.parameters is None:
        param_combis = [('N', 'delta_lth'), ('lth_ref', 'delta_lth'), ('N', 'lth_ref')]
    else:
        param_combis = [(args.parameters[0], args.parameters[1])]

    x_ran = parse_range(args.xrange)
    y_ran = parse_range(args.yrange)

    for par_comb in param_combis:
        make_two_var_plot(chi2_vals, par_comb, args.outbase,
                          max_dchi2=args.chi2range, xran=x_ran, yran=y_ran)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for making plots of a chi2 scan')
    parser.add_argument('chi2file', help='file containing the dataframe with the chi2 values.')
    parser.add_argument('outbase', help='Basename for created output plots and files.')
    parser.add_argument('-p', '--parameters', nargs=2, default=None,
                        help='Make 2D plot of passed parameters')
    parser.add_argument('-xr', '--xrange', type=str, default='None,None',
                        help='Plotting range for x parameter')
    parser.add_argument('-yr', '--yrange', type=str, default='None,None',
                        help='Plotting range for y parameter')
    parser.add_argument('-cr', '--chi2range', type=float, default=25,
                        help='Max delta chi2 value in plots')

    clargs = parser.parse_args()
    r.gROOT.SetBatch()

    main(clargs)

#!/usr/bin/env python
"""
Script that runs a scan of the chi2 values
"""

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import numpy as np
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.miscHelpers import parseVarBinning

from run_ratio_fit_new import scan_param_space


def get_ratio_hist(dataf, reff, treen, n_bins):
    """Get the 'cleaned up' ratio hist"""
    from run_ratio_histo_fit import get_histos, divide, set_bins_to_zero

    logging.debug('Creating histogram from datafile \'{}\' and '
                 'reffile \'{}\' using {} bins'.format(dataf, reff, n_bins))
    datah, refh = get_histos(dataf, reff, treen, n_bins)

    set_bins_to_zero(datah, 0, True)
    set_bins_to_zero(refh, 0, True)

    return divide(datah, refh)


def get_norm_scan_values(norm, hist, func, fix_norm=False):
    """
    Get the normalization parameter value(s).

    Args:
        norm: normalization settings to be parsed by parseVarBinning
        (i.e. list of strings)
        hist (ROOT.TH1D): histogram containing the data to fit if
        fix_norm is True
        func (ROOT.TF1): function to be fitted to the data if fix_norm
        is True
        fix_norm (bool): if true the normalization will be fixed to a
        free fit regardless of other settings from parsing

    Returns:
        norm_vals (numpy.array): the values to be used for fitting
    """
    logging.debug('Getting scan values for parameter N')
    if not fix_norm:
        vals = parseVarBinning(norm)

        # do some cleanup
        if np.any(vals <= 0):
            logging.info('Removing negative values from normalization scan params')
            vals = vals[vals > 0]
        if len(vals) == 0:
            logging.warning('After cleanup no normalization scan parameters are left.'
                           ' Fixing it to free fit.')
        else:
            return vals

    # if we are still here, we have to get the normalization from a free fit
    logging.info('Running a free fit to get the normalization')
    res = hist.Fit(func, 'S')
    logging.info('Fixing normalization to {}'.format(res.Parameter(0)))
    vals = np.array([res.Parameter(0)])

    return vals


def store_dataframe(df, outfile):
    """
    Store the dataframe either into a pkl file or into a root file via
    root_pandas.
    """
    logging.debug('Storing DataFrame to {}'.format(outfile))
    if not outfile.endswith('.pkl') and not outfile.endswith('.root'):
        logging.warning('Output file doesnot have .root or .pkl format. '
                        'Creating a .pkl file instead')
        logging.debug('Output filename before substitution: {}'.format(outfile))
        import re
        outfile = re.sub(r'(.*\.)(\w*)$', r'\1pkl', outfile)
        logging.debug('Output filename after substitution: {}'.format(outfile))

    logging.info('Writing resulting DataFramet to: {}'.format(outfile))
    # if .root is requested check if root_pandas is here, otherwise go to .pkl
    if outfile.endswith('.root'):
        try:
            from root_pandas import to_root
            to_root(df, outfile, 'chi2_values')
        except ImportError:
            logging.warning('Output to .root file was requested, but root_pandas'
                            ' was not found. Creating a .pkl file instead')
            outfile = outfile.replace('.pkl', '.root')

    if outfile.endswith('.pkl'):
        df.to_pickle(outfile)


def main(args):
    """Main"""
    fit_func = r.TF1('fit_func',
                     '[0] * (3 + [1]) / (3 + [1] + [2]) * (1 + ([1] + [2]) * x[0]*x[0]) / (1 + [1] * x[0]*x[0])',
                     -1, 1)

    # get the values that don't change with the binning here
    # the normalization is depending on the binning if it is fixed
    logging.debug('Getting parameter values for dlam = {} and lref = {}'
                  .format(args.dlam_range, args.lref_range))
    lref_values = parseVarBinning(args.lref_range)
    dlam_values = parseVarBinning(args.dlam_range)


    for n_bins in (int(v) for v in args.binnings.split(',')):
        ratioh = get_ratio_hist(args.datafile, args.reffile, args.treename, n_bins)
        norm_values = get_norm_scan_values(args.norm_range, ratioh,
                                           fit_func, args.fix_norm)
        scan_res = scan_param_space(ratioh, fit_func, {'N': norm_values,
                                                       'lref': lref_values,
                                                       'dlam': dlam_values})
        store_dataframe(scan_res, args.outfile.replace('XXX', str(n_bins)))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for scanning the parameter '
                                     'space and producing chi2 values or each point '
                                     'in the parameter space.'
                                     'The configurable values will be parsed by '
                                     'parseVarBinning')
    parser.add_argument('datafile', help='file containing the data tuple')
    parser.add_argument('reffile', help='file containing the reference tuple')
    parser.add_argument('outfile', help='output file to store the created data frame in.'
                        ' XXX will be replaced by the number of bins')
    parser.add_argument('-t', '--treename', default='chic_tuple',
                        help='name of the tree containting the tuple')
    parser.add_argument('-fn', '--fix-norm', default=False, action='store_true',
                        help='fix the normalization to the value from a free fit first')
    parser.add_argument('-b', '--binnings',default='48',
                        help='comma-separated list of bin numbers to use')
    parser.add_argument('-n', '--norm-range', default=['1'], nargs='+',
                        help='values to use for the normalization (only used if '
                        'fix-norm is not used!). If set to negative value (or zero) '
                        ' fix-norm will be assumed.')
    parser.add_argument('-l', '--lref-range', default=['0'], nargs='+',
                        help='values to use for the lambda_ref scan')
    parser.add_argument('-d', '--dlam-range', default=['0'], nargs='+',
                        help='values to use for the delta_lambda scan')

    clargs = parser.parse_args()

    r.gROOT.SetBatch()
    main(clargs)

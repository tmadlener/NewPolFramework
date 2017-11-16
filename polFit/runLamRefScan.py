#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
from ROOT import gROOT

from run_ratio_histo_fit import run
from runToyMCFits import read_gen_config, get_combinations
from utils.miscHelpers import condMkDir
from utils.plotHelpers import mkplot
from utils.TGraph_utils import createGraphSym
# import json

def run_scan(datafn, reffn, ref_range=[-0.3, 1], n_points=250):
    """
    Run a scan using n_points in ref_range
    """
    ref_scan_points = np.linspace(ref_range[0], ref_range[1], n_points)
    fit_results = []
    for lth_ref in ref_scan_points:
        fit_results.append(run(datafn, reffn, 'foo', 'genData',
                               chic1_limits=False, fix_ref=lth_ref,
                               fit_range=False, save=False))

    return fit_results


def make_graphs(lth_ref, delta_lth, delta_lth_err, chi2):
    """Create the delta labmda and chi2 graph from the passed fit_results"""
    lth_ref_err = np.zeros(len(lth_ref))
    delta_graph = createGraphSym(lth_ref, delta_lth, lth_ref_err, delta_lth_err)
    chi2_graph = createGraphSym(lth_ref, chi2, lth_ref_err, lth_ref_err) # since ndf is const, chi2 is enough

    return delta_graph, chi2_graph


def eval_fits(fit_results, datafn, plot=False, outdir='.'):
    """Evaluate the fit results"""
    lth_ref = np.array([r['lth_ref'][0] for r in fit_results])
    delta_lth = np.array([r['delta_lth'][0] for r in fit_results])
    delta_lth_err = np.array([r['delta_lth'][1] for r in fit_results])
    chi2 = np.array([r['chi2'] for r in fit_results])

    min_chi2_idx = np.argmin(chi2)
    lth_data_in = read_gen_config(datafn)['lthsig']

    if plot:
        delta_graph, chi2_graph = make_graphs(lth_ref, delta_lth + lth_ref, delta_lth_err, chi2)
        condMkDir(outdir)
        mkplot(delta_graph, saveAs='/'.join([outdir, 'delta_lth_v_lth_ref.pdf']), drawOpt='PE')
        mkplot(chi2_graph, saveAs='/'.join([outdir, 'chi2_v_lth_ref.pdf']), drawOpt='PE')

        return [lth_ref[min_chi2_idx], delta_lth[min_chi2_idx],
                delta_lth_err[min_chi2_idx], lth_data_in, np.min(chi2)]


def run_all_gens(gendir, outdir):
    data_ref_combis = get_combinations(gendir, outdir, True, ra_lth_ref=[0.5,0.5], ra_lth_data=[-0.5,-0.5])

    results = []

    for dataf, reff, outf in data_ref_combis:
        scan_res = run_scan(dataf, reff)
        results.append(np.array(eval_fits(scan_res, dataf, True, outf.split('.')[0])))

    df = pd.DataFrame(results, columns=['lth_ref', 'dlth', 'dlth_err', 'lth_in', 'min_chi2'])
    df.to_pickle('fit_results_df.pkl')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for running a scan of different'
                                     ' (fixed) lambda reference values over the ToyMC data')
    # parser.add_argument('dataFileName', help='filename of the data tuple')
    # parser.add_argument('refFileName', help='filename of the reference tuple')
    parser.add_argument('gendir', help='gen data directory')
    parser.add_argument('-o', '--outdir', help='output directory name', default='fit_results')

    args = parser.parse_args()

    gROOT.SetBatch()

    run_all_gens(args.gendir, args.outdir)

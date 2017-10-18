#!/usr/bin/env python

import glob
import os
import json
import pandas as pd
import numpy as np
from utils.miscHelpers import flatten

def get_matching_files(basedir, matchstr):
    """
    Get all files from (immediate) subdirectories of basedir that match the matchstr.
    """
    return glob.glob('/'.join([basedir, matchstr]))


def get_result_files_to_input(inputfile):
    """
    Get all result json files from the same directory as th inputfile is in
    """
    resdir = os.path.dirname(inputfile)
    return glob.glob('/'.join([resdir, 'fit_results_datagen*refgen*.json']))


def get_val(nested_dict, key):
    """
    Get the value to the key from the (possibly) nested dict values
    NOTE: this only works for the context we currently have, as it only descends down
    the first subdict in each dictionary recursively, so that it is possible that the key
    is in one of the dictionaiers, that gets never visited
    """
    if not key in nested_dict:
        for subdict in nested_dict.values():
            if not hasattr(subdict, '__iter__'):
                continue
            return get_val(subdict, key)
    else:
        return nested_dict[key]

    # the recursion has fallen through here, the key was nowhere to be found in
    # the nested dict (resp. not in the ones we looked at)
    return None


def get_all_vals(dict_list, key):
    """
    Get all values matching key from the list of (equally structured) dicts
    """
    vals = []
    for entry in dict_list:
        vals.append(get_val(entry, key))

    return vals

def get_final_results(resultjson):
    """
    Get the final results from the result json.
    For lph and ltp this is simply the last value, for lth it is the sum of all values
    added to the final result (from which also the error is taken)
    """
    with open(resultjson, 'r') as f:
        data = json.load(f)

    vals = {}

    # assuming here that the last result is stored in the json at the last position
    final_vals = data[-1]
    final_lth = get_val(final_vals, 'lth')
    # for lph and ltp simply return the result from the last fit
    vals['iteration'] = final_vals['iteration']
    vals['lph'] = get_val(final_vals, 'lph')
    vals['ltp'] = get_val(final_vals, 'ltp')
    vals['lth'] = final_lth

    all_lth = get_all_vals(data, 'lth')

    sum_ref_lth = np.sum(all_lth, axis=0)[0] # only need the central value

    res_lth = [final_lth[0] + sum_ref_lth, final_lth[1]]
    vals['lth_res'] = res_lth

    return vals


def get_input_lambdas(inputjson):
    """
    Get the input lambdas from the passed json and return them as a list, where
    the order is as follows: first data then reference: lth, lph, ltp
    """
    with open(inputjson, 'r') as f:
        data = json.load(f)
    return [data['data']['lthsig'], data['data']['lphsig'], data['data']['ltpsig'],
            data['ref']['lthsig'], data['ref']['lphsig'], data['ref']['ltpsig']]


def calc_mean_std_min_max(vals):
    """
    Calculate the (column-wise) mean, std deviation, min and max of the passed
    values
    """
    # put values into numpy array for easier computation of column wise values
    vals = np.array(vals)
    mean = vals.mean(axis=0)
    std = vals.std(axis=0)
    vmax = vals.max(axis=0)
    vmin = vals.min(axis=0)

    return np.array(list(flatten([mean, std, vmin, vmax])))


def check_converged(lth_vals, significance=0.1):
    """
    Check how many of the fits converged by checking the significance of the lth results
    """
    conv = 0
    for lth, err in lth_vals:
        if any(np.isnan([lth, err])):
            continue
        if abs(lth / err) < significance:
            conv += 1

    return conv


def get_result_lambdas(inputjson):
    """
    Get the resulting lambdas to the corresponding input json file
    Caluclates the mean, std, min and max value for the lth, lph, ltp and iterations
    Returns:
    - the mean and the std for lth, lph and ltp (as well as for their errors)
    - the min and max values for lth (but not for their errors)
    - the same values for lth_final, which is the result of the last fit
    - the mean, min and max values for iteraions
    - the number of converged fits
    - the number of total fits
    """
    resfiles = get_result_files_to_input(inputjson)
    raw_results = [get_final_results(rf) for rf in resfiles]

    lth_vals = get_all_vals(raw_results, 'lth_res')
    lth_finals = get_all_vals(raw_results, 'lth')
    lph_vals = get_all_vals(raw_results, 'lph')
    ltp_vals = get_all_vals(raw_results, 'ltp')
    iterations = get_all_vals(raw_results, 'iteration')
    # indices for the return values
    lambda_idcs = [0, 1, 2, 3, 4, 6]
    iter_idcs = [0, 2, 3]

    lth_res = calc_mean_std_min_max(lth_vals)[lambda_idcs]
    lph_res = calc_mean_std_min_max(lph_vals)[lambda_idcs]
    ltp_res = calc_mean_std_min_max(ltp_vals)[lambda_idcs]
    lth_final_res = calc_mean_std_min_max(lth_finals)[lambda_idcs]
    iterations = calc_mean_std_min_max(iterations)[iter_idcs]

    converged = check_converged(lth_finals)

    return [lth_res, lph_res, ltp_res, lth_final_res, iterations, converged, len(resfiles)]


def create_dataframe(basedir):
    """
    Create a dataframe collecting all the results from a given base directory.
    The data-frame will only contain the last (i.e. final) results for each iterative
    fit.
    There are no filters applied. nans in the dataframe are possible and probable!
    """
    # collect everything in an array of arrays and in the end construct the dataframe
    data_array = []

    in_columns = ['lth_data', 'lph_data', 'ltp_data', 'lth_ref', 'lph_ref', 'ltp_ref']
    res_columns = [
        'lth_mean', 'lth_err_mean', 'lth_std', 'lth_err_std', 'lth_min', 'lth_max',
        'lth_f_mean', 'lth_f_err_mean', 'lth_f_std', 'lth_f_err_std', 'lth_f_min', 'lth_f_max',
        'lph_mean', 'lph_err_mean', 'lph_std', 'lph_err_std', 'lph_min', 'lph_max',
        'ltp_mean', 'ltp_err_mean', 'ltp_std', 'ltp_err_std', 'ltp_min', 'ltp_max',
        'iter_mean', 'iter_min', 'iter_max', 'n_conv', 'n_fits'
    ]

    input_files = get_matching_files(basedir, '/*/input_lambdas.json')
    for infile in input_files:
        in_lambdas = get_input_lambdas(infile)
        res_lambdas = get_result_lambdas(infile)

        data_array.append(np.array(in_lambdas + list(flatten(res_lambdas))))

    return pd.DataFrame(data_array, columns=in_columns + res_columns)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script for analyzing the ToyMC fits')
    parser.add_argument('fitBaseDir', help='base directory where all the fits are stored')
    parser.add_argument('-o', '--outputfile', help='output (.pkl) file to which the data '
                        'frame will be stored.', dest='outfile', default='fit_results.pkl')

    args = parser.parse_args()

    frame = create_dataframe(args.fitBaseDir)
    frame.to_pickle(args.outfile)

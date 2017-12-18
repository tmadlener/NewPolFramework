#!/usr/bin/env python
"""
Script to average several chi2 scans
"""

import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s - %(funcName)s: %(message)s')
import pandas as pd
import numpy as np


from plot_chi2 import get_dataframe
from chi2_scan import store_dataframe

def check_compatible(df1, df2, pars):
    """
    Check if the two DataFrames can easily (i.e. trivially) be combined
    """
    # get the views into the parameters
    logging.debug('Checking if the two dataframes are compatible in the variables {}'
                  .format(pars))

    comp_dfs = df1[pars] == df2[pars]
    if np.all(comp_dfs):
        logging.debug('All values agree')
        return True

    logging.debug('There are incompatible values in the two dataframe')
    # COULDDO: some more diagnostics
    return False

def calc_exp_dchi2(chi2_vals):
    """Calculate the exponential of the negative delta chi2"""
    min_chi2 = chi2_vals.chi2.min()

    return np.exp(min_chi2 - chi2_vals.chi2)


def main(args):
    """Main"""
    comb_df = get_dataframe(args.chi2files.pop(0))
    comb_df['chi2'] = calc_exp_dchi2(comb_df)

    n_valid_files = 0

    for chi2file in args.chi2files:
        chi2_vals = get_dataframe(chi2file)
        if check_compatible(comb_df, chi2_vals, args.parameters):
            comb_df['chi2'] = comb_df['chi2'] + calc_exp_dchi2(chi2_vals)
            n_valid_files += 1

    comb_df['chi2'] = comb_df['chi2'] / n_valid_files
    comb_df['chi2'] = -np.log(comb_df['chi2'])
    store_dataframe(comb_df, args.outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Script to average several chi2 scans.')
    parser.add_argument('outfile', help='file to which the output dataframe will be stored')
    parser.add_argument('chi2files', nargs='+',
                        help='files containing dataframes fromc chi2 scans.')
    parser.add_argument('-p', '--parameters', nargs='+', default=['lth_ref', 'delta_lth'],
                        help='Parameters to use in compatibility check.')


    clargs = parser.parse_args()

    main(clargs)

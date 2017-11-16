#!/usr/bin/env python

import json
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from helpers import get_val
from utils.miscHelpers import flatten

result_jsons = glob.glob('lth_ref*/*/fit_results*.json')

def get_generations(filename):
    """Get the generation indices form the file names"""
    rgx = r'datagen_([0-9]+)_refgen_([0-9]+)'
    match = re.search(rgx, filename)
    if match:
        return [int(v) for v in match.groups()]


def create_df(jsons):
    """Create the dataframe containing the important data for all results"""
    columns = ['lth_in', 'lth_res', 'lth_err',
               'lph', 'lph_err', 'ltp', 'ltp_err', 'gen1', 'gen2', 'max_logL']

    data_array = []

    for res_json in jsons:
        with open(res_json, 'r') as f:
            data = json.load(f)[0]

            lth_in = data['lth_ref_input']
            lth_res = get_val(data, 'lth')
            lph = get_val(data, 'lph')
            ltp = get_val(data, 'ltp')
            max_logL = data['max_logL']
            gen1, gen2 = get_generations(res_json)

            data_array.append(np.array(list(flatten([lth_in] + lth_res + lph + ltp + [gen1, gen2, max_logL]))))

    return pd.DataFrame(data_array, columns=columns)


def make_overview_plot(df, bin_var, val_var, savename):
    """plot an overview of the val_var in bins of the bin_var in the dataframe df """
    group = df.groupby(bin_var)
    x_vals = group[val_var].mean().index
    mean = group[val_var].mean().values
    std = group[val_var].std().values
    v_max = group[val_var].max().values
    v_min = group[val_var].min().values

    fig = plt.figure()
    plt.errorbar(x_vals, mean, yerr=std, fmt='o', label='std')
    plt.errorbar(x_vals, mean, yerr=[mean - v_min, v_max - mean], fmt='o', label='min:max')

    plt.legend(loc='best')
    plt.xlabel(bin_var)
    plt.ylabel(val_var)
    plt.grid()
    fig.savefig(savename)


df = create_df(result_jsons)
df['pulls'] = (df['lth_res'] - (-0.5)) / df['lth_err']

df.to_pickle('result_df.pkl')
# df = pd.read_pickle('result_df.pkl')


## first make plots for each generation (yields 10 plots)
## sort dataframe by gen1
# df.sort_values(by=['gen1'])
# gen_vals = df['gen1'].unique().tolist()

make_overview_plot(df, 'lth_in', 'pulls', 'overview_pulls.pdf')
make_overview_plot(df, 'lth_in', 'lth_res', 'overview_lth_res.pdf')
make_overview_plot(df, 'lth_in', 'max_logL', 'overview_max_logl.pdf')
## print some summary information
# print('Pulls')
# for lth_in in df['lth_in'].unique().tolist():
#     lth_df = df.loc[df.lth_in == lth_in]

#     print('lth_ref = {}, mean = {}, std = {}'.format(lth_in, lth_df.pulls.mean(), lth_df.pulls.std()))

# print('log likelhood')
# for lth_in in df['lth_in'].unique().tolist():
#     lth_df = df.loc[df.lth_in == lth_in]

#     print('lth_ref = {}, mean = {}, std = {}, max = {}, min = {}'.format(lth_in, lth_df.max_logL.mean(), lth_df.max_logL.std(), lth_df.max_logL.max(), lth_df.max_logL.min()))

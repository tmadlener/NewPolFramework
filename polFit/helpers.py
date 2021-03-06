import logging

import pandas as pd
import numpy as np

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


def flat_list_tuple(tup_list):
    """
    From a list of tuples containing two elements, return a list of
    tuples, where the first elements is the tuple containing all first
    elemnts of the input list and the second contains one of all
    possible combinations of the items in the second elements of the
    tuple.
    """
    from itertools import product
    fix_funcs, fix_vals = zip(*tup_list) # split apart the fix_funcs and the values

    val_comb_gen = product(*[vals.tolist() for vals in fix_vals])
    return ((fix_funcs, val_comb) for val_comb in val_comb_gen)


def get_bin_centers(x_min, x_max, n_bins):
    """Get the values for the evaluation of the parameters in the bin centers"""
    bin_bounds = np.linspace(x_min, x_max, n_bins + 1)
    return 0.5 * (bin_bounds[1:] + bin_bounds[:-1])


def get_n_scan_points(scan_params):
    """Get the total number of scanning points to scan."""
    n = 1
    for _, vals in scan_params:
        # print(vals)
        n *= vals.shape[0]

    return n


def scan_chi2_params(hist, fit_func, scan_params, par_names):
    """
    Calculate the chisquare value between a histogram and a function
    for specific parameter sets.

    Args:
        hist (ROOT.TH1): The histogram containing the data to be fitted
        func (ROOT.TF1): The function to be fitted to the data.
        scan_params: A list of tuples, where each tuple corresponds to
        one paramter of func. Each tuple has to hold the following
        entries:
            [0] - A function taking a TF1 and a numeric value x, that
            is used to fixed the given parameter to the value x
            (e.g. func.FixParameter(x))
            [1] - All the values that should be tested for this
            parameter as a np.array
        par_names: A dictionary mapping integers to parameter
        names. Must have the same length as scan_params and start at
        0. The names will be used as column names for the returned
        DataFrame
    Returns:
        pandas.DataFrame: A DataFrame containing the values of all
        parameters and the corresponding chi2 value for each possible
        combination of parameter values in the input
    """
    from tqdm import tqdm

    n_params, n_names = len(scan_params), len(par_names)
    if n_params != n_names:
        if n_params < n_names:
            logging.warning('Got {} scan_params but {} par_names. Will ignore '
                           'the superfluous par_names'.format(n_params, n_names))
        else:
            logging.error('Got {} scan_params but only {} par_names. Need the '
                         'same number!'.format(n_params, n_names))

    n_total = get_n_scan_points(scan_params)
    logging.debug('Allocating results array of size {}'.format(n_total))

    # define the maximum number of scan point for which the (faster)
    # list of list approach is used before switching to a more memory
    # efficient numpy.array
    n_max_size = 30 * 1000**2
    if n_total < n_max_size:
        logging.debug('Using list of lists since size is smaller than {}'
                      .format(n_max_size))
        results = [[-1,-1,-1,-1] for _ in xrange(n_total)]
    else:
        logging.debug('Using numpy array, since size is greater than {}'
                      .format(n_max_size))
        results = np.zeros((n_total, 4), dtype=np.float64)

    logging.debug('Scanning {} points in total'.format(n_total))
    for j, (fix_funcs, val_comb) in enumerate(tqdm(flat_list_tuple(scan_params),
                                                   total=n_total, ncols=80)):
        for (i, func) in enumerate(fix_funcs):
            func(fit_func, val_comb[i])
            results[j][i] = val_comb[i]

        results[j][3] = hist.Chisquare(fit_func)

    logging.debug('Scanning done. Creating DataFrame')
    return pd.DataFrame(results, columns=par_names.values() + ['chi2'])

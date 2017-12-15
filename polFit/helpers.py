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
    n_params, n_names = len(scan_params), len(par_names)
    if n_params != n_names:
        if n_params < n_names:
            logger.warning('Got {} scan_params but {} par_names. Will ignore '
                           'the superfluous par_names'.format(n_params, n_names))
        else:
            logger.error('Got {} scan_params but only {} par_names. Need the '
                         'same number!'.format(n_params, n_names))

    results = []

    for fix_funcs, val_comb in flat_list_tuple(scan_params):
        pars = {}
        for (i, func) in enumerate(fix_funcs):
            func(fit_func, val_comb[i])
            pars[par_names[i]] = val_comb[i]

        pars['chi2'] = hist.Chisquare(fit_func)
        results.append(pars)

    return pd.DataFrame(results)

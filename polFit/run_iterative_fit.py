#!/usr/bin/env/ python

import argparse
import json
import os
import re
import numpy as np
from utils.TGraph_utils import createGraphSym
from utils.miscHelpers import tail

def parse_output(output):
    """
    Parse the output of one run (output is a list of tuples with only one entries each)

    Returns the final values of lth, lph and ltp (at the support points)
    """
    float_rgx = r'((-?[0-9].+[0-9]+)|(nan))' # floating point regex that also matches nans
    support_rgx = r'x = ' + float_rgx + r':\n'
    lth_rgx = r'\s*lth = ' + float_rgx + r' \+/- ' + float_rgx
    lph_rgx = r'\s*lph = ' + float_rgx + r' \+/- ' + float_rgx
    ltp_rgx = r'\s*ltp = ' + float_rgx + r' \+/- ' + float_rgx

    def match_and_extract_floats(string, rgx):
        """Match the string against the regex and extract all floating numbers from it"""
        m = re.match(rgx, string)
        if m:
            # with the regex defined as above we have for each match 3 groups, where
            # 1 + i is the index of the value that we want to convert to a float for the
            # i-th match
            group_idcs = xrange(1, len(m.groups()) + 1, 3)
            return [float(m.group(i)) for i in group_idcs]
        else:
            return None

    def cond_insert(d, curr_key, val, val_key):
        """conditionally insert val into d under val_key"""
        if curr_key is None or val is None: return
        d[curr_key][val_key] = val


    results = {}
    current_support = None

    for output_line in output:
        line = output_line[0]
        # print the output to stdout but avoid double newlines by stripping
        print(line.strip())

        support_try = match_and_extract_floats(line, support_rgx)
        if support_try is not None:
            current_support = support_try[0] # only one element in list in this case
            results[current_support] = {}

        cond_insert(results, current_support, match_and_extract_floats(line, lth_rgx), 'lth')
        cond_insert(results, current_support, match_and_extract_floats(line, lph_rgx), 'lph')
        cond_insert(results, current_support, match_and_extract_floats(line, ltp_rgx), 'ltp')

    return results


def run_fit(ref_lth, datafile, reffile, outfile, treename):
    """
    Run PoLPPD, with setting a reference lth, computed from the sum of the lths in the list
    """
    executable = '/'.join([os.environ['WORK'], 'NewPolMethodTests/Framework/polFit/runPolPPD'])

    fit_cmd = [executable, '--reffile', reffile, '--sigfile', datafile,
               '--sigtree', treename, '--reftree', treename, '--outfile', outfile,
               '--reflth']

    fit_cmd.append(str(ref_lth))

    print('========== Starting with lth_ref = {} ==========='.format(ref_lth))
    return parse_output(tail(fit_cmd))


def check_stop(lambdas, sigmas=1.0):
    """
    Check if lth of the results is close to zero (measured in significances)

    NOTE: currently only checks at one support point
    """
    # get any support point (they all have the same values currently with a constant fit)
    lth = lambdas.values()[0]['lth']

    # stop if there is a nan in the results
    if any(np.isnan(lth)): return True
    return abs(lth[0] / lth[1]) < sigmas


def calc_ref_lth(iteration_results):
    """
    Calculate the reference lambda from the previous iteration results

    NOTE: currently only constant lth is supported (no parametrization)
    """
    if len(iteration_results) == 0: return 0 # no iteration yet, so start at 0
    ref_lth = 0
    for res in iteration_results:
        # very ugly currently but that is the format we have and for now it works
        # currently seems to work also with added "iteration" key-value pair
        ref_lth -= res.values()[0]['lth'][0]

    return ref_lth


def run(datafile, reffile, outfile, treename, sigmaStop, maxIterations):
    """
    Run the iterative scheme, either until it converges (determined by sigmaStop) or until
    the maximum number of iterations is reached
    """
    iter_results = []
    it = 0
    while maxIterations < 0 or it < maxIterations:
        curr_ref = calc_ref_lth(iter_results)
        results = run_fit(curr_ref, datafile, reffile, outfile, treename)


        results["iteration"] = it # embed information of iteration into the results
        iter_results.append(results)

        if check_stop(results): break

        it += 1

    with open(outfile.replace('.root', '.json'), 'w') as f:
        json.dump(iter_results, f, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='script to run the iterative fitting procedure')
    parser.add_argument('dataFileName', help='filename of data tuple')
    parser.add_argument('refFileName', help='filename of ref tuple')
    parser.add_argument('outFileName', help='output file name', default='polPPD.root')
    parser.add_argument('-t', '--treename', help='treename', dest='tree', default='genData')
    parser.add_argument('-n', '--nMaxIterations', help='maximum number of iterations',
                        dest='maxIterations', default=-1, type=int)
    parser.add_argument('-s', '--stopSignificance', help='significance to stop the iteration',
                        default=0.25, type=float, dest='sigmas')

    args = parser.parse_args()

    run(args.dataFileName, args.refFileName, args.outFileName, args.tree,
        args.sigmas, args.maxIterations)

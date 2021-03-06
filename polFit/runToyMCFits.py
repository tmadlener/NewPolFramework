#!/usr/bin/env python

import argparse
import os
import json
import glob
import re
import datetime as dt
from utils.batch_utils import get_job_id, append_to_json
from utils.miscHelpers import condMkDirFile, tail, getBinIdx


def read_gen_config(genfile):
    """Read the generator settings from the generator dir"""
    gendir = os.path.dirname(genfile)
    with open('/'.join([gendir, 'input_lambdas.json']), 'r') as f:
        data = json.load(f)
    return data


def submit_job(datafile, reffile, outfile, maxIterations, significance, checkExact=False,
               checkExactDiff=False, lthStart = 0.0):
    """Submit a batch job that runs one fit using the passed files"""
    condMkDirFile(outfile)
    outdir = os.path.dirname(outfile)

    run_mode = 'normal' # for storage in the json file

    # combine data and refernce generator settings into one json and dump it to the output directory
    comb_config = {'data': read_gen_config(datafile),
                   'ref': read_gen_config(reffile)}
    with open('/'.join([outdir, 'input_lambdas.json']), 'w') as f:
        json.dump(comb_config, f, indent=2)

    if checkExact or checkExactDiff:
        max_iter = 1
        if checkExact and not checkExactDiff:
            lth_start = comb_config['ref']['lthsig']
            run_mode = 'exact_ref'
        if checkExactDiff and not checkExact:
            lth_start = comb_config['data']['lthsig'] - comb_config['ref']['lthsig']
            run_mode = 'exact_diff'
    else:
        max_iter = maxIterations
        lth_start = lthStart

    batch_script = '/'.join([os.environ['WORK'], 'NewPolMethodTests/Framework/polFit/batchRunIterative.sh'])
    cmd_list = ['sbatch', batch_script, datafile, reffile, 'genData', outfile,
                str(max_iter), str(significance), str(lth_start)]

    # batch submission yields only 1 element in the list, from which the job id will be
    # retrieved. second index is for getting first (and only) element of tuple returned
    # by tail
    job_id = get_job_id(list(tail(cmd_list))[0][0])

    append_to_json('/'.join([outdir, 'batch_job_info.json']),
                   {'job_id': job_id,
                    'sub_time': str(dt.datetime.now()),
                    'datafile': datafile,
                    'reffile': reffile,
                    'outfile': os.path.basename(outfile),
                    'lth_start': lth_start,
                    'run_mode': run_mode
                   })


def get_io_combi(datafile, reffile, outbase, only_same_gen=False,
                 ra_lth_ref=[-1,1], ra_lth_data=[-1,1], shift_gen=False):
    """
    Get the three tuple of datafile, inputfile, outputfile
    Default arguments are set, such that apriori every passed combination is accepted
    - shift_gen pairs generations which are off by one
    - only_same_gen pairs only same generations
    - they cannot be set simultaneously (will always return false that way)
    """
    def get_lth(filename):
        m = re.search('lth_(-?[0-9]+p[0-9]{2})', filename)
        if m:
            return m.group(1)
        else:
            return None

    def check_range(valstr, range_vals):
        """Check if the valstr (where the comma dot is a 'p' is in range_vals)"""
        val = float(valstr.replace('p', '.'))
        return val >= range_vals[0] and val <= range_vals[1]

    data_gen = getBinIdx(datafile, 'gen_')
    ref_gen = getBinIdx(reffile, 'gen_')
    if only_same_gen and ref_gen != data_gen:
        return None

    # if we want shifted (by 1) generations, check if data_gen == ref_gen + 1
    # the % 10 is to wrap around at 10. (So that 10 and 1 are paired)
    # NOTE: If there are more then 10 generations this will break (i.e. lead to less
    # thane xpected combinations)
    if shift_gen:
        if data_gen % 10 != (ref_gen + 1) % 10:
            return None

    data_lth = get_lth(datafile)
    ref_lth = get_lth(reffile)

    if not check_range(data_lth, ra_lth_data) or not check_range(ref_lth, ra_lth_ref):
        return None

    outdir = '_'.join(['data', 'lth', data_lth, 'ref', 'lth', ref_lth])
    ofn = '_'.join(['fit', 'results', 'datagen', str(data_gen), 'refgen', str(ref_gen)])
    outfile = '/'.join([outbase, outdir, ofn + '.root'])

    return (datafile, reffile, outfile)


def get_combinations(inputbase, outbase, only_same_gen=True,
                     ra_lth_ref=[-1,1], ra_lth_data=[-1,1]):
    """
    Get the combinations of data and reference files and produce the appropriate output
    file for it.
    """
    combis = []
    gen_files = glob.glob('/'.join([inputbase, 'lth_*_lph_*_ltp_*/genData_gen_*.root']))

    for df in gen_files:
        for rf in gen_files:
            # use shift_gen if only_same_gen is specified, since it is here mainly to
            # avoid a combinatoric explosion
            # NOTE: this currently assumes that there are only 10 generations
            combi = get_io_combi(df, rf, outbase, False,
                                 ra_lth_ref, ra_lth_data, only_same_gen)
            if combi is not None:
                combis.append(combi)

    return combis


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='script for submitting fitting jobs to '
                                     'the batch system.')
    parser.add_argument('genDataDir', help='directory that holds the generated ToyMC data'
                        ' in the format as it is generated by \'generateToyMC.py\'')
    parser.add_argument('outdir', help='output base directory for the fit results')
    parser.add_argument('-a', '--allCombis', action='store_true', dest='allCombis',
                        default=False, help='run all possible combinations of generations,'
                        ' instead of only running same generation data and reference')
    parser.add_argument('-c', '--combinations', help='run specific combination of passed files',
                        nargs='+', dest='combifiles', default=[])
    parser.add_argument('-n', '--nMaxIterations', help='maximum number of iterations',
                        dest='maxIterations', default=-1, type=int)
    parser.add_argument('-s', '--stopSignificance', help='significance to stop the iteration',
                        default=0.25, type=float, dest='sigmas')
    parser.add_argument('-x', '--checkExact', default=False, dest='checkExact',
                        action='store_true', help='Run the fit only once, but use the '
                        'exact lambda_theta reference as input to check if the output '
                        'is lambda_theta data')
    parser.add_argument('-d', '--checkExactDiff', default=False, action='store_true',
                        dest='checkExactDiff', help='Run the fit only once but use the '
                        'exact lambda difference (data - reference) as input to check if '
                        'the output is zero')
    parser.add_argument('-l', '--lthRefStart', type=float, default=0.0, dest='lthRefStart',
                        help='Starting value for lth_ref in the iterative fit')
    parser.add_argument('-rlr', '--rangeLthRef', default='-1,1', dest='rangeLthRef',
                        help='Range of lth_ref to consider in the combinations')
    parser.add_argument('-rld', '--rangeLthData', default='-1,1', dest='rangeLthData',
                        help='Range of lth_data to consider in the combinations')


    args = parser.parse_args()

    if len(args.combifiles) not in (0, 2):
        parser.error('Need non or exactly two arguments for --combinations')

    if args.checkExact and args.checkExactDiff:
        parser.error('Only one of checkExact and checkExactDiff can be true at the same time')

    if len(args.combifiles) > 0:
        data_ref_combis = [get_io_combi(args.combifiles[0], args.combifiles[1], args.outdir)]
    else:
        rangeLthRef = sorted([float(v) for v in args.rangeLthRef.split(',')])
        rangeLthData = sorted([float(v) for v in args.rangeLthData.split(',')])
        if len(rangeLthRef) != 2 or len(rangeLthData) != 2:
            parser.error('rangeLthRef and rangeLthData need to have two (comma-separated)'
                         ' elements that can be converted to floats')

        data_ref_combis = get_combinations(args.genDataDir, args.outdir, not args.allCombis,
                                           rangeLthRef, rangeLthData)

    for (datan, refn, outn) in data_ref_combis:
        submit_job(datan, refn, outn, args.maxIterations, args.sigmas,
                   args.checkExact, args.checkExactDiff, args.lthRefStart)

    print('Submitted {} jobs to the batch system'.format(len(data_ref_combis)))

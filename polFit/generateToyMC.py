#!/usr/bin/env python

import os
import json
import datetime as dt
import argparse
from utils.miscHelpers import condMkDir, tail, parseVarBinning
from utils.batch_utils import get_job_id, append_to_json


def submit_job(config, output_base_dir, N_start=1, N_end=1):
    """Submit a batch array job using the passed config. N specifies how many jobs should be launched"""

    rep_p = lambda x : '{:.2f}'.format(x).replace('.', 'p')
    lambda_args = [str(config[x]) for x in ['lthsig', 'lthbkg', 'lphsig', 'lphbkg', 'ltpsig', 'ltpbkg']]
    output_spec = '_'.join(['lth', rep_p(config['lthsig']), 'lph', rep_p(config['lphsig']),
                            'ltp', rep_p(config['ltpsig'])])

    output_dir = '/'.join([output_base_dir, output_spec])
    condMkDir(output_dir)

    # write the used config to the output directory for later retrieval if necessary
    with open('/'.join([output_dir, 'input_lambdas.json']), 'w') as of:
        json.dump(config, of, indent=2, sort_keys=True)

    batch_script = '/'.join([os.environ['WORK'], 'NewPolMethodTests/Framework/polFit/batchRunPolGen.sh'])
    outfile_base = 'genData'

    command_list = ['sbatch', '--array={}-{}'.format(N_start, N_end), batch_script, outfile_base]
    command_list += lambda_args + [str(config['nEvents']), output_dir]

    # the command yields only 1 element in the list, from which the job id will be retrieved
    # have to get the firt element of the tuple as well since tail returns a list of tuples
    job_id = get_job_id(list(tail(command_list))[0][0])

    append_to_json('/'.join([output_dir, 'batch_job_info.json']),
                   {'job_id': job_id, 'N_start': N_start, 'N_end': N_end,
                    'sub_time': str(dt.datetime.now())})


def create_configs(args):
    """
    Create the configurations needed by submit_job from the args parsed
    by argparse.
    """
    configs = []

    for lthsig in parseVarBinning(args.lthsig):
        for lthbkg in parseVarBinning(args.lthbkg):
            for lphsig in parseVarBinning(args.lphsig):
                for lphbkg in parseVarBinning(args.lphbkg):
                    for ltpsig in parseVarBinning(args.ltpsig):
                        for ltpbkg in parseVarBinning(args.ltpbkg):
                            configs.append(
                                {
                                    'lthsig': lthsig, 'lthbkg': lthbkg,
                                    'lphsig': lphsig, 'lphbkg': lphbkg,
                                    'ltpsig': ltpsig, 'ltpbkg': ltpbkg,
                                    'nEvents': args.nEvents
                                })

    return configs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='script for submitting jobs to the batch'
                                     ' system for generating ToyMC samples. The lambda '
                                     'input values understand the following syntax:\n'
                                     ' - list of numbers to use\n'
                                     ' - X:Y:Z linearly spaced array of numbers, where Y'
                                     ' is the spacing used to generate numbers between X'
                                     ' and Z\n'
                                     ' - X:Z,N linearly spaced array including X and Z'
                                     ' containing N entries in total',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--lthsig', help='values for signal lth to use', nargs='+',
                        dest='lthsig', default=['0'])
    parser.add_argument('--lthbkg', help='values for bkgnal lth to use', nargs='+',
                        dest='lthbkg', default=['0'])
    parser.add_argument('--lphsig', help='values for signal lph to use', nargs='+',
                        dest='lphsig', default=['0'])
    parser.add_argument('--lphbkg', help='values for bkgnal lph to use', nargs='+',
                        dest='lphbkg', default=['0'])
    parser.add_argument('--ltpsig', help='values for signal ltp to use', nargs='+',
                        dest='ltpsig', default=['0'])
    parser.add_argument('--ltpbkg', help='values for bkgnal ltp to use', nargs='+',
                        dest='ltpbkg', default=['0'])
    parser.add_argument('--nEvents', help='Number of events to generate', type=int,
                        dest='nEvents', default=10000)
    parser.add_argument('--nJobs', help='Number of generations per setting', type=int,
                        dest='nJobs', default=1)
    parser.add_argument('-o', '--outbasedir', help='output base directory (necessary '
                        'for checking if jobs were successful', dest='outbasedir',
                        default='genData')

    args = parser.parse_args()

    for c in create_configs(args):
        submit_job(c, args.outbasedir, 1, args.nJobs)

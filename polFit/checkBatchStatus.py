#!/usr/bin/env python

import argparse
import os
from utils.batch_utils import check_batch_file
import sys

def get_batch_info_files(basedir, jsonname):
    """
    Get all json file (names) containing the batch submission infos from the directories
    in the basedir.

    """
    # import glob
    # return glob.glob('/'.join([basedir, '*', 'batch_job_info.json']))

    info_files = []
    for directory,_,files in os.walk(basedir):
        for f in files:
            if f == jsonname:
                info_files.append(os.path.join(directory, f))
                break # there will only be one, so don't bother to search further

    return info_files

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Check the status of all batch jobs in'
                                     ' the passed subdirectory')
    parser.add_argument('checkdir', help='directory to check')
    parser.add_argument('-i', '--infofilename', help='name of the json in which the batch '
                        'submission info is stored in each subdirectory', dest='jsonname',
                        default='batch_job_info.json')

    args = parser.parse_args()

    failure_dirs = []
    for bf in get_batch_info_files(args.checkdir, args.jsonname):
        print('Checking jobs in {} ...'.format(bf))
        if check_batch_file(bf):
            print('All jobs completed successfully')
        else:
            failure_dirs.append(os.path.dirname(bf))
            print('Not all jobs completed successfully')

    if failure_dirs:
        print('The following directories had jobs that were not succesfully completed:')
        for fdir in failure_dirs:
            print(fdir)

        sys.exit(1)

#!/usr/bin/env python

import sys

import logging
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s - %(funcName)s: %(message)s')

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True

from utils.recurse import collectGraphs
from utils.plotHelpers import mkplot,  _defaultColors

# global array of rootfiles, to avoid getting them closed on going out of scope
_open_files = []

def collect_graphs(inputfiles):
    """Collect all the graphs from the input files"""
    graphs = {}
    global _open_files
    for inf in inputfiles:
        logging.debug('Opening file {}'.format(inf))
        f = r.TFile.Open(inf)
        _open_files.append(f)

        graphs[inf] = collectGraphs(f)

        logging.debug('Collected {} graphs'.format(len(graphs[inf])))

    return graphs


def parse_file_name(filename):
    """Get the bin threshold and the number of bins from the filename"""
    import re
    rgx = r'bin_thresh_([0-9]+).*n_bins_([0-9]+)'
    m = re.search(rgx, filename)
    if m:
        logging.debug('Matching \'{}\' to \'{}\' worked: {}'.format(rgx, filename, m.groups()))
        return [int(m.group(i)) for i in [1,2]]

    logging.warning('Could not match \'{}\' to \'{}\''.format(rgx, filename))
    return -1,-1


def select_graphs(graphs, sel_str):
    """From the graphs of all files select the ones matching"""
    from utils.miscHelpers import filterDict

    logging.debug('Selecting graphs matching \'{}\''.format(sel_str))
    sel_graphs = {}
    for inf in graphs:
        sel_cands = filterDict(graphs[inf], sel_str).values()
        logging.debug('Found {} matching graphs in \'{}\''.format(len(sel_cands), inf))
        if len(sel_cands) == 0:
            logging.warning('Could not get graphs matching \'{}\' '
                           'from file {}'.format(sel_str, inf))
            continue
        if len(sel_cands) > 1:
            # sort by length, since '_scan' is appended to the files
            sel_cands.sort(key=lambda x: len(x.GetName()))
            logging.info('Found {} graphs matching \'{}\' '
                        'in file {}, selecting {}'.format(len(sel_cands), sel_str,
                                                          inf, sel_cands[0]))
            logging.debug('Other candidates were {}'.format(sel_cands[1:]))

        sel_graphs[inf] = sel_cands[0]

    logging.debug('Selected {} graphs for {} files'.format(len(sel_graphs), len(graphs)))
    return sel_graphs


def sort_graphs(graphs):
    """Sort the graphs according to the bin threshold and then number of bins"""
    # first create a list of tuples so that keys and values stick together but
    # are orderable
    from operator import itemgetter

    logging.debug('Sorting graphs. Creating sorting list of tuple')
    gtuples = zip((parse_file_name(k) for k in graphs.keys()), graphs.values())

    gtuples.sort(key=itemgetter(0,1))
    logging.debug('Done sorting graphs')

    return (g[0] for g in gtuples), (g[1] for g in gtuples)


def set_marker_styles(graphs):
    """
    Set different marker styles, such that each marker - color
    combination is present only once
    """
    n_colors = len(_defaultColors())
    marker_styles = (20, 21, 22, 34, 47, 23, 33)

    for i, g in enumerate(graphs):
        i_col = (i + 1) / n_colors
        g.SetMarkerStyle(marker_styles[i_col])
        g.SetMarkerSize(1.5)


def set_line_styles(graphs):
    """
    Set different line styles, such that each line - color combination
    is present only once
    """
    n_colors = len(_defaultColors())
    line_styles = (1, 2, 8, 6)

    for i, g in enumerate(graphs):
        i_col = (i + 1) / n_colors
        g.SetLineStyle(line_styles[i_col])
        g.SetLineWidth(2)


def plot_graphs(graphs, sel_str, plotname, contour, ranges):
    """Plot all graphs matching the sel_str"""
    logging.debug('Creating plot for sel_str = \'{}\''.format(sel_str))
    sel_graphs = select_graphs(graphs, sel_str)
    logging.debug('Got {} graphs'.format(len(sel_graphs)))

    binnings, graphs = sort_graphs(sel_graphs)
    graphs = list(graphs) # need a list in any case

    if not contour:
        set_marker_styles(graphs)
        draw_opt = 'PE'
    else:
        set_line_styles(graphs)
        draw_opt='L'

    make_leg = lambda x: 'thresh = {}, bins = {}'.format(x[0], x[1])

    mkplot(graphs,
           legEntries=[make_leg(b) for b in binnings],
           xRange=ranges[0], yRange=ranges[1],
           saveAs=plotname, grid=True, drawOpt=draw_opt,
           legPos='botleft', yLabel='#Delta_{#lambda}', xLabel='#lambda_{ref}')


def main(inputfiles, contour, errorbars, outbase, errlevel, axis_ran_str=''):
    """Main"""
    # collect all graphs and then handle them according to args
    logging.info('Collecting graphs from {} files'.format(len(inputfiles)))
    all_graphs = collect_graphs(inputfiles)

    err_rgx = r'errlevel_' + errlevel.replace('.', '\.')

    if axis_ran_str:
        logging.debug('Getting axis ranges from {}'.format(axis_ran_str))
        ranges = [float(v) for v in axis_ran_str.split(',')]
        if len(ranges) < 4:
            print('Could not get 4 values from {}'.format(axis_ran_str))
            sys.exit(1)

        axis_ranges = []
        axis_ranges.append([ranges[0], ranges[1]])
        axis_ranges.append([ranges[2], ranges[3]])
        logging.debug('Axis ranges are: {} and {}'.format(axis_ranges[0],
                                                         axis_ranges[1]))

    else:
        axis_ranges = [None, None]

    if contour:
        plotname = outbase + '_contour.pdf'
        logging.info('Making contour plot \'{}\''.format(plotname))
        plot_graphs(all_graphs, r'contour.*' + err_rgx, plotname, True, axis_ranges)

    if errorbars:
        plotname = outbase +'_central.pdf'
        logging.info('Making central results plot \'{}\''.format(plotname))
        plot_graphs(all_graphs, r'fit.*' + err_rgx, plotname, False, axis_ranges)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Script for creating plots comparing'
                                     ' different uncertainty shapes')
    parser.add_argument('inputfiles', nargs='+', help='input files to use')
    parser.add_argument('-o', '--outbase', default='fit_results',
                        help='base name for the output plots')
    parser.add_argument('-c', '--contour', default=True, action='store_true',
                        help='make plot with contours')
    parser.add_argument('-nc', '--nocontour', dest='contour', action='store_false',
                        help='do not make plot with contours')
    parser.add_argument('-e', '--errorbars', default=False, action='store_true',
                        help='make plot with errorbars')
    parser.add_argument('-ne', '--noerrorbars', dest='errorbars', action='store_false',
                        help='do not make plot with errorbars')
    parser.add_argument('-l', '--errlevel', type=str, default='1.00',
                        help='desired error level (has to be present in files)')
    parser.add_argument('-z', '--zoom', help='zoom into desired region', default='')

    args = parser.parse_args()

    r.gROOT.SetBatch()
    r.gROOT.ProcessLine('gErrorIgnoreLevel = 1001')

    main(args.inputfiles, args.contour, args.errorbars, args.outbase, args.errlevel,
         args.zoom)

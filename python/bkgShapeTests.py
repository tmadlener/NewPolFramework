import argparse

parser = argparse.ArgumentParser(description='script to make plots for checking dependenc of weight on costh')
parser.add_argument('file', help='file containing all the information from the fit')
parser.add_argument('outdir', help='output_directory')

args = parser.parse_args()

from utils.miscHelpers import condMkDir
condMkDir(args.outdir)

import ROOT as r
import ROOT.RooFit as rf

r.gROOT.SetBatch()

"""
Helper functions
"""
def getVar(ws, name):
    var = ws.var(name)
    if not var: var = ws.function(name)

    return var


def getSRBoundaries(ws, state, nSig=3):
    """
    Get the boundaries of the signal region, defined by the state (1S, 2S or 3S)
    """
    mean = getVar(ws, 'mean' + state).getVal()
    sigma = getVar(ws, 'sigma' + state).getVal()

    return [mean - nSig*sigma, mean + nSig*sigma]


def makePlot(ws, raprange, SRs=None):
    """
    Plot the data and the fitted model
    """
    from utils.miscHelpers import createRandomString
    mass = ws.var('mass')
    model = ws.pdf('fullModel')
    data = ws.data('_'.join(['data', 'absCosth', raprange]))

    frame = mass.frame(rf.Range('fitRange'))
    ws.loadSnapshot('_'.join(['snap', 'absCosth', raprange]))

    data.plotOn(frame)
    model.plotOn(frame)

    name = createRandomString(8)
    c = r.TCanvas(name, 'c', 1000, 1000)
    c.cd()
    frame.Draw()

    if SRs is not None:
        from utils.plotHelpers import _defaultColors
        colors = {str(i) + 'S': col for (i, col) in enumerate(_defaultColors())}

        ymin = 0
        ymax = frame.GetMaximum()
        lines = [] # keep list, otherwise root will only plot one object due to gc
        for (s, vals) in SRs.iteritems():
            for v in vals:
                lines.append(r.TLine(v, ymin, v, ymax));
                lines[-1].SetLineColor(colors[s])
                lines[-1].Draw('same')


    c.SaveAs(args.outdir + '/plot_fit_' + raprange + '.pdf')


def calcFractions(ws, raprange):
    """
    Calculate the ratio of background events in the total range to the background
    events in the sideband regions
    """
    from math import sqrt
    states = ['1S', '2S', '3S']

    ws.loadSnapshot('_'.join(['snap', 'absCosth', raprange]))

    SRs = {}
    for state in states: SRs[state] = getSRBoundaries(ws, state, 3.0)
    # cleanup possible overlap between SR2S and SR3S
    if SRs['2S'][1] > SRs['3S'][0]:
        middleVal = (SRs['2S'][1] + SRs['3S'][0]) * 0.5
        SRs['2S'][1] = middleVal
        SRs['3S'][0] = middleVal

    makePlot(ws, raprange, SRs)

    # set ranges in mass variable
    mass = ws.var('mass')
    mass.setRange('_'.join(['RSB', raprange]), 8.6, SRs['1S'][0])
    mass.setRange('_'.join(['LSB', raprange]), SRs['3S'][1], 11.4)
    mass.setRange('_'.join(['IB', raprange]), SRs['1S'][1], SRs['2S'][0])
    for (s, vals) in SRs.iteritems():
        mass.setRange('_'.join([s, raprange]), vals[0], vals[1])

    bkgModel = ws.pdf('bkgPoly')
    bkgInt = bkgModel.createIntegral(r.RooArgSet(mass), rf.NormSet(r.RooArgSet(mass)),
                                     rf.Range('_'.join(['RSB', raprange])),
                                     rf.Range('_'.join(['LSB', raprange])),
                                     rf.Range('_'.join(['IB', raprange]))).getVal()


    srInt = {s: bkgModel.createIntegral(r.RooArgSet(mass), rf.NormSet(r.RooArgSet(mass)),
                                        rf.Range('_'.join([s, raprange]))).getVal()
             for s in SRs}

    return {s: (bkgInt + srI) / bkgInt for (s, srI) in srInt.iteritems()}


def createCosThBins(vals):
    low = [str(v) for v in vals[:-1]]
    high = [str(v) for v in vals[1:]]

    merge = [low[i] + 'to' + high[i] for i in range(len(low))]

    return [s.replace('.', 'p') for s in merge]


def makeGraphs(fracs, basename):
    from utils.TGraph_utils import createGraph
    from utils.miscHelpers import getBinEdgeVals
    from numpy import diff

    costh = [sum(getBinEdgeVals(n))*0.5 for n in fracs]
    costhErr = [diff(getBinEdgeVals(n))*0.5 for n in fracs]


    graphs = []
    for s in ['1S', '2S', '3S']:
        fs = [f[s] for (_,f) in fracs.iteritems()]

        graphs.append(createGraph(costh, fs, costhErr, costhErr, [0]*len(fs), [0]*len(fs)))
        graphs[-1].SetName('_'.join([basename, s]))

    return graphs


def getNEvents(ws, costhRanges):
    """foo"""
    from utils.miscHelpers import getBinEdgeVals

    print('number of events in |costh| bins:')
    for rg in costhRanges:
        ds = ws.data('_'.join(['data', 'absCosth', rg]))
        print('{}: {}'.format(getBinEdgeVals(rg), ds.numEntries()))


def printVar(ws, name, costhRanges):
    """
    Print variable from workspace for all costhRanges
    """
    from utils.miscHelpers import getBinEdgeVals

    print('values of {}:'.format(name))
    var = getVar(ws, name)
    for rg in costhRanges:
        ws.loadSnapshot('_'.join(['snap', 'absCosth', rg]))
        print('{}: {:.3f} +/- {:.3f}'.format(getBinEdgeVals(rg), var.getVal(), var.getError()))


"""
Script
"""
abscosthRan = createCosThBins([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1])

ff = r.TFile.Open(args.file)
ws = ff.Get('workspace')

fractions = {}
for rg in abscosthRan:
    fractions[rg] = calcFractions(ws, rg)

gFile = r.TFile(args.outdir + '/graphs.root', 'recreate')
graphs = makeGraphs(fractions, 'fracs')
for g in graphs:
    g.Write()


from utils.plotHelpers import mkplot
can = mkplot(graphs, xRange=[0,1.0], xLabel='|cos#theta|', yLabel='f_{total over SB}', ret=True,
             leg=True, drawOpt='P')
can.SaveAs(args.outdir + '/frac_v_abscosth.pdf')
gFile.Close()


for v in ['fBkg', 'f1S', 'f2S']:
    printVar(ws, v, abscosthRan)


getNEvents(ws, abscosthRan)

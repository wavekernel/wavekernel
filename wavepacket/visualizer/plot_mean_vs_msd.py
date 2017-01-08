# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, pylab
kAuPerAngstrom = 1.8897259885789
kAngstromPerAu = kAuPerAngstrom ** -1.0
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0

def plot(wavepacket_out, wavepacket_out_path, highlight_i, mean_min, mean_max, msd_min, msd_max, energy_min, energy_max, is_log, to_label, to_print, title, out_filename):
    condition = wavepacket_out['condition']
    fst_filter = wavepacket_out['setting']['fst_filter']
    end_filter = wavepacket_out['setting']['end_filter']
    indices = range(fst_filter, end_filter + 1)
    eigenvalues = condition['eigenvalues']
    eigenstate_msds = condition['eigenstate_msd_total']
    eigenstate_msds = map(lambda x: x * kAngstrom2PerAu2, eigenstate_msds)
    eigenstate_means = condition['eigenstate_mean_z']
    eigenstate_means = map(lambda x: x * kAngstromPerAu, eigenstate_means)

    if energy_min is None:
        energy_min = min(eigenvalues)
    if energy_max is None:
        energy_max = max(eigenvalues)
    if is_log:
        pylab.yscale('log')

    data_preplot = []
    for (energy, mean, msd) in zip(eigenvalues, eigenstate_means, eigenstate_msds):
        if energy_min <= energy and energy <= energy_max:
            data_preplot.append((mean, msd))
    data_preplot_unzip = zip(*data_preplot)
    eigenstate_means_preplot = data_preplot_unzip[0]
    eigenstate_msds_preplot = data_preplot_unzip[1]

    # Get limit values.
    pylab.plot(eigenstate_means_preplot, eigenstate_msds_preplot, 'o')
    if mean_min is None:
        mean_min = pylab.xlim()[0]
    if mean_max is None:
        mean_max = pylab.xlim()[1]
    if msd_min is None:
        msd_min = pylab.ylim()[0]
    if msd_max is None:
        msd_max = pylab.ylim()[1]

    # Recollect data in the new limit and reset the graph.
    data_mainplot = []
    if to_print:
        print '# Index, Energy, Mean Z, MSD'
    for (i, energy, mean, msd) in zip(indices, eigenvalues, eigenstate_means, eigenstate_msds):
        if energy_min <= energy and energy <= energy_max and \
           mean_min <= mean and mean <= mean_max and \
           msd_min <= msd and msd <= msd_max:
            data_mainplot.append((i, energy, mean, msd))
            if to_print:
                print i, energy, mean, msd
    assert(data_mainplot != [])
    pylab.clf()

    if is_log:
        pylab.yscale('log')
    pylab.xlim(mean_min, mean_max)
    pylab.ylim(msd_min, msd_max)
    xticks_new = list(pylab.xticks()[0])
    xticks_new.extend([mean_min, mean_max])
    pylab.xticks(xticks_new)
    pylab.xlim(mean_min, mean_max)  # limit setting again is needed.
    yticks_new = list(pylab.yticks()[0])
    yticks_new.extend([msd_min, msd_max])
    pylab.yticks(yticks_new)
    pylab.ylim(msd_min, msd_max)  # limit setting again is needed.

    data_mainplot_unzip = zip(*data_mainplot)
    eigenstate_means_mainplot = data_mainplot_unzip[2]
    eigenstate_msds_mainplot = data_mainplot_unzip[3]
    pylab.plot(eigenstate_means_mainplot, eigenstate_msds_mainplot, 'o', label=str(len(data_mainplot)) + ' points')

    if to_label:
        for (i, energy, mean, msd) in data_mainplot:
            pylab.text(mean, msd, str(i))

    if not highlight_i is None:
        data_mainplot_highlighted = filter(lambda d: d[0] == highlight_i, data_mainplot)
        means_highlighted = map(lambda d: d[2], data_mainplot_highlighted)
        msds_highlighted = map(lambda d: d[3], data_mainplot_highlighted)
        pylab.plot(means_highlighted, msds_highlighted, 'o', color='red', label=str(highlight_i), markersize=10)
    pylab.legend(numpoints=1)

    if title is None:
        title = wavepacket_out_path
    pylab.xlabel('Mean Z [$\AA$]')
    pylab.ylabel('MSD [$\AA^2$]')
    pylab.title(title)

    if out_filename is None:
        out_filename = re.sub("\.[^.]+$", "", wavepacket_out_path) + "_mean_vs_msd.png"
    pylab.savefig(out_filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavepacket_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-i', metavar='NUM', dest='highlight_i', type=int, default=None,
                        help='')
    parser.add_argument('-l', action='store_true', dest='to_label',
                        default=False, help='')
    parser.add_argument('-p', action='store_true', dest='to_print',
                        default=False, help='')
    parser.add_argument('--mean-min', metavar='MIN', dest='mean_min', type=float, default=None,
                        help='')
    parser.add_argument('--mean-max', metavar='MAX', dest='mean_max', type=float, default=None,
                        help='')
    parser.add_argument('--msd-min', metavar='MIN', dest='msd_min', type=float, default=None,
                        help='')
    parser.add_argument('--msd-max', metavar='MAX', dest='msd_max', type=float, default=None,
                        help='')
    parser.add_argument('--energy-min', metavar='MIN', dest='energy_min', type=float, default=None,
                        help='')
    parser.add_argument('--energy-max', metavar='MAX', dest='energy_max', type=float, default=None,
                        help='')
    parser.add_argument('--linear', action='store_false', dest='is_log',
                        default=True, help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', default=None,
                        help='')
    parser.add_argument('-o', metavar='FILE', dest='out_filename', default=None,
                        help='')
    args = parser.parse_args()

    with open(args.wavepacket_out_path, 'r') as fp:
        wavepacket_out = json.load(fp)
    plot(wavepacket_out, args.wavepacket_out_path, args.highlight_i,
         args.mean_min, args.mean_max, args.msd_min, args.msd_max, args.energy_min, args.energy_max,
         args.is_log, args.to_label, args.to_print, args.title, args.out_filename)

# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, pylab
kAuPerAngstrom = 1.8897259885789
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0
kAngstromPerAu = kAuPerAngstrom ** -1.0

def plot(wavekernel_out, wavekernel_out_path, to_show_msd, highlight_i, energy_min, energy_max, ymin, ymax, is_log, to_label, title, out_filename):
    condition = wavekernel_out['condition']
    eigenvalues = condition['eigenvalues']
    eigenstate_msd = condition['eigenstate_mean_z']
    eigenstate_msd = map(lambda x: x * kAngstromPerAu, eigenstate_msd)
    fst_filter = wavekernel_out['setting']['fst_filter']

    pylab.figure(figsize=(16, 12))

    if is_log:
        pylab.yscale('log')
    pylab.plot(eigenvalues, eigenstate_msd, 'o')

    if not highlight_i is None:
        j = highlight_i - fst_filter
        xs = [eigenvalues[j]]
        ys = [eigenstate_msd[j]]
        pylab.plot(xs, ys, 'o', color='red', label=str(highlight_i), markersize=10)
        pylab.legend(numpoints=1)

    if energy_min is None:
        energy_min = pylab.xlim()[0]
    if energy_max is None:
        energy_max = pylab.xlim()[1]
    if ymin is None:
        ymin = pylab.ylim()[0]
    if ymax is None:
        ymax = pylab.ylim()[1]
    pylab.xlim(energy_min, energy_max)
    pylab.ylim(ymin, ymax)

    xticks_new = list(pylab.xticks()[0])
    xticks_new.extend([energy_min, energy_max])
    pylab.xticks(xticks_new)
    pylab.xlim(energy_min, energy_max)  # limit setting again is needed.
    yticks_new = list(pylab.yticks()[0])
    yticks_new.extend([ymin, ymax])
    pylab.yticks(yticks_new)
    pylab.ylim(ymin, ymax)  # limit setting again is needed.

    pylab.rcParams.update({'font.size': 10})
    if to_label:
        j = fst_filter
        for x, y in zip(eigenvalues, eigenstate_msd):
            if energy_min <= x <= energy_max and ymin <= y <= ymax:
                pylab.text(x, y, str(j))
            j += 1

#    num_points = len(filter(
#        lambda (x, y): energy_min <= x and x <= energy_max and ymin <= y and y <= ymax,
#        zip(eigenvalues, eigenstate_msd)))
#
    if title is None:
        title = wavekernel_out_path
    pylab.xlabel('Energy [a.u.]')
    pylab.ylabel('Mean x [$\AA$]')
    pylab.title(title)

    if out_filename is None:
        out_filename = re.sub("\.[^.]+$", "", wavekernel_out_path) + "_eigenvalue_vs_mean.png"
    pylab.savefig(out_filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-m', action='store_true', dest='to_show_msd',
                        default=True, help='')
    parser.add_argument('-x', action='store_false', dest='to_show_msd',
                        default=True, help='')
    parser.add_argument('-i', metavar='NUM', dest='highlight_i', type=int, default=None,
                        help='')
    parser.add_argument('-l', action='store_true', dest='to_label',
                        default=False, help='')
    parser.add_argument('--energy-min', metavar='MIN', dest='energy_min', type=float, default=None,
                        help='')
    parser.add_argument('--energy-max', metavar='MAX', dest='energy_max', type=float, default=None,
                        help='')
    parser.add_argument('--ymin', metavar='MIN', dest='ymin', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='MAX', dest='ymax', type=float, default=None,
                        help='')
    parser.add_argument('--linear', action='store_false', dest='is_log',
                        default=True, help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', default=None,
                        help='')
    parser.add_argument('-o', metavar='FILE', dest='out_filename', default=None,
                        help='')
    args = parser.parse_args()

    if not args.to_show_msd:
        args.is_log = False

    with open(args.wavekernel_out_path, 'r') as fp:
        wavekernel_out = json.load(fp)
    plot(wavekernel_out, args.wavekernel_out_path, args.to_show_msd, args.highlight_i, args.energy_min, args.energy_max, args.ymin, args.ymax, args.is_log, args.to_label, args.title, args.out_filename)

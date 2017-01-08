# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, copy, pylab, math

kAuPerAngstrom = 1.8897259885789
kPsecPerAu = 2.418884326505e-5  # Time
kJoulePerEh = 4.35974434e-18  # [J / a.u. (energy)]
kBoltzmannInSI = 1.3806488e-23  # [J / K]
kAuPerKelvin = kBoltzmannInSI / kJoulePerEh

def convert(vss, temperature):  # temperature: float in Kelvein or None
    ts = []
    ys = []
    if temperature is None:
        def get_factor(eigenvalue_diff):
            return 1.
    else:
        def get_factor(eigenvalue_diff):
            return math.exp(- eigenvalue_diff / (temperature * kAuPerKelvin))
    for vs in vss:
        ts.append(vs['time'] * kPsecPerAu)
        eigenvalue_homo = vs['eigenvalues'][-1]
        sum_pr = 0.
        sum_factor = 0.
        for eigenvalue, ipratio in zip(vs['eigenvalues'], vs['ipratios']):
            factor = get_factor(eigenvalue_homo - eigenvalue)
            sum_pr += factor / ipratio
            sum_factor += factor
        ys.append(sum_pr / sum_factor)
    return ts, ys

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('eigenvalue_vs_ipratio_paths', metavar='JSON', type=str, nargs='+',
                        help='')
    #parser.add_argument('--xmin', metavar='XMIN', dest='x_min', type=float, default=None,
    #                    help='')
    #parser.add_argument('--xmax', metavar='XMAX', dest='x_max', type=float, default=None,
    #                    help='')
    parser.add_argument('--ymin', metavar='YMIN', dest='y_min', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='YMAX', dest='y_max', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default='',
                        help='')
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    parser.add_argument('-l', metavar='LABELS', dest='labels', type=str, default=None,
                        help='')
    parser.add_argument('--temperature', metavar='TEMPERATURE', dest='temperature', type=float, default=None,
                        help='')  # in Kelvin.
    args = parser.parse_args()

    if args.title == '':
        title = args.eigenvalue_vs_ipratio_paths[0]
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub(r'(_eigenvalue_vs_ipratio)?\.[^.]+$', '', args.eigenvalue_vs_ipratio_paths[0]) + \
                   '_eigenvalue_vs_ipratio_temporal.png'
    else:
        fig_path = args.fig_path

    if args.labels is None:
        labels = map(lambda i: 'sample %d' % i, range(len(args.eigenvalue_vs_ipratio_paths)))
    else:
        labels = args.labels.split(',')

    vsss = []
    for path in args.eigenvalue_vs_ipratio_paths:
        with open(path, 'r') as fp:
            vsss.append(json.load(fp))

    header = re.sub('\.[^.]+$', '', fig_path)

    min_eigenvalue = 1e100
    min_ipratio = 1e100
    max_eigenvalue = -1e100
    max_ipratio = -1e100
    for vss in vsss:
        for vs in vss:
            min_eigenvalue = min(min_eigenvalue, min(vs['eigenvalues']))
            min_ipratio = min(min_ipratio, min(vs['ipratios']))
            max_eigenvalue = max(max_eigenvalue, max(vs['eigenvalues']))
            max_ipratio = max(max_ipratio, max(vs['ipratios']))

    if args.time_start is None:
        time_start = None #vsss[0][0]['time']
    else:
        time_start = args.time_start
    if args.time_end is None:
        time_end = None #vsss[0][-1]['time']
    else:
        time_end = args.time_end

    if args.y_min is None:
        y_min = None
    else:
        y_min = args.y_min
    if args.y_max is None:
        y_max = None
    else:
        y_max = args.y_max

    output_filename = '%s.png' % (header)
    pylab.title(title)
    pylab.xlabel('Time [ps]')
    pylab.ylabel('Weighted average pratio in group')
    pylab.grid()
    pylab.xlim(time_start, time_end)
    pylab.ylim(y_min, y_max)        
    for i, vss in enumerate(vsss):
        ts, ys = convert(vss, args.temperature)
        pylab.plot(ts, ys, '+-', label=labels[i])
    pylab.legend(loc='upper right')
    pylab.savefig(output_filename)


        

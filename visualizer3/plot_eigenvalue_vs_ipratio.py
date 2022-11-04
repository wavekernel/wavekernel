# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, copy, pylab

kAuPerAngstrom = 1.8897259885789
kPsecPerAu = 2.418884326505e-5  # Time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('eigenvalue_vs_ipratio_paths', metavar='JSON', type=str, nargs='+',
                        help='')
    parser.add_argument('--xmin', metavar='XMIN', dest='x_min', type=float, default=None,
                        help='')
    parser.add_argument('--xmax', metavar='XMAX', dest='x_max', type=float, default=None,
                        help='')
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
    args = parser.parse_args()

    if args.title == '':
        title = args.eigenvalue_vs_ipratio_paths[0]
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub(r'(_eigenvalue_vs_ipratio)?\.[^.]+$', '', args.eigenvalue_vs_ipratio_paths[0]) + \
                   '_eigenvalue_vs_ipratio.png'
    else:
        fig_path = args.fig_path

    if args.labels is None:
        labels = ['sample %d' % i for i in range(len(args.eigenvalue_vs_ipratio_paths))]
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

    if args.x_min is None:
        x_min = min_eigenvalue
    else:
        x_min = args.x_min
    if args.x_max is None:
        x_max = max_eigenvalue
    else:
        x_max = args.x_max
    if args.y_min is None:
        y_min = 0.
    else:
        y_min = args.y_min
    if args.y_max is None:
        y_max = 1. / min_ipratio
    else:
        y_max = args.y_max

    for step in range(len(vsss[0])):
        output_filename = '%s_%06d.png' % (header, step)
        pylab.title(title + ' %.4f [ps]' % (vsss[0][step]['time'] * kPsecPerAu))
        pylab.xlabel('Eigenvalue [au]')
        pylab.ylabel('Pratio')
        pylab.grid()
        pylab.xlim(x_min, x_max)
        pylab.ylim(y_min, y_max)        
        for sample in range(len(vsss)):
            pratios = [1. / ipratio for ipratio in vsss[sample][step]['ipratios']]
            pylab.plot(vsss[sample][step]['eigenvalues'], pratios, '+-', label=labels[sample])
        pylab.legend(loc='upper left')
        pylab.savefig(output_filename)
        pylab.clf()

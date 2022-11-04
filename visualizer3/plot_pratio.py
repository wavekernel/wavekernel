# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def plot(pratio_extract, printed_value, ymin, ymax, time_start, time_end, title, fig_path):
    pylab.rcParams['font.size'] = 16
    fig = pylab.figure()
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)

    ts = [t * kPsecPerAu for t in pratio_extract['ts']]

    if time_start is None:
        time_start = min(ts)
    if time_end is None:
        time_end = max(ts)
    pylab.xlim(time_start, time_end)

    pylab.xlabel('Time [ps]')
    pylab.ylabel('Pratio')

    pylab.grid(True)

    mark_size = 7
    if printed_value == 'alpha':
        ys = pratio_extract['alpha_pratios']
    elif printed_value == 'psi':
        ys = pratio_extract['psi_pratios']
    else:
        assert(False)

    if ymin is None:
        ymin = 1.0
    if ymax is None:
        ymax = max(ys)
    pylab.ylim(ymin, ymax)

    pylab.plot(ts, ys, '+-', label=printed_value, ms=mark_size)

    pylab.legend(loc='upper left')
    pylab.title(title)
    pylab.savefig(fig_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('pratio_extract_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-p', metavar='VALUE', dest='printed_value', type=str, default='alpha',
                        help='')
    parser.add_argument('--ymin', metavar='Y_MIN', dest='ymin', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='Y_MAX', dest='ymax', type=float, default=None,
                        help='')
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default='',
                        help='')
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    args = parser.parse_args()

    if args.title == "":
        title = args.pratio_extract_path
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub('(_pratio)?\.[^.]+$', '', args.pratio_extract_path) + '_pratio.png'
    else:
        fig_path = args.fig_path

    with open(args.pratio_extract_path, 'r') as fp:
        pratio_extract = json.load(fp)
    plot(pratio_extract, args.printed_value, args.ymin, args.ymax, args.time_start, args.time_end, title, fig_path)

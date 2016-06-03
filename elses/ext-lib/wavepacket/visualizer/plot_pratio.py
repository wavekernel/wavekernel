# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def plot(pratio_extract, time_start, time_end, title, fig_path):
    pylab.rcParams['font.size'] = 16

    ts = map(lambda t: t * kPsecPerAu, pratio_extract['ts'])

    if time_start is None:
        time_start = min(ts)
    if time_end is None:
        time_end = max(ts)
    pylab.xlim(time_start, time_end)

    pylab.xlabel('Time [ps]')
    pylab.ylabel('Pratio')

    mark_size = 7
    #pylab.plot(ts, pratio_extract['psi_pratios'], '+', label='psi', ms=mark_size)
    pylab.plot(ts, pratio_extract['alpha_pratios'], '+', label='alpha', ms=mark_size)

    pylab.legend(loc='upper left')
    pylab.title(title)
    pylab.savefig(fig_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('pratio_extract_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
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
    plot(pratio_extract, args.time_start, args.time_end, title, fig_path)

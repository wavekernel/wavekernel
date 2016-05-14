# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def plot(pratio_extract, title, fig_path):
    pylab.rcParams['font.size'] = 16

    ts = map(lambda t: t * kPsecPerAu, pratio_extract['ts'])
    
    pylab.xlabel('Time [ps]')
    pylab.ylabel('Pratio')

    mark_size = 7
    pylab.plot(ts, pratio_extract['psi_pratios'], '+', label='psi', ms=mark_size)
    pylab.plot(ts, pratio_extract['alpha_pratios'], '+', label='alpha', ms=mark_size)

    pylab.legend(loc='upper left')
    pylab.title(title)
    pylab.savefig(fig_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('pratio_extract_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')    
    args = parser.parse_args()

    if args.title == "":
        title = args.pratio_extract_path
    else:
        title = args.title

    fig_path = re.sub('(_pratio)?\.[^.]+$', '', args.pratio_extract_path) + '_pratio.png'
    with open(args.pratio_extract_path, 'r') as fp:
        pratio_extract = json.load(fp)
    plot(pratio_extract, title, fig_path)

# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def plot(energy_calc, title, fig_path):
    ts = [s['time'] * kPsecPerAu for s in energy_calc['states']]
    energy_devs = [s['TB_energy_deviation'] for s in energy_calc['states']]

    mark_size = 5
    pylab.plot(ts, energy_devs, '+', label='energy_dev', ms=mark_size)

    pylab.legend(loc='upper left')
    #pylab.subplots_adjust(right=0.65)
    #pylab.legend(loc='upper right')
    pylab.title(title)
    pylab.xlabel('time [ps]')
    pylab.ylabel('energy deviation [a.u.]')
    pylab.savefig(fig_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('energy_calc_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('--energy-min', metavar='ENERGY_MIN', dest='energy_min',
                        type=float, default=None, help='')
    parser.add_argument('--energy-max', metavar='ENERGY_MAX', dest='energy_max',
                        type=float, default=None, help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    args = parser.parse_args()

    fig_path = re.sub('(_energy_dev)?\.[^.]+$', '', args.energy_calc_path) + '_energy_dev.png'
    with open(args.energy_calc_path, 'r') as fp:
        energy_calc = json.load(fp)

    title = args.energy_calc_path
    plot(energy_calc, title, fig_path)

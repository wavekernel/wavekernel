# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def read_plain_extracted(fp):
    rf = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    ts = []
    isamrs = []
    tb_energy = []
    nl_energy = []
    total_energy = []
    for line in fp:
        m = re.search(r'(?P<t>%s)\s+(?P<i>\d+)\s+(?P<isamr>True|False)\s+(?P<e1>%s)\s(?P<e2>%s)\s(?P<e3>%s)\s(?P<mean>%s)\s+(?P<msd>%s)\s+(?P<pip>%s)\s+(?P<aip>%s)' %
                      (rf, rf, rf, rf, rf, rf, rf, rf), line)
        if m:
            t = float(m.group('t'))
            isamr = m.group('isamr') == 'True'
            tbe = float(m.group('e1'))
            nle = float(m.group('e2'))
            totale = float(m.group('e3'))
            ts.append(t)
            isamrs.append(isamr)
            tb_energy.append(tbe)
            nl_energy.append(nle)
            total_energy.append(totale)

    #tb_energy = map(lambda e: e - tb_energy[0], tb_energy)
    #nl_energy = map(lambda e: e - nl_energy[0], nl_energy)
    #total_energy = map(lambda e: e - total_energy[0], total_energy)

    minimum = min(tb_energy)
    minimum = min(minimum, min(nl_energy))
    minimum = min(minimum, min(total_energy))
    maximum = max(tb_energy)
    maximum = max(maximum, max(nl_energy))
    maximum = max(maximum, max(total_energy))
    return {'ts': ts, 'isamrs': isamrs, 'max_eigenvalue': maximum, 'min_eigenvalue': minimum, 'tb_energy': tb_energy, 'nl_energy': nl_energy, 'total_energy': total_energy}

def plot(energy_calc, energy_min, energy_max, title, fig_path):
    if max(energy_calc['ts']) * kPsecPerAu < 0.001:
        ts = map(lambda t: t * kPsecPerAu * 1000, energy_calc['ts'])
        pylab.xlabel('Time [fs]')
    else:
        ts = map(lambda t: t * kPsecPerAu, energy_calc['ts'])
        pylab.xlabel('Time [ps]')

    min_eigenvalue = energy_calc['min_eigenvalue']
    max_eigenvalue = energy_calc['max_eigenvalue']
    diff_eigenvalue = max_eigenvalue - min_eigenvalue
    if energy_min is None:
        energy_min = min_eigenvalue - diff_eigenvalue * 0.01
    if energy_max is None:
        energy_max = max_eigenvalue + diff_eigenvalue * 0.01
    pylab.ylim(energy_min, energy_max)
    yticks_new = list(pylab.yticks()[0])
    yticks_new.extend([energy_min, energy_max])
    pylab.yticks(yticks_new)
    pylab.ylim(energy_min, energy_max)  # limit setting again is needed.
    pylab.ylabel('Energy [a.u.]')

    mark_size = 3
    pylab.plot(ts, energy_calc['tb_energy'], '+', label='TB_energy (H0)', ms=mark_size)
    pylab.plot(ts, energy_calc['nl_energy'], '+', label='NL_energy (H1)', ms=mark_size)
    pylab.plot(ts, energy_calc['total_energy'], '+', label='total_energy\n(H0 + H1)', ms=mark_size)

    line_min_eigenvalue = [energy_calc['min_eigenvalue']] * len(ts)
    line_max_eigenvalue = [energy_calc['max_eigenvalue']] * len(ts)
    pylab.plot(ts, line_min_eigenvalue, label='min energy')
    pylab.plot(ts, line_max_eigenvalue, label='max energy')

    pylab.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    pylab.subplots_adjust(right=0.65)
    #pylab.legend(loc='upper right')
    pylab.title(title)
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
    parser.add_argument('--plain', action='store_true', dest='is_plain_extracted_mode',
                        default=False, help='')
    args = parser.parse_args()

    if args.title == "":
        title = args.energy_calc_path
    else:
        title = args.title

    fig_path = re.sub('(_energy)?\.[^.]+$', '', args.energy_calc_path) + '_energy.png'
    with open(args.energy_calc_path, 'r') as fp:
        if args.is_plain_extracted_mode:
            energy_calc = read_plain_extracted(fp)
        else:
            energy_calc = json.load(fp)
    plot(energy_calc, args.energy_min, args.energy_max, title, fig_path)

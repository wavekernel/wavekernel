# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time

def print_energy(ts, xs, filename):
    with open(filename, 'w') as fp:
        fp.write('# time[ps]\tenergy[au]\n')
        for t, x in zip(ts, xs):
            fp.write('%f\t%f\n' % (t, x))

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

def diff_list(xs):
    n = len(xs)
    ys = [0.] * n
    for i in range(n - 1):
        ys[i + 1] = xs[i + 1] - xs[i]
    return ys

def plot(energy_calc, time_start, time_end, energy_min, energy_max, is_diff_mode, title, fig_path):
    if max(energy_calc['ts']) * kPsecPerAu < 0.001:
        ts = map(lambda t: t * kPsecPerAu * 1000, energy_calc['ts'])
        pylab.xlabel('Time [fs]')
    else:
        ts = map(lambda t: t * kPsecPerAu, energy_calc['ts'])
        pylab.xlabel('Time [ps]')

    if time_start is None:
        time_start = min(ts)
    if time_end is None:
        time_end = max(ts)
    pylab.xlim(time_start, time_end)

    pylab.ylabel('Energy [a.u.]')

    mark_size = 5
    ys = energy_calc['tb_energy']
    if is_diff_mode:
        ys = diff_list(ys)
    print_energy(ts, ys, 'tb_energy_%s.txt' % fig_path)
    pylab.plot(ts, ys, '+-', label='TB_energy (H0)', ms=mark_size)
    pylab.plot(ts, energy_calc['nl_energy'], '+', label='NL_energy (H1)', ms=mark_size)
    pylab.plot(ts, energy_calc['total_energy'], '+', label='total_energy\n(H0 + H1)', ms=mark_size)
    if 'eigenvalues_log' in energy_calc:
        ys = map(lambda l: l['eigenvalues'][0], energy_calc['eigenvalues_log'])
        if is_diff_mode:
            ys = diff_list(ys)
        print_energy(ts, ys, 'homo_energy_%s.txt' % fig_path)
        pylab.plot(ts, ys, 'x-', label='HOMO', ms=mark_size, color='red')

    pylab.ylim(energy_min, energy_max)
    pylab.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    pylab.grid(True)
    pylab.subplots_adjust(right=0.65)
    #pylab.legend(loc='upper right')
    pylab.title(title)
    pylab.savefig(fig_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('energy_calc_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')
    parser.add_argument('--energy-min', metavar='ENERGY_MIN', dest='energy_min',
                        type=float, default=None, help='')
    parser.add_argument('--energy-max', metavar='ENERGY_MAX', dest='energy_max',
                        type=float, default=None, help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    parser.add_argument('--plain', action='store_true', dest='is_plain_extracted_mode',
                        default=False, help='')
    parser.add_argument('--diff', action='store_true', dest='is_diff_mode',
                        default=False, help='')
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    args = parser.parse_args()

    if args.title == "":
        title = args.energy_calc_path
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub('(_energy)?\.[^.]+$', '', args.energy_calc_path)
        if args.is_diff_mode:
            fig_path += '_energy_diff.png'
        else:
            fig_path += '_energy.png'
    else:
        fig_path = args.fig_path

    with open(args.energy_calc_path, 'r') as fp:
        if args.is_plain_extracted_mode:
            energy_calc = read_plain_extracted(fp)
        else:
            energy_calc = json.load(fp)
    plot(energy_calc, args.time_start, args.time_end, args.energy_min, args.energy_max,
         args.is_diff_mode, title, fig_path)

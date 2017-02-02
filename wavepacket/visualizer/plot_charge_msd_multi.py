# -*- coding: utf-8 -*-
import argparse, json, sys, re, pylab, os.path, matplotlib

kAngstromPerAu = 0.529177249  # Length
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0  # Length
kPsecPerAu = 2.418884326505e-5  # Time
kBoltzmann = 8.6173303e-5  # [eV / K]

def avg(xss):
    m = len(xss)
    n = len(xss[0])
    xs = [0.] * n
    for j in range(n):
        for i in range(m):
            xs[j] += xss[i][j]
        xs[j] /= m
    return xs

def print_avged(ts, xs, filename):
    with open(filename, 'w') as fp:
        fp.write('# time[ps]\tMSD[Angstrom^2]\n')
        for t, x in zip(ts, xs):
            fp.write('%f\t%f\n' % (t, x))

def read_plain_extracted(fp):
    rf = r'[+-]?(\d+\.)?\d+([deE][+-]?\d+)?'
    ts = []
    means = []
    msds = []
    for line in fp:
        m = re.search(r'(?P<t>%s)\s+(?P<i>\d+)\s+(?P<isamr>True|False)\s+(?P<e1>%s)\s(?P<e2>%s)\s(?P<e3>%s)\s(?P<mean>%s)\s+(?P<msd>%s)\s+(?P<pip>%s)\s+(?P<aip>%s)' %
                      (rf, rf, rf, rf, rf, rf, rf, rf), line)
        if m:
            t = float(m.group('t'))
            mean = float(m.group('mean'))
            msd = float(m.group('msd'))
            ts.append(t)
            msds.append([0.0, 0.0, 0.0, msd])  # x, y, z, total
            means.append([mean, 0.0, 0.0])  # x, y, z
    return {'ts': ts, 'msds': msds, 'means': means}  # Note: TB_energy_deviations is not supported.

def plot_charge_moment(grouped_charge_moments,
                       msd_axis, msd_min, msd_max,
                       energy_min, energy_max, to_plot_tb_energy_deviation,
                       time_start, time_end,
                       time_start_diffusion, time_end_diffusion,
                       grouped_labels, group_to_avg_label, to_plot_avg_only,
                       title, fig_path):
    #font = {'size': 20}
    #matplotlib.rc('font', **font)
    fig = pylab.figure(figsize=(10, 7.5))
    # Cancel axis offset.
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)
    axis_name_to_num = {'x': 0, 'y': 1, 'z': 2, 'total': 3}
    msd_axis_num = axis_name_to_num[msd_axis]
    #mean_axis_num = axis_name_to_num[mean_axis]
    to_show_diffusion_coef = True

    group_to_tss = map(lambda cms:
                       map(lambda cm: map(lambda t: t * kPsecPerAu, cm['ts']), cms),
                       grouped_charge_moments)
    group_to_msdss = map(lambda cms:
                        map(lambda cm: map(lambda x: x[msd_axis_num] * kAngstrom2PerAu2, cm['msds']), cms),
                        grouped_charge_moments)
    #means = map(lambda m: m[mean_axis_num] * kAngstromPerAu, charge_moment['means'])
    #tb_energy_deviations = charge_moment['tb_energy_deviations']

    if time_start_diffusion is None:
        time_start_diffusion = min(group_to_tss[0][0])
    if time_end_diffusion is None:
        time_end_diffusion = max(group_to_tss[0][0])

    group_to_avg_msds = map(avg, group_to_msdss)
    #for g, avg_msds in enumerate(group_to_avg_msds):
    #    print_avged(ts, avg_msds, 'avg_%02d.txt' % g)  # ps, Angstrom^2

    group_to_diffusion_coef_info = []
    if to_show_diffusion_coef:
        for tss, avg_msds, avg_label in zip(group_to_tss, group_to_avg_msds, group_to_avg_label):
            xys = zip(tss[0], avg_msds)
            xys = filter(lambda xy: time_start_diffusion <= xy[0] and xy[0] < time_end_diffusion, xys)
            n = len(xys)
            x_bar = y_bar = xy_bar = xx_bar = 0.0
            for (x, y) in xys:
                x_bar += x
                y_bar += y
                xy_bar += x * y
                xx_bar += x ** 2.0
            x_bar /= float(n)
            y_bar /= float(n)
            xy_bar /= float(n)
            xx_bar /= float(n)
            a = (xy_bar - x_bar * y_bar) / (xx_bar - x_bar ** 2.0)
            b = (-x_bar * xy_bar + xx_bar * y_bar) / (xx_bar - x_bar ** 2.0)
            ts_new = map(lambda xy: xy[0], xys)
            ys_hat = map(lambda t: a * t + b, ts_new)
            intercept_relative_error_left = abs((b - xys[0][1]) / xys[0][1])
            intercept_relative_error_right = abs((a * xys[-1][0] + b - xys[-1][1]) / xys[-1][1])
            rmse = pylab.sqrt(sum(map(lambda xy: (xy[1] - a * xy[0] - b) ** 2.0, xys)) / float(n))
            print 'T0, T1: ', xys[0][0], xys[-1][0]
            print 'a T0 + b, y(T0), IRE(0): ', b, xys[0][1], intercept_relative_error_left
            print 'a T1 + b, y(T1), IRE(1): ', \
                a * xys[-1][0] + b, xys[0][1], intercept_relative_error_right
            print 'RMSE [angstrom^2]: ', rmse
            print 'diffusion coefficient [angstrom^2 / ps]: ', a / 2.0
            print 'diffusion coefficient [cm^2 / s]: ', a / 2.0 * 1e-4
            # '1.0' = [e]
            # 'kBoltzmann * 1.0' = [V / K]
            mbt = a / 2.0 * 1e-4 / (kBoltzmann * 1.0)
            print 'mobility * temperature [cm^2 K / V s]', mbt
            group_to_diffusion_coef_info.append({'mbt': mbt,
                                                 'IRE0': intercept_relative_error_left,
                                                 'IRE1': intercept_relative_error_right,
                                                 'RMSE': rmse})
            #k, s, d = fig_path_to_params(fig_path)
            #print 'ZZZ', fig_path, k, s, d, mbt
            pylab.plot(ts_new[0:n:n-1], ys_hat[0:n:n-1], 'x-', markersize=10, label='%s fit' % avg_label)

    pylab.xlabel('Time [ps]')
    pylab.ylabel('MSD ' + msd_axis + ' [$\AA^2$]', color='blue')
    pylab.grid(True)

    if not to_plot_avg_only:
        for tss, msdss, labels in zip(group_to_tss, group_to_msdss, grouped_labels):
            for ts, msds, label in zip(tss, msdss, labels):
                pylab.plot(ts, msds, '+-', label=label, markerfacecolor='none')

    for tss, avg_msds, avg_label in zip(group_to_tss, group_to_avg_msds, group_to_avg_label):
        pylab.plot(tss[0], avg_msds, '+-', label=avg_label,
                   markerfacecolor='none', markersize=10)  #, markeredgecolor='red'
    pylab.legend()

    if msd_min is None:
        msd_min = pylab.ylim()[0]
    if msd_max is None:
        msd_max = pylab.ylim()[1]
    pylab.ylim(msd_min, msd_max)
    #yticks_new = list(pylab.yticks()[0])
    #yticks_new.extend([msd_min, msd_max])
    #pylab.yticks(yticks_new)
    #pylab.ylim(msd_min, msd_max)  # limit setting again is needed.

    if time_start is None:
        time_start = pylab.xlim()[0]
    if time_end is None:
        time_end = pylab.xlim()[1]

    if to_show_diffusion_coef:
        for g, diffusion_coef_info in enumerate(group_to_diffusion_coef_info):
            pylab.text(time_start + (time_end - time_start) * 0.01,
                       msd_min + (msd_max - msd_min) * (0.94 - 0.04 * g),
                       str(diffusion_coef_info['mbt']) + ' [cm^2 K / V s]')
            #pylab.text(time_start + (time_end - time_start) * 0.01,
            #           msd_min + (msd_max - msd_min) * 0.9,
            #           'IRE(0) ' + str(intercept_relative_error_left))
            #pylab.text(time_start + (time_end - time_start) * 0.01,
            #           msd_min + (msd_max - msd_min) * 0.86,
            #           'IRE(1) ' + str(intercept_relative_error_right))
            #pylab.text(time_start + (time_end - time_start) * 0.01,
            #           msd_min + (msd_max - msd_min) * 0.82,
            #           'RMSE ' + str(rmse) + ' [$\AA^2$]')

    #pylab.twinx()
    ## Cancel axis offset.
    #ax = fig.gca()
    #ax.ticklabel_format(useOffset=False)

    pylab.xlim(time_start, time_end)
    #xticks_new = list(pylab.xticks()[0])
    #xticks_new.extend([time_start, time_end])
    #pylab.xticks(xticks_new)
    #pylab.xlim(time_start, time_end)  # limit setting again is needed.

    pylab.title(title)
    pylab.savefig(fig_path, dpi=80)  # dpi=80 correspond to figsize=(10, 7.5).

def parse_charge_moment_calc_paths(ss):
    # ['a1.json', 'a2.json', ':', 'b.json'] -> [['a1.json', 'a2.json'], ['b.json']]
    num_groups = len(filter(lambda s: s == ':', ss)) + 1
    grouped_paths = []
    for g in range(num_groups):
        grouped_paths.append([])
    g = 0
    num_paths = 0
    for s in ss:
        if s == ':':
            g += 1
        else:
            grouped_paths[g].append(s)
            num_paths += 1
    return num_paths, grouped_paths

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('charge_moment_calc_paths', metavar='JSON', type=str, nargs='+',
                        help='')
    parser.add_argument('--msd-axis', metavar='MSD_AXIS', dest='msd_axis', type=str, default='total',
                        help='')
    parser.add_argument('--msd-min', metavar='MSD_MIN', dest='msd_min', type=float, default=None,
                        help='')
    parser.add_argument('--msd-max', metavar='MSD_MAX', dest='msd_max', type=float, default=None,
                        help='')
    parser.add_argument('--energy-min', metavar='ENERGY_MIN', dest='energy_min', type=float, default=None,
                        help='')
    parser.add_argument('--energy-max', metavar='ENERGY_MAX', dest='energy_max', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default='',
                        help='')
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-d', action='store_true', dest='to_plot_tb_energy_deviation',
                        default=False, help='')
    parser.add_argument('-l', metavar='LABELS', dest='labels', type=str, default=None,
                        help='labels for each charge moment files separated by comma')
    parser.add_argument('--gl', metavar='GROUP_TO_AVG_LABEL', dest='group_to_avg_label', type=str, default=None,
                        help='labels for each averaged MSD in group separated by colon')
    parser.add_argument('--diffuse-start', metavar='TIME_START_DIFFUSION', dest='time_start_diffusion', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('--diffuse-end', metavar='TIME_END_DIFFUSION', dest='time_end_diffusion', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('--avg-only', action='store_true', dest='to_plot_avg_only',
                        default=False, help='')
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    args = parser.parse_args()

    num_paths, grouped_paths = parse_charge_moment_calc_paths(args.charge_moment_calc_paths)

    if args.title == '':
        title = args.charge_moment_calc_paths[0]
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub(r'(_charge_moment)?\.[^.]+$', '', args.charge_moment_calc_paths[0]) + \
                   '_charge_moment_msd_multi.png'
    else:
        fig_path = args.fig_path

    if args.labels is None:
        grouped_labels = grouped_paths
    else:
        grouped_labels = map(lambda g_str: g_str.split(','), args.labels.split(':'))
        assert(len(grouped_paths) == len(grouped_labels))  # The number of groups must be the same.
        for ps, ls in zip(grouped_paths, grouped_labels):
            assert(len(ps) == len(ls))  # The number of entries must be the same in all the groups.

    if args.group_to_avg_label is None:
        group_to_avg_label = map(lambda i: 'group %d avg' % (i + 1), range(len(grouped_paths)))
    else:
        group_to_avg_label = args.group_to_avg_label.split(':')
        assert(len(grouped_paths) == len(group_to_avg_label))  # The number of groups must be the same.

    grouped_charge_moments = []
    for ps in grouped_paths:
        charge_moments_in_group = []
        for p in ps:
            with open(p, 'r') as fp:
                charge_moments_in_group.append(json.load(fp))
        grouped_charge_moments.append(charge_moments_in_group)

    plot_charge_moment(grouped_charge_moments,
                       args.msd_axis, args.msd_min, args.msd_max,
                       args.energy_min, args.energy_max, args.to_plot_tb_energy_deviation,
                       args.time_start, args.time_end,
                       args.time_start_diffusion, args.time_end_diffusion,
                       grouped_labels, group_to_avg_label, args.to_plot_avg_only, title, fig_path)

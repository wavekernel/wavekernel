# -*- coding: utf-8 -*-
import argparse, json, sys, re, pylab, os.path, matplotlib

kAngstromPerAu = 0.529177249  # Length
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0  # Length
kPsecPerAu = 2.418884326505e-5  # Time
kBoltzmann = 8.6173303e-5  # [eV / K]

def fig_path_to_params(s):  # Temporary information retrieve function.
    rf = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    m = re.search(r'out_(?P<k>(para|meta))_447i_(?P<s>%s)s_(?P<d>%s)d' % (rf, rf), s)
    if m:
        return m.group('k'), float(m.group('s')), float(m.group('d'))
    else:
        assert(False)

def read_plain_extracted(fp):
    rf = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
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

def plot_charge_moment(charge_moment,
                       msd_axis, msd_min, msd_max,
                       mean_axis, mean_min, mean_max,
                       energy_min, energy_max, to_plot_tb_energy_deviation,
                       time_start, time_end,
                       time_start_diffusion, time_end_diffusion,
                       is_raw_mode, title, fig_path):
    #font = {'size': 20}
    #matplotlib.rc('font', **font)
    fig = pylab.figure(figsize=(10, 7.5))
    # Cancel axis offset.
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)
    axis_name_to_num = {'x': 0, 'y': 1, 'z': 2, 'total': 3}
    msd_axis_num = axis_name_to_num[msd_axis]
    mean_axis_num = axis_name_to_num[mean_axis]
    is_fs_mode = max(charge_moment["ts"]) * kPsecPerAu < 0.001
    to_show_diffusion_coef = not is_fs_mode  # fs mode is not supported now.

    if is_fs_mode:
        ts = map(lambda t: t * kPsecPerAu * 1000, charge_moment["ts"])
        pylab.xlabel("Time [fs]")
    else:
        ts = map(lambda t: t * kPsecPerAu, charge_moment["ts"])
        pylab.xlabel("Time [ps]")

    msds = map(lambda x: x[msd_axis_num] * kAngstrom2PerAu2, charge_moment['msds'])
    means = map(lambda m: m[mean_axis_num] * kAngstromPerAu, charge_moment['means'])
    tb_energy_deviations = charge_moment['tb_energy_deviations']

    if time_start_diffusion is None:
        time_start_diffusion = min(ts)
    if time_end_diffusion is None:
        time_end_diffusion = max(ts)

    if to_show_diffusion_coef and not is_fs_mode:
        xys = zip(ts, msds)
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
        print 'RMSE [Å^2]: ', rmse
        print 'diffusion coefficient [Å^2 / ps]: ', a / 2.0
        print 'diffusion coefficient [cm^2 / s]: ', a / 2.0 * 1e-4
        # '1.0' = [e]
        # 'kBoltzmann * 1.0' = [V / K]
        mbt = a / 2.0 * 1e-4 / (kBoltzmann * 1.0)
        print 'mobility * temperature [cm^2 K / V s]', mbt
        #k, s, d = fig_path_to_params(fig_path)
        #print 'ZZZ', fig_path, k, s, d, mbt
        pylab.plot(ts_new[0:n:n-1], ys_hat[0:n:n-1], 'x-', color='green', markersize=10)

    pylab.ylabel('MSD ' + msd_axis + ' [$\AA^2$]', color='blue')
    pylab.grid(True)
    pylab.plot(ts, msds, '+-', markeredgecolor='blue', label='MSD ' + msd_axis,
               markerfacecolor='none')

    if msd_min is None:
        msd_min = pylab.ylim()[0]
    if msd_max is None:
        msd_max = pylab.ylim()[1]
    pylab.ylim(msd_min, msd_max)
    yticks_new = list(pylab.yticks()[0])
    yticks_new.extend([msd_min, msd_max])
    pylab.yticks(yticks_new)
    pylab.ylim(msd_min, msd_max)  # limit setting again is needed.

    if time_start is None:
        time_start = pylab.xlim()[0]
    if time_end is None:
        time_end = pylab.xlim()[1]

    if to_show_diffusion_coef:
        pylab.text(time_start + (time_end - time_start) * 0.01,
                   msd_min + (msd_max - msd_min) * 0.94,
                   str(mbt) + ' [cm^2 K / V s]')
        pylab.text(time_start + (time_end - time_start) * 0.01,
                   msd_min + (msd_max - msd_min) * 0.9,
                   'IRE(0) ' + str(intercept_relative_error_left))
        pylab.text(time_start + (time_end - time_start) * 0.01,
                   msd_min + (msd_max - msd_min) * 0.86,
                   'IRE(1) ' + str(intercept_relative_error_right))
        pylab.text(time_start + (time_end - time_start) * 0.01,
                   msd_min + (msd_max - msd_min) * 0.82,
                   'RMSE ' + str(rmse) + ' [$\AA^2$]')

    pylab.twinx()
    # Cancel axis offset.
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)

    if to_plot_tb_energy_deviation:
        pylab.ylabel('TB energy deviation [a.u.]', color='red')
        pylab.plot(ts, tb_energy_deviations, '+-', color='red', label='TB energy dev')
        if energy_min is None:
            energy_min = pylab.ylim()[0]
        if energy_max is None:
            energy_max = pylab.ylim()[1]
        pylab.ylim(energy_min, energy_max)
        yticks_new = list(pylab.yticks()[0])
        yticks_new.extend([energy_min, energy_max])
        pylab.yticks(yticks_new)
        pylab.ylim(energy_min, energy_max)  # limit setting again is needed.
    else:
        pylab.ylabel('Mean ' + mean_axis + ' [$\AA$]', color='red')
        pylab.plot(ts, means, '+-', color='red', label='Mean ' + mean_axis)
        if mean_min is None:
            mean_min = pylab.ylim()[0]
        if mean_max is None:
            mean_max = pylab.ylim()[1]
        pylab.ylim(mean_min, mean_max)
        yticks_new = list(pylab.yticks()[0])
        yticks_new.extend([mean_min, mean_max])
        pylab.yticks(yticks_new)
        pylab.ylim(mean_min, mean_max)  # limit setting again is needed.

    pylab.xlim(time_start, time_end)
    xticks_new = list(pylab.xticks()[0])
    xticks_new.extend([time_start, time_end])
    pylab.xticks(xticks_new)
    pylab.xlim(time_start, time_end)  # limit setting again is needed.

    pylab.title(title)
    pylab.savefig(fig_path, dpi=80)  # dpi=80 correspond to figsize=(10, 7.5).

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('charge_moment_calc_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('--msd-axis', metavar='MSD_AXIS', dest='msd_axis', type=str, default='total',
                        help='')
    parser.add_argument('--msd-min', metavar='MSD_MIN', dest='msd_min', type=float, default=None,
                        help='')
    parser.add_argument('--msd-max', metavar='MSD_MAX', dest='msd_max', type=float, default=None,
                        help='')
    parser.add_argument('--mean-axis', metavar='MEAN_AXIS', dest='mean_axis', type=str, default='x',
                        help='')
    parser.add_argument('--mean-min', metavar='MEAN_MIN', dest='mean_min', type=float, default=None,
                        help='')
    parser.add_argument('--mean-max', metavar='MEAN_MAX', dest='mean_max', type=float, default=None,
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
    parser.add_argument('--diffuse-start', metavar='TIME_START_DIFFUSION', dest='time_start_diffusion', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('--diffuse-end', metavar='TIME_END_DIFFUSION', dest='time_end_diffusion', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('--plain', action='store_true', dest='is_plain_extracted_mode',
                        default=False, help='')
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    args = parser.parse_args()

    assert(not (args.is_plain_extracted_mode and args.msd_axis != 'total'))

    if args.title == '':
        title = args.charge_moment_calc_path
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub(r'(_charge_moment)?\.[^.]+$', '', args.charge_moment_calc_path) + \
                   '_charge_moment_msd.png'
    else:
        fig_path = args.fig_path

    with open(args.charge_moment_calc_path, 'r') as fp:
        if args.is_plain_extracted_mode:
            charge_moment = read_plain_extracted(fp)
        else:
            charge_moment = json.load(fp)
    plot_charge_moment(charge_moment,
                       args.msd_axis, args.msd_min, args.msd_max,
                       args.mean_axis, args.mean_min, args.mean_max,
                       args.energy_min, args.energy_max, args.to_plot_tb_energy_deviation,
                       args.time_start, args.time_end,
                       args.time_start_diffusion, args.time_end_diffusion,
                       args.is_plain_extracted_mode, title, fig_path)

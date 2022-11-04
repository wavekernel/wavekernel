# -*- coding: utf-8 -*-
import argparse, json, sys, re, pylab, os.path, matplotlib

kAngstromPerAu = 0.529177249  # Length
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0  # Length
kPsecPerAu = 2.418884326505e-5  # Time
kBoltzmann = 8.6173303e-5  # [eV / K]

def fig_path_to_params(s):  # Temporary information retrieve function.
    rf = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    m = re.search(r'out_(?P<k>(para|meta))_(?P<e>\d+)seed_(?P<t>%s)t_.*_(?P<s>%s)s_(?P<d>%s)d' % (rf, rf, rf), s)
    #m = re.search(r'out_(?P<k>(para|meta))_.*_(?P<s>%s)s_(?P<d>%s)d_(?P<g>%s)g' % (rf, rf, rf), s)
    if m:
        return m.group('k'), int(m.group('e')), float(m.group('t')), float(m.group('s')), float(m.group('d'))
        #return m.group('k'), float(m.group('s')), float(m.group('d')), float(m.group('g'))
    else:
        return None

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

def plot_diffusion_coef(ts, msds, label, time_start_diffusion, time_end_diffusion):
    if time_start_diffusion is None:
        time_start_diffusion = min(ts)
    if time_end_diffusion is None:
        time_end_diffusion = max(ts)
    xys = list(zip(ts, msds))
    xys = [xy for xy in xys if time_start_diffusion <= xy[0] and xy[0] < time_end_diffusion]
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
    ts_new = [xy[0] for xy in xys]
    ys_hat = [a * t + b for t in ts_new]
    intercept_relative_error_left = abs((b - xys[0][1]) / xys[0][1])
    intercept_relative_error_right = abs((a * xys[-1][0] + b - xys[-1][1]) / xys[-1][1])
    rmse = pylab.sqrt(sum([(xy[1] - a * xy[0] - b) ** 2.0 for xy in xys]) / float(n))
    print('T0, T1: ', xys[0][0], xys[-1][0])
    print('a T0 + b, y(T0), IRE(0): ', b, xys[0][1], intercept_relative_error_left)
    print('a T1 + b, y(T1), IRE(1): ', \
        a * xys[-1][0] + b, xys[0][1], intercept_relative_error_right)
    print('RMSE [angstrom^2]: ', rmse)
    print('diffusion coefficient [angstrom^2 / ps]: ', a / 2.0)
    print('diffusion coefficient [cm^2 / s]: ', a / 2.0 * 1e-4)
    # '1.0' = [e]
    # 'kBoltzmann * 1.0' = [V / K]
    mbt = a / 2.0 * 1e-4 / (kBoltzmann * 1.0)
    print('mobility * temperature [cm^2 K / V s]', mbt)
    if fig_path_to_params(fig_path) is not None:
        k, e, t, s, d = fig_path_to_params(fig_path)
        print('ZZZ', fig_path, k, e, t, s, d, mbt)
        #print 'ZZZ', fig_path, k, s, d, g, mbt
    pylab.plot(ts_new[0:n:n-1], ys_hat[0:n:n-1], 'x-', label=label, markersize=10)
    return mbt, rmse

def get_window_averaged_msds(ts, msds, window_width):
    new_ts = []
    new_msds = []
    # Remove dups first.
    t_prev = -1e100
    for t, msd in zip(ts, msds):
        if t != t_prev:
            new_msds.append(msd)
            new_ts.append(t)
            t_prev = t
    n = len(new_msds)
    m = n - window_width + 1
    msds_avg = [0.] * m
    for i in range(m):
        for j in range(window_width):
            msds_avg[i] += new_msds[i + j]
        msds_avg[i] /= window_width
    #print 'time range', max(new_ts[: m]) - min(new_ts[: m])
    return new_ts[: m], msds_avg

def plot_charge_moment(charge_moment,
                       msd_axis, msd_min, msd_max,
                       mean_axis, mean_min, mean_max,
                       energy_min, energy_max, to_plot_tb_energy_deviation,
                       time_start, time_end,
                       time_start_diffusion, time_end_diffusion, window_width,
                       is_raw_mode, title, fig_path):
    #font = {'size': 20}
    #matplotlib.rc('font', **font)
    fig = pylab.figure(figsize=(10, 7.5))
    # Cancel axis offset.
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)
    pylab.grid()

    axis_name_to_num = {'x': 0, 'y': 1, 'z': 2, 'total': 3}
    msd_axis_num = axis_name_to_num[msd_axis]
    mean_axis_num = axis_name_to_num[mean_axis]
    #is_fs_mode = max(charge_moment["ts"]) * kPsecPerAu < 0.001
    to_show_diffusion_coef = True #not is_fs_mode  # fs mode is not supported now.

    #if is_fs_mode:
    #    ts = map(lambda t: t * kPsecPerAu * 1000, charge_moment["ts"])
    #    pylab.xlabel("Time [fs]")
    #else:
    ts = [t * kPsecPerAu for t in charge_moment["ts"]]
    pylab.xlabel("Time [ps]")

    msds = [x[msd_axis_num] * kAngstrom2PerAu2 for x in charge_moment['msds']]
    means = [m[mean_axis_num] * kAngstromPerAu for m in charge_moment['means']]
    tb_energy_deviations = charge_moment['tb_energy_deviations']

    pylab.ylabel('MSD ' + msd_axis + ' [$\AA^2$]', color='blue')
    pylab.plot(ts, msds, '+-', markeredgecolor='blue', label=('MSD %s raw' % msd_axis),
               markerfacecolor='none')

    window_tss = []
    window_avg_msdss = []
    if window_widths is not None:
        for w in window_widths:
            window_ts, window_avg_msds = get_window_averaged_msds(ts, msds, w)
            window_tss.append(window_ts)
            window_avg_msdss.append(window_avg_msds)
            pylab.plot(window_ts, window_avg_msds, '+-', label=('MSD %s window %d' % (msd_axis, w)),
                       markerfacecolor='none')

    mbts = []
    rmses = []
    labels = []
    if to_show_diffusion_coef:  # and not is_fs_mode:
        mbt, rmse = plot_diffusion_coef(ts, msds, 'coef raw', time_start_diffusion, time_end_diffusion)
        mbts.append(mbt)
        rmses.append(rmse)
        labels.append('raw')
        if window_widths is not None:
            for w, window_ts, window_avg_msds in zip(window_widths, window_tss, window_avg_msdss):
                mbt, rmse = plot_diffusion_coef(window_ts, window_avg_msds, 'coef window %d' % w,
                                                time_start_diffusion, time_end_diffusion)
                mbts.append(mbt)
                rmses.append(rmse)
                labels.append('window %d' % w)

    #if msd_min is None:
    #    msd_min = pylab.ylim()[0]
    #if msd_max is None:
    #    msd_max = pylab.ylim()[1]
    #pylab.ylim(msd_min, msd_max)
    #yticks_new = list(pylab.yticks()[0])
    #yticks_new.extend([msd_min, msd_max])
    #pylab.yticks(yticks_new)
    pylab.ylim(msd_min, msd_max)  # limit setting again is needed (?).

    time_start_ = pylab.xlim()[0]
    time_end_ = pylab.xlim()[1]
    msd_min_ = pylab.ylim()[0]
    msd_max_ = pylab.ylim()[1]
    if to_show_diffusion_coef:
        for i, (mbt, rmse, label) in enumerate(zip(mbts, rmses, labels)):
            y_ratio = 0.94 - 0.04 * i
            pylab.text(time_start_ + (time_end_ - time_start_) * 0.01,
                       msd_min_ + (msd_max_ - msd_min_) * y_ratio,
                       '%s: %.2f [cm^2 K / V s], %.2f [$\AA^2$]' % (label, mbt, rmse))
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
   #
   #if to_plot_tb_energy_deviation:
   #    pylab.ylabel('TB energy deviation [a.u.]', color='red')
   #    pylab.plot(ts, tb_energy_deviations, '+-', color='red', label='TB energy dev')
   #    if energy_min is None:
   #        energy_min = pylab.ylim()[0]
   #    if energy_max is None:
   #        energy_max = pylab.ylim()[1]
   #    pylab.ylim(energy_min, energy_max)
   #    yticks_new = list(pylab.yticks()[0])
   #    yticks_new.extend([energy_min, energy_max])
   #    pylab.yticks(yticks_new)
   #    pylab.ylim(energy_min, energy_max)  # limit setting again is needed.
   #else:
   #    pylab.ylabel('Mean ' + mean_axis + ' [$\AA$]', color='red')
   #    pylab.plot(ts, means, '+-', color='red', label='Mean ' + mean_axis)
   #    if mean_min is None:
   #        mean_min = pylab.ylim()[0]
   #    if mean_max is None:
   #        mean_max = pylab.ylim()[1]
   #    pylab.ylim(mean_min, mean_max)
   #    yticks_new = list(pylab.yticks()[0])
   #    yticks_new.extend([mean_min, mean_max])
   #    pylab.yticks(yticks_new)
   #    pylab.ylim(mean_min, mean_max)  # limit setting again is needed.

    pylab.xlim(time_start, time_end)
    #xticks_new = list(pylab.xticks()[0])
    #xticks_new.extend([time_start, time_end])
    #pylab.xticks(xticks_new)
    #pylab.xlim(time_start, time_end)  # limit setting again is needed.

    pylab.legend()
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
    parser.add_argument('-w', metavar='WINDOW_WIDTHS', dest='window_widths_str', type=str, default=None,
                        help='')  # in comma separated list of window step size.
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

    if args.window_widths_str is None:
        window_widths = None
    else:
        window_widths = [int(s) for s in args.window_widths_str.split(',')]

    plot_charge_moment(charge_moment,
                       args.msd_axis, args.msd_min, args.msd_max,
                       args.mean_axis, args.mean_min, args.mean_max,
                       args.energy_min, args.energy_max, args.to_plot_tb_energy_deviation,
                       args.time_start, args.time_end,
                       args.time_start_diffusion, args.time_end_diffusion, window_widths,
                       args.is_plain_extracted_mode, title, fig_path)

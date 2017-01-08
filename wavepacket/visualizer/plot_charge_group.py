# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, re, os.path, matplotlib, math

kPsecPerAu = 2.418884326505e-5  # Time
kAuPerAngstrom = 1.8897259885789  # Length
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0  # Length

def plot_charge_group(charge_group, to_plot_msd_contribution, to_plot_diag,
                      time_start, time_end, group_start, group_end,
                      is_log_mode, log_min_exp, max_val_, stride, title, out_dir):
    num_groups = charge_group["num_groups"]
    if time_start is None:
        time_start = charge_group["charges_on_groups_all"][0]["t"]
    else:
        time_start /= kPsecPerAu
    if time_end is None:
        time_end = charge_group["charges_on_groups_all"][-1]["t"]
    else:
        time_end /= kPsecPerAu
    if group_start is None:
        group_start = 1
    if group_end is None:
        group_end = num_groups
    if max_val_ is None:
        if to_plot_msd_contribution:
            max_val = charge_group["msd_contributions_on_groups_max"] * kAngstrom2PerAu2
        else:
            max_val = charge_group["charges_on_groups_max"]
    else:
        max_val = max_val_
    if is_log_mode:
        min_val = 10.0 ** -log_min_exp
    else:
        min_val = 0

    font = {'size': 20}
    matplotlib.rc('font', **font)
    pylab.figure(figsize=(12, 7))
    i = j = 0
    for charge in charge_group["charges_on_groups_all"]:
        if time_start <= charge["t"] and charge["t"] <= time_end:
            if i % stride == 0:
                sys.stderr.write("plotting step number: " + str(charge["n"]) + "\n")
                pylab.title(title)
                pylab.xlabel("Group index", fontsize=20)
                if to_plot_msd_contribution:
                    pylab.ylabel("Contribution to MSD [$\AA^2$]", fontsize=20)
                else:
                    pylab.ylabel("Weight of charge", fontsize=20)
                pylab.xlim([group_start - 0.5, group_end + 0.5])
                pylab.ylim(min_val, max_val)
                xticks_new = list(pylab.xticks()[0])
                xticks_new.extend([group_start, group_end])
                pylab.xticks(xticks_new)
                pylab.xlim([group_start - 0.5, group_end + 0.5])  # limit setting again is needed.
                if not (max_val_ is None):
                    yticks_new = list(pylab.yticks()[0])
                    yticks_new.append(max_val)
                    pylab.yticks(yticks_new)
                    pylab.ylim(min_val, max_val)
                if to_plot_msd_contribution:
                    key = "co"
                    charge[key] = map(lambda x: x * kAngstrom2PerAu2, charge[key])  #
                else:
                    key = "ch"
                pylab.bar(range(group_start, group_end + 1),
                          charge[key][group_start - 1 : group_end], align="center", log=is_log_mode, color='black')
                if is_log_mode:
                    exp_min = math.log10(min_val)
                    exp_max = math.log10(max_val)
                    exp_diff = exp_max - exp_min
                    text_ys = map(lambda x: 10.0 ** (exp_max - exp_diff * x), [0.06, 0.14, 0.22])
                else:
                    diff = max_val - min_val
                    text_ys = map(lambda x: diff * x, [0.94, 0.86, 0.78])
                pylab.text(group_start, text_ys[0], "time [ps]: " + str(charge["t"] * kPsecPerAu))#, fontsize=24)
                pylab.text(group_start, text_ys[1], "MSD [$\AA^2$]: " + str(charge["m"] * kAngstrom2PerAu2))#, fontsize=24)
                pylab.text(group_start, text_ys[2], "IPRatio: " + str(charge["i"]))#, fontsize=24)
                if 'd' in charge and to_plot_diag:
                    pylab.twinx()
                    pylab.ylabel('diag(H)')
                    pylab.bar(range(group_start, group_end + 1),
                              charge['d'][group_start - 1 : group_end], align="center", log=is_log_mode, color='red', alpha=0.3)
                str_j = str(j)
                pylab.savefig(os.path.join(out_dir,
                                           "0" * (6 - len(str_j)) + str_j + ".png"))
                pylab.clf()
                j += 1
            i += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('charge_group_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('--stride', metavar='N', dest='stride', type=int, default=1,
                        help='')
    parser.add_argument('-c', action='store_true', dest='to_plot_msd_contribution',
                        default=False, help='')
    # start and end tiems: in picoseconds.
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')
    # Indices are inclusive (default: from 1 to num_groups).
    parser.add_argument('-l', metavar='LEFT', dest='group_start', type=int, default=None,
                        help='')
    parser.add_argument('-r', metavar='RIGHT', dest='group_end', type=int, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    parser.add_argument('--max', metavar='MAX', dest='max_val',
                        type=float, default=None, help='')
    parser.add_argument('--log', action='store_true', dest='is_log_mode',
                        default=False, help='')
    parser.add_argument('--log-min-exp', metavar='LOG_MIN_EXP', dest='log_min_exp',
                        type=float, default=8.0, help='')
    parser.add_argument('--no-diag', action='store_false', dest='to_plot_diag',
                        default=True, help='')

    args = parser.parse_args()

    out_dir = re.sub("\.[^.]+$", "", args.charge_group_path)
    if args.to_plot_msd_contribution:
        out_dir += "_msd_contrib"
    if args.is_log_mode:
        out_dir += "_log"
    if not os.path.isfile(args.charge_group_path):
        sys.stderr.write("[Error] file " + args.charge_group_path + " does not exist\n")
        sys.exit(1)
    if os.path.isdir(out_dir):
        sys.stderr.write("[Warning] output directory " + out_dir + " already exists\n")
    else:
        os.mkdir(out_dir)

    with open(args.charge_group_path, "r") as fp:
        charge_group = json.load(fp)

    if args.title == "":
        title = args.charge_group_path
    else:
        title = args.title

    plot_charge_group(charge_group, args.to_plot_msd_contribution, args.to_plot_diag,
                      args.time_start, args.time_end,
                      args.group_start, args.group_end,
                      args.is_log_mode, args.log_min_exp, args.max_val, args.stride, title, out_dir)

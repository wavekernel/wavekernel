# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re
kPsecPerAu = 2.418884326505e-5  # Time
kAuPerAngstrom = 1.8897259885789  # Length
kAngstromPerAu = kAuPerAngstrom ** -1.0  # Length
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0  # Length

def plot(alpha_calc, weight_max, time_start, time_end,
         num_printed_alpha_initial, num_printed_alpha_total,
         to_print_eigenstates, title, fig_path):
    alphas = alpha_calc["alphas"]
    fst_filter = alpha_calc["fst_filter"]
    num_filter = alpha_calc["num_filter"]

    if time_start is None:
        time_start = min(alpha_calc["ts"])
    if time_end is None:
        time_end = max(alpha_calc["ts"])
    # Convert limit times: picoseconds -> a.u.
    time_start = time_start / kPsecPerAu
    time_end = time_end / kPsecPerAu
    start_index = 0
    end_index = len(alpha_calc["ts"])
    for i in range(len(alpha_calc["ts"])):
        if alpha_calc["ts"][i] >= time_start:
            start_index = i
            break
    for i in range(len(alpha_calc["ts"])):
        if alpha_calc["ts"][i] >= time_end:
            end_index = i
            break
    ts = alpha_calc["ts"][start_index : end_index]

    pylab.subplot(1,2,1)
    if weight_max is None:
        ylim_max = max(map(lambda alpha:max(alpha["weights"][start_index : end_index]), alphas)) * 1.05
        pylab.ylim(0.0, ylim_max)
    else:
        ylim_max = weight_max
        pylab.ylim(0.0, ylim_max)
        yticks_new = list(pylab.yticks()[0])
        yticks_new.append(ylim_max)
        pylab.yticks(yticks_new)
        pylab.ylim(0.0, ylim_max)  # limit setting again is needed.

    # Convert data times: a.u. -> picoseconds.
    if max(ts) * kPsecPerAu < 0.001:
        ts = map(lambda t: t * kPsecPerAu * 1000, ts)
        pylab.xlabel("Time [fs]")
    else:
        ts = map(lambda t: t * kPsecPerAu, ts)
        pylab.xlabel("Time [ps]")
    pylab.ylabel("Weight")
    pylab.title(title)

    plotted_indices = set()
    for alpha in alphas:
        alpha["total_weight"] = sum(alpha["weights"][start_index : end_index])
        alpha["initial_weight"] = alpha["weights"][0]
    alphas.sort(key=lambda alpha: -alpha["total_weight"])
    plotted_indices.update(map(lambda alpha: alpha["i"], alphas[: num_printed_alpha_total]))
    alphas.sort(key=lambda alpha: -alpha["initial_weight"])
    plotted_indices.update(map(lambda alpha: alpha["i"], alphas[: num_printed_alpha_initial]))
    alphas_max_is_over_threshold = filter(lambda alpha: max(alpha["weights"]) >= 0.1, alphas)
    plotted_indices.update(map(lambda alpha: alpha["i"], alphas_max_is_over_threshold))

    #max_eigenvalue = max(map(lambda alpha: alpha["eigenvalue"], alphas))
    pylab.plot([], [], label="index(energy [a.u.], Mean X [$\AA$], MSD [$\AA^2$])", color="white")
    if to_print_eigenstates:
        print "# index, energy [a.u.], Mean X [angstrom], MSD [angstrom^2]"
    for alpha in alphas:
        if alpha["i"] in plotted_indices:
            plotted = (alpha["i"],
                       alpha["eigenvalue"],
                       alpha["mean"][0] * kAngstromPerAu,  # Mean of X axis.
                       alpha["msd"][3] * kAngstrom2PerAu2)  # Total MSD.
            pylab.plot(ts,
                       alpha["weights"][start_index : end_index],
                       "-",
                       label="%d (%.10f, %.1f, %.1f)" % plotted)
            if to_print_eigenstates:
                print "%d %f %f %f" % plotted
    pylab.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., mode="expand")
    pylab.savefig(fig_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('alpha_calc_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('--weight-max', metavar='WEIGHT_MAX', dest='weight_max', type=float, default=None,
                        help='')
    # start and end tiems: in picoseconds.
    parser.add_argument('-s', metavar='TIME_START', dest='time_start', type=float, default=None,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='')
    parser.add_argument('--ninit', metavar='NUM_PRINTED_ALPHA_INITIAL', dest='num_printed_alpha_initial',
                        type=int, default=10, help='')
    parser.add_argument('--ntotal', metavar='NUM_PRINTED_ALPHA_TOTAL', dest='num_printed_alpha_total',
                        type=int, default=10, help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    parser.add_argument('-p', action='store_true', dest='to_print_eigenstates',
                        default=False, help='')

    args = parser.parse_args()

    fig_path = re.sub("(_alpha)?\.[^.]+$", "", args.alpha_calc_path) + "_alpha.png"
    if args.title == "":
        title = args.alpha_calc_path
    else:
        title = args.title

    with open(args.alpha_calc_path, "r") as fp:
        alpha_calc = json.load(fp)
    plot(alpha_calc, args.weight_max, args.time_start, args.time_end,
         args.num_printed_alpha_initial, args.num_printed_alpha_total,
         args.to_print_eigenstates, title, fig_path)

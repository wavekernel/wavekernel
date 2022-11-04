# -*- coding: utf-8 -*-

import json, sys, re, pylab, os.path, matplotlib, argparse

kAngstromPerAu = 0.529177249
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0
kPsecPerAu = 2.418884326505e-5

def plot_charge_multi(charge_moments_with_labels, fig_path):
    #font = {'size': 20}
    #matplotlib.rc('font', **font)
    #pylab.figure(figsize=(12, 12))
    pylab.subplot(3, 1, 1)
    pylab.ylabel("MSD x [$\AA^2$]")
    pylab.subplot(3, 1, 2)
    pylab.ylabel("Mean x [$\AA$]")
    pylab.subplot(3, 1, 3)
    pylab.ylabel("ipratio")
    for entry in charge_moments_with_labels:
        charge_moment, label = entry[0], entry[1]
        times = [t * kPsecPerAu for t in charge_moment["ts"]]
        msds_x = [msd[0] for msd in charge_moment["msds"]]
        msds_x_ang = [x * kAngstrom2PerAu2 for x in msds_x]
        ipratios = charge_moment["ipratios"]
        means_x = [mean[0] for mean in charge_moment["means"]]
        means_x_ang = [x * kAngstromPerAu for x in means_x]
        label_folded = ""
        max_len_in_line = 8
        for i_begin in range(0, len(label), max_len_in_line):
            i_end = min(i_begin + max_len_in_line, len(label))
            label_folded += label[i_begin : i_end]
            if i_end != len(label):
                label_folded += "\n"
        pylab.subplot(3, 1, 1)
        pylab.plot(times, msds_x_ang, label=label_folded)
        pylab.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        pylab.subplot(3, 1, 2)
        pylab.plot(times, means_x_ang, label=label_folded)
        pylab.subplot(3, 1, 3)
        pylab.plot(times, ipratios, label=label_folded)
    #pylab.legend(loc='lower right')
    pylab.xlabel("time [ps]")#, fontsize=20)
    pylab.subplots_adjust(right=0.75)
    pylab.savefig(fig_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('charge_moment_paths', metavar='JSON', type=str, nargs='+',
                        help='plural *_charge_moment.json')
    parser.add_argument('-l', dest='legend_labels', default='',
                        help='set labels for the graph delimited with comma (default: input JSON paths)')
    parser.add_argument('-o', dest='fig_path', default='charge_multi.png',
                        help='set output figure path')
    args = parser.parse_args()

    charge_moment_paths = args.charge_moment_paths
    charge_moments = []
    for path in charge_moment_paths:
        with open(path, "r") as fp:
            charge_moments.append(json.load(fp))

    if args.legend_labels == '':
        legend_labels = charge_moment_paths
    else:
        legend_labels = args.legend_labels.split(",")
        assert(len(charge_moments) == len(legend_labels))

    plot_charge_multi(list(zip(charge_moments, legend_labels)), args.fig_path)

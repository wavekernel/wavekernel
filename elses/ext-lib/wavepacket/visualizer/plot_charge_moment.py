# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re

kPsecPerAu = 2.418884326505e-5  # Time
kAngstromPerAu = 0.529177249  # Length

def plot(charge_moment, ymin, ymax, title, fig_path):
    means = map(list, zip(*charge_moment["means"]))
    msds = map(list, zip(*charge_moment["msds"]))
    means_ang = map(lambda xs: map(lambda x: x * kAngstromPerAu, xs), means)
    variances_ang = map(lambda xs: map(lambda x: math.sqrt(x) * kAngstromPerAu, xs), msds)

    if max(charge_moment["ts"]) * kPsecPerAu < 0.001:
        ts = map(lambda t: t * kPsecPerAu * 1000, charge_moment["ts"])
        pylab.xlabel("Time [fs]")
    else:
        ts = map(lambda t: t * kPsecPerAu, charge_moment["ts"])
        pylab.xlabel("Time [ps]")

    minval = min(map(lambda xs: min(xs), means_ang + variances_ang))
    maxval = max(map(lambda xs: max(xs), means_ang + variances_ang))
    if ymin is None:
        ymin = minval - (maxval - minval) * 0.1
    if ymax is None:
        ymax = maxval + (maxval - minval) * 0.5
    pylab.ylim(ymin, ymax)
    yticks_new = list(pylab.yticks()[0])
    yticks_new.extend([ymin, ymax])
    pylab.yticks(yticks_new)
    pylab.ylim(ymin, ymax)  # limit setting again is needed.
    pylab.ylabel("Mean / Variance [$\AA$]", fontsize=20)

    pylab.title(title)
    pylab.plot(ts, means_ang[0], color="blue", label="X mean")
    pylab.plot(ts, means_ang[1], color="green", label="Y mean")
    pylab.plot(ts, means_ang[2], color="red", label="Z mean")
    pylab.plot(ts, variances_ang[0], color="blue",
               linestyle="--", label="X variance")
    pylab.plot(ts, variances_ang[1], color="green",
               linestyle="--", label="Y variance")
    pylab.plot(ts, variances_ang[2], color="red",
               linestyle="--", label="Z variance")
    pylab.plot(ts, variances_ang[3], color="black",
               linestyle="--", label="total variance")
    pylab.legend(ncol=2)
    pylab.savefig(fig_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('charge_moment', metavar='JSON', type=str,
                        help='')
    parser.add_argument('--ymin', metavar='YMIN', dest='ymin', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='YMAX', dest='ymax', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    args = parser.parse_args()

    with open(args.charge_moment, "r") as fp:
        charge_moment = json.load(fp)
    fig_path = re.sub("(_charge_moment)?\.[^.]+$", "", args.charge_moment) + "_charge_moment.png"
    if args.title == "":
        title = args.charge_moment
    else:
        title = args.title
    plot(charge_moment, args.ymin, args.ymax, title, fig_path)

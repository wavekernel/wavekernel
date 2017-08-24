# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, re, os.path, matplotlib, math

kPsecPerAu = 2.418884326505e-5  # Time
kAuPerAngstrom = 1.8897259885789  # Length
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0  # Length

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('charge_group_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default="",
                        help='')
    parser.add_argument('--log', action='store_true', dest='is_log_mode',
                        default=False, help='')
    args = parser.parse_args()

    with open(args.charge_group_path, "r") as fp:
        charge_group = json.load(fp)
    
    if args.title == "":
        title = args.charge_group_path
    else:
        title = args.title

    charges = []
    time = map(lambda f:f["t"]*kPsecPerAu,charge_group["charges_on_groups_all"])
    for i in xrange(charge_group["num_groups"]) :
        charge = map(lambda f:f["ch"][i],charge_group["charges_on_groups_all"])
        pylab.plot(time,charge,label="group%d"%(i+1),linewidth=3)
        charges.append(charge)
    if args.is_log_mode :
        pylab.yscale('log')
        output = os.path.join(args.charge_group_path[:-5]+"_log.png")
    else :
        output = os.path.join(args.charge_group_path[:-5]+".png")

    pylab.title(title)
    #pylab.ylim(0.0,0.017)
    pylab.grid()
    pylab.legend()
    pylab.xlabel("Time [ps]")
    pylab.ylabel("Weight")
    pylab.savefig(output)

    with open(args.charge_group_path[:-5]+".dat","w") as fp :
        fp.write("# %18s  "%("time[ps]"))
        for i in xrange(len(charges)) :
            string = "group-%d "%(i+1)
            fp.write("%20s "%string)
        fp.write("\n")
        for data in zip(time,*charges) :
            for value in data :
                fp.write("%20s "%(str(value)))
            fp.write("\n")

# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, pylab
kAuPerAngstrom = 1.8897259885789
kAngstromPerAu = kAuPerAngstrom ** -1.0
kAngstrom2PerAu2 = kAuPerAngstrom ** -2.0

def select(wavekernel_out, num_states, min_energy):
    condition = wavekernel_out['condition']
    fst_filter = wavekernel_out['setting']['fst_filter']
    end_filter = wavekernel_out['setting']['end_filter']
    indices = range(fst_filter, end_filter + 1)
    eigenvalues = condition['eigenvalues']
    eigenstate_msds = condition['eigenstate_msd_total']
    eigenstate_msds = map(lambda x: x * kAngstrom2PerAu2, eigenstate_msds)
    eigenstate_means = condition['eigenstate_mean_x']
    eigenstate_means = map(lambda x: x * kAngstromPerAu, eigenstate_means)

    dim = condition['dim']
    homo = dim / 2
    left = min(eigenstate_means)
    right = max(eigenstate_means)
    center = left + (right - left) / 2.0
    assert(homo * 2 == dim and fst_filter <= homo and homo <= end_filter)

    print '# Index, Energy, Mean X, MSD, dist from center'
    states = []
    for i in range(end_filter - fst_filter + 1):
        energy = eigenvalues[i]
        if energy >= min_energy:
            mean_x = eigenstate_means[i]
            msd = eigenstate_msds[i]
            states.append((i + fst_filter, energy, mean_x, msd))
    states.sort(key=lambda s: abs(s[2] - center))
    for i in range(min(num_states, len(states))):
        (j, energy, mean_x, msd) = states[i]
        print j, energy, mean_x, msd, abs(mean_x - center)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-n', metavar='NUM', dest='num_states', type=int, default=30,
                        help='')
    parser.add_argument('-e', metavar='MIN', dest='min_energy', type=float, default=-1e10,
                        help='')
    args = parser.parse_args()

    with open(args.wavekernel_out_path, 'r') as fp:
        wavekernel_out = json.load(fp)
    select(wavekernel_out, args.num_states, args.min_energy)

# Old codes.
#def main():
#    filename = sys.argv[1]
#    all_list = eigenlevel_mean_msd_list(filename)
#    # input parameter------------------------------
#    benzene = 100      # meta and para for 1200atom
#    mean_center = 350  # angstrom for para 1200atom
#    #mean_center = 300 # angstrom for meta 1200atom
#    mean_side   = 10   # angstrom
#    msd_min = 1e+2
#    msd_max = 1e+4
#    mean_side_limit = 300
#    #-----------------------------------------------
#    for name in all_list:
#        all_list[name] = all_list[name][len(all_list[name])-benzene:]
#
#    while True:
#        select_eigenlevel = []
#        #print "mean_side=",mean_side
#        for i in range(len(all_list["mean_x"])):
#            judge_num = abs(all_list["mean_x"][i] - mean_center)
#            if judge_num < mean_side:
#                if all_list["msd"][i] <= msd_max and msd_min <= all_list["msd"][i]:
#                    select_eigenlevel.insert(0,all_list["index"][i])
#                    #print all_list["index"][i]
#                else:
#                    pass
#            else:
#                pass
#        if len(select_eigenlevel) >= 10 or mean_side >= mean_side_limit:
#            break
#        else:
#            mean_side += 1
#
#    print 'select eigenlevel-> "',
#    for i in range(len(select_eigenlevel)):
#        print select_eigenlevel[i],
#    print '"'
#    print "mean_side =",mean_side
#
#def eigenlevel_mean_msd_list(filename):
#    all_list={"index" :[],"energy":[],
#              "mean_x":[],"msd"   :[]}
#    for line in open(filename,"r"):
#        if line.find("#")>=0:
#            pass
#        else:
#            line = str2list(line)
#            all_list["index"].insert(0,int(line[0]))
#            all_list["energy"].insert(0,float(line[1]))
#            all_list["mean_x"].insert(0,float(line[2]))
#            all_list["msd"].insert(0,float(line[3]))
#    return all_list
#
#def str2list(string):
#    string = string.strip()
#    string = re.sub("  *"," ",string)
#    return re.split(" ",string)
#
#if __name__=="__main__":
#    main()

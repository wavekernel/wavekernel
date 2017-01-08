# -*- coding: utf-8 -*-
import argparse, json, sys, re, pylab, os.path, matplotlib
kPsecPerAu = 2.418884326505e-5  # Time

def read_raw(fp):
    rf = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    ts = []
    aps = []
    for line in fp:
        m = re.search(r'(?P<t>%s)\s+(?P<i>\d+)\s+(?P<isamr>True|False)\s+(?P<e1>%s)\s(?P<e2>%s)\s(?P<e3>%s)\s(?P<mean>%s)\s+(?P<msd>%s)\s+(?P<pip>%s)\s+(?P<aip>%s)' %
                      (rf, rf, rf, rf, rf, rf, rf, rf), line)
        if m:
            t = float(m.group('t')) * kPsecPerAu
            aip = float(m.group('aip'))
            ts.append(t)
            aps.append(1.0 / aip)
    return {'ts': ts, 'aps': aps}

if __name__ == '__main__':
    with open(sys.argv[1]) as fp:
        data = read_raw(fp)

    fig_path = re.sub(r'\.[^.]+$', '', sys.argv[1]) + '_alpha_pratio.png'        
    pylab.plot(data['ts'], data['aps'], '+')
    pylab.xlabel('time [ps]')
    pylab.ylabel('PR')
    #pylab.ylim(1.0, 600.0)
    pylab.title(sys.argv[1])
    pylab.savefig(fig_path)

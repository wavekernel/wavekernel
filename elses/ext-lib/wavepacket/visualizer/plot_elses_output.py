# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, copy, pylab, math

kAuPerAngstrom = 1.8897259885789
kPsecPerAu = 2.418884326505e-5  # Time
kJoulePerEh = 4.35974434e-18  # [J / a.u. (energy)]
kBoltzmannInSI = 1.3806488e-23  # [J / K]
kAuPerKelvin = kBoltzmannInSI / kJoulePerEh

def read_elses_output(to_read_from_stdout, fp):
    steps = []
    potentials = []
    for line in fp:
        if to_read_from_stdout:
            m = re.search(r'IONMOV3\s+\d+\s+(\d+)\s+([-\d\.]+)', line)
            if m:
                print line.strip(), '->', m.group(1), m.group(2)
                step = int(m.group(1))
                potential = float(m.group(2))            
                steps.append(step)
                potentials.append(potential)            
        else:
            m = re.search(r'Energy summary \(eV/atom\):\s*(\d+)\s+([-\d\.]+)', line)
            if m:
                print line.strip(), '->', m.group(1), m.group(2)
                step = int(m.group(1))
                potential = float(m.group(2))            
                steps.append(step)
                potentials.append(potential)            
    return (steps, potentials)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('elses_outputs', metavar='ELSES_OUTPUTS', type=str, nargs='+',
                        help='')
    parser.add_argument('--ymin', metavar='YMIN', dest='y_min', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='YMAX', dest='y_max', type=float, default=None,
                        help='')
    parser.add_argument('-t', metavar='TITLE', dest='title', type=str, default=None,
                        help='')
    parser.add_argument('-s', metavar='STEP_START', dest='step_start', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-e', metavar='STEP_END', dest='step_end', type=float, default=None,
                        help='')  # in ps.
    parser.add_argument('-o', metavar='OUT', dest='fig_path', type=str, default=None,
                        help='')
    parser.add_argument('-l', metavar='LABELS', dest='labels', type=str, default=None,
                        help='')
    parser.add_argument('--stdout', action='store_true', dest='to_read_from_stdout',
                        default=False, help='')    
    args = parser.parse_args()

    if args.title is None:
        title = args.elses_outputs[0]
    else:
        title = args.title

    if args.fig_path is None:
        fig_path = re.sub(r'\.[^.]+$', '', args.elses_outputs[0]) + '.png'
    else:
        fig_path = args.fig_path

    if args.labels is None:
        labels = map(lambda i: 'sample %d' % i, range(len(args.elses_outputs)))
    else:
        labels = args.labels.split(',')
    assert(len(labels) == len(args.elses_outputs))

    outputs = []
    for path in args.elses_outputs:
        with open(path, 'r') as fp:
            outputs.append(read_elses_output(args.to_read_from_stdout, fp))

    header = re.sub('\.[^.]+$', '', fig_path)

    if args.step_start is None:
        step_start = min(outputs[0][0])
    else:
        step_start = args.step_start
    if args.step_end is None:
        step_end = max(outputs[0][0])
    else:
        step_end = args.step_end

    if args.y_min is None:
        y_min = min(outputs[0][1])
    else:
        y_min = args.y_min
    if args.y_max is None:
        y_max = max(outputs[0][1])
    else:
        y_max = args.y_max    

    output_filename = '%s.png' % header
    fig = pylab.figure()
    ax = fig.gca()
    ax.ticklabel_format(useOffset=False)    
    pylab.title(title)
    pylab.xlabel('Step')
    if args.to_read_from_stdout:
        pylab.ylabel('Potential [au]')
    else:
        pylab.ylabel('Potential [eV/atom]')
    pylab.grid()
    pylab.xlim(step_start, step_end)
    pylab.ylim(y_min, y_max)
    for i, output in enumerate(outputs):
        pylab.plot(output[0], output[1], '+-', label=labels[i])
    pylab.legend(loc='upper right')
    pylab.savefig(output_filename)

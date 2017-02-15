# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, numpy
kAuPerAngstrom = 1.8897259885789
kKelvinPerAu = 3.157746e5
kSizeOfReal = 8

def get_real_array(split_dir, is_little_endian, element):
    if element == []:
        return []
    elif isinstance(element[0], basestring):  # Supposed to be binary output mode.
        first = element[1]
        last = element[2]
        count = last - first + 1
        with open(os.path.join(split_dir, element[0]), 'rb') as fp:
            fp.seek(kSizeOfReal * (first - 1), os.SEEK_SET)
            xs_str = fp.read(kSizeOfReal * count)
        format_char_endian = '<' if is_little_endian else '>'
        return struct.unpack(format_char_endian + str(count) + 'd', xs_str)
    else:
        return element

def read_xyz(fp):
    regexp_float = r'[+-]?(\d+\.)?\d+([deE][+-]?\d+)?'
    xyz = []
    is_first_line = True
    for line in fp:
        if is_first_line:
            m = re.search(r'(\d+)', line)
            if m:
                num_atoms = int(m.group(1))
                is_first_line = False
            else:
                sys.exit(1)
        else:
            m = re.match(r'(?P<a>\w+)\s+(?P<x>%s)\s+(?P<y>%s)\s+(?P<z>%s)' %
                         (regexp_float, regexp_float, regexp_float), line)
            if m:
                atom = m.group('a')
                pos = [float(m.group('x')) * kAuPerAngstrom,
                       float(m.group('y')) * kAuPerAngstrom,
                       float(m.group('z')) * kAuPerAngstrom]
                xyz.append((atom, numpy.array(pos)))
        if len(xyz) >= num_atoms:
            break
    return xyz

def get_moments(xs):
    mean = sum(xs) / len(xs)
    variance = sum(map(lambda x: (x - mean) ** 2.0, xs)) / len(xs)
    return (mean, variance)

def add_step_maxwell(state, split_dir, is_little_endian, xyz):
    element_to_mass = {'H': 1.8371526e3, 'C': 2.1874661e4}
    atom_mass = map(lambda atom: element_to_mass[atom[0]], xyz)
    atom_speed = get_real_array(split_dir, is_little_endian, state['atom_speed'])
    atom_energy = map(lambda (mass, speed): 0.5 * mass * (speed ** 2.0), zip(atom_mass, atom_speed))
    (mean, variance) = get_moments(atom_energy)
    print mean * kKelvinPerAu, numpy.sqrt(variance) * kKelvinPerAu, \
        max(map(lambda e: abs(e - mean), atom_energy)) * kKelvinPerAu

def add_step_harmonic(state, split_dir, is_little_endian, xyz):
    atom_energy = get_real_array(split_dir, is_little_endian, state['atom_perturb'])
    (mean, variance) = get_moments(atom_energy)
    print mean * kKelvinPerAu, numpy.sqrt(variance) * kKelvinPerAu, \
        max(map(lambda e: abs(e - mean), atom_energy)) * kKelvinPerAu

def calc(wavepacket_out, wavepacket_out_path, xyz, stride, is_little_endian, start_time):
    if wavepacket_out['setting']['h1_type'] == 'maxwell':
        is_maxwell = True
    elif wavepacket_out['setting']['h1_type'] == 'harmonic':
        is_maxwell = False
    else:
        assert(False)
    print 'mean [K], deviation [K], largest diff from mean [K]'
    if wavepacket_out['setting']['is_output_split']:
        split_dir = os.path.dirname(wavepacket_out_path)
        for meta in wavepacket_out['split_files_metadata']:
            path = os.path.join(split_dir, meta['filename'])
            with open(path, 'r') as fp:
                diff = datetime.datetime.now() - start_time
                sys.stderr.write(str(diff) + " reading: " + path + "\n")
                states_split = json.load(fp)
            for state in states_split:
                if state['step_num'] % stride == 0:
                    if is_maxwell:
                        add_step_maxwell(state, split_dir, is_little_endian, xyz)
                    else:
                        add_step_harmonic(state, split_dir, is_little_endian, xyz)
    else:
        split_dir = ''
        for state in wavepacket_out['states']:
            if state['step_num'] % stride == 0:
                if is_maxwell:
                    add_step_maxwell(state, split_dir, is_little_endian, xyz)
                else:
                    add_step_harmonic(state, split_dir, is_little_endian, xyz)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavepacket_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('xyz_path', metavar='XYZ', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    with open(args.wavepacket_out_path, 'r') as fp:
        wavepacket_out = json.load(fp)
    with open(args.xyz_path, 'r') as fp:
        xyz = read_xyz(fp)
    calc(wavepacket_out, args.wavepacket_out_path, xyz, args.skip_stride_num, args.is_little_endian, start_time)

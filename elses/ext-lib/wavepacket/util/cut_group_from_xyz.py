# -*- coding: utf-8 -*-
import sys, re

def read_xyz(fp):
    xyz = []
    for i, line in enumerate(fp):
        line = line.strip()
        if i > 1:
            xyz.append(line)
    return xyz

def read_group_id(fp):
    is_nums_read = False
    for line in fp:
        if is_nums_read:
            m = re.search('(\d+)\s+(\d+)', line)
            if m:
                atom_to_group[int(m.group(1)) - 1] = int(m.group(2)) - 1
        else:
            m = re.search('(\d+)\s+(\d+)\s+(\d+)', line)
            if m:
                is_nums_read = True
                num_atoms = int(m.group(1))
                num_groups = int(m.group(2))
                num_max_group_size = int(m.group(3))
                atom_to_group = [0] * num_atoms
    return {'num_groups': num_groups,
            'atom_to_group': atom_to_group}

def print_xyz(xyz, fp):
    num_atoms = len(xyz)
    fp.write('%d\n' % num_atoms)
    fp.write('#\n')
    for line in xyz:
        fp.write(line + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'Usage: python cut_group_from_xyz.py <xyz> <group id> <group start> <group end>'
        sys.exit(0)
    xyz_path = sys.argv[1]
    group_id_path = sys.argv[2]
    # group id file format is 1-origin.
    cut_group_start = int(sys.argv[3]) - 1
    cut_group_final = int(sys.argv[4]) - 1
    with open(xyz_path, 'r') as fp:
        xyz = read_xyz(fp)
    with open(group_id_path, 'r') as fp:
        group_id = read_group_id(fp)
    xyz_cut = []
    for i, line in enumerate(xyz):
        group = group_id['atom_to_group'][i]
        if cut_group_start <= group and group <= cut_group_final:
            xyz_cut.append(line)
    out_filename = re.sub('\.[^.]+$', '', xyz_path) + \
                   ('_group%d-%d.xyz' % (cut_group_start + 1, cut_group_final + 1))
    with open(out_filename, 'w') as fp:
        print_xyz(xyz_cut, fp)

# -*- coding: utf-8 -*-
# π共役鎖系のベンゼン環にランダムなひねりを与える.
import argparse, sys, re, numpy, copy
kAuPerAngstrom = 1.8897259885789
kRadianPerDegree = numpy.pi / 180.0

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

def make_group_to_atoms(xyz, group_id):
    group_to_atoms = []
    for i in range(group_id['num_groups']):
        group_to_atoms.append([])
    for i in range(len(xyz)):
        g = group_id['atom_to_group'][i]
        group_to_atoms[g].append(xyz[i])
    return group_to_atoms

#               C3 - C5
#              /       \
# C0 ≡ C1 - C2         C7
#              \       /
#               C4 - C6

# [numpy.array] -> numpy.array
# Get vector (C0 => C1).
def get_triple_bond_vector(carbons):
    c0 = get_rotate_origin(carbons)
    remains = filter(lambda c: (c != c0).any(), carbons)
    # C1 is the nearest carbon to C0.
    c1 = min(remains, key=lambda c: numpy.linalg.norm(c - c0))
    return c1 - c0

# [numpy.array] -> numpy.array
# Get point C0.
def get_rotate_origin(carbons):
    n = len(carbons)
    assert(n == 8)
    center = numpy.sum(carbons, axis=0) / n
    # C0 is the most distant carbon from the center of gravity.
    c0 = max(carbons, key=lambda c: numpy.linalg.norm(c - center))
    return c0

# (numpy.array, [atom]) -> [atom]
# atom = (element, numpy.array)
def shift_atoms(origin, atoms):
    return map(lambda (e, x): (e, x - origin), atoms)

# (numpy.array, float, numpy.array) -> numpy.array
def rotate_around(rotate_axis, angle, position):
    rotate_axis /= numpy.linalg.norm(rotate_axis)
    nx, ny, nz = rotate_axis[0], rotate_axis[1], rotate_axis[2]
    c = numpy.cos(angle)
    cm = 1.0 - c
    s = numpy.sin(angle)
    rotator = numpy.array(
        [[nx * nx * cm + c, nx * ny * cm - nz * s, nz * nx * cm + ny * s],
         [nx * ny * cm + nz * s, ny * ny * cm + c, ny * nz * cm - nx * s],
         [nz * nx * cm - ny * s, ny * nz * cm + nx * s, nz * nz * cm + c]])
    return numpy.dot(rotator, position)

def sample_angle(deviation):
    if deviation == 0.0:
        return 0.0
    else:
        return numpy.random.normal(0.0, deviation)

def is_rotatable_group(atoms_in_group):
    carbons = filter(lambda a: a[0] == 'C', atoms_in_group)
    hydrogens = filter(lambda a: a[0] == 'H', atoms_in_group)
    return (len(carbons) == 8 and (len(hydrogens) == 4 or len(hydrogens) == 5))

def rotate_group(atoms_in_group, angle):
    carbons = filter(lambda a: a[0] == 'C', atoms_in_group)
    carbons_coord = map(lambda a: a[1], carbons)
    rotate_axis = get_triple_bond_vector(carbons_coord)
    origin = get_rotate_origin(carbons_coord)
    shifted = shift_atoms(origin, atoms_in_group)
    rotated = map(lambda a: (a[0], rotate_around(rotate_axis, angle, a[1])), shifted)
    return shift_atoms(-origin, rotated)

def print_xyz(xyz, fp):
    num_atoms = len(xyz)
    fp.write('%d\n' % num_atoms)
    fp.write('#\n')
    for atom in xyz:
        s = '%s %.6f %.6f %.6f\n' % (atom[0],
                                     atom[1][0] / kAuPerAngstrom,
                                     atom[1][1] / kAuPerAngstrom,
                                     atom[1][2] / kAuPerAngstrom)
        fp.write(s)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('xyz_path', metavar='XYZ', type=str,
                        help='')
    parser.add_argument('group_id_path', metavar='GROUP_ID', type=str,
                        help='')
    parser.add_argument('angle_deviation', metavar='ANGLE', type=float,
                        help='')
    parser.add_argument('-l', metavar='LOG', dest='log_filename', type=str,
                        default='rotate_piconjpolymer_log.txt', help='')
    args = parser.parse_args()

    angle_deviation = args.angle_deviation * kRadianPerDegree  # degree -> radian.
    with open(args.xyz_path, 'r') as fp:
        xyz = read_xyz(fp)
    with open(args.group_id_path, 'r') as fp:
        group_id = read_group_id(fp)
    group_to_atoms = make_group_to_atoms(xyz, group_id)

    out_filename = re.sub('\.[^.]+$', '_tav%d.xyz' %
                          int(angle_deviation / kRadianPerDegree), args.xyz_path)

    rotated_atoms = []
    with open(args.log_filename, 'w') as fp:
        fp.write('filename: %s\n' % out_filename)
        fp.write('input degree deviation: %f\n' % (angle_deviation / kRadianPerDegree))
        for group in range(group_id['num_groups']):
            atoms_in_group = group_to_atoms[group]
            if is_rotatable_group(atoms_in_group):
                angle = sample_angle(angle_deviation)
                fp.write('%d %f\n' % (group + 1, angle / kRadianPerDegree))
                rotated_group = rotate_group(atoms_in_group, angle)
            else:
                fp.write('%d skip\n' % (group + 1))
                rotated_group = atoms_in_group
            rotated_atoms.extend(rotated_group)

    with open(out_filename, 'w') as fp:
        print_xyz(rotated_atoms, fp)

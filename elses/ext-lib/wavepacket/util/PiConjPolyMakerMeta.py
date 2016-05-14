# coding: utf-8
import numpy as np
import math
import sys

def orthogonal_to_spherical(vector):
    r = np.linalg.norm(vector)
    if r == 0.0:
        return np.array([0.0, 0.0, 0.0])
    theta = np.arccos(vector[2] / r)
    phi = np.arctan2(vector[1], vector[0])
    return np.array([r, theta, phi])

def spherical_to_orthogonal(vector):
    r, theta, phi = vector[0], vector[1], vector[2]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])

def rotate_around(vector, axis, theta):
    axis /= np.linalg.norm(axis)
    nx, ny, nz = axis[0], axis[1], axis[2]
    c = np.cos(theta)
    s = np.sin(theta)
    c1 = 1.0 - c
    rotator = np.array([[nx * nx * c1 + c, nx * ny * c1 - nz * s, nz * nx * c1 + ny * s],
                        [nx * ny * c1 + nz * s, ny * ny * c1 + c, ny * nz * c1 - nx * s],
                        [nz * nx * c1 - ny * s, ny * nz * c1 + nx * s, nz * nz * c1 + c]])
    return np.dot(rotator, vector)

def add_terminal_hydrogen(joint, bond_length):
    start = joint['end_point']
    di = joint['direction']
    di /= np.linalg.norm(di)
    structure = [(start + di * bond_length, 'H')]
    next_joint = {'end_point': structure[0][0], 'direction': di}
    return (structure, next_joint)

def add_carbons_triple_bond(joint, bond_length):
    carbons_length = 1.18559
    start = joint['end_point']
    di = joint['direction']
    di /= np.linalg.norm(di)
    positions = [start + di * bond_length]
    positions.append(positions[0] + di * carbons_length)
    next_joint = {'end_point': positions[1], 'direction': di}
    elements = ['C', 'C']
    return (zip(positions, elements), next_joint)

def add_benzene(joint, bond_length, plane_normal, next_carbon_num):
    assert(1 <= next_carbon_num and next_carbon_num <= 5)
    a = 1.42514  # C-C length in benzene.
    b = 1.06300  # C-H length in benzene.
    start = joint['end_point']
    di = joint['direction']
    di /= np.linalg.norm(di)
    di1 = rotate_around(di, plane_normal, math.pi / 3.0)
    di2 = rotate_around(di, plane_normal, - math.pi / 3.0)

    position_carbon_0 = start + di * bond_length
    positions_template = [
        di1 * a - di2 * b, # H1
        di2 * a - di1 * b, # H2
        (di1 + di) * a + di1 * b, # H3
        (di2 + di) * a + di2 * b, # H4
        di * (a * 2.0 + b), # (H5)
        np.array([0.0, 0.0, 0.0]), # C5
        di1 * a, # C6
        di2 * a, # C7
        (di1 + di) * a, # C8
        (di2 + di) * a, # C9
        di * a * 2.0] # C10
    def translate(position):
        return position_carbon_0 + position
    positions = map(translate, positions_template)

    end_H_index, end_C_index = next_carbon_num - 1, next_carbon_num + 5
    end_point = positions[end_C_index]
    end_direction = positions[end_H_index] - end_point
    end_direction /= np.linalg.norm(end_direction)
    next_joint = {'end_point': end_point, 'direction': end_direction}
    positions.pop(end_H_index)

    elements = ['H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C']
    return (zip(positions, elements), next_joint)

def make_conj_poly(num_monomers, meta_length):
    plane_normal = np.array([0.0, 1.0, 0.0])
    ch_length = 1.06300  # C-H length in benzene.
    cc_length = 1.46400  # C-C length between benzene and C-C triple bond.

    structure = [(np.array([0.0, 0.0, 0.0]), 'H')]
    first_joint = {'end_point': structure[0][0], 'direction': np.array([-np.sin(math.pi / 6.0), 0.0, np.cos(math.pi / 6.0)])}
    benzene, joint = add_benzene(first_joint, ch_length, plane_normal, 5)
    structure.extend(benzene)
    group_ids = [0] * 11  # H(1) + benzene(10)

    for monomer in range(1, num_monomers):
        cc, joint = add_carbons_triple_bond(joint, cc_length)
        structure.extend(cc)
        if monomer == 0 or meta_length == 0:
            next_carbon_num = 5
        else:
            if monomer % (meta_length * 2) == 0:
                next_carbon_num = 4
            elif monomer % meta_length == 0:
                next_carbon_num = 3
            else:
                next_carbon_num = 5
        benzene, joint = add_benzene(joint, cc_length, plane_normal, next_carbon_num)
        structure.extend(benzene)
        group_ids.extend([max(group_ids) + 1] * 12)  # C-C(2) + benzene(10)

    hydrogen, joint = add_terminal_hydrogen(joint, ch_length)
    structure.extend(hydrogen)
    group_ids.append(max(group_ids))
    return structure, group_ids

def get_max_group_size(group_id):
    num_groups = max(group_id) + 1
    group_sizes = [0] * num_groups
    for i in group_id:
        group_sizes[i] += 1
    return max(group_sizes)

def print_group_ids(group_ids, filename):
    num_atoms = len(structure)
    num_groups = max(group_ids) + 1
    max_group_size = get_max_group_size(group_ids)
    with open(filename, 'w') as fp:
        fp.write('# number of atoms, number of groups, maximum number of group atoms\n')
        fp.write('%d %d %d\n' % (num_atoms, num_groups, max_group_size))
        for i in range(num_atoms):
            fp.write('%d %d\n' % (i + 1, group_ids[i] + 1))

def print_structure(structure, filename):
    template = '%s  %+10.6f  %+10.6f  %+10.6f\n'
    with open(filename, 'w') as fp:
        fp.write('%d\n' % len(structure))
        fp.write('# generated by %s\n' % ' '.join(sys.argv))
        for atom in structure:
            fp.write(template % (atom[1], atom[0][0], atom[0][1], atom[0][2]))

def get_filenames(num_monomers, meta_length):
    str_nums = ("meta_nm%d_m%d") % (num_monomers, meta_length)
    return (str_nums + ".xyz", str_nums + "_group_id.txt")

if __name__=='__main__':
    if len(sys.argv) <= 1:
        print 'Usage: python PiConjPolyMakerMeta.py <num_monomers> [<meta_length>]'
        sys.exit(0)

    num_monomers = int(sys.argv[1])
    if len(sys.argv) >= 3:
        meta_length = int(sys.argv[2])
    else:
        meta_length = 1  # meta_length = 0 means para structure.

    if meta_length != 0 and (num_monomers - 1) % (2 * meta_length) != 0:
        num_monomers_example_up = (num_monomers - 1) / (2 * meta_length) * (2 * meta_length) + 1
        num_monomers_example_down = ((num_monomers - 1) / (2 * meta_length) + 1) * (2 * meta_length) + 1
        print "[Warning] meta structure is truncated. The number of monomers must be " + \
        "2 * meta_length + 1 (For example, %d or %d)." % (num_monomers_example_up, num_monomers_example_down)

    structure, group_ids = make_conj_poly(num_monomers, meta_length)
    filename_xyz, filename_group_id = get_filenames(num_monomers, meta_length)
    print "filenames:", filename_xyz, filename_group_id
    print_structure(structure, filename_xyz)
    print_group_ids(group_ids, filename_group_id)

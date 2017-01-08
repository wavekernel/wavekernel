# -*- coding: utf-8 -*-
import json, sys, re, os.path, math
kAuPerAngstrom = 1.8897259885789

def read_xyz(fp):
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    xyz = []
    is_first_line = True
    for line in fp:
        if is_first_line:
            m = re.search(r"(\d+)", line)
            if m:
                num_atoms = int(m.group(1))
                is_first_line = False
            else:
                sys.exit(1)
        else:
            m = re.match(r"(?P<a>\w+)\s+(?P<x>%s)\s+(?P<y>%s)\s+(?P<z>%s)" %
                         (regexp_float, regexp_float, regexp_float), line)
            if m:
                atom = m.group("a")
                pos = [float(m.group("x")) * kAuPerAngstrom,
                       float(m.group("y")) * kAuPerAngstrom,
                       float(m.group("z")) * kAuPerAngstrom]
                xyz.append((atom, pos))
        if len(xyz) >= num_atoms:
            break
    return xyz

def read_group_id(fp):
    is_nums_read = False
    for line in fp:
        if is_nums_read:
            m = re.search("(\d+)\s+(\d+)", line)
            if m:
                atom_to_group[int(m.group(1)) - 1] = int(m.group(2)) - 1
        else:
            m = re.search("(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                is_nums_read = True
                num_atoms = int(m.group(1))
                num_groups = int(m.group(2))
                num_max_group_size = int(m.group(3))
                atom_to_group = [0] * num_atoms
    return {"num_groups": num_groups,
            "atom_to_group": atom_to_group}

def get_normal(r1, r2, r3):
    (x1, y1, z1) = (r2[0] - r1[0], r2[1] - r1[1], r2[2] - r1[2])
    (x2, y2, z2) = (r3[0] - r1[0], r3[1] - r1[1], r3[2] - r1[2])
    assert(x1 != 0.0 and x1 * y2 != x2 * y1)
    zn = 1.0
    yn = (x2 * z1 - x1 * z2) / (x1 * y2 - x2 * y1) * zn
    xn = - (y1 * yn + z1 * zn) / x1
    rn = math.sqrt(xn * xn + yn * yn + zn * zn)
    return [xn / rn, yn / rn, zn / rn]

# 0 <= theta <= pi / 2  ->   1 >= cos(theta) >= 0
def get_cos_dihedral_angle(carbons1, carbons2):
    n1 = get_normal(carbons1[0], carbons1[1], carbons1[2])
    n2 = get_normal(carbons2[0], carbons2[1], carbons2[2])
    return abs(n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2])

def make_group_to_atoms(xyz, group_id):
    group_to_atoms = []
    for i in range(group_id["num_groups"]):
        group_to_atoms.append([])
    for i in range(len(xyz)):
        g = group_id["atom_to_group"][i]
        group_to_atoms[g].append(xyz[i])
    return group_to_atoms

def get_three_carbon_positions(atoms):
    carbons = filter(lambda a: a[0] == "C", atoms)
    carbons.sort(key=lambda a: a[1][0])
    return (carbons[0][1], carbons[3][1], carbons[4][1])

if __name__ == "__main__":
    xyz_path = sys.argv[1]
    group_id_path = sys.argv[2]
    with open(xyz_path, "r") as fp:
        xyz = read_xyz(fp)
    with open(group_id_path, "r") as fp:
        group_id = read_group_id(fp)
    group_to_atoms = make_group_to_atoms(xyz, group_id)
    group_to_carbon_positions = map(get_three_carbon_positions, group_to_atoms)
    cos_dihedral_angles = []
    for i in range(len(group_to_carbon_positions) - 1):
        carbons1 = group_to_carbon_positions[i]
        carbons2 = group_to_carbon_positions[i + 1]
        c = get_cos_dihedral_angle(carbons1, carbons2)
        cos_dihedral_angles.append(c)
    print sum(cos_dihedral_angles) / len(cos_dihedral_angles)

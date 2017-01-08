# -*- coding: utf-8 -*-
# XYZファイルからX軸で等間隔に区切ったグループIDファイルを生成する.
import json, pylab, sys, re, os

regexp_float = r"-?(\d+\.)?\d+([deE][+-]?\d+)?"

def read_xyz(fp):
    i = 0
    line_to_x = []
    for line in fp:
        m = re.match(r"(\w+)\s+(%s)" % regexp_float, line)
        if m:
            line_to_x.append(float(m.group(2)))
    return line_to_x

def group_xs(line_to_x, num_groups):
    min_x = min(line_to_x)
    max_x = max(line_to_x) * 1.000001
    def x_to_group(x):
        return int(num_groups * (x - min_x) / (max_x - min_x))
    return map(x_to_group, line_to_x)

def print_group(atom_to_group):
    num_atoms = len(atom_to_group)
    print "# Made by set_group_id_divide.py"
    print "# number of atoms, number of groups, maximum number of group atoms"
    print "%d %d %d" % (num_atoms, max(atom_to_group) + 1, 0)
    print "#"
    for i in range(0, num_atoms):
        print "%d %d" % (i + 1, atom_to_group[i] + 1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python set_group_id_divide.py <xyz file> <n>"
    else:
        num_groups = int(sys.argv[2])
        with open(sys.argv[1], "r") as fp:
            atom_to_group = group_xs(read_xyz(fp), num_groups)
            print_group(atom_to_group)

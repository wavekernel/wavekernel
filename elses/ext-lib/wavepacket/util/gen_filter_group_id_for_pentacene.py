# -*- coding: utf-8 -*-
import sys, re

if __name__ == '__main__':
    num_atoms_all = 2592
    num_atoms_per_molecule = 36
    assert(num_atoms_all % num_atoms_per_molecule == 0)
    num_real_groups = num_atoms_all / num_atoms_per_molecule
    num_groups = int(sys.argv[1])
    assert(num_real_groups % num_groups == 0)
    filename = 'filter_group_id_group%d' % num_groups
    print filename
    with open(filename, 'w') as fp:
        fp.write('# pentacene film in %d groups\n' % num_groups)
        fp.write('# number of atoms, number of groups, maximum number of group atoms\n')
        fp.write('%d %d %d\n#\n' % (num_atoms_all, num_groups, num_atoms_all / num_groups))
        atom = 1
        for g in range(num_groups):
            for a in range(num_atoms_all / num_groups):
                fp.write('%d %d\n' % (atom, g + 1))
                atom += 1

import sys, re

def read_vector(fp):
    xs = []
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    for line in fp:
        m = re.search(r"(\d+\s+)?(?P<x>%s)" % regexp_float, line)
        if m:
            xs.append(float(m.group("x")))
    return xs

def read_and_print_gro(fp_in, num_atoms_in_molecule, indices, fp_out):
    line_num = 0
    index_molecule = 0
    is_end = indices == []
    for line in fp_in:
        if line_num == 0:
            fp_out.write(line)
        elif line_num == 1:
            num_atoms_in = int(line)
            fp_out.write('%d\n' % (num_atoms_in_molecule * len(indices)))
        elif not is_end and 2 + num_atoms_in_molecule * indices[index_molecule] <= line_num and line_num < 2 + num_atoms_in_molecule * (indices[index_molecule] + 1):
            fp_out.write(line)
            if line_num == 1 + num_atoms_in_molecule * (indices[index_molecule] + 1):
                index_molecule += 1
                if index_molecule >= len(indices):
                    is_end = True
        elif line_num == num_atoms_in + 2:  # Box.
            fp_out.write(line)
        line_num += 1

if __name__ == '__main__':
    gro_in_filename = sys.argv[1]
    num_atoms_in_molecule = int(sys.argv[2])
    eigenvector_filename = sys.argv[3]
    threshold = float(sys.argv[4])
    gro_out_filename = re.sub('\.[^.]+$', ('_cut_threshold%.3f.gro' % threshold), gro_in_filename)

    with open(eigenvector_filename) as fp:
        eigenvector = read_vector(fp)
    indexed_eigenvector = zip(range(len(eigenvector)), eigenvector)
    elements_above_threshold = filter(lambda (i, e): abs(e) >= threshold, indexed_eigenvector)
    indices_above_threshold = map(lambda (i, e): i, elements_above_threshold)
    print '%d molecules above threshold' % len(elements_above_threshold)
    for (i, e) in elements_above_threshold:
        print i + 1, e  # 0-origin to 1-origin.

    with open(gro_in_filename) as fp_in:
        with open(gro_out_filename, 'w') as fp_out:
            read_and_print_gro(fp_in, num_atoms_in_molecule, indices_above_threshold, fp_out)

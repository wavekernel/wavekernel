import sys, re

def read_vector(fp):
    xs = []
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    for line in fp:
        m = re.search(r"(\d+\s+)?(\d+\s+)?(?P<x>%s)" % regexp_float, line)
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

def read_and_print_xyz(fp_in, num_atoms_in_molecule, indices, fp_out):
    is_first_line = True
    atom_line_num = 0
    index_molecule = 0
    for line in fp_in:
        if line[0] == '#':
            fp_out.write(line)
        elif is_first_line:
            ss = line.split()
            num_atoms_in = int(ss[0])
            ss[0] = str(num_atoms_in_molecule * len(indices))
            fp_out.write(' '.join(ss) + '\n')
            is_first_line = False
            if indices == []:
                return
        else:
            if num_atoms_in_molecule * indices[index_molecule] <= atom_line_num and \
               atom_line_num < num_atoms_in_molecule * (indices[index_molecule] + 1):
                fp_out.write(line)
                if atom_line_num == num_atoms_in_molecule * (indices[index_molecule] + 1) - 1:
                    index_molecule += 1
                    if index_molecule >= len(indices):
                        return
            atom_line_num += 1
            if atom_line_num >= num_atoms_in:
                return

def get_pratio(xs):
    sum2 = 0.0
    sum4 = 0.0
    for x in xs:
        sum2 += x ** 2.0
        sum4 += x ** 4.0
    return sum2 ** 2.0 / sum4

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print 'usage: python cut_molecules_on_eigenvector.py <gro/xyz> <#atoms in molecule> <eigenvector> <threshold>'
        sys.exit(0)
    input_filename = sys.argv[1]
    file_type = input_filename[-3 :]
    num_atoms_in_molecule = int(sys.argv[2])
    eigenvector_filename = sys.argv[3]
    threshold = float(sys.argv[4])

    assert(file_type == 'xyz' or file_type == 'gro')
    output_filename = re.sub('\.[^.]+$', ('_cut_threshold%.4f.%s' % (threshold, file_type)), input_filename)
    print 'output:', output_filename

    with open(eigenvector_filename) as fp:
        eigenvector = read_vector(fp)
        print 'pratio: ', get_pratio(eigenvector)
    indexed_eigenvector = zip(range(len(eigenvector)), eigenvector)
    elements_above_threshold = filter(lambda (i, e): abs(e) >= threshold, indexed_eigenvector)
    square_sum_above_threshold = sum(map(lambda (i, e): e ** 2.0, elements_above_threshold))
    indices_above_threshold = map(lambda (i, e): i, elements_above_threshold)
    print '%d molecules above threshold' % len(elements_above_threshold)
    print 'square sum of the eigenvector elements above threshold: ', square_sum_above_threshold
    for (i, e) in elements_above_threshold:
        print i + 1, e  # 0-origin to 1-origin.

    with open(input_filename) as fp_in:
        with open(output_filename, 'w') as fp_out:
            if file_type == 'gro':
                read_and_print_gro(fp_in, num_atoms_in_molecule, indices_above_threshold, fp_out)
            elif file_type == 'xyz':
                read_and_print_xyz(fp_in, num_atoms_in_molecule, indices_above_threshold, fp_out)

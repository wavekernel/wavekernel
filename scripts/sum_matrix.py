import re, sys, os, argparse

def read_matrix(fp):  # Note: do not concern symmetricity.
    elements = []
    regexp_float = r'[+-]?(\d+\.)?\d+([deE][+-]?\d+)?'
    is_first_line = True
    num_elements = 0
    for line in fp:
        if line == '' or line[0] == '%':
            continue
        elif is_first_line:
            m = re.search(r'(\d+)\s+(\d+)\s+(\d+)', line)
            if m:
                rows = int(m.group(1))
                cols = int(m.group(2))
                nnz = int(m.group(3))
                print 'header: ', rows, cols, nnz
                is_first_line = False
            else:
                assert(False)
        else:
            m = re.search(r'(\d+)\s+(\d+)\s+(?P<x>%s)' % regexp_float, line)
            i = int(m.group(1))
            j = int(m.group(2))
            x = float(m.group('x'))
            elements.append(((i, j), x))
            num_elements += 1
            if num_elements >= nnz:
                break
    return {'rows': rows, 'cols': cols, 'elements': elements}

def mmwrite(matrix, filename, symmetricity, comment):
    with open(filename, 'w') as fp:
        fp.write('%%%%MatrixMarket matrix coordinate real %s\n' % symmetricity)
        fp.write('%%%s\n' % comment)
        num_nonzeros = len(matrix['elements'])
        fp.write('%7d %7d %12d\n' % (matrix['rows'], matrix['cols'], num_nonzeros))
        indices = matrix['elements'].keys()
        indices.sort()
        for index in indices:
            x = matrix['elements'][index]
            fp.write('%7d %7d % .16e\n' % (index[0], index[1], x))  # The space between % and . is intended.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('matrices', metavar='MATRICES', type=str, nargs='+',
                        help='Input matrices in MatrixMarket Format')
    parser.add_argument('-o', metavar='OUTPUT', type=str, dest='output_path',
                        default='sum.mtx', help='')
    parser.add_argument('-s', action='store_true', dest='is_symmetric',
                        default=False, help='')
    args = parser.parse_args()

    assert(not args.output_path in args.matrices)

    summed_elements = {}
    for i, f in enumerate(args.matrices):
        print 'start reading %d / %d: %s' % (i + 1, len(sys.argv) - 1, f)
        with open(f) as fp:
            matrix = read_matrix(fp)
            if i == 0:
                rows = matrix['rows']
                cols = matrix['cols']
            else:
                assert(matrix['rows'] == rows)
                assert(matrix['cols'] == cols)
            for index, x in matrix['elements']:
                if not index in summed_elements:
                    summed_elements[index] = 0.0
                summed_elements[index] += x
    summed_matrix = {'rows': rows, 'cols': cols, 'elements': summed_elements}
    symmetricity = 'symmetric' if args.is_symmetric else 'general'
    mmwrite(summed_matrix, args.output_path, symmetricity, 'generated by sum_matrix.py')

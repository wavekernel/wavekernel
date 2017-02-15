import re, sys

def read(fp):
    data = []
    max_i = 0
    max_j = 0
    for line in fp:
        m = re.match('\w*\(\s*(\d+),\s*(\d+)\)=\s+([-+.DEde\d]+)', line)
        if m:
            i = int(m.group(1))
            j = int(m.group(2))
            max_i = max(i, max_i)
            max_j = max(j, max_j)
            data.append((i, j, m.group(3)))
    return data, max_i, max_j

if __name__ == '__main__':
    with open(sys.argv[1]) as fp:
        data, max_i, max_j = read(fp)
    header = '%%MatrixMarket matrix coordinate real general'
    print header
    print '%d %d %d' % (max_i, max_j, len(data))
    for (i, j, x) in data:
        print i, j, x
    


import sys

def is_init_procedure(line):
    return line.find('set_structure_data') >= 0 or \
        False

if __name__ == '__main__':
    with open(sys.argv[1]) as fp:  # Output of fipppx -A -I cpu,call,hwm
        i = -1
        for line in fp:
            if line.find('Procedures profile (Total thread cost basis)') >= 0:
                i = 0
            elif i >= 0:
                i += 1
            if  8 <= i and i <= 21:
                line = line.rstrip()
                if not is_init_procedure(line):
                    print line

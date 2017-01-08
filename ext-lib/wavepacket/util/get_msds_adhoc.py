# -*- coding: utf-8 -*-
import json, sys, re, os.path, struct
regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
if __name__ == '__main__':
    for f in sys.argv[1:]:
        sys.stderr.write(f + '\n')
        with open(f, 'r') as fp:
            num_line = 0
            count = -1
            for line in fp:
                num_line += 1
                if (num_line % 1 == 0):
                    sys.stderr.write('line ' + str(num_line))
                m = re.search(r'"time":\s*(%s)' % regexp_float, line)
                if m:
                    time = float(m.group(1))
                else:
                    m = re.search(r'charge_coordinate_msd', line)
                    if m:
                        count = 0
                    elif 0 <= count and count < 3:
                        count += 1
                    elif count == 3:
                        msd = float(line)
                        sys.stdout.write('%f %f\n' % (time, msd))
                        sys.stdout.flush()
                        count = -1

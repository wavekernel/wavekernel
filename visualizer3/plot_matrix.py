# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re
import numpy as np

def read_mat(fp, i1, i2):
    dim = i2 - i1 + 1
    mat = np.zeros((dim, dim))
    for line in fp:
        m = re.match(r'\w+\(\s*(\d+),\s*(\d+)\)=\s*([-+.\dED]+)', line)
        if m:
            i_str = m.group(1)
            j_str = m.group(2)
            x_str = m.group(3)
        else:
            line_split = line.split(' ')
            i_str = line_split[0]
            j_str = line_split[1]
            x_str = line_split[2]
        i = int(i_str) - 1
        j = int(j_str) - 1
        #if abs(float(line_split[2].replace('D', 'E'))) > 0.0001:
        #    print i+1, j+1, float(line_split[2].replace('D', 'E'))
        #if i > i2 and j > i2:
        #    return mat
        if (i1 - 1 <= i and i <= i2 - 1) and (i1 - 1 <= j and j <= i2 - 1):
            x_str = x_str.replace('D', 'E')
            try:
                x = float(x_str)
            except:
                x = 0.0
            mat[i - i1 + 1][j - i1 + 1] = abs(x)
    return mat

if __name__ == '__main__':
    i1 = 1
    i2 = 33
    dim = i2 - i1 + 1
    for input_filename in sys.argv[1 :]:
        with open(input_filename) as fp:
            mat = read_mat(fp, i1, i2)
        print(input_filename, ' read end')
        x = np.linspace(i1 - 1, i2, dim + 1)
        y = np.linspace(i1 - 1, i2, dim + 1)
        X, Y = np.meshgrid(x, y)
        pylab.pcolor(X, Y, mat, vmin=0.0, vmax=1.0)
        pylab.xlabel('row index')
        pylab.ylabel('col index')
        pylab.title("%s abs((Y'S'Y)_{i,j}) %d - %d" % (input_filename, i1, i2))
        pylab.xlim(i1 - 1, i2)
        pylab.ylim(i1 - 1, i2)
        pylab.colorbar()
        output_filename = re.sub('\.[^.]+$', '.png', input_filename)
        pylab.savefig(output_filename)
        pylab.clf()

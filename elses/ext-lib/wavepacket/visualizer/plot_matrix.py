# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re
import numpy as np

def read_mat(fp, i1, i2):
    dim = i2 - i1 + 1
    mat = np.zeros((dim, dim))
    for line in fp:
        line_split = line.split(' ')
        i = int(line_split[0]) - 1
        j = int(line_split[1]) - 1
        #if abs(float(line_split[2].replace('D', 'E'))) > 0.0001:
        #    print i+1, j+1, float(line_split[2].replace('D', 'E'))
        if i >= i2 and j >= i2:
            return mat
        if (i1 - 1 <= i and i < i2) and (i1 - 1 <= j and j < i2):
            x_str = line_split[2].replace('D', 'E')
            x = float(x_str)
            mat[i - i1 + 1][j - i1 + 1] = abs(x)
    return mat

if __name__ == '__main__':
    i1 = 1750
    i2 = 1797
    dim = i2 - i1 + 1
    with open(sys.argv[1]) as fp:
        mat = read_mat(fp, i1, i2)
    print 'read end'
    x = np.linspace(i1 - 1, i2, dim)
    y = np.linspace(i1 - 1, i2, dim)
    X, Y = np.meshgrid(x, y)
    pylab.pcolor(X, Y, mat)
    pylab.xlabel('row index')
    pylab.ylabel('col index')
    pylab.title("map of abs((Y'S'Y)_{i,j}) for index %d - %d" % (i1, i2))
    pylab.xlim(i1 - 1, i2)
    pylab.ylim(i1 - 1, i2)
    pylab.colorbar()
    pylab.savefig('mat.png')

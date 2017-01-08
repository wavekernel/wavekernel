# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re
import numpy as np
import scipy as sp
import scipy.io

if __name__ == '__main__':
    #i1 = 1750
    #i2 = 1797
    #dim = i2 - i1 + 1
    mat = abs(sp.io.mmread(sys.argv[1]).todense().A)
    dim = mat.shape[0]
    for i in range(dim):
        mat[i, i] = 0.0
    i1 = 0
    i2 = 300
    x = np.linspace(i1 - 0.5, i2 - 0.5, i2 - i1 + 1)
    y = np.linspace(i1 - 0.5, i2 - 0.5, i2 - i1 + 1)
    X, Y = np.meshgrid(x, y)
    pylab.pcolor(X, Y, mat[i1 : i2, i1 : i2], cmap="Greys")  # matplotlib.org/examples/color/colormaps_reference.html
    pylab.xlabel('row index')
    pylab.ylabel('col index')
    pylab.title("map of abs(X) for index %d - %d" % (i1, i2 - 1))
    pylab.xlim(i1 - 0.5, i2 - 0.5)
    pylab.ylim(i1 - 0.5, i2 - 0.5)
    cbar = pylab.colorbar()
    pylab.savefig('mat.png')

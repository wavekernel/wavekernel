# -*- coding: utf-8 -*-
import sys
import numpy as np
import numpy.random
import scipy as sp
import scipy.io
import scipy.sparse as sps

def gen_almost_eye(dim):
    mat = sps.lil_matrix(sps.eye(dim))
    offdiag = np.random.rand(dim - 1) / 100.0
    mat.setdiag(offdiag, k=1)
    mat.setdiag(offdiag, k=-1)
    filename = 'almost_eye_dim%d.mtx' % dim
    return mat, filename

if __name__ == '__main__':
    dim = int(sys.argv[1])
    mat, filename = gen_almost_eye(dim)
    print mat
    scipy.io.mmwrite(filename, mat)

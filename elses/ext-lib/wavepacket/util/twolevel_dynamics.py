# -*- coding: utf-8 -*-
import numpy as np
import scipy as sp

if __name__ == "__main__":
    l1 = -0.4800133228958759
    l2 = -0.4882793176400123
    eps = 0.0200505186338377148

    l0 = (l1 + l2) / 2.0
    ld = (l2 - l1) / 2.0
    lp = l0 + np.sqrt(ld ** 2.0 + eps ** 2.0)
    lm = l0 - np.sqrt(ld ** 2.0 + eps ** 2.0)
    vp = np.array([eps, lp - l1])
    vm = np.array([eps, lm - l1])

    c1 = 1.0
    c2 = 0.0

    a = np.array([[eps, eps], [lp - l1, lm - l1]])
    b = np.array([c1, c2])
    x = np.linalg.solve(a, b)

    dt = 1.0

    flag = False

    for i in range(0, 100):
        t = dt * i
        c = x[0] * np.exp(-1.0j * lp * t) * vp + x[1] * np.exp(-1.0j * lm * t) * vm
        c_to_a = np.array([[np.exp(1.0j * l1 * t), 0.0], [0.0, np.exp(1.0j * l2 * t)]])
        a = np.dot(c_to_a, c)
        if flag:
            print 't: ', t
            print 'real', a.real
            print 'imag', a.imag
        else:
            print t, a[0].real, a[1].real

# -*- coding: utf-8 -*-
import numpy, sys
from scipy.stats import maxwell, bernoulli
import scipy.constants as const

mass_carbon = 2.1874661e4  # [a.u. (mass)]
mass_hydrogen = 1.8371526e3  # [a.u. (mass)]
kBoltzmann = 1.0  # Boltzmann's constant [a.u.]
kJoulePerEh = const.physical_constants['hartree-joule relationship'][0]  # [J / a.u. (energy)]
kBoltzmannInSI = const.physical_constants['Boltzmann constant'][0]  # [J / K]
kAuPerKelvin = kBoltzmannInSI / kJoulePerEh  # Temperature
kAccelRatio = 0.1
kDeltaTime = 1.0

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print 'Usage: python test_maxwell_dist.py <temperature (in Kelvin)> <atom (H or C)> <num of random variables>'
        sys.exit(0)
    temperature = float(sys.argv[1]) * kAuPerKelvin
    atom = sys.argv[2]  # 'H' or 'C'
    n = int(sys.argv[3])  # The number of generated random numbers.
    if atom == 'H':
        mass = mass_hydrogen
    elif atom == 'C':
        mass = mass_carbon
    else:
        assert(False)
    scale = numpy.sqrt(kBoltzmann * temperature / mass)
    speeds = [0.0] * n
    cumulations = [0.0] * n

    speeds[0] = maxwell.rvs(scale=scale)
    cumulations[0] = maxwell.cdf(speeds[0], scale=scale)
    for i in range(1, n):
        # bernoulli.rvs(p) = 1 in probability p, 0 in probability (1 - p).
        accelerate = (bernoulli.rvs(cumulations[i - 1]) == 0)
        if accelerate:
            speeds[i] = speeds[i - 1] * (1.0 + kDeltaTime * kAccelRatio)
        else:
            speeds[i] = max(0.0, speeds[i - 1] * (1.0 - kDeltaTime * kAccelRatio))
        # Cumulative probability P(x < speeds[i]) where x is a random variable
        # from the Maxwell distribution.
        cumulations[i] = maxwell.cdf(speeds[i], scale=scale)
    energies = map(lambda s: 0.5 * mass * (s ** 2.0), speeds)
    print 'speed [a.u.]    cumulation   energy [a.u.]'
    for i in range(n):
        print speeds[i], cumulations[i], energies[i]
    print 'temperature [a.u.]:', temperature
    print 'average energy [a.u.]:', sum(energies) / n

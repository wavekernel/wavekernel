import json, sys, pylab, re
import numpy as np
kPsecPerAu = 2.418884326505e-5  # Time.
kAngstromPerAu = 0.529177249  # Length
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0  # Length

def eval_distribution(z_angstrom, d):
    acc = 0.0
    z = z_angstrom / kAngstromPerAu
    for mean, msd, alpha_weight in zip(d['z_means'], d['msds'], d['alpha_weights']):
        acc += alpha_weight / np.sqrt(2 * np.pi * msd) * np.exp(-(z - mean) ** 2.0 / (2.0 * msd))
    return acc

def get_mean_and_msd_of_distribution(d):
    mean_distribution = 0.0
    msd_distribution = 0.0
    for mean, alpha_weight in zip(d['z_means'], d['alpha_weights']):
        mean_distribution += alpha_weight * mean
    for mean, msd, alpha_weight in zip(d['z_means'], d['msds'], d['alpha_weights']):
        msd_distribution += alpha_weight * ((mean - mean_distribution) ** 2.0 + msd)
    return (mean_distribution, msd_distribution)

def write_texts(d, z_min, y_max):
    mean_distribution, msd_distribution = get_mean_and_msd_of_distribution(d)
    pylab.text(z_min, y_max * 0.95, 'time: ' + str(alpha_distribution['time'] * kPsecPerAu) + '[ps]')
    pylab.text(z_min, y_max * 0.9, 'mean(dist): ' + str(mean_distribution * kAngstromPerAu) + '[A]')
    pylab.text(z_min, y_max * 0.85, 'msd(dist): ' + str(msd_distribution * kAngstrom2PerAu2) + '[A^2]')
    pylab.text(z_min, y_max * 0.8, 'msd(actual): ' + str(d['actual_msd'] * kAngstrom2PerAu2) + '[A^2]')

def write_point_texts(d):
    for eigenvalue, mean, msd, alpha_weight in zip(d['eigenvalues'], d['z_means'], d['msds'], d['alpha_weights']):
        if alpha_weight > 5e-2:
            pylab.text(mean * kAngstromPerAu, alpha_weight, '%.5f, %.2f' % (eigenvalue, msd * kAngstrom2PerAu2))

def plot_points(d):
    for m, w, i in zip(d['z_means'], d['alpha_weights'], d['indices']):
        m = m * kAngstromPerAu
        if i == d['fst_filter'] + d['num_filter'] - 1:
            color = 'red'
            print(i, color)
        else:
            color = 'green'
        pylab.plot([m, m], [0.0, w], 'x-', color=color)

for filename in sys.argv[1 :]:
    print('plotting', filename)
    with open(filename) as fp:
        alpha_distribution = json.load(fp)
    z_min = 0.0  # Angstrom.
    z_max = 700.0  # Angstrom.
    zs = np.linspace(z_min, z_max, num=700)
    ys = [eval_distribution(z, alpha_distribution) for z in zs]
    pylab.title(filename)
    #pylab.xlim(-0.473, -0.47)
    y_max = 0.02
    pylab.ylim(0.0, y_max)
    pylab.plot(zs, ys, '+')
    pylab.xlabel('Position [Angstrom]')
    pylab.ylabel('Probability density (dist)')
    write_texts(alpha_distribution, z_min, y_max)

    pylab.twinx()
    pylab.ylabel('Probability density (point)')
    pylab.ylim(0.0, 1.0)
    plot_points(alpha_distribution)
    write_point_texts(alpha_distribution)

    pylab.savefig(re.sub("\.[^.]+$", ".png", filename))
    pylab.clf()

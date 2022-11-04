import json, sys, pylab, re
import numpy as np
import matplotlib
kPsecPerAu = 2.418884326505e-5  # Time.
kAngstromPerAu = 0.529177249  # Length
kAngstrom2PerAu2 = kAngstromPerAu ** 2.0  # Length
kCoefWeightAmplitude = 10.
kAlphaWeightThreshold = 1e-8

#def eval_distribution(z_angstrom, d):
#    acc = 0.0
#    z = z_angstrom / kAngstromPerAu
#    for mean, msd, alpha_weight in zip(d['means'], d['msds'], d['alpha_weights']):
#        acc += alpha_weight / np.sqrt(2 * np.pi * msd) * np.exp(-(z - mean) ** 2.0 / (2.0 * msd))
#    return acc
#
#def get_mean_and_msd_of_distribution(d):
#    mean_distribution = 0.0
#    msd_distribution = 0.0
#    for mean, alpha_weight in zip(d['means'], d['alpha_weights']):
#        mean_distribution += alpha_weight * mean
#    for mean, msd, alpha_weight in zip(d['means'], d['msds'], d['alpha_weights']):
#        msd_distribution += alpha_weight * ((mean - mean_distribution) ** 2.0 + msd)
#    return (mean_distribution, msd_distribution)

def write_texts(d, x_min, y_min, y_max):
    #mean_distribution, msd_distribution = get_mean_and_msd_of_distribution(d)
    w = y_max - y_min
    pylab.text(x_min, y_min + w * 0.95, 'time: ' + str(d['time'] * kPsecPerAu) + '[ps]')
    pylab.text(x_min, y_min + w * 0.9, 'MSD(actual): ' + str(d['actual_msd'] * kAngstrom2PerAu2) + '[A^2]')

def write_point_texts(d):
    for eigenvalue, mean, msd, alpha_weight in zip(d['eigenvalues'], d['z_means'], d['msds'], d['alpha_weights']):
        if alpha_weight > 1e-2 or alpha_weight <= kAlphaWeightThreshold:
            pylab.text(eigenvalue, mean * kAngstromPerAu, '%.2f,%.1f' % (alpha_weight, msd * kAngstrom2PerAu2), fontsize=8)

#def plot_points(d):
#    for m, w, i in zip(d['means'], d['alpha_weights'], ):
#        m = m * kAngstromPerAu
#        if i == :
#            color = 'red'
#            print i, color
#        else:
#            color = 'green'
#        pylab.plot([m, m], [0.0, w], 'x-', color=color)

min_eigenvalue = 1e100
max_eigenvalue = -1e100
min_mean = 1e100
max_mean = -1e100
min_msd = 1e100
max_msd = -1e100
for filename in sys.argv[1 :]:
    print('reading', filename)
    with open(filename) as fp:
        d = json.load(fp)
    min_eigenvalue = min(min_eigenvalue, min(d['eigenvalues']))
    max_eigenvalue = max(max_eigenvalue, max(d['eigenvalues']))
    min_mean = min(min_mean, min(d['z_means']))
    max_mean = max(max_mean, max(d['z_means']))
    min_msd = min(min_msd, min(d['msds']))
    max_msd = max(max_msd, max(d['msds']))
min_mean = min_mean * kAngstromPerAu
max_mean = max_mean * kAngstromPerAu
min_msd = min_msd * kAngstrom2PerAu2
max_msd = max_msd * kAngstrom2PerAu2
width_log_msd = np.log(max_msd / min_msd)

for filename in sys.argv[1 :]:
    print('plotting', filename)
    with open(filename) as fp:
        d = json.load(fp)
    #z_min = 0.0  # Angstrom.
    #z_max = 700.0  # Angstrom.
    #zs = np.linspace(z_min, z_max, num=700)
    #ys = map(lambda z: eval_distribution(z, alpha_distribution), zs)
    pylab.title(filename)
    for e, z_mean, msd, w, i in zip(d['eigenvalues'], d['z_means'], d['msds'], d['alpha_weights'], d['indices']):
        ratio_log_msd = np.log(msd / min_msd) / width_log_msd
        if i == d['fst_filter'] + d['num_filter'] - 1:
            color = 'red'
        else:
            color = 'blue' #matplotlib.colors.hsv_to_rgb([ratio_log_msd, 1., 1.])
        if w > kAlphaWeightThreshold:
            fill = 'full'
            alpha = 1. - np.exp(- w * kCoefWeightAmplitude)
        else:
            fill = 'none'
            alpha = 1.
        pylab.plot([e], [z_mean * kAngstromPerAu], 'o', alpha=alpha, color=color, fillstyle=fill)
    pylab.xlabel('Eigenvalue [a.u.]')
    pylab.ylabel('Mean Z [A]')
    pylab.xlim(min_eigenvalue, max_eigenvalue)
    pylab.ylim(min_mean, max_mean)
    pylab.grid()
    write_texts(d, min_eigenvalue, min_mean, max_mean)
    write_point_texts(d)
    #pylab.xlim(-0.473, -0.47)
    #y_max = 0.1
    #pylab.ylim(0.0, y_max)
    #pylab.plot(zs, ys, '+')
    #pylab.xlabel('Position [Angstrom]')
    #pylab.ylabel('Probability density (dist)')
    #write_texts(alpha_distribution, z_min, y_max)

    #pylab.twinx()
    #pylab.ylabel('Probability density (point)')
    #pylab.ylim(0.0, 1.0)
    #plot_points(alpha_distribution)
    #write_point_texts(alpha_distribution)

    pylab.savefig(re.sub('\.[^.]+$', '', filename) + '.png')
    pylab.clf()

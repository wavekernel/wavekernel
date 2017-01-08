import json, sys, pylab, re
kPsecPerAu = 2.418884326505e-5  # Time
for filename in sys.argv[1 :]:
    with open(filename) as fp:
        data = json.load(fp)
    pylab.plot(data['eigenvalues'], data['alpha_nrm2'], 'x-')
    pylab.title('%s %.3f [ps]' % (filename, data['time'] * kPsecPerAu))
    pylab.xlim(-0.473, -0.47)
    pylab.ylim(0.0, 1.0)
    pylab.xlabel('energy [a.u.]')
    pylab.ylabel('amplitude')
    pylab.savefig(re.sub("\.[^.]+$", ".png", filename))
    pylab.clf()

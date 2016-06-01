import sys, re, pylab, math, matplotlib, json

#evPerAu = 27.211  # Energy.

def read_vector(fp):
    xs = []
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    for line in fp:
        m = re.search(r"\d+\s+(?P<x>%s)" % regexp_float, line)
        if m:
            xs.append(float(m.group("x")))
    return xs

if __name__=="__main__":
    #eigenvalues_path = sys.argv[1]
    #ipratios_path = sys.argv[2]
    #title = sys.argv[3]

    json_path = sys.argv[1]
    title = sys.argv[2]

    with open(json_path) as fp:
        data = json.load(fp)
    eigenvalues = data['condition']['eigenvalues']
    ipratios = data['condition']['eigenstate_ipratio']

    #with open(eigenvalues_path, "r") as fp:
    #    eigenvalues = read_vector(fp)
    #with open(ipratios_path, "r") as fp:
    #    ipratios = read_vector(fp)
    #homo_index = 1266  #len(ipratios) / 2 - 1
    #log_pratios = map(lambda x: -math.log10(x), ipratios)
    pratios = map(lambda x: 1.0 / x, ipratios)
    marker_style = matplotlib.markers.MarkerStyle(marker=None, fillstyle=u'full')
    pylab.xlabel("Eigenvalue [a.u.]")
    #pylab.ylabel("Log10 (Participation Ratio)")
    pylab.ylabel("Participation Ratio")
    #pylab.scatter(eigenvalues, log_pratios, facecolors='none', edgecolors='black', linewidth=0.5, label="Eigenstate")
    pylab.scatter(eigenvalues, pratios, facecolors='none', edgecolors='black', linewidth=0.5, label="Eigenstate")
    #pylab.scatter(eigenvalues[homo_index], log_pratios[homo_index], facecolors='none', edgecolors='b', linewidth=1.5, label="HOMO")
    #i = 1266
    #pylab.scatter(eigenvalues[i], pratios[i], facecolors='none', edgecolors='b', linewidth=1.5, label=str(i+1))
    #i = 1255
    #pylab.scatter(eigenvalues[i], pratios[i], facecolors='none', edgecolors='r', linewidth=1.5, label=str(i+1))
    pylab.legend(loc="upper left")
    pylab.savefig(title + "_pr.png")

# -*- coding: utf-8 -*-
# psiが書かれたファイルの集合から, 山﨑さんのAutoWaveで可視化を行える
# 単一のJSON形式ファイルに変換する.

import json, sys, re, os.path

def read(path):
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    xs = []
    with open(path, "r") as fp:
        for line in fp:
            m = re.search(r"\(\s*\d+\s*,\s*\d+\s*\)=\s*(?P<x>%s)" % regexp_float, line)
            if m:
                xs.append(float(m.group("x")))
            else:
                m = re.search(r"\s*\d+\s+\d+\s*(?P<x>%s)" % regexp_float, line)
                if m:
                    xs.append(float(m.group("x")))
    return xs

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python dats_to_json.py [dat]..."
        sys.exit(0)
    dat_paths = sys.argv[1:]
    dats = map(read, dat_paths)
    dat_dims = map(len, dats)
    assert(min(dat_dims) == max(dat_dims))
    dim = dat_dims[0]
    states = []
    for dat in dats:
        psi = {"real": dat, "imag": [0.0] * dim}
        states.append({"psi": psi})
    out = {"condition": {"dim": dim},
           "states": states}
    with open("out.json", "w") as fp:
        json.dump(out, fp)

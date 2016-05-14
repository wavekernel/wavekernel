# -*- coding: utf-8 -*-
import json, sys, re, os.path

def ipr(xs):
    sum_power4 = 0.0
    sum_power2 = 0.0
    for x in xs:
        sum_power4 += x ** 4.0
        sum_power2 += x ** 2.0
    return sum_power4 / (sum_power2 ** 2.0)

if __name__ == "__main__":
    regexp_float = r"[+-]?(\d+\.)?\d+([deE][+-]?\d+)?"
    xs = []
    with open(sys.argv[1], "r") as fp:
        for line in fp:
            #print line
            m = re.search(r"(\(?\s*\d+[\s,]+\d+\s*(\)=)?\s*)?(?P<x>%s)" % regexp_float, line)
            if m:
                xs.append(float(m.group("x")))
        #print xs
    print ipr(xs)

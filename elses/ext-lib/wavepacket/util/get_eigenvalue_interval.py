# -*- coding: utf-8 -*-
import json, sys

def get_interval(xs):
    assert(len(xs) >= 2)
    min_interval = max_interval = xs[1] - xs[0]
    for i in range(1, len(xs) - 1):
        interval = xs[i + 1] - xs[i]
        min_interval = min(min_interval, interval)
        max_interval = max(max_interval, interval)
    return (min_interval, max_interval)


if __name__ == "__main__":
    out_json_path = sys.argv[1]
    with open(out_json_path, "r") as fp:
        out = json.load(fp)
    eigenvalues = out["condition"]["eigenvalues"]
    print get_interval(eigenvalues)

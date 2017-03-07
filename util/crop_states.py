# -*- coding: utf-8 -*-

import json, sys, os.path, re, copy

def crop(wavekernel_out, stride):
    wavekernel_out_copy = copy.deepcopy(wavekernel_out)
    num_ticks = len(wavekernel_out["states"])
    wavekernel_out_copy["states"] = []
    for i in range(0, num_ticks, stride):
        wavekernel_out_copy["states"].append(wavekernel_out["states"][i])
    return wavekernel_out_copy

def crop_from_splits(out, out_dir, stride, step_num_first, step_num_last):
    states_joined = []
    for meta in out["split_files_metadata"]:
        fp = open(os.path.join(out_dir, meta["filename"]), "r")
        states_split = json.load(fp)
        fp.close()
        match = re.search("(\d{6,6})-(\d{6,6})", meta["filename"])
        first = max(int(match.group(1)), step_num_first)
        last = min(int(match.group(2)), step_num_last)
        for step in range(first, last + 1):
            if (step - step_num_first) % stride == 0:  # step number is 1-origin.
                states_joined.append(states_split[step - first])
    out["setting"]["is_output_split"] = False
    out.pop("split_files_metadata")
    out["states"] = states_joined
    return out

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: python crop_states.py <JSON file> <crop stride>"
        sys.exit(0)
    fp = open(sys.argv[1], 'r')
    wavekernel_out = json.load(fp)
    fp.close()
    stride = int(sys.argv[2])
    step_num_first = int(sys.argv[3])
    step_num_last = int(sys.argv[4])

    if wavekernel_out["setting"]["is_output_split"]:
        wavekernel_out_cropped = crop_from_splits(wavekernel_out,
                                                  os.path.dirname(sys.argv[1]),
                                                  stride, step_num_first, step_num_last)
    else:
        wavekernel_out_cropped = crop(wavekernel_out, stride)
    print json.dumps(wavekernel_out_cropped, indent=2)

# -*- coding: utf-8 -*-

import json, sys, os.path, re, copy

def print_times(wavekernel_out):
    timer_events = wavekernel_out["events"]
    max_event_name_length = max(map(lambda event:len(event["name"]), timer_events))
    timer_events.sort(key=lambda event:event["val"], reverse=True)
    print "event_name%s\tnum_repeated\ttime" % (" " * (max_event_name_length - len("event_name")))
    for event in timer_events:
        space_len = max_event_name_length - len(event["name"])
        print "%s%s\t%d\t%f" % (event["name"], " " * space_len, event["num_repeated"], event["val"])

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "Usage: python print_times.py <JSON file>"
        sys.exit(0)
    fp = open(sys.argv[1], 'r')
    wavekernel_out = json.load(fp)
    fp.close()
    print_times(wavekernel_out)

# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct
kAuPerAngstrom = 1.8897259885789
kSizeOfReal = 8

def get_real_array(split_dir, is_little_endian, element):
    if element == []:
        return []
    elif isinstance(element[0], basestring):  # Supposed to be binary output mode.
        first = element[1]
        last = element[2]
        count = last - first + 1
        with open(os.path.join(split_dir, element[0]), 'rb') as fp:
            fp.seek(kSizeOfReal * (first - 1), os.SEEK_SET)
            xs_str = fp.read(kSizeOfReal * count)
        format_char_endian = '<' if is_little_endian else '>'
        return struct.unpack(format_char_endian + str(count) + 'd', xs_str)
    else:
        return element

def add_step(state, extracted_types, split_dir, is_little_endian,
             ts, zzzs):
    # Common.
    ts.append(state["time"])
    # Alpha
    if "zzz" in extracted_types:
        psi_real = get_real_array(split_dir, is_little_endian, state["psi"]["real"])
        psi_imag = get_real_array(split_dir, is_little_endian, state["psi"]["imag"])
        zzs = []
        print 'start', state["time"]
        for i, (re, im) in enumerate(zip(psi_real, psi_imag)):
            j = i / 3594
            if len(zzs) <= j:
                zzs.append(0.0)
            zzs[j] += re ** 2.0 + im ** 2.0
        #print zzs
        sum2_2 = sum(zzs) ** 2.0
        sum4 = sum(map(lambda zz: zz ** 2.0, zzs))
        print len(zzs), sum2_2, sum4, sum4 / sum2_2
        zzzs.append(sum4 / sum2_2)

def calc(wavepacket_out, extracted_types, stride, wavepacket_out_path, is_little_endian, start_time):
    cond = wavepacket_out["condition"]
    # Common.
    dim = cond["dim"]
    ts = []
    zzzs = []

    assert(wavepacket_out["setting"]["is_output_split"])
    split_dir = os.path.dirname(wavepacket_out_path)
    for meta in wavepacket_out["split_files_metadata"]:
        path = os.path.join(split_dir, meta["filename"])
        with open(path, "r") as fp:
            diff = datetime.datetime.now() - start_time
            sys.stderr.write(str(diff) + " reading: " + path + "\n")
            states_split = json.load(fp)
        for state in states_split: #["states"]:
            if state["step_num"] % stride == 0:
                add_step(state, extracted_types, split_dir, is_little_endian,
                         ts, zzzs)

    header = re.sub("\.[^.]+$", "", wavepacket_out_path)
    if "zzz" in extracted_types:
        result_zzz = {"ts": ts,
                      "zzzs": zzzs}
        filename_zzz = header + "_zzz.json"
        with open(filename_zzz, "w") as fp:
            json.dump(result_zzz, fp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavepacket_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    parser.add_argument('--type', dest='extracted_types_comma_separated', type=str,
                        default='zzz', help='')
    args = parser.parse_args()

    extracted_types = args.extracted_types_comma_separated.split(',')
    assert(all(map(lambda t: t == "zzz", extracted_types)))

    start_time = datetime.datetime.now()

    if not os.path.isfile(args.wavepacket_out_path):
        sys.stderr.write("file " + args.wavepacket_out_path + " does not exist\n")
        sys.exit(1)

    with open(args.wavepacket_out_path, "r") as fp:
        wavepacket_out = json.load(fp)
    calc(wavepacket_out, extracted_types, args.skip_stride_num, args.wavepacket_out_path,
         args.is_little_endian, start_time)

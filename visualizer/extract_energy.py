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

def add_step(state, ts, tb_energy, nl_energy, total_energy):
    # Energies
    ts.append(state["time"])
    tb_energy.append(state["TB_energy"])
    nl_energy.append(state["NL_energy"])
    total_energy.append(state["total_energy"])

def calc(wavekernel_out, stride, wavekernel_out_path, is_little_endian, start_time):
    cond = wavekernel_out["condition"]
    # Energy
    dim = cond["dim"]
    eigenvalues = cond["eigenvalues"]
    ts = []
    tb_energy = []
    nl_energy = []
    total_energy = []

    if wavekernel_out["setting"]["is_output_split"]:
        split_dir = os.path.dirname(wavekernel_out_path)
        for meta in wavekernel_out["split_files_metadata"]:
            path = os.path.join(split_dir, meta["filename"])
            with open(path, "r") as fp:
                diff = datetime.datetime.now() - start_time
                sys.stderr.write(str(diff) + " reading: " + path + "\n")
                states_split = json.load(fp)
            for state in states_split:
                if state["step_num"] % stride == 0:
                    add_step(state, ts, tb_energy, nl_energy, total_energy)
    else:
        for state in wavekernel_out["states"]:
            if state["step_num"] % stride == 0:
                add_step(state, ts, tb_energy, nl_energy, total_energy)

    header = re.sub("\.[^.]+$", "", wavekernel_out_path)
    filename_energy = header + "_energy.json"
    result_energy = {"min_eigenvalue": min(eigenvalues),
                     "max_eigenvalue": max(eigenvalues),
                     "ts": ts,
                     "tb_energy": tb_energy,
                     "nl_energy": nl_energy,
                     "total_energy": total_energy}
    with open(filename_energy, "w") as fp:
        json.dump(result_energy, fp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    if not os.path.isfile(args.wavekernel_out_path):
        sys.stderr.write("file " + args.wavekernel_out_path + " does not exist\n")
        sys.exit(1)

    with open(args.wavekernel_out_path, "r") as fp:
        wavekernel_out = json.load(fp)

    calc(wavekernel_out, args.skip_stride_num, args.wavekernel_out_path,
         args.is_little_endian, start_time)

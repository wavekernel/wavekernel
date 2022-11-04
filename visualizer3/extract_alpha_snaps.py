# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct
kAuPerAngstrom = 1.8897259885789
kSizeOfReal = 8

def get_real_array(split_dir, is_little_endian, element):
    if element == []:
        return []
    elif isinstance(element[0], str):  # Supposed to be binary output mode.
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


def get_eigenvalues(split_dir, is_little_endian, out, input_step):
    for s in out["structures"]:
        if s["input_step"] == input_step:
            return get_real_array(split_dir, is_little_endian, s["eigenvalues"])
    assert(False)  # Specified input_step not found.

def read_step(state, split_dir, is_little_endian,
              num_filter, eigenvalues, header):
    t = state["time"]
    s = state["step_num"]
    alpha_real = get_real_array(split_dir, is_little_endian, state["alpha"]["real"])
    alpha_imag = get_real_array(split_dir, is_little_endian, state["alpha"]["imag"])
    alpha_nrm2 = [r_i[0] ** 2.0 + r_i[1] ** 2.0 for r_i in zip(alpha_real, alpha_imag)]
    assert(len(eigenvalues) == len(alpha_nrm2))
    output = {"time": t, "step_num": s, "eigenvalues": eigenvalues, "alpha_nrm2": alpha_nrm2}
    with open('%s_%06d_alpha.json' % (header, s), 'w') as fp:
        json.dump(output, fp)

def calc(wavekernel_out, stride, wavekernel_out_path, is_little_endian, start_time):
    cond = wavekernel_out["condition"]
    # Common.
    dim = cond["dim"]
    ts = []
    # Alpha
    fst_filter = wavekernel_out["setting"]["fst_filter"]
    num_filter = wavekernel_out["setting"]["end_filter"] - fst_filter + 1
    last_input_step = 0

    assert(wavekernel_out["setting"]["is_output_split"])
    split_dir = os.path.dirname(wavekernel_out_path)
    header = re.sub("\.[^.]+$", "", wavekernel_out_path)
    for meta in wavekernel_out["split_files_metadata"]:
        path = os.path.join(split_dir, meta["filename"])
        with open(path, "r") as fp:
            diff = datetime.datetime.now() - start_time
            sys.stderr.write(str(diff) + " reading: " + path + "\n")
            states_split = json.load(fp)
        for state in states_split["states"]:
            if state["step_num"] % stride == 0:
                eigenvalues = get_eigenvalues(split_dir, is_little_endian, states_split, state["input_step"])
                if state["input_step"] > last_input_step:
                    last_input_step += 1
                    read_step(state, split_dir, is_little_endian,
                              num_filter, eigenvalues, header)


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

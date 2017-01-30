# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, os, os.path
kAuPerAngstrom = 1.8897259885789  # Length.
kPsecPerAu = 2.418884326505e-5  # Time.
kSizeOfReal = 8
kAlphaWeightThreshold = 1e-8

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

def get_input_step_info(split_dir, is_little_endian, out, input_step):
    for s in out['structures']:
        if s['input_step'] == input_step:
            eigenvalues = get_real_array(split_dir, is_little_endian, s['eigenvalues'])
            z_means = get_real_array(split_dir, is_little_endian, s['eigenstate_mean_z'])
            msds = get_real_array(split_dir, is_little_endian, s['eigenstate_msd_total'])
            return (eigenvalues, z_means, msds)
    assert(False)  # Specified input_step not found.

def read_and_write_step(state, split_dir, is_little_endian,
                        fst_filter, num_filter, input_step_info, header, i, out_dir):
    t = state['time']
    s = state['step_num']
    actual_msd = state['charge_coordinate_msd'][3]
    alpha_real = get_real_array(split_dir, is_little_endian, state['alpha']['real'])
    alpha_imag = get_real_array(split_dir, is_little_endian, state['alpha']['imag'])
    alpha_weights = map(lambda (r, i): r ** 2.0 + i ** 2.0, zip(alpha_real, alpha_imag))
    (eigenvalues, z_means, msds) = input_step_info
    assert(len(eigenvalues) == len(alpha_weights) == num_filter)
    state_zipped = zip(range(fst_filter, fst_filter + num_filter), eigenvalues, z_means, msds, alpha_weights)
    def is_innegligible(state):
        is_highest = state[0] == fst_filter + num_filter - 1  # Save highest eigenstate forcibly.
        is_weight_large = state[4] > kAlphaWeightThreshold
        return is_highest or is_weight_large
    state_zipped_innegligible = filter(is_innegligible, state_zipped)

    (indices, eigenvalues, z_means, msds, alpha_weights) = map(list, zip(*state_zipped_innegligible))
    output = {'time': t, 'step_num': s, 'actual_msd': actual_msd,
              'indices': indices, 'fst_filter': fst_filter, 'num_filter': num_filter,
              'eigenvalues': eigenvalues, 'z_means': z_means, 'msds': msds, 'alpha_weights': alpha_weights}
    out_filename = '%06d.json' % i
    out_path = os.path.join(out_dir, out_filename)
    with open(out_path, 'w') as fp:
        json.dump(output, fp, indent=2)

def calc(wavepacket_out, stride, wavepacket_out_path, is_little_endian, out_dir, start_time, time_end):
    cond = wavepacket_out['condition']
    # Common.
    dim = cond['dim']
    ts = []
    # Alpha
    fst_filter = wavepacket_out['setting']['fst_filter']
    num_filter = wavepacket_out['setting']['end_filter'] - fst_filter + 1

    assert(wavepacket_out['setting']['is_output_split'])
    split_dir = os.path.dirname(wavepacket_out_path)
    header = re.sub('\.[^.]+$', '', wavepacket_out_path)
    i = 0
    for meta in wavepacket_out['split_files_metadata']:
        path = os.path.join(split_dir, meta['filename'])
        with open(path, 'r') as fp:
            diff = datetime.datetime.now() - start_time
            sys.stderr.write(str(diff) + ' reading: ' + path + '\n')
            states_split = json.load(fp)
        for state in states_split['states']:
            if state['step_num'] % stride == 0:
                input_step_info = get_input_step_info(split_dir, is_little_endian, states_split, state['input_step'])
                read_and_write_step(state, split_dir, is_little_endian,
                                    fst_filter, num_filter, input_step_info, header, i, out_dir)
                i += 1
                if (not time_end is None) and (state['time'] * kPsecPerAu >= time_end):
                    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavepacket_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('-e', metavar='TIME_END', dest='time_end', type=float, default=None,
                        help='[ps]')
    parser.add_argument('-d', metavar='OUT_DIR', dest='out_dir', type=str, default=None,
                        help='output directory')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    if not os.path.isfile(args.wavepacket_out_path):
        sys.stderr.write('file ' + args.wavepacket_out_path + ' does not exist\n')
        sys.exit(1)

    if args.out_dir is None:
        out_dir = re.sub('\.[^.]+$', '', args.wavepacket_out_path) + '_alpha_distribution'
    else:
        out_dir = args.out_dir
    if os.path.isdir(out_dir):
        print '[Warn] directory %s already exists' % out_dir
    else:
        os.mkdir(out_dir)

    with open(args.wavepacket_out_path, 'r') as fp:
        wavepacket_out = json.load(fp)
    calc(wavepacket_out, args.skip_stride_num, args.wavepacket_out_path,
         args.is_little_endian, out_dir, start_time, args.time_end)

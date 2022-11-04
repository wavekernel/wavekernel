# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, copy, pylab
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

def get_eigenvalues(split_dir, is_little_endian, out, input_step):  # Can read both from single or split output.
    for s in out['structures']:
        if s['input_step'] == input_step:
            return {'input_step': input_step,
                    'time': s['time'],
                    'eigenvalues': get_real_array(split_dir, is_little_endian, s['eigenvalues'])}
    assert(False)  # Specified input_step not found.

def calc(wavekernel_out, stride, wavekernel_out_path, is_little_endian, start_time):
    setting = wavekernel_out['setting']
    cond = wavekernel_out['condition']
    # Common.
    dim = cond['dim']
    # xyz
    eigenvalues_acc = []
    last_input_step = 0

    if wavekernel_out['setting']['is_output_split']:
        split_dir = os.path.dirname(wavekernel_out_path)
        for meta in wavekernel_out['split_files_metadata']:
            path = os.path.join(split_dir, meta['filename'])
            with open(path, 'r') as fp:
                diff = datetime.datetime.now() - start_time
                sys.stderr.write(str(diff) + ' reading: ' + path + '\n')
                states_split = json.load(fp)
            for state in states_split['states']:
                if (state['input_step'] - 1) % stride == 0 and state['input_step'] > last_input_step:
                    eigenvalues = get_eigenvalues(split_dir, is_little_endian, states_split, state['input_step'])
                    eigenvalues_acc.append(eigenvalues)
                    last_input_step = state['input_step']
    else:
        for state in wavekernel_out['states']:
            if (state['input_step'] - 1) % stride == 0 and state['input_step'] > last_input_step:
                eigenvalues = get_eigenvalues(split_dir, is_little_endian, states_split, state['input_step'])
                eigenvalues_acc.append(eigenvalues)
                last_input_step = state['input_step']
    return eigenvalues_acc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('-n', metavar='BINS', dest='num_bins', type=int, default=100,
                        help='')
    parser.add_argument('--xmin', metavar='XMIN', dest='x_min', type=float, default=None,
                        help='')
    parser.add_argument('--xmax', metavar='XMAX', dest='x_max', type=float, default=None,
                        help='')
    parser.add_argument('--ymin', metavar='YMIN', dest='y_min', type=float, default=None,
                        help='')
    parser.add_argument('--ymax', metavar='YMAX', dest='y_max', type=float, default=None,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    if not os.path.isfile(args.wavekernel_out_path):
        sys.stderr.write('file ' + args.wavekernel_out_path + ' does not exist\n')
        sys.exit(1)

    with open(args.wavekernel_out_path, 'r') as fp:
        wavekernel_out = json.load(fp)
    eigenvalues_acc = calc(wavekernel_out, args.skip_stride_num, args.wavekernel_out_path,
                           args.is_little_endian, start_time)
    acc = []
    acc_homo = []
    for e in eigenvalues_acc:
        acc.extend(e['eigenvalues'])
        acc_homo.append(e['eigenvalues'][-1])  # Suppose the highest eigenvalue is HOMO.
    min_eigenvalue = min(acc)
    max_eigenvalue = max(acc) + 1e-12
    width = (max_eigenvalue - min_eigenvalue) / args.num_bins
    header = re.sub('\.[^.]+$', '', args.wavekernel_out_path)
    output_filename = header + '_eigenhist.png'
    print('min_eigenvalue: ', min_eigenvalue)
    print('max_eigenvalue: ', max_eigenvalue)
    print('width: ', width)
    print('output_filename: ', output_filename)

    if args.x_min is None:
        x_min = min_eigenvalue
    else:
        x_min = args.x_min
    if args.x_max is None:
        x_max = max_eigenvalue
    else:
        x_max = args.x_max
    pylab.title(header)
    pylab.xlabel('Eigenvalue [au]')
    pylab.ylabel('Density of states')
    pylab.grid()
    hist_range = (x_min, x_max)
    pylab.hist(acc, bins=args.num_bins, normed=True, range=hist_range, label='All')
    pylab.hist(acc_homo, bins=args.num_bins, normed=True, range=hist_range, color='red', alpha=0.5, label='HOMO')
    pylab.xlim(x_min, x_max)
    pylab.ylim(args.y_min, args.y_max)
    pylab.legend(loc='best')
    pylab.savefig(output_filename)

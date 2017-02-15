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

def get_eigenvalue_vs_ipratio(split_dir, is_little_endian, out, input_step):  # Can read both from single or split output.
    for s in out['structures']:
        if s['input_step'] == input_step:
            return {'input_step': input_step,
                    'time': s['time'],
                    'eigenvalues': get_real_array(split_dir, is_little_endian, s['eigenvalues']),
                    'ipratios': get_real_array(split_dir, is_little_endian, s['eigenstate_ipratio_on_groups'])
            }
    assert(False)  # Specified input_step not found.

def calc(wavepacket_out, stride, wavepacket_out_path, is_little_endian, start_time):
    setting = wavepacket_out['setting']
    cond = wavepacket_out['condition']
    # Common.
    dim = cond['dim']
    # xyz
    vss_acc = []
    last_input_step = 0

    if wavepacket_out['setting']['is_output_split']:
        split_dir = os.path.dirname(wavepacket_out_path)
        for meta in wavepacket_out['split_files_metadata']:
            path = os.path.join(split_dir, meta['filename'])
            with open(path, 'r') as fp:
                diff = datetime.datetime.now() - start_time
                sys.stderr.write(str(diff) + ' reading: ' + path + '\n')
                states_split = json.load(fp)
            for state in states_split['states']:
                if (state['input_step'] - 1) % stride == 0 and state['input_step'] > last_input_step:
                    vs = get_eigenvalue_vs_ipratio(split_dir, is_little_endian, states_split, state['input_step'])
                    vss_acc.append(vs)
                    last_input_step = state['input_step']
    else:
        for state in wavepacket_out['states']:
            if (state['input_step'] - 1) % stride == 0 and state['input_step'] > last_input_step:
                vs = get_eigenvalue_vs_ipratio(split_dir, is_little_endian, states_split, state['input_step'])
                vss_acc.append(vs)
                last_input_step = state['input_step']
    return vss_acc

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavepacket_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    start_time = datetime.datetime.now()

    if not os.path.isfile(args.wavepacket_out_path):
        sys.stderr.write('file ' + args.wavepacket_out_path + ' does not exist\n')
        sys.exit(1)

    with open(args.wavepacket_out_path, 'r') as fp:
        wavepacket_out = json.load(fp)
    vss_acc = calc(wavepacket_out, args.skip_stride_num, args.wavepacket_out_path,
                   args.is_little_endian, start_time)
    output_path = re.sub('\.[^.]+$', '', args.wavepacket_out_path) + '_eigenvalue_vs_ipratio.json'
    with open(output_path, 'w') as fp:
        json.dump(vss_acc, fp, indent=2)

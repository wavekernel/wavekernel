# -*- coding: utf-8 -*-
import json, sys, re, os.path, struct
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

def make_restarter(wavekernel_out, wavekernel_out_path, is_little_endian):
    if wavekernel_out['setting']['is_output_split']:
        # Delete elements. They are discarded by wavekernel in restart mode if exist.
        wavekernel_out.pop('condition')
        wavekernel_out.pop('events')
        uninherited_settings = ['version', 'command', 'limit_t',
                                'num_mpi_processes', 'num_mpi_processes_row', 'num_mpi_processes_col',
                                'num_omp_threads']
        for key in uninherited_settings:
            wavekernel_out['setting'].pop(key)
        # Get the last state.
        split_dir = os.path.dirname(wavekernel_out_path)
        path_last = os.path.join(split_dir, wavekernel_out['split_files_metadata'][-1]['filename'])
        with open(path_last, 'r') as fp:
            states_split_last = json.load(fp)
        last_state = states_split_last[-1]
        uninherited_state_elements = ['alpha', 'charges_on_basis', 'charges_on_atoms',
                                      'charge_coordinate_mean', 'charge_coordinate_msd',
                                      'NL_energy', 'TB_energy', 'total_energy']
        for key in uninherited_state_elements:
            last_state.pop(key)
        last_state_psi_real = get_real_array(split_dir, is_little_endian, last_state['psi']['real'])
        last_state_psi_imag = get_real_array(split_dir, is_little_endian, last_state['psi']['imag'])
        last_state['psi'] = {'real': last_state_psi_real, 'imag': last_state_psi_imag}
        if wavekernel_out['setting']['h1_type'] == 'maxwell':
            last_state_atom_speed = get_real_array(split_dir, is_little_endian, last_state['atom_speed'])
            last_state['atom_speed'] = last_state_atom_speed
        elif wavekernel_out['setting']['h1_type'] == 'harmonic':
            last_state_atom_perturb = get_real_array(split_dir, is_little_endian, last_state['atom_perturb'])
            last_state['atom_perturb'] = last_state_atom_perturb
        elif wavekernel_out['setting']['h1_type'] == 'zero':  # Nothing to inherit.
            ()
        else:
            assert(False)
        wavekernel_out['states'] = [last_state]
        #
        match = re.search("(\d{6,6})-(\d{6,6})", path_last)
        restart_total_states_count = int(match.group(2))
        wavekernel_out['setting']['restart_total_states_count'] = restart_total_states_count
        wavekernel_out['setting']['restart_input_step'] = last_state['input_step']
    else:
        assert(False)

if __name__ == '__main__':
    wavekernel_out_path = sys.argv[1]
    with open(wavekernel_out_path, 'r') as fp:
        wavekernel_out = json.load(fp)

    if (len(sys.argv) >= 3 and sys.argv[2] == '--big-endian'):
        is_little_endian = False
    else:
        is_little_endian = True

    make_restarter(wavekernel_out, wavekernel_out_path, is_little_endian)
    with open('restart.json', 'w') as fp:
        json.dump(wavekernel_out, fp, indent=2)

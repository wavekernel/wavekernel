# -*- coding: utf-8 -*-
import argparse, json, sys, re, os, datetime, struct, copy
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

def get_charges_on_groups(group_info, charges_on_atoms):
    charges_on_groups = [0.0] * group_info["num_groups"]
    for atom in range(0, len(charges_on_atoms)):
        charges_on_groups[group_info["atom_to_group"][atom]] += charges_on_atoms[atom]
    return charges_on_groups

def make_dummy_group_id(num_atoms):
    return [[a] for a in range(1, num_atoms + 1)]  # One atom, one group.

def read_group_id(wavekernel_out):
    if "group_id" in wavekernel_out["condition"]:
        group_id = wavekernel_out["condition"]["group_id"]
    else:
        group_id = make_dummy_group_id(wavekernel_out["condition"]["num_atoms"])
    num_groups = len(group_id)
    num_atoms = len(wavekernel_out["condition"]["elements"])
    atom_to_group = [0] * num_atoms
    for g in range(num_groups):
        # Atom number in group_id is 1-oriented, while internal representation of visualiers are 0-oriented.
        for a in group_id[g]:
            atom_to_group[a - 1] = g
    return {"num_groups": num_groups,
            "atom_to_group": atom_to_group}

def get_xyz(split_dir, is_little_endian, out, input_step):  # Can read both from single or split output.
    for s in out["structures"]:
        if s["input_step"] == input_step:
            xs = get_real_array(split_dir, is_little_endian, s["coordinates_x"])
            ys = get_real_array(split_dir, is_little_endian, s["coordinates_y"])
            zs = get_real_array(split_dir, is_little_endian, s["coordinates_z"])
            return list(zip(xs, ys, zs))
    assert(False)  # Specified input_step not found.

def get_eigenvalues(split_dir, is_little_endian, out, input_step, first, last):  # Can read both from single or split output.
    for s in out['structures']:
        if s['input_step'] == input_step:
            return {'input_step': input_step,
                    'time': s['time'],
                    'eigenvalues': get_real_array(split_dir, is_little_endian, s['eigenvalues'])[first : last]}
    assert(False)  # Specified input_step not found.

def get_msd(charges, xyz, mean):
    assert(len(charges) == len(xyz))
    normalizer = 1.0 / sum(charges)
    msd = 0.0
    for atom in range(0, len(charges)):
        msd_atom = 0.0
        for axis in range(3):  # 0: x, 1: y, 2: z
            msd_atom += (xyz[atom][axis] - mean[axis]) ** 2.0
        msd += msd_atom * charges[atom] * normalizer
    return msd

def get_msd_contributions_on_groups(charges, xyz, mean, group_info, state):
    contributions = [0.0] * group_info["num_groups"]
    assert(len(charges) == len(xyz))
    normalizer = 1.0 / sum(charges)
    contributions = [0.0] * group_info["num_groups"]
    for atom in range(len(xyz)):
        msd_atom = 0.0
        for axis in range(3):  # 0: x, 1: y, 2: z
            msd_atom += (xyz[atom][axis] - mean[axis]) ** 2.0
        contributions[group_info["atom_to_group"][atom]] += msd_atom * charges[atom] * normalizer
    return contributions

def write_xyz(elements, xyz_coordinates, fp):
    ANGSTROM_PER_AU = 1.0 / kAuPerAngstrom
    to_write_header = False
    if to_write_header:
        header = [0.0, 0.0, 0.0]  # TODO: read box from ELSES.
    num_atoms = len(elements)
    x_min = y_min = z_min = 1e100
    for i in range(num_atoms):
        x = xyz_coordinates[0][i][0]
        y = xyz_coordinates[0][i][1]
        z = xyz_coordinates[0][i][2]
        x_min = min(x, x_min)
        y_min = min(y, y_min)
        z_min = min(z, z_min)
    for input_step, coord in enumerate(xyz_coordinates):
        assert(len(coord) == num_atoms)
        if to_write_header:
            fp.write('%d %.6f %.6f %.6f\n' % (num_atoms, header[0], header[1], header[2]))
        else:
            fp.write('%d\n' % num_atoms)
        fp.write('# input_step: %d\n' % input_step)
        for e, (x, y, z) in zip(elements, coord):
            s = '%s %.6f %.6f %.6f\n' % (e, x * ANGSTROM_PER_AU, y * ANGSTROM_PER_AU, z * ANGSTROM_PER_AU)
            fp.write(s)

def add_step(setting, state, extracted_types, split_dir, is_little_endian,
             num_filter, xyz, eigenvalues, group_info,
             ts, eigenvalues_log, tb_energy, nl_energy, total_energy,
             means, msds, ipratios, tb_energy_deviations,
             psi_pratios, alpha_pratios,
             charges_on_groups_all, charges_on_groups_max, msd_contributions_on_groups_max, alphas):
    # Common.
    ts.append(state["time"])
    # Energies
    if "energy" in extracted_types:
        eigenvalues_log.append(eigenvalues)
        tb_energy.append(state["TB_energy"])
        nl_energy.append(state["NL_energy"])
        total_energy.append(state["total_energy"])
    # Common for MSD and charge on groups.
    if "msd" in extracted_types or "group" in extracted_types:
        charges = get_real_array(split_dir, is_little_endian, state["charges_on_atoms"])
        mean = get_real_array(split_dir, is_little_endian, state["charge_coordinate_mean"])
        msd = get_real_array(split_dir, is_little_endian, state["charge_coordinate_msd"])
    # MSD
    if "msd" in extracted_types:
        msd_total_recalculated = get_msd(charges, xyz, mean)
        if abs(msd[3] - msd_total_recalculated) / msd[3] > 1e-4:
            print("msd error 1: ", msd[3], msd_total_recalculated)
        means.append(mean)
        msds.append(msd)
        tb_energy_deviations.append(state["TB_energy_deviation"])
        ipratios.append(state["psi"]["ipratio"])
    # Charge on groups
    if "group" in extracted_types:
        charges_on_groups = get_charges_on_groups(group_info, charges)
        msd_contributions_on_groups = get_msd_contributions_on_groups(charges, xyz, mean, group_info, state)
        charges_on_groups_max["val"] = max(charges_on_groups_max["val"], max(charges_on_groups))
        msd_contributions_on_groups_max["val"] = max(msd_contributions_on_groups_max["val"], max(msd_contributions_on_groups))
        if abs(msd[3] - sum(msd_contributions_on_groups)) / msd[3] > 1e-4:
            print("msd error 2: ", msd[3], sum(msd_contributions_on_groups))
        if setting["h1_type"] == "harmonic_for_nn_exciton":
            diag = copy.deepcopy(state["H_diag"])
            num_sites = len(diag) / 3
            for site in range(num_sites):
                diag[site * 3] += state["atom_perturb"][site * 2] + state["atom_perturb"][site * 2 + 1]
                if site < num_sites - 1:
                    diag[site * 3 + 1] += state["atom_perturb"][site * 2] + state["atom_perturb"][(site + 1) * 2 + 1]
                    diag[site * 3 + 2] += state["atom_perturb"][(site + 1) * 2] + state["atom_perturb"][site * 2 + 1]
            charges_on_groups_all.append({"n": state["step_num"],
                                          "t": state["time"],
                                          "ch": charges_on_groups,
                                          "co": msd_contributions_on_groups,
                                          "m": msd[3],
                                          "i": state["psi"]["ipratio"],
                                          "d": diag})
        else:
            charges_on_groups_all.append({"n": state["step_num"],
                                          "t": state["time"],
                                          "ch": charges_on_groups,
                                          "co": msd_contributions_on_groups,
                                          "m": msd[3],
                                          "i": state["psi"]["ipratio"]})
    # Alpha
    if "alpha" in extracted_types:
        alpha_real = get_real_array(split_dir, is_little_endian, state["alpha"]["real"])
        alpha_imag = get_real_array(split_dir, is_little_endian, state["alpha"]["imag"])
        for i in range(0, num_filter):
            real_part = alpha_real[i]
            imag_part = alpha_imag[i]
            if abs(real_part) < 1e-20:  # Alleviate underflow error.
                real_part = 0.
            if abs(imag_part) < 1e-20:
                imag_part = 0.
            absval = real_part ** 2 + imag_part ** 2
            if absval < 1e-6:
                absval = 0.0  # Reduce file size.
            alphas[i]["weights"].append(absval)
    # pratio
    if "pratio" in extracted_types:
        psi_pratios.append(1.0 / state["psi"]["ipratio"])
        alpha_pratios.append(1.0 / state["alpha"]["ipratio"])

def extract_main(wavekernel_out, extracted_types, stride, wavekernel_out_path, is_little_endian, out_dir, start_time):
    setting = wavekernel_out["setting"]
    cond = wavekernel_out["condition"]
    # Common.
    dim = cond["dim"]
    initial_eigenvalues = cond["eigenvalues"]
    # Select filtered eigenvalues from the largest one to max 100-th under it.
    last_eigenvalue_index = len(initial_eigenvalues)
    first_eigenvalue_index = max(0, last_eigenvalue_index - 100)
    ts = []
    # Charge on groups.
    charges_on_groups_all = []
    charges_on_groups_max = {"val": 0.0}
    msd_contributions_on_groups_max = {"val": 0.0}
    # Energy
    eigenvalues_log = []
    tb_energy = []
    nl_energy = []
    total_energy = []
    # MSD
    means = []
    msds = []
    ipratios = []
    tb_energy_deviations = []
    # Alpha
    fst_filter = wavekernel_out["setting"]["fst_filter"]
    num_filter = wavekernel_out["setting"]["end_filter"] - fst_filter + 1
    alphas = []
    # pratio
    psi_pratios = []
    alpha_pratios = []
    # xyz
    xyz_coordinates = []
    last_input_step = 0
    for i in range(num_filter):
        alphas.append({"i": i + fst_filter,
                       "eigenvalue": initial_eigenvalues[i],
                       "msd": [cond["eigenstate_msd_x"][i],
                               cond["eigenstate_msd_y"][i],
                               cond["eigenstate_msd_z"][i],
                               cond["eigenstate_msd_total"][i]],
                       "mean": [cond["eigenstate_mean_x"][i],
                                cond["eigenstate_mean_y"][i],
                                cond["eigenstate_mean_z"][i]],
                       "weights": []})

    group_info = read_group_id(wavekernel_out)
    if wavekernel_out["setting"]["is_output_split"]:
        split_dir = os.path.dirname(wavekernel_out_path)
        for meta in wavekernel_out["split_files_metadata"]:
            path = os.path.join(split_dir, meta["filename"])
            with open(path, "r") as fp:
                diff = datetime.datetime.now() - start_time
                sys.stderr.write(str(diff) + " reading: " + path + "\n")
                states_split = json.load(fp)
            for state in states_split["states"]:
                if state["step_num"] % stride == 0:
                    xyz = get_xyz(split_dir, is_little_endian, states_split, state["input_step"])
                    eigenvalues = get_eigenvalues(split_dir, is_little_endian, states_split, state['input_step'],
                                                  first_eigenvalue_index, last_eigenvalue_index)
                    if state["input_step"] > last_input_step:
                        xyz_coordinates.append(xyz)
                        last_input_step += 1
                    add_step(setting, state, extracted_types, split_dir, is_little_endian,
                             num_filter, xyz, eigenvalues, group_info,
                             ts, eigenvalues_log, tb_energy, nl_energy, total_energy,
                             means, msds, ipratios, tb_energy_deviations,
                             psi_pratios, alpha_pratios,
                             charges_on_groups_all, charges_on_groups_max, msd_contributions_on_groups_max,
                             alphas)
    else:
        for state in wavekernel_out["states"]:
            if state["step_num"] % stride == 0:
                xyz = get_xyz('', is_little_endian, wavekernel_out, state["input_step"])
                eigenvalues = get_eigenvalues('', is_little_endian, states_split, state['input_step'],
                                              first_eigenvalue_index, last_eigenvalue_index)
                if state["input_step"] > last_input_step:
                    xyz_coordinates.append(xyz)
                    last_input_step += 1
                add_step(setting, state, extracted_types, '', is_little_endian,
                         num_filter, xyz, eigenvalues, group_info,
                         ts, eigenvalues_log, tb_energy, nl_energy, total_energy,
                         means, msds, ipratios, psi_pratios, alpha_pratios,
                         charges_on_groups_all, charges_on_groups_max, msd_contributions_on_groups_max,
                         alphas)

    basedir, tail = os.path.split(wavekernel_out_path)
    if out_dir is not None:
        basedir = out_dir
    header = os.path.join(basedir, re.sub("\.[^.]+$", "", tail))

    if "group" in extracted_types:
        result_charge_group = {"num_groups": group_info["num_groups"],
                               "charges_on_groups_max": charges_on_groups_max["val"],
                               "charges_on_groups_all": charges_on_groups_all,
                               "msd_contributions_on_groups_max": msd_contributions_on_groups_max["val"]}
        filename_charge_group = header + "_charge_group.json"
        with open(filename_charge_group, "w") as fp:
            json.dump(result_charge_group, fp, indent=2)

    if "energy" in extracted_types:
        result_energy = {"min_eigenvalue": min(initial_eigenvalues),
                         "max_eigenvalue": max(initial_eigenvalues),
                         "ts": ts,
                         "eigenvalues_log": eigenvalues_log,
                         "tb_energy": tb_energy,
                         "nl_energy": nl_energy,
                         "total_energy": total_energy}
        filename_energy = header + "_energy.json"
        with open(filename_energy, "w") as fp:
            json.dump(result_energy, fp, indent=2)

    if "msd" in extracted_types:
        result_charge_moment = {"ts": ts,
                                "means": means,
                                "msds": msds,
                                "ipratios": ipratios,
                                "tb_energy_deviations": tb_energy_deviations}
        filename_charge_moment = header + "_charge_moment.json"
        with open(filename_charge_moment, "w") as fp:
            json.dump(result_charge_moment, fp, indent=2)

    if "alpha" in extracted_types:
        result_alpha = {"ts": ts,
                        "fst_filter": fst_filter,
                        "num_filter": num_filter,
                        "alphas": alphas}
        filename_alpha = header + "_alpha.json"
        with open(filename_alpha, "w") as fp:
            json.dump(result_alpha, fp, indent=2)

    if "pratio" in extracted_types:
        result_pratio = {"ts": ts,
                         "psi_pratios": psi_pratios,
                         "alpha_pratios": alpha_pratios}
        filename_pratio = header + "_pratio.json"
        with open(filename_pratio, "w") as fp:
            json.dump(result_pratio, fp, indent=2)

    if "xyz" in extracted_types:
        filename_xyz = header + "_position.xyz"
        with open(filename_xyz, "w") as fp:
            write_xyz(cond["elements"], xyz_coordinates, fp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('-o', metavar='OUTDIR', dest='out_dir', type=str, default='.',
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    parser.add_argument('--type', dest='extracted_types_comma_separated', type=str,
                        default='group,energy,msd,alpha,xyz,pratio', help='')
    args = parser.parse_args()

    extracted_types = args.extracted_types_comma_separated.split(',')
    assert(all([t == "group" or t == "energy" or t == "msd" or t == "alpha" or t == "xyz" or t == "pratio" for t in extracted_types]))

    start_time = datetime.datetime.now()

    # Argument checking.
    if not os.path.isfile(args.wavekernel_out_path):
        sys.stderr.write("File %s does not exist\n" % args.wavekernel_out_path)
        sys.exit(1)
    if not os.path.isdir(args.out_dir):
        sys.stderr.write("Directory %s does not exist\n" % args.out_dir)
        sys.exit(1)

    with open(args.wavekernel_out_path, "r") as fp:
        wavekernel_out = json.load(fp)
    extract_main(wavekernel_out, extracted_types, args.skip_stride_num, args.wavekernel_out_path,
                 args.is_little_endian, args.out_dir, start_time)

#===section0===#
import os,re,json,os.path,copy,sys,math, struct
from multiprocessing import Pool, cpu_count, current_process
import time
from itertools import chain

kPsecPerAu = 2.418884326505e-5  # Time
kFsecPerAu = 2.418884326505e-2  # Time
AuPerOng = 0.5291772489999785 
OngPerAu = 1.889725988579
kSizeOfReal = 8

def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)

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

def crop(wavepacket_out, stride,is_little_endian,max_num,min_num,target):
    wavepacket_out_copy = copy.deepcopy(wavepacket_out)
    num_ticks = len(wavepacket_out["states"])
    wavepacket_out_copy["states"] = []
    for states in wavepacket_out["states"]:
        if target == False:
            if not states["is_after_matrix_replace"] and \
               states["step_num"]%stride == 0 and \
               min_num <= states["step_num"] and states["step_num"] <= max_num:
                wavepacket_out_copy["states"].append(
                    get_elements(states,is_little_endian,""))
        else :
            for i in target:
                if not states["is_after_matrix_replace"] \
                   and states["step_num"] == i:
                    wavepacket_out_copy["states"].append(
                        get_elements(states,is_little_endian,""))
    return wavepacket_out_copy

def get_elements(states,is_little_endian,split_dir):
    elements = {}
    elements["time"] = states["time"]
    #sys.stderr.write(str(elements["time"])+ "\n")
    elements["step_num"] = states["step_num"]
    psi = {}
    psi["real"] = get_real_array(split_dir, is_little_endian,
                                 states["psi"]["real"])
    psi["imag"] = get_real_array(split_dir, is_little_endian,
                                 states["psi"]["imag"])
    elements["psi"] = psi
    elements["charge_coordinate_mean"] = states["charge_coordinate_mean"]
    elements["charge_coordinate_msd"] = states["charge_coordinate_msd"]
    elements["input_step"] = states["input_step"]
    return elements

def read_json(metadata):
    filename = metadata["filename"]
    first = metadata["min_step_num"]
    last = metadata["max_step_num"]
    states_joined = []
    if last < m_num or first > M_num:
        sys.stderr.write("skipping: " + filename + "\n")
    else:
        sys.stderr.write("reading: " + filename + "\n")
        with open(os.path.join(Split_dir, filename), "r") as fp:
            states_split = json.load(fp)
        for states in states_split["states"]:
            step = states["step_num"]
            if step % Stride == 0 and m_num <= step and step <= M_num and \
               not states["is_after_matrix_replace"]:
                #states_joined.append(states_split[step - first])
                states_joined += [get_elements(
                    states,Is_little_endian,Split_dir)]
    return states_joined

def read_json_target(metadata):
    filename = metadata["filename"]
    first = metadata["min_step_num"]
    last = metadata["max_step_num"]
    test = []
    for i in Target:
        if first <= i and i <= last :
            test += [i]
    states_joined = []
    if test == []:
        sys.stderr.write("skipping: " + filename + "\n")
    else:
        sys.stderr.write("reading: " + filename + "\n")
        fp = open(os.path.join(Split_dir, filename), "r")
        states_split = json.load(fp)
        fp.close()
        for states in states_split["states"]:
            step = states["step_num"]
            for i in test:
                if step == i :
                    states_joined += [get_elements(
                        states,Is_little_endian,Split_dir)]
    return states_joined

def crop_from_splits( out , stride , max_num , min_num ,
                     parallel_flag , target , core_num ,
                      is_little_endian , split_dir ):
    filenames = out["split_files_metadata"]
    states_joined = []
    global Stride,M_num,m_num,Target,Is_little_endian,Split_dir
    Stride = stride
    M_num = max_num
    m_num = min_num
    Target = target
    Is_little_endian = is_little_endian
    Split_dir = split_dir
    if parallel_flag :
        p = Pool(core_num)
        if Target == False:
            state_joined = p.map(read_json, filenames)
        else:
            state_joined = p.map(read_json_target, filenames)
    else :
        state_joined = []
        if Target == False:
            for i in filenames:
                state_joined += [read_json(i)]
        else:
            for i in filenames:
                state_joined += [read_json_target(i)]
    #print len(state_joined),type(state_joined)
    states_joined = []
    #states_joined = chain.from_iterable(state_joined)
    for states in state_joined:
        states_joined += states
    out["setting"]["is_output_split"] = False
    out.pop("split_files_metadata")
    out["states"] = states_joined
    return out

def print_wavefunction(LCAO_header, position_data,
                       states, out_dir,
                       name,multistep_input_mode,matrix_num):  # name = real/imag
    #m = 200
    out_dir_part = os.path.join(out_dir, name)
    if os.path.isdir(out_dir_part):
        pass
    else :
        os.mkdir(out_dir_part)

    out_paths = []
    t_index = 0
    for t in states:
        time = t["time"]
        out_path = os.path.join(
            out_dir_part, "psi_" + name + "_" + "%08d" % t_index +
            "_step%08d"%t["step_num"] +
            "_t%010.2ffs"%(t["time"]*kFsecPerAu) + ".txt")
        out_paths.append(out_path)
        out = open(out_path, 'w')
        out.write(LCAO_header[0])
        for i in range(len(LCAO_header)-2):
            if multistep_input_mode:
                mat_step = t["input_step"] - 1
                out.write(LCAO_header[i+1]%(
                    (position_data["xyz_t"]["x"][mat_step][i],
                     position_data["xyz_t"]["y"][mat_step][i],
                     position_data["xyz_t"]["z"][mat_step][i])))
            else:
                out.write(LCAO_header[i+1]%(
                    (position_data["xyz_t"]["x"][0][i],
                     position_data["xyz_t"]["y"][0][i],
                     position_data["xyz_t"]["z"][0][i])))
        out.write(LCAO_header[-1])
        out.write("         1\n")
        for j in range(matrix_num):
            out.write("        %d        %1.15f   k=         1\n" %
                      (j + 1, t["psi"][name][j]))
        t_index += 1
        out.close()
    return out_paths

def LCAO_converter(wave_out, position_data, out_dir,periodic,basis_path):
    header,matrix_num = read_basis(basis_path)#LCAO_header(position_data,periodic)
    states = wave_out["states"]
    base_number = wave_out["condition"]["dim"]
    if wave_out["setting"]["is_multistep_input_mode"]:
        multistep_input_mode = True
        sys.stderr.write("mode multi step %s [a.u.]\n"%\
                         wave_out["setting"]["multistep_input_read_interval"])
    else:
        multistep_input_mode = False
        sys.stderr.write("mode one step\n")
    out_paths = print_wavefunction(
        header,position_data, states, out_dir, "real",multistep_input_mode,matrix_num)
    out_paths.extend(print_wavefunction(
        header, position_data, states, out_dir, "imag",multistep_input_mode,matrix_num))
    return out_paths

def read_basis(filename):
    header = []
    with open(filename,"r") as f:
        lines = ""
        line = f.readline()# file_format= v0.04.05
        lines += line
        line = f.readline()# atom_num matrix_num
        lines += line
        line = str2list(line)
        atom_num = int(line[0])
        matrix_num = int(line[1])

        for i in range(3):# skip unit cell
            line = f.readline()
            lines += line
        header += [lines]

        atom_info = []
        for i in range(atom_num):# skip atomic orbital
            lines = ""
            line = f.readline() # atom info num
            lines += line
            line = str2list(line)
            atom = line[0]
            line = f.readline() # atom atomic orbital num
            lines += line
            line = str2list(line)
            orbital = int(line[0])
            atom_info += [(atom,orbital)]
            for j in range(6):
                line = f.readline()
                lines += line
            line = f.readline()
            lines += "      %f     %f    %f\n"
            header += [lines]

        line = f.readline() # LCAO coefficients
        header += [line]
        """
        line = f.readline() # vector num
        line = str2list(line)
        vector_num = int(line[0])

        matrix = make_zero_matrix(matrix_num)
        for i in xrange(vector_num):
            counter = 0
            for j in xrange(atom_num):
                for k in xrange(atom_info[j][1]):
                    line = f.readline()
                    line = str2list(line)
                    matrix[counter][i] = float(line[1])
                    counter += 1
        orbital_group = []
        for j in xrange(atom_num):
            for k in xrange(atom_info[j][1]):
                orbital_group += [j]
        """
        print(len(header)) 
    return header,matrix_num

def LCAO_header(position_data,periodic):
    header = []
    line = ""
    line +=  "# file_format= v0.04.05\n" 
    line +=  "        %d       %d\n" %\
            (position_data["atom_num"],position_data["matrix_num"])
    if periodic["x_mode"]:
        x = periodic["x"]
    else :
        x = position_data["input_mesh_grid"]["x_max"] - position_data["input_mesh_grid"]["x_min"]
    if periodic["y_mode"]:
        y = periodic["y"]
    else:
        y = position_data["input_mesh_grid"]["y_max"] - position_data["input_mesh_grid"]["y_min"]
    if periodic["z_mode"]:
        z = periodic["z"]
    else:
        z = position_data["input_mesh_grid"]["z_max"] - position_data["input_mesh_grid"]["z_min"]

    line +=  "     %f      0.00000000      0.00000000    %d\n" % (x,periodic["x_mode"]) 
    line +=  "      0.00000000     %f      0.00000000    %d\n" % (y,periodic["y_mode"])
    line +=  "      0.00000000      0.00000000     %f    %d\n" % (z,periodic["z_mode"])
    header += [line]
    for atom in position_data["atom"]:
        line = ""
        if atom == "H":
            line += " 1\n"
            line += "       1\n"
            line += "     1\n"
            line += "     0\n"
            line += " 1.300E+00\n"
            line += " 0.000E+00\n"
            line += " 0.000E+00\n"
            line += " 0.000E+00\n"
            line += "      %f     %f    %f\n"
        elif atom == "C":
            line += " 6\n"
            line += "       4\n"
            line += "     2     2     2     2\n"
            line += "     0     1     1     1\n"
            line += " 1.710E+00 1.625E+00 1.625E+00 1.625E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += "      %f     %f    %f\n"
        elif atom == "O":
            line += " 8\n"
            line += "       4\n"
            line += "     2     2     2     2\n"
            line += "     0     1     1     1\n"
            line += " 2.575E+00 2.275E+00 2.275E+00 2.275E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += " 0.000E+00 0.000E+00 0.000E+00 0.000E+00\n"
            line += "      %f     %f    %f\n"
        header += [line]
    header += ["# LCAO coefficients\n"]
    return header

if __name__ == "__main__":
    wave_out_path = sys.argv[1]
    output_wavefunction_path = sys.argv[2]
    out_dir = sys.argv[3]

    if len(sys.argv) > 4:
        stride = int(sys.argv[4])
    else:
        stride = 1

    wave_out_fp = open(wave_out_path, "r")
    wave_out = json.load(wave_out_fp)
    wave_out_fp.close()

    if wave_out["setting"]["is_output_split"]:
        print("In is_output_split")
        wave_out_cropped = crop_from_splits(
            wave_out,os.path.dirname(sys.argv[1]),stride)
    else:
        wave_out_cropped   = crop(wave_out, stride)
 
    LCAO_converter(wave_out_cropped, output_wavefunction_path, out_dir)

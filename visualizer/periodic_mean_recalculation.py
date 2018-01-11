# coding: utf-8
import sys, re, json, os.path, struct, argparse
import numpy as math
kSizeOfReal = 8
kAngstromPerAu = 0.529177249  # Length

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

def get_charges_on_atoms(wavepacket_out, split_dir, is_little_endian, stride, wavepacket_out_path):
    charges_on_atoms = []
    input_steps = []
    step_num = [] #
    readed_step = [] #
    structures_list = [] #
    needed_read_step = [] #
    ipratios = []
    tb_energy_deviations = []
    times = []
    
    if wavepacket_out["setting"]["is_output_split"]:
        split_dir = os.path.dirname(wavepacket_out_path)
        for meta in wavepacket_out["split_files_metadata"]:
            path = os.path.join(split_dir, meta["filename"])
            with open(path, "r") as fp:
                sys.stderr.write("reading: " + path + "\n")
                states_split = json.load(fp)
            for state in states_split["states"]:
                if state["step_num"] % stride == 0:
                    times.append(state["time"])
                    cs = get_real_array(split_dir, is_little_endian, state["charges_on_atoms"])
                    charges_on_atoms.append(cs)
                    ipratios.append(state["psi"]["ipratio"])
                    tb_energy_deviations.append(state["TB_energy_deviation"])
                    input_steps.append(state["input_step"])
                    step_num.append(state["step_num"])
                    if not state["input_step"] in needed_read_step :
                        needed_read_step.append(state["input_step"])
            
            for structures in states_split["structures"] :
                if not structures["input_step"] in readed_step :
                    if structures["input_step"] in needed_read_step :
                        dic = {}
                        x = get_real_array(split_dir, is_little_endian, structures["coordinates_x"])
                        y = get_real_array(split_dir, is_little_endian, structures["coordinates_y"])
                        z = get_real_array(split_dir, is_little_endian, structures["coordinates_z"])
                        dic["xyz"] = zip(wavepacket_out["condition"]["elements"],zip(x,y,z))
                        dic["input_step"] = structures["input_step"]
                        readed_step.append(structures["input_step"])
                        structures_list.append(dic)
                    else :
                        dic = {"xyz":None,"input_step":structures["input_step"]}
                        readed_step.append(structures["input_step"])
                        structures_list.append(dic)
    else:
        for state in wavepacket_out["states"]:
            if state["step_num"] % stride == 0:
                times.append(state["time"])
                cs = get_real_array(split_dir, is_little_endian, state["charges_on_atoms"])
                charges_on_atoms.append(cs)
                ipratios.append(state["psi"]["ipratio"])
                tb_energy_deviations.append(state["TB_energy_deviation"])
                input_steps.append(state["input_step"])
                step_num.append(state["step_num"])
                if not state["input_step"] in needed_read_step :
                    needed_read_step.append(state["input_step"])

        for structures in wavepacket_out["structures"] :
            if not structures["input_step"] in readed_step :
                if structures["input_step"] in needed_read_step :
                    dic = {}
                    x = get_real_array(split_dir, is_little_endian, structures["coordinates_x"])
                    y = get_real_array(split_dir, is_little_endian, structures["coordinates_y"])
                    z = get_real_array(split_dir, is_little_endian, structures["coordinates_z"])
                    dic["xyz"] = zip(wavepacket_out["condition"]["elements"],zip(x,y,z))
                    dic["input_step"] = structures["input_step"]
                    readed_step.append(structures["input_step"])
                    structures_list.append(dic)
                else :
                    dic = {"xyz":None,"input_step":structures["input_step"]}
                    readed_step.append(structures["input_step"])
                    structures_list.append(dic)

    return  charges_on_atoms, input_steps, structures_list, step_num, ipratios, tb_energy_deviations, times

def boundry(structures_list,Ls) :
    # 0~L の間に入るように周期境界条件をかける
    # Ls : 周期セル (x,y,z)
    # Ls = [8.0,8.0,2.5]
    new_structures_list = []
    for dic in structures_list :
        xyz = []
        for atom,(x,y,z) in dic["xyz"] :
            while Ls[0] < x   : x -= Ls[0] ;
            while x     < 0.0 : x += Ls[0] ;
            while Ls[1] < y   : y -= Ls[1] ;
            while y     < 0.0 : y += Ls[1] ;
            while Ls[2] < z   : z -= Ls[2] ;
            while z     < 0.0 : z += Ls[2] ;
            xyz.append( (atom,(x,y,z)) )
        new_dic = {"xyz":xyz,"input_step":dic["input_step"]}
        new_structures_list.append(dic)
    return new_structures_list

def calculation(charges_on_atoms,input_step,structures_list,step_num,Ls) :
    means = []
    msds = []

    old_mean_x = 0.5*Ls[0]
    old_mean_y = 0.5*Ls[1]
    old_mean_z = 0.5*Ls[2]
    
    for charge_on_atom,i in zip(charges_on_atoms,input_step) :
        xyz = structures_list[i-1]["xyz"]
        mean_x = 0.0
        mean_y = 0.0
        mean_z = 0.0
        for charge,(atom,(x,y,z)) in zip(charge_on_atom,xyz) :
            mean_x += charge*math.exp(1.0j*2.0*math.pi*x/Ls[0])
            mean_y += charge*math.exp(1.0j*2.0*math.pi*y/Ls[1])
            mean_z += charge*math.exp(1.0j*2.0*math.pi*z/Ls[2])
        mean_x = math.log(mean_x).imag * Ls[0]/(2.0*math.pi)
        mean_y = math.log(mean_y).imag * Ls[1]/(2.0*math.pi)
        mean_z = math.log(mean_z).imag * Ls[2]/(2.0*math.pi)

        # numpy.log は -0.5*L < sita < 0.5*L で出力されるため, 0 < sita < L になおす
        if mean_x < 0.0 : mean_x += Ls[0] ;
        if mean_y < 0.0 : mean_y += Ls[1] ;
        if mean_z < 0.0 : mean_z += Ls[2] ;

        # 一つ前の mean に近い値を mean にする
        while old_mean_x < mean_x : x_dummy = mean_x ; mean_x -= Ls[0] ;
        while mean_x < old_mean_x : x_dummy = mean_x ; mean_x += Ls[0] ;
        if abs(old_mean_x - mean_x) < abs(old_mean_x - x_dummy) :
            mean_x = mean_x
            old_mean_x = mean_x
        else :
            mean_x = x_dummy
            old_mean_x = x_dummy
        
        while old_mean_y < mean_y : y_dummy = mean_y ; mean_y -= Ls[1] ;
        while mean_y < old_mean_y : y_dummy = mean_y ; mean_y += Ls[1] ;
        if abs(old_mean_y - mean_y) < abs(old_mean_y - y_dummy) :
            mean_y = mean_y
            old_mean_y = mean_y
        else :
            mean_y = y_dummy
            old_mean_y = y_dummy
        
        while old_mean_z < mean_z : z_dummy = mean_z ; mean_z -= Ls[2] ;
        while mean_z < old_mean_z : z_dummy = mean_z ; mean_z += Ls[2] ;
        if abs(old_mean_z - mean_z) < abs(old_mean_z - z_dummy) :
            mean_z = mean_z
            old_mean_z = mean_z
        else :
            mean_z = z_dummy
            old_mean_z = z_dummy

        # charge の msd の計算
        msd_x = 0.0
        msd_y = 0.0
        msd_z = 0.0
        for charge,(atom,(x,y,z)) in zip(charge_on_atom,xyz) :
            # mean に近い値を使用 
            while mean_x < x : x_dummy = x ; x -= Ls[0] ;
            while x < mean_x : x_dummy = x ; x += Ls[0] ;
            x = min( [abs(mean_x - x),abs(mean_x - x_dummy)] )
            while mean_y < y : y_dummy = y ; y -= Ls[1] ;
            while y < mean_y : y_dummy = y ; y += Ls[1] ;
            y = min( [abs(mean_y - y),abs(mean_y - y_dummy)] )
            while mean_z < z : z_dummy = z ; z -= Ls[2] ;
            while z < mean_z : z_dummy = z ; z += Ls[2] ;
            z = min( [abs(mean_z - z),abs(mean_z - z_dummy)] )
            
            msd_x += charge*x**2
            msd_y += charge*y**2
            msd_z += charge*x**2
        msd = msd_x + msd_y + msd_z
        means.append( [mean_x,mean_y,mean_z] )
        msds.append( [msd_x,msd_y,msd_z,msd] )
    return means,msds

def wave_center_msd(means) :
    msds = [[0.0,0.0,0.0,0.0]]
    x_0 , y_0 , z_0 = means[0]
    for x,y,z in means[1:] :
        msd_x = (x - x_0)**2
        msd_y = (y - y_0)**2
        msd_z = (z - z_0)**2
        msd = msd_x + msd_y + msd_z
        msds.append( [msd_x,msd_y,msd_z,msd] )
    return msds

def output_json(wavepacket_out_path, charges_on_atoms, input_steps, structures, step_num):
    header = re.sub("\.[^.]+$", "", wavepacket_out_path)
    filename_charges_on_atoms = header+"_boundry_charge_moment.json"
    result_charges_on_atoms = {"msds":msds,
                               "ts":ts,
                               "ipratios":ipratios,
                               "means":means}
    with open(filename_charges_on_atoms, "w") as fp:
        json.dump(result_charges_on_atoms, fp, indent=2)

def read_boundry_L(filename) :
    filetype = filename[-4:]
    if filetype == ".xyz" :
        with open(filename) as fp :
            line = fp.readline()
        line = line.split()
        Ls = [float(x)/kAngstromPerAu for x in line[1:]]
    elif filetype == ".xml" :
        from xml.etree import ElementTree
        tree = ElementTree.parse(filename)
        structure = tree.getroot()
        unitcell = structure.find("unitcell")

        Ls = []
        i = 0
        for vector in unitcell.findall("vector") :
            if vector.attrib["unit"] == "a.u." :
                Ls.append( float(vector.text.split()[i]) )
            elif vector.attrib["unit"] == "angstrom" :
                Ls.append( float(vector.text.split()[i])/kAngstromPerAu )
            else :
                print "[error] %s is not suport" % vector.attrib["unit"]
                sys.exit()
            i += 1
    else :
        print "[error] %s is not suport"%(filetype)
    print Ls
    return Ls

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('wavekernel_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('xyz_path', metavar='[XYZ/XML]', type=str,
                        help='need boundry cell')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num', type=int, default=1,
                        help='')
    parser.add_argument('--big-endian', action='store_false', dest='is_little_endian',
                        default=True, help='')
    args = parser.parse_args()

    wavepacket_out_path = args.wavekernel_out_path
    with open(wavepacket_out_path, "r") as fp:
        wavepacket_out = json.load(fp)

    Ls = read_boundry_L(args.xyz_path)

    stride = args.skip_stride_num
    split_dir = os.path.dirname(wavepacket_out_path)
    is_little_endian = args.is_little_endian
    charges_on_atoms, input_steps, structures_list, step_num, ipratios, tb_energy_deviations, times = \
                      get_charges_on_atoms(wavepacket_out, split_dir, is_little_endian, stride, wavepacket_out_path)

    #structures_list = boundry(structures_list,Ls)
    means,msds = calculation(charges_on_atoms,input_steps,structures_list,step_num,Ls)
    wave_center_msds = wave_center_msd(means)

    header = re.sub("\.[^.]+$", "", wavepacket_out_path)
    filename_charges_on_atoms = header+"_periodic_charge_moment.json"
    result_charges_on_atoms = {"tb_energy_deviations":tb_energy_deviations,
                               "msds":msds,
                               "ts":times,
                               "ipratios":ipratios,
                               "means":means,
                               "wave_center_msds":wave_center_msds}
    with open(filename_charges_on_atoms, "w") as fp:
        json.dump(result_charges_on_atoms, fp, indent=2)
    

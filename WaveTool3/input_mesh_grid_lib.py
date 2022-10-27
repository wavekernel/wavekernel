# coding: utf-8
import re,sys,os,os.path,math

def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)

def read_xyz_t(filename,periodic):# for "step > 1"
    xyz_t = {}
    x_t = []
    y_t = []
    z_t = []
    with open(filename,"r") as fp:
        while True:
            x = []
            y = []
            z = []
            atom_num = fp.readline() # atom num
            if not atom_num:
                break
            atom_num = str2list(atom_num)
            if periodic["flag"]  != False :
                periodic["x_t"] += [float(atom_num[1])*1.889725988579]
                periodic["y_t"] += [float(atom_num[2])*1.889725988579]
                periodic["z_t"] += [float(atom_num[3])*1.889725988579]
            atom_num = int(atom_num[0])
            line = fp.readline() # comment
            for i in range(atom_num):
                line = fp.readline()
                line = str2list(line)
                if periodic["x_mode"] :
                    x += [math.fmod(float(line[1].replace("D","E"))*1.889725988579,periodic["x_t"][-1])]
                    if x[-1] < 0 :
                        x[-1] += periodic["x_t"][-1]
                else :
                    x += [float(line[1].replace("D","E"))*1.889725988579]
                if periodic["y_mode"] :
                    y += [math.fmod(float(line[2].replace("D","E"))*1.889725988579,periodic["y_t"][-1])]
                    if y[-1] < 0 :
                        y[-1] += periodic["y_t"][-1]
                else :
                    y += [float(line[2].replace("D","E"))*1.889725988579]
                if periodic["z_mode"] :
                    z += [math.fmod(float(line[3].replace("D","E"))*1.889725988579,periodic["z_t"][-1])]
                    if z[-1] < 0 :
                        z[-1] += periodic["z_t"][-1]
                else :
                    z += [float(line[3].replace("D","E"))*1.889725988579]
            x_t += [x]
            y_t += [y]
            z_t += [z]
    xyz_t["x"] = x_t
    xyz_t["y"] = y_t
    xyz_t["z"] = z_t
    print("step num",len(xyz_t["x"]))
    return xyz_t

def read_xyz(filename,periodic):
    xyz = {}
    x = []
    y = []
    z = []
    atom = []
    fp = open(filename,"r")
    atom_num = fp.readline() # atom num
    atom_num = int(str2list(atom_num)[0])
    line = fp.readline() # comment
    counter = 0
    H_counter = 0
    C_counter = 0
    O_counter = 0
    #print type(atom_num)
    for i in range(atom_num):
        line = fp.readline()
        if atom_num <= counter :
            break
        line = str2list(line)
        if line[0] == "H":
            H_counter += 1
        elif line[0] == "C":
            C_counter += 1
        elif line[0] == "O":
            O_counter += 1
        if periodic["x_mode"] :
            x += [math.fmod(float(line[1].replace("D","E"))*1.889725988579,periodic["x"])]
        else :
            x += [float(line[1].replace("D","E"))*1.889725988579]
        if periodic["y_mode"] :
            y += [math.fmod(float(line[2].replace("D","E"))*1.889725988579,periodic["y"])]
        else :
            y += [float(line[2].replace("D","E"))*1.889725988579]
        if periodic["z_mode"] :
            z += [math.fmod(float(line[3].replace("D","E"))*1.889725988579,periodic["z"])]
        else :
            z += [float(line[3].replace("D","E"))*1.889725988579]
        atom += [line[0]]
        counter += 1
    print("atom num" ,counter)
    fp.close()
    xyz["x"] = x
    xyz["y"] = y
    xyz["z"] = z
    data = {}
    data["xyz"] = xyz # xyz data
    data["xyz_t"] = read_xyz_t(filename,periodic) # xyz data (step > 1)
    data["C_num"] = C_counter # C atom num
    data["H_num"] = H_counter # H atom num
    data["O_num"] = O_counter # O atom num
    data["atom"] = atom # all atom data
    data["atom_num"] = atom_num # all atom num
    data["matrix_num"] = data["C_num"]*4 + data["H_num"]*1 + data["O_num"]*4
    #print "matrix num",data["matrix_num"]
    return data

def input_mesh_grid_t(xyz_t,cutoff):
    input_mesh_grid = {}
    input_mesh_grid["x_min"] = min(xyz_t["x"][0]) - cutoff
    input_mesh_grid["x_max"] = max(xyz_t["x"][0]) + cutoff
    input_mesh_grid["y_min"] = min(xyz_t["y"][0]) - cutoff
    input_mesh_grid["y_max"] = max(xyz_t["y"][0]) + cutoff
    input_mesh_grid["z_min"] = min(xyz_t["z"][0]) - cutoff
    input_mesh_grid["z_max"] = max(xyz_t["z"][0]) + cutoff
    for t in range(len(xyz_t["x"])):
        input_mesh_grid["x_min"] = min([input_mesh_grid["x_min"],
                                       min(xyz_t["x"][t]) - cutoff])
        input_mesh_grid["x_max"] = max([input_mesh_grid["x_max"],
                                       max(xyz_t["x"][t]) + cutoff])
        input_mesh_grid["y_min"] = min([input_mesh_grid["y_min"],
                                       min(xyz_t["y"][t]) - cutoff])
        input_mesh_grid["y_max"] = max([input_mesh_grid["y_max"],
                                       max(xyz_t["y"][t]) + cutoff])
        input_mesh_grid["z_min"] = min([input_mesh_grid["z_min"],
                                       min(xyz_t["z"][t]) - cutoff])
        input_mesh_grid["z_max"] = max([input_mesh_grid["z_max"],
                                       max(xyz_t["z"][t]) + cutoff])
    return input_mesh_grid

def main(a,xyz_path,mesh,periodic):
    data = read_xyz(xyz_path,periodic)

    #x_min = min(data["xyz"]["x"])
    #x_max = max(data["xyz"]["x"])
    #y_min = min(data["xyz"]["y"])
    #y_max = max(data["xyz"]["y"])
    #z_min = min(data["xyz"]["z"])
    #z_max = max(data["xyz"]["z"])

    data["mesh"] = mesh 
    data["xyz_path"] = xyz_path # xyz file name
    data["cutoff_au"] = a
    data["input_mesh_grid"] = input_mesh_grid_t(data["xyz_t"],a)
    #print_input_mesh_grid(x_min,x_max,y_min,y_max,z_min,z_max,mesh,a)
    #print data
    return data

def make_position_alpha(wave_out_cropped,position_data,alpha):
    mean = wave_out_cropped["states"][0]["charge_coordinate_mean"]
    msd = wave_out_cropped["states"][-1]["charge_coordinate_msd"]
    msd = math.sqrt(msd[3])*alpha
    input_mesh_grid = position_data["input_mesh_grid"]
    if mean[0] - msd < input_mesh_grid["x_min"]:
        x_min = input_mesh_grid["x_min"]
    else:
        x_min = mean[0] - msd
    if mean[0] + msd > input_mesh_grid["x_max"]:
        x_max = input_mesh_grid["x_max"]
    else:
        x_max = mean[0] + msd
    if mean[1] - msd < input_mesh_grid["y_min"]:
        y_min = input_mesh_grid["y_min"]
    else:
        y_min = mean[1] - msd
    if mean[1] + msd > input_mesh_grid["y_max"]:
        y_max = input_mesh_grid["y_max"]
    else:
        y_max = mean[1] + msd
    if mean[2] - msd < input_mesh_grid["z_min"]:
        z_min = input_mesh_grid["z_min"]
    else:
        z_min = mean[2] - msd
    if mean[2] + msd > input_mesh_grid["z_max"]:
        z_max = input_mesh_grid["z_max"]
    else:
        z_max = mean[2] + msd
    a = position_data["cutoff_au"]
    mesh = position_data["mesh"]
    print_input_mesh_grid(x_min,x_max,y_min,y_max,z_min,z_max,mesh,a)

def make_alpha(wave_out_cropped,alpha,mesh,a):
    mean = wave_out_cropped["states"][0]["charge_coordinate_mean"]
    msd = wave_out_cropped["states"][-1]["charge_coordinate_msd"]
    msd = math.sqrt(msd[3])*alpha
    x_min = mean[0] - msd
    x_max = mean[0] + msd
    y_min = mean[1] - msd
    y_max = mean[1] + msd
    z_min = mean[2] - msd
    z_max = mean[2] + msd
    print_input_mesh_grid(x_min,x_max,y_min,y_max,z_min,z_max,mesh,a)

def print_input_mesh_grid(x_min,x_max,y_min,y_max,z_min,z_max,mesh,a):
    with open("input_mesh_grid.txt","w") as f:
        f.write('%d %d %d\n'%(int((x_max - x_min+a*2)*mesh)+10,
                              int((y_max - y_min+a*2)*mesh)+10,
                              int((z_max - z_min+a*2)*mesh)+10))
        f.write('%f %f %f au\n'%(x_min-a,y_min-a,z_min-a))
        f.write('%f %f %f au\n'%(x_max - x_min+a*2,
                                 y_max - y_min+a*2,
                                 z_max - z_min+a*2))
    
if __name__=="__main__":
    if len(sys.argv) < 2:
        print("Usage: python make_input_mesh_grid.py <xyz file> <mesh> <a>")
        sys.exit(0)
    xyz_path = sys.argv[1]
    mesh = float(sys.argv[2])
    a = float(sys.argv[3])
    main(a,xyz_path,mesh)

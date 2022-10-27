#  import module and file open
import math
import re
import sys,os
from multiprocessing import Pool, cpu_count, current_process

def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)

def charge(cube_path):
    if cube_path.find("psi_real") >= 0:
        cube_imag_path = re.sub("real", "imag", cube_path)
        cube_char_path = re.sub("real", "char", cube_path)
        Charge_make(cube_path, cube_imag_path, cube_char_path)

def Charge_make(real_path, imag_path, out_path):
    with open(real_path, "r") as real_fp :
        with open(imag_path, "r") as imag_fp :

            set_list1 = []
            #  import module and file open
import math
import re
import sys
def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)
def Charge_make(real_path, imag_path, out_path):
    real_fp = open(real_path, "r")
    imag_fp = open(imag_path, "r")

    set_list1 = []
    for m in range(6):
        set_read1 = real_fp.readline()
        set_read1 = " "+set_read1
        set_read1 = set_read1.rstrip()
        set_read1 = re.sub("\n", "", set_read1)
        set_read1 = set_read1.split()
        set_list1.append(set_read1)
    #print(set_list1)
    atom_number = abs(int(set_list1[2][0]))
    x_meshsize = int(set_list1[3][0])
    y_meshsize = int(set_list1[4][0])
    z_meshsize = int(set_list1[5][0])
    num_cube_lines = (int(z_meshsize / 6) + 1) * y_meshsize * x_meshsize
    real_fp.seek(0)

    out_fp = open(out_path, "w")

    for j in range(0, atom_number + 7):
        line = real_fp.readline()
        line_dummy = imag_fp.readline()
        out_fp.write(line)

    for i in range(num_cube_lines):
        line_real = real_fp.readline()
        line_imag = imag_fp.readline()
        if line_real != "" and line_imag != "":
            line_real = str2list(line_real)
            line_imag = str2list(line_imag)
            for j in range(len(line_real)):
                charge = float(line_real[j]) ** 2 + float(line_imag[j]) ** 2
                charge = math.sqrt(charge)
                out_fp.write("  " + str(charge))
        out_fp.write("\n")

    real_fp.close()
    imag_fp.close()
    out_fp.close()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python charge_makes_cube.py <real> <imege> <output>")
        sys.exit(0)
    real_path = sys.argv[1]
    imag_path = sys.argv[2]
    out_path = sys.argv[3]
    Charge_make(real_path, imag_path, out_path)
            num_cube_lines = (int(z_meshsize / 6) + 1) * y_meshsize * x_meshsize
            real_fp.seek(0)

            with open(out_path, "w") as out_fp :
                for j in range(0, atom_number + 7):
                    line = real_fp.readline()
                    line_dummy = imag_fp.readline()
                    out_fp.write(line)
                for i in range(num_cube_lines):
                    line_real = real_fp.readline()
                    line_imag = imag_fp.readline()
                    if line_real != "" and line_imag != "":
                        line_real = str2list(line_real)
                        line_imag = str2list(line_imag)
                        for j in range(len(line_real)):
                            charge = float(line_real[j]) ** 2 + float(line_imag[j]) ** 2
                            charge = math.sqrt(charge)
                            out_fp.write("  " + str(charge))
                        if i != num_cube_lines:
                            out_fp.write("\n")

    del set_list1,atom_number,x_meshsize,z_meshsize,num_cube_lines,line,line_dummy,line_real,line_imag,charge

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python charge_makes_cube.py <dir> [parallel]")
        sys.exit(0)
    out_path = sys.argv[1]
    file_list = os.listdir(out_path + "/real/")
    out_paths = []
    for files in file_list:
        if files.find(".cube")>=0:
            out_paths.append(out_path + "/real/" + files)
    sortKey = lambda f: f if not f.startswith('.') else f[1:]
    out_paths.sort(key=sortKey)

    if len(sys.argv) == 3:
        p = Pool(16)
        p.map(charge, out_paths)
    else :
        for i in out_paths:
            charge(i)


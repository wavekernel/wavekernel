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

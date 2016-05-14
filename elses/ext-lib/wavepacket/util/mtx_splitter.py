# coding: utf-8
# Author: Yukiya Abe
import re, sys

def insert_num(filename, n):
    return "%s_%06d%s" % (filename[: -4], n, filename[-4 :])

def insert_str(filename, s):
    return "%s_%s%s" % (filename[: -4], s, filename[-4 :])

def mtx_split_end(filename,end):
    counter = 1
    with open(filename,"r") as f:
        while True:
            matinfo = f.readline()
            if not matinfo: break
            comment = f.readline()
            atom_num = f.readline()
            atom_num_int = str2list(atom_num)
            comment_int = str2list(comment)
            step = int(comment_int[-1])
            if step > end :
                break
            with open(insert_num(filename, counter), "w") as g:
                g.write(matinfo)
                g.write(comment)
                g.write(atom_num)
                atom_num = int(atom_num_int[-1])
                for i in range(atom_num):
                    line = f.readline()
                    g.write(line)
            counter += 1

def mtx_split_normal(filename):
    with open(filename,"r") as f:
        while True:
            matinfo = f.readline()
            if not matinfo: break
            comment = f.readline()
            atom_num = f.readline()
            atom_num_int = str2list(atom_num)
            comment_int = str2list(comment)
            with open(insert_num(filename, int(comment_int[-1]) + 1), "w") as g:
                g.write(matinfo)
                g.write(comment)
                g.write(atom_num)
                atom_num = int(atom_num_int[-1])
                for i in range(atom_num):
                    line = f.readline()
                    g.write(line)

def mtx_split_last(filename):
    with open(filename,"r") as f:
        while True:
            matinfo = f.readline()
            if not matinfo: break
            comment = f.readline()
            atom_num = f.readline()
            atom_num_int = str2list(atom_num)
            comment_int = str2list(comment)
            filelist = []
            filelist += [matinfo]
            filelist += [comment]
            filelist += [atom_num]
            atom_num = int(atom_num_int[-1])
            for i in range(atom_num):
                line = f.readline()
                filelist += [line]
        with open(insert_num(filename, int(comment_int[-1]) + 1), "w") as g:
            for line in filelist:
                g.write(line)

def mtx_split_step(filename,step):
    with open(filename,"r") as f:
        while True:
            matinfo = f.readline()
            if not matinfo: break
            comment = f.readline()
            atom_num = f.readline()
            atom_num_int = str2list(atom_num)
            comment_int = str2list(comment)
            filelist = []
            filelist += [matinfo]
            filelist += [comment]
            filelist += [atom_num]
            atom_num = int(atom_num_int[-1])
            for i in range(atom_num):
                line = f.readline()
                filelist += [line]
            if int(comment_int[-1]) == step:
                break
        with open(insert_num(filename, int(comment_int[-1]) + 1), "w") as g:
            for line in filelist:
                g.write(line)

def mtx_split_stride(filename,stride):
    with open(filename,"r") as f:
        with open(insert_str(filename, "stride"), "w") as g:
            while True:
                matinfo = f.readline()
                if not matinfo: break
                comment = f.readline()
                atom_num = f.readline()
                atom_num_int = str2list(atom_num)
                comment_int = str2list(comment)
                filelist = []
                filelist += [matinfo]
                filelist += [comment]
                filelist += [atom_num]
                atom_num = int(atom_num_int[-1])
                for i in range(atom_num):
                    line = f.readline()
                    filelist += [line]
                if int(comment_int[2]) % stride == 0:
                    for line in filelist:
                        g.write(line)
                        #break

def mtx_split_stride_first_end(filename,stride,first,end):
    with open(filename,"r") as f:
        with open(insert_str(filename, "stride_%06d-%06d" % (first, end)), "w") as g:
            while True:
                matinfo = f.readline()
                if not matinfo: break
                comment = f.readline()
                atom_num = f.readline()
                atom_num_int = str2list(atom_num)
                comment_int = str2list(comment)
                filelist = []
                filelist += [matinfo]
                filelist += [comment]
                filelist += [atom_num]
                atom_num = int(atom_num_int[-1])
                for i in range(atom_num):
                    line = f.readline()
                    filelist += [line]
                if int(comment_int[2]) % stride == 0 and int(comment_int[2]) >= first :
                    for line in filelist:
                        g.write(line)
                if int(comment_int[2]) > end:
                    break

def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    string = re.split(" ",string)
    return string


def main(argv):
    argc = len(argv)
    filename = argv[1]

    if argc >= 3:
        mode = argv[2]
    else :
        mode = "last"

    if mode == "normal":
        mtx_split_normal(filename)
    elif mode == "last":
        mtx_split_last(filename)
    elif mode == "step":
        if argc == 4:
            step = int(argv[3]);
        else:
            step = 0;
        mtx_split_step(filename,step)
    elif mode == "stride":
        if argc == 4 :
            stride = int(argv[3]);
        else:
            stride = 10;
        mtx_split_stride(filename,stride)
    elif mode == "stride-first":
        if argc == 5 :
            stride = int(argv[3])
            first = int(argv[4])
        else :
            stride = 10
            first = 0
        mtx_split_stride_first(filename,stride,first)
    elif mode == "stride-first-end":
        if argc == 6:
            stride = int(argv[3])
            first = int(argv[4])
            end = int(argv[5])
        else:
            stride = 10
            first = 0
            end = 1
        mtx_split_stride_first_end(filename,stride,first,end)
    elif mode == "split-end":
        if argc == 4:
            end = int(argv[3])
        mtx_split_end(filename, end)
    else:
        print "unknown mode"

usage = """Usage: python mtx_splitter_abe.py <mtx file> [<mode>]
normal
last
step <step>
stride <stride>
stride-first <stride> <first>
stride-first-end <stride> <first> <end>
split-end <end>
"""

if __name__=="__main__":
    if len(sys.argv) < 2:
        print usage
    else:
        main(sys.argv)

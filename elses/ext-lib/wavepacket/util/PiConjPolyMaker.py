# coding: utf-8
import numpy as np
import math
import sys

def get_max_group_size(group_id):
    num_groups = max(group_id) + 1
    group_sizes = [0] * num_groups
    for i in group_id:
        group_sizes[i] += 1
    return max(group_sizes)

def main_GroupID_on(num_monomers, torsion_angle_mean, torsion_angle_variance, filenames):
    pos, unit_cell, atom_list = template("para")
    rot = 0.0
    temp = "  %+10.6f  %+10.6f  %+10.6f"
    position = [0.0, 0.0, -1.0]
    group_id = []
    with open(filenames[0], "w") as fp_xyz:
        fp_xyz.write(str(num_monomers * 12) + "\n#\n")
        fp_xyz.write(("H" + temp + "\n") % (position[0], position[1], position[2]))
        group_id.append(0)
        for i in range(num_monomers):
            if i == num_monomers - 1:
                for j in range(len(atom_list) - 2):
                    position = pos[j]
                    position = rotate(position, rot, "z")
                    fp_xyz.write((atom_list[j] + temp + "\n") %
                                 (position[0] + unit_cell[0] * i,
                                  position[1] + unit_cell[1] * i,
                                  position[2] + unit_cell[2] * i))
                    group_id.append(i)
                fp_xyz.write(("H" + temp + "\n") %
                             (position[0] + unit_cell[0] * i,
                              position[1] + unit_cell[1] * i,
                              position[2] + unit_cell[2] * i + 1.0))
                group_id.append(i)
            else:
                for j in range(len(atom_list)):
                    position = pos[j]
                    position = rotate(position, rot, "z")
                    fp_xyz.write((atom_list[j] + temp + "\n") %
                                 (position[0] + unit_cell[0] * i,
                                  position[1] + unit_cell[1] * i,
                                  position[2] + unit_cell[2] * i))
                    group_id.append(i)
            if torsion_angle_variance == 0.0:
                angle = torsion_angle_mean
            else:
                angle = np.random.normal(torsion_angle_mean, torsion_angle_variance)
            rot += angle
    num_atoms = len(group_id)
    num_groups = max(group_id) + 1
    max_group_size = get_max_group_size(group_id)

    with open(filenames[1], "w") as fp_group_id:
        fp_group_id.write(
            "# number of atoms, number of groups, maximum number of group atoms\n")
        fp_group_id.write("%d %d %d\n" % (num_atoms, num_groups, max_group_size))
        for i in range(num_atoms):
            fp_group_id.write("%d %d\n" % (i + 1, group_id[i] + 1))

def main_GroupID_off(num_monomers, torsion_angle_mean):
    pos,unit_cell,atom_list = template("para")
    argvs = sys.argv
    argc  = len(argvs)
    if argc == 2:
        torsion_angle_mean = float(argvs[1])
    else:
        torsion_angle_mean = 0.0
    rot = 0.0
    temp = "  %+10.6f  %+10.6f  %+10.6f"
    position=[0.0, 0.0, -1.0]
    print "H",temp%(position[0],position[1],position[2])
    for i in range(N):
        if i == num_monomers - 1:
            for j in range(len(atom_list)-2):
                position = pos[j]
                position = rotate(position,rot,"z")
                print atom_list[j],temp%(position[0]+unit_cell[0]*i,
                                         position[1]+unit_cell[1]*i,
                                         position[2]+unit_cell[2]*i)
            print "H",temp%(position[0]+unit_cell[0]*i+1.0,
                            position[1]+unit_cell[1]*i,
                            position[2]+unit_cell[2]*i)
        else:
            for j in range(len(atom_list)):
                position = pos[j]
                position = rotate(position,rot,"z")
                print atom_list[j],temp%(position[0]+unit_cell[0]*i,
                                         position[1]+unit_cell[1]*i,
                                         position[2]+unit_cell[2]*i)
        rot += torsion_angle_mean
    print "  "
    print "  "
    print "  "

def gaussian_header_Opt():
    print "#p Opt, B3LYP/STO-3G formcheck pop=full punch=archive scf=Tight nosymm"
    print " "
    print "comment"
    print " "
    print "0 1"

def gaussian_header_Sp():
    print "#p Sp, B3LYP/STO-3G formcheck pop=full punch=archive scf=Tight nosymm"
    print " "
    print "comment"
    print " "
    print "0 1"

def rotate(position, theta, axis):
    theta = theta*math.pi/180.0
    if axis == "x":
        rot = [[1.0, 0.0, 0.0],
               [0.0, np.cos(theta), -np.sin(theta)],
               [0.0, np.sin(theta), np.cos(theta)]]
        position = np.dot(rot, position)
        return position
    elif axis == "y":
        rot = [[np.cos(theta), 0.0, -np.sin(theta)],
               [0.0          , 1.0,  0.0          ],
               [np.sin(theta), 0.0,  np.cos(theta)]]
        position = np.dot(rot,position)
        return position
    elif axis == "z":
        rot = [[np.cos(theta), -np.sin(theta), 0.0],
               [np.sin(theta),  np.cos(theta), 0.0],
               [0.0          ,  0.0          , 1.0]]
        position = np.dot(rot,position)
        return position

def template(mol):
    if mol == "para":
        global a,b,cc1,cc2,cc3
        # For Gaussian ***************************
        #a   = 1.42494
        #b   = 1.09828
        #cc1 = a
        #cc2 = 1.45247
        #cc3 = 1.21536

        # For ELSES ******************************
        a   = 1.42514
        b   = 1.06300
        cc1 = a
        cc2 = 1.46400
        cc3 = 1.18559

        unit_cell = [0.0, 0.0, 4*a*np.cos(60*math.pi/180.0)+2*cc2+cc3]

        s = np.sin(math.pi / 3.0)
        c = np.cos(math.pi / 3.0)
        sa = s * a
        sb = s * b
        ca = c * a
        cb = c * b
        pos = [[  sa + sb, 0.0, ca - cb], # H
               [- sa - sb, 0.0, ca - cb], # H
               [  sa + sb, 0.0, 3.0 * ca + cb], # H
               [- sa - sb, 0.0, 3.0 * ca + cb], # H
               [0.0, 0.0, 0.0], # C
               [  sa, 0.0, ca], # C
               [- sa, 0.0, ca], # C
               [  sa, 0.0, 3.0 * ca], # C
               [- sa, 0.0, 3.0 * ca], # C
               [0.0, 0.0, 4.0 * ca], # C
               [0.0, 0.0, 4.0 * ca + cc2], # C
               [0.0, 0.0, 4.0 * ca + cc2 + cc3]] # C

        atom_list = ["H","H","H","H","C","C","C","C","C","C","C","C"]
    else:
        pos = [0.0,0.0,0.0]
        unit_cell = [1.0 ,1.0 ,1.0]
        atom_list = ["H"]

    return pos,unit_cell,atom_list

def get_filenames(num_monomers, torsion_angle_mean, torsion_angle_variance, is_straight_mode):
    if is_straight_mode:
        str_nums = "para_nm%d" % num_monomers
    else:
        str_nums = ("para_nm%d_ta%d_tav%d") % (num_monomers, int(torsion_angle_mean), int(torsion_angle_variance))
    return (str_nums + ".xyz", str_nums + "_group_id.txt")

if __name__=="__main__":
    if len(sys.argv) <= 1:
        print "Usage: python PiConjPolyMaker.py <num_monomers> [<torsion_angle_mean>] [<torsion_angle_variance>]"
        sys.exit(0)

    num_monomers = int(sys.argv[1])
    is_straight_mode = len(sys.argv) == 2

    if len(sys.argv) >= 3:
        torsion_angle_mean = float(sys.argv[2])
    else:
        torsion_angle_mean = 0.0
    if len(sys.argv) >= 4:
        torsion_angle_variance = float(sys.argv[3])
    else:
        torsion_angle_variance = 0.0

    filenames = get_filenames(num_monomers, torsion_angle_mean, torsion_angle_variance, is_straight_mode)
    print "filenames:", filenames[0], filenames[1]

    # For ELSES_setting!!################
    # Choose (1) or (2)
    main_GroupID_on(num_monomers, torsion_angle_mean, torsion_angle_variance, filenames)  # (1)
    #main_GroupID_off(num_monomers, torsion_angle_mean) # (2)

    # For Gaussian_setting!!#############
    # Choose (1) or (2)
    #gaussian_header_Opt() # (1)
    #gaussian_header_Sp()  # (2)
    # For Gaussian Input File create.
    #main_GroupID_off(num_monomers, torsion_angle_mean)

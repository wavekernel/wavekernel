# coding: utf-8
import argparse,re,sys,os,os.path

def read_group_id(group_id_name):
    group_id = []
    with open(group_id_name,"r") as f:
        line = f.readline()
        line = f.readline()
        a = str2list(line)
        for i in range(int(a[0])):
            line = f.readline()
            n = str2list(line)
            group_id += [int(n[1])]
    return group_id
def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)
def read_xyz(filename):
    xyz = {}
    x = []
    y = []
    z = []
    fp = open(filename,"r")
    atom_num = fp.readline() # atom num
    atom_num = int(str2list(atom_num)[0])
    line = fp.readline() # comment
    counter = 0
    for line in fp:
        if atom_num <= counter :
            break
        line = str2list(line)
        x += [float(line[1].replace("D","E"))*1.889725988579]
        y += [float(line[2].replace("D","E"))*1.889725988579]
        z += [float(line[3].replace("D","E"))*1.889725988579]
        counter += 1
    print(counter ,"atom")
    xyz["x"] = x
    xyz["y"] = y
    xyz["z"] = z
    return xyz
def main(xyz_path,group_path,group_min,group_max,mesh,cutoff):
    xyz = read_xyz(xyz_path)
    group_id = list(range(len(xyz["x"])))
    if group_path != False :
        group_id = read_group_id(group_path)
    x = []
    y = []
    z = []
    #print group_id
    for index, item in enumerate(group_id):
        if group_min <= item and item <= group_max:
            x += [xyz["x"][index]]
            y += [xyz["y"][index]]
            z += [xyz["z"][index]]
    x_min = min(x)
    x_max = max(x)
    y_min = min(y)
    y_max = max(y)
    z_min = min(z)
    z_max = max(z)
    a = cutoff
    with open("input_mesh_grid.txt","w") as f:
        f.write('%d %d %d\n'%(int((x_max - x_min+a*2)*mesh)+10,int((y_max - y_min+a*2)*mesh)+10,int((z_max - z_min+a*2)*mesh)+10))
        f.write('%f %f %f au\n'%(x_min-a,y_min-a,z_min-a))
        f.write('%f %f %f au\n'%(x_max - x_min+a*2,y_max - y_min+a*2,z_max - z_min+a*2))
    print((int(x_max - x_min+a*2)*mesh+10)*(int(y_max - y_min+a*2)*mesh+10)*(int(z_max - z_min+a*2)*mesh+10)*1.0e-6*24.7 ,"MB")
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('xyz_path', metavar='XYZ', type=str,
                        help='')
    parser.add_argument('-group', metavar='GROUP', dest='group_path', type=str, default=False,
                        help='')
    parser.add_argument('-group-min', metavar='GROUP-MIN', dest='group_min', type=int, default=1,
                        help='')
    parser.add_argument('-group-max', metavar='GROUP-MAX', dest='group_max', type=int, default=10000000,
                        help='')
    parser.add_argument('-cutoff-au', metavar='CUTOFF', dest='cutoff', type=float, default=8.0,
                        help='')
    parser.add_argument('-mesh', metavar='MESH', dest='mesh', type=float, default=1.0,
                        help='')
    args = parser.parse_args()
    
    main(args.xyz_path,args.group_path,args.group_min,args.group_max,args.mesh,args.cutoff)

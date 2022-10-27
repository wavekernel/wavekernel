import argparse,json, datetime, struct
import os, re, shutil, sys ,math
import subprocess
from multiprocessing import Pool, cpu_count, current_process
import time
try:
    from mpi4py import MPI
    MPI_flag = True
except:
    MPI_flag = False
    rank = 0

import charge_makes_cube
import LCAO_converter
import input_mesh_grid_lib

def str2list(string):
    string = string.strip()
    string = re.sub("  *"," ",string)
    return re.split(" ",string)

def mkdir2(outdir):
    if os.path.isdir(outdir):
        sys.stderr.write("output directory " + outdir + " already exists\n")
    else :
        os.mkdir(outdir)
def check_file_exists(path):
    if not os.path.isfile(path):
        sys.stderr.write("file " + path + " does not exist\n")
        sys.exit(1)
def charge(cube_path):
    if cube_path.find("psi_real") >= 0:
        cube_imag_path = re.sub("real", "imag", cube_path)
        cube_char_path = re.sub("real", "char", cube_path)
        charge_makes_cube.Charge_make(cube_path, cube_imag_path, cube_char_path)
        del cube_imag_path,cube_char_path
def target_get(target):
    if target != False:
        target_str = target.split(",")
        target = []
        for i in target_str:
            target += [int(i)]
    else:
        target = False
    #del target_str
    return target
def get_OMP():
    try:
        Parallel_core = int(os.environ.get("OMP_NUM_THREADS"))
    except:
        Parallel_core = cpu_count()
    return Parallel_core
def get_periodic(periodic,xyz_file_name):
    periodic_mode = {}
    periodic_mode["x_mode"] = False
    periodic_mode["y_mode"] = False
    periodic_mode["z_mode"] = False
    periodic_mode["flag"] = False
    if len(periodic) >= 1:
        with open(xyz_file_name) as f:
            line = f.readline()
        line = str2list(line)
        if len(line) < 4:
            sys.stderr.write("not found periodic infomation\n")
            sys.exit(1)
        else :
            periodic_mode["x"] = float(line[1])*1.889725988579
            periodic_mode["y"] = float(line[2])*1.889725988579
            periodic_mode["z"] = float(line[3])*1.889725988579
            periodic_mode["x_t"] = []
            periodic_mode["y_t"] = []
            periodic_mode["z_t"] = []
    if len(periodic) == 1:
        periodic_mode["flag"] = 1
        periodic_mode[periodic+"_mode"] = True
    elif len(periodic) == 2:
        periodic_mode["flag"] = 2
        periodic_mode[periodic[0]+"_mode"] = True
        periodic_mode[periodic[1]+"_mode"] = True
    elif len(periodic) == 3:
        periodic_mode["flag"] = 3
        periodic_mode[periodic[0]+"_mode"] = True
        periodic_mode[periodic[1]+"_mode"] = True
        periodic_mode[periodic[2]+"_mode"] = True
    return periodic_mode
def make_cube(arg):
    out_path , elses_generate_cube_path ,cutoff= arg
    dirname = os.path.dirname(out_path)
    basename = os.path.basename(out_path)
    shutil.copy(out_path, os.getcwd())
    command = [elses_generate_cube_path,
               "-lcao_coef_file=%s" % (basename),
                "-cube_filename_header=%s_" % (basename[:17]),
                "-cutoff_au=%f" % (cutoff)]
    subp = subprocess.Popen(" ".join(command), shell=True)
    os.waitpid(subp.pid, 0)
    os.remove(basename)
    match = re.search("(.+).txt$", basename)
    cube_path = os.path.join(dirname, match.group(1) + ".cube")
    #cube_paths.append(cube_path)
    # tmp_000001.cube is a default output name of elses-generate-cubefile
    shutil.move("%s_000001.cube"%(basename[:17]), cube_path)
    return cube_path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('elses_generate_cube_path',
                        metavar='elses-generate-cubefile',
                        type=str,help='')
    parser.add_argument('wave_out_path', metavar='JSON', type=str,
                        help='')
    parser.add_argument('position_path', metavar='XYZ', type=str,
                        help='')
    parser.add_argument('basis_path', metavar='BASIS', type=str,
                        help='')
    parser.add_argument('-s', metavar='STRIDE', dest='skip_stride_num',
                        type=int, default=1,help='')
    parser.add_argument('-load-min', metavar='LOAD-MIN',
                        dest='load_min_num', type=int, default=0,help='')
    parser.add_argument('-load-max', metavar='LOAD-MAX',
                        dest='load_max_num', type=int,
                        default=10000000,help='')
    parser.add_argument('-cutoff-au', metavar='CUTOFF', dest='cutoff',
                        type=float, default=8,help='')
    parser.add_argument('-input_mesh_grid', metavar='PATH',
                        dest='input_mesh_grid_path', type=str,
                        default=False,help='')
    parser.add_argument('-alpha', metavar='ALPHA', dest='alpha',
                        type=float, default=False,help='')
    parser.add_argument('-mesh', metavar='MESH', dest='mesh',
                        type=float, default=1.0,help='')
    parser.add_argument('-target', metavar='target', dest='target',
                        type=str, default=False,help='ex) 12,24,15')
    parser.add_argument('-core', metavar='CORE', dest='core_num',
                        default=get_OMP(),type=int,help='')
    parser.add_argument('-periodic', metavar='periodic', dest='periodic',
                        default="",type=str,help='ex) x or zy or xyz')
    parser.add_argument('--parallel', action='store_true',
                        dest='parallel_flag',default=False, help='')
    parser.add_argument('--position', action='store_true',
                        dest='position_flag',default=False,help='')
    parser.add_argument('--LCAO', action='store_true',
                        dest='lcao_flag',default=False,help='')
    parser.add_argument('--big-endian', action='store_false',
                        dest='is_little_endian',default=True, help='')
    args = parser.parse_args()

    if MPI_flag == True:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    if rank == 0 :
        print("MPI_flag",MPI_flag)
        if MPI_flag :
            print("node = ",size)
        start_time = time.time()
        # set parameter 
        target = target_get(args.target)
        print("core = ",args.core_num)
        elses_generate_cube_path = args.elses_generate_cube_path
        wave_out_path = args.wave_out_path
        is_little_endian = args.is_little_endian
        stride = args.skip_stride_num
        max_num = args.load_max_num
        min_num = args.load_min_num

        # check input file
        check_file_exists(elses_generate_cube_path)
        check_file_exists(wave_out_path)
        check_file_exists(args.position_path)
        to_use_input_mesh_grid = False
        if args.input_mesh_grid_path != False and\
           args.position_flag == False and args.alpha == False :
            to_use_input_mesh_grid = True
            input_mesh_grid_path = args.input_mesh_grid_path
            check_file_exists(input_mesh_grid_path)

        # read xyz file
        periodic = get_periodic(args.periodic,args.position_path)
        xyz_load_time_start = time.time()
        position_data = input_mesh_grid_lib.main(
            args.cutoff,args.position_path,args.mesh,periodic)
        xyz_load_time_end = time.time()

        # make output dir
        outdir = wave_out_path + ".cube"
        mkdir2(outdir)

        # read json file
        json_load_time_start = time.time()
        with open(wave_out_path, "r") as fp:
            wave_out = json.load(fp)
        # split json
        if wave_out["setting"]["is_output_split"]:
            split_dir = os.path.dirname(wave_out_path)
            wave_out_cropped = LCAO_converter.crop_from_splits(
                wave_out,stride,max_num,min_num,args.parallel_flag,
                target,args.core_num,is_little_endian,split_dir)
        else:
            wave_out_cropped = LCAO_converter.crop(
                wave_out, stride,is_little_endian,max_num,min_num,target)
        json_load_time_end = time.time()

        # LCAO convert
        make_LCAO_time_start = time.time()
        out_paths = LCAO_converter.LCAO_converter(
            wave_out_cropped, position_data, outdir,periodic,args.basis_path)
        make_LCAO_time_end = time.time()
        del wave_out_cropped

        # make input_mesh_grid (periodic xyz)
        if periodic["flag"] != False:
            if periodic["x_mode"]:
                x_min , x_max = 0.0,periodic["x"]
            else:
                x_min,x_max = position_data["input_mesh_grid"]["x_min"],position_data["input_mesh_grid"]["x_max"]
            if periodic["y_mode"]:
                y_min , y_max = 0.0,periodic["y"]
            else:
                y_min,y_max = position_data["input_mesh_grid"]["y_min"],position_data["input_mesh_grid"]["y_max"]
            if periodic["z_mode"]:
                z_min , z_max = 0.0,periodic["z"]
            else:
                z_min,z_max = position_data["input_mesh_grid"]["z_min"],position_data["input_mesh_grid"]["z_max"]
            #input_mesh_grid_lib.print_input_mesh_grid(0.0,periodic["x"],0.0,periodic["y"],0.0,periodic["z"],args.mesh,0)
            input_mesh_grid_lib.print_input_mesh_grid(x_min,x_max,y_min,y_max,z_min,z_max,args.mesh,0)
        # make input_mesh_grid (position file)
        if args.position_flag == True:
            input_mesh_grid_lib.print_input_mesh_grid(
                position_data["input_mesh_grid"]["x_min"],
                position_data["input_mesh_grid"]["x_max"],
                position_data["input_mesh_grid"]["y_min"],
                position_data["input_mesh_grid"]["y_max"],
                position_data["input_mesh_grid"]["z_min"],
                position_data["input_mesh_grid"]["z_max"],
                args.mesh,args.cutoff)

        # make alpha input_mesh_grid
        if args.position_flag != False and args.alpha != False :
            input_mesh_grid_lib.make_position_alpha(
                wave_out_cropped,position_data,args.alpha)
        elif to_use_input_mesh_grid == False and args.alpha != False:
            input_mesh_grid_lib.make_alpha(
                wave_out_cropped,args.alpha,args.mesh,args.cutoff)

        # move input_mesh_grid
        if to_use_input_mesh_grid:
            if os.path.dirname(input_mesh_grid_path) != os.getcwd():
                if os.path.isfile(os.getcwd()+"/input_mesh_grid.txt"):
                    sys.stderr.write("input_mesh_grid is already exists\n")
                else:
                    shutil.copy(input_mesh_grid_path, os.getcwd())
        # Remove the copied input_mesh_grid.txt
        #if to_use_input_mesh_grid:
            #if os.path.dirname(input_mesh_grid_path) != os.getcwd():
                #os.remove(os.path.basename(input_mesh_grid_path))

    if args.lcao_flag : # only LCAO mode 
        if rank == 0:
            print("only LCAO mode complate")
        sys.exit()

    if MPI_flag : 
        if rank == 0 :
            args_list = []
            for i in range(size):
                args_mpi = []
                for j in range(i,len(out_paths),size):
                    args_mpi += [ (out_paths[j],
                                   elses_generate_cube_path ,args.cutoff) ]
                args_list += [args_mpi]
        else:
            args_list = []
        comm.Barrier()
        args_list = comm.scatter(args_list, root=0)
    else :
        args_list = []
        for out_path in out_paths:
            args_list += [ (out_path , elses_generate_cube_path ,args.cutoff )]

    # used elses-generate-cubefile
    if rank ==0 :
        egc_time_start = time.time()
    cube_paths = []
    for arg in args_list:
        cube_paths += [make_cube(arg)]# elses-generate-cubefile
    if MPI_flag :
        comm.Barrier()
    if rank ==0 :
        egc_time_end = time.time()

    # Make cubefiles of charge from psi_real_* and psi_imag_*
    if rank == 0:
        make_charge_time_start = time.time()
        print("char start")
        mkdir2(os.path.join(outdir, "char"))

    # It is supposed that paths of cubefiles for real part are form of
    # "hoge/real/psi_real_fuga" and these for imaginary part are form of
    # "hoge/imag/psi_imag_fuga".
    if args.parallel_flag :
        p = Pool(args.core_num)
        p.map(charge, cube_paths)
    else :
        for i in cube_paths:
            charge(i)
    if MPI_flag:
        comm.Barrier()
    if rank == 0:
        make_charge_time_end = time.time()
        end_time = time.time()
        print("WaveTool complete")
        if MPI_flag:
            file_name = "time_MPI_node_%d_core_%d.txt"%(size,args.core_num)
        else :
            file_name = "time_OMP_core_%d.txt"%(args.core_num)
        with open(file_name,"w") as f:
            f.write("all time " + str(end_time - start_time) + " \n"  )
            f.write("xyz load time " +
                    str(xyz_load_time_end - xyz_load_time_start) + " \n" )
            f.write("json load time " +
                    str(json_load_time_end - json_load_time_start) + " \n" )
            f.write("make LCAO time " +
                    str(make_LCAO_time_end - make_LCAO_time_start) + " \n" )
            f.write("elses-generate-cubefile time " +
                    str(egc_time_end - egc_time_start) + " \n")
            f.write("make charge time " +
                    str(make_charge_time_end - make_charge_time_start) + " \n")

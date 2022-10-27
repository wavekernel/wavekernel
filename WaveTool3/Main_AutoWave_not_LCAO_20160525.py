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

def get_OMP():
    try:
        Parallel_core = int(os.environ.get("OMP_NUM_THREADS"))
    except:
        Parallel_core = cpu_count()
    return Parallel_core

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
                        metavar='<elses-generate-cubefile>',
                        type=str,help='')
    parser.add_argument('wave_out_path', metavar='<JSON.CUBE>', type=str,
                        help='')
    parser.add_argument('-core', metavar='CORE', dest='core_num',
                        default=get_OMP(),type=int,help='')
    parser.add_argument('--parallel', action='store_true',
                        dest='parallel_flag',default=False, help='')
    parser.add_argument('-cutoff-au', metavar='CUTOFF', dest='cutoff',
                        type=float, default=8,help='')
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
        print("core = ",args.core_num)
        elses_generate_cube_path = args.elses_generate_cube_path
        wave_out_path = args.wave_out_path

        # check input file
        #check_file_exists(elses_generate_cube_path)
        #check_file_exists(wave_out_path)

    ####################
    # get output_wavefunction
    if rank == 0:
        out_paths = []
        # real part
        file_list = os.listdir(wave_out_path + "/real/")
        for f in file_list:
            if f.find(".txt")>=0:
                out_paths.append(wave_out_path + "/real/" + f)
        # imag part
        file_list = os.listdir(wave_out_path + "/imag/")
        for f in file_list:
            if f.find(".txt")>=0:
                out_paths.append(wave_out_path + "/imag/" + f)        
    else :
        out_paths = []

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
        mkdir2(os.path.join(wave_out_path, "char"))

    # It is supposed that paths of cubefiles for real part are form of
    # "hoge/real/psi_real_fuga" and these for imaginary part are form of
    # "hoge/imag/psi_imag_fuga".
    if args.parallel_flag :
        p = Pool(args.core_num)
        p.map(charge, cube_paths)
    else :
        for i in cube_paths:
            charge(i)
    """
    if MPI_flag:
        comm.Barrier()
    if rank == 0:
        make_charge_time_end = time.time()
        end_time = time.time()
        print "WaveTool complete"
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
    """
    

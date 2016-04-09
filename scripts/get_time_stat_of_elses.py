import re, sys, os.path

def avg(xs):
    return sum(xs) / len(xs)

if __name__ == '__main__':
    dirname = sys.argv[1]
    assert(os.path.isdir(dirname))
    filenames = os.listdir(dirname)
    log_node_paths = []
    i = 0
    max_i = 100000000
    for f in filenames:
        if i >= max_i:
            break
        if f.find('log-node') >= 0:
            log_node_paths.append(os.path.join(dirname, f))
            i += 1
    befor_times = []
    body_total_times = []
    body_mpi_times = []
    body_wait_times = []  

    for i, f in enumerate(log_node_paths):
        body_total_times_in_node = []  # [node]
        body_mpi_times_in_node = []   
        body_wait_times_in_node = []  
        if i % 500 == 0:
            print 'read file %d / %d: %s' % (i + 1, len(log_node_paths), f)
        with open(f) as fp:
            is_reading_num_nodes = True
            for line in fp:
                if is_reading_num_nodes:
                    m = re.search(r'INFO-MPI-OMP: P_MPI, P_OMP=\s*(\d+)\s+\d+', line)
                    if m:
                        num_nodes = int(m.group(1))
                        is_reading_num_nodes = False
                else:
                    m = re.search(r'elaps-time\((lap|befor)\s+MDloop\)=\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
                    if m:
                        mode = m.group(1)
                        #step = int(m.group(2))
                        time_total = float(m.group(3))
                        time_mpi = float(m.group(4))
                        time_wait = float(m.group(5))
                        if mode == 'lap':
                            body_total_times_in_node.append(time_total)
                            body_mpi_times_in_node.append(time_mpi)
                            body_wait_times_in_node.append(time_wait)
                        else:  # mode == 'befor'
                            befor_times.append(time_total)
        for lap, (tt, tm, tw) in enumerate(zip(body_total_times_in_node,
                                               body_mpi_times_in_node,
                                               body_wait_times_in_node)):
            if len(body_total_times) <= lap:
                body_total_times.append([])
            if len(body_mpi_times) <= lap:
                body_mpi_times.append([])
            if len(body_wait_times) <= lap:
                body_wait_times.append([])
            body_total_times[lap].append(tt)
            body_mpi_times[lap].append(tm)
            body_wait_times[lap].append(tw)

    # befor_times[node]        
    # body_total_times[lap][node]
    # body_mpi_times[lap][node]
    # body_wait_times[lap][node]

    befor_total_node0 = befor_times[0]
    befor_total_node1 = befor_times[1]
    befor_total_nodeavg = avg(befor_times)

    total_node0 = map(lambda t: t[0], body_total_times)  # [lap]
    total_node1 = map(lambda t: t[1], body_total_times)  # [lap]
    total_node0_lapavg = avg(total_node0)
    total_node1_lapavg = avg(total_node1)
    total_nodeavg = map(lambda ts: sum(ts) / len(ts), body_total_times)  # [lap]
    print 'total_nodeavg (per lap): ', total_nodeavg
    total_nodeavg_lapavg = avg(total_nodeavg)

    mpi_node0 = map(lambda t: t[0], body_mpi_times)  # [lap]
    mpi_node1 = map(lambda t: t[1], body_mpi_times)  # [lap]
    mpi_node0_lapavg = avg(mpi_node0)
    mpi_node1_lapavg = avg(mpi_node1)
    mpi_nodeavg = map(avg, body_mpi_times)  # [lap]
    print 'mpi_nodeavg (per lap): ', mpi_nodeavg
    mpi_nodeavg_lapavg = avg(mpi_nodeavg)

    wait_node0 = map(lambda t: t[0], body_wait_times)  # [lap]
    wait_node1 = map(lambda t: t[1], body_wait_times)  # [lap]    
    wait_node0_lapavg = avg(wait_node0)
    wait_node1_lapavg = avg(wait_node1)
    wait_nodemax = map(max, body_wait_times)  # [lap]
    print 'wait_nodemax (per lap): ', wait_nodemax
    wait_nodemax_lapavg = avg(wait_nodemax)

    print '# num_nodes befor(0) befor(1) avg.befor avg.lap(0) avg.lap(1) avg.avg.lap avg.mpi(0) avg.mpi(1) avg.avg.mpi avg.wait(0) avg.wait(1) avg.max.wait'
    print '%d %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f' % (num_nodes,
                                                                              befor_total_node0,
                                                                              befor_total_node1,
                                                                              befor_total_nodeavg,
                                                                              total_node0_lapavg,
                                                                              total_node1_lapavg,
                                                                              total_nodeavg_lapavg,
                                                                              mpi_node0_lapavg,
                                                                              mpi_node1_lapavg,
                                                                              mpi_nodeavg_lapavg,
                                                                              wait_node0_lapavg,
                                                                              wait_node1_lapavg,
                                                                              wait_nodemax_lapavg)


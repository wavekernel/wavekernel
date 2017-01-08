import re, sys, os.path

if __name__ == '__main__':
    dirname = sys.argv[1]  # Output directory of ELSES.
    assert(os.path.isdir(dirname))
    log_master_path0 = os.path.join(dirname, 'log-node000000.txt')
    log_master_path1 = os.path.join(dirname, 'log-node000001.txt')
    data = []
    steps = None
    for f in [log_master_path0, log_master_path1]:
        lap_times_total = []
        lap_times_mpi = []
        lap_times_wait = []
        with open(f) as fp:
            for line in fp:
                m = re.search(r'elaps-time\((lap|befor)\s+MDloop\)=\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
                if m:
                    mode = m.group(1)
                    step = int(m.group(2))
                    #if step == 1 or step == 6:  # Skip when the first and last step consume long time due to output of restart.xml
                    #    continue
                    print step
                    time_total = float(m.group(3))
                    time_mpi = float(m.group(4))
                    time_wait = float(m.group(5))
                    if mode == 'lap':
                        lap_times_total.append(time_total)
                        lap_times_mpi.append(time_mpi)
                        lap_times_wait.append(time_wait)
                    else:  # mode == 'befor'
                        befor_time = time_total
        if steps is None:
            steps = len(lap_times_total)
        else:
            assert(steps == len(lap_times_total))
            assert(steps == len(lap_times_mpi))
            assert(steps == len(lap_times_wait))
        data.append((befor_time, sum(lap_times_total) / steps, sum(lap_times_mpi) / steps, sum(lap_times_wait) / steps))
    print 'steps: ', steps
    print '# befor, avg_time_total, avg_time_mpi, avg_time_wait'
    for (befor, avg_time_total, avg_time_mpi, avg_time_wait) in data:
        print '%.4f %.4f %.4f %.4f' % (befor, avg_time_total, avg_time_mpi, avg_time_wait)
    if len(data) == 2:
        print '%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f' % (data[0][0], data[1][0], data[0][1], data[1][1],
                                                           data[0][2], data[1][2], data[0][3], data[1][3])

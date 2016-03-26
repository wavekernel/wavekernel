import re, sys, os.path

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
            if i > 0 and i % 500 == 0:
                print 'file %d read' % i
            i += 1
    data = []
    for f in log_node_paths:
        lap_times_total = []
        lap_times_mpi = []
        lap_times_wait = []
        with open(f) as fp:
            for line in fp:
                m = re.search(r'elaps-time\((lap|befor)\s+MDloop\)=\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
                if m:
                    mode = m.group(1)
                    step = int(m.group(2))
                    time_total = float(m.group(3))
                    time_mpi = float(m.group(4))
                    time_wait = float(m.group(5))
                    if mode == 'lap':
                        lap_times_total.append(time_total)
                        lap_times_mpi.append(time_mpi)
                        lap_times_wait.append(time_wait)
                    else:  # mode == 'befor'
                        befor_time = time_total
        for i, t in zip(range(len(lap_times_wait)), lap_times_wait):
            if len(data) <= i:
                data.append([])
            data[i].append(t)
    times_wait_for_laps = data  #map(list, zip(*data))

    max_wait = 0.0
    for lap in range(len(times_wait_for_laps)):
        print '%d %.4f' % (lap, max(times_wait_for_laps[lap]))
        max_wait = max(max_wait, max(times_wait_for_laps[lap]))
    print 'max %.4f' % max_wait

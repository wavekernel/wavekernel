import re, sys, os.path

def get_avg(xs):
    return sum(xs) / len(xs)

if __name__ == '__main__':
    filename = sys.argv[1]
    totals = {}
    mpis = {}
    waits = {}
    befors = {}
    nodes = {}
    steps = {}
#split82944_co10.0_prj100_mtcrs_lpg32MB_j5133430/log-node010313.txt:elaps-time(lap    MDloop)=              2      32.640000       1.656767       3.519759
    with open(filename) as fp:
        for line in fp:
            m = re.search(r'log-node(\d+)\.txt:elaps-time\((lap|befor)\s+MDloop\)=\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
            if m:
                node = int(m.group(1))
                mode = m.group(2)
                step = int(m.group(3))
                time_total = float(m.group(4))
                time_mpi = float(m.group(5))
                time_wait = float(m.group(6))
                if mode == 'lap':
                    nodes[node] = True
                    steps[step] = True
                    totals[(node, step)] = time_total
                    mpis[(node, step)] = time_mpi
                    waits[(node, step)] = time_wait
                else:  # mode == 'befor'
                    befors[(node, step)] = time_total
    nodes = sorted(nodes.keys())
    num_nodes = len(nodes)
    steps = sorted(steps.keys())
    steps = filter(lambda s: s > 0, steps)
    num_steps = len(steps)
    print 'num_nodes: ', num_nodes
    print 'steps: ', steps
    
    avgs_total = []
    avgs_mpi = []
    maxs_wait = []
    maxs_wait_node = []
    mins_wait = []
    mins_wait_node = []    
    for s in steps:
        avg_total = 0.0
        avg_mpi = 0.0
        max_wait = 0.0
        max_wait_node = None
        min_wait = 100000.0
        min_wait_node = None
        for n in nodes:
            avg_total += totals[(n, s)]
            avg_mpi += mpis[(n, s)]
            w = waits[(n, s)]
            if max_wait < w:
                max_wait = w
                max_wait_node = n
            if min_wait > w:
                min_wait = w
                min_wait_node = n
        avg_total /= num_nodes
        avg_mpi /= num_nodes
        avgs_total.append(avg_total)
        avgs_mpi.append(avg_mpi)
        maxs_wait.append(max_wait)
        maxs_wait_node.append(max_wait_node)
        mins_wait.append(min_wait)
        mins_wait_node.append(min_wait_node)
    avg_avg_total = get_avg(avgs_total)
    avg_avg_mpi = get_avg(avgs_mpi)
    avg_max_wait = get_avg(maxs_wait)
    avg_min_wait = get_avg(mins_wait)

    print ''
    print 'max wait.\nlap node wait'
    for s, n, t in zip(steps, maxs_wait_node, maxs_wait):
        print '%d %5d %.4f' % (s, n, t)
    print ''
    print 'min wait.\nlap node wait'
    for s, n, t in zip(steps, mins_wait_node, mins_wait):
        print '%d %5d %.4f' % (s, n, t)

    print ''
    print 'avg.avg.total, avg.avg.mpi, avg.max.wait'
    print '%.4f %.4f %.4f' % (avg_avg_total, avg_avg_mpi, avg_max_wait)

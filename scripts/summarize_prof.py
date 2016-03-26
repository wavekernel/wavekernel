import re, sys

def line_to_data(line):
    r = r'[+-]?(\d+\.)?\d+([deE][+-]?\d+)?'
    #       Elapsed(s)       MFLOPS MFLOPS/PEAK(%)           MIPS   MIPS/PEAK(%) 
    m = re.match(r'\s+(?P<elapsed>%s)\s+(?P<mflops>%s)\s+(?P<mflops_eff>%s)\s+(?P<mips>%s)\s+(?P<mips_eff>%s)\s+Application' % 
                 (r, r, r, r, r), line)
    if m:
        elapsed = float(m.group('elapsed'))
        mflops = float(m.group('mflops'))
        mflops_eff = float(m.group('mflops_eff'))
        return (elapsed, mflops, mflops_eff)
    else:
        return None

if __name__ == '__main__':
    if len(sys.argv) == 2:
        with open(sys.argv[1]) as fp:
            for line in fp:
                data = line_to_data(line)
                if data:
                    print '# elapsed(s)    MFLOPS    MFLOPS/PEAK(%)'
                    print data[0], data[1], data[2]
                    break
    elif len(sys.argv) == 5:
        with open(sys.argv[1]) as fp:
            for line in fp:
                data1 = line_to_data(line)
                if data1:
                    break
        with open(sys.argv[2]) as fp:
            for line in fp:
                data2 = line_to_data(line)
                if data2:
                    break
        flop1 = data1[0] * data1[1]
        flop2 = data2[0] * data2[1]
        main_flop = flop1 - flop2
        main_avg_elapsed_per_step = float(sys.argv[3])
        main_num_steps = int(sys.argv[4])
        main_elapsed = main_avg_elapsed_per_step * main_num_steps
        main_flops = main_flop / main_elapsed
        peak_flops = data1[1] / (data1[2] / 100.0)
        main_eff = main_flops / peak_flops * 100.0
        print '## elapsed(s)    MFLOPS    MFLOPS/PEAK(%)'
        print '# Total'
        print data1[0], data1[1], data1[2]
        print '# Initial Only'
        print data2[0], data2[1], data2[2]
        print '# Main Loop Only. flop_main, flop_initial: ', flop1, flop2
        print main_elapsed, main_flops, main_eff
        

      

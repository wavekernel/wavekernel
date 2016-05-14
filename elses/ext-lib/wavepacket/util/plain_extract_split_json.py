import json, sys

if __name__ == '__main__':
    print '# %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('time', 'input_step', 'is_after_matrix_replace',
                                                        'tb_energy', 'nl_energy', 'total_energy',
                                                        'charge_mean', 'charge_msd', 'psi_ipratio', 'alpha_ipratio')
    for f in sys.argv[1 :]:
        with open(f) as fp:
            out = json.load(fp)
        for s in out['states']:
            t = s['time']
            i = s['input_step']
            isamr = s['is_after_matrix_replace']
            etb = s['TB_energy']
            ebl = s['NL_energy']
            etotal = s['total_energy']
            mean = s['charge_coordinate_mean'][0]  # x
            msd = s['charge_coordinate_msd'][3]  # total
            pp = s['psi']['ipratio']
            ap = s['alpha']['ipratio']
            print '%.16e\t%d\t%r\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e\t%.16e' % (t, i, isamr, etb, ebl, etotal, mean, msd, pp, ap)

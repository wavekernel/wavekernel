#!/bin/sh
d=`dirname $0`
f=$1
python ${d}/plot_pratio.py ${f}_pratio.json
python ${d}/plot_alpha.py ${f}_alpha.json
python ${d}/plot_energy.py ${f}_energy.json
python ${d}/plot_charge_msd.py ${f}_charge_moment.json
convert +append ${f}_alpha.png ${f}_pratio.png ${f}_charge_moment_msd.png ${f}_energy.png ${f}_append.png

#!/bin/sh -

#########################################
## script for evaluating difference between Gaussian and ELSES with paramset.

# include settings
if [ ! -r ../setting.sh ]; then
    echo "Error: not found setting.sh in $PWD/.."
    exit 1
else
    . ../setting.sh
fi

# include variables

$elses_readcf ../param.cfg > ../config.sh
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses_readcf"
    exit 1
fi

if [ ! -r ../config.sh ]; then
    echo "Error: not found config.sh in $PWD/.."
    exit 1
else
    . ../config.sh
fi


if [ "$spin" == "" ]; then
    spin="nospin"
elif [ "$spin" == "off" ]; then
    spin="nospin"
elif [ "$spin" == "on" ]; then
    spin="spin"
fi


# output files for this script.
output_data_optimized="data_optimized.dat"
output_data_conformation="data_conformation.dat"
output_energy_eigen_conformation="energy_eigen_conformation.dat"
output_energy_total_conformation="energy_total_conformation.dat"
output_energy_difference="energy_difference.dat"
output_MO_match="MO_match_information.dat"

# gnuplot script file to output
gnuplot_energy_eigen="gnuplot_eigen_conformation.gpl"
gnuplot_energy_total="gnuplot_total_conformation.gpl"

# check paramset
if [ $# != 15 ]; then
    echo "run.sh Ds Dp Dd Z1s Z1p Z1d Z2s Z2p Z2d CAs CAp CAd Rs Rp Rd"
#                1  2  3  4   5   6   7   8   9   10  11  12  13 14 15
    exit 1
fi

Ds=$1
Dp=$2
Dd=$3
Z1s=$4
Z1p=$5
Z1d=$6
Z2s=$7
Z2p=$8
Z2d=$9
CAs=${10}
CAp=${11}
CAd=${12}
Rs=${13}
Rp=${14}
Rd=${15}

## calculate additional parameters
# chemical hardness is opposite to Dp
hardness=`echo "-($Dp)" | bc -l`

if [ "$doublezeta_n_s" == "" ]; then
    doublezeta_n_s="0"
fi
if [ "$doublezeta_n_p" == "" ]; then
    doublezeta_n_p="0"
fi
if [ "$doublezeta_n_d" == "" ]; then
    doublezeta_n_d="0"
fi

# c1,c2 for doublezeta
$elses_convdz $Z1s $Z2s $doublezeta_n_s $CAs > _tmp; read C1s C2s < _tmp; rm _tmp
$elses_convdz $Z1p $Z2p $doublezeta_n_p $CAp > _tmp; read C1p C2p < _tmp; rm _tmp
$elses_convdz $Z1d $Z2d $doublezeta_n_d $CAd > _tmp; read C1d C2d < _tmp; rm _tmp

# create element file with given paramset.
sed -e "s/VARIABLE_DIAG_S/$Ds/g" \
    -e "s/VARIABLE_DIAG_P/$Dp/g" \
    -e "s/VARIABLE_DIAG_D/$Dd/g" \
    -e "s/VARIABLE_ZETA_S/$Z1s/g" \
    -e "s/VARIABLE_ZETA_P/$Z1p/g" \
    -e "s/VARIABLE_ZETA_D/$Z1d/g" \
    -e "s/VARIABLE_ZETA2_S/$Z2s/g" \
    -e "s/VARIABLE_ZETA2_P/$Z2p/g" \
    -e "s/VARIABLE_ZETA2_D/$Z2d/g" \
    -e "s/VARIABLE_C1_S/$C1s/g" \
    -e "s/VARIABLE_C1_P/$C1p/g" \
    -e "s/VARIABLE_C1_D/$C1d/g" \
    -e "s/VARIABLE_C2_S/$C2s/g" \
    -e "s/VARIABLE_C2_P/$C2p/g" \
    -e "s/VARIABLE_C2_D/$C2d/g" \
    -e "s/VARIABLE_REPULSIVE_S/$Rs/g" \
    -e "s/VARIABLE_REPULSIVE_P/$Rp/g" \
    -e "s/VARIABLE_REPULSIVE_D/$Rd/g" \
    -e "s/VARIABLE_HARDNESS/$hardness/g" \
    template_${element}.xml > ${element}.xml

#########################################
## calculation of single point enegies with observed optimized structure

# generate target files
#echo molecule:${molecule}
#echo basetyep:${basetype}
#echo bond_optimized:${bond_optimized}
#echo angle_optimized:${angle_optimized}
#echo spin:${spin}

#exit

$elses_genmol $molecule elses $basetype $bond_optimized $angle_optimized $spin >  molecule.xml
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses_genmol"
    exit 1
fi

# run elses main program
#$elses config.xml > /dev/null >& /dev/null
$elses config.xml > elses.log
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses"
    exit 1
fi
if [ ! -r output_eigen_levels.txt ]; then
    echo "Error: failed in program `basename $elses`, parameter was not good"
    exit 1
fi
if [ ! -r output_wavefunction.txt ]; then
    echo "Error: failed in program `basename $elses`, parameter was not good"
    exit 1
fi
if [ ! -r Output.txt ]; then
    echo "Error: failed in program `basename $elses`, parameter was not good"
    exit 1
fi

# obtain eigen energies and their AO coefficients and total energy
# from gaussian output file and save these date to a file

$elses_obtdat elses molecule.xml output_eigen_levels.txt output_wavefunction.txt Output.txt > $output_data_optimized

# find the heightest eigen energy of orbitals that are occupied
energy_homo_optimized=`grep HOMO $output_data_optimized | awk '{ print $2;}'`
# obtain total energy
energy_total_optimized=`grep TOTAL $output_data_optimized | awk '{print $2}'`

rm -f output_eigen_levels.txt Output.txt
#rm -f log-node* output_tdos.txt restart.xml molecule.xml

#########################################
## calculation of some point enegies with some bonds/angles

# loop of some bonds/angles
echo "# bond, angle, eigen energies [a.u.]" > $output_energy_eigen_conformation
echo "# bond, angle, total energies [a.u.]" > $output_energy_total_conformation
echo "# MO match information" > ${output_MO_match}

#echo "ELSES calculation of bonds: $bonds_conformation, angles: $angles_conformation"
for bond in $bonds_conformation; do
    for angle in $angles_conformation; do
        # generate target files
	$elses_genmol $molecule elses $basetype $bond $angle $spin >  molecule.xml
	if [ ! $? == 0 ]; then
	    echo "Error: failed in program $elses_genmol"
	    exit 1
	fi

        # run elses main program
	$elses config.xml > /dev/null >& /dev/null
	if [ ! $? == 0 ]; then
	    echo "Error: failed in program $elses"
	    exit 1
	fi

	if [ ! -r output_eigen_levels.txt ]; then
	    echo "Error: failed in program `basename $elses`, parameter was not good"
	    exit 1
	fi

        # obtain eigen energies and their AO coefficients and total energy
        # from gaussian output file and save these date to a file
	$elses_obtdat elses molecule.xml output_eigen_levels.txt output_wavefunction.txt Output.txt > $output_data_conformation
	$elses_obtdat elses-sort $output_data_conformation ../gaussian/$output_data_optimized > _tmp

	# obtain the eigen energies and subtract them by the optimized homo energy
	energy_eigen=`grep EIGEN _tmp \
	    | awk "{ printf \" %+f\", \\$2-($energy_homo_optimized);}"`

        # obtain the total energy and subtract it by the optimized total energy
	energy_total=`grep TOTAL _tmp \
	    | awk "{print \\$2-($energy_total_optimized)}"`

	# obtain the MO matching information
	cp _tmp _tmp_MO_match
	#MO_mathcing=`grep 'MO matching' _tmp`
	#EIGEN_ENERGY_STATE=`grep EIGEN _tmp`

	rm _tmp

        # save eigen energies
	echo $bond $angle $energy_eigen >> $output_energy_eigen_conformation

        # save total energy
	echo $bond $angle $energy_total >> $output_energy_total_conformation

	# save MO matching information
	echo " " >> ${output_MO_match}
        echo ${bond} ${angle} >> ${output_MO_match}
	cat ${output_MO_match} _tmp_MO_match > _tmp
	cp _tmp ${output_MO_match}

	rm _tmp_MO_match
	rm _tmp



    done

    if [ `echo $angles_conformation | wc -w` -gt 1 ]; then
	echo "" >> $output_energy_eigen_conformation
	echo "" >> $output_energy_total_conformation
    fi
done

#echo "eigen energies are saved in file " $output_energy_eigen_conformation
#echo "total energies are saved in file " $output_energy_total_conformation



if [ ! -r ../gaussian/$output_energy_eigen_conformation ]; then
    echo "Error: not found data file", ../gaussian/$output_energy_eigen_conformation
    exit 1
fi
#echo -n "diff in eigen energies: "
$elses_comdat $output_energy_eigen_conformation ../gaussian/$output_energy_eigen_conformation ${nbands}  > $output_energy_difference
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses_comdat"
    exit 1
fi

if [ ! -r ../gaussian/$output_energy_total_conformation ]; then
    echo "Error: not found data file", ../gaussian/$output_energy_total_conformation
    exit 1
fi
#echo -n "diff in total energies: "
$elses_comdat $output_energy_total_conformation ../gaussian/$output_energy_total_conformation 1 >>  $output_energy_difference
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses_comdat"
    exit 1
fi

num_bands=`grep NUM ../gaussian/$output_data_optimized | awk '{print $2}'`

# create gnuplot script files
if [ `echo $bonds_conformation | wc -w` -gt 1 ]; then
    if [ `echo $angles_conformation | wc -w` -gt 1 ]; then
#--------------
	echo -e \
	    "\nN=$num_bands" \
	    "\nset grid" \
	    "\nset xlabel \"Bond length [angstrom]\"" \
	    "\nset ylabel \"Bond angle  [degree]\"" \
	    "\nset zlabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Eigen energies by ELSES and Gaussian\"\n" \
	    > $gnuplot_energy_eigen

	for i in `seq 1 $num_bands`; do
	    if [ $i == "1" ]; then
		echo " splot \"$output_energy_eigen_conformation\"             u 1:2:2+$i  w lp pt 1 ps 2 lw 2 lt 1 t \"ELSES\""
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 1:2:2+$i  w lp pt 2 ps 2 lw 2 lt 2 t \"Gaussian\""
	    else
		echo "replot \"$output_energy_eigen_conformation\"             u 1:2:2+$i  w lp pt 1 ps 2 lw 2 lt 1 notitle"
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 1:2:2+$i  w lp pt 2 ps 2 lw 2 lt 2 notitle"
	    fi
	done >> $gnuplot_energy_eigen

	echo -e \
	    "\npause -1 \"hit return key\"" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_eigen%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    >> $gnuplot_energy_eigen

	echo -e \
	    "\nset grid" \
	    "\nset xlabel \"Bond length [angstrom]\"" \
	    "\nset ylabel \"Bond angle  [degree]\"" \
	    "\nset zlabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Total energies by ELSES and Gaussian\"" \
	    "\n" \
	    "\n splot \"$output_energy_total_conformation\"             u 1:2:3 w lp ps 2 lw 2 t \"ELSES\"" \
	    "\nreplot \"../gaussian/$output_energy_total_conformation\" u 1:2:3 w lp ps 2 lw 2 t \"Gaussian\"" \
	    "\npause -1 \"hit return key\"" \
	    "\n" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_total%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    > $gnuplot_energy_total

#--------------
    else
#--------------
	echo -e \
	    "\nN=$num_bands" \
	    "\nset grid" \
	    "\nset xlabel \"Bond length [angstrom]\"" \
	    "\nset ylabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Eigen energies by ELSES and Gaussian\"\n" \
	    > $gnuplot_energy_eigen

	for i in `seq 1 $num_bands`; do
	    if [ $i == "1" ]; then
		echo "  plot \"$output_energy_eigen_conformation\"             u 1:2+$i  w lp pt 1 ps 2 lw 2 lt 1 t \"ELSES\""
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 1:2+$i  w lp pt 2 ps 2 lw 2 lt 2 t \"Gaussian\""
	    else
		echo "replot \"$output_energy_eigen_conformation\"             u 1:2+$i  w lp pt 1 ps 2 lw 2 lt 1 notitle"
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 1:2+$i  w lp pt 2 ps 2 lw 2 lt 2 notitle"
	    fi
	done >> $gnuplot_energy_eigen

	echo -e \
	    "\npause -1 \"hit return key\"" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_eigen%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    >> $gnuplot_energy_eigen

	echo -e \
	    "\nset grid" \
	    "\nset xlabel \"Bond length [angstrom]\"" \
	    "\nset ylabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Total energies by ELSES and Gaussian\"" \
	    "\n" \
	    "\n  plot \"$output_energy_total_conformation\"             u 1:3 w lp ps 2 lw 2 t \"ELSES\"" \
	    "\nreplot \"../gaussian/$output_energy_total_conformation\" u 1:3 w lp ps 2 lw 2 t \"Gaussian\"" \
	    "\npause -1 \"hit return key\"" \
	    "\n" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_total%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    > $gnuplot_energy_total
#--------------
    fi
else
    if [ `echo $angles_conformation | wc -w` -gt 1 ]; then
#--------------
	echo -e \
	    "\nN=$num_bands" \
	    "\nset grid" \
	    "\nset xlabel \"Bond angle  [degree]\"" \
	    "\nset ylabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Eigen energies by ELSES and Gaussian\"\n" \
	    > $gnuplot_energy_eigen

	for i in `seq 1 $num_bands`; do
	    if [ $i == "1" ]; then
		echo "  plot \"$output_energy_eigen_conformation\"             u 2:2+$i  w lp pt 1 ps 2 lw 2 lt 1 t \"ELSES\""
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 2:2+$i  w lp pt 2 ps 2 lw 2 lt 2 t \"Gaussian\""
	    else
		echo "replot \"$output_energy_eigen_conformation\"             u 2:2+$i  w lp pt 1 ps 2 lw 2 lt 1 notitle"
		echo "replot \"../gaussian/$output_energy_eigen_conformation\" u 2:2+$i  w lp pt 2 ps 2 lw 2 lt 2 notitle"
	    fi
	done >> $gnuplot_energy_eigen

	echo -e \
	    "\npause -1 \"hit return key\"" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_eigen%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    >> $gnuplot_energy_eigen

	echo -e \
	    "\nset grid" \
	    "\nset xlabel \"Bond angle  [degree]\"" \
	    "\nset ylabel \"Eigen energy [a.u.]\"" \
	    "\nset title  \"Total energies by ELSES and Gaussian\"" \
	    "\n" \
	    "\n  plot \"$output_energy_total_conformation\"             u 2:3 w lp lw 2 ps 2 t \"ELSES\"" \
	    "\nreplot \"../gaussian/$output_energy_total_conformation\" u 2:3 w lp lw 2 ps 2 t \"Gaussian\"" \
	    "\npause -1 \"hit return key\"" \
	    "\n" \
	    "\nset terminal post color eps \"Arial\"" \
	    "\nset output \"${gnuplot_energy_total%.gpl}.eps\"" \
	    "\nreplot" \
	    "\nset ter x11\n" \
	    > $gnuplot_energy_total
#--------------
    else
    echo "only one data, no plot"
    fi
fi

#rm -f output_eigen_levels.txt Output.txt molecule.xml log-node* output_tdos.txt

exit 0

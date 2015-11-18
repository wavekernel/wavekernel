#!/bin/sh -

#########################################
## script for evaluating Gaussian with paramset.

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

if [ "$basetype" == "" ]; then
    basetype="STO-3G"
fi

# output of Gaussian
output_data_optimized="data_optimized.dat"
output_data_conformation="data_conformation.dat"
output_energy_eigen_conformation="energy_eigen_conformation.dat"
output_energy_total_conformation="energy_total_conformation.dat"

#########################################
## calculation of single point enegies with observed optimized structure

# generate target files
$elses_genmol $molecule gaussian $basetype $bond_optimized $angle_optimized $spin > molecule.com
if [ ! $? == 0 ]; then
    echo "Error: failed in program $elses_genmol"
    exit 1
fi

# run gaussian
$gaussian < molecule.com > molecule.out
if [ ! $? == 0 ]; then
    echo "Error: failed in program $gaussian"
    exit 1
fi

# obtain eigen energies and their AO coefficients and total energy
# from gaussian output file and save these date to a file
$elses_obtdat gaussian molecule.out $basetype > $output_data_optimized

# obtain the HOMO energy
energy_homo_optimized=`grep HOMO $output_data_optimized | awk '{ print $2;}'`
# obtain the total energy
energy_total_optimized=`grep TOTAL $output_data_optimized | awk '{print $2}'`

#rm -f molecule.com molecule.out *.FChk fort.7

#########################################
## calculation of some point enegies with some bonds/angles

# loop of some bonds/angles
echo "# bond, angle, eigen energies [a.u.]" > $output_energy_eigen_conformation
echo "# bond, angle, total energies [a.u.]" > $output_energy_total_conformation

# echo "Gaussian calculation of bonds: $bonds_conformation, angles: $angles_conformation"
for bond in $bonds_conformation; do
    for angle in $angles_conformation; do
        # generate target files
	$elses_genmol $molecule gaussian $basetype $bond $angle $spin  > molecule.com
	if [ ! $? == 0 ]; then
	    echo "Error: failed in program $elses_genmol"
	    exit 1
	fi

        # run gaussian
	$gaussian < molecule.com > molecule.out
	if [ ! $? == 0 ]; then
	    echo "Error: failed in program $gaussian"
	    exit 1
	fi

        # convert eigen energies
        # obtain eigen energies and their AO coefficients from gaussian output file
        # and save these date to a file
        $elses_obtdat gaussian molecule.out $basetype > $output_data_conformation

	# obtain the eigen energies and subtract them by the optimized homo energy
	energy_eigen=`grep EIGEN $output_data_conformation \
	    | awk "{ printf \" %+f\", \\$2-($energy_homo_optimized);}"`

        # obtain the total energy and subtract it by the optimized total energy
	energy_total=`grep TOTAL $output_data_conformation \
	    | awk "{print \\$2-($energy_total_optimized)}"`

        # save eigen energies
	echo $bond $angle $energy_eigen >> $output_energy_eigen_conformation

        # save total energy
	echo $bond $angle $energy_total >> $output_energy_total_conformation
    done

    if [ `echo $angles_conformation | wc -w` -gt 1 ]; then
	echo "" >> $output_energy_eigen_conformation
	echo "" >> $output_energy_total_conformation
    fi
done

#echo "eigen energies are saved in file " $output_energy_eigen_conformation
#echo "total energies are saved in file " $output_energy_total_conformation

#rm -f molecule.com molecule.out *.FChk fort.7

exit 0

#!/bin/sh
echo 'Sample test'
echo 'rm output*.txt'
rm output*.txt
echo 'rm Output.txt'
rm Output.txt
echo 'rm log*.txt'
rm log*.txt
echo 'rm *.cube*'
rm *.cube*
echo '../../../bin/elses-xml-generate generate.xml C6H6.xml > log-xml-gen.txt'
../../../bin/elses-xml-generate generate.xml C6H6.xml > log-xml-gen.txt
echo '../../shell_scripts/elses_test.sh 0'
../../shell_scripts/elses_test.sh 0
echo '../../../bin/elses-generate-cubefile 1 3 > log-gen-cube.txt'
../../../bin/elses-generate-cubefile 1 3 > log-gen-cube.txt
echo 'head -n 100 eigen_state_000001.cube > eigen_state_000001.cube.header'
head -n 100 eigen_state_000001.cube > eigen_state_000001.cube.header
echo 'head -n 100 eigen_state_000002.cube > eigen_state_000002.cube.header'
head -n 100 eigen_state_000002.cube > eigen_state_000002.cube.header
echo 'head -n 100 eigen_state_000003.cube > eigen_state_000003.cube.header'
head -n 100 eigen_state_000003.cube > eigen_state_000003.cube.header
echo 'gzip *.cube'
gzip *.cube
echo 'Copy the result files into ./result_present directory'
cp -p Output.txt ./result_present
cp -p *.cube.header ./result_present
cp -p log*.txt ./result_present
cp -p output*.txt ./result_present
cp -p position.xyz ./result_present
cp -p restart.xml ./result_present

#!/bin/sh -

echo "#### Updating utility programs"
(cd util; make)

#echo "#### Running Gaussian"
#(cd gaussian; ./run.sh)

echo "#### Running ELSES"
(cd elses; ./optimize.sh)

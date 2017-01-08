#!/bin/sh -

# include settings
if [ ! -r ../setting.sh ]; then
    echo "Error: not found setting.sh in $PWD/.."
    exit 1
else
    . ../setting.sh
fi

$elses_optpar | tee optimize.log


#! /bin/bash

rm *_
rm *.csv
rm plot*
rm taurenmd.log
rm traj_OUTPUT.dcd

source activate taurenmd

python ../taurenmd.py -c test_1.json

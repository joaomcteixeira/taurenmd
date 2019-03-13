#! /bin/bash

rm _*
rm *.csv
rm plot*
rm taurenmd.log
rm traj_OUTPUT.dcd

source activate taurenmd

# basic test with RMSDS using MDTraj
python ../taurenmd.py -c test_1.json
# basic test with RMSDs using MDAnalysis
python ../taurenmd.py -c test_2.json
# basic test with atom_selection selecting chain 1 - MDTraj
python ../taurenmd.py -c test_3.json
# basic test with atom_selection selecting chain 1 - MDAnalysis
python ../taurenmd.py -c test_4.json
# savetraj with atom and frame selection - MDTraj
python ../taurenmd.py -c test_5.json
# savetraj with atom and frame selection - MDAnalysis
python ../taurenmd.py -c test_6.json
# combined RMSDs with atom selection (chain) - MDTraj
python ../taurenmd.py -c test_7.json
# same as above with empty selection
python ../taurenmd.py -c test_8.json
# combined RMSDs with atom selection (chain) - MDAnalysis
python ../taurenmd.py -c test_9.json
# combined RMSDs with atom selection (resid) - MDAnalysis
python ../taurenmd.py -c test_11.json

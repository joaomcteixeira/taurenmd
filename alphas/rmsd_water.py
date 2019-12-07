import sys

import mdtraj as md
from taurenmd.libs import libmdt
from bioplottemplates.plots import param


t1 = libmdt.mdtraj_load_traj(sys.argv[1], sys.argv[2])
t2 = libmdt.mdtraj_load_traj(sys.argv[1], sys.argv[3])

print(t1.n_frames)
t3 = t1.join(t2)
print(t3.n_frames)

atom_select = t3.top.select('water and not name DUM')
t3.atom_slice(atom_select, inplace=True)

rmsds = md.rmsd(t3, t3, frame=547)

param.plot(
    list(range(t3.n_frames)),
    rmsds)

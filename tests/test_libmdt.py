"""Test libmdt."""
import mdtraj as md

from taurenmd.libs import libmdt

from . import toptest, toptest_cif, trajtest


def test_load_traj_pdb():
    """Test loading traj."""
    traj = libmdt.load_traj(toptest, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10


def test_load_traj_cif():
    """Test loading traj."""
    traj = libmdt.load_traj(toptest_cif, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10

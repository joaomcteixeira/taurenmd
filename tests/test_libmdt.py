"""Test libmdt."""
import mdtraj as md

from taurenmd.libs import libmdt

from . import toptest, toptest_cif, trajtest


def test_load_traj_pdb():
    """Test loading traj."""
    traj = libmdt.load_traj(toptest, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10


def test_load_traj_pdb_2():
    """Test loading traj."""
    traj = libmdt.load_traj(toptest, [trajtest, trajtest])
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 20


def test_load_traj_pdb_3():
    """Test loading traj."""
    traj = libmdt.load_traj(toptest, [trajtest, trajtest], insort=True)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 20


def test_load_traj_cif():
    """Test loading traj."""
    assert toptest_cif.suffix == '.cif'
    traj = libmdt.load_traj(toptest_cif, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10

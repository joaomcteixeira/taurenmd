"""Test libmdt."""
import mdtraj as md
import pytest

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
    traj = libmdt.load_traj(toptest_cif, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10


def test_attempt_load_cif_SIMTK_1():
    """Test attempt load CIF."""
    libmdt.attempt_to_load_top_from_simtk(toptest_cif)


def test_attempt_load_cif_SIMTK_2():
    """Test attempt load CIF."""
    assert toptest.str() == libmdt.attempt_to_load_top_from_simtk(toptest)


def test_load_traj_cif_import_error():
    """Test loading traj."""
    libmdt.SIMTK = False
    with pytest.raises(SystemExit) as err:
        libmdt.load_traj(toptest_cif, trajtest)
    assert err.type == SystemExit
    assert err.value.code == 0


def test_simtk_import_error():
    """Test import error message."""
    with pytest.raises(SystemExit) as err:
        libmdt._log_simtkimport_error()
    
    assert err.type == SystemExit
    assert err.value.code == 0

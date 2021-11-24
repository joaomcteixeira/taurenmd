"""Test libmdt."""
import importlib
import sys

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
    assert toptest_cif.suffix == '.cif'
    traj = libmdt.load_traj(toptest_cif, trajtest)
    assert isinstance(traj, md.Trajectory)
    assert len(traj) == 10


def test_load_traj_cif_import_error():
    """Test loadtraj_cif import error as if simtk was not installed."""
    sys.modules['openmm.app.pdbxfile'] = None
    importlib.reload(libmdt)
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

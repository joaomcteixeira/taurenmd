"""Test libopenmm."""
import importlib
import sys

import pytest
from openmm.app import pdbxfile

from taurenmd.libs import libopenmm

from . import toptest, toptest_cif, trajtest


def test_simtk():
    """Test simtk is installed."""
    assert libopenmm.SIMTK


def test_attempt_load_cif_SIMTK_1():
    """Test attempt load CIF."""
    mol = libopenmm.attempt_to_load_top_from_simtk(toptest_cif)
    assert isinstance(mol, pdbxfile.PDBxFile)


def test_attempt_load_cif_SIMTK_2():
    """Test attempt load CIF."""
    with pytest.raises(ValueError):
        libopenmm.attempt_to_load_top_from_simtk(toptest)


def test_simtk_import_error():
    """Test import error message."""
    with pytest.raises(SystemExit) as err:
        libopenmm._log_simtkimport_error()
    assert err.type == SystemExit
    assert err.value.code == 0


def test_without_simtk():
    """Test without simtk installed."""
    sys.modules['openmm.app.pdbxfile'] = None
    importlib.reload(libopenmm)
    with pytest.raises(SystemExit) as err:
        libopenmm.attempt_to_load_top_from_simtk(toptest_cif)
    assert err.type == SystemExit
    assert err.value.code == 0


def test_load_traj_cif_import_error():
    """Test loadtraj_cif import error as if simtk was not installed."""
    sys.modules['openmm.app.pdbxfile'] = None
    from taurenmd.libs import libmdt
    importlib.reload(libmdt)
    with pytest.raises(SystemExit) as err:
        libmdt.load_traj(toptest_cif, trajtest)
    assert err.type == SystemExit
    assert err.value.code == 0

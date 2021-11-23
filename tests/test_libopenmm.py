"""Test libopenmm."""
import pytest
from openmm.app import pdbxfile

from taurenmd.libs import libopenmm

from . import toptest, toptest_cif


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

"""Test libmda."""
import copy

import MDAnalysis as mda
import numpy as np
import pytest

from taurenmd.libs import libmda as la

from . import toptest, trajtest


def test_load_universe_1():
    """Test load MDA Universe."""
    universe = la.load_universe(
        toptest,
        trajtest,
        )
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 10


def test_load_universe_2():
    """Test load MDA Universe multiple trajs."""
    universe = la.load_universe(
        toptest,
        trajtest,
        trajtest,
        trajtest,
        )
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 30


def test_load_universe_3():
    """Test load MDA Universe multiple trajs."""
    mapp = map(str, [trajtest] * 3)
    universe = la.load_universe(toptest, *mapp)
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 30


def test_report_universe():
    """
    Test universe report.

    Only interface, does not test return value.
    """
    universe = mda.Universe(
        toptest.str(),
        trajtest.str(),
        )
    la.report(universe)


def test_mdaalignto():
    """Test MDA alignto."""
    universe = mda.Universe(toptest.str(), trajtest.str())
    reference = mda.Universe(toptest.str())
    ugroup = universe.select_atoms('all')

    for _ts in universe.trajectory[1:]:
        prev = copy.copy(ugroup.positions)
        assert np.all(np.equal(prev, ugroup.positions))
        la.mdaalignto(universe, reference, 'name CA')
        assert np.all(np.not_equal(prev, ugroup.positions))


def test_mdaalignto_exit():
    """Test exit on ZeroDivision."""
    universe = mda.Universe(toptest.str(), trajtest.str())
    reference = mda.Universe(toptest.str())
    with pytest.raises((ValueError, ZeroDivisionError)):
        la.mdaalignto(universe, reference, 'name ZZZZZ')


def test_draw_atom_label():
    """Test label maker."""
    universe = mda.Universe(
        toptest.str(),
        trajtest.str(),
        )

    atom_group = universe.select_atoms('resnum 10-12 and name CA')

    labels = la.draw_atom_label_from_atom_group(atom_group)
    assert isinstance(labels, list)
    assert labels == ['A.10Ser.CA', 'A.11Ile.CA', 'A.12Leu.CA']

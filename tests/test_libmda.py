"""Test libmda."""
import copy

import MDAnalysis as mda
import numpy as np
import pytest

from taurenmd.libs import libmda as la

from . import toptest, toptest_cif, trajtest


def test_load_universe_1():
    """Test load MDA Universe."""
    universe = la.load_universe(
        toptest,
        trajtest,
        )
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 10


def test_load_universe_1_cif():
    """Test load MDA Universe."""
    universe = la.load_universe(
        toptest_cif,
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


@pytest.mark.parametrize(
    'x, frame',
     (('1', 1),
      ('-1', -1),
      ('1ns', 1000),
      ('1.5ns', 1500),
      ('12e3', 12e3),
      ('12e3ps', 12e3)),
    )
def test_convert_str_time(x, frame):
    """
    Test convert string to time.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """
    assert frame == la.convert_time_to_frame(x, dt=1)


def test_convert_str_time_dt():
    """
    Test convert string to time in ps.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """
    assert 1 == la.convert_time_to_frame("10ps", dt=10)


@pytest.mark.parametrize(
    'value',
    [
        '0.1',
        'nothing',
        ],
    )
def test_convert_str_time_raise(value):
    """
    Test convert string to time ValueError.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """
    with pytest.raises(ValueError):
        la.convert_time_to_frame(value, dt=1)

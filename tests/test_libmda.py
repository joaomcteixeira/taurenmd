"""Test libmda."""
import copy
from math import isclose

import MDAnalysis as mda
import numpy as np
import pytest

from taurenmd.libs import libmda as la

from . import toptest, toptest_cif, trajtest


@pytest.fixture
def universe():
    """Return a Universe."""
    return la.load_universe(toptest, trajtest)


def test_load_universe_1(universe):
    """Test load MDA Universe."""
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 10


def test_load_universe_1_cif():
    """Test load MDA Universe."""
    universe = la.load_universe(toptest_cif, trajtest)
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


def test_load_universe_2_insort():
    """Test load MDA Universe multiple trajs."""
    universe = la.load_universe(
        toptest,
        trajtest,
        trajtest,
        trajtest,
        insort=True,
        )
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 30


def test_load_universe_3():
    """Test load MDA Universe multiple trajs."""
    mapp = map(str, [trajtest] * 3)
    universe = la.load_universe(toptest, *mapp)
    assert isinstance(universe, mda.Universe)
    assert len(universe.trajectory) == 30


def test_report_universe(universe):
    """
    Test universe report.

    Only interface, does not test return value.
    """
    la.report(universe)


def test_mdaalignto(universe):
    """Test MDA alignto."""
    reference = mda.Universe(toptest.str())
    ugroup = universe.select_atoms('all')

    for _ts in universe.trajectory[1:]:
        prev = copy.copy(ugroup.positions)
        assert np.all(np.equal(prev, ugroup.positions))
        la.mdaalignto(universe, reference, 'name CA')
        assert np.all(np.not_equal(prev, ugroup.positions))


def test_mdaalignto_exit(universe):
    """Test exit on ZeroDivision."""
    reference = mda.Universe(toptest.str())
    with pytest.raises((ValueError, ZeroDivisionError)):
        la.mdaalignto(universe, reference, 'name ZZZZZ')


def test_draw_atom_label(universe):
    """Test label maker."""
    atom_group = universe.select_atoms('resnum 10-12 and name CA')
    labels = la.draw_atom_label_from_atom_group(atom_group)
    assert isinstance(labels, list)
    assert labels == ['A.10Ser.CA', 'A.11Ile.CA', 'A.12Leu.CA']


@pytest.mark.parametrize(
    'x, frame',
    (
        ('1', 1),
        ('-1', -1),
        ('1ns', 1000),
        ('1.5ns', 1500),
        ('12e3', 12e3),
        ('12e3ps', 12e3),
        ),
    )
def test_convert_str_time(x, frame):
    """
    Test convert string to time.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """  # noqa: E501
    assert frame == la.convert_time_to_frame(x, dt=1)


@pytest.mark.parametrize(
    'time_,dt,expected',
    [
        ('10ps', 10, 1),
        ('0.05ns', 50, 1),
        ]
    )
def test_convert_str_time_dt(time_, dt, expected):
    """
    Test convert string to time in ps.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """  # noqa: E501
    assert expected == la.convert_time_to_frame(time_, dt=dt)


@pytest.mark.parametrize(
    'value',
    [
        'nothing',
        ''
        ],
    )
def test_convert_str_time_raise(value):
    """
    Test convert string to time ValueError.

    Tests taken from: https://github.com/MDAnalysis/mdacli/blob/15f6981df7b14ef5d52d64b56953d276291068ab/tests/test_utils.py#L22-L41
    """  # noqa: E501
    with pytest.raises(ValueError):
        la.convert_time_to_frame(value, dt=1)


@pytest.mark.parametrize(
    'timestep,expected',
    [
        ('ns', [0.001, 0.002, 0.003, 0.004, 0.005]),
        ('ps', [1, 2, 3, 4, 5]),
        ]
    )
def test_create_x_data(universe, timestep, expected):
    """Test create x_data."""
    xdata, xlabel = la.create_x_data(universe, timestep, list(range(1, 6)))
    assert xlabel == f'Time ({timestep})'
    assert all(isclose(i, j, rel_tol=0.00001) for i, j in zip(xdata, expected))


def test_create_x_data_w_frames(universe):
    """Test create x_data."""
    xdata, xlabel = la.create_x_data(universe, False, list(range(10)))
    assert xlabel == 'Frames'
    expected = list(range(10))
    assert all(isclose(i, j, rel_tol=0.001) for i, j in zip(xdata, expected))


@pytest.mark.parametrize(
    'f1, f2, expected',
    [
        (0, -1, 9),
        (0, 1, 1),
        (0, 0, 0),
        ]
    )
def test_get_time(universe, f1, f2, expected):
    """Test get time from Univerise."""
    result = la.get_timestep(universe, i0=f1, i1=f2)
    assert isclose(result, expected, rel_tol=0.001)


@pytest.mark.parametrize(
    'slc,expected',
    [
        (slice(3, None), [3, 4, 5, 6, 7, 8, 9]),
        (slice(None, None), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        (slice(None, 3), [0, 1, 2]),
        (slice(3, None, 3), [3, 6, 9]),
        ]
    )
def test_get_frame_list_from_slice(universe, slc, expected):
    """Test get a frame list from slice."""
    result = la.get_frame_list_from_slice(universe, slc)
    assert result == expected


@pytest.mark.parametrize(
    's,expected',
    [
        ('0ps', 0),
        ('1ps', 1),
        ('2ps', 2),
        ('3ps', 3),
        ('4ps', 4),
        ('5ps', 5),
        ('6ps', 6),
        ('7ps', 7),
        ('8ps', 8),
        ('9ps', 9),
        ('10ps', 10),
        ('0.005ns', 5),
        ('4', 4),
        (3, 3),
        (0, 0),
        ]
    )
def test_get_time_or_frame_to_frame(universe, s, expected):
    """Test get frame."""
    result = la.convert_time_or_frame_to_frame(s, universe)
    assert result == expected


@pytest.mark.parametrize(
    'start,stop,step,expected',
    [
        (None, None, 1, slice(None, None, 1)),
        (2, 5, 1, slice(2, 5, 1)),
        (8, 15, 1, slice(8, 15, 1)),
        ('2ps', None, None, slice(2, None, None)),
        ('0.004ns', '8ps', None, slice(4, 8, None)),
        ]
    )
def test_get_frame_slices(universe, start, stop, step, expected):
    """Test get frame slices."""
    result = la.get_frame_slices(
        universe,
        start=start,
        stop=stop,
        step=step,
        )

    assert result == expected

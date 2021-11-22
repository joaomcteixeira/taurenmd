"""Test libcalc."""
from math import sqrt

import MDAnalysis as mda
import numpy as np
import pytest

from taurenmd.libs import libcalc as lc

from . import toptest, trajtest


def test_mda_rmsd():
    """Test MDA RMSDS."""
    u = mda.Universe(toptest.str(), trajtest.str())
    result = lc.mda_rmsd(u)
    assert result.shape == (10,)
    assert result.dtype == np.float


def test_mda_rmsd_2():
    """Test MDA RMSDS."""
    u = mda.Universe(toptest.str(), trajtest.str())
    result = lc.mda_rmsd(u, frame_slice='::2')
    assert result.shape == (5,)
    assert result.dtype == np.float


def test_mda_rmsf():
    """Test MDA RMSF."""
    u = mda.Universe(toptest.str(), trajtest.str())
    CA = u.select_atoms('name CA')
    result = lc.mda_rmsf(CA, frame_slice=None)
    assert isinstance(result, np.ndarray)
    assert result.size == len(CA)
    assert result.dtype == np.float


@pytest.mark.parametrize(
    'p1,p2,p3,expected',
    [
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
            ),
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            np.array([-1.0, 0.0, 0.0])
            ),
        ],
    )
def test_calc_plane_normal(p1, p2, p3, expected):
    """Test calc plane normal."""
    result = lc.calc_plane_normal(p1, p2, p3)
    assert np.all(np.equal(result, expected))


@pytest.mark.parametrize(
    'p1,p2,p3,expected',
    [
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            (0.0, 0.0, 1.0, 0.0),
            ),
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            (-1.0, 0.0, 0.0, 0.0)
            ),
        ],
    )
def test_calc_plane_eq(p1, p2, p3, expected):
    """Test calc plane normal."""
    result = lc.calc_plane_eq(p1, p2, p3)
    print(result)
    assert all(abs(a - b) < 0.00001 for a, b in zip(result, expected))


@pytest.mark.parametrize(
    'a1,b1,c1,a2,b2,c2,aunit,expected',
    [
        (0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 'degrees', 90.0),
        (0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 'radians', 1.5708),
        (0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 'degrees', 90.0),
        (1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 'degrees', 0.0),
        (1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 'radians', 0.0),
        ],
    )
def test_calc_planes_angle(a1, b1, c1, a2, b2, c2, aunit, expected):
    """Test angle between planes."""
    result = lc.calc_planes_angle(a1, b1, c1, a2, b2, c2, aunit)
    assert abs(result - expected) < 0.00001


@pytest.mark.parametrize(
    'coords, expected',
    [
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1]]), np.array([90], dtype=float)),  # noqa: E501
        (np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]), np.array([90, -90], dtype=float)),  # noqa: E501
        (
            np.array((
                (1, 0, 0.000),
                (0, 0, 0.000),
                (0, 1, 0.000),
                (-1, 1, 0.000),
                (-1, 1, 1),
                (0, 1, 1),
                (0, 0, 0)
                )),
            np.array([-180, 90, -0, -45], dtype=float),
            ),
        (
            np.array((
                (1, 0, 0.000),
                (0, 0, 0.000),
                (0, 1, 0.000),
                (-1, 1, 0.000),
                (-1, 1, 1),
                )),
            np.array([-180, 90], dtype=float),
            ),
        ]
    )
def test_calculate_torsions(coords, expected):
    """Tests torsion calculation."""
    result = np.degrees(lc.calc_torsion_angles(coords))
    assert np.all(np.equal(result, expected))


@pytest.mark.parametrize(
    'p1,p2,p3,p4,expected',
    [
        (
            np.array([1, 0, 0]),
            np.array([0, 0, 0]),
            np.array([0, 1, 0]),
            np.array([
                [1, 1, 0],  # 0
                [0, 1, 1],
                [sqrt(1), 1, sqrt(1)],
                ]),
            np.array([0, -90, -45]),
            )
        ]
    )
def test_calculate_torsions_set(p1, p2, p3, p4, expected):
    """Tests torsion calculation."""
    result = np.degrees(lc.torsion_set(p1, p2, p3, p4))
    assert np.all(np.equal(result, expected))

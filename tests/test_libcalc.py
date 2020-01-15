"""Test libcalc."""
import MDAnalysis as mda
import numpy as np
import pytest
from pyquaternion import Quaternion as Q

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
            np.array([0.0, 0.0, -1.0])
            ),
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            np.array([1.0, 0.0, 0.0])
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
            (0.0, 0.0, -1.0, 0.0),
            ),
        (
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, -1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
            (1.0, 0.0, 0.0, 0.0)
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


def test_gen_quaternion_rot():
    """Test quaternion rotation."""
    result = lc.generate_quaternion_rotations(
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        num=360
        )
    assert isinstance(result, list)
    assert len(result) == 360
    assert isinstance(result[0], tuple)
    assert len(result[0]) == 2
    assert isinstance(result[0][0], Q)
    assert isinstance(result[0][1], Q)


def test_sort_by_minimum():
    """Test sort quaternion rotations."""
    rotations = lc.generate_quaternion_rotations(
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        num=10,
        endpoint=False,
        )
    for i in range(len(rotations)):
        print(rotations[i][0].degrees)
    result = lc.sort_by_minimum_Qdistances(
        rotations,
        np.array([1.0, 0.0, 0.0]),
        )

    for i in range(len(result)):
        print(result[i][0].degrees)

    assert isinstance(result, list)
    assert isinstance(result[0], tuple)
    assert len(result[0]) == 2
    assert isinstance(result[0][0], Q)
    assert isinstance(result[0][1], Q)
    assert result[0][0].degrees < 0.001
    assert abs(result[0][-1].degrees - 180) < 0.001

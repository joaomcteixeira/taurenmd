"""Test libutil."""
from math import isclose

import pytest

from taurenmd.libs import libutil


@pytest.mark.parametrize(
    's,expected',
    [
        ['00', (0.0, '')],
        ['0', (0.0, '')],
        ['-0', (0.0, '')],
        ['0.', (0.0, '')],
        ['.0ns', (0.0, 'ns')],
        ['12.3ps', (12.3, 'ps')],
        ['12.34minute', (12.34, 'minute')],
        ['-12.34', (-12.34, '')],
        ['10E40sd', (float('10E40'), 'sd')],
        ['10e40sd', (float('10e40'), 'sd')],
        ['1E4', (float('1E4'), '')],
        ['1E-4hours', (float('1E-4'), 'hours')],
        ['-10E40nanosecond', (float('-10E40'), 'nanosecond')],
        ['-10E-40', (float('-10E-40'), '')],
        ['-10e-40', (float('-10e-40'), '')],
        ['-.1E-4', (float('-.1E-4'), '')],
        ['10.2E30', (float('10.2E30'), '')],
        ['10.2e30', (float('10.2e30'), '')],
        ['.10E30', (float('.10E30'), '')],
        ['-10.2E30', (float('-10.2E30'), '')],
        ['-.10E30', (float('-.10E30'), '')],
        ]
    )
def test_string_to_timestep(s, expected):
    """Test correct conversion."""
    time_, unit_ = libutil.split_time_unit(s)
    assert isclose(time_, expected[0])
    assert unit_ == expected[1]


@pytest.mark.parametrize(
    's',
    [
        '12.34 # with comment',
        '12.34.12fdsfdsfds',
        '1E4.4',
        '.10E30E',
        '10.2.E30',
        '123#with wrong comment',
        'E10',
        'e10',
        '.E10',
        '-e',
        '-E10',
        '10-10E19',
        ]
    )
def test_string_to_timestep_wrong(s):
    """Test correct conversion."""
    with pytest.raises(ValueError):
        libutil.split_time_unit(s)

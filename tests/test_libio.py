"""Test libio."""
import pytest

from taurenmd.core import Path
from taurenmd.libs import libio as io

from . import datap, export_data_expected, export_data_expected_2

@pytest.mark.parametrize(
    'in1,expected',
    [
        ('traj_1.dcd', 1),
        ('traj2.dcd', 2),
        ('tr4j2.dcd', 2),
        ('traj_3.dcd', 3),
        ('traj_1231.dcd', 1231),
        ('traj_0011.dcd', 11),
        ('traj_1_.dcd', 1),
        ('traj_20200101_1.dcd', 1),
        ],
    )
def test_get_number(in1, expected):
    """Test get number from path."""
    result = io.get_number(in1)
    assert result == expected


@pytest.mark.parametrize(
    'in1,expected',
    [
        (['t_1.dcd', 't_11.dcd', 't_2.dcd'], ['t_1.dcd', 't_2.dcd', 't_11.dcd']),
        ],
    )
def test_sort_numbered_input_1(in1, expected):
    """Test sort numbered inputs."""
    result = io.sort_numbered_input(*in1)
    assert result == expected

@pytest.mark.parametrize(
    'in1,error',
    [
        (['t_1.dcd', 't_11.dcd', 't_2.dcd'], ValueError),
        ]
    )
def test_sort_numbered_inputs_error(in1, error):
    """Test sort numbered inputs raised Errors."""
    with pytest.raises(error):
        io.sort_numbered_input(in1)


def test_report_input():
    """
    Test report_input.

    Tests only the interface, there is no return value.
    """
    io.report_input('topology.pdb', 'traj1.dcd')


@pytest.mark.parametrize(
    'ipath,prefix,expected',
    [
        ('traj_output.xtc', 'my_prefix', 'my_prefixtraj_output.xtc'),
        ('traj_output.xtc', 'my_prefix_', 'my_prefix_traj_output.xtc'),
        (None, 'my_prefix', 'my_prefix'),
        (None, 'my_prefix.pdb', 'my_prefix.pdb'),
        ],
    )
def test_add_prefix_to_path(ipath, prefix, expected):
    """Test adding prefix to path."""
    result = io.add_prefix_to_path(ipath, prefix)
    assert isinstance(result, Path)
    assert result.name == expected


@pytest.mark.parametrize(
    'ipath,suffix,expected',
    [
        ('traj_output.xtc', 'my_suffix', 'traj_outputmy_suffix.xtc'),
        ('traj_output.xtc', '_my_suffix', 'traj_output_my_suffix.xtc'),
        ('traj_output.xtc', '_my_suffix.pdb', 'traj_output_my_suffix.pdb'),
        ('traj_output', '_my_suffix', 'traj_output_my_suffix'),
        (None, '_my_suffix', '_my_suffix'),
        (None, '_my_suffix.pdb', '_my_suffix.pdb'),
        ],
    )
def test_add_suffix_to_path(ipath, suffix, expected):
    """Test adding prefix to path."""
    result = io.add_suffix_to_path(ipath, suffix)
    assert isinstance(result, Path)
    assert result.name == expected


@pytest.mark.parametrize(
    'inputs,expected',
    [
        (('traj',), 'traj_frame0.pdb'),
        (('traj', 1), 'traj_frame1.pdb'),
        (('traj', 401), 'traj_frame401.pdb'),
        (('traj', 4, 'dcd'), 'traj_frame4.dcd'),
        (('traj', 4, 'dcd', 4), 'traj_frame0004.dcd'),
        (('traj', 4, 'dcd', 4, '_suffixed.tpr'), 'traj_suffixed.tpr'),
        ],
    )
def test_mk_frame_path(inputs, expected):
    """
    Test mk frame path.

    Notice ``inputs`` is a tuple to unpack.
    """
    result = io.mk_frame_path(*inputs)
    assert isinstance(result, Path)
    assert result.name == expected


@pytest.mark.parametrize(
    'in1,expected',
    [
        ('file', Path('file')),
        (None, Path()),
        ],
    )
def test_get_ipath(in1, expected):
    result = io._get_ipath(in1)
    assert isinstance(result, Path)
    assert result == expected


@pytest.mark.parametrize(
    'top,traj,expected',
    [
        ('top.pdb', 'traj.xtc', 'top.pdb'),
        (Path('top.pdb'), 'traj.xtc', 'top.pdb'),
        ('_top.pdb', 'traj.xtc', 'traj_top.pdb'),
        ('_top', 'traj.xtc', 'traj_top.xtc'),
        ('top_', 'traj.xtc', 'top_traj.xtc'),
        ],
    )
def test_parse_top_output(top, traj, expected):
    """Test parse topology output."""
    result = io.parse_top_output(top, traj_output=traj)
    assert isinstance(result, Path)
    assert result.name == expected


def test_export_data_to_file_1():
    """Test export data to file."""
    xdata = list(range(10))
    ydata = list(range(0, 20, 2))
    header="#this is an header"
    fout = Path(datap, 'exporttest.csv')
    io.export_data_to_file(
        xdata,
        ydata,
        fname=fout,
        header=header,
        delimiter='-',
        fmt='{:.1f}',
        )

    fres = fout.open().read()
    fexp = export_data_expected.open().read()
    assert fres == fexp
    fout.unlink()


def test_export_data_to_file_2():
    """Test export data to file several series."""
    xdata = list(range(10))
    ydata = list(range(0, 20, 2))
    ydata = [ydata, ydata]
    assert isinstance(ydata, list)
    assert len(ydata) == 2
    header="#this is an header"
    fout = Path(datap, 'exporttest.csv')
    io.export_data_to_file(
        xdata,
        *ydata,
        fname=fout,
        header=header,
        delimiter='-',
        fmt='{:.1f}',
        )

    fres = fout.open().read()
    fexp = export_data_expected_2.open().read()
    assert fres == fexp
    fout.unlink()

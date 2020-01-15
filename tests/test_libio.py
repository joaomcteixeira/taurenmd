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
        (
            ['t_1.dcd', 't_11.dcd', 't_2.dcd'],
            ['t_1.dcd', 't_2.dcd', 't_11.dcd'],
            ),
        (['b.dcd', 'c.dcd', 'a.dcd'], ['a.dcd', 'b.dcd', 'c.dcd']),
        ],
    )
def test_sort_numbered_input_1(in1, expected):
    """Test sort numbered inputs."""
    result = io.sort_numbered_input(*in1)
    assert result == expected


@pytest.mark.parametrize(
    'in1,error',
    [
        (['t_1.dcd', 't_11.dcd', 't_2.dcd'], TypeError),
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
    """Test get ipath private function."""
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
    header = "#this is an header\n"
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
    header = "#this is an header\n"
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


@pytest.mark.parametrize(
    'inputs,expected',
    [
        ((10,), list(range(10))),
        ((10, 2, 5), list(range(2, 5))),
        ((10, 2, 50), list(range(2, 10))),
        ((10, None, None, 2), list(range(0, 10, 2))),
        ((10, None, None, None, '1,2,45,65'), [1, 2, 45, 65]),
        ((10, None, None, None, [1, 2, 45, 65]), [1, 2, 45, 65]),
        ((10, None, None, None, ['1', '2', '45', '65']), [1, 2, 45, 65]),
        ((None, None, None, None, ['1', '2', '45', '65']), [1, 2, 45, 65]),
        ],
    )
def test_frame_list(inputs, expected):
    """Test make frame list."""
    result = io.frame_list(*inputs)
    assert result == expected


@pytest.mark.parametrize(
    'lent,flist,error',
    [
        (10, 123, ValueError),
        (None, None, TypeError),
        ]
    )
def test_frame_list_error(lent, flist, error):
    """Test frame list error."""
    with pytest.raises(error):
        io.frame_list(lent, flist=flist)


def test_frame_slice():
    """Test frame_slice."""
    result = io.frame_slice(1, 100, 2)
    assert result == slice(1, 100, 2)


@pytest.mark.parametrize(
    'value,start,stop,step,expected',
    [
        ('1,100,2', None, None, None, slice(1, 100, 2)),
        ('::2', None, None, None, slice(None, None, 2)),
        ('1:100:2', None, None, None, slice(1, 100, 2)),
        ('1:None:2', None, None, None, slice(1, None, 2)),
        (None, 10, None, None, slice(10, None, None)),
        (slice(None), None, None, None, slice(None, None, None)),
        (None, 10, 100, None, slice(10, 100, None)),
        (None, 10, 100, 3, slice(10, 100, 3)),
        ((0, 50, 3), None, None, None, slice(0, 50, 3)),
        ((0, '50', 3), None, None, None, slice(0, 50, 3)),
        ((None, 100, None), None, None, None, slice(None, 100, None)),
        ('10', None, None, None, slice(10, None, None)),
        (10, None, None, None, slice(None, 10, None)),
        (None, None, None, None, slice(None, None, None)),
        (None, {}, None, None, slice(None, None, None)),
        ],
    )
def test_evaluate_to_slice(value, start, stop, step, expected):
    """Test evaluate to slice."""
    result = io.evaluate_to_slice(
        value=value,
        start=start,
        stop=stop,
        step=step,
        )
    assert result == expected


def test_evaluate_to_slice_error():
    """Test evaluate to slice ValueError."""
    with pytest.raises(ValueError):
        io.evaluate_to_slice(value={})

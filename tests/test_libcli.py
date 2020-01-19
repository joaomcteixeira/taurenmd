"""Test libcli."""
import argparse

import pytest

from taurenmd import Path, references
from taurenmd.libs import libcli as lc
from taurenmd.logger import CMDFILE


def test_load_args():
    """Test load_args func."""
    parser = argparse.ArgumentParser()
    # this has to --cov because otherwise argparse catches the
    # pytest commands
    parser.add_argument('--cov', nargs=argparse.REMAINDER)
    result = lc.load_args(parser)
    assert isinstance(result, argparse.Namespace)


def test_maincli():
    """Test maincli."""
    def myfunc(cov, *args, **kwargs):
        return cov
    parser = argparse.ArgumentParser()
    parser.add_argument('--cov', nargs=argparse.REMAINDER)
    result = lc.maincli(parser, myfunc)
    # since we have it, lets play it with and close the circle
    assert result == (
        '--cov-report=term-missing '
        '--cov-append --cov-config=.coveragerc -vv tests'
        ).split()


def test_save_refs():
    """Test save references."""
    references.add("zello world")
    lc.save_references()
    p1 = Path(CMDFILE)
    s = p1.open().readlines()
    assert s[-1].split(':')[1][1:] == 'zello world'
    p1.unlink()


def test_CustomParser():
    """Test Custom Parser."""
    assert issubclass(lc.CustomParser, argparse.ArgumentParser)
    assert hasattr(lc.CustomParser, 'error')


def test_CustomParser_error():
    """Test CP error."""
    with pytest.raises(SystemExit) as error:
        lc.CustomParser().error('my error')
    assert error.type == SystemExit
    assert error.value.code == 2


@pytest.mark.parametrize(
    'key,value,expected',
    [
        ('tilte=', 'my title', 'my title'),
        ('colors=', 'k,b,y,m', ('k', 'b', 'y', 'm')),
        ('xrange=', '1.5,100', (1.5, 100)),
        ('xmin=', '1.5', 1.5),
        ('xmin=', '15', 15),
        ('var1', '', True),
        ],
    )
def test_ParamsToDick(key, value, expected):
    """Test params to dict parsing."""
    p = lc.ParamsToDict(dest='plot', option_strings='--plot')
    namespace = argparse.Namespace
    parser = argparse.ArgumentParser()
    p(
        parser,
        namespace,
        [f'{key}{value}'],
        )
    v = vars(namespace)

    assert 'plotvars' in v
    assert v['plotvars'][key.rstrip('=')] == expected


def test_save_command():
    """Tests only the interface."""
    lc.save_command('testcommandsave', 1, 2, 3, 4)
    Path('testcommandsave').unlink()


def test_add_subparser():
    """Test adds subparser."""

    def mainfunc():
        return 'this is main func'

    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    parserm = argparse.ArgumentParser()
    parserm.add_argument('foo')
    parserm.add_argument('--bar', nargs=1, default=False)
    pseudomodule = argparse.Namespace(
        _name='mycmd',
        __doc__='Module documentation',
        _help='short help',
        ap=parserm,
        main=mainfunc,
        )

    lc.add_subparser(subparser, pseudomodule)
    v = vars(parser.parse_args('mycmd this_foo --bar b'.split()))
    assert v['foo'] == 'this_foo'
    assert v['bar'] == ['b']


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-a degrees', 'degrees'),
        ('--aunit degrees', 'degrees'),
        ('-a radians', 'radians'),
        ('--aunit radians', 'radians'),
        ('', 'degrees'),
        ]
    )
def test_angle_unit(cmd, expected):
    """Test angle unit argument."""
    parser = argparse.ArgumentParser()
    lc.add_angle_unit_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['aunit'] == expected


def test_angle_unit_error():
    """Test sysexit error."""
    parser = argparse.ArgumentParser()
    lc.add_angle_unit_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args('-a'.split())


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-i', True),
        ('--insort', True),
        ('', False),
        ]
    )
def test_insort(cmd, expected):
    """Test angle unit argument."""
    parser = argparse.ArgumentParser()
    lc.add_insort_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['insort'] == expected


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('topology.pdb', 'topology.pdb'),
        ]
    )
def test_topology_input(cmd, expected):
    """Test topology argument."""
    parser = argparse.ArgumentParser()
    lc.add_topology_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['topology'] == expected


def test_topology_input_error():
    """Test topology input error."""
    parser = argparse.ArgumentParser()
    lc.add_topology_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([])


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('traj1.dcd', ['traj1.dcd']),
        ('traj1.dcd traj2.dcd', ['traj1.dcd', 'traj2.dcd']),
        ]
    )
def test_trajectories(cmd, expected):
    """Test trajectories argument."""
    parser = argparse.ArgumentParser()
    lc.add_trajectories_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['trajectories'] == expected


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('traj1.dcd', 'traj1.dcd'),
        ]
    )
def test_trajectory(cmd, expected):
    """Test trajectory argument."""
    parser = argparse.ArgumentParser()
    lc.add_trajectory_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['trajectory'] == expected


def test_trajectory_error():
    """Test error in trajectory input."""
    parser = argparse.ArgumentParser()
    lc.add_trajectory_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args('traj1.dcd traj2.dcd'.split())


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-s 1 -e 100 -p 4', (1, 100, 4)),
        ('-e 100 -p 4', (None, 100, 4)),
        ('-s 100 -p 4', (100, None, 4)),
        ('-s 100 -e 4', (100, 4, None)),
        ]
    )
def test_slice_arg(cmd, expected):
    """Test slice arguments."""
    parser = argparse.ArgumentParser()
    lc.add_slice_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['start'] == expected[0]
    assert v['stop'] == expected[1]
    assert v['step'] == expected[2]


@pytest.mark.parametrize(
    'cmd,expected',
    [
        (['-l', "name CA"], 'name CA'),
        (['--selection', "name CA and segid A"], 'name CA and segid A'),
        ([], 'all'),
        ],
    )
def test_atom_selection_args(cmd, expected):
    """Test atom selection arg."""
    parser = argparse.ArgumentParser()
    lc.add_atom_selection_arg(parser)
    v = vars(parser.parse_args(cmd))
    assert v['selection'] == expected


@pytest.mark.parametrize(
    'cmd',
    [
        ('-l selA selB'),
        ('--selection selB selC selD'),
        ],
    )
def test_atom_selection_error(cmd):
    """Test atom selection error."""
    parser = argparse.ArgumentParser()
    lc.add_atom_selection_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(cmd.split())


@pytest.mark.parametrize(
    'cmd,expected',
    [
        (['-g', "name CA"], ['name CA']),
        (['--selections', "name CA", "segid A"], ['name CA', 'segid A']),
        ([], None),
        ],
    )
def test_atom_selections_args(cmd, expected):
    """Test atom selection arg."""
    parser = argparse.ArgumentParser()
    lc.add_atom_selections_arg(parser)
    v = vars(parser.parse_args(cmd))
    assert v['selections'] == expected


@pytest.mark.parametrize(
    'cmd',
    [
        ('-g'),
        ('--selections'),
        ],
    )
def test_atom_selections_error(cmd):
    """Test atom selection error."""
    parser = argparse.ArgumentParser()
    lc.add_atom_selections_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(cmd.split())


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-t 1', [1]),
        ('--flist 1 2', [1, 2]),
        ],
    )
def test_framelist_arg(cmd, expected):
    """Test frame list."""
    parser = argparse.ArgumentParser()
    lc.add_frame_list_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['flist'] == expected


def test_framelist_error():
    """Test framelist error."""
    parser = argparse.ArgumentParser()
    lc.add_frame_list_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(['-t'])


@pytest.mark.parametrize(
    'cmd,expected',
    [
        (
            ['-z', 'segid A', 'segid B', 'segid C'],
            ['segid A', 'segid B', 'segid C']
            ),
        (
            ['--plane-selection', 'segid A', 'segid B', 'segid C'],
            ['segid A', 'segid B', 'segid C']
            ),
        ],
    )
def test_plane_selection(cmd, expected):
    """Test plane selection."""
    parser = argparse.ArgumentParser()
    lc.add_plane_selection_arg(parser)
    v = vars(parser.parse_args(cmd))
    assert v['plane_selection'] == expected


@pytest.mark.parametrize(
    'cmd',
    [
        (['-z']),
        (['-z', 'A']),
        (['-z', 'A', 'B']),
        (['-z', 'A', 'B', 'D', 'E']),
        ],
    )
def test_plane_selection_error(cmd):
    """Test plane selection error."""
    parser = argparse.ArgumentParser()
    lc.add_plane_selection_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(cmd)


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-r 1', 1),
        ('--ref-frame 2', 2),
        ('', 0),
        ],
    )
def test_reference_frame(cmd, expected):
    """Test reference frame."""
    parser = argparse.ArgumentParser()
    lc.add_reference_frame_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['ref_frame'] == expected


@pytest.mark.parametrize(
    'cmd',
    [
        (['-r', 'A', 'B', 'D', 'E']),
        ],
    )
def test_reference_frame_error(cmd):
    """Test plane selection error."""
    parser = argparse.ArgumentParser()
    lc.add_reference_frame_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(cmd)


@pytest.mark.parametrize(
    'cmd,expected',
    [
        (['--plot', 'title=mytitle'], {'title': 'mytitle'}),
        (
            ['--plot', 'title=mytitle', 'colors=a,b,c'],
            {
                'title': 'mytitle',
                'colors': ('a', 'b', 'c')
                }),
        ],
    )
def test_plot_arg(cmd, expected):
    """Test plot args."""
    parser = argparse.ArgumentParser()
    lc.add_plot_arg(parser)
    v = vars(parser.parse_args(cmd))
    assert v['plot'] is True
    assert 'plotvars' in v
    assert v['plotvars'] == expected


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-o', '_frame0.pdb'),
        ('', False),
        ('--top-output', '_frame0.pdb'),
        ('--top-output frame0_', 'frame0_'),
        ],
    )
def test_topoutput_args(cmd, expected):
    """Test top output."""
    parser = argparse.ArgumentParser()
    lc.add_top_output_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['top_output'] == expected
    

@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-d traj.pdb', 'traj.pdb'),
        ('', 'traj_out.dcd'),
        ('--traj-output traj_out.pdb', 'traj_out.pdb'),
        ('--traj-output newtraj.xtc', 'newtraj.xtc'),
        ],
    )
def test_trajoutput_args(cmd, expected):
    """Test traj output."""
    parser = argparse.ArgumentParser()
    lc.add_traj_output_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['traj_output'] == expected


@pytest.mark.parametrize(
    'cmd',
    [
        ('-t'),
        ('--traj-output'),
        ],
    )
def test_traj_output_arg_error(cmd):
    """Test top output error."""
    parser = argparse.ArgumentParser()
    lc.add_traj_output_arg(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([cmd])


@pytest.mark.parametrize(
    'cmd,expected',
    [
        ('-x', 'results.csv'),
        ('--export', 'results.csv'),
        ('-x data.csv', 'data.csv'),
        ('--export data.csv', 'data.csv'),
        ('', False),
        ],
    )
def test_data_export_arg(cmd, expected):
    """Test data export arg."""
    parser = argparse.ArgumentParser()
    lc.add_data_export_arg(parser)
    v = vars(parser.parse_args(cmd.split()))
    assert v['export'] == expected

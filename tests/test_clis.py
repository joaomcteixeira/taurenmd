"""Test taurenmd client interfaces."""
import argparse
import copy
import functools
import inspect
import sys

import pytest

# import bellow by alphabetical order the cli interfaces you are
# implementing
from taurenmd import __main__ as main_module
from taurenmd import (  # noqa: E133
    cli,
    cli_distances,
    cli_fext,
    cli_imagemol,
    cli_nosol,
    cli_pangle,
    cli_report,
    cli_rmsd,
    cli_rmsf,
    cli_rotations,
    cli_trajedit,
)
from taurenmd.libs import libcli as lc

from . import Path, toptest, trajtest


subclients = [
    # add your new client to this list.
    cli_distances,
    cli_fext,
    cli_imagemol,
    cli_nosol,
    cli_pangle,
    cli_report,
    cli_rmsd,
    cli_rmsf,
    cli_rotations,
    cli_trajedit,
    ]


@pytest.fixture(params=subclients)
def client(request):
    """Loop all individual clients."""
    return request.param


@pytest.fixture(params=[cli] + subclients)
def all_clients(request):
    """Consider individual clientes plus main cli."""
    return request.param


def test_main_module():
    """Test main module entry point."""
    assert hasattr(main_module, 'maincli')


def test_main_module_if():
    """Test cli main entry point."""
    lines = inspect.getsource(main_module).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''


def test_cli__version(all_clients):
    """Test clients have version parameter."""
    with pytest.raises(SystemExit) as err:
        all_clients.ap.parse_args(['-v'])
    assert err.value.code == 0


def test_cli_script_1():
    """Test cli main entry point."""
    assert hasattr(cli, 'maincli')


def test_cli_script_2():
    """Test cli has load_args function."""
    assert hasattr(cli, 'load_args')


def test_cli_script_3():
    """Test cli has _ap function."""
    assert hasattr(cli, '_ap')


def test_cli__ap_returns():
    """Test cli _ap returns ap."""
    assert isinstance(cli._ap(), argparse.ArgumentParser)


def test_cli_script_4():
    """Test cli main if __name__ part."""
    lines = inspect.getsource(cli).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''


def test_cli_load_args():
    """Test cli load_args functions."""
    backup = copy.deepcopy(sys.argv)
    sys.argv = ['taurenmd', 'report', toptest.str(), trajtest.str()]
    result = cli.load_args()
    assert isinstance(result, argparse.Namespace)
    sys.argv = copy.deepcopy(backup)


def test_cli_maincli_execution():
    """Test cli maincli execution."""
    backup = copy.deepcopy(sys.argv)
    sys.argv = ['taurenmd', 'report', toptest.str(), trajtest.str()]
    cli.maincli()
    sys.argv = copy.deepcopy(backup)


def test_cli_maincli_SystemExit():
    """Test cli.py maincli() SystemExit."""
    backup = copy.deepcopy(sys.argv)
    sys.argv = ['taurenmd']
    with pytest.raises(SystemExit):
        cli.maincli()
    sys.argv = copy.deepcopy(backup)


def test_clients_have_maincli(client):
    """Test all modules have maincli."""
    assert hasattr(client, 'maincli')


def test_clientes_have__ap(client):
    """Test all modules have _ap."""
    assert hasattr(client, '_ap')


def test_clients_have__ap_returns(client):
    """Test all modules _ap() return ap."""
    assert isinstance(client._ap(), argparse.ArgumentParser)


def test_clients_have_help(client):
    """Test all modules have help."""
    assert hasattr(client, '_help')


def test_clients_have_main(client):
    """Test all modules have main."""
    assert hasattr(client, 'main')


def test_clients_if_main_code(client):
    """Test if __name__ == '__main__' code lines."""
    lines = inspect.getsource(client).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''


def test_clients_have_ap(client):
    """Test all modules have ap."""
    assert hasattr(client, 'ap')


def test_clients_main_partial_type(client):
    """Test maincli partial type."""
    assert isinstance(client.maincli, functools.partial)


def test_clients_maincli_partial_func(client):
    """Test maincli partial args."""
    assert client.maincli.func == (lc.maincli)


def test_clients_maincli_partial_args(client):
    """Test maincli partial args."""
    assert client.maincli.args == (client.ap, client.main)


@pytest.mark.parametrize(
    'module,name',
    [
        # add here the name of your client, this ensures that changes
        # in the command interface get noticed by the tests.
        # (cli_NAME, 'NAME'),
        (cli_distances, 'dist'),
        (cli_fext, 'fext'),
        (cli_imagemol, 'imagemol'),
        (cli_nosol, 'nosol'),
        (cli_pangle, 'pangle'),
        (cli_report, 'report'),
        (cli_rmsd, 'rmsd'),
        (cli_rmsf, 'rmsf'),
        (cli_rotations, 'rotations'),
        (cli_trajedit, 'trajedit'),
        ],
    )
def test_cli__name(module, name):
    """Test name messages."""
    assert module._name == name


# write a test function that calls the main function for your client.
# and add any additional modification or output confirmation that you
# find appropriate so that test coverage is 100% and you asure the
# correct output is generate. Read through the other client for examples.
# there is a final test to add at the end of the file, do not miss it.

# def test_cli_NAME_1():
#     """Test cli NAME."""
#     cli_NAME.main(
#         #params...
#         )


def test_cli_distances_1():
    """Test cli distances."""
    cli_distances.main(
        toptest,
        [trajtest],
        sel1='resnum 10',
        sel2='resnum 20',
        start=5,
        end=None,
        step=None,
        export='dist.csv',
        plot=True,
        plotvars={'filename': 'dist.pdf'},
        )
    Path('dist.csv').unlink()
    Path('dist.pdf').unlink()


def test_cli_distances_2():
    """Test cli distances."""
    cli_distances.main(
        toptest,
        [trajtest],
        sel1='resnum 10',
        sel2='resnum 20',
        start=None,
        end=5,
        step=None,
        export=False,
        plot=False,
        )
    assert not Path('dist.csv').exists()
    assert not Path('dist.pdf').exists()


def test_cli_fext_1():
    """Test client frame extract."""
    cli_fext.main(
        toptest,
        [trajtest],
        start=4,
        stop=6,
        step=None,
        prefix='timestep_',
        ext='pdb',
        selection='name CA',
        )
    p1 = Path('timestep_04.pdb')
    p2 = Path('timestep_05.pdb')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_imagemol():
    """Test imagemol."""
    cli_imagemol.main(
        toptest,
        trajtest,
        traj_output='testimage.dcd',
        top_output='testtop.pdb',
        protocol=1,
        )
    p1 = Path('testimage.dcd')
    p2 = Path('testtop.pdb')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_imagemol_2():
    """Test imagemol."""
    cli_imagemol.main(
        toptest,
        trajtest,
        traj_output='testimage.dcd',
        top_output=False,
        protocol=2,
        )
    p1 = Path('testimage.dcd')
    assert p1.exists()
    p1.unlink()


def test_cli_nosol_1():
    """Test cli nosol."""
    cli_nosol.main(
        toptest,
        trajtest,
        selection='name CA',
        maintain=None,
        top_output='_f0.pdb',
        traj_output='nosol.dcd',
        )
    p1 = Path('nosol.dcd')
    p2 = Path('nosol_f0.pdb')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_nosol_2():
    """Test cli nosol."""
    cli_nosol.main(
        toptest,
        trajtest,
        selection=None,
        maintain=['CL'],
        top_output='_f1.pdb',
        traj_output='nosol.dcd',
        )
    p1 = Path('nosol.dcd')
    p2 = Path('nosol_f1.pdb')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_nosol_3():
    """Test cli nosol."""
    cli_nosol.main(
        toptest,
        trajtest,
        selection=None,
        maintain=['CL', 'NA'],
        top_output=None,
        traj_output='nosol.dcd',
        )
    p1 = Path('nosol.dcd')
    p2 = Path('nosol_f1.pdb')
    assert p1.exists()
    assert not p2.exists()
    p1.unlink()


def test_cli_pangle_1():
    """Test client plane angle."""
    cli_pangle.main(
        toptest,
        [trajtest],
        ['resnum 10', 'resnum 20', 'resnum 30'],
        aunit='degrees',
        export='angles.csv',
        plot=True,
        plotvars={'filename': 'angles.pdf'},
        )
    p1 = Path('angles.csv')
    p2 = Path('angles.pdf')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_pangle_2():
    """Test client plane angle."""
    cli_pangle.main(
        toptest,
        [trajtest],
        ['resnum 10', 'resnum 20', 'resnum 30'],
        aunit='radians',
        export=False,
        plot=False,
        )


def test_cli_report_1():
    """Test cli report."""
    cli_report.main(toptest, [trajtest])


def test_cli_rmsd_1():
    """Test cli rmsd."""
    cli_rmsd.main(
        toptest,
        [trajtest],
        selections=['name CA'],
        )


def test_cli_rmsd_2():
    """Test cli rmsd."""
    cli_rmsd.main(
        toptest,
        [trajtest],
        selections=None,
        export='rmsds.csv',
        plot=True,
        plotvars={'filename': 'rmsds.pdf'},
        )
    p1 = Path('rmsds.csv')
    p2 = Path('rmsds.pdf')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_rmsf_1():
    """Test cli rmsf."""
    cli_rmsf.main(
        toptest,
        [trajtest],
        selections=['name CA'],
        )


def test_cli_rmsf_2():
    """Test cli rmsd."""
    cli_rmsf.main(
        toptest,
        [trajtest],
        selections=None,
        export='rmsfs.csv',
        plot=True,
        plotvars={'filename': 'rmsfs.pdf'},
        )
    p1 = Path('rmsfs.csv')
    p2 = Path('rmsfs.pdf')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_rmsf_error1():
    """Test cli rmsd."""
    with pytest.raises(TypeError):
        cli_rmsf.main(
            toptest,
            [trajtest],
            selections='name CA',
            )


def test_cli_rmsf_error2():
    """Test cli rmsd."""
    with pytest.raises(TypeError):
        cli_rmsf.main(
            toptest,
            [trajtest],
            selections=[],
            )


def test_cli_rotations_1():
    """Test cli rotations."""
    cli_rotations.main(
        toptest,
        [trajtest],
        plane_selection=['resnum 10', 'resnum 20', 'resnum 30'],
        aunit='radians',
        export='rots.csv',
        )

    p1 = Path('roll_angles_rots.csv')
    p2 = Path('pitch_angles_rots.csv')
    p3 = Path('yaw_angles_rots.csv')
    assert p1.exists()
    assert p2.exists()
    assert p3.exists()
    p1.unlink()
    p2.unlink()
    p3.unlink()


def test_cli_rotations_2():
    """Test cli rotations."""
    cli_rotations.main(
        toptest,
        [trajtest],
        plane_selection=['resnum 10', 'resnum 20', 'resnum 30'],
        )


def test_cli_trajedit_1():
    """Test trajedit 1."""
    cli_trajedit.main(
        toptest,
        [trajtest],
        insort=True,
        selection='name CA',
        traj_output='traj.dcd',
        top_output='topooo.pdb',
        unwrap=True,
        align='name CA',
        )
    p1 = Path('traj.dcd')
    p2 = Path('topooo.pdb')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()


def test_cli_trajedit_2():
    """Test trajedit 2."""
    cli_trajedit.main(
        toptest,
        [trajtest],
        )
    p1 = Path('traj_out.dcd')
    assert p1.exists()
    p1.unlink()


def test_cli_trajedit_3():
    """Test trajedit 3."""
    cli_trajedit.main(
        toptest,
        [trajtest],
        top_output='topooo.pdb',
        )
    p2 = Path('topooo.pdb')
    p1 = Path('traj_out.dcd')
    assert p1.exists()
    assert p2.exists()
    p1.unlink()
    p2.unlink()

# test argparse


@pytest.mark.parametrize(
    'module,cmd',
    [
        # add command examples to this list
        (cli, 'fext top.pdb traj1.xtc traj2.xtc'),
        (cli, 'dist top.pdb traj1.xtc traj2.xtc -l1 selA -l2 selB -s 1 -e 2 -p 3 -x --plot'),  # noqa: E501
        (cli, 'imagemol top.pdb traj.xtc -d tout.xtc -o'),
        (cli, 'imagemol top.pdb traj.xtc -d tout.xtc -o out.pdb'),
        (cli, 'pangle top.pdb traj_1.xtc traj_2.xtc -x --plot title=title -s 1 -e 2 -p 3 -z selA selB selC -a radians'),  # noqa: E501
        (cli, 'trajedit top.pdb traj1.xtc traj2.xtc -i -l segA -s 1 -e 10 -p 2 -d tout.xtc -o'),  # noqa: E501
        (cli, 'nosol top.pdb traj.xtc -d tout.xtc -o out.pdb -m NA'),
        (cli, 'rmsd top.pdb traj1.xtc traj2.xtc -g segA segB -s 1 -e 10 -p 2 -r 10 -x data.csv --plot'),  # noqa: E501
        (cli, 'rmsd top.pdb traj1.xtc traj2.xtc -g segA segB -s 1 -e 10 -p 2 -r 10 -x data.csv --plot title=1'),  # noqa: E501
        (cli, 'rmsf top.pdb traj1.xtc traj2.xtc -g segA segB -s 1 -e 10 -p 2 -x data.csv --plot title=1'),  # noqa: E501
        (cli, 'rotations top.pdb traj1.xtc traj2.xtc -z segA segB segC -s 1 -e 10 -p 2 -a radians -x data.csv'),  # noqa: E501
        (cli, 'report top.pdb traj1.xtc traj2.xtc'),
        ],
    )
def test_ap_interface(module, cmd):
    """Test argparse interface of command lines."""
    module.ap.parse_args(cmd.split())

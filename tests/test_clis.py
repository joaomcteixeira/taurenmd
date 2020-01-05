"""Test taurenmd client interfaces."""
import argparse
import functools
import inspect

import pytest

from taurenmd import (
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

from . import *

@pytest.fixture(
    params=[
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
    )
def client(request):
    return request.param


def test_has_maincli(client):
    """Test all modules have maincli."""
    assert hasattr(client, 'maincli')


def test_has__ap(client):
    """Test all modules have _ap."""
    assert hasattr(client, '_ap')


def test_has__ap_returns(client):
    """Test all modules have _ap."""
    assert isinstance(client._ap(), argparse.ArgumentParser)


def test_has_help(client):
    """Test all modules have help."""
    assert hasattr(client, '_help')


def test_has_main(client):
    """Test all modules have main."""
    assert hasattr(client, 'main')


def test_has_ap(client):
    """Test all modules have ap."""
    assert hasattr(client, 'ap')


def test_main_partial_type(client):
    """Test maincli partial type."""
    assert isinstance(client.maincli, functools.partial)


def test_maincli_partial_func(client):
    """Test maincli partial args."""
    assert client.maincli.func == (lc.maincli)


def test_maincli_partial_args(client):
    """Test maincli partial args."""
    assert client.maincli.args == (client.ap, client.main)



@pytest.mark.parametrize(
    'module,name',
    [
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
def test_help(module, name):
    """Test name messages."""
    assert module._name == name


def test_if_main_code(client):
    """Test if __name__ == '__main__' code lines."""
    lines = inspect.getsource(client).split('\n')
    assert lines[-3] == "if __name__ == '__main__':"
    assert lines[-2] == "    maincli()"
    assert lines[-1] == ''


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


def test_cli_nosol_2():
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

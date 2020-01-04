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

"""Test libcli."""
import argparse

import pytest

from taurenmd import Path
from taurenmd.libs import libcli as lc


def test_references():
    assert isinstance(lc.ref_mdt, str)
    assert isinstance(lc.ref_mda, str)
    assert isinstance(lc.ref_mda_selection, str)
    assert isinstance(lc.ref_mda_unwrap, str)
    assert isinstance(lc.ref_mda_alignto, str)
    assert isinstance(lc.ref_plottemplates_param, str)
    assert isinstance(lc.ref_plottemplates_labeldots, str)
    assert isinstance(lc.ref_pyquaternion, str)


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
    assert result == '--cov-report=term-missing -vv tests'.split()


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

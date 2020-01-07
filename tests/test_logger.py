"""Test logger module."""

import taurenmd.logger as tl


def test_has_T():
    """Test has T."""
    assert hasattr(tl, 'T')


def test_has_S():
    """Test has S."""
    assert hasattr(tl, 'S')


def test_TitleLog():
    """Test TitleLog formatter."""
    t = tl.TitleLog('section title: {}', 'chapter 1')
    assert '\n* Section Title: chapter 1 ...' == str(t)


def test_SubLog():
    """Test SubLog formatter."""
    s = tl.SubLog('some args: {}, {}', '1', '2')
    assert '    some args: 1, 2' == str(s)


def test_SubLog_2():
    """Test SubLog formatter."""
    s = tl.SubLog('some args: {}, {}', '1', '2', indent=1, spacer='*')
    assert '****some args: 1, 2' == str(s)

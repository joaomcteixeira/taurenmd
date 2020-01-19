"""Test libplot."""

from taurenmd import Path
from taurenmd.libs import libplot as lplt


def test_param():
    """Test param plot."""
    p = Path('param.pdf')
    lplt.param(
        list(range(10)),
        list(range(10)),
        filename='param.pdf',
        )
    assert p.exists()
    p.unlink()


def test_label_dots():
    """Test label_dots plot."""
    p = Path('label_dots.pdf')
    lplt.param(
        ['a', 'b', 'c'],
        list(range(3)),
        filename='label_dots.pdf',
        )
    assert p.exists()
    p.unlink()

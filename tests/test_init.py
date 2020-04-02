"""Test things in init."""

import pytest

from taurenmd import _controlled_exit


def test_controlled_exit():
    """Test controlled exit function."""
    with pytest.raises(SystemExit) as exit:
        _controlled_exit()
    assert exit.type == SystemExit
    assert exit.value.code == 126

"""Test core definitions."""
import taurenmd.core as tcore


def test_references():
    """Test scientific references to third party libraries."""
    assert isinstance(tcore.ref_mdt, str)
    assert isinstance(tcore.ref_mda, str)
    assert isinstance(tcore.ref_mda_selection, str)
    assert isinstance(tcore.ref_mda_unwrap, str)
    assert isinstance(tcore.ref_mda_alignto, str)
    assert isinstance(tcore.ref_plottemplates_param, str)
    assert isinstance(tcore.ref_plottemplates_labeldots, str)
    assert isinstance(tcore.ref_pyquaternion, str)
    assert isinstance(tcore.ref_numpy, str)

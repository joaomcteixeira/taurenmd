from taurenmd.libs import libutil


class TestFrameSlice:

    def test_frame_slice_all(self):
        """Test 'all'."""
        result = libutil.frame_slice('all')
        assert result == slice(None, None, None)

    def test_frame_slice_none(self):
        """Test None."""
        result = libutil.frame_slice(None)
        assert result == slice(None, None, None)

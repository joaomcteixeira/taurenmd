from taurenmd.libs import libutil


class TestFrameSlice:

    def test_frame_(self):
        """Test None."""
        result = libutil.frame_slice()
        assert result == slice(None, None, None)

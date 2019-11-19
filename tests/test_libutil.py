from taurenmd.libs import libutil


class TestFrameSlice:

    def test_frame_(self):
        """Test None."""
        result = libutil.frame_slice()
        assert result == slice(None, None, None)

    def test_tuple1(self):
        result = libutil.frame_slice()

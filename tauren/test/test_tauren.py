import pytest

from tauren.tauren import TaurenTraj

def test_string_to_slice_1():
    """single digit '1'"""
    
    frames500 = list(range(1, 501))
    
    s1 = slice(0, 1, 1)
    s2 = TaurenTraj._get_frame_slicer_from_string("1")
    
    assert frames500[s1] == frames500[s2]

def test_string_to_slice_2():
    """range '10:250'"""
    
    frames500 = list(range(1, 501))
    
    s1 = slice(9, 250, 1)
    s2 = TaurenTraj._get_frame_slicer_from_string("10:250")
    
    assert frames500[s1] == frames500[s2]

def test_string_to_slice_3():
    """range '200:'"""
    
    frames500 = list(range(1, 501))
    
    s1 = slice(199, 500, 1)
    s2 = TaurenTraj._get_frame_slicer_from_string("200:")
    
    assert frames500[s1] == frames500[s2]

def test_string_to_slice_4():
    """range ':300'"""
    
    frames500 = list(range(1, 501))
    
    s1 = slice(0, 300, 1)
    s2 = TaurenTraj._get_frame_slicer_from_string(":300")
    
    assert frames500[s1] == frames[s2]


def test_string_to_slice_5():
    """range with step 1:500:100"""
    
    frames500 = list(range(1, 501))
    
    s1 = slice(0, 500, 100)
    s2 = TaurenTraj._get_frame_slicer_from_string("1:500:100")
    
    assert frames500[s1] == frames500[s2]
    

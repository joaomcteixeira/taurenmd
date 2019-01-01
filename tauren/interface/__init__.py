import sys
import os

sys.path.append(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))

__all__ = [
    "interface",
    ]

from tauren.interface import interface as interface

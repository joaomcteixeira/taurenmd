"""General utilities."""
import os
import re
from itertools import zip_longest


def make_list(*items):
    """Make a list out of objects."""
    new = []

    for item in items:
        if isinstance(item, str):
            new.append(item)
        else:
            new.extend(item)
    return new


def split_time_unit(s):
    """
    Split time and units.

    Follows the regex: https://regex101.com/r/LZAbil/2

    Returns
    -------
    tuple (float, str)
        Value as float, units as str.

    Raises
    ------
    IndexError
        Tuple could not be found. This happens when a number is not
        present in the start of the string.
    """
    type_regex = re.compile(r'^(\-?\d+\.?\d*|\-?\.\d+|\-?\.?\d+[eE]\-?\d+|-?\d+\.?\d*[eE]\d+)($|[a-z]*$)')  # noqa: E501
    try:
        value, unit = type_regex.findall(s)[0]
    except IndexError as err:
        raise ValueError(f'Time string not appropriate: {s!r}') from err
    return float(value), unit


def make_csv_lines_in_interleaved_manner(data, labels):
    """
    Create CSV lines from a list of lists in interleaved manner.

    Labels and data do not need to be of the same size. Resulting
    sizes will be completed with empty strings for those series with
    less values.

    Parameters
    ----------
    data : list of lists of numbers
        The data series.

    labels : list of lists of str
        The labels refering to the data series.

    Returns
    -------
    str
        A string merging the lines.
    """
    # require
    assert all(len(i) == len(j) for i, j in zip(data, labels)), \
        'data and labels size do not match'

    dataz = zip_longest(*data, fillvalue='')
    labelz = zip_longest(*labels, fillvalue='')

    lines = []
    for di, li in zip(dataz, labelz):
        _lines = (f'{_l},{_d}' for _d, _l in zip(di, li) if _d)
        lines.append(','.join(_lines))
    return os.linesep.join(lines)

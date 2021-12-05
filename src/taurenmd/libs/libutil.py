"""General utilities."""
import re


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

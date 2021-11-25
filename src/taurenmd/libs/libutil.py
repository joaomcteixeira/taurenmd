"""General utilities."""


def make_list(*items):
    """
    Makes a list out of objects.
    """
    new = []

    for item in items:
        if isinstance(item, str):
            new.append(item)
        else:
            new.extend(item)
    return new

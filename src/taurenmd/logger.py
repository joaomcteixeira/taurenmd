"""Manages operations with logging."""


class TitleLog:
    """Format string to title."""
    def __init__(self, msg):
        self.msg = msg.title()
    
    def __str__(self):
        return '\n* {} ...'.format(self.msg)


class SubLog:
    """
    Format string to bullet point like structure.
    
    This format performs nicely under the `TitleLog` formatting.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return '    {}'.format(self.msg)


T = TitleLog
S = SubLog

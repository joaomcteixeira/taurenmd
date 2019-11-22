"""Manages operations with logging."""


class TitleLog:
    """Format string to title."""
    def __init__(self, msg, *args):
        self.msg = msg.title()
        self.args = args
    
    def __str__(self):
        output = self.msg.format(*self.args)
        return '\n* {} ...'.format(output)


class SubLog:
    """
    Format string to bullet point like structure.
    
    This format performs nicely under the `TitleLog` formatting.
    """
    
    def __init__(self, msg, *args, indent=1, spacer=' '):
        self.msg = msg
        self.args = args
        self.indent = spacer * 4 * indent
    
    def __str__(self):
        output = self.msg.format(*self.args)
        return '{}{}'.format(self.indent, output)


T = TitleLog
S = SubLog

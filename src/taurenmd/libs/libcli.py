"""Common operations for client interfaces."""
import argparse
import sys
from datetime import datetime


# https://stackoverflow.com/questions/4042452
class CustomParser(argparse.ArgumentParser):
    """Custom Parser class."""
    
    def error(self, message):
        """Present error message."""
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


class ParamsToDict(argparse.Action):
    """Convert command-line parameters to dict."""
    
    def __call__(self, parser, namespace, values, option_string=None):
        """Executes."""

        bool_value = {
            'true': True,
            'false': False,
            }

        param_dict = {}
        for kv in values:
            # print(param_dict, kv)
            try:
                k, v = kv.split('=')
            except ValueError:
                param_dict[kv] = True
            else:
                if ',' in v:
                    vs = v.split(',')
                    try:
                        param_dict[k] = tuple(float(i) for i in vs)
                    except (ValueError, TypeError):
                        param_dict[k] = tuple(i for i in vs)

                else:
                    try:
                        param_dict[k] = float(v)
                    except (ValueError, TypeError):  # is string or list
                        param_dict[k] = bool_value.get(v.lower(), v)
            
        setattr(namespace, self.dest, param_dict)


def save_command(fname, *args):
    """Append the input command to a file."""
    with open(fname, 'a') as fh:
        fh.write(
            '[{}] {}\n'.format(
                datetime.now().strftime("%d/%m/%Y, %H:%M:%S"),
                ' '.join(args),
                )
            )

def add_subparser(parser, module):
    """
    Adds a subcommand to a parser.
    """
    new_ap = parser.add_parser(
        module._name,
        description=module.__doc__,
        help=module._help,
        parents=[module.ap],
        add_help=False,
        )

    new_ap.set_defaults(func=module.main)


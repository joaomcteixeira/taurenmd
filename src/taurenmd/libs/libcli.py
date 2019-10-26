"""Common operations for client interfaces."""
import argparse
import sys


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
        param_dict = {}
        for kv in values:
            # print(param_dict, kv)
            try:
                k, v = kv.split('=')
            except ValueError:
                param_dict[kv] = True
            else:
                if ',' in v:
                    v = v.split(',')
                    
                param_dict[k] = v
        
        setattr(namespace, self.dest, param_dict)

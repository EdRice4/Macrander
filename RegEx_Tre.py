# {{{ Imports
from os import getcwd, listdir
import re
from numpy import loadtxt
import argparse
# }}}


# {{{ FileIO
class FileIO(object):

    """ {{{ Docstrings
    A class in which file input/output is handled.

    Namely:

        1.) A new tree is written to a file.

    }}} """

    def write_tre_file(self, sub_tre_file):

        """ {{{ Docstrings
        Writes a new tree to a file.
        }}} """

        with open(self.new_tre_file, 'w') as ntf:
            ntf.write(sub_tre_file)


# {{{ RegEx
class RegEx(FileIO):

    """ {{{ Docstrings
    A class in which all regular expression functionality is stored.
    }}} """

    def substitute(self):
        replace = dict((re.escape(k), v) for k, v in self.dfile)
        pattern = re.compile('|'.join(replace.iterkeys()))
        sub_tre_file = pattern.sub(lambda m: replace[re.escape(m.group(0))],
                                   self.tre)
        return sub_tre_file


class IterRegistry(type):

    """A class allowing for interation over instances of TreFile class."""

    def __iter__(cls):
        return iter(cls.registry)


class TreFile(RegEx):

    """A class in which all the necessary parameters corresponding to each
       respective tre file are stored."""

    __metaclass__ = IterRegistry
    registry = []

    def __init__(self, tre_file, dict_file):
        with open(tre_file, 'r') as tre:
            self.tre = tre.read()
        with open(dict_file, 'r') as dfile:
            self.dfile = loadtxt(dfile, dtype=str, delimiter='\t')
        self.new_tre_file = tre_file.replace('.tre', '_subbed.tre')
        self.registry.append(self)

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('tre_file', type=str, help=('Name of tre file to be '
                                                    'substituted. Specify '
                                                    'None if running in '
                                                    'batch mode.')
                        )
arg_parser.add_argument('dict_file', type=str, help=('Name of dictionary '
                                                     'file. Specify None if '
                                                     'running in batch mode.')
                        )
arg_parser.add_argument('-b', '--batch', help=('Run script in batch mode. '
                                               'i.e. perform substition on '
                                               'all tre files directory with '
                                               'their respective dictionary '
                                               'files.'),
                        action='store_true'
                        )
args = arg_parser.parse_args()

if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    tre_files = sorted(filter(lambda x: '.tre' in x, fid))
    dict_files = sorted(filter(lambda x: '.txt' in x, fid))
    for i, j in zip(tre_files, dict_files):
        TreFile(i, j)
else:
    TreFile(args.tre_file, args.dict_file)

for tre in TreFile:
    sub_tre_file = tre.substitute()
    tre.write_tre_file(sub_tre_file)

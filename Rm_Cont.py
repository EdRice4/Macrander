# {{{ Imports
from os import getcwd, listdir
import re
import argparse
# }}}


# {{{ RegEx
class RegEx(object):

    """ {{{ Docstrings
    A class in which regular expression functionality is stored.
    }}} """

    # {{{ strip_and_escape
    def strip_and_escape(self):
        rm_strip_esc = map(lambda x: re.escape(x.strip()), self.rm)
        return rm_strip_esc
    # }}}

    # {{{ compile_rm
    def compile_rm(self):
        pattern = self.strip_and_escape()
        pattern = re.compile('|'.join(pattern))
        return pattern
    # }}}

    # {{{ remove_cont
    def remove_cont(self):
        pattern = self.compile_rm()
        for i in self.fas:
            if bool(pattern.match(i)):
                self.fas.remove(i)
    # }}}
# }}}


# TODO: What sort of files am I working with here?
# {{{ FileIO
class FileIO(RegEx):

    """ {{{ Docstrings
    A class in which file input/output is stored.
    }}} """

    # {{{ write_new_output
    def write_new_output(self):
        new_fasta_seq = '>'.join(self.fas)
        with open(self.new_fasta, 'w') as fas:
            fas.write('>')
            fas.write(new_fasta_seq)
    # }}}
# }}}


# {{{ IterRegistry
class IterRegistry(type):

    """ {{{ Docstrings
    A metaclass allowing for iterations over instances of the FastaFile
    class.
    }}} """

    # {{{ __iter__
    def __iter__(cls):
        return iter(cls.registry)
    # }}}
# }}}


# {{{ FastaFile
class FastaFile(FileIO):

    """ {{{ Docstrings
    A class in which all the necessary parameters corresponding to each
    respective fasta file are stored.
    }}} """

    # {{{ __metaclass__
    __metaclass__ = IterRegistry
    registry = []
    # }}}

    # {{{ __init__
    def __init__(self, fasta_file, rem_file):
        with open(fasta_file, 'r') as fas:
            fas = fas.read()
            self.fas = (fas[1:]).split('>')
        # TODO: What sort of files am I working with here?
        self.decont = fasta_file.replace('.fasta', '_new.fasta')
        with open(rem_file, 'r') as rem:
            self.rm = rem.readlines()
        self.registry.append(self)
    # }}}
# }}}

# {{{ ArgParser
arg_parser = argparse.ArgumentParser(
        prog='Remove_Contaminants.py',
        description=(
            'Intended to be used in sequence with Fasta2Phylip.py. Given an '
            'output file \'contaminated\' with unique IDs, will substitute '
            'these IDs with the original IDs from the fasta file.'
            ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
arg_parser.add_argument(
        'cont_file', type=str,
        help=(
                'Name of output file containing \'contaminating\' unique '
                'sequences which are to be substituted with original IDs.'
        ),
        default=None
        )
arg_parser.add_argument(
        'dict_file', type=str,
        help=(
                'Name of dictionary file (from Fasta2Phylip.py) containing '
                'dictionary of unique and original identities.'
                ),
        default=None
        )
arg_parser.add_argument(
        '-b', '--batch',
        help=(
                'Run script in batch mode. i.e. perform substitution with all '
                'respective \'contaminated\' files and their corresponding '
                'dictionary files.'
                ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Batch
if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    # TODO: What sort of files am I working with here?
    fasta_files = sorted([x for x in fid if '.fasta' in x])
    rem_files = sorted([x for x in fid if '.txt' in x])
    for i, j in zip(fasta_files, rem_files):
        FastaFile(i, j)
else:
    FastaFile(args.fasta_file, args.rem_file)
# }}}


# {{{ Run
for fas in FastaFile:
    fas.remove_cont()
    fas.write_new_fasta()
# }}}

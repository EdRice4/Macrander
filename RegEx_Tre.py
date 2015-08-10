# {{{ Header
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code written by: Edwin Rice
# email: edwinricethe4th@gmail.com
# phone: +1 (513) 426-4187
# github: https://github.com/EdRice4
#
# Initiation date|time: 08/10/2015|12:59:38
# }}}


# {{{ Imports
from os import getcwd, listdir
import re
import argparse
from sys import exit
# }}}


# {{{ FileIO
class FileIO(object):

    """ {{{ Docstrings

    A class in which file input/output is handled.

    }}} """

    # {{{ read_tree_file
    def read_tree_file(self, tree_file):

        """ {{{ Docstrings

        Reads a tree file into a string.

        }}} """

        # Open file in read mode
        with open(tree_file, 'r') as tree:
            # Read in as one, long string
            tree = tree.read()
        return tree
    # }}}

    # {{{ read_dict_file
    def read_dict_file(self, dict_file):

        """ {{{ Docstrings

        Reads a dictionary file into a dictionary with the format:

            dictionary = {

                'Substitute ID' : 'Original ID'

                }

        """

        # Initiate empty dictionary to store IDs
        IDs = {}
        # Open file in read mode
        with open(dict_file, 'r') as dictionary:
            dictionary = dict_file.readlines()
        for line in dictionary:
            # Partition line into two values, based on location of tab "\t"
            # character
            line = line.split('\t')
            # Set values for readability
            Substitute_ID = line[0]
            Original_ID = line[1]
            # Add to dictionary
            IDs[Substitute_ID] = Original_ID
        return IDs
    # }}}

    # {{{ write_tree_file
    def write_tree_file(self, sub_tre_file):

        """ {{{ Docstrings

        Writes a new tree to a file.

        }}} """

        with open(self.new_tre_file, 'w') as new_tree_file:
            new_tree_file.write(sub_tre_file)
    # }}}


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
# }}}


class TreFile(RegEx):

    """ {{{ Docstrings

    A class in which all the necessary parameters corresponding to each
    respective tre file are stored.

    }}} """

    def __init__(self, tree_file, dict_file):
        self._tree = self.read_tree_file(tree_file)
        self._dict = self.read_dict_file(dict_file)
        self._sub_tree_file = self._tree.replace('.tre', '_subbed.tre')
        self._sub_tree = self.substitute()
        self.write_tree_file()

# {{{ ArgParse
arg_parser = argparse.ArgumentParser(
        prog='RegEx_Tre.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            'A script intended to be used as the final component in the '
            'workflow of F2P.py. F2P.py produces three output files: 1.) '
            'A phylip file containing "original" IDs and sequences; 2.) '
            'A phylip file containing "substitue" IDs and seqeunces (i.e. '
            '~10 digit integers); 3.) A dictionary file containing each '
            'original ID and its repective substitute ID. This script accepts '
            'a tree ouput file (or any output file) containing the substitute '
            'IDs and replaces each, utilizing regular expressions, with its '
            'respective original ID as defined by the dictionary file.'
            ),
        epilog=(
            'Note: When running in batch mode, all tre files you wish to run '
            'should contain the string ".tre" in the file name while all '
            'dictionary files you wish to run should contain the string '
            '".txt" in the file name. Furthermore, each tree file and its '
            'respective dictionary file should be named analogously. For '
            'instance, "Grammostola_rosea.tre" and "Grammostola_rosea.txt."'
            )
        )
arg_parser.add_argument(
        'tre_file', type=str, help='Name of tre file to be substituted.',
        default=None
        )
arg_parser.add_argument(
        'dict_file', type=str,
        help=(
                'Name of corresponding dictionary file containing original '
                'IDs and respective substitute IDs.'
                ),
        default=None
        )
arg_parser.add_argument(
        '-b', '--batch',
        help=(
                'Run script in batch mode. i.e. perform substition on all tre '
                'files withint directory and their respective dictionary '
                'files.'
                ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Batch
# Instantiate instance of TreFile class with each tre file within directory and
# its corresponding dictionary file
if args.batch:
    # Get current working directory
    cwd = getcwd()
    # Get all files in cwd
    fid = listdir(cwd)
    # Filter out tre files
    # NOTE: All tree files you wish to run should contain the string ".tre"
    # and all files within the directory containing this string will be run
    tre_files = sorted(filter(lambda x: '.tre' in x, fid))
    # Filter out dictionary files
    # NOTE: All dictionary files you wish to run should contain the string
    # ".txt" and all files within the directory containing this string will be
    # run
    dict_files = sorted(filter(lambda x: '.txt' in x, fid))
    # Instantiate instance of TreFile class for each tree file and its
    # respective dictionary file found
    for i, j in zip(tre_files, dict_files):
        TreFile(i, j)
# Else, utilize user-specifed strings to instantiate single instance of PepFile
# class
else:
    # Exit if no proper tree/dictionary file specified.
    if not args.tree_file or args.dict_file:
        exit(
                'You did not specify a tree/dictionary file to run nor did '
                'you run the script in batch mode. Try again.'
                )
    TreFile(args.tre_file, args.dict_file)
# }}}

for tre in TreFile:
    sub_tre_file = tre.substitute()
    tre.write_tre_file(sub_tre_file)

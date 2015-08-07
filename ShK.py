#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code written by: Edwin Rice
# email: edwinricethe4th@gmail.com
# phone: +1 (513) 426-4187
# github: https://github.com/EdRice4
#
# Initiation date|time: 08/05/2015|18:16:26

# {{{ Imports
from os import getcwd, listdir
from re import match
import argparse
# }}}


# {{{ PepFile
class PepFile(object):

    """ {{{ Docstrings

    Class in which pep file functionality is stored.

    }}} """

    # {{{ read_pep_into_dict
    def read_pep_into_dict(self, pep_file):

        """ {{{ Docstrings

        Reads the pep file into a dictionary in the format of:

            dictionary = {

                    'Pep_ID' : 'Pep_seq'

                    }

        }}} """

        # Read in pep file
        with open(pep_file, 'r') as pep:
                pep = pep.readlines()
        # Filter out IDs
        ids = [line for line in pep if line[0] == '>']
        # Filter out peptide sequences
        seqs = [line for line in pep if line[0] != '>']
        # Concatenate into dictionary
        pep_dict = dict(zip(ids, seqs))
        return pep_dict
    # }}}
# }}}


# {{{ IterRegistry
class IterRegistry(type):

    """ {{{ Docstrings

    A metaclass allowing for iterations over the Data class.

    }}} """

    # {{{ __iter__
    def __iter__(cls):
        return iter(cls.registry)
    # }}}
# }}}


# TODO:Where is this inheriting from?
# {{{ Data
class Data(object):

    """ {{{ Docstrings

    A class in which all pertinent data for each repsective peptide file
    is stored.

    }}} """

    # {{{ __metaclass__
    __metaclass__ = IterRegistry
    registry = []
    # }}}

    # {{{ __init__
    def __init__(self, pep_file):
        self._pep_dict = self.read_pep_into_file(pep_file)
        self.registry.append(self)
    # }}}
# }}}


# {{{ Batch
if args.batch:
    # Get current working directory
    cwd = getcwd()
    # Get all files in cwd
    files = listdir(cwd)
    # Filter out pep files
    # NOTE: All files you wish to run should contain the string '.pep' in the
    # name.
    pep_files = [x for x in files if '.pep' in x]
    # Instantiate intances of Data calss for all pep files found
    for i in pep_files:
        Data(i)
else:
    Data(args.pep_file)
# }}}


# {{{ ArgParse
arg_parser = argparse.ArgumentParser(
        prog='ShK.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
                'Given a list of peptide sequences, putatively containing '
                'the ShK domain (as described in Rangaraju et al., "Potassium'
                'Channel Modulation by a Toxin Domain in Matrix '
                'Metalloprotease 23." DOI: 10.1074/jbc.M109.071266 and in the'
                'Simple Modular Architecture Research Tool database:'
                'http://smart.embl.de/smart/do_annotation.pl?DOMAIN=SM00254),'
                'will robustly filter out those peptide sequences and '
                'characterize the nature of each sequence.'
                ),
        epilog=(
                'Note: When running in batch mode, all files you wish to run '
                'should contain the string ".pep" in the file name.'
                )
        )
arg_parser.add_argument(
        '-pep', type=str, help='Name of pep file you wish to run.',
        default=None
        )
arg_parser.add_argument(
        '-b', '--batch',
        help=(
                'Run script in batch mode, i.e. run every pep file in '
                'directory.'
                ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Run
for instance in Data:
# }}}

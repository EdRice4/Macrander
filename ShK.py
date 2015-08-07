# {{{ Header
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code written by: Edwin Rice
# email: edwinricethe4th@gmail.com
# phone: +1 (513) 426-4187
# github: https://github.com/EdRice4
#
# Initiation date|time: 08/05/2015|18:16:26
# }}}

# {{{ Imports
from os import getcwd, listdir
import re
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

    # {{{ write_pep_dict_to_file
    def write_pep_dict_to_file(self):

        """ {{{ Docstrings

        Writes the filtered pep dict into a pep file with the format of:

            '>Pep_ID\n'
            'Pep_seq\n'

        }}} """

        # Open filtered pep file in write mode
        with open(self._filtered_pep_file, 'w') as pep_file:
            for k, v in self._filtered_pep_dict:
                # Format ID and sequence to print correclty
                Pep_ID = '>{0}\n'.format(k)
                Pep_seq = '{0}\n'.format(v)
                # Write to file
                pep_file.write(Pep_ID)
                pep_file.write(Pep_seq)
    # }}}
# }}}


# {{{ SearchParse
class SearchParse(PepFile):

    """ {{{ Docstrings

    A class in which all regular expression searching and subsequent data
    parsing is stored.

    }}} """

    # {{{ SearchPattern
    # First, must define serach pattern
    amino_acids = [
            'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
            'F', 'P', 'S', 'T', 'W', 'Y', 'V'
            ]
    # "|" in regular expression speak means "or"
    amino_acids = '|'.join(amino_acids)
    # "{m}" in regular expression speak means "match exactly this number of
    # preceding expression"
    C3C2C = (
            'C[' + amino_acids + ']{3}C[' + amino_acids + ']{2}C'
            )
    # "+" in regular expression speak means "match 1 or more repetitions of
    # preceding expression"
    CNCNC = (
            'C[' + amino_acids + ']+C[' + amino_acids + ']+C'
            )
    # Compile into pattern re.search can utilize
    C3C2C = re.compile(C3C2C)
    CNCNC = re.compile(CNCNC)
    # }}}

    # {{{ search_the_6_Cs (Arr, matey)
    def search_the_6_Cs(self, pep_seq):

        """ {{{ Docstrings

        Utilizing regular expressions (as provided by the "re" python module),
        searches a peptide sequence for the pattern: C???C??C where ? = any
        amino acid.

        }}} """

        # Does peptide sequence contain C3C2C pattern?
        contains_C3C2C = re.search(C3C2C, pep_seq)
        # If not, return FALSE
        if not contains_C3C2C:
            return FALSE
        # If so, continue
        else:
            # Get index of pattern
            C3C2C_position = pep.index(
                    contains_shk.group()
                    )
            # Truncate peptide sequence so only includes 50 amino acids that
            # precede occurence of pattern
            truncated_pep_seq = pep_seq[
                    C3C2C_position - 50:C3C2C_position
                    ]
            # Does truncated peptide sequence contain 3 preceding Cs?
            contains_CNCNC = re.seach(CNCNC, truncated_pep_seq)
            # If so, return whole peptide sequence
            if contains_CNCNC:
                return pep_seq
            # If not, return FALSE
            else:
                return FALSE
    # }}}

    # {{{ filter_pep_dict
    def filter_pep_dict(self):

        """ {{{ Docstrings

        Utilizing previous function, search_the_6_Cs, filters peptide
        sequences, only keeping those which return a TRUE value (i.e. contain
        the aforementioned pattern).

        }}} """

        # Initialze empty dictionary where peptide sequences containing the
        # aforementioned pattern will be stored
        filtered_pep_dict = {}
        # Itereate over keys, values in pep_dict
        for k, v in self._pep_dict.iteritems():
            if search_the_6_Cs(v):
                filtered_pep_dict[k] = v
        # Return filtered_pep_dict
        return filtered_pep_dict
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
        self._pep_dict = self.read_pep_into_dict(pep_file)
        self._filtered_pep_dict = self.filter_pep_dict()
        self._filtered_pep_file = pep_file.replace('.pep', '_filtered.pep')
        self.write_pep_dict_to_file()
        self.registry.append(self)
    # }}}
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


# {{{ Run
for instance in Data:
# }}}

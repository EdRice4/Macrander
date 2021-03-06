#!/usr/bin/env python
# -*- coding: utf-8 -*-


# {{{ Header
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
from sys import exit
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
        # Filter out IDs and strip each entry of leading/trailing whitespace
        Pep_IDs = [line.strip() for line in pep if line[0] == '>']
        # Filter out peptide sequences and strip each entry of
        # leading/trailing whitespace
        Pep_seqs = [line.strip() for line in pep if line[0] != '>']
        # Concatenate into dictionary
        pep_dict = dict(zip(Pep_IDs, Pep_seqs))
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
            for k, v in self._filtered_pep_dict.iteritems():
                # Set values for readability
                Pep_ID = k
                # Handle addition of newline character in manner contingent
                # upon length; improve readability when writing
                # If only one match returned, append newline "\n" character
                if len(v) == 1:
                    Pep_seq = v[0] + '\n'
                # Else, join each item in list with newline character
                else:
                    Pep_seq = '\n'.join(v)
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

    # {{{ search_the_6_Cs (Arr, matey)
    def search_the_6_Cs(self, pep_seq):

        """ {{{ Docstrings

        Utilizing regular expressions (as provided by the "re" python module),
        searches a peptide sequence for occurrences of the ShK domain as
        defined below.

        }}} """

        # {{{ SearchPattern
        # First, must define serach pattern
        # All possible values in peptide sequence; a little bit more specific
        # than "\w"
        amino_acids = [
                'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K',
                'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
                ]
        # "|" in regular expression speak means "or"
        amino_acids = '|'.join(amino_acids)
        # In regular expression speak, "{m}" means "match exactly this number
        # of the preceding expression", "+" means "match 1 or more repetitions
        # of preceding expression", and {m,n} means "match m through n
        # repetitions of preceding expression, attempting to match as many as
        # possible"
        # ::MODIFIABLE::
        # NOTE: If you would like to search further upstream for preceding
        # cysteines, simply modify second integer in {1,50}
        ShK_domain = (
                'C[' + amino_acids + ']+C[' + amino_acids + ']+C[' +
                amino_acids + ']{1,50}C[' + amino_acids + ']{3}C[' +
                amino_acids + ']{2}C'
                )
        # Compile into pattern re.findall can utilize
        ShK_domain = re.compile(ShK_domain)
        # }}}

        # Does peptide sequence contain ShK domain?
        putative_ShK_domains = re.findall(ShK_domain, pep_seq)
        # If matches found, return list of all matches (as strings)
        if putative_ShK_domains:
            return putative_ShK_domains
        # Else, return false
        else:
            return False
    # }}}

    # {{{ filter_pep_dict
    def filter_pep_dict(self):

        """ {{{ Docstrings

        Utilizing previous function, search_the_6_Cs, filters peptide
        sequences, only keeping those which return a TRUE value (i.e. contain
        at least one occurrence the ShK domain).

        }}} """

        # Initialze empty dictionary where peptide sequences containing at
        # at least one occurrence of the ShK domain will be stored
        filtered_pep_dict = {}
        # Itereate over keys, values in pep_dict
        for k, v in self._pep_dict.iteritems():
            # Set values for readability
            Pep_ID = k
            Pep_seq = self.search_the_6_Cs(v)
            # If it contains the ShK domain, include it in filtered pep dict,
            # else, discard
            if Pep_seq:
                filtered_pep_dict[Pep_ID] = Pep_seq
        # Return filtered_pep_dict
        return filtered_pep_dict
    # }}}
# }}}


# {{{ Data
class Data(SearchParse):

    """ {{{ Docstrings

    A class in which all pertinent data for each repsective peptide file
    is stored.

    }}} """

    # {{{ __init__
    def __init__(self, pep_file):
        self._pep_dict = self.read_pep_into_dict(pep_file)
        self._filtered_pep_dict = self.filter_pep_dict()
        self._filtered_pep_file = pep_file.replace('.pep', '_filtered.pep')
        self.write_pep_dict_to_file()
    # }}}
# }}}


# {{{ Main argument parser
arg_parser = argparse.ArgumentParser(
        prog='ShK.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
                'Given a list of peptide sequences, putatively containing '
                'the ShK domain (as described in Rangaraju et al., "Potassium '
                'Channel Modulation by a Toxin Domain in Matrix '
                'Metalloprotease 23." DOI: 10.1074/jbc.M109.071266 and in the '
                'Simple Modular Architecture Research Tool database: '
                'http://smart.embl.de/smart/do_annotation.pl?DOMAIN=SM00254), '
                'will robustly filter out those peptide sequences and '
                'characterize the nature of each sequence.'
                ),
        epilog=(
                'Note: When running in batch mode, all files you wish to run '
                'should contain the string ".pep" in the file name and all '
                'files containing this string in the name will be run. For '
                'explanations on potentially useful, user configurable '
                'variables, search for ::MODIFIABLE:: in file.'
                )
        )
arg_parser.add_argument(
        '-pep_file', type=str, help=(
                'Name of pep file you wish to run.'
                ),
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
# Instantiate instance of Data class with every pep file in directory if user
# specified batch flag
if args.batch:
    # Get current working directory
    cwd = getcwd()
    # Get all files in cwd
    files = listdir(cwd)
    # Filter out pep files
    # NOTE: All files you wish to run should contain the string '.pep' in the
    # name and all files containing this string will be run if you specify
    # the batch flag
    pep_files = [x for x in files if '.pep' in x]
    # Instantiate intances of Data class for all pep files found
    for i in pep_files:
        Data(i)
# Else, utilize user-specified string to instantiate single instance of Data
# class
else:
    # Exit if no proper file specified
    if not args.pep_file:
            exit(
                'You did not specify a pep file to run nor did you run the '
                'script in batch mode. Try again.'
                )
    # Else, continue
    Data(args.pep_file)
# }}}

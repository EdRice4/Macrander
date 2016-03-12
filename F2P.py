#!/usr/bin/env python
# -*- coding: utf-8 -*-


# {{{ Header
# Code written by: Edwin Rice
# email: edwinricethe4th@gmail.com
# phone: +1 (513) 426-4187
# github: https://github.com/EdRice4
#
# Initiation date|time: 08/10/2015|14:55:58
# }}}


# {{{ Imports
from os import getcwd, listdir
from random import randrange
import argparse
# }}}


# {{{ FastaFile
class FastaFile(object):

    """ {{{ Docstrings

    A class in which the fasta file is read.

    }}} """

    # {{{ get_original_data
    def get_original_data(self, fasta_file):

        """ {{{ Docstrings

        Reads fasta file into dictionary, whose keys are the original sequence
        identifiers and values are the corresponding sequences, given the
        name of the fasta file as a string.

        }}} """

        # Open fasta file in read mode
        with open(fasta_file, 'r') as fasta:
            # Read into list
            fasta = fasta.readlines()
        # Filter out sequence IDs and strip each entry of leading/trailing
        # whitespace
        # TODO: Remove greater than (">") symbol?
        seq_IDs = [line.strip() for line in fasta if line[0] == '>']
        # Filter out nucleotide sequences and strip each entry of
        # leading/trailing whitespace
        seqs = [line.strip() for line in fasta if line[0] != '>']
        # Concatenate into dictionary
        seq_dict = dict(zip(seq_IDs, seqs))
        return seq_dict
    # }}}

    # {{{ generate_unique_ids
    def generate_unique_ids(self, original_seq_dict):

        """ {{{ Docstrings

        Randomly generate unique identifier (integer from 0 through 9999999999)
        and use as key in two dictionaries. The first contains each unique ID
        as keys and corresponding sequences as values. The second contains each
        unique ID as keys and each corresponding original ID as values.

        }}} """

        # Initiate two empy dictionaries in which pertinent variables are
        # stored
        unique_seq_dict = {}
        unique_ID_dict = {}
        # Iterate over key, value pairs
        for ID, seq in original_seq_dict.iteritems():
            # Define variables for readability
            original_ID = ID
            # Generate unique ID
            unique_ID = randrange(0, 9999999999)
            # Define key, value paris in novel dictionaries
            unique_seq_dict[unique_ID] = seq
            unique_ID_dict[unique_ID] = original_ID
        return(unique_seq_dict, unique_ID_dict)
    # }}}
# }}}


# {{{ PhylipFile
class PhylipFile(FastaFile):

    """ {{{ Docstrings

    A class in which phylip file output functionality is stored.

    }}} """

    # {{{ add_args
    @staticmethod
    def add_args():

        """ {{{ Docstrings

        Add argument group "phylip" to namespace and subsequent pertinent
        arguments to aforementioned group.

        }}} """

        args_phylip = arg_parser.add_argument_group(
                'phylip', 'Arguments pertaining to phylip file output.'
                )
        args_phylip.add_argument(
                '-s', '--sequential', help=(
                        'Write phylip file in sequential format.'
                        ),
                action='store_true'
                )
    # }}}

    # {{{ __init__
    def __init__(self, fasta_file):

        """ {{{ Docstrings

        Upon instantiation of PhylipFile instance, run all pertinent
        functions.

        }}} """

        original_phylip_name = self._path.replace(
                '.fasta', '_original.phylip'
                )
        unique_phylip_name = self._path.replace(
                '.fasta', '_unique.phylip'
                )
        dictionary_name = self._path.replace('.fasta', '.txt')
        original_seq_dict = self.get_original_data(fasta_file)
        unique_seq_dict, unique_ID_dict = self.generate_unique_ids(
                original_seq_dict
                )
        if args.sequential:
            self.write_phylip_file_sequential(
                    original_seq_dict, original_phylip_name
                    )
            self.write_phylip_file_sequential(
                    unique_seq_dict, unique_phylip_name
                    )
        else:
            self.write_phylip_file_interleaved(
                    original_seq_dict, original_phylip_name
                    )
            self.write_phylip_file_interleaved(
                    unique_seq_dict, unique_phylip_name
                    )
        self.write_dictionary_file(unique_ID_dict, dictionary_name)
    # }}}

    # {{{ write_phylip_file_sequential
    def write_phylip_file_sequential(self, seq_dictionary, phylip_file):

        """ {{{ Docstrings

        A phylip file is written in sequential format, given a sequence
        dictionary and name of novel phylip file to be written to.

        }}} """

        # Open phylip file in write mode
        with open(phylip_file, 'w') as phy:
            # Iterate over key, value pairs
            for ID, seq in seq_dictionary.iteritems():
                # Define line
                line = '{0}\t{1}\t'.format(ID, seq)
                # Write to file
                phy.write(line)
    # }}}

    # {{{ write_phylip_file_interleaved
    def write_phylip_file_interleaved(self, seq_dictionary, phylip_file):

        """ {{{ Docstrings

        A phylip file is written in interleaved format, given a sequence
        dictionary and name of novel phylip file to be written to.

        }}} """

        # Open phylip file in write mode
        with open(phylip_file, 'w') as phy:
            # Iterate over key, value pairs
            for ID, seq in seq_dictionary.iteritems():
                # Define line
                line = '{0}\t{1}\n'.format(ID, seq[0:51])
                # Write to file, only including first 50 characters of
                # "seq" value
                phy.write(line)
            # Write empty line to file
            phy.write('\n')
            # Iterate over values
            for seq in seq_dictionary.itervalues():
                # Iterate over sequence, starting from nucleotide at position
                # 51, in increments of 50 nucleotides
                for position in range(51, len(seq), 50):
                    # Define line
                    line = '{0}\n'.format(seq[position:position + 50])
                    # Write nucleotides to file
                    phy.write(line)
                # Write empty line to file
                phy.write('\n')
        # }}}

    # {{{ write_dictionary
    def write_dictionary_file(self, unique_ID_dict, dictionary_name):

        """ {{{ Docstrings

        Writes unique and original IDs to tab-delimited file.

        }}} """

        # Open dictionary in write mode
        with open(dictionary_name, 'w') as dfile:
            # Iterate over key, value pairs
            for unique, orig in unique_ID_dict.iteritems():
                # Define line
                line = '{0}\t{1}\n'.format(unique, orig)
                # Write to file
                dfile.write(line)
    # }}}
# }}}


# {{{ FastaFile
class FastaFile(PhylipFile):

    """ {{{ Docstrings

    A class in which all of the necessary parameters corresponding to each
    respective fasta file are stored.

    }}} """

    # {{{ add_args
    @staticmethod
    def add_args():

        """ {{{ Docstrings

        Add mutually exclusive argument group "pep_batch" to namespace and
        subsequent pertinent arguments to aforementioned group.

        }}} """

        pep_batch = arg_parser.add_mutually_exclusive_group()
        pep_batch.add_argument(
                '-f', '--fasta_file', type=str, help=(
                        'Name of fasta sequence to be converted.'
                        ),
                default=None
                )
        pep_batch.add_argument(
                '-b', '--batch', help=(
                        'Run script in batch mode. i.e. convert all fasta'
                        'files in directory to phylip.'
                        ),
                action='store_true'
                )
    # }}}

    # {{{ __init__
    def __init__(self, fasta_file):

        """ {{{ Docstrings

        Upon instantiation of FastaFile instance, run all pertinent
        functions.

        }}} """

        self._path = fasta_file
        PhylipFile.__init__(self, fasta_file)
    # }}}

# }}}


# {{{ Main argument parser
arg_parser = argparse.ArgumentParser(
        prog='Fasta2Phylip.py',
        description=(
                'Converts a fasta file to phylip format, either writing it in '
                'sequential or interleaved format, depending on user '
                'specification. In total, 3 files will be produced: 1.) A '
                'phylip file with original IDs. 2.) A phylip file with unique '
                'IDs. 3.) A dictionary file with unique and original IDs.'
                'The intended workflow is analgous to the following: '
                'fasta_file >F2P.py> phylip_file_unique >Analysis> '
                'output_unique >Rm_Cont.py> output_original. Designed to '
                'circumvent the fact that many analyses requiring phylip '
                'files as input truncate IDs at 10 chracters.'
                ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
# Run add_args for each class when passing "-h" flag and prior to instantiating
# instances of any class
if __name__ == '__main__':
    PhylipFile.add_args()
    FastaFile.add_args()
# Parse args into namespace objects
args = arg_parser.parse_args()
# Print informative output
if args.sequential:
    print(
            'You have chosen to write the phylip file in sequential '
            'format.'
            )
else:
    print(
            'You have chosen to write the phylip file in interleaved '
            'format.'
            )
# }}}


# {{{ Batch
# Instantiate instance of FastaFile class with every fasta file in directory if
# user specified batch flag
if args.batch:
    # Get current working directory
    cwd = getcwd()
    # Get all files in cwd
    files = listdir(cwd)
    # Filter out fasta files
    # NOTE: All files you wish to run should contain the string '.fasta' in the
    # name and all files containing this string will be run if you specify
    # the batch flag
    fasta_files = [x for x in files if '.fasta' in x]
    # Instantiate intances of FastaFile class for all fasta files found
    for fasta_file in fasta_files:
        FastaFile(fasta_file)
# Else, utilize user-specified string to instantiate single instance of Data
# class
else:
    FastaFile(args.fasta_file)
# }}}

# {{{ Imports
from os import getcwd, listdir
from random import randrange
import argparse
# }}}


# {{{ IDSeq
class IDSeq(object):

    """ {{{ Docstrings

    A class in which the fasta file is read.

    Namely:

            1.) Two lists are generated; one containing each identifier and
                the other containing each corresponding sequence.

    }}} """

    def generate_original_ids_and_seqs(self):

        """ {{{ Docstrings
        Reads fasta file into list.
        }}} """

        for i in self.fasta_file:
            if i[0] == '>':
                self.identities.append(i.strip('>\r\n'))
            else:
                self.sequences.append(i.strip())
# }}}


# {{{ PhylipFile
class PhylipFile(IDSeq):

    """ {{{ Docstrings

    A class in which phylip file output is stored.

    Namely:

        1.) The original (i.e. containing the sequence IDs from the fasta file)
            phylip file is written in sequential format.

        2.) The original phylip file is written in interleaved format.

        3.) The unique (i.e. containing the randomly generated [see
            generate_unique_ids]) phylip file is written in sequential format.

        4.) The unique phylip file is written in interleaved format.

    }}} """

    # {{{ write_original_phylip_file_sequential
    def write_original_phylip_file_sequential(self):

        """ {{{ Docstrings
        The original phylip file is written in sequential format.
        }}} """

        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j + '\n')
    # }}}

    # {{{ write_original_phylip_file_interleaved
    def write_original_phylip_file_interleaved(self):

        """ {{{ Docstrings
        The original phylip file is writtin in interleaved format.
        }}} """

        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j[0:50] + '\n')
            phy.write('\n')
            for i in range(50, len(self.sequences[0]), 50):
                for j in self.sequences:
                    phy.write(j[i:i + 50] + '\n')
                phy.write('\n')
    # }}}

    # {{{ write_unique_phylip_file_sequential
    def write_unique_phylip_file_sequential(self):

        """ {{{ Docstrings
        The unique phylip file is written in sequential format.
        }}} """

        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j + '\n')
    # }}}

    # {{{ write_unique_phylip_file_interleaved
    def write_unique_phylip_file_interleaved(self):

        """ {{{ Docstrings
        The unique phylip file is written in interleaved format.
        }}} """

        with open(self.unique_phylip_name, 'w') as phy:
            for i, j in zip(self.unique_identities.itervalues(),
                            self.sequences):
                phy.write(str(i) + '\t' + j[0:50] + '\n')
            phy.write('\n')
            for i in range(50, len(self.sequences[0]), 50):
                for j in self.sequences:
                    phy.write(j[i:i + 50] + '\n')
                phy.write('\n')
    # }}}
# }}}


# {{{ DictionaryFile
class DictionaryFile(PhylipFile):

    """ {{{ Docstrings

    A class in which dictionary file ouput is stored.

    Namely:

        1.) The unique IDs are generated and stored in a dictionary in the
            format of:

                dictionary = {
                        'XXXXXXXXX' : 'Original_ID',
                        'XXXXXXXXX' : 'Original_ID'
                        }

        2.) The dictionary file is written in the format of:

                Unique_ID\tOriginal_ID
                'XXXXXXXXX'\t'Original_ID'

            Where XXXXXXXXX is a randomly generated integer between the
            integers 0 and 999999999.

    }}} """

    # {{{ generate_unique_ids
    def generate_unique_ids(self):

        """ {{{ Docstrings
        The unique IDs are generated and stored in a dictionary.
        }}} """

        for i in self.identities:
            self.unique_identities[randrange(0, 9999999)] = i
    # }}}

    # {{{ write_dictionary
    def write_dictionary_file(self):

        """ {{{ Docstrings
        The unique and original IDs are written to a file.
        }}} """

        with open(self.dictionary_name, 'w') as dfile:
            for key, value in self.unique_identities.iteritems():
                dfile.write(str(key) + '\t' + value + '\n')
    # }}}
# }}}


# {{{ IterRegistry
class IterRegistry(type):

    """ {{{ Docstrings
    A metaclass to allow for iteration of instances of the FastaFile class.
    }}} """

    # {{{ __iter__
    def __iter__(cls):
        return iter(cls.registry)
    # }}}
# }}}


# {{{ FastaFile
class FastaFile(DictionaryFile):

    """ {{{ Docstrings
    A class in which all of the necessary parameters corresponding to each
    respective fasta file are stored.


    Namely:

        1.) The sequences IDs are stored in a list, self.identities.

        2.) The sequences themselves are stored in a list, self.sequences.

        3.) Randomly generated unique identities (see generate_unique_ids) are
            stored in a dictionary, self.unique_identities.

    }}} """

    # {{{ __metaclass__
    __metaclass__ = IterRegistry
    registry = []
    # }}}

    # {{{ __init__
    def __init__(self, fasta_file):
        self.path = fasta_file
        with open(self.path, 'r') as phy:
            self.fasta_file = phy.readlines()
        self.original_phylip_name = self.path.replace('.fasta',
                                                      '_original.phylip')
        self.unique_phylip_name = self.path.replace('.fasta',
                                                    '_unique.phylip')
        self.dictionary_name = self.path.replace('.fasta', '.txt')
        self.identities = []
        self.sequences = []
        self.unique_identities = {}
        self.registry.append(self)
    # }}}
# }}}


# {{{ ArgParser
arg_parser = argparse.ArgumentParser(
        prog='Fasta2Phylip.py',
        description=(
                'Converts a fasta file to phylip format, either writing it in'
                'sequential or interleaved format, depending on user'
                'specification. In total, 3 files will be produced: 1.) A '
                'phylip file with original IDs. 2.) A phylip file with unique '
                'IDs. 3.) A dictionary file with unique and original IDs.'
                'The intended workflow is analgous to the following: '
                'fasta_file >F2P.py> phylip_file_unique >Analysis> '
                'output_unique >Rm_Cont.py> output_original. Designed to '
                'circumvent the fact that many analysis requiring phylip '
                'files as input truncate IDs at 10 chracters.'
                ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
arg_parser.add_argument(
        'fasta_file', type=str, help='Name of fasta sequence to be converted.',
        default=None
        )
arg_parser.add_argument(
        '-s', '--sequential', help='Write phylip file in sequential format.',
        action='store_true'
        )
arg_parser.add_argument(
        '-b', '--batch', help=(
                'Run script in batch mode; i.e. convert all fasta files in '
                'directory to phylip.'
                ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Batch
if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    fasta_files = [x for x in fid if '.fasta' in x]
    for i in fasta_files:
        FastaFile(i)
else:
    FastaFile(args.fasta_file)
# }}}


# {{{ Run
for fasta in FastaFile:
    fasta.generate_original_ids_and_seqs()
    fasta.generate_unique_ids()
    if args.sequential:
        print(
                'You have chosen to write the phylip file in sequential '
                'format.'
                )
        fasta.write_original_phylip_file_sequential()
        fasta.write_unique_phylip_file_sequential()
    else:
        print(
                'You have chosen to write the phylip file in interleaved '
                'format.'
                )
        fasta.write_original_phylip_file_interleaved()
        fasta.write_unique_phylip_file_interleaved()
    fasta.write_dictionary_file()
# }}}

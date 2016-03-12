# {{{ Imports
from os import getcwd, listdir
from string import maketrans
from re import search, sub
import argparse
# }}}

# {{{ ExtractData
class ExtractData(object):

    """ {{{ Docstrings
    A class in which all data extraction functionality, i.e. data is simply
    read, functionality is stored.

    Namely:

        1.) All pertinent parameters, including the sequence ID, nucleotides,
            and strand, are extracted from the TransDecoder output and stored
            in a dictionary with the form of:

                dictionary = {
                        'Seq_ID' : [
                                (start, end, 'strand'),
                                (start, end, 'strand')
                                ]
                        }

            Note that the dictionary supports multiple hits for the same
            sequence ID and that the parameters are stored as a list of
            tuples.

        2.) All pertinent parameters, including sequence ID and sequence, are
            extracted from the fasta file and stored in a dictionary with the
            form of:

                dictionary = {
                        'Seq_ID' : 'sequence'
                        }

            Also of note, all new line characters (\n) are removed to prevent
            them from interfering with counting the nucleotides.

    }}} """

    # {{{ extract_params_from_pep
    def extract_params_from_pep(self):

        """ {{{ Docstrings
        Reads TransDecoder output into dictionary.
        }}} """

        seq_of_int = {}
        # Peptide sequence information not pertinent
        pep = filter(lambda x: x[0] == '>', self.pep)
        for i in pep:
            # Sequence ID, nucleotides, and strand, respectively
            data = search('(\w+):(\d+)-(\d+)\((\+|-)\)', i)
            diff = int(data.group(2)) - int(data.group(3))
            # If forward strand (+)
            if diff < 0:
                # Python begins indexing at 0 while TD begins at 1
                start = int(data.group(2)) - 1
                # TODO: WHY AM I DOING THIS?
                end = start - (diff - 1)
            # If reverse strand (-)
            else:
                start = int(data.group(3)) - 1
                # TODO: WHY AM I DOING THIS?
                end = start + (diff + 1)
            # If already present in dictionary
            if seq_of_int.get(data.group(1)):
                seq_of_int[data.group(1)].append((start, end, data.group(4)))
            # If not already present in dictionary
            else:
                seq_of_int[data.group(1)] = [(start, end, data.group(4))]
        return seq_of_int
    # }}}

    # {{{ extract_params_from_fas
    def extract_params_from_fas(self):

        """ {{{ Docstrings
        Reads fasta file into dictionary. Note: Please ensure fasta file
        contains no empty lines. Use the command:

            sed '/^$/d/ [file] > [output]

        }}} """

        fasta_dict = {}
        # Python runs the split function first which generates an empty line
        # at the beginning
        fasta = self.fas[1:].split('>')
        # Use for loop as multiple manipulations of 'i'
        for i in fasta:
            # {{{ Modifiable
            # ::Modifiable::
            # The first argument in the split function is the character that
            # separates each sequence ID with its respective sequence. This
            # character should not be present in the sequence ID itself.
            # The script assumes that each ID is on a separate line so that
            # the fasta file looks like:
            #
            #       >Seq_ID\n
            #       sequence\n
            #       sequence\n ...
            #
            # }}}
            i = i.split('\n', 1)
            seq_id = i[0].strip()
            seq = i[1].replace('\n', '')
            fasta_dict[seq_id] = seq
        return fasta_dict
    # }}}
# }}}


# {{{ DataParse
class DataParse(ExtractData):

    """ {{{ Docstrings
    A class in which all data parsing, i.e. the boundary between disparate file
    types, functionality is stored.

    Namely:

        1.) If applicable (if reverse), the reverse complement of the strand
            is generated.

        2.) Data from TransDecoder is used to extract the pertinent nucleotides
            from the fasta file, and if applicable, the reverse complement of
            the sequence is generated (utilizing the previous function), which
            is/are stored in a filtered fasta dictionary in the form of:

            dictionary = {
                    'Seq_ID' : [
                            'sequence',
                            'sequence'
                            ]
                    }

    }}} """

    # {{{ rev_compl
    def rev_compl(params, seq):

        """ {{{ Docstrings
        Returns the reverse complement of a given sequence, if applicable.
        }}} """

        # Define variables
        start, end, strand = params
        tbl = maketrans('ATCG', 'TAGC')
        # If forward strand, do nothing
        if strand == '+':
            return seq[start:end]
        # If reverse, return reverse complement
        else:
            seq = seq[start:end]
            # Reverse
            rev_seq = seq[::-1]
            # Reverse complement
            rev_compl_seq = rev_seq.translate(tbl)
            return rev_compl_seq
    # }}}

    # {{{ extract_pertinent_seq
    def extract_pertinent_seq(self, fd, soi):

        """ {{{ Docstrings
        Returns pertinent sequences, and if applicable, reverse complement of
        such, and stores in dictionary.
        }}} """

        ffd = {}
        for i in soi.iteritems():
            # Get sequence ID
            seq_id = i[0]
            # Get sequence from fasta dictionary
            seq = fd[i[0]]
            # Build filtered fasta dictionary
            ffd[seq_id] = map(
                    lambda x: DataParse.rev_compl(x, seq), i[1]
                    )
        return ffd
    # }}}
# }}}


# {{{ FileIO
class FileIO(ExtractData):

    """ {{{ Docstrings
    A class in which all file input/output functionality is stored.

    Namely:

        1.) The resulting sequence IDs and their respective sequences are
            manipulated in order to improve readability when written. That is
            to say that a '\n' (newline) character is inserted at every 50
            nucleotides.

        2.) The 'pretty' dictionary generated by the previous function is
            written to a fasta file.

    }}} """

    # {{{ make_pretty
    def make_pretty(self, rcfd):

        """ {{{ Docstrings
        The filtered fasta dictionary is manipulated to improve readability
        when written.
        }}} """

        pretty_dict = {}
        repl_funct = lambda m: m.group(0) + '\n'
        for i in rcfd.iteritems():
            pretty_dict[i[0]] = map(lambda x: sub('[A-Z]{50}', repl_funct,
                                    x), i[1])
        return pretty_dict
    # }}}

    # {{{ write_dict
    def write_dict(self, pd):

        """ {{{ Docstrings
        The pretty dictionary is written to a fasta file.
        }}} """

        with open(self.fas_new, 'w') as fas:
            for i in pd.iteritems():
                fas.write('>%s\n' % i[0])
                fas.write(('\n>%s\n' % i[0]).join(i[1]))
                fas.write('\n')
    # }}}
# }}}


# {{{ IterRegistry
class IterRegistry(type):

    """ {{{
    A metaclass allowing for iterations over the PepFile and FastaFile
    class.
    }}} """

    # {{{ __iter__
    def __iter__(cls):
        return iter(cls.registry)
    # }}}
# }}}


# {{{ PepFastaFile
class PepFastaFile(FileIO):

    """ {{{ Docstrings
    A class in which all the necessary parameters corresponding to each
    respective peptide file are stored.
    }}} """

    # {{{ __metaclass__
    __metaclass__ = IterRegistry
    registry = []
    # }}}

    # {{{ __init__
    def __init__(self, pep_file, fasta_file):
        with open(pep_file, 'r') as pep:
            self.pep = pep.readlines()
        with open(fasta_file, 'r') as fas:
            self.fas = fas.read()
        self.fas_new = fasta_file.replace('.fasta', '_new.fasta')
        self.registry.append(self)
    # }}}
# }}}


# {{{ ArgParser
arg_parser = argparse.ArgumentParser(
        prog='Peptide2Nucleotide.py',
        description=(
                'Converts peptide sequence output from TransDecoder into '
                'pertinent sequences, and if applicable, generates reverse '
                'complement of aforementioned sequences.'
                ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
arg_parser.add_argument(
        'pep_file', type=str,
        help=(
                'Name of TransDecoder output file containing peptide '
                'sequences of interest.'
                ),
        default=None
        )
arg_parser.add_argument(
        'fasta_file', type=str,
        help=(
                'Name of fasta file containing nucleotide sequences of '
                'interest. See docstrings for extract_params_from_fas and '
                '\'::Modifiable::\' for notes on the necessary format for '
                'fasta files.'
        ),
        default=None
        )
arg_parser.add_argument(
        '-b', '--batch',
        help=(
                'Run script in batch mode. i.e. perform sequence selection '
                'with all peptide sequences in directory and their respective '
                'fasta files.'
                ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Batch
if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    pep_files = sorted(filter(lambda x: '.pep.' in x, fid))
    fasta_files = sorted(filter(lambda x: '.fasta' in x, fid))
    for i, j in zip(pep_files, fasta_files):
        PepFastaFile(i, j)
else:
    PepFastaFile(args.pep_file, args.fasta_file)
# }}}


# {{{ Run
for f in PepFastaFile:
    seqs_of_int = f.extract_params_from_pep()
    fasta_dict = f.extract_params_from_fas()
    filt_fasta_dict = f.extract_pertinent_seq(fasta_dict, seqs_of_int)
    pretty_dict = f.make_pretty(filt_fasta_dict)
    f.write_dict(pretty_dict)
# }}}

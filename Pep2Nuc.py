# {{{ Imports
from os import getcwd, listdir
from string import maketrans
from re import search, sub
import argparse
# }}}

# {{{ DataParse
class DataParse(object):

    """ {{{ Docstrings
    A class in which all data parsing functionality is stored.

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


class ExtractData(DataParse):

    """A class in which all data extraction, i.e. the boundary between
       disparate file types, functionality is stored."""

    def extract_pertinent_seq(self, fd, soi):
        ffd = {}
        for i in fd.iteritems():
            if soi.get(i[0], False):
                ffd[i[0]] = [i[1][soi[i[0]][j]:soi[i[0]][j + 1]]
                             for j in range(0, len(soi[i[0]]), 3)]
        return ffd

    def rev_compl(self, ffd, soi):
        tbl = maketrans('ATCG', 'TAGC')
        rev_compl_fasta_dict = {}
        for i in ffd.iteritems():
            if soi[i[0]][2] == '-':
                rev_compl_fasta_dict[i[0]] = map(lambda x:
                                                 x[::-1].translate(tbl), i[1])
                #rev_compl_fasta_dict[i[0]] = map(lambda x: x.translate(tbl),
                #                                 i[1])
            else:
                rev_compl_fasta_dict[i[0]] = ffd[i[0]]
        return rev_compl_fasta_dict


class FileIO(ExtractData):

    """A class in which all file input/output functionality is stored."""

    def make_pretty(self, rcfd):
        pretty_dict = {}
        repl_funct = lambda m: m.group(0) + '\n'
        for i in rcfd.iteritems():
            pretty_dict[i[0]] = map(lambda x: sub('[A-Z]{50}', repl_funct,
                                    x), i[1])
        return pretty_dict

    def write_dict(self, pd):
        with open(self.fas_new, 'w') as fas:
            for i in pd.iteritems():
                fas.write('>%s\n' % i[0])
                fas.write(('\n>%s\n' % i[0]).join(i[1]))
                fas.write('\n')


class IterRegistry(type):

    """A metaclass allowing for iterations over the PepFile and FastaFile
       class."""

    def __iter__(cls):
        return iter(cls.registry)


class PepFastaFile(FileIO):

    """A class in which all the necessary parameters corresponding to each
       respective peptide file are stored."""

    __metaclass__ = IterRegistry
    registry = []

    def __init__(self, pep_file, fasta_file):
        with open(pep_file, 'r') as pep:
            self.pep = pep.readlines()
        with open(fasta_file, 'r') as fas:
            self.fas = fas.read()
        self.fas_new = fasta_file.replace('.fasta', '_new.fasta')
        self.registry.append(self)


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('pep_file', type=str, help=('Name of TransDecoder '
                                                    'output file containing '
                                                    'peptide sequences of '
                                                    'interest. Specify None '
                                                    'if running script in '
                                                    'batch mode.')
                        )
arg_parser.add_argument('fasta_file', type=str, help=('Name of fasta file '
                                                      'containing nucleotide '
                                                      'sequences of '
                                                      'interest. Specify '
                                                      'None if running '
                                                      'script in batch mode.')
                        )
arg_parser.add_argument('-b', '--batch', help=('Run script in batch mode. '
                                               'i.e. perform sequence '
                                               'selection with all peptide '
                                               'sequences in directory and '
                                               'their respective fasta '
                                               'files.'),
                        action='store_true'
                        )
args = arg_parser.parse_args()

if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    pep_files = sorted(filter(lambda x: '.pep.' in x, fid))
    fasta_files = sorted(filter(lambda x: '.fasta' in x, fid))
    for i, j in zip(pep_files, fasta_files):
        PepFastaFile(i, j)
else:
    PepFastaFile(args.pep_file, args.fasta_file)

for f in PepFastaFile:
    seqs_of_int = f.extract_params_from_pep()
    fasta_dict = f.extract_params_from_fas()
    filt_fasta_dict = f.extract_pertinent_seq(fasta_dict, seqs_of_int)
    rev_compl_fasta_dict = f.rev_compl(filt_fasta_dict, seqs_of_int)
    pretty_dict = f.make_pretty(rev_compl_fasta_dict)
    f.write_dict(pretty_dict)

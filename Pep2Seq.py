from os import getcwd, listdir
import argparse


class IterRegistry(type):

    """A metaclass allowing for iterations over the PepFile and FastaFile
       class."""

    def __iter__(cls):
        return iter(cls.registry)


class PepFile(object):

    """A class in which all the necessary parameters corresponding to each
       respective peptide file are stored."""

    __metaclass__ = IterRegistry
    registry = []

    def __init__(self, pep_file):
        self.title = pep_file
        self.registry.append(self)


class FastaFile(object):

    """A class in which all the necessary parameters corresponding to each
       respective fasta file are stored."""

    __metaclass__ = IterRegistry
    registry = []

    def __init__(self, fasta_file):
        self.title = fasta_file
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
    pep_files = sorted(map(lambda x: '.pep.' in x, fid))
    fasta_files = sorted(map(lambda x: '.fasta' in x, fid))
    for i, j in zip(pep_files, fasta_files):
        PepFile(i)
        FastaFile(j)
else:
    PepFile(args.pep_file)
    FastaFile(args.fasta_file)

for pep in PepFile:
    print(pep.title)
for fas in FastaFile:
    print(fas.title)

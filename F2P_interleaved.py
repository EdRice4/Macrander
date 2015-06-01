from os import getcwd, listdir
from random import randrange
import argparse


class IDSeq(object):

    """Generate two ordered lists; one containing each identifier and the
    other containing each corresponding sequence."""

    def generate_original_ids_and_seqs(self):
        for i in self.fasta_file:
            if i[0] == '>':
                self.identities.append(i.strip('>\r\n'))
            else:
                self.sequences.append(i.strip())


class OriginalPhylipFile(IDSeq):

    """Write phylip file in interleaved format containing original identities
       and their corresponding sequences."""

    def write_original_phylip_file_sequential(self):
        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j + '\n')

    def write_original_phylip_file_interleaved(self):
        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j[0:50] + '\n')
            phy.write('\n')
            for i in range(50, len(self.sequences[0]), 50):
                for j in self.sequences:
                    phy.write(j[i:i + 50] + '\n')
                phy.write('\n')


class UniqueID(OriginalPhylipFile):

    """Pair each original identity with randomly generated string of
       10 digits."""

    def generate_unique_ids(self):
        for i in self.identities:
            self.unique_identities[randrange(0, 9999999)] = i


class UniquePhylipFile(UniqueID):

    """Write phylip file in interleaved format containing original identities
       and their corresponding sequences."""

    def write_unique_phylip_file_sequential(self):
        with open(self.original_phylip_name, 'w') as phy:
            for i, j in zip(self.identities, self.sequences):
                phy.write(i + '\t' + j + '\n')

    def write_unique_phylip_file_interleaved(self):
        with open(self.unique_phylip_name, 'w') as phy:
            for i, j in zip(self.unique_identities.itervalues(),
                            self.sequences):
                phy.write(str(i) + '\t' + j[0:50] + '\n')
            phy.write('\n')
            for i in range(50, len(self.sequences[0]), 50):
                for j in self.sequences:
                    phy.write(j[i:i + 50] + '\n')
                phy.write('\n')


class DictionaryFile(UniquePhylipFile):

    """Write original and corresponding unique identities to tab delimited
       file."""

    def write_dictionary_file(self):
        with open(self.dictionary_name, 'w') as dfile:
            for key, value in self.unique_identities.iteritems():
                dfile.write(str(key) + '\t' + value + '\n')


class IterRegistry(type):

    """A metaclass to allow for iteration of instances of the FastaFile
       class."""

    def __iter__(cls):
        return iter(cls.registry)


class FastaFile(DictionaryFile):

    """A class in which all of the necessary parameters corresponding to each
       respective fasta file are stored."""

    __metaclass__ = IterRegistry
    registry = []

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


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('fasta_file', type=str, help=('Name of fasta '
                                                      'sequence to be '
                                                      'converted. Specify '
                                                      'None if running in '
                                                      'batch mode.')
                        )
arg_parser.add_argument('-s', '--sequential', help=('Write phylip file in '
                                                    'sequential format. By '
                                                    'default, script writes '
                                                    'in interleaved format.'),
                        action='store_true'
                        )
arg_parser.add_argument('-b', '--batch', help=('Run script in batch mode; '
                                               'i.e. convert all fasta '
                                               'files in directory to '
                                               'phylip.'),
                        action='store_true'
                        )
args = arg_parser.parse_args()

if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    fasta_files = filter(lambda x: '.fasta' in x, fid)
    for i in fasta_files:
        FastaFile(i)
else:
    FastaFile(args.fasta_file)

for fasta in FastaFile:
    fasta.generate_original_ids_and_seqs()
    fasta.generate_unique_ids()
    if args.sequential:
        print('sequential')
        fasta.write_original_phylip_file_sequential()
        fasta.write_unique_phylip_file_sequential()
    else:
        print('interleaved')
        fasta.write_original_phylip_file_interleaved()
        fasta.write_unique_phylip_file_interleaved()
    fasta.write_dictionary_file()

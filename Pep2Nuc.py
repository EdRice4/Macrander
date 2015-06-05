from os import getcwd, listdir
from string import maketrans
from re import search, sub
import argparse


class DataParse(object):

    """A class in which all data parsing functionality is stored."""

    def extract_params_from_pep(self):
        seq_of_int = {}
        pep = filter(lambda x: x[0] == '>', self.pep)
        for i in pep:
            data = search('(\w+):(\d+)-(\d+)\((\+|-)\)', i)
            diff = int(data.group(2)) - int(data.group(3))
            if diff < 0:
                start = int(data.group(2)) - 1
                end = start - (diff + 2)
            else:
                start = int(data.group(3)) - 1
                end = start + (diff - 2)
            if seq_of_int.get(data.group(1), False):
                seq_of_int[data.group(1)].extend([start, end, data.group(4)])
            else:
                seq_of_int[data.group(1)] = [start, end, data.group(4)]
        return seq_of_int

    def split_and_build_dict(self):
        fasta_dict = {}
        fasta = self.fas[1:].split('>')
        for i in fasta:
            seq_id = i.split(' ', 1)[0]
            seq = i.split('\n', 1)[-1]
            fasta_dict[seq_id.replace('\n', '')] = seq.replace('\n', '')
        return fasta_dict


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
                rev_compl_fasta_dict[i[0]] = map(lambda x: x[::-1], ffd[i[0]])
                rev_compl_fasta_dict[i[0]] = map(lambda x: x.translate(tbl),
                                                 ffd[i[0]])
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
    fasta_dict = f.split_and_build_dict()
    filt_fasta_dict = f.extract_pertinent_seq(fasta_dict, seqs_of_int)
    rev_compl_fasta_dict = f.rev_compl(filt_fasta_dict, seqs_of_int)
    pretty_dict = f.make_pretty(rev_compl_fasta_dict)
    f.write_dict(pretty_dict)

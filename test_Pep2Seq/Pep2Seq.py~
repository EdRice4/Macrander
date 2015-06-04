from os import getcwd, listdir
from string import maketrans
import re
import argparse


class DataParse(object):

    """A class in which all data parsing functionality is stored."""

    def extract_params_from_pep(self):
        seq_of_int = {}
        pep = filter(lambda x: x[0] == '>', self.pep)
        for i in pep:
            data = re.search('(\w+):(\d+)-(\d+)\((\+|-)\)', i)
            diff = int(data.group(2)) - int(data.group(3))
            if diff < 0:
                position = data.group(2) + ':' + str(int(data.group(2)) - diff)
            if diff > 0:
                position = data.group(3) + ':' + str(int(data.group(3)) + diff)
            seq_of_int[data.group(1)] = [position, data.group(4)]
        return seq_of_int

    def split_and_build_dict(self):
        fasta_dict = {}
        fasta = self.fas[1:].split('>')
        for i in fasta:
            seq_id = i.split(' ', 1)[0]
            seq = i.split('\n', 1)[-1]
            fasta_dict[seq_id] = seq.strip()
        return fasta_dict


class RegEx(DataParse):

    """A class in which all regular expression functionality is stored."""

    def compile_pattern(self, parameters):
        seq_ids = zip(*parameters)[0]
        pattern = map(lambda x: re.escape(x), seq_ids)
        pattern = re.compile('|'.join(pattern))
        return pattern


class ExtractData(RegEx):

    """A class in which all data extraction, i.e. the boundary between
       disparate file types, functionality is stored."""

    def extract_pertinent_seq(self, fasta_dict, seq_of_int, match_object):
        filt_fasta_dict = {}
        for i, j in zip(fasta_dict.iteritems(), seq_of_int.iteritems()):
            if bool(match_object.match(i[0])):
                filt_fasta_dict[i[0]] = i[1][int(j[i[0]][0])]
        return filt_fasta_dict

    def rev_compl(self, filt_fasta_dict, seq_of_int):
        tbl = maketrans('ATCG', 'TAGC')
        rev_compl_fasta_dict = {}
        for i, j in zip(filt_fasta_dict.iteritems(), seq_of_int.iteritems()):
            if j[i[0]][1] == '-':
                rev_compl_fasta_dict[i[0]] = (i[1][::-1]).translate(tbl)
        return rev_compl_fasta_dict


class FileIO(ExtractData):

    """A class in which all file input/output functionality is stored."""

    def write_dict(self, rev_compl_fasta_dict):
        with open(self.fas_new, 'w') as fas:
            for i, j in zip(rev_compl_fasta_dict.iterkeys(),
                            rev_compl_fasta_dict.itervalues()):
                fas.write(i + '\t' + j + '\n')


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
    seq_id_match_pattern = f.compile_pattern(seqs_of_int)
    filt_fasta_dict = f.extract_pertinent_seq(fasta_dict, seqs_of_int,
                                              seq_id_match_pattern)
    print(filt_fasta_dict)

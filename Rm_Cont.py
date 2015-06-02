from os import getcwd, listdir
import re
import argparse


class RegEx(object):

    """A class in which regular expression functionality is stored."""

    def strip_and_escape(self):
        rm_strip_esc = map(lambda x: re.escape(x.strip()), self.rm)
        return rm_strip_esc

    def compile_rm(self):
        pattern = self.strip_and_escape()
        pattern = re.compile('|'.join(pattern))
        return pattern

    def remove_cont(self):
        pattern = self.compile_rm()
        for i in self.rm:
            if bool(pattern.match(i)):
                self.rm.remove(i)


class FileIO(RegEx):

    """A class in which file input/output is stored."""

    def write_new_fasta(self):
        with open(self.new_fasta, 'w') as fas:
            fas.write('>')
            fas.write(self.rm.join('>'))


class IterRegistry(type):

    """A metaclass allowing for iterations over instances of the FastaFile
       class."""

    def __iter__(cls):
        return iter(cls.registry)


class FastaFile(FileIO):

    """A class in which all the necessary parameters corresponding to each
       respective fasta file."""

    def __init__(self, fasta_file, rem_file):
        with open(fasta_file, 'r') as fas:
            self.fasta = fas.read()
        self.new_fasta = fasta_file.replace('.fasta', '_new.fasta')
        with open(rem_file, 'r') as rem:
            self.rm = rem.readlines()
        self.registry.append(self)

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('fasta_file', type=str, help=('Name of fasta sequence '
                                                      'containing sequences '
                                                      'to be removed. Specify '
                                                      'None if running in '
                                                      'batch mode.')
                        )
arg_parser.add_argument('rem_file', type=str, help=('Name of file containing '
                                                    'identities of sequences '
                                                    'to be removed. Specify '
                                                    'None if running in batch '
                                                    'mode.')
                        )
arg_parser.add_argument('-b', '--batch', help=('Run script in batch mode. '
                                               'i.e. perform removal with all '
                                               'fasta sequences in directory '
                                               'and their corresponding txt '
                                               'files.'),
                        action='store_true'
                        )
args = arg_parser.parse_args()

if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    fasta_files = sorted(filter(lambda x: '.fasta' in x, fid))
    rem_files = sorted(filter(lambda x: '.txt' in x, fid))
    for i, j in zip(fasta_files, rem_files):
        FastaFile(i, j)
else:
    FastaFile(args.fasta_file, args.rem_file)

for fas in FastaFile:
    fas.remove_cont()
    fas.write_new_fasta()

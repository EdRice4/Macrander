from re import search
from os import getcwd, listdir
import argparse


# {{{ FileIO class
class FileIO(object):

    """ {{{ Docstrings
    A class in which all file input/output functionality is stored.

    Namely:

        1.) A GO.txt/TPM.txt file is read in the form of a list with each
            item (line) stripped of trailing white space.
        2.) A concatenated GO/TPM dictionary is written to a file in the form
            of:

                transcript_id\tgene_ontology\n
                Seq_ID\tGO:GO_ID^GO^data:TPM\tGO:GO_ID^GO^data:TPM\n

            See docstrings in DataParse for explanation of TPM values.
        3.) Files output by write_concat_GO_dicts are read into list, later to
            be concatenated into a massive string so that the GO IDs and
            corresponding TPM values can be extracted.
        4.) A concatenated GO/TPM dictionary is written to a file in the form
            of:

                GO_id\tCumulative_TPM\n
                GO_ID\tTPM
            See docstrings in DataParse for explanation of cumulative TPM.
    }}} """

    # {{{ read_GO_TPM_file
    def read_GO_TPM_file(self, *GO_TPM_file):

        """ {{{ Docstrings
        Reads GO and TPM files into dictionary.
        }}} """

        dictionary = {}
        for i in GO_TPM_file:
            with open(i, 'r') as f:
                f = f.readlines()[1:]
            f = map(lambda x: x.strip(), f)
            if '_GO.txt' in i:
                dictionary['GO'] = f
            else:
                dictionary['TPM'] = f
        return(dictionary['GO'], dictionary['TPM'])
    # }}}

    # {{{ write_concat_GO_dicts
    def write_concat_GO_dicts(self, *GO_dict):

        """ {{{ Docstrings
        Writes concatenated GO/TPM dictionary to file.
        }}} """

        for i, j in zip(self.IDs[0:3], GO_dict):
            with open(i, 'w') as f:
                f.write('transcript_id\tgene_ontology\n')
                for k, v in j.iteritems():
                    f.write(k + '\t' + '\t'.join(v) + '\n')
    # }}}

    # {{{ read_concat_file
    def read_concat_file(self):

        """ {{{ Docstrings
        Reads files output by write_concat_GO_dicts into list.
        }}} """

        file_list = []
        for i in self.IDs[0:3]:
            with open(i, 'r') as cf:
                cf = cf.read()
            file_list.append(cf)
        return file_list
    # }}}

    # {{{ write_cum_file
    def write_cum_file(self, cum_data):

        """ {{{ Docstrings
        Writes the cumulative GO/TPM data to file.
        }}} """

        with open(self.IDs[3], 'w') as cum:
            cum.write('GO_id\tCumulative_TPM\n')
            for k, v in cum_data.iteritems():
                cum.write(k + '\t' + str(v) + '\n')
# }}}


# {{{ DataParse class
class DataParse(FileIO):

    """ {{{ Docstrings
    A class in which all data parsing functionality is stored.

    Namely:

        1.) A GO dictionary is constructed in the form of:

                dictionary = {
                        'Seq_ID' : [
                                'GO:GO_ID^GO^data',
                                'GO:GO_ID^GO^data'
                                ]
                        }
            Alternatively, a TPM dictionary is constructed in the form of:

                dictionary = {'Seq_ID' : TPM}

        2.) Utilizing this dictionary, the values are filtered to construct
            novel dictionaries categorized by the primary GO category
            (cellular component, biological process, or molecular function).
            i.e. only GO hits corresponding to the pertinent primary GO
            category are included.
        3.) A concatenated GO/TPM dictionary is constructed in the form of:

                dictionary = {
                        'Seq_ID' : [
                                'GO:GO_ID^GO^data:TPM',
                                'GO:GO_ID^GO^data:TPM'
                                ]
                        }

            Where the TPM value is the original TPM value divided by the number
            of hits for the corresponding sequence for that primary GO
            category.
        4.) See 5.
        5.) A novel concatenated GO/TPM dictionary is constructed in the form
            of:

                dictionary = {'GO_ID' : cumulative TPM}

            Where the value of cumulative TPM is the sum off all TPM values
            across all primary GO categories for that corresponding GO_ID.
    }}} """

    # {{{ build_GO_TPM_dict
    def build_GO_TPM_dict(self, GO_TPM_file_list, Is_GO):

        """ {{{ Docstrings
        Constructs dictionary, seperating sequence IDs from data.
        }}} """

        f = map(lambda x: x.split('\t'), GO_TPM_file_list)
        if Is_GO:
            go_dict = dict((i[0], i[1].split('`')) for i in f)
            return go_dict
        else:
            tmp_dict = dict((i[0], int(i[1])) for i in f)
            return tmp_dict
    # }}}

    # {{{ filter_GO_dict
    def filter_GO_dict(self, GO_dict, *GO_category):

        """ {{{ Docstrings
        Filters GO data based on the primary GO category.
        }}} """

        dictionary = {}
        for i in GO_category:
            tmp_dict = {}
            for k, v in GO_dict.iteritems():
                tmp_value = filter(lambda x: i in x, v)
                if tmp_value:
                    tmp_dict[k] = tmp_value
            dictionary[i] = tmp_dict
                # {{{ return empty
                # uncomment this block of code if you would like script to
                # print 'empty' in gene_ontology column if no GO hits in
                # pertinent GO category
                # else:
                #     dictionary[k] = 'empty'
                # }}}
        return(
                dictionary['cellular_component'],
                dictionary['biological_process'],
                dictionary['molecular_function']
                )
        # }}}

    # {{{ concatenate_GO_TPM_data
    def concatenate_GO_TPM_data(self, TPM_dict, *filtered_GO_dicts):

        """ {{{ Docstrings
        Concatenates GO and TPM for subsequent writing to file.
        }}} """

        dictionary = {}
        for i in filtered_GO_dicts:
            tmp_dict = {}
            for k, v in i.iteritems():
                tmp_dict[k] = map(
                        lambda x: x + ':{0}'.format(TPM_dict[k] / len(v)), v
                        )
            if i == go_cc:
                dictionary['go_tpm_cc'] = tmp_dict
            elif i == go_bp:
                dictionary['go_tpm_bp'] = tmp_dict
            else:
                dictionary['go_tpm_mf'] = tmp_dict
        return(
                dictionary['go_tpm_cc'], dictionary['go_tpm_bp'],
                dictionary['go_tpm_mf']
                )
    # }}}

    # {{{ prep_cum_data
    def prep_cum_data(self, list_of_concat_files):

        """ {{{ Docstrings
        Filters and further concatenates files for subsequent parsing.
        }}} """

        cf = ''.join(list_of_concat_files)
        cf = cf.replace('\n', '\t')
        cf = cf.split('\t')
        cf = filter(lambda x: 'GO:' in x, cf)
        return cf
    # }}}

    # {{{ build_cum_dict
    def build_cum_dict(self, prepped_cum_data):

        """ {{{ Docstrings
        Utilizing re module, extracts pertinent data and constructs
        dictionary.
        }}} """

        dictionary = {}
        match_objects = map(
                lambda x: search('GO:(\d+).*:(\d+)', x), prepped_cum_data
                )
        for i in match_objects:
            if not dictionary.get(i.group(1), False):
                dictionary[i.group(1)] = int(i.group(2))
            else:
                dictionary[i.group(1)] = (
                        dictionary[i.group(1)] + int(i.group(2))
                        )
        return dictionary
    # }}}
# }}}


# {{{ IterRegistry metaclass
class IterRegistry(type):

    """ {{{ Docstrings
    A metaclass allowing for iterations over instances of the GO_TPM class.
    }}} """

    # {{{ __iter__
    def __iter__(cls):
        return iter(cls.registry)
    # }}}
# }}}


# {{{ # GO_TPM class
class GO_TPM(DataParse):

    # {{{ __metaclass__
    __metaclass__ = IterRegistry
    registry = []
    # }}}

    """ {{{ Docstrings
    A class in which all pertinent data and parameters correspnding to
    each each Seq and their respective GO and TMP files are stored.
    }}} """

    # {{{ __init__
    def __init__(self, seq, GO_file, TPM_file):
        self.go, self.tpm = self.read_GO_TPM_file(GO_file, TPM_file)
        self.IDs = (
                '{0}_CC{1},{0}_BP{1},{0}_MF{1},{0}_CUM{1}'.format(seq, '.txt')
                ).split(',')
        self.registry.append(self)
    # }}}
# }}}


# {{{ ArgParse
arg_parser = argparse.ArgumentParser(
        prog='GO+TPM',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
                'Concatencates gene ontology and quantification data from '
                'RNA-seq (TPM) in a more informative manner.'
                ),
        epilog=('Note: All GO and TPM files must contain \'_GO\' and \'_TPM\' '
                'in name, respectively, and must be named in analogous '
                'manner. e.g. Seq_One_GO.txt, Seq_One_TPM.txt'
                )
        )
arg_parser.add_argument(
        '-GO', type=str, help='Name of GO file to be concatenated.',
        default='None'
        )
arg_parser.add_argument(
        '-TPM', type=str, help='Name of TPM file to be concatenated.',
        default='None'
        )
arg_parser.add_argument(
        '-b', '--batch', help=('Run script in batch mode. i.e. Perform '
                               'concatenation with all respectie GO files and '
                               'their corresponding TPM files.'
                               ),
        action='store_true'
        )
args = arg_parser.parse_args()
# }}}


# {{{ Instantiate instances of GO_TPM class
if args.batch:
    cwd = getcwd()
    fid = listdir(cwd)
    GO_files = sorted(filter(lambda x: '_GO.txt' in x, fid))
    TPM_files = sorted(filter(lambda x: '_TPM.txt' in x, fid))
    for i, j in zip(GO_files, TPM_files):
        GO_TPM(i.replace('_GO.txt', ''), i, j)
else:
    GO_TPM(args.GO, args.TPM)
# }}}


# {{{ Run functions for instances of GO_TPM class
for data in GO_TPM:
    go = data.build_GO_TPM_dict(data.go, True)
    tpm = data.build_GO_TPM_dict(data.tpm, False)
    go_cc, go_bp, go_mf = data.filter_GO_dict(
            go, 'cellular_component', 'biological_process',
            'molecular_function'
            )
    go_tpm_cc, go_tpm_bp, go_tpm_mf = data.concatenate_GO_TPM_data(
            tpm, go_cc, go_bp, go_mf
            )
    data.write_concat_GO_dicts(go_tpm_cc, go_tpm_bp, go_tpm_mf)
    concat_file_list = data.read_concat_file()
    prepped_cum_data = data.prep_cum_data(concat_file_list)
    cum_data = data.build_cum_dict(prepped_cum_data)
    data.write_cum_file(cum_data)
# }}}

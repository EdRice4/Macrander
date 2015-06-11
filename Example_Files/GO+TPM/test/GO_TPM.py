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

                transcript_id\tgene_ontology
                Seq_ID\tGO:GO_ID^GO^data:TPM\tGO:GO_ID^GO^data:TPM\n

        3.) The concatenated GO/TPM files are read to be parsed for subsequent
            data.
    }}} """

    # {{{ read_GO_TPM_file
    def read_GO_TPM_file(self, *GO_TPM_file):
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
        for i, j in zip(self.IDs, GO_dict):
            with open(i, 'w') as f:
                f.write('transcript_id\tgene_ontology\n')
                for k, v in j.iteritems():
                    f.write(k + '\t' + '\t'.join(v) + '\n')
    # }}}

    # {{{ read_concat_file
    def read_concat_file(self):

        """ {{{ Docstrings
        Reads files output by write_concat_GO_dicts into list so that the GO
        IDs and corresponding TPM values can be extracted.
        and corresponding TPM values
        }}} """

        file_list = []
        for i in self.IDs:
            with open(i, 'r') as cf:
                cf = cf.read()
            file_list.append(cf)
        return file_list
    # }}}
# }}}


# {{{ DataParse class
class DataParse(FileIO):

    """ {{{ Docstrings
    A class in which all data parsing functionality is stored.

    Namely:

        1.) A dictionary with sequences IDs as keys and corresponding GO terms
            as values is constructed.
        2.) Utilizing this dictionary, the values are filtered to construct
            novel dictionaries categorized by the primary GO category.
    }}} """

    # {{{ build_GO_TPM_dict
    def build_GO_TPM_dict(self, GO_TPM_file_list, Is_GO):

        """ {{{ Docstrings
        Builds dictionary, seperating sequence IDs from data.
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
        Filters GO data based on GO category (e.g. Cellular Component,
        Biological Process, Molecular Function)
        }}} """

        filtered_GO_dict = {}
        for i in GO_category:
            for k, v in GO_dict.iteritems():
                tmp_value = filter(lambda x: GO_category in x, v)
                if tmp_value:
                    filtered_GO_dict[GO_category] = dict(k=tmp_value)
                # {{{ return empty
                # uncomment this block of code if you would like script to
                # print 'empty' in gene_ontology column if no GO hits in
                # pertinent GO category
                # else:
                #     filtered_GO_dict[k] = 'empty'
                # }}}
        return(
                filtered_GO_dict['cellular_component'],
                filtered_GO_dict['biological_process'],
                filtered_GO_dict['molecular_function']
                )
    # }}}

    # {{{ concatenate_GO_TPM_data
    def concatenate_GO_TPM_data(self, filtered_GO_dict, TPM_dict):
        GO_TPM_dict = {}
        for k, v in filtered_GO_dict.iteritems():
            GO_TPM_dict[k] = map(
                    lambda x: x + ':{0}'.format(TPM_dict[k] / len(v)), v
                    )
        return GO_TPM_dict
    # }}}

    # {{{ prep_concat_data
    def prep_concat_data(self, list_of_concat_files):

        """ {{{ Docstrings
        Prepares and further concatenates files for subsequent parsing.
        }}} """

        cf = ''.join(list_of_concat_files)
        cf = cf.repalce('\n', '\t')
        cf = cf.split('\t')
        cf = filter(lambda x: 'GO:' in x, cf)
        return cf
    # }}}

    # {{{ build_concat_dict
    def build_concat_dict(self, prepped_concat_data):

        """ {{{ Docstrings
        Utilizing re module, constrcuts dictionary in the form of:

            dictionary = {'GO_ID' : cumulative TPM value}

        }}} """

        dictionary = {}
        match_objects = map(
                lambda x: search('GO:(\d+).*:(\d+)', x), prepped_concat_data
                )
        for i in match_objects:
            if not dictionary.get(i.group(0), False):
                dictionary[i.group(0)] = int(i.group(2))
            else:
                dictionary[i.group(0)] = dictionary[i.group(0)] + i.group(2)
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
        #self.IDs = (
                #'{0}_CC{1},{0}_BP{1},{0}_MF{1}'.format(seq, '.txt')
        #        ).split(',')
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
    print(go_cc, go_bp, go_mf)
    #go_bp = data.filter_GO_dict(go, 'biological_process')
    #go_mf = data.filter_GO_dict(go, 'molecular_function')
    #go_tpm_cc = data.concatenate_GO_TPM_data(go_cc, tpm)
    #data.write_concat_GO_dicts(go_tpm_cc)
# }}}

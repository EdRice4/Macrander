from numpy import genfromtxt
from os import getcwd, listdir


class FileIO(object):

    """
    A class in which all file input/output functionality is stored.

    Namely:

        1.) A GO.txt/TPM.txt file is read in the form of a list with each
            item (line) stripped of trailing white space.
        2.) A filtered GO dictionary is written to a file in the form of:
            transcript_id\tgene_ontology.
    """

    def read_GO_TPM_file(self, GO_TPM_file):
        with open(GO_TPM_file, 'r') as f:
            f = f.readlines()[1:]
        f = map(lambda x: x.strip(), f)
        return f

    def write_filtered_GO_dict(self, GO_dict):
        with open(self.ID + '.txt', 'w') as filtered:
            filtered.write('transcript_id\tgene_ontology\n')
            for k, v in GO_dict.iteritems():
                filtered.write(
                        k + '\t' + '\t'.join(v) + '\n'
                        )


class DataParse(FileIO):

    """
    A class in which all data parsing functionality is stored.

    Namely:

        1.) A dictionary with sequences IDs as keys and corresponding GO terms
            as values is constructed.
        2.) Utilizing this dictionary, the values are filtered to construct
            novel dictionaries categorized by the primary GO category.
    """

    def build_GO_dict(self, GO_TPM_file_list, Is_GO):
        f = map(lambda x: x.split('\t'), GO_TPM_file_list)
        if Is_GO:
            go_dict = dict((i[0], i[1].split('`')) for i in f)
            return go_dict
        else:
            tmp_dict = dict((i[0], i[1]) for i in f)
            return tmp_dict

    def filter_GO_dict(self, GO_dict, GO_category):
        filtered_GO_dict = {}
        for k, v in GO_dict.iteritems():
            tmp_value = filter(lambda x: GO_category in x, v)
            if tmp_value:
                filtered_GO_dict[k] = tmp_value
            # uncomment this block of code if you would like script to print
            # 'empty' in gene_ontology column if no GO hits in pertinent
            # GO category
            # else:
            #     filtered_GO_dict[k] = 'empty'
        return filtered_GO_dict

    def concatenate_GO_TPM_data(self, filtered_GO_dict, TPM_dict):
        GO_TPM_dict = {}
        for k, v in filtered_GO_dict.iteritems():
            GO_TPM_dict[k] = [
                    map(lambda x: x + ':%s' % (TPM_dict[k] / len(x)), v)
                    ]
        return GO_TPM_dict

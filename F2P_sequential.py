import argparse

def convert_fasta_to_phylip(fasta_seq, name_of_new_phy_file):

    """Take a sequence file in fasta format and convert to phylip format
       without truncation of identifiers."""

    identities = []
    sequences = []
    with open(fasta_seq, 'r') as fasta:
        fasta = fasta.readlines()
    for i in fasta:
        if i[0] == '>':
            identities.append(i)
        else:
            sequences.append(i)
    with open(name_of_new_phy_file, 'w') as phy:
        for i, j in zip(identities, sequences):
            phy.write(i.strip('>\r\n') + '\t' + j)

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('fasta_sequence', type=str, help=('Name of fasta '
                        'sequence to be converted.'))
arg_parser.add_argument('phy_file', type=str, help=('Name of new phylip file.'))
args = arg_parser.parse_args()

convert_fasta_to_phylip(args.fasta_sequence, args.phy_file)

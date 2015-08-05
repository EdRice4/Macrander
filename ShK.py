#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code written by: Edwin Rice
# email: edwinricethe4th@gmail.com
# phone: +1 (513) 426-4187
# github: https://github.com/EdRice4
#
# Initiation date|time: 08/05/2015|18:16:26

# {{{ Imports
from re import match
import argparse
# }}}


# {{{ FastaFile
class FastaFile(object):

    """ {{{ Docstrings

    Class in which fasta file functionality is stored.

    }}} """

    # {{{ read_pep_into_dict
    def read_pep_into_dict(self, pep_file):

        """ {{{ Docstrings

        Reads the pep file into a dictionary in the format of:

            dictionary = {

                    'Pep_ID' : 'Pep_seq'

                    }

        }}} """

        with open(pep_file, 'r') as pep:
                pep = pep.readlines()
    # }}}
# }}}

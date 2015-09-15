import os
import logging
import unittest

from riboplot import ribocore, riboplot, ribocount

# use testing configuration
CFG = ribocount.CONFIG = riboplot.config.TestingConfig()
logging.disable(logging.CRITICAL)


class RiboCountTestCase(unittest.TestCase):

    def test_unrelated_fasta_file(self):
        """If an unrelated fasta file is used, raise an error"""
        parser = ribocount.create_parser()
        args = parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.UNRELATED_FASTA])
        self.assertRaises(ribocore.ArgumentError, ribocount.main, args)

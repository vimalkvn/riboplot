import os
import logging
import unittest

from riboplot import ribocore, riboplot, ribocount

# use testing configuration
CONFIG = ribocount.CONFIG = riboplot.config.TestingConfig()
logging.disable(logging.CRITICAL)

RIBO_FILE = os.path.join(CONFIG.TEST_DATA_DIR, '5hRPFsorted.bam')
RNA_FILE = os.path.join(CONFIG.TEST_DATA_DIR, '5hmRNAsorted.bam')
TRANSCRIPT_NAME = 'gi|148357119|ref|NM_001098396.1|'
TRANSCRIPTOME_FASTA = os.path.join(CONFIG.TEST_DATA_DIR, 'zebrafish.fna')
TRANSCRIPTOME_FASTA_MINUS1 = os.path.join(CONFIG.TEST_DATA_DIR, 'zebrafish_minus1.fna')
UNRELATED_FASTA = os.path.join(CONFIG.TEST_DATA_DIR, 'unrelated.fna')


class RiboCountTestCase(unittest.TestCase):

    def test_unrelated_fasta_file(self):
        """If an unrelated fasta file is used, raise an error"""
        parser = ribocount.create_parser()
        args = parser.parse_args(['-b', RIBO_FILE, '-f', UNRELATED_FASTA])
        self.assertRaises(ribocore.ArgumentError, ribocount.main, args)

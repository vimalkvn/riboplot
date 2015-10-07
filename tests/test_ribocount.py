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

    def test_get_ribo_counts(self):
        """Get read counts upstream of the longest ORF"""
        with ribocore.open_pysam_file(CFG.RIBO_FILE, ftype='bam') as f:
            counts, reads = ribocore.get_ribo_counts(ribo_fileobj=f, transcript_name=CFG.TRANSCRIPT_NAME)

            # 529 is the start position of the longest orf for the test transcript
            counts_upstream, reads_upstream = ribocore.filter_ribo_counts(counts=counts, orf_start=529)

        print '\nTotal read counts: {}\nUpstream read counts: {}'.format(reads, reads_upstream)
        self.assertTrue(reads_upstream > 1, msg='There should be reads upstream for this BAM file')
        self.assertTrue(reads > 1, msg='There should be Ribo-Seq reads for this BAM file')
        self.assertTrue(len(counts) > len(counts_upstream), msg='Total read counts should be higher than upstream read counts')

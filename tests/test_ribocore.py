# -*- coding: utf-8 -*-
import os
import logging
import unittest

from riboplot import config, ribocore

# use testing configuration
CFG = config.TestingConfig()
logging.disable(logging.CRITICAL)


# reference_orfs were obtained manually (UGENE, ORF finder)
REFERENCE_ORFS = [{'start': 45, 'stop': 203, 'sequence': 'ATGATTGAGAGGGACGGCCGGGGGCATTCGTATTGCGCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACGAAAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAAGAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAG'},
                  {'start': 219, 'stop': 254, 'sequence': 'ATGCCGACCCGCGATCCGGCGGCGTTATTCCCATGA'},
                  {'start': 251, 'stop': 328, 'sequence': 'ATGACCCGCCGGGCAGCGTGCGGGAAACCACGAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTGAAACTTAA'},
                  {'start': 306, 'stop': 374, 'sequence': 'ATGGTTGCAAAGCTGAAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAA'},
                  {'start': 465, 'stop': 512, 'sequence': 'ATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGTTCATTCCGATAA'},
                  {'start': 529, 'stop': 945, 'sequence': 'ATGCTAAATAGTTACGCGGCCCCGCGCGGTCGGCGTCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACGCGAGATGGAGCAATAACAGGTCTGTGATGCCCTTAGATGTCCGGGGCTGCACGCGCGCCACAATGGGCGGATCAACGTGTGCCTACCCTGCGCCGAGAGGCGCGGGTAACCCGTTGAACCCCGCTCGTGATTGGGACTGGGGCTTGAAACTGTTTCCCATCAACGAGGAATTCCCAGTAAGCGCAGGTCATAAGCTTGCGTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGTCCTCGGATCGGCCCCGCCGGGGCTCCTCGCCGGGCCCTGGCGGAGCGCCGAGAAGACGATCAAACTTGATCCTCTAG'},
                  {'start': 606, 'stop': 617, 'sequence': 'ATGGAGCAATAA'},
                  {'start': 638, 'stop': 718, 'sequence': 'ATGTCCGGGGCTGCACGCGCGCCACAATGGGCGGATCAACGTGTGCCTACCCTGCGCCGAGAGGCGCGGGTAACCCGTTGA'},
                  {'start': 854, 'stop': 862, 'sequence': 'ATGGTTTAG'}]


class BamTestCase(unittest.TestCase):
    """Check if all arguments sent on the command line are valid."""

    def test_is_bam_valid(self):
        """Test if BAM file is valid."""
        valid = ribocore.is_bam_valid(CFG.RIBO_FILE)
        self.assertTrue(valid)

        # test with a FASTA file (which is not BAM)
        self.assertRaises(ValueError, ribocore.is_bam_valid, CFG.TRANSCRIPTOME_FASTA)

    def test_bam_has_index(self):
        """Check if BAM file has an index."""
        # RPF file has an index
        has_index = ribocore.bam_has_index(CFG.RIBO_FILE)
        self.assertTrue(has_index)

        # RNA file doesn't have an index
        has_index = ribocore.bam_has_index(CFG.RNA_FILE)
        self.assertFalse(has_index)

    def test_create_bam_index(self):
        """Index a BAM file."""
        ribocore.create_bam_index(CFG.RNA_FILE)

        # check if index exists
        has_index = ribocore.bam_has_index(CFG.RNA_FILE)
        self.assertTrue(has_index)

        # remove index
        os.remove('{}.bai'.format(CFG.RNA_FILE))


class FastaTestCase(unittest.TestCase):

    def test_is_fasta_valid(self):
        """A valid FASTA file can be opened with pysam.FastaFile."""
        self.assertTrue(ribocore.is_fasta_valid(CFG.TRANSCRIPTOME_FASTA))

    def test_get_fasta_records(self):
        """Given a transcriptome fasta file and a transcript name, it should \
        be possible to get the sequence and length of a given transcript.

        """
        record = ribocore.get_fasta_records(CFG.TRANSCRIPTOME_FASTA, [CFG.TRANSCRIPT_NAME])
        self.assertEqual(len(record[CFG.TRANSCRIPT_NAME]['sequence']), CFG.TRANSCRIPT_LENGTH)

    def test_get_fasta_record(self):
        """Get a single FASTA record from a transcriptome."""
        record = ribocore.get_fasta_record(fasta_file=CFG.TRANSCRIPTOME_FASTA,
                                           transcript_name=CFG.TRANSCRIPT_NAME)
        self.assertEqual(record[CFG.TRANSCRIPT_NAME], CFG.TRANSCRIPT_SEQUENCE)
        self.assertEqual(len(record[CFG.TRANSCRIPT_NAME]), CFG.TRANSCRIPT_LENGTH)


class RiboCoreTestCase(unittest.TestCase):

    def test_get_three_frame_orfs(self):
        """Get ORFs in frames 1, 2 and 3."""
        orfs = ribocore.get_three_frame_orfs(sequence=CFG.TRANSCRIPT_SEQUENCE,
                                             starts=['ATG'], stops=['TAG', 'TAA', 'TGA'])

        # function should return the same ORFs as reference
        self.assertEqual(len(orfs), len(REFERENCE_ORFS))
        for item in orfs:
            self.assertTrue(item in REFERENCE_ORFS)

    def test_get_longest_orf(self):
        """Get the longest ORF from a list."""
        # longest ORF in reference is sequence with start 529
        longest_orf = ribocore.get_longest_orf(REFERENCE_ORFS)
        self.assertEqual(longest_orf['start'], 529)

import os
import shutil
import logging
import unittest
import tempfile

from riboplot import ribocore, riboplot

# use testing configuration
CFG = riboplot.CONFIG = riboplot.config.TestingConfig()
logging.disable(logging.CRITICAL)


class CheckArgumentsTestCase(unittest.TestCase):
    """Check if all arguments sent on the command line are valid."""
    parser = riboplot.create_parser()

    def test_bedtools_missing(self):
        """If bedtools is not in PATH, raise an error."""
        args = self.parser.parse_args(
            ['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME, '-n', CFG.RNA_FILE])
        save_path = os.environ['PATH']
        os.environ['PATH'] = ''
        self.assertRaises(OSError, ribocore.check_optional_arguments, ribo_file=args.ribo_file, rna_file=args.rna_file)
        os.environ['PATH'] = save_path

    def test_valid_read_length(self):
        """Read length should be a valid integer."""
        args = self.parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA,
                                       '-t', CFG.TRANSCRIPT_NAME, '-l', '28'])
        ribocore.check_optional_arguments(ribo_file=args.ribo_file, read_length=args.read_length)

    def test_invalid_read_length(self):
        """An error is raised if an invalid read length is used."""
        args = self.parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME,
                                       '-l', '-1'])  # invalid read length -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_length=args.read_length)

        args = self.parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME,
                                       '-l', '100'])  # invalid read length 100
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_length=args.read_length)

    def test_valid_read_offset(self):
        """Read offset should be positive."""
        args = self.parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME,
                                       '-s', '-1'])  # invalid read offset -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_offset=args.read_offset)

    def test_missing_transcript_in_fasta(self):
        """If a transcript is missing in FASTA, an error is raised."""
        args = self.parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME])  # invalid read offset -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_required_arguments,
                          args.ribo_file, args.transcriptome_fasta, 'hello')

    def test_missing_transcript_in_bam(self):
        """If a transcript is missing in BAM, an error is raised."""
        # testing with an unrelated BAM file
        args = self.parser.parse_args(['-b', CFG.UNRELATED_BAM, '-f', CFG.TRANSCRIPTOME_FASTA,
                                       '-t', CFG.TRANSCRIPT_NAME])
        self.assertRaises(ribocore.ArgumentError, ribocore.check_required_arguments, args.ribo_file,
                          args.transcriptome_fasta, args.transcript_name)


class RNACountsTestCase(unittest.TestCase):

    def test_get_rna_counts(self):
        """Test get RNA counts for transcript from RNA-Seq BAM file. Assumes bedtools is installed."""
        counts = riboplot.get_rna_counts(CFG.RNA_FILE, CFG.TRANSCRIPT_NAME)
        self.assertIsInstance(counts, dict)
        self.assertTrue(len(counts) > 0)

    def test_invalid_rna_file(self):
        """If an invalid RNA file is provided, generate an error message"""
        # using transcriptome FASTA file as the invalid RNA file for test
        parser = riboplot.create_parser()
        args = parser.parse_args(['-b', CFG.RIBO_FILE,  '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME, '-n', CFG.TRANSCRIPTOME_FASTA])
        self.assertRaises(ValueError, ribocore.check_optional_arguments, ribo_file=args.ribo_file, rna_file=args.rna_file)


class RiboPlotTestCase(unittest.TestCase):

    def test_get_codon_positions(self):
        """Get codon positions in all frames given a sequence."""
        # the positions on this sequence were calculated manually.
        fasta = ('AACCGGAGCACCCAGAGAAAACCCACGCAAACGCAGGGAGAATTTGCAAACTCCACACA'
                 'GAAATGCCAGCTGATCCAGCCGAGCCTCGAGTCAGCATCCTTGCTTGTTGGATGCCTGA'
                 'TTGCAGTTCAACTCCAAACTCAGTTGGACCAGCTGATCAGTG')
        codon_positions = riboplot.get_start_stops(fasta)
        expected = {1: {'starts': [], 'stops': []},
                    2: {'starts': [], 'stops': [71, 116, 152]},
                    3: {'starts': [63, 111], 'stops': []}}
        self.assertEqual(codon_positions, expected)

    def test_valid_riboplot_run(self):
        """A good riboplot run"""
        output_dir = tempfile.mkdtemp()
        print 'Output path is {}'.format(output_dir)
        parser = riboplot.create_parser()
        args = parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', CFG.TRANSCRIPT_NAME,
                                  '-o', output_dir])
        riboplot.main(args)
        for fname in ('riboplot.png', 'riboplot.svg', 'RiboCounts.csv'):
            self.assertTrue(os.path.exists(os.path.join(output_dir, fname)))
        shutil.rmtree(output_dir)

    def test_transcript_with_no_counts(self):
        """If the transcript has no ribocounts, no plot should be produced."""
        transcript = 'gi|62955616|ref|NM_001017822.1|'  # has no reads
        output_dir = tempfile.mkdtemp()
        parser = riboplot.create_parser()
        args = parser.parse_args(['-b', CFG.RIBO_FILE, '-f', CFG.TRANSCRIPTOME_FASTA, '-t', transcript, '-o', output_dir])
        self.assertRaises(ribocore.RiboPlotError, riboplot.main, args)
        for fname in ('riboplot.png', 'riboplot.svg', 'RiboCounts.csv'):
            self.assertFalse(os.path.exists(os.path.join(output_dir, fname)))
        shutil.rmtree(output_dir)

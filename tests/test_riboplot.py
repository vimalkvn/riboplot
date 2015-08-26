import os
import shutil
import logging
import unittest
import tempfile

from riboplot import ribocore, riboplot

# use testing configuration
CONFIG = riboplot.CONFIG = riboplot.config.TestingConfig()
logging.disable(logging.CRITICAL)

RIBO_FILE = os.path.join(CONFIG.TEST_DATA_DIR, '5hRPFsorted.bam')
RNA_FILE = os.path.join(CONFIG.TEST_DATA_DIR, '5hmRNAsorted.bam')
TRANSCRIPT_NAME = 'gi|148357119|ref|NM_001098396.1|'
TRANSCRIPTOME_FASTA = os.path.join(CONFIG.TEST_DATA_DIR, 'zebrafish.fna')
UNRELATED_BAM = os.path.join(CONFIG.TEST_DATA_DIR, 'unrelated.bam')


class CheckArgumentsTestCase(unittest.TestCase):
    """Check if all arguments sent on the command line are valid."""
    parser = riboplot.create_parser()

    def test_bedtools_missing(self):
        """If bedtools is not in PATH, raise an error."""
        args = self.parser.parse_args(
            ['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME, '-n', RNA_FILE])
        save_path = os.environ['PATH']
        os.environ['PATH'] = ''
        self.assertRaises(OSError, ribocore.check_optional_arguments, ribo_file=args.ribo_file, rna_file=args.rna_file)
        os.environ['PATH'] = save_path

    def test_is_bam_valid(self):
        """Test if BAM file is valid."""
        valid = ribocore.is_bam_valid(RIBO_FILE)
        self.assertTrue(valid)

        # test with a FASTA file (which is not BAM)
        self.assertRaises(ValueError, ribocore.is_bam_valid, TRANSCRIPTOME_FASTA)

    def test_bam_has_index(self):
        """Check if BAM file has an index."""
        # RPF file has an index
        has_index = ribocore.bam_has_index(RIBO_FILE)
        self.assertTrue(has_index)

        # RNA file doesn't have an index
        has_index = ribocore.bam_has_index(RNA_FILE)
        self.assertFalse(has_index)

    def test_create_bam_index(self):
        """Index a BAM file."""
        ribocore.create_bam_index(RNA_FILE)

        # check if index exists
        has_index = ribocore.bam_has_index(RNA_FILE)
        self.assertTrue(has_index)

        # remove index
        os.remove('{}.bai'.format(RNA_FILE))

    def test_valid_read_length(self):
        """Read length should be a valid integer."""
        args = self.parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA,
                                       '-t', TRANSCRIPT_NAME, '-l', '28'])
        ribocore.check_optional_arguments(ribo_file=args.ribo_file, read_length=args.read_length)

    def test_invalid_read_length(self):
        """An error is raised if an invalid read length is used."""
        args = self.parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME,
                                       '-l', '-1'])  # invalid read length -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_length=args.read_length)

        args = self.parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME,
                                       '-l', '100'])  # invalid read length 100
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_length=args.read_length)

    def test_valid_read_offset(self):
        """Read offset should be positive."""
        args = self.parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME,
                                       '-s', '-1'])  # invalid read offset -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_optional_arguments,
                          ribo_file=args.ribo_file, read_offset=args.read_offset)

    def test_is_fasta_valid(self):
        """A valid FASTA file can be opened with pysam.FastaFile."""
        self.assertTrue(ribocore.is_fasta_valid(TRANSCRIPTOME_FASTA))

    def test_missing_transcript_in_fasta(self):
        """If a transcript is missing in FASTA, an error is raised."""
        args = self.parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME])  # invalid read offset -1
        self.assertRaises(ribocore.ArgumentError, ribocore.check_required_arguments,
                          args.ribo_file, args.transcriptome_fasta, 'hello')

    def test_missing_transcript_in_bam(self):
        """If a transcript is missing in BAM, an error is raised."""
        # testing with an unrelated BAM file
        args = self.parser.parse_args(['-b', UNRELATED_BAM, '-f', TRANSCRIPTOME_FASTA,
                                       '-t', TRANSCRIPT_NAME])
        self.assertRaises(ribocore.ArgumentError, ribocore.check_required_arguments, args.ribo_file,
                          args.transcriptome_fasta, args.transcript_name)


class RNACountsTestCase(unittest.TestCase):

    def test_get_rna_counts(self):
        """Test get RNA counts for transcript from RNA-Seq BAM file. Assumes bedtools is installed."""
        counts = riboplot.get_rna_counts(RNA_FILE, TRANSCRIPT_NAME)
        self.assertIsInstance(counts, dict)
        self.assertTrue(len(counts) > 0)

    def test_invalid_rna_file(self):
        """If an invalid RNA file is provided, generate an error message"""
        # using transcriptome FASTA file as the invalid RNA file for test
        parser = riboplot.create_parser()
        args = parser.parse_args(['-b', RIBO_FILE,  '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME, '-n', TRANSCRIPTOME_FASTA])
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
        args = parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', TRANSCRIPT_NAME,
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
        args = parser.parse_args(['-b', RIBO_FILE, '-f', TRANSCRIPTOME_FASTA, '-t', transcript, '-o', output_dir])
        self.assertRaises(ribocore.RiboPlotError, riboplot.main, args)
        for fname in ('riboplot.png', 'riboplot.svg', 'RiboCounts.csv'):
            self.assertFalse(os.path.exists(os.path.join(output_dir, fname)))
        shutil.rmtree(output_dir)

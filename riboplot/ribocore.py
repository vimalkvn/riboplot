"""Common functions. """
import pysam
import logging
import subprocess


# create logger for the entire program
log = logging.getLogger('riboplot')
log.setLevel(logging.DEBUG)

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(module)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
ch.setFormatter(formatter)
log.addHandler(ch)


class ArgumentError(Exception):
    """Raised when invalid arguments are sent in the command line."""
    pass


class RiboPlotError(Exception):
    """General errors relating to riboplot."""
    pass


class RiboCountError(Exception):
    """General errors relating to ribocount."""
    pass


class RNACountsError(Exception):
    """For errors related to RNA Coverage generation using bedtools. """
    pass


def is_bam_valid(bam_file):
    """Check if bam file is valid. Raises a ValueError if pysam cannot read the file.

    #TODO: pysam  does not differentiate between BAM and SAM
    """
    try:
        f = pysam.AlignmentFile(bam_file)
    except ValueError:
        raise
    except:
        raise
    else:
        f.close()

    return True


def bam_has_index(bam_file):
    """Check if bam file has an index. Returns True/False."""
    has_index = None
    with pysam.AlignmentFile(bam_file, 'rb') as bam_fileobj:
        try:
            bam_fileobj.fetch(bam_fileobj.references[0])
        except ValueError as err:
            if err.message == 'fetch called on bamfile without index':
                bam_fileobj.close()
                has_index = False
        else:
            has_index = True
    return has_index


def create_bam_index(bam_file):
    """Create an index for the given BAM file."""
    pysam.index(bam_file)


def get_bam_read_lengths(bam_file):
    """For a given BAM file, return an unique list of read lengths present."""
    read_lengths = []
    with pysam.AlignmentFile(bam_file, 'rb') as bam_fileobj:
        for record in bam_fileobj:
            read_lengths.append(record.query_length)
    return set(read_lengths)


def is_fasta_valid(fasta_file):
    """Check if fasta file is valid. Raises a ValueError if pysam cannot read the file.

    #TODO: pysam  does not differentiate between BAM and SAM
    """
    try:
        f = pysam.FastaFile(fasta_file)
    except IOError:
        raise
    else:
        f.close()
    return True


def get_fasta_record(fasta):
    """Return a single transcript name from a valid fasta file"""
    f = pysam.FastaFile(fasta)
    transcript = f.references[0]
    f.close()
    return transcript


def get_fasta_records(fasta, transcripts):
    """Return list of transcript records from the given fasta file.
    Each record will be of the form {'sequence_id': {'sequence': 'AAA', 'length': 3}}

    trascripts should be provided as a list of sequence id's.

    """
    records = {}
    f = pysam.FastaFile(fasta)
    for transcript in transcripts:
        try:
            sequence, length = f.fetch(transcript), f.get_reference_length(transcript)
        except KeyError:
            msg = 'Transcript "{}" does not exist in transcriptome FASTA file'.format(transcript)
            log.error(msg)
            raise ArgumentError(msg)
        records[transcript] = {'sequence': sequence, 'length': length}
    f.close()
    return records


def get_ribo_counts(ribo_fileobj, transcript_name, read_length=0):
    """Get total reads and read counts in all 3 frames for the give transcript
    from input BAM file (indexed).

    """
    read_counts = {}
    total_reads = 0

    for record in ribo_fileobj.fetch(transcript_name):
        if read_length:
            if record.query_length != read_length:
                continue
        total_reads += 1
        position = record.pos + 1

        try:
            read_counts[position]
        except KeyError:
            read_counts[position] = {1: 0, 2: 0, 3: 0}

        rem = position % 3
        if rem == 0:
            read_counts[position][3] += 1
        else:
            read_counts[position][rem] += 1

    return read_counts, total_reads


def check_required_arguments(ribo_file, transcriptome_fasta, transcript_name=None):
    """Check required arguments of both riboplot and ribocount."""
    # Is this a valid BAM file? i.e., can pysam read it?
    try:
        is_bam_valid(ribo_file)
    except ValueError:
        log.error('The given RiboSeq BAM file is not valid')
        raise

    # Does the BAM file have an index? If not, create it.
    if not bam_has_index(ribo_file):
        log.info('Creating an index for the BAM file...')
        create_bam_index(ribo_file)

    # Is FASTA file valid?
    fasta_valid = False
    try:
        fasta_valid = is_fasta_valid(transcriptome_fasta)
    except IOError:
        log.error('Transcriptome FASTA file is not valid')
        raise

    if fasta_valid:
        if transcript_name:
            try:
                get_fasta_records(transcriptome_fasta, [transcript_name])
            except IOError:
                log.error('Could not get FASTA sequence of "{}" from transcriptome FASTA file'.format(transcript_name))
                raise
        else:
            # ribocount doesn't have a transcript option so we get the first
            # sequence name from the fasta file
            transcript_name =  get_fasta_record(transcriptome_fasta)

        # check if transcript also exists in BAM
        with pysam.AlignmentFile(ribo_file, 'rb') as bam_file:
            if transcript_name not in bam_file.references:
                msg = 'Transcript "{}" does not exist in BAM file'.format(transcript_name)
                log.error(msg)
                raise ArgumentError(msg)


def check_optional_arguments(ribo_file, read_length=None, read_offset=None, rna_file=None):
    """Check all optional arguments."""
    if rna_file:
        try:
            subprocess.check_output(['bedtools', '--version'])
        except OSError:
            log.error('Could not find bedtools in PATH. bedtools is required'
                      'for generating RNA coverage plot.')
            raise
        # Is this a valid BAM file? i.e., can pysam read it?
        try:
            is_bam_valid(rna_file)
        except ValueError:
            log.error('The given RNASeq BAM file is not valid')
            raise

    # If read_length is given, it must be a positive integer or reads of that
    # length must exist in the BAM file
    if read_length:
        if read_length < 0:
            msg = 'Read length must be a positive value'
            log.error(msg)
            raise ArgumentError(msg)

        bam_read_lengths = get_bam_read_lengths(ribo_file)
        if read_length not in bam_read_lengths:
            msg = 'Reads of the length "{}" does not exist in the BAM file'.format(read_length)
            log.error(msg)
            raise ArgumentError(msg)

    # If read_offset is given, it must be a positive integer
    if read_offset:
        if read_offset < 0:
            msg = 'Read offset must be 0 or greater'
            log.error(msg)
            raise ArgumentError(msg)

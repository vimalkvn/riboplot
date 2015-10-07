"""Common functions. """
import pysam
import logging
import subprocess
from contextlib import contextmanager


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


class BamFileError(Exception):
    """Errors related to BAM file"""
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


@contextmanager
def open_pysam_file(fname, ftype):
    """Open a BAM or FASTA file with pysam (for use with "with" statement)"""
    try:
        if ftype == 'bam':
            fpysam = pysam.AlignmentFile(fname, 'rb')
        elif ftype == 'fasta':
            fpysam = pysam.FastaFile(fname)
        yield fpysam
    except:
        raise
    else:
        fpysam.close()


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


def get_first_transcript_name(fasta_file):
    """Return the first FASTA sequence from the given FASTA file.

    Keyword arguments:
    fasta_file -- FASTA format file of the transcriptome

    """
    with open_pysam_file(fname=fasta_file, ftype='fasta') as f:
        transcript_name = f.references[0]
    return transcript_name


def get_fasta_record(fasta_file, transcript_name):
    """Return a single transcript from a valid fasta file as a record.

    record[transcript_name] = sequence

    Keyword arguments:
    fasta_file -- FASTA format file of the transcriptome
    transcript_name -- Name of the transcript as in the FASTA header

    """
    with open_pysam_file(fname=fasta_file, ftype='fasta') as f:
        sequence = f.fetch(transcript_name)
    return {transcript_name: sequence}


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


def get_three_frame_orfs(sequence, starts=None, stops=None):
    """Find ORF's in frames 1, 2 and 3 for the given sequence.

    Positions returned are 1-based (not 0)

    Return format [{'start': start_position, 'stop': stop_position, 'sequence': sequence}, ]

    Keyword arguments:
    sequence -- sequence for the transcript
    starts -- List of codons to be considered as start (Default: ['ATG'])
    stops -- List of codons to be considered as stop (Default: ['TAG', 'TGA', 'TAA'])

    """
    if not starts:
        starts = ['ATG']

    if not stops:
        stops = ['TAG', 'TGA', 'TAA']

    # Find ORFs in 3 frames
    orfs = []
    for frame in range(3):
        start_codon = None
        orf = ''
        for position in range(frame, len(sequence), 3):
            codon = sequence[position:position + 3]
            if codon in starts:
                # We have found a start already, so add codon to orf and
                # continue. This is an internal MET
                if start_codon is not None:
                    orf += codon
                    continue

                # New orf start
                start_codon = position
                orf = codon
            else:
                # if sequence starts with ATG, start_codon will be 0
                if start_codon is None:
                    # We haven't found a start codon yet
                    continue
                orf += codon
                if codon in stops:
                    # orfs[start_codon + 1] = orf
                    orfs.append({'start': start_codon + 1, 'stop': position + 3, 'sequence': orf})

                    # Reset
                    start_codon = None
                    orf = ''
    return orfs


def get_longest_orf(orfs):
    """Find longest ORF from the given list of ORFs."""
    sorted_orf = sorted(orfs, key=lambda x: len(x['sequence']), reverse=True)[0]
    return sorted_orf


def filter_ribo_counts(counts, orf_start=None, orf_stop=None):
    """Filter read counts  and return only upstream of orf_start or downstream
    of orf_stop.

    Keyword arguments:
    counts -- Ribo-Seq read counts obtained from get_ribo_counts.
    orf_start -- Start position of the longest ORF.
    orf_stop -- Stop position of the longest ORF.

    """
    filtered_counts = dict.copy(counts)
    for position in counts:
        if orf_start and orf_stop:
            # if only upstream and downstream reads are required, check if
            # current position is upstream or downstream of the ORF start/stop
            # if not, remove from counts
            if (position > orf_start and position < orf_stop):
                filtered_counts.pop(position)
        elif orf_start:
            # check if current position is upstream of ORF start. if not, remove
            if position > orf_start:
                filtered_counts.pop(position)
        elif orf_stop:
            # check if current position is downstream of ORF stop. If not,
            # remove
            if position < orf_stop:
                filtered_counts.pop(position)

    # calculate total reads for this transcript
    total_reads = sum(sum(item.values()) for item in filtered_counts.values())
    return filtered_counts, total_reads


def get_ribo_counts(ribo_fileobj, transcript_name, read_length=0):
    """Get read counts in 3 frames for the given transcript from input BAM file
    (indexed).

    Keyword arguments:
    ribo_fileobj -- file object - BAM file opened using pysam AlignmentFile
    transcript_name -- Name of transcript to get counts for
    read_length (optional) -- If provided, get counts only for reads of this length.

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

        if not bam_has_index(ribo_file):
            msg = ('Could not create an index for this BAM file. Is this a valid BAM file '
                   'and/or is the BAM file sorted by chromosomal coordinates?')
            log.error(msg)
            raise BamFileError(msg)

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
            transcript_name = get_first_transcript_name(transcriptome_fasta)

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

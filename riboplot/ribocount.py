# -*- coding: utf-8 -*-
"""ribocount"""
import sys

# check dependencies
try:
    import pysam
except ImportError as e:
    sys.exit('Could not import the "pysam" module\n\nImporting failed with'
             '{}\n\n'.format(e))

import os
# import timeit
import shutil
import zipfile
import logging
import argparse

import ribocore
import config

# Default is production
CONFIG = config.ProductionConfig()

log = logging.getLogger('riboplot')


class ErrorLogFormatter(logging.Formatter):
    """Custom error log format for the HTML file"""

    def format(self, record):
        return '<h2>RiboCount Error</h2><p>{}</p>'.format(record.msg)


def create_parser():
    """Argument parser. """
    parser = argparse.ArgumentParser(
        prog='ribocount.py', description='Output read counts for all transcripts')

    # required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-b', '--ribo_file', help='Ribo-Seq alignment file in BAM format', required=True)
    required.add_argument('-f', '--transcriptome_fasta', help='FASTA format file of the transcriptome', required=True)

    # optional arguments
    parser.add_argument('-l', '--read_length', help='Read length to consider (default: %(default)s)',
                        metavar='INTEGER', type=int)
    parser.add_argument('-s', '--read_offset', help='Read offset (default: %(default)s)',
                        metavar='INTEGER', type=int, default=0)
    parser.add_argument('-m', '--html_file', help='Output file for results (HTML)', default='ribocount.html')
    parser.add_argument('-o', '--output_path', help='Files are saved in this directory', default='output')
    parser.add_argument('-d', '--debug', help='Flag. Produce debug output', action='store_true')

    return parser


def main(args):
    """Main program"""
    (ribo_file, transcriptome_fasta, read_length, read_offset, output_path, html_file) = (
        args.ribo_file, args.transcriptome_fasta, args.read_length,
        args.read_offset, args.output_path, args.html_file)

    log.debug('Supplied arguments\n{}'.format(
        '\n'.join(['{:<20}: {}'.format(k, v) for k, v in vars(args).items()])))

    # error messages (simple format) are written to html file
    fh = logging.FileHandler(html_file)
    fh.setLevel(logging.ERROR)
    fh.setFormatter(ErrorLogFormatter('%(message)s'))
    log.addHandler(fh)

    log.info('Checking if required arguments are valid...')
    ribocore.check_required_arguments(ribo_file=ribo_file, transcriptome_fasta=transcriptome_fasta)
    log.info('Done')

    log.info('Checking if optional arguments are valid...')
    ribocore.check_optional_arguments(ribo_file=ribo_file, read_length=read_length, read_offset=read_offset)
    log.info('Done')

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # zip_dir contents will be written here and a zip archive will be created
    # from this directory
    zip_dir = os.path.join(output_path, 'ribocount_output')
    if not os.path.exists(zip_dir):
        os.mkdir(zip_dir)

    csv_dir = os.path.join(zip_dir, 'csv')
    if not os.path.exists(csv_dir):
        os.mkdir(csv_dir)

    count = 0
    table_content = ''
    bam_fileobj = pysam.AlignmentFile(ribo_file, 'rb')
    fasta_file = pysam.FastaFile(transcriptome_fasta)
    log.info('Get RiboSeq read counts for all transcripts in FASTA')
    for transcript in fasta_file.references:
        rp_counts, rp_reads = ribocore.get_ribo_counts(bam_fileobj, transcript, read_length)
        if not rp_reads:  # no reads for this transcript. skip.
            continue

        log.debug('Writing counts for {}'.format(transcript))
        count += 1
        csv_file = 'RiboCounts{}.csv'.format(count)
        with open(os.path.join(csv_dir, csv_file), 'w') as f:
            f.write('"Position","Frame 1","Frame 2","Frame 3"\n')

            for pos in range(1, len(fasta_file[transcript]) + 1):
                if pos in rp_counts:
                    f.write('{0},{1},{2},{3}\n'.format(
                        pos, rp_counts[pos][1], rp_counts[pos][2], rp_counts[pos][3]))
                else:
                    f.write('{0},{1},{2},{3}\n'.format(pos, 0, 0, 0))
        table_content += ('<tr><td>{0}</td><td>{1}</td><td>{2}</td><td><a href="csv/{3}">'
                          '{3}</a></td></tr>'.format(count, transcript, rp_reads, csv_file))
    fasta_file.close()
    bam_fileobj.close()

    # only for display in HTML
    if not read_length:
        read_length = 'All'

    log.info('Done')

    if not count:
        if read_length:
            log.info('No transcripts found for read length {}'.format(read_length))
        else:
            log.info('No transcripts found')
    else:
        with open(os.path.join(CONFIG.PKG_DATA_DIR, 'ribocount.html')) as g,\
                open(os.path.join(zip_dir, 'index.html'), 'w') as h:
            h.write(g.read().format(count=count, length=read_length, table_content=table_content))

        for asset in ('css', 'js'):
            asset_dir = os.path.join(zip_dir, asset)
            if not os.path.exists(asset_dir):
                os.mkdir(asset_dir)
            asset_data_dir = os.path.join(CONFIG.PKG_DATA_DIR, asset)
            for fname in os.listdir(asset_data_dir):
                shutil.copy(os.path.join(asset_data_dir, fname),
                            os.path.join(zip_dir, asset, fname))

        log.info('Creating zip file')
        os.chdir(output_path)
        with zipfile.ZipFile('ribocount_output.zip', 'w') as zipf:
            for root, d, f in os.walk('ribocount_output'):
                for name in f:
                    zipf.write(os.path.join(root, name))
        shutil.rmtree('ribocount_output')
        os.chdir('../')
        log.debug('Writing HTML report')

        with open(os.path.join(CONFIG.PKG_DATA_DIR, 'ribocount_index.html')) as j, open(args.html_file, 'w') as k:
            k.write(j.read().format(count=count, read_length=read_length))
    log.info('Finished')


def run():
    """Run program"""
    parsed = create_parser()
    args = parsed.parse_args()
    main(args)


if __name__ == '__main__':
    run()

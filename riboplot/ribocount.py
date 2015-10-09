# -*- coding: utf-8 -*-
"""ribocount"""
import os
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

    count_group = parser.add_mutually_exclusive_group()
    count_group.add_argument('-v', '--count_five', help='Flag. Output reads in 5\' region', action='store_true')
    count_group.add_argument('-r', '--count_three', help='Flag. Output reads in 3\' region', action='store_true')
    parser.add_argument('-m', '--html_file', help='Output file for results (HTML)', default='ribocount.html')
    parser.add_argument('-o', '--output_path', help='Files are saved in this directory', default='output')
    parser.add_argument('-d', '--debug', help='Flag. Produce debug output', action='store_true')

    return parser


def main(args):
    """Main program"""
    (ribo_file, transcriptome_fasta, read_length, read_offset, count_five, count_three,
     output_path, html_file) = \
        (args.ribo_file, args.transcriptome_fasta, args.read_length, args.read_offset,
         args.count_five, args.count_three, args.output_path, args.html_file)

    log.debug('Supplied arguments\n{}'.format(
        '\n'.join(['{:<20}: {}'.format(k, v) for k, v in vars(args).items()])))

    # error messages (simple format) are written to html file
    fh = logging.FileHandler(html_file)
    fh.setLevel(logging.ERROR)
    fh.setFormatter(ErrorLogFormatter('%(message)s'))
    log.addHandler(fh)

    log.info('Checking if provided arguments are valid...')
    ribocore.check_required_arguments(ribo_file=ribo_file, transcriptome_fasta=transcriptome_fasta)
    ribocore.check_optional_arguments(ribo_file=ribo_file, read_length=read_length, read_offset=read_offset)
    log.info('Done')

    with ribocore.open_pysam_file(fname=ribo_file, ftype='bam') as b, ribocore.open_pysam_file(fname=transcriptome_fasta, ftype='fasta') as f:
        # Total valid transcript count (ones with reads)
        count = 0
        prime = None
        table_body = ''  # HTML table body content
        if count_five:
            log.info('Only 5\' read counts requested')
            prime = '5'
        elif count_three:
            log.info('Only 3\' read counts requested')
            prime = '3'

        # create output directories
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

        log.info('Get RiboSeq read counts for all transcripts in FASTA')
        for transcript in f.references:
            ribo_counts, ribo_reads = ribocore.get_ribo_counts(b, transcript, read_length)
            if not ribo_reads:  # no reads for this transcript. skip.
                continue

            # By default, all counts will be written (ribo_counts)
            # If 5' or 3' counts requested, filter and use
            # those counts for printing instead
            write_counts = ribo_counts
            log.debug('Total read counts {}'.format(ribo_reads))

            # find longest ORF and filter counts based on whether 5' or 3' is
            # requested
            longest_orf = {}
            if count_five or count_three:
                # use default start and stop codons and find ORFs in all 3
                # frames (+)
                orfs = ribocore.get_three_frame_orfs(sequence=f[transcript])
                if not len(orfs):
                    log.debug('No ORFs for transcript {0}'.format(transcript))
                    continue
                longest_orf = ribocore.get_longest_orf(orfs=orfs)
                orf_start, orf_stop = longest_orf['start'], longest_orf['stop']
                log.info('Transcript: {0} Longest ORF Start: {1}, Stop: {2}'.format(transcript, orf_start, orf_stop))

                if count_five:
                    write_counts, five_reads = ribocore.filter_ribo_counts(counts=ribo_counts, orf_start=orf_start)
                    log.debug('5\' region read counts: {}'.format(five_reads))
                elif count_three:
                    write_counts, three_reads = ribocore.filter_ribo_counts(counts=ribo_counts, orf_stop=orf_stop)
                    log.debug('3\' region read counts: {}'.format(three_reads))

            if not len(write_counts):
                # no counts for transcript
                continue

            log.debug('Writing counts to CSV file for transcript {}'.format(transcript))
            count += 1
            csv_file = 'RiboCounts{}.csv'.format(count)
            with open(os.path.join(csv_dir, csv_file), 'w') as cw:
                cw.write('"Position","Frame 1","Frame 2","Frame 3"\n')
                for pos in range(1, len(f[transcript]) + 1):
                    if pos in write_counts:
                        cw.write('{0},{1},{2},{3}\n'.format(
                            pos, write_counts[pos][1], write_counts[pos][2], write_counts[pos][3]))
                    else:
                        cw.write('{0},{1},{2},{3}\n'.format(pos, 0, 0, 0))
            # HTML table
            table_body += '<tr><td>{0}</td><td>{1}</td>'.format(transcript, ribo_reads)
            if count_five:
                table_body += '<td>{0}</td>'.format(five_reads)
            elif count_three:
                table_body += '<td>{0}</td>'.format(three_reads)
            table_body += '<td><a href="csv/{0}">{0}</a></td></tr>'.format(csv_file)
        table_body += '</tbody>'

    # only for display in HTML
    if not read_length:
        read_length = 'All'

    if not count:
        if read_length:
            log.info('No transcripts found for read length {}'.format(read_length))
        else:
            log.info('No transcripts found')
    else:
        if prime:
            template = 'ribocount_prime.html'
        else:
            template = 'ribocount.html'
        with open(os.path.join(CONFIG.PKG_DATA_DIR, template)) as g,\
                open(os.path.join(zip_dir, 'index.html'), 'w') as h:
            h.write(g.read().format(count=count, length=read_length, prime=prime, table_body=table_body))

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

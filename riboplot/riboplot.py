# -*- coding: utf-8 -*-
import sys

# check dependencies
module_errors = {}
try:
    import pysam
except ImportError as e:
    module_errors['pysam'] = str(e)

try:
    from matplotlib import pyplot as plt
except ImportError as e:
    module_errors['matplotlib'] = str(e)
else:
    from matplotlib import gridspec
    from matplotlib.font_manager import FontProperties

if len(module_errors):
    msg = 'Could not import all required modules:\n\n'
    for module, error in module_errors.items():
        msg += 'Importing "{0}" failed with \n{1}\n\n'.format(module, error)
    sys.exit(msg)

import os
import re
import config
# import timeit
import shutil
import logging
import tempfile
import argparse
import subprocess

import ribocore
# Default is production
CONFIG = config.ProductionConfig()

# create logger
log = logging.getLogger('riboplot')


class ErrorLogFormatter(logging.Formatter):
    """Custom error log format for the HTML file"""

    def format(self, record):
        return '<h2>RiboPlot Error</h2><p>{}</p>'.format(record.msg)


def get_start_stops(transcript_sequence, start_codons=None, stop_codons=None):
        """Return start and stop positions for all frames in the given
        transcript.

        """
        if not start_codons:
            start_codons = ['ATG']
        if not stop_codons:
            stop_codons = ['TAA', 'TAG', 'TGA']

        seq_frames = {1: {'starts': [], 'stops': []},
                      2: {'starts': [], 'stops': []},
                      3: {'starts': [], 'stops': []}}

        for codons, positions in ((start_codons, 'starts'),
                                  (stop_codons, 'stops')):
            if len(codons) > 1:
                pat = re.compile('|'.join(codons))
            else:
                pat = re.compile(codons[0])

            for m in re.finditer(pat, transcript_sequence):
                # Increment position by 1, Frame 1 starts at position 1 not 0
                start = m.start() + 1
                rem = start % 3
                if rem == 1:  # frame 1
                    seq_frames[1][positions].append(start)
                elif rem == 2:  # frame 2
                    seq_frames[2][positions].append(start)
                elif rem == 0:  # frame 3
                    seq_frames[3][positions].append(start)
        return seq_frames


def get_rna_counts(rna_file, transcript_name):
    """Get coverage for a given RNA BAM file, return read counts. """
    # check if the RNA file exists
    if not os.path.exists(rna_file):
        msg = 'RNA-Seq BAM file "{}" does not exist'.format(rna_file)
        logging.error(msg)
        raise OSError(msg)
    rna_counts = {}

    cov_file = tempfile.NamedTemporaryFile(delete=False)
    try:
        subprocess.check_call(
            ['bedtools', 'genomecov', '-ibam', rna_file,
             '-bg'], stdout=cov_file)
    except subprocess.CalledProcessError as e:
        # needs testing
        raise ribocore.RNACountsError('Could not generate coverage for RNA BAM file. \n{}\n'.format(e))
    for line in open(cov_file.name):
        line = line.split()
        if line[0] == transcript_name:
            position, count = int(line[1]) + 1, int(line[3])
            rna_counts[position] = count
    cov_file.close()
    os.unlink(cov_file.name)
    return rna_counts


def set_axis_color(axis, color):
    """Sets the spine color of all sides of an axis (top, right, bottom, left)."""
    axis.spines['top'].set_color(color)
    axis.spines['right'].set_color(color)
    axis.spines['bottom'].set_color(color)
    axis.spines['left'].set_color(color)


def plot_profile(ribo_counts, transcript_name, transcript_length,
                 start_stops, read_length=None, read_offset=None, rna_counts=None,
                 html_file='index.html', output_path='output'):
    """Plot read counts (in all 3 frames) and RNA coverage if provided for a
    single transcript.

    """
    gs = gridspec.GridSpec(6, 1, height_ratios=[0.1, 6.8, 0.1, 1, 1, 1], hspace=0.35)
    font_xsmall = {'family': 'sans-serif', 'color': '#555555', 'weight': 'normal', 'size': 'x-small'}

    # plot for frames legend
    ax1 = plt.subplot(gs[0], axisbg='white')
    ax1.text(0.95, 0.1, "frame 1", size=6, ha="right", va="center", color='white',
             bbox=dict(boxstyle="square", color='tomato'))
    ax1.text(0.97, 0.1, "2", size=6, ha="right", va="center", color='white',
             bbox=dict(boxstyle="square", color='limegreen'))
    ax1.text(0.99, 0.1, "3", size=6, ha="right", va="center", color='white',
             bbox=dict(boxstyle="square", color='deepskyblue'))

    # riboseq bar plots
    ax2 = plt.subplot(gs[1])
    label = 'Ribo-Seq count'
    if read_length:
        label = 'Ribo-Seq count ({}-mer)'.format(read_length)
    ax2.set_ylabel(label, fontdict=font_xsmall, labelpad=10)

    # rna coverage if available
    ax_rna = None
    if rna_counts:
        ax_rna = ax2.twinx()
        ax_rna.set_ylabel('RNA-Seq count', fontdict=font_xsmall, labelpad=10)
        ax_rna.bar(rna_counts.keys(), rna_counts.values(), facecolor='#e6e6e6', edgecolor='#e6e6e6', label='RNA', linewidth=0.5)
        ax_rna.set_zorder(1)

    frame_counts = {1: {}, 2: {}, 3: {}}
    for k, v in ribo_counts.iteritems():
        for fr in (1, 2, 3):
            if v[fr] > 0:
                frame_counts[fr][k] = v[fr]
                break

    cnts = []
    [cnts.extend(item.values()) for item in frame_counts.values()]
    y_max = int(round(max(cnts) * 1.2))
    ax2.set_ylim(0.0, y_max)
    ax2.set_zorder(2)
    ax2.patch.set_facecolor('none')

    for frame, color in ((1, 'tomato'), (2, 'limegreen'), (3, 'deepskyblue')):
        if read_offset:
            x_vals = [pos + read_offset for pos in frame_counts[frame].keys()]
        else:
            x_vals = frame_counts[frame].keys()
        ax2.bar(x_vals, frame_counts[frame].values(), color=color, facecolor=color, edgecolor=color, linewidth=0.5)

    ax3 = plt.subplot(gs[2], axisbg='white')
    ax3.text(0.90, 0.1, "start codon", size=6, ha="right", va="center", color='#555555',
             bbox=dict(boxstyle="square", facecolor='white', edgecolor='#c6c6c6'))
    ax3.text(0.99, 0.1, "stop codon", size=6, ha="right", va="center", color='white',
             bbox=dict(boxstyle="square", color='#777777'))

    ax4 = plt.subplot(gs[3], sharex=ax2, axisbg='#c6c6c6')
    ax5 = plt.subplot(gs[4], sharex=ax2, axisbg='#c6c6c6')
    ax6 = plt.subplot(gs[5], sharex=ax2, axisbg='#c6c6c6')

    for axis in (ax1, ax3):
        axis.tick_params(top=False, left=False, right=False, bottom=False, labeltop=False,
                         labelleft=False, labelright=False, labelbottom=False)
        set_axis_color(axis, 'white')

    axes = [ax2]
    if ax_rna:
        axes.append(ax_rna)

    fp = FontProperties(size='5')
    for axis in axes:
        set_axis_color(axis, '#f7f7f7')
        for item in (axis.get_xticklabels() + axis.get_yticklabels()):
            item.set_fontproperties(fp)
            item.set_color('#555555')

    for axis, frame in ((ax4, 1), (ax5, 2), (ax6, 3)):
        set_axis_color(axis, '#C6C6C6')
        for item in (axis.get_xticklabels()):
            item.set_fontproperties(fp)
            item.set_color('#555555')
        axis.set_ylim(0, 0.2)
        axis.set_xlim(0, transcript_length)
        starts = [(item, 1) for item in start_stops[frame]['starts']]
        stops = [(item, 1) for item in start_stops[frame]['stops']]
        start_colors = ['#ffffff' for item in starts]
        axis.broken_barh(starts, (0.11, 0.2),
                         facecolors=start_colors, edgecolors=start_colors, label='start', zorder=5, linewidth=0.5)
        stop_colors = ['#777777' for item in stops]
        axis.broken_barh(stops, (0, 0.2), facecolors=stop_colors,
                         edgecolors=stop_colors, label='stop', zorder=5, linewidth=0.5)
        axis.set_ylabel('{}'.format(frame),
                        fontdict={'family': 'sans-serif', 'color': '#555555',
                                  'weight': 'normal', 'size': '6'},
                        rotation='horizontal', labelpad=10, verticalalignment='center')
        axis.tick_params(top=False, left=False, right=False, labeltop=False,
                         labelleft=False, labelright=False, direction='out')
    plt.xlabel('Transcript length ({} nt)'.format(transcript_length), fontdict=font_xsmall, labelpad=10)
    plt.title('Transcript {}'.format(transcript_name),
              fontdict={'family': 'sans-serif', 'color': '#222222',
                        'weight': 'normal', 'size': 'x-small'}, y=13.0)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    plt.savefig(os.path.join(output_path, 'riboplot.svg'))
    plt.savefig(os.path.join(output_path, 'riboplot.png'), dpi=300)

    with open(os.path.join(CONFIG.PKG_DATA_DIR, 'riboplot.html')) as g, open(os.path.join(output_path, html_file), 'w') as h:
        h.write(g.read().format(transcript_name=transcript_name))

    css_dir = os.path.join(output_path, 'css')
    if not os.path.exists(css_dir):
        os.mkdir(css_dir)

    css_data_dir = os.path.join(CONFIG.PKG_DATA_DIR, 'css')
    for fname in os.listdir(css_data_dir):
        shutil.copy(os.path.join(css_data_dir, fname), os.path.join(output_path, 'css', fname))


def create_parser():
    """Argument parser. """
    parser = argparse.ArgumentParser(
        prog='riboplot.py', description='Plot and output read counts for a single transcript')

    required = parser.add_argument_group('required arguments')
    required.add_argument('-b', '--ribo_file', help='Ribo-Seq alignment file in BAM format', required=True)
    required.add_argument('-f', '--transcriptome_fasta', help='FASTA format file of the transcriptome', required=True)
    required.add_argument('-t', '--transcript_name', help='Transcript name', metavar='TEXT', required=True)

    # plot function - optional arguments
    parser.add_argument('-n', '--rna_file', help='RNA-Seq alignment file (BAM)')
    parser.add_argument('-l', '--read_length', help='Read length to consider (default: %(default)s)',
                        metavar='INTEGER', type=int)
    parser.add_argument('-s', '--read_offset', help='Read offset (default: %(default)s)',
                        metavar='INTEGER', type=int, default=0)
    parser.add_argument('-m', '--html_file', help='Output file for results (HTML)', default='riboplot.html')
    parser.add_argument('-o', '--output_path', help='Files are saved in this directory', default='output')
    parser.add_argument('-d', '--debug', help='Flag. Produce debug output', action='store_true')

    return parser


def main(args):
    """Main program"""
    (ribo_file, rna_file, transcript_name, transcriptome_fasta, read_length,
     read_offset, output_path, html_file) = (
         args.ribo_file, args.rna_file, args.transcript_name, args.transcriptome_fasta,
         args.read_length, args.read_offset, args.output_path, args.html_file)

    log.debug('Supplied arguments\n{}'.format(
        '\n'.join(['{:<20}: {}'.format(k, v) for k, v in vars(args).items()])))

    # error messages (simple format) are written to html file
    fh = logging.FileHandler(html_file)
    fh.setLevel(logging.ERROR)
    fh.setFormatter(ErrorLogFormatter('%(message)s'))
    log.addHandler(fh)

    log.info('Checking if required arguments are valid...')
    ribocore.check_required_arguments(ribo_file=ribo_file, transcriptome_fasta=transcriptome_fasta,
                                      transcript_name=transcript_name)
    log.info('Done')

    log.info('Checking if optional arguments are valid...')
    ribocore.check_optional_arguments(ribo_file=ribo_file, read_length=read_length, read_offset=read_offset,
                                      rna_file=rna_file)
    log.info('Done')

    log.info('Get ribo-seq read counts and total reads in Ribo-Seq...')
    bam_fileobj = pysam.AlignmentFile(ribo_file, 'rb')
    ribo_counts, total_reads = ribocore.get_ribo_counts(bam_fileobj, transcript_name, read_length)
    bam_fileobj.close()

    if not ribo_counts:
        msg = ('No RiboSeq read counts for transcript {}. No plot will be '
               'generated!'.format(transcript_name))
        log.error(msg)
        raise ribocore.RiboPlotError(msg)
    else:
        log.info('Get RNA counts for the given transcript...')
        mrna_counts = {}
        if rna_file:
            try:
                mrna_counts = get_rna_counts(rna_file, transcript_name)
            except OSError as e:
                log.error(e)
                raise

            if not mrna_counts:
                log.warn('No RNA counts for this transcript from the given RNA Seq file. '
                         'RNA-Seq coverage will not be generated')
        else:
            log.debug('No RNA-Seq data provided. Not generating coverage')

        log.info('Get sequence and length of the given transcripts from FASTA file...')
        fasta_records = ribocore.get_fasta_records(transcriptome_fasta, [transcript_name])
        transcript_seq, transcript_length = (fasta_records[transcript_name]['sequence'],
                                             fasta_records[transcript_name]['length'])

        log.info('Get start/stop positions in transcript sequence (3 frames)...')
        codon_positions = get_start_stops(transcript_seq)

        if not os.path.exists(output_path):
            os.mkdir(output_path)

        log.info('Writing RiboSeq read counts for {}'.format(transcript_name))
        with open(os.path.join(output_path, 'RiboCounts.csv'), 'w') as f:
            f.write('"Position","Frame 1","Frame 2","Frame 3"\n')

            for pos in range(1, transcript_length + 1):
                if pos in ribo_counts:
                    f.write('{0},{1},{2},{3}\n'.format(
                        pos, ribo_counts[pos][1], ribo_counts[pos][2], ribo_counts[pos][3]))
                else:
                    f.write('{0},{1},{2},{3}\n'.format(pos, 0, 0, 0))

        log.info('Generating RiboPlot...')
        plot_profile(ribo_counts, transcript_name, transcript_length,
                     codon_positions, read_length, read_offset, mrna_counts,
                     html_file=args.html_file, output_path=args.output_path)
    log.info('Finished!')


def run():
    """Run program"""
    parsed = create_parser()
    args = parsed.parse_args()
    main(args)


if __name__ == '__main__':
    run()

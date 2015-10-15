# -*- coding: utf-8 -*-
import sys

# check dependencies
module_errors = {}

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
        transcript_sequence = transcript_sequence.upper()  # for comparison with codons below
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


def set_axis_color(axis, color, alpha=None):
    """Sets the spine color of all sides of an axis (top, right, bottom, left)."""
    for side in ('top', 'right',  'bottom', 'left'):
        spine = axis.spines[side]
        spine.set_color(color)
        if alpha is not None:
            spine.set_alpha(alpha)


def get_color_palette(scheme):
    """Return colors for a given scheme. Default colors are returned for an item
    if undefined in scheme.

    """
    color_schemes = {
        'default': {
            'frames': ['tomato', 'limegreen', 'deepskyblue'], 'background': '#fafafa',
            'color': '#616161', 'ticks': '#757575', 'start': '#ffffff', 'stop': '#919191',
            'rna': '#e0e0e0', 'axis': '#e0e0e0', 'grey': '#bdbdbd'
        },
        'colorbrewer': {
            'frames': ['#fc8d62', '#66c2a5', '#8da0cb']
        },
        'rgb': {
            'frames': ['red', 'green', 'blue']
        },
        'greyorfs': {}
    }

    colors = {}
    for k, v in color_schemes['default'].items():
        try:
            vals = color_schemes[scheme][k]
        except KeyError:
            vals = v
        colors[k] = vals
    return colors


def plot_profile(ribo_counts, transcript_name, transcript_length,
                 start_stops, read_length=None, read_offset=None, rna_counts=None,
                 color_scheme='default', html_file='index.html', output_path='output'):
    """Plot read counts (in all 3 frames) and RNA coverage if provided for a
    single transcript.

    """
    colors = get_color_palette(scheme=color_scheme)
    gs = gridspec.GridSpec(3, 1, height_ratios=[6, 1.3, 0.5], hspace=0.35)
    font_axis = {'family': 'sans-serif', 'color': colors['color'], 'weight': 'bold', 'size': 7}

    # riboseq bar plots
    gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    ax2 = plt.subplot(gs2[0])
    label = 'Ribo-Seq count'

    if read_length:
        label = 'Ribo-Seq count ({}-mer)'.format(read_length)
    ax2.set_ylabel(label, fontdict=font_axis, labelpad=10)

    # rna coverage if available
    ax_rna = None
    if rna_counts:
        ax_rna = ax2.twinx()
        ax_rna.set_ylabel('RNA-Seq count', fontdict=font_axis, labelpad=10)
        ax_rna.bar(rna_counts.keys(), rna_counts.values(), facecolor=colors['rna'],
                   edgecolor=colors['rna'], label='RNA', linewidth=0)
        ax_rna.set_zorder(1)

    frame_counts = {1: {}, 2: {}, 3: {}}
    for k, v in ribo_counts.iteritems():
        for fr in (1, 2, 3):
            if v[fr] > 0:
                frame_counts[fr][k] = v[fr]
                break

    cnts = []
    [cnts.extend(item.values()) for item in frame_counts.values()]
    y_max = float(max(cnts) * 1.25)
    ax2.set_ylim(0.0, y_max)
    ax2.set_zorder(2)
    ax2.patch.set_facecolor('none')

    for frame in (1, 2, 3):
        color = colors['frames'][frame - 1]
        if read_offset:
            x_vals = [pos + read_offset for pos in frame_counts[frame].keys()]
        else:
            x_vals = frame_counts[frame].keys()
        ax2.bar(x_vals, frame_counts[frame].values(), color=color, facecolor=color, edgecolor=color, linewidth=0)

    # ORF architecture
    gs3 = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1], hspace=0.1)
    if color_scheme == 'greyorfs':
        axisbg = [colors['grey'] for i in range(3)]
    else:
        axisbg = colors['frames']

    ax4 = plt.subplot(gs3[0], sharex=ax2, axisbg=axisbg[0])
    ax5 = plt.subplot(gs3[1], sharex=ax2, axisbg=axisbg[1])
    ax6 = plt.subplot(gs3[2], sharex=ax2, axisbg=axisbg[2])
    ax6.set_xlabel('Transcript length ({} nt)'.format(transcript_length), fontdict=font_axis, labelpad=6)

    # Legend
    gs4 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[2], hspace=0.1)
    ax7 = plt.subplot(gs4[0], axisbg=colors['background'])
    set_axis_color(ax7, colors['background'])

    ax7.text(0.03, 0.1, "Codons: ", size=5, ha='center', va='center', color=colors['color'],
             fontdict={'weight': 'bold'})
    ax7.text(0.10, 0.1, "START", size=5, ha="center", va="center", color=colors['color'],
             bbox=dict(boxstyle="square", facecolor=colors['start'], edgecolor=colors['color'], linewidth=0.2))
    ax7.text(0.15, 0.1, "STOP", size=5, ha="center", va="center", color='white',
             bbox=dict(boxstyle="square", color=colors['stop']))
    ax7.text(0.22, 0.1, "Frames: ", size=5, ha='center', va='center', color=colors['color'],
             fontdict={'weight': 'bold'})
    ax7.text(0.27, 0.1, "1", size=5, ha="center", va="center", color='white',
             bbox=dict(boxstyle="square", color=colors['frames'][0]))
    ax7.text(0.29, 0.1, "2", size=5, ha="center", va="center", color='white',
             bbox=dict(boxstyle="square", color=colors['frames'][1]))
    ax7.text(0.31, 0.1, "3", size=5, ha="center", va="center", color='white',
             bbox=dict(boxstyle="square", color=colors['frames'][2]))

    # No ticks or labels for ORF 1, 2 and Legend
    for axis in (ax4, ax5, ax7):
        axis.tick_params(top=False, left=False, right=False, bottom=False, labeltop=False,
                         labelleft=False, labelright=False, labelbottom=False)

    axes = [ax2]
    if ax_rna:
        axes.append(ax_rna)

    fp = FontProperties(size='5')
    for axis in axes:
        set_axis_color(axis, colors['axis'])
        axis.tick_params(colors=colors['ticks'])
        for item in (axis.get_xticklabels() + axis.get_yticklabels()):
            item.set_fontproperties(fp)
            item.set_color(colors['color'])

    for axis, frame in ((ax4, 1), (ax5, 2), (ax6, 3)):
        if color_scheme == 'greyorfs':
            color = colors['grey']
        else:
            color = colors['frames'][frame - 1]
        set_axis_color(axis, color, alpha=0.05)
        axis.patch.set_alpha(0.3)  # opacity of ORF architecture
        for item in (axis.get_xticklabels()):
            item.set_fontproperties(fp)
            item.set_color(colors['color'])
        axis.set_ylim(0, 0.2)
        axis.set_xlim(0, transcript_length)
        starts = [(item, 1) for item in start_stops[frame]['starts']]
        stops = [(item, 1) for item in start_stops[frame]['stops']]
        start_colors = [colors['start'] for item in starts]
        axis.broken_barh(starts, (0.11, 0.2),
                         facecolors=start_colors, edgecolors=start_colors, label='start', zorder=5, linewidth=0)
        stop_colors = [colors['stop'] for item in stops]
        axis.broken_barh(stops, (0, 0.2), facecolors=stop_colors,
                         edgecolors=stop_colors, label='stop', zorder=5, linewidth=0)
        axis.set_ylabel('{}'.format(frame),
                        fontdict={'family': 'sans-serif', 'color': colors['color'],
                                  'weight': 'normal', 'size': '6'},
                        rotation='horizontal', labelpad=10, verticalalignment='center')
        axis.tick_params(top=False, left=False, right=False, labeltop=False,
                         labelleft=False, labelright=False, direction='out', colors=colors['ticks'])
    plt.title('{}'.format(transcript_name),
              fontdict={'family': 'sans-serif', 'color': colors['color'],
                        'weight': 'bold', 'size': 8, 'y': 20})
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    plt.savefig(os.path.join(output_path, 'riboplot.svg'), facecolor=colors['background'])
    plt.savefig(os.path.join(output_path, 'riboplot.png'), dpi=600, facecolor=colors['background'])

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
    parser.add_argument('-c', '--color_scheme', help='Color scheme to use (default: %(default)s)',
                        choices=['default', 'colorbrewer', 'rgb', 'greyorfs'], default='default')
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
    ribocore.check_required_arguments(
        ribo_file=ribo_file, transcriptome_fasta=transcriptome_fasta, transcript_name=transcript_name)
    log.info('Done')

    log.info('Checking if optional arguments are valid...')
    ribocore.check_optional_arguments(
        ribo_file=ribo_file, read_length=read_length, read_offset=read_offset, rna_file=rna_file)
    log.info('Done')

    log.info('Get sequence and length of the given transcript from FASTA file...')
    record = ribocore.get_fasta_record(transcriptome_fasta, transcript_name)
    transcript_sequence = record[transcript_name]
    transcript_length = len(transcript_sequence)

    log.info('Get ribo-seq read counts and total reads in Ribo-Seq...')
    with ribocore.open_pysam_file(fname=ribo_file, ftype='bam') as bam_fileobj:
        ribo_counts, total_reads = ribocore.get_ribo_counts(
            ribo_fileobj=bam_fileobj, transcript_name=transcript_name, read_length=read_length)

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

        log.info('Get start/stop positions in transcript sequence (3 frames)...')
        codon_positions = get_start_stops(transcript_sequence)

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
                     color_scheme=args.color_scheme,
                     html_file=args.html_file, output_path=args.output_path)
    log.info('Finished!')


def run():
    """Run program"""
    parsed = create_parser()
    args = parsed.parse_args()
    main(args)


if __name__ == '__main__':
    run()

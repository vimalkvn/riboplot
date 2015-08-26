.. _usage:

=====
Usage
=====

RiboPlot
--------
Plot and output Ribo-Seq read counts of a single transcript from an alignment file (sorted BAM).

Parameters
..........

1. Ribo-Seq alignment file (Sorted BAM file)
++++++++++++++++++++++++++++++++++++++++++++
A Bowtie 1 output (BAM) from an alignment of Ribo-Seq data to the transcriptome. This BAM
file should be sorted. This can be done using one of the following methods.

1. RiboGalaxy_ -> Sort Data -> Sort BAM dataset.
2. ``samtools sort input.bam inputsorted``

2. Transcriptome (FASTA)
++++++++++++++++++++++++
A FASTA format file with sequences of the transcripts.

3. Name of the transcript to plot (Text)
++++++++++++++++++++++++++++++++++++++++
The name of the transcript to plot **should** match the name in the transcriptome (FASTA)
and the Ribo-Seq/RNA-Seq alignment (BAM).

4. RNA coverage [optional] (Sorted BAM file)
++++++++++++++++++++++++++++++++++++++++++++
If you have RNA-Seq data (sorted BAM), you can select the option to plot RNA coverage.

5. Read lengths to consider [Optional] (Integer - 0 or greater)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If this option is provided, only Ribo-Seq data of the given length is considered.

6. Offset [optional] (Integer - 0 or greater)
+++++++++++++++++++++++++++++++++++++++++++++
If this option is provided, this offset is added to the read alignment positions.

Output
......
1. Plots (PNG and SVG)
++++++++++++++++++++++
Ribo-Seq read counts as a bar plot in 3 frames (color codes: 1: red, 2: green, 3: blue)

RNA coverage as a gray background (if the RNA coverage option was selected).

The open reading frame architecture appears below the plot with start (ATG) and stop codons ('TAA', 'TAG', 'TGA') in all 3 frames.

The color codes are start (white) and stop (dark gray).

.. image:: ../images/riboplot.png

2. RiboSeq read counts (CSV)
++++++++++++++++++++++++++++
In 3 frames for each position in the transcript.


Command line
............
``riboplot`` can also be run on the command line. The usage is ::

    usage: riboplot [-h] -b RIBO_FILE -f TRANSCRIPTOME_FASTA -t TEXT
                    [-n RNA_FILE] [-l INTEGER] [-s INTEGER] [-m HTML_FILE]
                    [-o OUTPUT_PATH] [-d]

    Plot and output read counts for a single transcript

    optional arguments:
    -h, --help
        show this help message and exit

    -n RNA_FILE, --rna_file RNA_FILE
        RNA-Seq alignment file (BAM)

    -l INTEGER, --read_length INTEGER
        Read length to consider (default: None)

    -s INTEGER, --read_offset INTEGER
        Read offset (default: 0)

    -m HTML_FILE, --html_file HTML_FILE
        Output file for results (HTML)

    -o OUTPUT_PATH, --output_path OUTPUT_PATH
        Files are saved in this directory

    -d, --debug
        Flag. Produce debug output

    required arguments:
    -b RIBO_FILE, --ribo_file RIBO_FILE
        Ribo-Seq alignment file in BAM format

    -f TRANSCRIPTOME_FASTA, --transcriptome_fasta TRANSCRIPTOME_FASTA
        FASTA format file of the transcriptome

    -t TEXT, --transcript_name TEXT
        Transcript name

RiboCount
---------
Output read counts for all transcripts in an alignment.

Parameters
..........
1. Ribo-Seq alignment file (Sorted BAM file)
++++++++++++++++++++++++++++++++++++++++++++
A Bowtie 1 output (BAM) from an alignment of Ribo-Seq data to the transcriptome. This BAM
file should be sorted. This can be done using one of the following methods.

1. RiboGalaxy_ -> Sort Data -> Sort BAM dataset.
2. ``samtools sort input.bam inputsorted``

2. Transcriptome (FASTA)
++++++++++++++++++++++++
A FASTA format file with sequences of the transcripts.

3. Read lengths to consider [optional] (Integer - 0 or greater)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If this option is provided, only Ribo-Seq data of the given length is considered.

4. Offset [optional] (Integer - 0 or greater)
+++++++++++++++++++++++++++++++++++++++++++++
If this option is provided, this offset is added to the read alignment positions.

Output
......
Read counts for all transcripts in the alignment (ZIP)
++++++++++++++++++++++++++++++++++++++++++++++++++++++
The output file ``ribocount_output.zip`` should first be uncompressed. This will generate
a folder called ``ribocount_output``. Open ``index.html`` in a web browser to view the results of ribocount.

Total reads for each transcript will be displayed in a table along with the name of the transcript and a link
to the CSV file containing the read counts in 3 frames for each position in the transcript.

.. image:: ../images/ribocount.png

Command line
............
``ribocount`` can also be run on the command line. The usage is ::

    usage: ribocount [-h] -b RIBO_FILE -f TRANSCRIPTOME_FASTA [-l INTEGER]
    [-s INTEGER] [-m HTML_FILE] [-o OUTPUT_PATH] [-d]

    Output read counts for all transcripts

    optional arguments:

        -h, --help            show this help message and exit

        -l INTEGER, --read_length INTEGER
            Read length to consider (default: None)

        -s INTEGER, --read_offset INTEGER
            Read offset (default: 0)

        -m HTML_FILE, --html_file HTML_FILE

            Output file for results (HTML)

        -o OUTPUT_PATH, --output_path OUTPUT_PATH
            Files are saved in this directory

        -d, --debug           Flag. Produce debug output

    required arguments:

        -b RIBO_FILE, --ribo_file RIBO_FILE
            Ribo-Seq alignment file in BAM format

        -f TRANSCRIPTOME_FASTA, --transcriptome_fasta TRANSCRIPTOME_FASTA
            FASTA format file of the transcriptome

.. links
.. _RiboGalaxy: http://ribogalaxy.ucc.ie


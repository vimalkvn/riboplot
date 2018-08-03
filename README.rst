========
RiboPlot
========

.. image:: https://img.shields.io/pypi/v/riboplot.svg
   :target: https://pypi.python.org/pypi/riboplot

RiboPlot includes programs to plot and output Ribo-Seq read counts from an alignment file (BAM format).

There are two programs in the package:

**riboplot**
  Plot and output read counts (csv) for a single transcript.

**ribocount**
  Output read counts for all transcripts in an alignment.

Quickstart
----------

Installation
............
riboplot can be installed easily with Conda_ using the command::

  conda create -n riboplot riboplot

This is the **recommended** method. Please check the
installation_ section of the documentation for alternative
installation methods.

Get sample data
...............
For convenience, this repository includes sample Ribo-Seq and
RNA-seq data (from `Bazzini et.al`_) that have been processed
through the RiboGalaxy pipeline - removal of adapters and rRNA
followed by aligning the reads to the transcriptome. The resulting
alignment (SAM) file is then converted to BAM format which is
then sorted. If you would like to learn more about this workflow,
please check the `published page`_ on RiboGalaxy.

To get sample data, clone this repository::

  git clone https://github.com/vimalkvn/riboplot

Activate the conda environment created above::

  source activate riboplot
	
Change to the directory containing sample data::

  cd tests/data

Plot and output Ribo-Seq read counts for a single transcript (riboplot)
.......................................................................
::

  riboplot -b 5hRPFsorted.bam -n 5hmRNAsorted.bam \
  -f zebrafish.fna -t "gi|41055123|ref|NM_201172.1|"

Here::

  -b Ribo-Seq alignment file in BAM format

  -n RNA-Seq alignment file (BAM) 
  (optional, for including RNA coverage in the plot)

  -f FASTA format file of the transcriptome

  -t Name of the transcript as in the FASTA file

Outputs will be in the *output* folder. Plot is saved as
*riboplot.png* and the read counts will be in *RiboCounts.csv*.

.. image:: docs/images/riboplot_default.png
   :scale: 60 %
   :alt: RiboPlot of a single transcript

Output read counts for all transcripts (ribocount)
..................................................
::

  ribocount -b 5hRPFsorted.bam -f zebrafish.fna

When the run is complete, a zip file called
*ribocount_output.zip* will be created in the *output* folder.
Uncompress this file and open *index.html* in a browser::

   cd output/
   unzip ribocount_output.zip 
   cd ribocount_output/
   xdg-open index.html 

.. image:: docs/images/ribocount.png
   :scale: 60 %
   :alt: RiboCount of all transcripts

   
For detailed usage, please check the online documentation
http://pythonhosted.org/riboplot

PyPI: https://pypi.python.org/pypi/riboplot

Github Repository: https://github.com/vimalkvn/riboplot

Online: http://ribogalaxy.ucc.ie (under RiboSeq Analysis/RiboPlot)

Free software: GPL license.

.. Links

.. _Bazzini et.al: https://www.ncbi.nlm.nih.gov/pubmed/24705786
.. _Conda: https://conda.io
.. _installation: https://pythonhosted.org/riboplot/installation.html
.. _published page: https://ribogalaxy.ucc.ie/u/vimalkumarvelayudhan/p/using-riboplot-and-ribocount



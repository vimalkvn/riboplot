.. :changelog:

History
=======
0.3 (2015-12-10)
----------------
* Support for multiple read lengths and corresponding offsets.
* Optional arguments check is now split into multiple steps.
* Remove the calculation of read lengths present in the BAM alignment. 

0.2.5 (2015-11-30)
------------------
* Bugfix: Use known versions of libraries as dependencies in setup.py.
  Otherwise, pip installs latest versions which haven't been tested.
* Remove Codons label, change START to AUG [riboplot].

0.2.4 (2015-11-25)
------------------
* Bugfix: ribocount now returns correct read counts if an offset is provided.
* Bugfix: Don't include read counts in the longest ORF start or stop positions
  i.e., only include reads upstream or downstream of the start or stop positions.

0.2.3 (2015-11-24)
------------------
* Bugfix: Recalculate read frame positions after applying offset.
* Add sequence to CSV output.

0.2.2 (2015-10-29)
------------------
* Use default linewidth (riboplot), Minor changes to default plot style.
* Use 'Agg' as the default matplotlib backend (prevent $DISPLAY errors).
* Use smaller images for the help section.
* Fix typo in HISTORY.rst.
* Update Github repository URL.

0.2.1 (2015-10-16)
------------------
Fix: Add mock again to setup requirements. Matplolib install under Galaxy fails otherwise.

0.2.0 (2015-10-15)
------------------
**riboplot**

* Support for color schemes - options: default, colorbrewer, rgb, greyorfs.
* Plot changes - colored, light background for ORF architecture, legend moved down.
* Image DPI change 300 |srarr| 600.
* Sample scripts to plot multiple transcripts and multiple color schemes.

**ribocount**

* Ability to output read counts in 5'(``-v``) and 3'(``-r``) leader region.
       
Bugfixes

1. Report error if a bam index could not be generated (ex. bam sorted using 
   reference names instead of chromosome coordinates).
2. No start/stop values when transcript sequences are in lower case.

Tests

* Split into three different files.
* Updated test configuration.

Minor

* XML wrapper changes - remove 'leader', add help text and output label change [ribocount].
* Remove unused pysam import [riboplot].
* Replace doc includes with sym links (breaks Galaxy toolshed otherwise).
* Updated CSS/Table styles (removed number column from table).
* y limit is now 1.25x [riboplot].
* Removed line width on bars [riboplot].

0.1.1 (2015-09-08)
------------------
* Remove invalid Python 3 version (whl).
* Remove unused test data files.
* Remove strict dependence on matplotlib and pysam versions.
* Cleanup MANIFEST.

0.1.0 (2015-08-24)
------------------
* First release on PyPI.

.. substitutions  
.. |srarr|    unicode:: U+02192 .. RIGHTWARDS ARROW

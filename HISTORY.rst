.. :changelog:

History
=======
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

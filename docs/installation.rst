============
Installation
============
RiboPlot is available as a standalone program or as a wrapper for `Galaxy <http://galaxyproject.org>`_.
Please follow the steps below depending on how you would like to use it.

.. note::

    1. RNA coverage plot requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

       This release of riboplot has been tested with bedtools versions ``2.17.0`` and ``2.24``.

       On Ubuntu and derivatives, bedtools can be installed from the repositories using::

          sudo apt-get install bedtools

    2. Matplotlib installation requires freetype headers and libraries installed
       (``libfreetype6-dev`` on Ubuntu 14.04).


Standalone
----------
At the command line::

    $ pip install riboplot

Or, if you have virtualenv and virtualenvwrapper installed::

    $ mkvirtualenv riboplot
    $ pip install riboplot

.. note::

    If you are using easy_install, please use a recent version of setuptools.

    On Ubuntu 12.04, you might need to tell virtualenv to use setuptools instead of
    distribute (default) like this::

        mkvirtualenv --setuptools riboplot

Galaxy wrapper
--------------
Login with an administrator account on your Galaxy instance. Search for **riboplot** 
from the Galaxy Main Toolshed and install it on your instance. Dependencies will be installed automatically. See note on freetype requirement for Matplotlib above. 



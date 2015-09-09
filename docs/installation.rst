============
Installation
============
RiboPlot is available as a standalone program or as a wrapper for `Galaxy <http://galaxyproject.org>`_.
Please follow the steps below depending on how you would like to use it.

.. note::

    RNA coverage plot requires `bedtools <https://github.com/arq5x/bedtools2>`_ to be installed.

    This release of riboplot has been tested with bedtools version ``2.17.0``.

    On Ubuntu and derivatives, bedtools can be installed from the repositories using::

        sudo apt-get install bedtools

Standalone
----------
At the command line::

    $ pip install riboplot

Or, if you have virtualenvwrapper installed::

    $ mkvirtualenv riboplot
    $ pip install riboplot

.. note::

    If you are using easy_install, please use a recent version of setuptools.

Galaxy wrapper
--------------
Login with an administrator account on your Galaxy instance. Search for **riboplot** from the Galaxy Main Toolshed and install it on your instance. Dependencies (except bedtools, see note below) will be installed automatically.

